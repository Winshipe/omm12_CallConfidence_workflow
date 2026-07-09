"""
rules/variant_calling.smk — align blended reads to the reference and call variants
====================================================================================

This module replaces the breseq-based caller with a BWA-MEM + GATK
HaplotypeCaller pipeline, which is the industry-standard approach for
short-read variant calling.

Overview of steps
-----------------
1. bwa_index         — index each unmutated reference FASTA so BWA can
                       align reads against it (only run once per reference,
                       not once per scenario).

2. bwa_align         — align the blended paired-end reads for a scenario to
                       the (merged) reference with BWA-MEM, sort the result,
                       and mark duplicate read pairs with samtools.

3. gatk_haplotype_caller — call variants with GATK HaplotypeCaller in
                           EMIT_ALL_SITES or standard mode, producing a VCF.

Why the unmutated reference?
-----------------------------
The pipeline measures whether the variant caller can *discover* the mutations
that were artificially introduced.  Aligning to the original, unmutated
sequence means every called variant is a genuine discovery (or a false
positive), not a pre-existing difference.

Multiple references
-------------------
When a scenario mixes reads from more than one reference genome, all unmutated
FASTAs are concatenated into a single combined reference before indexing and
alignment.  GATK and BWA both handle multi-contig references natively, so no
special flags are needed.

Output files (per scenario × replicate)
----------------------------------------
  results/variant_calling/{scenario}/{replicate}/output.vcf.gz        final variant calls
  results/variant_calling/{scenario}/{replicate}/output.vcf.gz.tbi    VCF index (tabix)
  results/variant_calling/{scenario}/{replicate}/aligned.bam          sorted, deduplicated BAM
  results/variant_calling/{scenario}/{replicate}/aligned.bam.bai      BAM index

Configuration keys read from config["variant_calling"]
-------------------------------------------------------
  threads          : int   — CPU threads for BWA and GATK (default 8)
  bwa_extra_flags  : str   — optional extra flags passed to bwa mem
  gatk_extra_flags : str   — optional extra flags passed to HaplotypeCaller
  min_base_quality : int   — minimum base quality score for GATK (default 20)
"""

import os


# ---------------------------------------------------------------------------
# Helper: collect the unmutated reference FASTAs for a given scenario
# ---------------------------------------------------------------------------

def scenario_references(wildcards):
    """
    Return a list of unmutated FASTA paths for every reference that
    contributes reads in this scenario.

    Duplicates are removed while preserving the order in which references
    first appear in the config.  This list is used both to build the merged
    reference and as an explicit Snakemake input dependency.
    """
    contribs = config["scenarios"][wildcards.scenario]
    seen = {}
    for c in contribs:
        ref = c["ref_id"]
        if ref not in seen:
            seen[ref] = config["references"][ref]
    return list(seen.values())


def merged_ref_path(wildcards):
    """
    Return the path of the merged (concatenated) reference FASTA for a
    scenario.  This is the single file that BWA and GATK will use.
    """
    return f"results/variant_calling/{wildcards.scenario}/{wildcards.replicate}/reference.fasta"


# ---------------------------------------------------------------------------
# Rule 1 — merge references and build BWA index
# ---------------------------------------------------------------------------

rule bwa_index:
    """
    Concatenate all per-scenario reference FASTAs into one combined FASTA,
    then build a BWA index and a samtools FASTA index (.fai) and a GATK
    sequence dictionary (.dict).

    Why merge?
    ----------
    BWA and GATK both expect a single reference file.  Concatenating the
    FASTAs gives each contig a unique name, so reads from different source
    genomes are mapped to the correct contig without any naming conflicts.

    Why the .fai and .dict?
    -----------------------
    GATK requires both a samtools FASTA index (.fai) and a sequence
    dictionary (.dict) to quickly look up contig lengths and offsets.
    """
    input:
        # All unmutated FASTAs for this scenario (resolved by the helper above)
        refs=scenario_references,
    output:
        # The merged reference FASTA
        ref="results/variant_calling/{scenario}/{replicate}/reference.fasta",
        # BWA index files (BWA creates these automatically alongside the FASTA)
        bwt="results/variant_calling/{scenario}/{replicate}/reference.fasta.bwt",
        pac="results/variant_calling/{scenario}/{replicate}/reference.fasta.pac",
        ann="results/variant_calling/{scenario}/{replicate}/reference.fasta.ann",
        amb="results/variant_calling/{scenario}/{replicate}/reference.fasta.amb",
        sa ="results/variant_calling/{scenario}/{replicate}/reference.fasta.sa",
        # samtools FASTA index (required by GATK)
        fai="results/variant_calling/{scenario}/{replicate}/reference.fasta.fai",
        # GATK sequence dictionary (required by GATK)
        dic="results/variant_calling/{scenario}/{replicate}/reference.dict",
    log:
        "logs/variant_calling/{scenario}/{replicate}.index.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        # ── Step 1: merge all reference FASTAs into one file ──────────────
        # 'cat' simply concatenates the files; because each FASTA has its
        # own '>' header lines, the result is a valid multi-contig FASTA.
        cat {input.refs} > {output.ref} 2>> {log}

        # ── Step 2: build BWA index ───────────────────────────────────────
        # BWA index creates five companion files (.bwt, .pac, .ann, .amb, .sa)
        # alongside the reference FASTA.  These are look-up structures BWA
        # uses to rapidly find where each read aligns.
        bwa index {output.ref} >> {log} 2>&1

        # ── Step 3: samtools FASTA index (.fai) ──────────────────────────
        # Needed by GATK to quickly look up any genomic position.
        samtools faidx {output.ref} >> {log} 2>&1

        # ── Step 4: GATK sequence dictionary (.dict) ──────────────────────
        # GATK requires this file to know the names and lengths of all contigs.
        gatk CreateSequenceDictionary \
            --REFERENCE {output.ref} \
            --OUTPUT    {output.dic} \
            >> {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule 2 — align reads with BWA-MEM and mark duplicates
# ---------------------------------------------------------------------------

rule bwa_align:
    """
    Align blended paired-end reads to the merged reference with BWA-MEM,
    then sort and mark PCR duplicates.

    Steps performed in a single shell block to avoid writing large
    intermediate files to disk:

    1. bwa mem    — align reads; output is unsorted SAM written to stdout
    2. samtools sort — sort alignments by coordinate; output is a BAM file
    3. samtools markdup — flag duplicate read pairs so GATK can ignore them
       (duplicates arise from PCR amplification and do not represent
        independent observations of the same variant)
    4. samtools index — create a .bai index so GATK can seek into the BAM

    Read-group tag (@RG)
    --------------------
    GATK requires every read in the BAM to have an @RG (read-group) tag
    that carries at minimum an ID, a sample name (SM), and a platform (PL).
    The -R flag to bwa mem embeds this information in the BAM header.
    """
    input:
        r1  ="results/blended/{scenario}/{replicate}/{scenario}_R1.fastq.gz",
        r2  ="results/blended/{scenario}/{replicate}/{scenario}_R2.fastq.gz",
        ref =rules.bwa_index.output.ref,
        bwt =rules.bwa_index.output.bwt,   # explicit dependency on index files
    output:
        bam ="results/variant_calling/{scenario}/{replicate}/aligned.bam",
        bai ="results/variant_calling/{scenario}/{replicate}/aligned.bam.bai",
        nc = temp("results/variant_calling/{scenario}/{replicate}/namecollate.bam.tmp"),
        fm = temp("results/variant_calling/{scenario}/{replicate}/fixmate.bam.tmp"),
        sorted = temp("results/variant_calling/{scenario}/{replicate}/sorted.bam.tmp"),
    params:
        # Read-group string embedded in the BAM header.
        # ID  : unique run identifier (scenario + replicate)
        # SM  : sample name shown in VCF genotype columns
        # PL  : sequencing platform (must be one of GATK's recognised values)
        # LB  : library name (used by markdup to identify duplicate pairs)
        rg=lambda wc: (
            f"@RG\\tID:{wc.scenario}_{wc.replicate}"
            f"\\tSM:{wc.scenario}"
            f"\\tPL:ILLUMINA"
            f"\\tLB:{wc.scenario}_{wc.replicate}"
        ),
        extra=config["variant_calling"].get("bwa_extra_flags", ""),
    threads:
        config["variant_calling"]["threads"]
    resources:
        # Reserve enough RAM for samtools sort (roughly 768 MB per thread)
        mem_mb=lambda wc, threads: threads * 768,
    log:
        "logs/variant_calling/{scenario}/{replicate}.align.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        # Pipe: bwa mem → samtools sort → write sorted BAM
        # -R adds the read-group tag required by GATK
        # -t sets the number of threads
        # 2>> {log} appends stderr (progress messages) to the log file
        bwa mem \
            -R '{params.rg}' \
            -t {threads} \
            -o {output.bam} \
            {params.extra} \
            {input.ref} \
            {input.r1} {input.r2} \
            2>> {log} \
            2>&1
        samtools collate -o {output.nc} {output.bam} >> {log} 2>&1
        samtools fixmate -m {output.nc} {output.fm} >> {log} 2>&1
        samtools sort -@ {threads} -o {output.sorted} {output.fm} >> {log} 2>&1
        # Mark PCR duplicate read pairs in the sorted BAM.
        # samtools markdup flags duplicates without removing them;
        # GATK will then skip flagged reads automatically.
        samtools markdup \
            -@ {threads} \
            {output.sorted} \
            {output.bam}.tmp \
            >> {log} 2>&1
        mv {output.bam}.tmp {output.bam}

        # Index the BAM so that GATK can jump to any genomic position quickly.
        samtools index {output.bam} >> {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule 3 — call variants with GATK HaplotypeCaller
# ---------------------------------------------------------------------------

rule gatk_haplotype_caller:
    """
    Call single-nucleotide variants (and small indels) with GATK
    HaplotypeCaller.

    What HaplotypeCaller does
    -------------------------
    For each genomic region with sufficient read coverage it:
      1. Locally re-assembles the reads into candidate haplotypes.
      2. Evaluates each haplotype against the reads using a probabilistic
         model.
      3. Reports the most likely genotype at every site where a variant is
         supported.

    Key flags used
    --------------
    --emit-ref-confidence NONE
        Output only variant sites (default behaviour); change to BP_RESOLUTION
        or GVCF to get per-base coverage information.
    --min-base-quality-score
        Ignore bases with a quality below this threshold when assembling
        haplotypes (read from config, default 20).
    --sample-ploidy 2
        Assume a diploid organism.  For haploid bacteria, set this to 1 in
        the config; for polyploids, increase accordingly.

    Output
    ------
    A bgzip-compressed, tabix-indexed VCF file.  These formats are required
    by many downstream tools and are more space-efficient than plain VCF.
    """
    input:
        bam =rules.bwa_align.output.bam,
        bai =rules.bwa_align.output.bai,
        ref =rules.bwa_index.output.ref,
        fai =rules.bwa_index.output.fai,
        dic =rules.bwa_index.output.dic,
    output:
           # results/variant_calling/scenario_equal_mix/rep1/output.vcf
        vcf="results/variant_calling/{scenario}/{replicate}/output.vcf"#.gz",
#        tbi="results/variant_calling/{scenario}/{replicate}/output.vcf.gz.tbi",
    params:
        min_base_quality=config["variant_calling"].get("min_base_quality", 20),
        ploidy          =config["variant_calling"].get("ploidy", 2),
        extra           =config["variant_calling"].get("gatk_extra_flags", ""),
    threads:
        config["variant_calling"]["threads"]
    resources:
        # GATK's default Java heap is 4 GB; scale with thread count
        mem_mb=lambda wc, threads: max(8000, threads * 1500),
    log:
        "logs/variant_calling/{scenario}/{replicate}.gatk.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        gatk HaplotypeCaller \
            --reference              {input.ref} \
            --input                  {input.bam} \
            --output                 {output.vcf} \
            --min-base-quality-score {params.min_base_quality} \
            --sample-ploidy          {params.ploidy} \
            --native-pair-hmm-threads {threads} \
            {params.extra} \
            >> {log} 2>&1

        # GATK writes the .tbi index automatically when the output filename
        # ends in .vcf.gz, but we declare it explicitly so Snakemake knows
        # it was produced and can use it as an input to downstream rules.
        """
