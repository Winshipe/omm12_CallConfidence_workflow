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

3. Variant calling, with the caller chosen by config["variant_calling"]["caller"]:

   "haplotypecaller" (default)
     gatk_haplotype_caller    — GATK HaplotypeCaller (GVCF mode) → haplotypecaller.vcf

   "mutect2"
     gatk_mutect2             — GATK Mutect2 in tumor-only mode → mutect2.unfiltered.vcf
                                (+ the .stats file FilterMutectCalls needs)
     gatk_filter_mutect_calls — GATK FilterMutectCalls → mutect2.filtered.vcf;
                                failing calls keep a non-PASS FILTER value and
                                are ignored by assess_variants.py

   Mutect2 models each variant's allele fraction directly instead of a fixed
   ploidy, so it suits the mixed (mutated + unmutated) samples simulated here,
   where mutations sit at the scenario's mutated_fraction.

4. select_variant_calls  — copy the selected caller's VCF to output.vcf, the
                           file assessment reads.  Each caller keeps its own
                           file, so switching caller does not overwrite the
                           other caller's results.

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
  results/variant_calling/{scenario}/{replicate}/haplotypecaller.vcf     HaplotypeCaller calls
  results/variant_calling/{scenario}/{replicate}/mutect2.filtered.vcf    Mutect2 calls with FILTER set
  results/variant_calling/{scenario}/{replicate}/output.vcf              selected caller's calls
  results/variant_calling/{scenario}/{replicate}/aligned.bam          sorted, deduplicated BAM
  results/variant_calling/{scenario}/{replicate}/aligned.bam.bai      BAM index

Configuration keys read from config["variant_calling"]
-------------------------------------------------------
  threads          : int   — CPU threads for BWA and GATK (default 8)
  bwa_extra_flags  : str   — optional extra flags passed to bwa mem
  gatk_extra_flags : str   — optional extra flags passed to HaplotypeCaller
  min_base_quality : int   — minimum base quality score for GATK (default 20)
  caller           : str   — "haplotypecaller" (default) or "mutect2"
  mutect2_extra_flags       : str — optional extra flags passed to Mutect2
  filter_mutect_extra_flags : str — optional extra flags passed to FilterMutectCalls
"""

import os


# ---------------------------------------------------------------------------
# Variant caller selection
# ---------------------------------------------------------------------------

# VCF written by each supported caller (relative to the scenario/replicate dir)
CALLER_VCF = {
    "haplotypecaller": "haplotypecaller.vcf",
    "mutect2":         "mutect2.filtered.vcf",
}
VARIANT_CALLER = str(config["variant_calling"].get("caller", "haplotypecaller")).lower()
if VARIANT_CALLER not in CALLER_VCF:
    raise ValueError(
        f"config variant_calling.caller must be one of {sorted(CALLER_VCF)}, "
        f"got {VARIANT_CALLER!r}"
    )


def selected_caller_vcf(wildcards):
    """VCF produced by the configured caller for this scenario × replicate."""
    return (f"results/variant_calling/{wildcards.scenario}/{wildcards.replicate}/"
            f"{CALLER_VCF[VARIANT_CALLER]}")


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
        sam ="results/variant_calling/{scenario}/{replicate}/aligned.sam",
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
        samtools view {output.bam} > {output.sam}
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
        # Copied to output.vcf by select_variant_calls when this caller is selected
        vcf="results/variant_calling/{scenario}/{replicate}/haplotypecaller.vcf",
    params:
        min_base_quality=config["variant_calling"].get("min_base_quality", 20),
        ploidy          =config["variant_calling"].get("ploidy", 2), 
        extra           =config["variant_calling"].get("gatk_extra_flags", ""),
    threads:
        config["variant_calling"]["threads"]
    resources:
        # GATK's default Java heap is 4 GB; scale with thread count
        mem_mb= config["variant_calling"].get("mem_mb",lambda wc, threads: max(8000, threads * 1500)),
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
            -ERC GVCF \
            {params.extra} \
            >> {log} 2>&1

        # GATK writes the .tbi index automatically when the output filename
        # ends in .vcf.gz, but we declare it explicitly so Snakemake knows
        # it was produced and can use it as an input to downstream rules.
        """


# ---------------------------------------------------------------------------
# Rule 3 (alternative) — call variants with GATK Mutect2 + FilterMutectCalls
# ---------------------------------------------------------------------------

rule gatk_mutect2:
    """
    Call somatic-style variants with GATK Mutect2 in tumor-only mode (no
    matched normal): every read in the BAM belongs to the one sample.

    Unlike HaplotypeCaller, Mutect2 does not assume a ploidy; it estimates
    each variant's allele fraction, so low-fraction mutations in a mixed
    sample can be called.  Mutect2 leaves QUAL as "." — call confidence is
    expressed through FilterMutectCalls' FILTER column instead.

    Mutect2 also writes <output>.stats, which FilterMutectCalls requires.
    """
    input:
        bam =rules.bwa_align.output.bam,
        bai =rules.bwa_align.output.bai,
        ref =rules.bwa_index.output.ref,
        fai =rules.bwa_index.output.fai,
        dic =rules.bwa_index.output.dic,
    output:
        vcf  ="results/variant_calling/{scenario}/{replicate}/mutect2.unfiltered.vcf",
        stats="results/variant_calling/{scenario}/{replicate}/mutect2.unfiltered.vcf.stats",
    params:
        min_base_quality=config["variant_calling"].get("min_base_quality", 20),
        extra           =config["variant_calling"].get("mutect2_extra_flags", ""),
    threads:
        config["variant_calling"]["threads"]
    resources:
        mem_mb=config["variant_calling"].get("mem_mb", 8000),
    log:
        "logs/variant_calling/{scenario}/{replicate}.mutect2.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        gatk Mutect2 \
            --reference              {input.ref} \
            --input                  {input.bam} \
            --output                 {output.vcf} \
            --min-base-quality-score {params.min_base_quality} \
            --native-pair-hmm-threads {threads} \
            {params.extra} \
            >> {log} 2>&1
        """


rule gatk_filter_mutect_calls:
    """
    Apply GATK FilterMutectCalls to the Mutect2 calls.  Every call is kept in
    the output; calls that fail a filter get the failing filter names in the
    FILTER column instead of PASS, and assess_variants.py ignores them.
    """
    input:
        vcf  =rules.gatk_mutect2.output.vcf,
        stats=rules.gatk_mutect2.output.stats,
        ref  =rules.bwa_index.output.ref,
        fai  =rules.bwa_index.output.fai,
        dic  =rules.bwa_index.output.dic,
    output:
        vcf="results/variant_calling/{scenario}/{replicate}/mutect2.filtered.vcf",
    params:
        extra=config["variant_calling"].get("filter_mutect_extra_flags", ""),
    resources:
        mem_mb=config["variant_calling"].get("mem_mb", 8000),
    log:
        "logs/variant_calling/{scenario}/{replicate}.filter_mutect.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        gatk FilterMutectCalls \
            --reference {input.ref} \
            --variant   {input.vcf} \
            --stats     {input.stats} \
            --output    {output.vcf} \
            {params.extra} \
            >> {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Rule 4 — hand the selected caller's VCF to assessment
# ---------------------------------------------------------------------------

rule select_variant_calls:
    """
    Copy the VCF of the caller chosen in config["variant_calling"]["caller"]
    to output.vcf, the file assess_variants reads.
    """
    input:
        vcf=selected_caller_vcf,
    output:
        vcf="results/variant_calling/{scenario}/{replicate}/output.vcf",
    shell:
        "cp {input.vcf} {output.vcf}"
