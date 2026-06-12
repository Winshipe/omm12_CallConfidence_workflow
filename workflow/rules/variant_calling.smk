"""
rules/variant_calling.smk
=========================
Align blended reads to unmutated reference genome(s) and call variants
using BWA-MEM (alignment) and GATK HaplotypeCaller (variant calling).

Pipeline overview for each scenario
------------------------------------
1. merge_references      – if a scenario uses more than one reference genome,
                           concatenate their FASTA files into a single combined
                           reference so BWA and GATK only need to handle one file.
2. bwa_index             – build the BWA and samtools index files for the
                           (possibly merged) reference.
3. bwa_align             – align the paired-end reads to the reference with
                           BWA-MEM.  Adds a read-group tag (required by GATK).
4. sort_and_index_bam    – sort the BAM by coordinate and build an index.
                           GATK requires a sorted, indexed BAM.
5. mark_duplicates       – flag PCR duplicates with GATK MarkDuplicates so
                           they do not inflate variant evidence.
6. gatk_haplotype_caller – call variants with GATK HaplotypeCaller in
                           GVCF mode (one intermediate file per scenario).
7. gatk_genotype_gvcfs   – convert the GVCF to a final VCF containing
                           only variant sites.
8. filter_variants       – apply hard quality filters to flag low-confidence
                           calls (FILTER column in the VCF).

Output that is consumed by assess.smk
--------------------------------------
  results/breseq/{scenario}/output/output.vcf

  The output path is intentionally kept identical to the breseq module so
  that assess.smk requires no changes when swapping callers.

Configuration keys used from config.yaml
-----------------------------------------
  config["references"]                      — dict {ref_id: fasta_path}
  config["scenarios"]                       — dict {scenario: [contributions]}
  config["variant_calling"]["threads"]      — CPU threads for BWA / GATK
  config["variant_calling"]["mem_mb"]       — RAM for GATK (megabytes)
  config["variant_calling"]["min_base_quality_score"]  — -mbq flag (default 20)
  config["variant_calling"]["gatk_extra_flags"]        — any extra GATK flags

  See the updated config.yaml snippet at the bottom of this file for the
  recommended variant_calling block.
"""


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def scenario_ref_ids(wildcards):
    """
    Return the list of unique reference IDs that contribute reads in
    this scenario (preserving the order they appear in the config).
    """
    seen = []
    for c in config["scenarios"][wildcards.scenario]:
        if c["ref_id"] not in seen:
            seen.append(c["ref_id"])
    return seen


def scenario_fasta_paths(wildcards):
    """
    Return the list of unmutated FASTA paths for all references in
    this scenario, in the same order as scenario_ref_ids().
    """
    return [config["references"][r] for r in scenario_ref_ids(wildcards)]


def merged_ref_path(wildcards):
    """
    Return the path to the (possibly merged) reference FASTA for a scenario.
    If only one reference is used the original FASTA is symlinked; if multiple
    references are used they are concatenated into a combined file.
    """
    return f"results/gatk/{{scenario}}/reference/combined_reference.fasta"


# ---------------------------------------------------------------------------
# Rule 1 – merge reference FASTAs for this scenario
# ---------------------------------------------------------------------------

rule merge_references:
    """
    Concatenate all reference FASTAs that contribute to this scenario into
    a single FASTA file.

    Why?  BWA and GATK each work against a single reference file.  When a
    scenario mixes reads from more than one genome, both genomes must be
    present in the reference so that reads align to the correct sequence.

    If only one reference is used, this step simply copies the file.
    The sequence IDs inside each FASTA are preserved unchanged, so
    downstream tools can still tell which genome a variant comes from.
    """
    input:
        fastas=scenario_fasta_paths,
    output:
        combined="results/gatk/{scenario}/reference/combined_reference.fasta",
    log:
        "logs/variant_calling/{scenario}.merge_refs.log",
    shell:
        # cat works for both single and multiple input files
        "cat {input.fastas} > {output.combined} 2> {log}"


# ---------------------------------------------------------------------------
# Rule 2 – index the combined reference (BWA + samtools dict/fai)
# ---------------------------------------------------------------------------

rule bwa_index:
    """
    Build the index files that BWA and GATK need before they can use a
    reference genome.

    BWA index  (.amb / .ann / .bwt / .pac / .sa)
      — used by BWA-MEM during alignment.

    samtools faidx (.fai)
      — a plain-text index that lets tools quickly look up any region of
        the genome without reading the whole file.

    GATK CreateSequenceDictionary (.dict)
      — a SAM-format header listing every contig name and length; GATK
        requires this to validate that the reference matches the BAM.
    """
    input:
        ref="results/gatk/{scenario}/reference/combined_reference.fasta",
    output:
        # BWA index files
        amb="results/gatk/{scenario}/reference/combined_reference.fasta.amb",
        ann="results/gatk/{scenario}/reference/combined_reference.fasta.ann",
        bwt="results/gatk/{scenario}/reference/combined_reference.fasta.bwt",
        pac="results/gatk/{scenario}/reference/combined_reference.fasta.pac",
        sa="results/gatk/{scenario}/reference/combined_reference.fasta.sa",
        # samtools index
        fai="results/gatk/{scenario}/reference/combined_reference.fasta.fai",
        # GATK sequence dictionary
        dict="results/gatk/{scenario}/reference/combined_reference.dict",
    log:
        "logs/variant_calling/{scenario}.bwa_index.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        # Build BWA index (creates .amb .ann .bwt .pac .sa alongside the FASTA)
        bwa index {input.ref} &>> {log}

        # Build samtools FASTA index (.fai)
        samtools faidx {input.ref} &>> {log}

        # Build GATK sequence dictionary (.dict)
        gatk CreateSequenceDictionary \
            --REFERENCE {input.ref} \
            --OUTPUT {output.dict} \
            &>> {log}
        """


# ---------------------------------------------------------------------------
# Rule 3 – align reads to the combined reference with BWA-MEM
# ---------------------------------------------------------------------------

rule bwa_align:
    """
    Align paired-end reads from a scenario to the combined reference genome
    using BWA-MEM, the standard short-read aligner for Illumina data.

    Key options used:
      -t   — number of CPU threads (set by config)
      -R   — read-group tag; GATK *requires* this field to be present.
             The tag records the sample name (SM:), library (LB:), and
             sequencing platform (PL:ILLUMINA).

    The alignment is piped directly to samtools to convert to BAM format,
    avoiding a large intermediate SAM file on disk.
    """
    input:
        r1="results/blended/{scenario}/{scenario}_R1.fastq.gz",
        r2="results/blended/{scenario}/{scenario}_R2.fastq.gz",
        ref="results/gatk/{scenario}/reference/combined_reference.fasta",
        # Explicit dependency on the index so Snakemake schedules correctly
        bwt="results/gatk/{scenario}/reference/combined_reference.fasta.bwt",
    output:
        bam=temp("results/gatk/{scenario}/aligned.bam"),
    params:
        # Read-group string.  SM (sample name) matches the scenario name so
        # GATK can track which sample produced each call.
        rg=r"@RG\tID:{scenario}\tSM:{scenario}\tLB:{scenario}\tPL:ILLUMINA",
    threads:
        config["variant_calling"]["threads"]
    log:
        "logs/variant_calling/{scenario}.bwa_align.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        bwa mem \
            -t {threads} \
            -R '{params.rg}' \
            {input.ref} \
            {input.r1} {input.r2} \
        | samtools view -b -o {output.bam} \
        &> {log}
        """


# ---------------------------------------------------------------------------
# Rule 4 – sort BAM by coordinate and build index
# ---------------------------------------------------------------------------

rule sort_and_index_bam:
    """
    Sort the aligned reads by their position in the reference genome and
    build a BAM index.

    Why sort?  GATK requires reads to be in coordinate order so it can
    process one genomic region at a time without loading the whole file.

    Why index?  The index (.bai file) lets tools jump directly to any
    position in the BAM without reading from the start — essential for
    performance on large files.
    """
    input:
        bam="results/gatk/{scenario}/aligned.bam",
    output:
        sorted_bam=temp("results/gatk/{scenario}/aligned.sorted.bam"),
        bai=temp("results/gatk/{scenario}/aligned.sorted.bam.bai"),
    threads:
        config["variant_calling"]["threads"]
    log:
        "logs/variant_calling/{scenario}.sort_bam.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        samtools sort \
            -@ {threads} \
            -o {output.sorted_bam} \
            {input.bam} \
            &> {log}

        samtools index {output.sorted_bam} &>> {log}
        """


# ---------------------------------------------------------------------------
# Rule 5 – mark PCR duplicates
# ---------------------------------------------------------------------------

rule mark_duplicates:
    """
    Identify and flag reads that are PCR duplicates using GATK
    MarkDuplicates.

    PCR duplicates are multiple reads that originate from the *same* DNA
    fragment amplified during library preparation.  They are not independent
    observations of the genome and would artificially inflate the apparent
    support for a variant.

    MarkDuplicates does NOT remove the reads; it sets a flag in the BAM so
    that GATK HaplotypeCaller knows to down-weight them.  A metrics file
    recording how many duplicates were found is written alongside the BAM.
    """
    input:
        bam="results/gatk/{scenario}/aligned.sorted.bam",
        bai="results/gatk/{scenario}/aligned.sorted.bam.bai",
    output:
        bam="results/gatk/{scenario}/deduped.bam",
        bai="results/gatk/{scenario}/deduped.bam.bai",
        metrics="results/gatk/{scenario}/duplicate_metrics.txt",
    resources:
        mem_mb=config["variant_calling"].get("mem_mb", 16000),
    log:
        "logs/variant_calling/{scenario}.markdup.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        gatk MarkDuplicates \
            --INPUT  {input.bam} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --CREATE_INDEX true \
            &> {log}
        """


# ---------------------------------------------------------------------------
# Rule 6 – call variants with GATK HaplotypeCaller (GVCF mode)
# ---------------------------------------------------------------------------

rule gatk_haplotype_caller:
    """
    Run GATK HaplotypeCaller to identify candidate variant sites.

    This rule produces a GVCF (Genomic VCF) rather than a standard VCF.
    A GVCF contains a record for *every* position in the genome — both
    variant and non-variant sites — which allows GATK to distinguish
    "no evidence for a variant" from "no data at all".

    Important flags used here:
      --ploidy 1
          Bacterial genomes are haploid; this tells GATK to expect only
          one copy of each chromosome.  Using the default (diploid) would
          produce incorrect genotype calls.

      --emit-ref-confidence GVCF
          Activates GVCF mode.

      --min-base-quality-score
          Ignore bases with a Phred quality score below this threshold
          (default 20, i.e. ≥99 % base-call accuracy).

    The GVCF is an intermediate file that is converted to a final VCF in
    the next rule.
    """
    input:
        bam="results/gatk/{scenario}/deduped.bam",
        bai="results/gatk/{scenario}/deduped.bam.bai",
        ref="results/gatk/{scenario}/reference/combined_reference.fasta",
        fai="results/gatk/{scenario}/reference/combined_reference.fasta.fai",
        dict="results/gatk/{scenario}/reference/combined_reference.dict",
    output:
        gvcf=temp("results/gatk/{scenario}/raw_calls.g.vcf.gz"),
    params:
        mbq=config["variant_calling"].get("min_base_quality_score", 20),
        extra=config["variant_calling"].get("gatk_extra_flags", ""),
    threads:
        config["variant_calling"]["threads"]
    resources:
        mem_mb=config["variant_calling"].get("mem_mb", 16000),
    log:
        "logs/variant_calling/{scenario}.haplotype_caller.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        gatk HaplotypeCaller \
            --reference {input.ref} \
            --input     {input.bam} \
            --output    {output.gvcf} \
            --sample-ploidy 1 \
            --emit-ref-confidence GVCF \
            --min-base-quality-score {params.mbq} \
            --native-pair-hmm-threads {threads} \
            {params.extra} \
            &> {log}
        """


# ---------------------------------------------------------------------------
# Rule 7 – genotype the GVCF to produce a standard VCF
# ---------------------------------------------------------------------------

rule gatk_genotype_gvcfs:
    """
    Convert the GVCF produced by HaplotypeCaller into a standard VCF that
    contains only variant sites.

    GenotypeGVCFs re-evaluates the raw likelihoods across all positions,
    applies genotype priors, and emits one line per variant.  The output
    is a compressed, indexed VCF (.vcf.gz).
    """
    input:
        gvcf="results/gatk/{scenario}/raw_calls.g.vcf.gz",
        ref="results/gatk/{scenario}/reference/combined_reference.fasta",
        fai="results/gatk/{scenario}/reference/combined_reference.fasta.fai",
        dict="results/gatk/{scenario}/reference/combined_reference.dict",
    output:
        vcf_gz=temp("results/gatk/{scenario}/genotyped.vcf.gz"),
    resources:
        mem_mb=config["variant_calling"].get("mem_mb", 16000),
    log:
        "logs/variant_calling/{scenario}.genotype_gvcfs.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        gatk GenotypeGVCFs \
            --reference {input.ref} \
            --variant   {input.gvcf} \
            --output    {output.vcf_gz} \
            --sample-ploidy 1 \
            &> {log}
        """


# ---------------------------------------------------------------------------
# Rule 8 – hard-filter variants and write the final VCF
# ---------------------------------------------------------------------------

rule filter_variants:
    """
    Apply hard quality filters to flag low-confidence variant calls.

    GATK's recommended approach for small datasets (where Variant Quality
    Score Recalibration, VQSR, cannot be used) is "hard filtering": apply
    simple threshold rules and mark calls that fail with a label in the
    FILTER column of the VCF.

    Filters applied (SNPs only):
      QD < 2.0   — QualByDepth: variant quality divided by read depth.
                   Low values suggest the site is noisy relative to
                   coverage.
      FS > 60.0  — FisherStrand: strand bias measured as a Phred score.
                   High values mean far more supporting reads come from
                   one strand, suggesting a sequencing artefact.
      MQ < 40.0  — RMSMappingQuality: average mapping quality of reads
                   covering the site.  Low values mean reads align poorly.
      MQRankSum < -12.5
                 — Difference in mapping quality between ref and alt reads.
      ReadPosRankSum < -8.0
                 — Whether the alt allele tends to appear only at the ends
                   of reads (a common artefact).

    Calls that *pass* all filters are marked PASS; calls that fail one or
    more are labelled with the filter name(s).  No calls are removed — the
    assess_variants script uses the QUAL score to decide what to trust.

    The output VCF is written to the path expected by assess.smk so that
    the rest of the pipeline works without modification.
    """
    input:
        vcf_gz="results/gatk/{scenario}/genotyped.vcf.gz",
        ref="results/gatk/{scenario}/reference/combined_reference.fasta",
        fai="results/gatk/{scenario}/reference/combined_reference.fasta.fai",
        dict="results/gatk/{scenario}/reference/combined_reference.dict",
    output:
        # *** This path matches what assess.smk expects ***
        vcf="results/breseq/{scenario}/output/output.vcf",
    log:
        "logs/variant_calling/{scenario}.filter_variants.log",
    conda:
        "../envs/bwa_gatk.yaml"
    shell:
        """
        gatk VariantFiltration \
            --reference {input.ref} \
            --variant   {input.vcf_gz} \
            --output    {output.vcf} \
            --filter-expression "QD < 2.0"            --filter-name "LowQD" \
            --filter-expression "FS > 60.0"           --filter-name "HighFS" \
            --filter-expression "MQ < 40.0"           --filter-name "LowMQ" \
            --filter-expression "MQRankSum < -12.5"   --filter-name "LowMQRankSum" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "LowReadPosRankSum" \
            &> {log}
        """


# ---------------------------------------------------------------------------
# Usage note
# ---------------------------------------------------------------------------
# To run only the variant-calling module for a specific scenario:
#
#   snakemake --use-conda \
#       results/breseq/scenario_equal_mix/output/output.vcf
#
# To run all scenarios defined in the config:
#
#   snakemake --use-conda results/breseq/{scenario}/output/output.vcf \
#       --wildcards scenario=$(python -c \
#           "import yaml; c=yaml.safe_load(open('config/config.yaml')); \
#            print(' '.join(c['scenarios']))")
#
# ---------------------------------------------------------------------------
#
# Recommended variant_calling block for config/config.yaml
# ---------------------------------------------------------
#
# variant_calling:
#   threads: 8          # CPU threads for BWA and GATK
#   mem_mb: 16000       # RAM for GATK steps (megabytes; 16 GB is usually enough)
#   min_base_quality_score: 20   # ignore bases with Phred quality below this
#   gatk_extra_flags: ""         # any additional GATK HaplotypeCaller flags
#
# ---------------------------------------------------------------------------
