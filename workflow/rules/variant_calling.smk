"""
rules/variant_calling.smk — align blended reads to the reference and call variants
===================================================================================

breseq is run independently for every (scenario × replicate) combination.
Each run aligns the blended reads produced by blend_reads.smk against the
unmutated reference genome(s) and calls single-nucleotide polymorphisms in
polymorphism-prediction mode.

Why the unmutated reference?
-----------------------------
The whole point of the benchmark is to see whether the variant caller can
*discover* mutations relative to the original reference.  Using the unmutated
reference ensures that any called variants reflect true introduced mutations
(or false positives from the caller), not pre-existing differences.

Multiple references
--------------------
When a scenario blends reads from more than one reference genome, all
unmutated reference FASTAs are passed to breseq with separate -r flags.
breseq handles multi-reference inputs natively.

Output
------
  results/breseq/{scenario}/{replicate}/output/output.vcf
  results/breseq/{scenario}/{replicate}/output/summary.html
"""


def scenario_references(wildcards):
    """
    Return a list of unmutated FASTA paths for every reference that
    contributes reads in this scenario.  Duplicates are removed while
    preserving the order in which references first appear.
    """
    contribs = config["scenarios"][wildcards.scenario]
    seen     = {}
    for c in contribs:
        ref = c["ref_id"]
        if ref not in seen:
            seen[ref] = config["references"][ref]
    return list(seen.values())


rule run_breseq:
    """
    Run breseq in polymorphism mode against the unmutated reference(s).

    The output directory structure follows breseq conventions:
      output/output.vcf  — called variants in VCF format
      output/summary.html — human-readable HTML summary

    One independent breseq run is executed per (scenario × replicate) so that
    replicates are not mixed.
    """
    input:
        r1   = "results/blended/{scenario}/{replicate}/{scenario}_R1.fastq.gz",
        r2   = "results/blended/{scenario}/{replicate}/{scenario}_R2.fastq.gz",
        refs = scenario_references,
    output:
        vcf     = "results/breseq/{scenario}/{replicate}/output/output.vcf",
        summary = "results/breseq/{scenario}/{replicate}/output/summary.html",
    params:
        outdir    = "results/breseq/{scenario}/{replicate}",
        extra     = config["variant_calling"]["breseq_extra_flags"],
        # Build a '-r <ref>' flag for each reference
        ref_flags = lambda wc, input: " ".join(f"-r {r}" for r in input.refs),
    threads:
        config["variant_calling"]["threads"]
    log:
        "logs/variant_calling/{scenario}/{replicate}.log",
    resources:
        mem_mb=32000,
    conda:
        "../envs/breseq.yaml"
    shell:
        """
        breseq \
            {params.ref_flags} \
            --num-processors {threads} \
            --polymorphism-prediction \
            --brief-html-output \
            --output {params.outdir} \
            {params.extra} \
            {input.r1} {input.r2} \
            &> {log}
        """
