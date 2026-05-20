"""
rules/variant_calling.smk  — align blended reads to unmutated references and
                              call variants with breseq
"""


def scenario_references(wildcards):
    """
    Collect the unmutated FASTA files for all references participating
    in this scenario.  breseq accepts multiple -r references.
    """
    contribs = config["scenarios"][wildcards.scenario]
    seen = {}
    for c in contribs:
        ref = c["ref_id"]
        if ref not in seen:
            seen[ref] = config["references"][ref]
    return list(seen.values())


rule run_breseq:
    """
    Run breseq in polymorphism mode against the unmutated reference(s).
    The output directory contains output/output.gd (Genome Diff format)
    with all called variants.
    """
    input:
        r1="results/blended/{scenario}/{scenario}_R1.fastq.gz",
        r2="results/blended/{scenario}/{scenario}_R2.fastq.gz",
        refs=scenario_references,
    output:
        gd="results/breseq/{scenario}/output/output.vcf",
        summary="results/breseq/{scenario}/output/summary.html",
    params:
        outdir="results/breseq/{scenario}",
        extra=config["variant_calling"]["breseq_extra_flags"],
        ref_flags=lambda wc, input: " ".join(f"-r {r}" for r in input.refs),
    threads:
        config["variant_calling"]["threads"]
    log:
        "logs/variant_calling/{scenario}.log",
    resources:
        mem_mb = 32000,
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
