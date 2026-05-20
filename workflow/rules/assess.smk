"""
rules/assess.smk  — compare breseq calls against ground-truth mutations

For each scenario, gather the ground-truth TSV files for every mutated
reference that contributes reads, then run the assessment script.
"""


def scenario_ground_truth_tsvs(wildcards):
    """
    Return a dict {ref_id: mutations_tsv} for all refs in this scenario
    that have a non-zero mutated_fraction.
    """
    contribs = config["scenarios"][wildcards.scenario]
    out = {}
    for c in contribs:
        ref = c["ref_id"]
        if c["mutated_fraction"] > 0 and ref not in out:
            out[ref] = f"results/mutated/{ref}/{ref}.mutations.tsv"
    return list(out.values())


rule assess_variants:
    """
    Compare breseq Genome Diff output to ground-truth mutation tables.
    Produces a TSV with one row per expected mutation, recording whether
    it was detected, the called ALT, and the breseq quality score.
    """
    input:
        vcf="results/breseq/{scenario}/output/output.vcf",
        ground_truth=scenario_ground_truth_tsvs,
    output:
        tsv="results/assessment/{scenario}_assessment.tsv",
    params:
        min_quality=config["assessment"]["min_quality"],
    log:
        "logs/assessment/{scenario}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assess_variants.py"
