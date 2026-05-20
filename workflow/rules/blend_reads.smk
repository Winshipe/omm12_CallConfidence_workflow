"""
rules/blend_reads.smk  — combine simulated reads at scenario-defined abundances

Each scenario in config["scenarios"] is a list of contributions:
  ref_id            : reference ID (must match a key in config["references"])
  mutated_fraction  : fraction of reads for this ref drawn from the mutated pool
  abundance         : relative abundance (normalised to sum to 1)

The output is a gzip-compressed FASTQ pair for each scenario.
"""

import math


def scenario_input_reads(wildcards):
    """
    Return the flat list of R1/R2 FASTQ files needed for this scenario.
    Each contribution requires both the unmutated and mutated simulated reads
    (the blending script picks the right proportions at runtime).
    """
    contribs = config["scenarios"][wildcards.scenario]
    inputs = {}
    for i, c in enumerate(contribs):
        ref = c["ref_id"]
        for mutated in ("unmutated", "mutated"):
            for read in ("R1", "R2"):
                key = f"contrib_{i}_{mutated}_{read}"
                inputs[key] = (
                    f"results/simulated/{ref}/{mutated}/{ref}_{read}.fastq.gz"
                )
    return inputs


rule blend_reads:
    """
    Randomly subsample and interleave reads from each reference contribution
    to produce the target total read count for the scenario.
    """
    input:
        unpack(scenario_input_reads),
    output:
        r1="results/blended/{scenario}/{scenario}_R1.fastq.gz",
        r2="results/blended/{scenario}/{scenario}_R2.fastq.gz",
    params:
        scenario_cfg=lambda wc: config["scenarios"][wc.scenario],
        total_reads=config["blend_total_reads"],
        seed=config["mutation"]["seed"],
    log:
        "logs/blend/{scenario}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/blend_reads.py"
