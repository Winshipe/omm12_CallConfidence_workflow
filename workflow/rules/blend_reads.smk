"""
rules/blend_reads.smk — combine simulated reads at scenario-defined abundances
===============================================================================

For each (scenario × replicate) combination, reads simulated from multiple
references are mixed at user-specified relative abundances and mutated
fractions to produce a single paired-end FASTQ pair.

Scenario configuration
-----------------------
Each scenario in config["scenarios"] is a list of contributions, e.g.:

    scenario_equal_mix:
      - ref_id: ref_A
        mutated_fraction: 0.5   # half the ref_A reads come from the mutated version
        abundance: 1.0
      - ref_id: ref_B
        mutated_fraction: 0.5
        abundance: 1.0

The same scenario definition is reused for every replicate; only the
underlying simulated reads (and their mutations) differ between replicates.

The blending script (blend_reads.py) samples exactly
config["blend_total_reads"] read pairs per output file.

Output
------
  results/blended/{scenario}/{replicate}/{scenario}_R1.fastq.gz
  results/blended/{scenario}/{replicate}/{scenario}_R2.fastq.gz
"""


def scenario_input_reads(wildcards):
    """
    Build the dictionary of FASTQ input paths that blend_reads.py needs.

    For each contribution i in the scenario, we need both the 'unmutated' and
    'mutated' simulated reads (the blending script picks the right proportions
    at runtime based on mutated_fraction).  Keys follow the naming convention
    expected by blend_reads.py:  contrib_{i}_{unmutated|mutated}_{R1|R2}
    """
    contribs = config["scenarios"][wildcards.scenario]
    inputs   = {}
    for i, c in enumerate(contribs):
        ref = c["ref_id"]
        for mutated in ("unmutated", "mutated"):
            for read in ("R1", "R2"):
                key = f"contrib_{i}_{mutated}_{read}"
                inputs[key] = (
                    f"results/simulated/{ref}/{wildcards.replicate}"
                    f"/{mutated}/{ref}_{read}.fastq.gz"
                )
    return inputs


def replicate_seed(wildcards):
    """
    Derive a per-replicate blending seed from the base seed.
    This ensures the random subsampling in blend_reads.py is also independent
    across replicates.
    """
    base  = int(config["mutation"]["base_seed"])
    index = int(wildcards.replicate.replace("rep", ""))
    return base + index


#rule blend_reads:
#    """
#    Randomly subsample and interleave reads from each reference contribution
#    to produce a FASTQ pair for a given scenario and replicate.
#
#    The total number of read pairs written is config["blend_total_reads"].
#    Reads are shuffled after blending so they are not ordered by source
#    reference or mutation status.
#    """
#    input:
#        unpack(scenario_input_reads),
#    output:
#        r1="results/blended/{scenario}/{replicate}/{scenario}_R1.fastq.gz",
#        r2="results/blended/{scenario}/{replicate}/{scenario}_R2.fastq.gz",
#    params:
#        scenario_cfg = lambda wc: config["scenarios"][wc.scenario],
#        #total_reads  = config["blend_total_reads"],
#        seed         = replicate_seed,
#    log:
#        "logs/blend/{scenario}/{replicate}.log",
#    conda:
#        "../envs/python.yaml"
#    script:
#        "../scripts/blend_reads.py"

def prepare_seqtk_sampling_script(wildcards):
    rel_abs_dict = dict()
    maf_abs_dict = dict()
    for genome_cfg in config["scenarios"][wildcards.scenario]:
        rel_abs_dict[genome_cfg["ref_id"]] = float(genome_cfg["abundance"])
        maf_abs_dict[genome_cfg["ref_id"]] = float(genome_cfg["mutated_fraction"])
    sum_abundances = sum(rel_abs_dict.values())
    output = ""
    template = "seqtk sample -s 42 results/simulated/{{ref_id}}/{{replicate}}/{{mut_or_ref}}/{{ref_id}}_R{{p}}.fastq.gz {{fract}} >> results/blended/{{scenario}}/{{replicate}}/{{scenario}}_R{{p}}.fastq"
    for ref_id, ra in rel_abs_dict.items():
        mut_fraction = ra * maf_abs_dict[ref_id] #relative abundance * mutant fraction
        ref_fraction = ra * (1 - maf_abs_dict[ref_id])
        for p in [1,2]:
            for mr in ["mutated","unmutated"]:
                output += temp.format(ref_id=ref_id, replicate=wildcards.replicate, mut_or_ref=mr, p=p, fract=mut_fraction if mr == "mutated" else ref_fraction, scenario=wildcards.scenario)
    output += "gzip results/blended/{{scenario}}/{{replicate}}/{{scenario}}_R1.fastq".format(scenario=wildcards.scenario, replicate=wildcards.replicate)
    output += "gzip results/blended/{{scenario}}/{{replicate}}/{{scenario}}_R2.fastq".format(scenario=wildcards.scenario, replicate=wildcards.replicate)

    return output
rule blend_reads:
    input:
        unpack(scenario_input_reads)
    output:
        r1="results/blended/{scenario}/{replicate}/{scenario}_R1.fastq.gz",
        r2="results/blended/{scenario}/{replicate}/{scenario}_R2.fastq.gz",
    conda:
        "../envs/seqtk.yaml"
    run:
        shell(prepare_seqtk_sampling_script(wildcards))
