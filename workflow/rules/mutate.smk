"""
rules/mutate.smk  — add mutations to reference sequences
"""


def get_ref_fasta(wildcards):
    return config["references"][wildcards.ref_id]


rule mutate_reference:
    """
    Apply a substitution model to a reference FASTA and produce:
      - a mutated FASTA
      - a TSV ground-truth table (position, ref_base, alt_base, mutation_type)
    """
    input:
        ref=get_ref_fasta,
    output:
        mutated_fasta="results/mutated/{ref_id}/{ref_id}.mutated.fasta",
        mutations_tsv="results/mutated/{ref_id}/{ref_id}.mutations.tsv",
    params:
        model=config["mutation"]["model"],
        rate=config["mutation"]["substitution_rate"],
        kappa=config["mutation"].get("kappa", 2.0),
        gc_freq=config["mutation"].get("gc_freq", 0.5),
        seed=sum([ord(c) for c in "results/mutated/{ref_id}/{ref_id}.mutated.fasta"]) #config["mutation"]["seed"],
    log:
        "logs/mutate/{ref_id}.log",
    conda:
        "../envs/mutate.yaml"
    script:
        "../scripts/mutate_reference.py"
