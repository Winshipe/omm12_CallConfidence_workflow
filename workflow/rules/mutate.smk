"""
rules/mutate.smk — introduce random substitutions into reference sequences
===========================================================================

One independent set of mutations is produced for every (ref_id × replicate)
combination.  Each replicate uses a different random seed (derived from
config["mutation"]["base_seed"] + the replicate index) so that mutations are
statistically independent while still being fully reproducible.

Output files
------------
  results/mutated/{ref_id}/{replicate}/{ref_id}.mutated.fasta
      The reference sequence with substitutions applied.

  results/mutated/{ref_id}/{replicate}/{ref_id}.mutations.tsv
      Ground-truth table of every mutation introduced:
        seq_id | position | ref_base | alt_base | mutation_type

These files are consumed by simulate_reads.smk (mutated FASTA) and
assess.smk (mutations TSV).
"""


# ---------------------------------------------------------------------------
# Helper: look up the original reference FASTA for a given ref_id wildcard
# ---------------------------------------------------------------------------

def get_ref_fasta(wildcards):
    """Return the path to the unmutated reference FASTA."""
    return config["references"][wildcards.ref_id]


# ---------------------------------------------------------------------------
# Helper: compute a per-replicate seed from the base seed in config
# ---------------------------------------------------------------------------

def replicate_seed(wildcards):
    """
    Derive a deterministic integer seed for this replicate.

    Replicates are named rep1, rep2, …  The index is parsed from the name and
    added to the base seed so that every replicate is independent yet
    reproducible (re-running with the same config always gives the same result).
    """
    base   = int(config["mutation"]["base_seed"])
    # wildcards.replicate is e.g. "rep1"; strip the "rep" prefix to get the index
    index  = int(wildcards.replicate.replace("rep", ""))*10
    return base + index


# ---------------------------------------------------------------------------
# Rule
# ---------------------------------------------------------------------------

rule mutate_reference:
    """
    Apply a substitution model to a reference FASTA and produce a mutated
    version plus a ground-truth TSV recording every change.

    Supported substitution models (set via config["mutation"]["model"]):
      jukes_cantor  — equal rates for all substitution types
      tamura_nei    — allows different transition / transversion rates and
                      arbitrary nucleotide frequencies

    The random seed is replicate-specific so that each replicate yields a
    different (but reproducible) set of mutations.
    """
    input:
        ref=get_ref_fasta,
    output:
        mutated_fasta="results/mutated/{ref_id}/{replicate}/{ref_id}.mutated.fasta",
        mutations_tsv="results/mutated/{ref_id}/{replicate}/{ref_id}.mutations.tsv",
    params:
        model   = config["mutation"]["model"],
        rate    = config["mutation"]["substitution_rate"],
        kappa   = config["mutation"].get("kappa", 2.0),
        gc_freq = config["mutation"].get("gc_freq", 0.5),
        # Seed is evaluated at runtime via the lambda so it has access to wildcards
        seed    = replicate_seed,
    log:
        "logs/mutate/{ref_id}/{replicate}.log",
    conda:
        "../envs/mutate.yaml"
    script:
        "../scripts/mutate_reference.py"
