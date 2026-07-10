"""
rules/simulate_reads.smk — simulate MiSeq paired-end reads with ART
====================================================================

Reads are simulated independently for every (ref_id × replicate) combination.
Each replicate uses a different ART random seed (base_seed + replicate_index)
so that sequencing errors and read positions vary between replicates even when
the underlying FASTA is the same.

Two versions of reads are produced for each reference × replicate:
  • unmutated — reads from the original reference FASTA
  • mutated   — reads from the per-replicate mutated FASTA

Both are needed by blend_reads.smk, which mixes them at the ratios defined in
config["scenarios"].

Optional empirical error profile
---------------------------------
If config["simulation"]["empirical_reads_R1"] and empirical_reads_R2 are set to
real FASTQ files, ART's profiler is used to learn a custom quality-score
distribution before simulation.  Leave both as null to use ART's built-in
MiSeq (MSv3) profile.
"""


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_ref_fasta(wildcards):
    """Return the unmutated reference FASTA for a given ref_id."""
    return config["references"][wildcards.ref_id]


def art_ref_fasta(wildcards):
    """
    Return the correct FASTA depending on whether we are simulating
    'mutated' or 'unmutated' reads.
    """
    if wildcards.mutated == "mutated":
        # Per-replicate mutated FASTA produced by mutate.smk
        return (
            f"results/mutated/{wildcards.ref_id}/"
            f"{wildcards.replicate}/{wildcards.ref_id}.mutated.fasta"
        )
    # Original reference
    return config["references"][wildcards.ref_id]


def replicate_seed(wildcards):
    """
    Derive a per-replicate ART seed from the base seed.
    This makes sequencing errors independent across replicates.
    """
    base  = int(config["mutation"]["base_seed"])
    index = int(wildcards.replicate.replace("rep", ""))
    return base + index


# ---------------------------------------------------------------------------
# Optional rule: learn an empirical error profile from real reads
# (only executed when empirical_reads_R1/R2 are set in config)
# ---------------------------------------------------------------------------

rule art_learn_profile:
    """
    Run ART's built-in profiler on real MiSeq reads to derive a custom
    quality-score model.  This rule is only triggered when empirical read
    paths are supplied in the config.
    """
    input:
        r1=config["simulation"]["empirical_reads_R1"] or [],
        r2=config["simulation"]["empirical_reads_R2"] or [],
    output:
        profile_r1="results/art_profile/empirical_R1.txt",
        profile_r2="results/art_profile/empirical_R2.txt",
    log:
        "logs/simulation/art_learn_profile.log",
    conda:
        "../envs/art.yaml"
    shell:
        """
        art_profiler_illumina {output.profile_r1} {input.r1} 150 &>> {log}
        art_profiler_illumina {output.profile_r2} {input.r2} 150 &>> {log}
        """


def art_profile_flags(wildcards):
    """
    Return the ART command-line flags that select the quality-score profile:
      - empirical profile flags if real reads were provided
      - built-in platform profile flag otherwise
    """
    if config["simulation"]["empirical_reads_R1"] is not None:
        return (
            "--qprof1 results/art_profile/empirical_R1.txt "
            "--qprof2 results/art_profile/empirical_R2.txt"
        )
    return f"-ss {config['simulation']['platform']}"


# ---------------------------------------------------------------------------
# Rule: simulate reads for one reference × {unmutated|mutated} × replicate
# ---------------------------------------------------------------------------

rule simulate_reads:
    """
    Run ART to produce paired-end FASTQ reads from a single reference FASTA.

    The {mutated} wildcard is either 'mutated' or 'unmutated'; the {replicate}
    wildcard is e.g. 'rep1', 'rep2', etc.  A unique random seed is derived per
    replicate so that each replicate's reads are independent.

    Output files are gzip-compressed to save disk space.
    """
    input:
        ref=art_ref_fasta,
        # If an empirical profile is configured, declare the profile files as
        # inputs so Snakemake runs art_learn_profile first.
        profile=lambda wc: (
            [
                "results/art_profile/empirical_R1.txt",
                "results/art_profile/empirical_R2.txt",
            ]
            if config["simulation"]["empirical_reads_R1"] is not None
            else []
        ),
    output:
        r1="results/simulated/{ref_id}/{replicate}/{mutated}/{ref_id}_R1.fastq.gz",
        r2="results/simulated/{ref_id}/{replicate}/{mutated}/{ref_id}_R2.fastq.gz",
        sam="results/simulated/{ref_id}/{replicate}/{mutated}/{ref_id}_.sam"
    params:
        profile_flags = art_profile_flags,
        read_len      = config["simulation"]["read_length"],
        mflen         = config["simulation"]["mean_fragment_length"],
        sdev          = config["simulation"]["std_fragment_length"],
        fcov          = config["simulation"]["coverage"],
        seed          = replicate_seed,
        out_prefix    = lambda wc: (
            f"results/simulated/{wc.ref_id}/{wc.replicate}/{wc.mutated}/{wc.ref_id}"
        ),
    threads:
        config["simulation"]["threads"]
    log:
        "logs/simulation/{ref_id}/{replicate}_{mutated}.log",
    conda:
        "../envs/art.yaml"
    shell:
        """
        art_illumina \
            {params.profile_flags} \
            --paired \
            --in {input.ref} \
            --out {params.out_prefix}_ \
            --len {params.read_len} \
            --mflen {params.mflen} \
            --sdev {params.sdev} \
            --fcov {params.fcov} \
            --rndSeed {params.seed} \
            -sam \
            --noALN \
            &> {log}

        # ART names outputs <prefix>_1.fq / _2.fq — rename to _R1/_R2 convention
        gzip {params.out_prefix}_1.fq
        gzip {params.out_prefix}_2.fq
        mv {params.out_prefix}_1.fq.gz {output.r1}
        mv {params.out_prefix}_2.fq.gz {output.r2}
        """
