"""
rules/simulate_reads.smk  — simulate MiSeq paired-end reads with ART

If empirical reads are supplied in the config, ART is first run in
profile-learning mode to build a custom error model; otherwise the
built-in MSv3 profile is used.
"""


def get_ref_fasta(wildcards):
    return config["references"][wildcards.ref_id]


def art_ref_fasta(wildcards):
    """Return mutated or unmutated FASTA depending on the mutated flag."""
    if wildcards.mutated == "mutated":
        return f"results/mutated/{wildcards.ref_id}/{wildcards.ref_id}.mutated.fasta"
    return config["references"][wildcards.ref_id]


# ------------------------------------------------------------------
# Optional: learn empirical error profile from real reads
# ------------------------------------------------------------------

rule art_learn_profile:
    """
    Use ART's built-in profiler to derive a custom quality-score
    profile from real MiSeq reads (only executed when empirical reads
    are configured).
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
    """Return ART flags selecting either the empirical or built-in profile."""
    if config["simulation"]["empirical_reads_R1"] is not None:
        return (
            "--qprof1 results/art_profile/empirical_R1.txt "
            "--qprof2 results/art_profile/empirical_R2.txt"
        )
    return f"-ss {config['simulation']['platform']}"


# ------------------------------------------------------------------
# Simulate reads for each reference × {unmutated, mutated}
# ------------------------------------------------------------------

rule simulate_reads:
    """
    Run ART to produce paired-end FASTQ reads for one reference
    (either unmutated or mutated version).
    """
    input:
        ref=art_ref_fasta,
        # Implicitly depend on the profile if empirical reads are given
        profile=lambda wc: (
            ["results/art_profile/empirical_R1.txt",
             "results/art_profile/empirical_R2.txt"]
            if config["simulation"]["empirical_reads_R1"] is not None
            else []
        ),
    output:
        r1="results/simulated/{ref_id}/{mutated}/{ref_id}_R1.fastq.gz",
        r2="results/simulated/{ref_id}/{mutated}/{ref_id}_R2.fastq.gz",
    params:
        profile_flags=art_profile_flags,
        read_len=config["simulation"]["read_length"],
        mflen=config["simulation"]["mean_fragment_length"],
        sdev=config["simulation"]["std_fragment_length"],
        fcov=config["simulation"]["coverage"],
        seed=config["simulation"]["seed"],
        out_prefix=lambda wc: (
            f"results/simulated/{wc.ref_id}/{wc.mutated}/{wc.ref_id}"
        ),
    threads:
        config["simulation"]["threads"]
    log:
        "logs/simulation/{ref_id}_{mutated}.log",
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
        # ART names files <prefix>_1.fastq.gz / _2.fastq.gz — rename to _R1/_R2
        gzip {params.out_prefix}_1.fq
        gzip {params.out_prefix}_2.fq
        mv {params.out_prefix}_1.fq.gz {output.r1}
        mv {params.out_prefix}_2.fq.gz {output.r2}
        """
