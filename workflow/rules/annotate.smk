"""
rules/annotate.smk
==================
Annotate reference genomes for transposable elements (TEs), insertion sequences
(IS elements), and other repeats using MMseqs2.

Pipeline overview
-----------------
1. get_amino_acid_seqs  – call genes with Prodigal (AA + NT FASTA outputs)
2. mmseqs_search        – easy-search each reference against every configured
                          repeat/TE database; produces one hit table per
                          (ref_id × db) pair
3. cluster_genes_by_ani – cluster all reference genes at 90 % ANI with
                          MMseqs2 easy-cluster to identify shared / homologous
                          sequences across references
4. find_homologs        – parse cluster output → BED-style TSV of genes that
                          appear in more than one reference ("challenging"
                          regions)
5. annotation_hits      – aggregate target that collects all per-db hit tables
                          and both TSV outputs for a given ref_id

Configuration keys expected under config["annotation"]
-------------------------------------------------------
  mmseqs_databases : dict[str, str]   # name → path to each MMseqs2 target DB
  split_mem_limit  : str  (optional)  # e.g. "164G"; defaults to "164G"
  max_memory       : str  (optional)  # e.g. "196G"; defaults to "196G"
  threads          : int  (optional)  # CPU threads; defaults to 16

config["references"] : dict[str, str]  # ref_id → path to reference FASTA
"""

import os


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_ref_fasta(wildcards):
    """Return the reference FASTA path for a given ref_id wildcard."""
    return config["references"][wildcards.ref_id]


# ---------------------------------------------------------------------------
# Aggregate rule – collect all annotation outputs for every reference
# ---------------------------------------------------------------------------

rule annotation_hits:
    """
    Pseudo-target that forces production of:
      - one mmseqs hit table per (ref_id × database) combination
      - the cross-reference homolog BED TSV
      - the 'challenging regions' TSV for the given ref_id
    """
    input:
        # One hit file per (db, ref_id) combination
        expand(
            "results/annotation/{ref_id}/{db}_hits",
            db=list(config["annotation"]["mmseqs_databases"].keys()),
            ref_id=list(config["references"].keys()),
        ),
        "results/annotation/{ref_id}/homologs.tsv",
        "results/annotation/{ref_id}/{ref_id}.challenging.tsv",


# ---------------------------------------------------------------------------
# Rule 1 – predict protein-coding genes with Prodigal
# ---------------------------------------------------------------------------

rule get_amino_acid_seqs:
    """
    Run Prodigal in metagenomic mode (-p meta) on a reference FASTA to produce:
      - {ref_id}.faa : translated amino-acid sequences (used by mmseqs_search)
      - {ref_id}.fna : nucleotide gene sequences
    """
    input:
        ref=get_ref_fasta,
    output:
        out_aa="results/annotation/{ref_id}/{ref_id}.faa",
        out_nt="results/annotation/{ref_id}/{ref_id}.fna",
        out_gff="results/annotation/{ref_id}/{ref_id}.gff",
    log:
        "logs/annotation/{ref_id}.prodigal.log",
    conda:
        "../envs/prodigal.yaml"
    shell:
        # -p meta  : metagenomic / anonymous mode (no training step required)
        # -d       : write nucleotide sequences to out_nt
        # -a       : write amino-acid sequences to out_aa
        # -i       : input FASTA
        """
        prodigal -p meta \
            -i {input.ref} \
            -a {output.out_aa} \
            -d {output.out_nt} \
            -f gff \
            -o {output.out_gff} \
            &> {log}
        """


# ---------------------------------------------------------------------------
# Rule 2 – search predicted proteins against repeat / TE databases
# ---------------------------------------------------------------------------

rule mmseqs_search:
    """
    Use MMseqs2 easy-search to compare predicted AA sequences for one reference
    against one repeat/TE database.  Outputs a tabular hit file in m8 format.

    Memory and thread parameters are read from config with safe defaults.
    """
    input:
        # Prodigal-predicted amino-acid sequences for this reference
        query=rules.get_amino_acid_seqs.output.out_aa,
        # Path to the pre-built MMseqs2 target database
        db_path=lambda wildcards: config["annotation"]["mmseqs_databases"][wildcards.db],
    output:
        hits="results/annotation/{ref_id}/{db}_hits",
    log:
        "logs/annotation/{ref_id}.{db}.mmseqs_search.log",
    conda:
        "../envs/mmseqs2.yaml"
    params:
        split_mem_limit=config["annotation"].get("split_mem_limit", "164G"),
    resources:
        # mem_mb should be an integer (MB); convert a string like "196G" → 200704
        mem_mb=int(config["annotation"].get("max_memory", "196").rstrip("G")) * 1024,
        cpus_per_task=int(config["annotation"].get("threads", 16)),
    shell:
        """
        dbpath={input.db_path}
        if [[ dbpath == *.f*a ]]; then
            dbpath_temp="${{filename%.*}}"
            mmseqs createdb $dbpath $dbpath_temp
            dbpath=$dbpath_temp
        fi
        mmseqs easy-search \
            --split-memory-limit {params.split_mem_limit} \
            --threads {resources.cpus_per_task} \
            {input.query} \
            $dbpath \
            {output.hits} \
            /tmp/$SLURM_JOB_ID
            &> {log}

        """


# ---------------------------------------------------------------------------
# Rule 3 – cluster all reference genes by ANI to find cross-reference homologs
# ---------------------------------------------------------------------------

rule cluster_genes_by_ani:
    """
    Concatenate all reference nucleotide gene FASTAs and cluster at 90 % ANI
    with MMseqs2 easy-cluster.  The resulting cluster TSV (two-column:
    representative → member) is parsed by find_homologs to identify genes
    shared across multiple references.
    """
    input:
        # Collect nucleotide gene FASTAs for every reference
        fa_files=expand(
            "results/annotation/{ref_id}/{ref_id}.fna",
            ref_id=list(config["references"].keys()),
        ),
    output:
        clusters="results/annotation/clusters.txt",
        # easy-cluster also writes _rep_seq.fasta and _all_seqs.fasta;
        # declare the prefix directory so Snakemake can track the tmp files
        concat_fa=temp("results/annotation/concat_refs.fa"),
    log:
        "logs/annotation/cluster_genes_by_ani.log",
    conda:
        "../envs/mmseqs2.yaml"
    params:
        out_prefix="results/annotation/clusters",
        min_seq_id=0.9,
        min_aln_len=50,
    resources:
        mem_mb=int(config["annotation"].get("max_memory", "196").rstrip("G")) * 1024,
    shell:
        """
        cat {input.fa_files} > {output.concat_fa}
        mmseqs easy-cluster \
            --min-seq-id {params.min_seq_id} \
            --min-aln-len {params.min_aln_len} \
            {output.concat_fa} \
            {params.out_prefix} \
            /tmp/$SLURM_JOB_ID \
            &> {log}
        mv {params.out_prefix}_cluster.tsv {output.clusters}
        """


# ---------------------------------------------------------------------------
# Rule 4 – identify homologous / challenging regions across references
# ---------------------------------------------------------------------------

rule find_homologs:
    """
    Parse the MMseqs2 cluster TSV to find genes present in more than one
    reference (homologs).  For each homologous gene, look up its genomic
    coordinates from the Prodigal FASTA headers and write a BED-style TSV:

        chrom  start  stop  name  score  strand

    Output is used downstream to flag 'challenging' regions that may confound
    short-read mapping or variant calling.
    """
    input:
        # Two-column cluster file: representative<TAB>member
        clusters="results/annotation/clusters.txt",
        # Nucleotide FASTAs for all references (headers carry coordinate info)
        fa_files = "results/annotation/{ref_id}/{ref_id}.fna",
#        fa_files=expand(
#            "results/annotation/{ref_id}/{ref_id}.fna",
#            ref_id=list(config["references"].keys()),
#        ),
    output:
        ref_hits="results/annotation/{ref_id}/{ref_id}.challenging.tsv",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/find_homologs.py"


