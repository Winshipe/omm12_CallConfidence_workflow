# CallConfidence - Variant-Calling Benchmarking Pipeline
# CallConfidence - Variant-Calling Benchmarking Pipeline

A modular Snakemake workflow for benchmarking variant callers on simulated
metagenomic reads from mixed, mutated reference sequences.

---

## Overview

```
references (FASTA)
       │
       ├─► [mutate]          → mutated FASTA + ground-truth TSV
       │
       ├─► [annotate]        → repeat / TE annotation TSV  (MMseqs2)
       │
       ├─► [simulate_reads]  → simulated MiSeq reads (ART)
       │       (unmutated & mutated versions of each reference)
       │
       ├─► [blend_reads]     → per-scenario blended FASTQ pairs
       │       (user-defined abundances & mutated fractions)
       │
       ├─► [variant_calling] → VCF file
       ├─► [variant_calling] → VCF file
       │
       └─► [assess]          → per-mutation detection TSV
```

---

## Quick start

### 1. Install Snakemake

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake>=7
conda activate snakemake
```

If you're running snakemake on a cluster with a scheduler, you will likely need to add an executor plugin to allow snakemake to talk to the cluster's scheduler.  https://snakemake.github.io/snakemake-plugin-catalog/index.html  If you're unsure which to add, talk to your cluster's support staff

### 2. Configure your run

Edit `config/config.yaml`:

- Add your reference FASTA paths under `references`.
- Point `annotation.mmseqs_databases` at your MMseqs2-formatted TE / IS databases.
- Define blend scenarios under `scenarios`.
- Optionally supply real MiSeq reads for empirical error-profile learning under
  `simulation.empirical_reads_R1/R2`.

### 3. Dry run

```bash
snakemake --snakefile workflow/Snakefile \
          --configfile config/config.yaml \
          --use-conda \
          --cores 1 \
          --dry-run
```

### 4. Run

```bash
snakemake --snakefile workflow/Snakefile \
          --configfile config/config.yaml \
          --use-conda \
          --cores 16
```

On a cluster with SLURM:

```bash
snakemake --snakefile workflow/Snakefile \
          --configfile config/config.yaml \
          --use-conda \
          --executor slurm \
          --jobs 50 \
          --default-resources slurm_partition=standard mem_mb=8000
```

---

## Directory structure

```
snakemake_workflow/
├── config/
│   └── config.yaml             ← edit this
├── resources/
│   ├── references/             ← put reference FASTAs here
│   └── databases/              ← put MMseqs2 DBs here
├── workflow/
│   ├── Snakefile               ← main entry point
│   ├── envs/
│   │   ├── art.yaml
│   │   ├── variant_calling.yaml
│   │   ├── mmseqs2.yaml
│   │   ├── mutate.yaml
│   │   └── python.yaml
│   ├── rules/
│   │   ├── annotate.smk
│   │   ├── assess.smk
│   │   ├── blend_reads.smk
│   │   ├── mutate.smk
│   │   ├── simulate_reads.smk
│   │   └── variant_calling.smk
│   └── scripts/
│       ├── annotate_repeats.py
│       ├── assess_variants.py
│       ├── blend_reads.py
│       ├── mmseqs_search.py
│       └── mutate_reference.py
└── results/                    ← created during the run
    ├── annotation/
    ├── assessment/
    ├── blended/
    ├── variant_calling/
    ├── mutated/
    └── simulated/
```

---

## Modules

### mutate

Applies a user-specified substitution model to each reference.

| Model | Config key | Notes |
|-------|-----------|-------|
| Jukes-Cantor | `jukes_cantor` | Equal rates among all substitution classes |
| Tamura & Nei 1993 | `tamura_nei` | Separate Ti/Tv rates, unequal base frequencies |

Key config keys: `mutation.model`, `mutation.substitution_rate`, `mutation.kappa`,
`mutation.gc_freq`, `mutation.seed`.

### annotate

This module uses mmseqs2 to search for hits in the databases provided by the user.
I use the Phrog and TnCentral+ISfinder databases.  The databases for mmseqs need to be 
downloaded from their respective websites and formatted but this is relatively simple using 
mmseqs createdb command (`mmseqs2 createdb sequences.fasta dbname`; described in the mmseqs2 user manual). This module also uses prodigal to 
predict genes and then uses mmseqs2's clustering to group the genes by 90% ANI to identify
close homologs.  This will then output a tsv (Tab Separated Values) with columns for chromosome, 
start, stop and information about the "dangerous" region (ie is it a TE, IS, phage or duplicated).

### simulate_reads

Uses `art_illumina` to simulate paired-end MiSeq reads.  If real reads are
supplied as fastqs with the paths placed in the config file, ART's `art_profiler_illumina` is invoked
first to build a custom quality-score profile.

### blend_reads

Combines simulated reads according to differerent user specified relative abundances and minor allele fractions (called "scenarios", also in the config file).  The user can specify multiple scenarios to assess in the config file.


### variant_calling

Maps the output from blend_reads to the reference genome(s) using BWA and then calls variants with GATK HaplotypeCaller 
Maps the output from blend_reads to the reference genome(s) using BWA and then calls variants with GATK HaplotypeCaller 

### assess

Parses the VCF output and cross-references each expected
Parses the VCF output and cross-references each expected
mutation with the ground truth (the record of which mutations were generated), recording detection status and quality scores.

---

## Outputs

| Path | Description |
|------|-------------|
| `results/mutated/{ref_id}/{ref_id}.mutated.fasta` | Mutated reference sequence |
| `results/mutated/{ref_id}/{ref_id}.mutations.tsv` | Ground-truth mutation table |
| `results/annotation/{ref_id}/{ref_id}.repeats.tsv` | Repeat/TE annotation |
| `results/simulated/{ref_id}/{mutated}/{ref_id}_R{1,2}.fastq.gz` | Simulated reads |
| `results/blended/{scenario}/{scenario}_R{1,2}.fastq.gz` | Blended scenario reads |
| `results/variant_calling/{scenario}/output/output.vcf` | Output from variant caller (GATK) |
| `results/assessment/{scenario}_assessment.tsv` | Per-mutation assessment |

The assessment TSV contains the following columns:

```
ref_id  seq_id  position  ref_base  alt_base  mutation_type
detected  vcf_alt  vcf_freq  vcf_quality  above_threshold
detected  vcf_alt  vcf_freq  vcf_quality  above_threshold
```
