#!/usr/bin/env python3
"""
scripts/blend_reads.py
───────────────────────
Blend simulated reads from multiple references at user-defined abundances
and mutated fractions into a single paired-end FASTQ pair.

Algorithm
─────────
For each contribution i:
  n_total_i   = round(total_reads * normalised_abundance_i)
  n_mutated_i = round(n_total_i   * mutated_fraction_i)
  n_wild_i    = n_total_i - n_mutated_i

The script randomly subsamples exactly that many read *pairs* from the
appropriate FASTQ.gz files and writes them to the output files.

Snakemake injects:
  snakemake.input            — flat dict of FASTQ paths (contrib_i_{un}mutated_{R1,R2})
  snakemake.output.r1 / .r2
  snakemake.params.scenario_cfg  — list of {ref_id, mutated_fraction, abundance}
  snakemake.params.total_reads
  snakemake.params.seed
  snakemake.log[0]
"""

import gzip
import logging
import math
import random
import sys
from itertools import islice
from pathlib import Path

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


# ── helpers ───────────────────────────────────────────────────────────────────

def fastq_records(path):
    """Yield (name, seq, plus, qual) tuples from a (possibly gzipped) FASTQ."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            name = fh.readline()
            if not name:
                break
            seq  = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            yield name, seq, plus, qual


def count_fastq_records(path):
    """Fast record count without parsing — just counts '@' header lines."""
    opener = gzip.open if str(path).endswith(".gz") else open
    count = 0
    with opener(path, "rt") as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            fh.readline()  # seq
            fh.readline()  # +
            fh.readline()  # qual
            count += 1
    return count


def reservoir_sample(path, k, rng):
    """Return k randomly chosen (name, seq, plus, qual) tuples."""
    reservoir = []
    for i, record in enumerate(fastq_records(path)):
        if i < k:
            reservoir.append(record)
        else:
            j = rng.randint(0, i)
            if j < k:
                reservoir[j] = record
    return reservoir


def write_fastq_records(records, fh):
    for name, seq, plus, qual in records:
        fh.write(name)
        fh.write(seq)
        fh.write(plus)
        fh.write(qual)


# ── main ─────────────────────────────────────────────────────────────────────

rng           = random.Random(snakemake.params.seed)
scenario_cfg  = snakemake.params.scenario_cfg
total_reads   = int(snakemake.params.total_reads)

# Normalise abundances
total_abundance = sum(c["abundance"] for c in scenario_cfg)
norm_abundances = [c["abundance"] / total_abundance for c in scenario_cfg]

# Compute per-contribution read counts
contribution_counts = []
running = 0
for i, (c, norm_ab) in enumerate(zip(scenario_cfg, norm_abundances)):
    if i < len(scenario_cfg) - 1:
        n = round(total_reads * norm_ab)
    else:
        n = total_reads - running     # absorb rounding remainder
    running += n
    n_mut   = round(n * c["mutated_fraction"])
    n_wild  = n - n_mut
    contribution_counts.append((n_wild, n_mut))
    log.info(
        f"Contribution {i} (ref={c['ref_id']}): "
        f"{n_wild} unmutated + {n_mut} mutated reads"
    )

# Sample reads from each contribution
all_r1, all_r2 = [], []

for i, (c, (n_wild, n_mut)) in enumerate(zip(scenario_cfg, contribution_counts)):
    ref = c["ref_id"]
    for mutated, n in (("unmutated", n_wild), ("mutated", n_mut)):
        if n == 0:
            continue
        key_r1 = f"contrib_{i}_{mutated}_R1"
        key_r2 = f"contrib_{i}_{mutated}_R2"
        path_r1 = snakemake.input[key_r1]
        path_r2 = snakemake.input[key_r2]

        log.info(f"  Sampling {n} pairs from {Path(path_r1).name} ({mutated})")
        sampled_r1 = reservoir_sample(path_r1, n, rng)
        sampled_r2 = reservoir_sample(path_r2, n, rng)

        if len(sampled_r1) < n:
            log.warning(
                f"  Requested {n} reads but only {len(sampled_r1)} available "
                f"in {path_r1}. Using all available."
            )

        all_r1.extend(sampled_r1)
        all_r2.extend(sampled_r2)

# Shuffle so reads are not ordered by source
combined = list(zip(all_r1, all_r2))
rng.shuffle(combined)
all_r1, all_r2 = zip(*combined) if combined else ([], [])

log.info(f"Writing {len(all_r1)} read pairs to output")

with gzip.open(snakemake.output.r1, "wt") as fh:
    write_fastq_records(all_r1, fh)

with gzip.open(snakemake.output.r2, "wt") as fh:
    write_fastq_records(all_r2, fh)

log.info("Done.")
