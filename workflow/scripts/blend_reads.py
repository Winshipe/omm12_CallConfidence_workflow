#!/usr/bin/env python3
"""
scripts/blend_reads.py
───────────────────────
Blend simulated reads from multiple references into a single paired-end
FASTQ pair.

Algorithm
─────────
Every read pair from every input FASTQ pair (the unmutated and mutated
simulations of each contribution) is used — no subsampling.  R1 and R2 records
are read in lockstep so mates stay together, then the pairs are shuffled
(seeded) so reads are not ordered by source, and written to the outputs.

Because all reads are used, the scenario's abundance and mutated_fraction
values do not change which reads are written; they are only logged.

Snakemake injects:
  snakemake.input            — flat dict of FASTQ paths (contrib_i_{un}mutated_{R1,R2})
  snakemake.output.r1 / .r2
  snakemake.params.scenario_cfg  — list of {ref_id, mutated_fraction, abundance}
  snakemake.params.seed
  snakemake.log[0]
"""

import gzip
import logging
import random
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


def read_pairs(path_r1, path_r2):
    """Return a list of (R1 record, R2 record) tuples, read in lockstep."""
    r1 = list(fastq_records(path_r1))
    r2 = list(fastq_records(path_r2))
    if len(r1) != len(r2):
        raise ValueError(
            f"{path_r1} has {len(r1)} reads but {path_r2} has {len(r2)}; "
            "paired FASTQ files must contain the same number of reads"
        )
    return list(zip(r1, r2))


def write_fastq_records(records, fh):
    for name, seq, plus, qual in records:
        fh.write(name)
        fh.write(seq)
        fh.write(plus)
        fh.write(qual)


# ── main ─────────────────────────────────────────────────────────────────────

rng          = random.Random(snakemake.params.seed)
scenario_cfg = snakemake.params.scenario_cfg

pairs = []
for i, c in enumerate(scenario_cfg):
    log.info(
        f"Contribution {i} (ref={c['ref_id']}, abundance={c['abundance']}, "
        f"mutated_fraction={c['mutated_fraction']}): using all reads"
    )
    for mutated in ("unmutated", "mutated"):
        path_r1 = snakemake.input[f"contrib_{i}_{mutated}_R1"]
        path_r2 = snakemake.input[f"contrib_{i}_{mutated}_R2"]
        contrib_pairs = read_pairs(path_r1, path_r2)
        log.info(f"  {len(contrib_pairs)} pairs from {Path(path_r1).name} ({mutated})")
        pairs.extend(contrib_pairs)

# Shuffle so reads are not ordered by source; mates move together
rng.shuffle(pairs)

log.info(f"Writing {len(pairs)} read pairs to output")

with gzip.open(snakemake.output.r1, "wt") as fh:
    write_fastq_records((r1 for r1, _ in pairs), fh)

with gzip.open(snakemake.output.r2, "wt") as fh:
    write_fastq_records((r2 for _, r2 in pairs), fh)

log.info("Done.")
