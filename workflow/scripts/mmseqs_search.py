#!/usr/bin/env python3
"""
scripts/mmseqs_search.py
─────────────────────────
Search a query MMseqs2 database against one or more target databases.
Results from all databases are merged into a single m8 TSV.

Snakemake injects:
  snakemake.input.query_db      path to the MMseqs2 query DB directory
  snakemake.output.hits         merged m8-format TSV output
  snakemake.params.databases    list of target DB paths
  snakemake.params.sensitivity
  snakemake.params.evalue
  snakemake.params.tmp
  snakemake.threads
  snakemake.log[0]
"""

import subprocess
import sys
import logging
import os
import tempfile
from pathlib import Path

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


def run(cmd, **kwargs):
    log.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        log.error(result.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    log.info(result.stdout)
    return result


query_db  = snakemake.input.query_db
databases = snakemake.params.databases
sens      = snakemake.params.sensitivity
evalue    = snakemake.params.evalue
tmp_dir   = snakemake.params.tmp
threads   = str(snakemake.threads)
out_hits  = snakemake.output.hits

os.makedirs(tmp_dir, exist_ok=True)
partial_files = []

for db in databases:
    db_name = Path(db).name
    result_db = os.path.join(tmp_dir, f"result_{db_name}")
    result_m8  = os.path.join(tmp_dir, f"result_{db_name}.m8")
    tmp_sub    = os.path.join(tmp_dir, db_name)
    os.makedirs(tmp_sub, exist_ok=True)

    # Search
    run([
        "mmseqs", "search",
        f"{query_db}/queryDB",
        db,
        result_db,
        tmp_sub,
        "-s", str(sens),
        "-e", str(evalue),
        "--threads", threads,
        "--search-type", "3",   # nucleotide vs nucleotide
    ])

    # Convert to m8
    run([
        "mmseqs", "convertalis",
        f"{query_db}/queryDB",
        db,
        result_db,
        result_m8,
        "--format-mode", "0",
    ])

    partial_files.append(result_m8)
    log.info(f"Finished search against {db}")

# Merge all per-database m8 files
with open(out_hits, "w") as fout:
    fout.write("query\ttarget\tpident\talign_len\tmismatches\tgap_opens\t"
               "q_start\tq_end\tt_start\tt_end\tevalue\tbitscore\tdatabase\n")
    for m8, db in zip(partial_files, databases):
        db_name = Path(db).name
        with open(m8) as fin:
            for line in fin:
                fout.write(line.rstrip() + f"\t{db_name}\n")

log.info(f"Merged results written to {out_hits}")
