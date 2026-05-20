#!/usr/bin/env python3
"""
scripts/annotate_repeats.py
────────────────────────────
Post-process raw MMseqs2 m8 hits into a clean annotation TSV.

Output columns:
  seq_id  start  end  strand  element_name  database  evalue  bitscore

Snakemake injects:
  snakemake.input.hits
  snakemake.output.tsv
  snakemake.wildcards.ref_id
  snakemake.log[0]
"""

import logging
import csv

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

EXPECTED_COLS = [
    "query", "target", "pident", "align_len", "mismatches", "gap_opens",
    "q_start", "q_end", "t_start", "t_end", "evalue", "bitscore", "database",
]


def parse_element_name(target_field):
    """Extract a human-readable element name from the MMseqs2 target string."""
    # Many TE databases encode the element name as the first token before '#'
    # or '|'.  Fall back to the full target string.
    for sep in ("#", "|", " "):
        if sep in target_field:
            return target_field.split(sep)[0].strip()
    return target_field


annotations = []

with open(snakemake.input.hits) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        q_start = int(row["q_start"])
        q_end   = int(row["q_end"])

        if q_start <= q_end:
            start, end, strand = q_start, q_end, "+"
        else:
            start, end, strand = q_end, q_start, "-"

        annotations.append({
            "seq_id":       row["query"],
            "start":        start,
            "end":          end,
            "strand":       strand,
            "element_name": parse_element_name(row["target"]),
            "database":     row["database"],
            "evalue":       row["evalue"],
            "bitscore":     row["bitscore"],
        })

# Sort by seq_id then start position
annotations.sort(key=lambda r: (r["seq_id"], r["start"]))

with open(snakemake.output.tsv, "w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=["seq_id", "start", "end", "strand",
                    "element_name", "database", "evalue", "bitscore"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(annotations)

log.info(f"Wrote {len(annotations)} annotations to {snakemake.output.tsv}")
