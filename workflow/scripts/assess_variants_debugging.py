#!/usr/bin/env python3
"""
scripts/assess_variants.py
───────────────────────────
Compare variant calls (VCF format) against ground-truth mutation TSVs
to produce a per-mutation assessment table.

Output columns
──────────────
  ref_id          reference sequence ID
  position        1-based position of the expected mutation
  ref_base        expected reference allele
  alt_base        expected alternative allele
  mutation_type   transition | transversion
  detected        True | False — whether the caller called a variant here
  vcf_alt         the ALT allele reported in the VCF (or "." if not called)
  vcf_freq        allele frequency from the VCF INFO/FORMAT field (or ".")
  vcf_quality     QUAL score from the VCF (or ".")
  above_threshold True | False — QUAL >= min_quality AND detected

Snakemake injects:
  snakemake.input.vcf             path to variant calls VCF
  snakemake.input.ground_truth    dict {ref_id: mutations_tsv_path}
  snakemake.output.tsv
  snakemake.params.min_quality
  snakemake.log[0]
"""

import csv
import logging
import sys
from collections import defaultdict

class Myinput:
    def __init__(self, first, second):
        self.vcf = first
        self.ground_truth = [second]
class Myoutput:
    def __init__(self, first):
        self.tsv = first
class Myparams:
    def __init__(self, first):
        self.min_quality = first

class Snakemake:
    def __init__(self, vcf_path, gt_path, outpath, minq, log_path):
        self.input = Myinput(vcf_path, gt_path)
        self.output = Myoutput(outpath)
        self.params = Myparams(float(minq))
        self.log = [log_path]

import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process a VCF file with ground truth TSV files and quality filtering."
    )

    # Required VCF file path
    parser.add_argument(
        "vcf_path",
        type=str,
        help="Path to the input VCF file"
    )

    # One or more ground truth TSV files
    parser.add_argument(
        "ground_truth",
        type=str,
        nargs="+",  # allows arbitrary number (at least 1)
        help="Paths to ground truth TSV files"
    )

    # Minimum quality (float)
    parser.add_argument(
        "--min-quality",
        type=float,
        required=True,
        help="Minimum quality threshold"
    )

    # Log file path
    parser.add_argument(
        "--log",
        type=str,
        required=True,
        help="Path to the log file"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    print("VCF file:", args.vcf_path)
    print("Ground truth TSVs:", args.ground_truth)
    print("Min quality:", args.min_quality)
    print("Log file:", args.log)
snakemake = Snakemake(args[],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


# ── Parse VCF ─────────────────────────────────────────────────────────────────

def parse_vcf(vcf_path):
    """
    Parse a VCF file and return a dict:
      {(chrom, position): {'alt': str, 'freq': str, 'quality': float | None}}

    Only SNP records (single-base REF and ALT) are considered; multi-allelic
    sites use the first ALT allele. Allele frequency is resolved by checking
    (in order): INFO/AF, INFO/FREQ, then the AF sub-field of the first sample
    FORMAT column.
    """
    calls = {}

    with open(vcf_path) as fh:
        format_keys = []  # populated when we hit the #CHROM header line

        for line in fh:
            line = line.rstrip()

            # Skip meta-information lines
            if line.startswith("##"):
                continue

            # Column-header line — nothing to parse, but marks end of header
            if line.startswith("#CHROM"):
                continue

            fields = line.split("\t")
            if len(fields) < 5:
                continue

            chrom   = fields[0]
            pos     = int(fields[1])
            ref     = fields[3]
            alt_raw = fields[4]
            qual    = fields[5] if len(fields) > 5 else "."
            info    = fields[7] if len(fields) > 7 else "."

            # Only handle SNPs (single-base substitutions)
            alt = alt_raw.split(",")[0]  # take first ALT for multi-allelic
            if len(ref) != 1 or len(alt) != 1 or alt in (".", "*"):
                continue

            # Parse QUAL
            try:
                quality = float(qual)
            except ValueError:
                quality = None

            # Resolve allele frequency -----------------------------------------
            # 1. Try INFO field (AF=... or FREQ=...)
            freq = "."
            info_dict = {}
            for token in info.split(";"):
                if "=" in token:
                    k, v = token.split("=", 1)
                    info_dict[k] = v

            if "AF" in info_dict:
                freq = info_dict["AF"].split(",")[0]  # first value for multi-allelic
            elif "FREQ" in info_dict:
                freq = info_dict["FREQ"].split(",")[0]

            # 2. Fall back to FORMAT/sample column (e.g. GATK-style GT:AD:AF)
            if freq == "." and len(fields) >= 10:
                fmt_keys  = fields[8].split(":")
                fmt_vals  = fields[9].split(":")
                fmt       = dict(zip(fmt_keys, fmt_vals))
                if "AF" in fmt:
                    freq = fmt["AF"].split(",")[0]

            calls[(chrom, pos)] = {
                "alt":     alt,
                "freq":    freq,
                "quality": quality,
            }

    log.info(f"Parsed {len(calls)} SNP calls from {vcf_path}")
    return calls


# ── Load ground truth ─────────────────────────────────────────────────────────

def load_ground_truth(tsv_path, ref_id):
    """Return list of mutation dicts from a ground-truth TSV."""
    mutations = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mutations.append({
                "ref_id":        ref_id,
                "seq_id":        row["seq_id"],
                "position":      int(row["position"]),
                "ref_base":      row["ref_base"],
                "alt_base":      row["alt_base"],
                "mutation_type": row["mutation_type"],
            })
    log.info(f"Loaded {len(mutations)} expected mutations from {tsv_path}")
    return mutations


# ── Main ──────────────────────────────────────────────────────────────────────

min_quality = float(snakemake.params.min_quality)
calls       = parse_vcf(snakemake.input.vcf)

expected_mutations = []
for tsv_path in snakemake.input.ground_truth:
    ref_id = tsv_path.split("/")[-2] #path should look like results/mutated/ref_A/ref_A.mutations.tsv
    expected_mutations.extend(load_ground_truth(tsv_path, ref_id))

rows = []
for m in expected_mutations:
    key  = (m["seq_id"], m["position"])
    call = calls.get(key)

    if call is None:
        detected    = False
        vcf_alt     = "."
        vcf_freq    = "."
        vcf_qual    = "."
        above_thresh = False
    else:
        detected    = True
        vcf_alt     = call["alt"]
        vcf_freq    = call["freq"]
        vcf_qual    = call["quality"] if call["quality"] is not None else "."

        above_thresh = (
            call["quality"] is not None
            and call["quality"] >= min_quality
        )

    rows.append({
        "ref_id":          m["ref_id"],
        "seq_id":          m["seq_id"],
        "position":        m["position"],
        "ref_base":        m["ref_base"],
        "alt_base":        m["alt_base"],
        "mutation_type":   m["mutation_type"],
        "detected":        detected,
        "vcf_alt":         vcf_alt,
        "vcf_freq":        vcf_freq,
        "vcf_quality":     vcf_qual,
        "above_threshold": above_thresh,
    })

n_detected = sum(r["detected"] for r in rows)
n_above    = sum(r["above_threshold"] for r in rows)
log.info(
    f"Summary: {n_detected}/{len(rows)} expected mutations detected; "
    f"{n_above} above quality threshold ({min_quality})"
)

fieldnames = [
    "ref_id", "seq_id", "position", "ref_base", "alt_base", "mutation_type",
    "detected", "vcf_alt", "vcf_freq", "vcf_quality", "above_threshold",
]

with open(snakemake.output.tsv, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

log.info(f"Assessment written to {snakemake.output.tsv}")
