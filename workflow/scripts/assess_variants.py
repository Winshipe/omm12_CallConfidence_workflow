#!/usr/bin/env python3
"""
scripts/assess_variants.py
───────────────────────────
Compare variant calls (VCF format) against ground-truth mutation TSVs to
produce a per-mutation assessment table for one (scenario × replicate).

This script is called by the 'assess_variants' rule in assess.smk.  It is
run once per replicate; the results from all replicates are later combined by
the 'aggregate_replicates' rule into a single cross-replicate summary.

Output columns
──────────────
  ref_id          — reference ID (e.g. ref_A)
  seq_id          — sequence / contig ID within the reference
  position        — 1-based genomic position of the expected mutation
  ref_base        — expected reference allele
  alt_base        — expected alternative allele
  mutation_type   — transition | transversion
  detected        — True | False — was any variant called at this position?
  vcf_alt         — the ALT allele reported by breseq (or "." if not called)
  vcf_freq        — allele frequency from the VCF INFO/FORMAT field (or ".")
  vcf_quality     — QUAL score from the VCF (or ".")
  above_threshold — True | False — detected AND QUAL >= min_quality
  replicate       — replicate label (e.g. "rep1"), added for traceability

Snakemake injects
─────────────────
  snakemake.input.vcf             path to the breseq VCF for this replicate
  snakemake.input.ground_truth    list of mutation TSV paths for this replicate
  snakemake.output.tsv            path for the output assessment TSV
  snakemake.params.min_quality    minimum QUAL score to count as a true positive
  snakemake.params.replicate      replicate label string (e.g. "rep1")
  snakemake.log[0]                path for the log file
"""

import csv
import logging
import pandas as pd

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

#read the vcf file
found_variants = pd.read_csv(
    snakemake.input.vcf,\
    sep="\t",\
    comment="#",\
    names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'scenario_equal_mix']\
)

#read the tsv containing the ground truth
ground_truth = pd.read_csv(snakemake.input.ground_truth,sep="\t")
log.info(f"Loaded {ground_truth.shape[0]} expected mutations from {snakemake.input.ground_truth}")
ground_truth.columns = ["CHROM","POS","REF","ALT","mutation_type"] # rename columns to match the VCF names

found_variants = found_variants[(found_variants["QUAL"] >= snakemake.params.min_quality) & (found_variants["CHROM"] == ground_truth.CHROM[0])]
log.info(f"Parsed {found_variants.shape[0]} SNP calls from {snakemake.input.vcf}")

#here we do an outer join which preserves values in both the left and right tables and add an indicator column
#if a variant is found only in the left table (ie ground_truth) its a false negative
#if a variant is found only in the right table (ie the vcf) its a false positive
#and if its in both its a true positive
combined_variants = pd.merge(\
    ground_truth,\
    found_variants,\
    how="outer",\
    on=["CHROM","POS","REF","ALT"],\
    indicator=True\
    )[["CHROM","POS","REF","ALT","_merge"]]

combined_variants["Truthiness"] = combined_variants["_merge"].apply(lambda value: {"both":"TP", "right_only": "FP", "left_only":"FN"}[value])
combined_variants = combined_variants.drop('_merge', axis=1)

combined_variants["replicate"] = snakemake.params.replicate

combined_variants.to_csv(snakemake.output.tsv, sep="\t", index=False)


#Claude's attempt based on my previous code and prompts below.  Looks nice but is overly verbose when a few lines of pandas will suffice :)

# ── Parse VCF ─────────────────────────────────────────────────────────────────

# def parse_vcf(vcf_path):
#     """
#     Parse a VCF file and return a dict keyed by (chrom, position):
#       {(chrom, pos): {'alt': str, 'freq': str, 'quality': float | None}}

#     Only SNP records (single-base REF and single-base ALT) are retained.
#     For multi-allelic sites only the first ALT allele is used.

#     Allele frequency is resolved by checking (in order):
#       1. INFO field — AF=... or FREQ=...
#       2. FORMAT/sample column — AF sub-field (GATK-style GT:AD:AF)
#     """
#     calls = {}

#     with open(vcf_path) as fh:
#         for line in fh:
#             line = line.rstrip()

#             if line.startswith("##"):   # VCF meta-information — skip
#                 continue
#             if line.startswith("#CHROM"):  # column header — skip
#                 continue

#             fields = line.split("\t")
#             if len(fields) < 5:
#                 continue

#             chrom   = fields[0]
#             pos     = int(fields[1])
#             ref     = fields[3]
#             alt_raw = fields[4]
#             qual    = fields[5] if len(fields) > 5 else "."
#             info    = fields[7] if len(fields) > 7 else "."

#             # Only handle SNPs (len 1 REF and ALT)
#             alt = alt_raw.split(",")[0]
#             if len(ref) != 1 or len(alt) != 1 or alt in (".", "*"):
#                 continue

#             # Parse QUAL score
#             try:
#                 quality = float(qual)
#             except ValueError:
#                 quality = None

#             # Resolve allele frequency
#             freq      = "."
#             info_dict = {}
#             for token in info.split(";"):
#                 if "=" in token:
#                     k, v = token.split("=", 1)
#                     info_dict[k] = v

#             if "AF" in info_dict:
#                 freq = info_dict["AF"].split(",")[0]
#             elif "FREQ" in info_dict:
#                 freq = info_dict["FREQ"].split(",")[0]

#             # Fall back to FORMAT/sample column
#             if freq == "." and len(fields) >= 10:
#                 fmt_keys = fields[8].split(":")
#                 fmt_vals = fields[9].split(":")
#                 fmt      = dict(zip(fmt_keys, fmt_vals))
#                 if "AF" in fmt:
#                     freq = fmt["AF"].split(",")[0]

#             calls[(chrom, pos)] = {
#                 "alt":     alt,
#                 "freq":    freq,
#                 "quality": quality,
#             }

#     log.info(f"Parsed {len(calls)} SNP calls from {vcf_path}")
#     return calls


# # ── Load ground truth ─────────────────────────────────────────────────────────

# def load_ground_truth(tsv_path, ref_id):
#     """
#     Read a ground-truth mutation TSV and return a list of mutation dicts.

#     Expected columns: seq_id | position | ref_base | alt_base | mutation_type
#     """
#     mutations = []
#     with open(tsv_path) as fh:
#         reader = csv.DictReader(fh, delimiter="\t")
#         for row in reader:
#             mutations.append({
#                 "ref_id":        ref_id,
#                 "seq_id":        row["seq_id"],
#                 "position":      int(row["position"]),
#                 "ref_base":      row["ref_base"],
#                 "alt_base":      row["alt_base"],
#                 #"mutation_type": row["mutation_type"],
#             })
#     log.info(f"Loaded {len(mutations)} expected mutations from {tsv_path}")
#     return mutations


# # ── Main ──────────────────────────────────────────────────────────────────────

# min_quality = float(snakemake.params.min_quality)
# replicate   = snakemake.params.replicate

# calls = parse_vcf(snakemake.input.vcf)

# # Load ground truth from every contributing reference for this replicate.
# # The TSV path encodes the replicate (e.g. results/mutated/ref_A/rep1/ref_A.mutations.tsv),
# # so the ref_id is extracted from the third-to-last path component.
# expected_mutations = []
# for tsv_path in snakemake.input.ground_truth:
#     # Path structure: results/mutated/{ref_id}/{replicate}/{ref_id}.mutations.tsv
#     parts  = tsv_path.replace("\\", "/").split("/")
#     ref_id = parts[-3]   # third from the end is the ref_id directory
#     expected_mutations.extend(load_ground_truth(tsv_path, ref_id))

# # Build one output row per expected mutation
# rows = []
# for m in expected_mutations:
#     key  = (m["seq_id"], m["position"])
#     call = calls.get(key)

#     if call is None:
#         detected     = False
#         vcf_alt      = "."
#         vcf_freq     = "."
#         vcf_qual     = "."
#         above_thresh = False
#     else:
#         detected     = True
#         vcf_alt      = call["alt"]
#         vcf_freq     = call["freq"]
#         vcf_qual     = call["quality"] if call["quality"] is not None else "."
#         above_thresh = (
#             call["quality"] is not None
#             and call["quality"] >= min_quality
#         )

#     rows.append({
#         "ref_id":          m["ref_id"],
#         "seq_id":          m["seq_id"],
#         "position":        m["position"],
#         "ref_base":        m["ref_base"],
#         "alt_base":        m["alt_base"],
#         #"mutation_type":   m["mutation_type"],
#         "detected":        detected,
#         "vcf_alt":         vcf_alt,
#         "vcf_freq":        vcf_freq,
#         "vcf_quality":     vcf_qual,
#         "above_threshold": above_thresh,
#         "replicate":       replicate,
#     })

# n_detected = sum(r["detected"] for r in rows)
# n_above    = sum(r["above_threshold"] for r in rows)
# log.info(
#     f"Replicate {replicate}: "
#     f"{n_detected}/{len(rows)} expected mutations detected; "
#     f"{n_above} above quality threshold ({min_quality})"
# )

# fieldnames = [
#     "ref_id", "seq_id", "position", "ref_base", "alt_base", #"mutation_type",
#     "detected", "vcf_alt", "vcf_freq", "vcf_quality", "above_threshold",
#     "replicate",
# ]

# with open(snakemake.output.tsv, "w", newline="") as fh:
#     writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
#     writer.writeheader()
#     writer.writerows(rows)

# log.info(f"Assessment written to {snakemake.output.tsv}")
