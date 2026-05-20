"""
rules/assess.smk — compare breseq calls to ground-truth mutations
==================================================================

This sub-workflow has two stages:

Stage 1 — per-replicate assessment
------------------------------------
For each (scenario × replicate), compare the breseq VCF output against the
ground-truth mutation TSVs for that specific replicate.  This produces one
TSV per (scenario × replicate):

  results/assessment/{scenario}/{replicate}_assessment.tsv

Columns:
  ref_id | seq_id | position | ref_base | alt_base | mutation_type |
  detected | vcf_alt | vcf_freq | vcf_quality | above_threshold | replicate

Stage 2 — cross-replicate aggregation
---------------------------------------
All per-replicate TSVs for a scenario are concatenated into a single summary
table, with a 'replicate' column added so results can be grouped and compared:

  results/assessment/{scenario}_all_replicates.tsv

This aggregated file is the primary output used for downstream statistical
analysis (e.g. sensitivity as a function of mutated_fraction, or consistency
of calls across replicates).
"""

import os


# ---------------------------------------------------------------------------
# Helper: derive the replicate list from config
# ---------------------------------------------------------------------------

REPLICATES = [f"rep{i+1}" for i in range(config["replicates"])]


# ---------------------------------------------------------------------------
# Helper: ground-truth TSVs for one (scenario × replicate) pair
# ---------------------------------------------------------------------------

def scenario_ground_truth_tsvs(wildcards):
    """
    Return a list of ground-truth mutation TSV paths for every reference that
    contributes reads to this scenario with a non-zero mutated_fraction.

    Each replicate has its own set of mutation TSVs (mutations are drawn
    independently per replicate), so the replicate wildcard is included in
    the path.
    """
    contribs = config["scenarios"][wildcards.scenario]
    seen     = {}
    for c in contribs:
        ref = c["ref_id"]
        if c["mutated_fraction"] > 0 and ref not in seen:
            seen[ref] = (
                f"results/mutated/{ref}/{wildcards.replicate}"
                f"/{ref}.mutations.tsv"
            )
    return list(seen.values())


# ---------------------------------------------------------------------------
# Stage 1 — assess one (scenario × replicate)
# ---------------------------------------------------------------------------

rule assess_variants:
    """
    Compare a single breseq VCF against the ground-truth mutation tables for
    the same replicate.

    For each expected mutation (from the ground-truth TSV), the rule checks
    whether breseq called a variant at that position and records:
      - whether it was detected at all
      - the called ALT allele
      - the allele frequency reported by breseq
      - the QUAL score
      - whether the QUAL score meets the minimum threshold

    One output row is written per expected mutation.
    """
    input:
        vcf          = "results/breseq/{scenario}/{replicate}/output/output.vcf",
        ground_truth = scenario_ground_truth_tsvs,
    output:
        tsv="results/assessment/{scenario}/{replicate}_assessment.tsv",
    params:
        min_quality = config["assessment"]["min_quality"],
        replicate   = lambda wc: wc.replicate,
    log:
        "logs/assessment/{scenario}/{replicate}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assess_variants.py"


# ---------------------------------------------------------------------------
# Stage 2 — aggregate all replicates for one scenario
# ---------------------------------------------------------------------------

def all_replicate_tsvs(wildcards):
    """Return the list of per-replicate assessment TSVs for this scenario."""
    return expand(
        "results/assessment/{scenario}/{replicate}_assessment.tsv",
        scenario=wildcards.scenario,
        replicate=REPLICATES,
    )


rule aggregate_replicates:
    """
    Concatenate per-replicate assessment TSVs into one table for a scenario.

    A 'replicate' column is added to each row so that downstream analysis can
    distinguish results from different replicates and compute statistics such
    as:
      - mean / variance of detection rate across replicates
      - sensitivity as a function of allele frequency or mutation type
      - consistency of quality scores for the same mutation across replicates

    The header line is written once; subsequent files are appended without
    their header to avoid duplicates.
    """
    input:
        tsvs=all_replicate_tsvs,
    output:
        tsv="results/assessment/{scenario}_all_replicates.tsv",
    log:
        "logs/assessment/{scenario}_aggregate.log",
    run:
        import csv, logging

        logging.basicConfig(
            filename=log[0],
            level=logging.INFO,
            format="%(asctime)s %(levelname)s %(message)s",
        )
        logger = logging.getLogger(__name__)

        header_written = False
        total_rows     = 0

        with open(output.tsv, "w", newline="") as out_fh:
            writer = None

            for tsv_path in sorted(input.tsvs):
                # Derive the replicate label from the filename
                # e.g. "results/assessment/scenario_equal_mix/rep2_assessment.tsv"
                #  → "rep2"
                rep_label = os.path.basename(tsv_path).replace("_assessment.tsv", "")

                with open(tsv_path, newline="") as in_fh:
                    reader = csv.DictReader(in_fh, delimiter="\t")

                    for row in reader:
                        # Attach the replicate label to every row
                        row["replicate"] = rep_label

                        if writer is None:
                            # Initialise writer with the fieldnames from the
                            # first file, plus the new 'replicate' column
                            fieldnames = list(reader.fieldnames) + ["replicate"]
                            writer = csv.DictWriter(
                                out_fh,
                                fieldnames=fieldnames,
                                delimiter="\t",
                            )
                            writer.writeheader()

                        writer.writerow(row)
                        total_rows += 1

                logger.info(f"Appended {rep_label} from {tsv_path}")

        logger.info(
            f"Aggregated {total_rows} rows across {len(input.tsvs)} replicates "
            f"→ {output.tsv}"
        )
