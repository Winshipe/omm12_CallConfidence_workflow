#!/usr/bin/env python3
"""
scripts/mutate_reference.py
───────────────────────────
Apply a substitution model to a reference FASTA and write:
  - a mutated FASTA
  - a ground-truth TSV:  seq_id  position  ref_base  alt_base  mutation_type

Snakemake injects:
  snakemake.input.ref
  snakemake.output.mutated_fasta
  snakemake.output.mutations_tsv
  snakemake.params.{model, rate, kappa, gc_freq, seed}
  snakemake.log[0]
"""

import sys
import logging
import random
import math
from pathlib import Path

# Route all output to the Snakemake log file
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


# ── FASTA I/O ────────────────────────────────────────────────────────────────

def read_fasta(path):
    """Return list of (header, sequence) tuples."""
    records = []
    header, seq_parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def write_fasta(path, records, line_width=70):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i : i + line_width] + "\n")


# ── Substitution models ───────────────────────────────────────────────────────

BASES = ("A", "C", "G", "T")
PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T"}


def jukes_cantor_probs(ref_base):
    """
    Equal rate substitution to any of the three alternative bases.
    Returns dict {base: prob} for the *alternative* bases only.
    """
    others = [b for b in BASES if b != ref_base]
    return {b: 1 / 3 for b in others}


def tamura_nei_probs(ref_base, kappa, pi):
    """
    Tamura & Nei (1993) substitution probabilities for *alternative* bases.
    pi = {"A": freq_A, "C": freq_C, "G": freq_G, "T": freq_T}
    kappa = transition / transversion rate ratio
    """
    probs = {}
    for alt in BASES:
        if alt == ref_base:
            continue
        if (ref_base in PURINES) == (alt in PURINES):
            # transition
            weight = kappa * pi[alt]
        else:
            # transversion
            weight = pi[alt]
        probs[alt] = weight

    # normalise
    total = sum(probs.values())
    return {b: v / total for b, v in probs.items()}


def mutate_sequence(seq, model, rate, kappa, gc_freq, rng):
    """
    Walk every position and stochastically introduce substitutions.
    Returns (mutated_seq_str, list_of_mutation_dicts).
    """
    at_freq = (1.0 - gc_freq) / 2.0
    gc_per  = gc_freq / 2.0
    pi = {"A": at_freq, "T": at_freq, "G": gc_per, "C": gc_per}

    seq = list(seq)
    mutations = []

    for pos, ref_base in enumerate(seq):
        if ref_base not in BASES:        # skip ambiguous / gap characters
            continue
        if rng.random() > rate:
            continue                     # no mutation at this site

        if model == "jukes_cantor":
            alt_probs = jukes_cantor_probs(ref_base)
        elif model == "tamura_nei":
            alt_probs = tamura_nei_probs(ref_base, kappa, pi)
        else:
            raise ValueError(f"Unknown model: {model}")

        alts, weights = zip(*alt_probs.items())
        alt_base = rng.choices(alts, weights=weights, k=1)[0]

        # classify
        if (ref_base in PURINES) == (alt_base in PURINES):
            mut_type = "transition"
        else:
            mut_type = "transversion"

        seq[pos] = alt_base
        mutations.append(
            {
                "position": pos + 1,   # 1-based
                "ref_base": ref_base,
                "alt_base": alt_base,
                #"mutation_type": mut_type,
            }
        )

    return "".join(seq), mutations


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    rng = random.Random(snakemake.params.seed)

    model    = snakemake.params.model
    rate     = float(snakemake.params.rate)
    kappa    = float(snakemake.params.kappa)
    gc_freq  = float(snakemake.params.gc_freq)

    log.info(f"Model: {model}  rate: {rate}  kappa: {kappa}  gc_freq: {gc_freq}")

    records = read_fasta(snakemake.input.ref)
    log.info(f"Read {len(records)} sequence(s) from {snakemake.input.ref}")

    mutated_records = []
    all_mutations = []

    for header, seq in records:
        seq_id = header.split()[0]
        mutated_seq, muts = mutate_sequence(seq, model, rate, kappa, gc_freq, rng)
        mutated_records.append((f"{header} [mutated]", mutated_seq))
        for m in muts:
            m["seq_id"] = seq_id
        all_mutations.extend(muts)
        log.info(f"  {seq_id}: {len(muts)} mutations introduced out of {len(seq)} bp")

    write_fasta(snakemake.output.mutated_fasta, mutated_records)
    log.info(f"Wrote mutated FASTA to {snakemake.output.mutated_fasta}")

    # Write TSV
    with open(snakemake.output.mutations_tsv, "w") as fh:
        fh.write("seq_id\tposition\tref_base\talt_base\tmutation_type\n")
        for m in all_mutations:
            fh.write(
                f"{m['seq_id']}\t{m['position']}\t{m['ref_base']}\t"
                f"{m['alt_base']}\n" #\t{m['mutation_type']}\n"
            )
    log.info(f"Wrote ground-truth TSV to {snakemake.output.mutations_tsv}")
    log.info("Done.")


main()
