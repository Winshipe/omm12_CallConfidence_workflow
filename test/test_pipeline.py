#!/usr/bin/env python3
"""
test_pipeline.py
================
Entry point for the CallConfidence workflow test suite.

No test framework, Snakemake, conda environment or bioinformatics tool is
needed: the suite uses only the Python standard library.  Tests that exercise
scripts importing pandas / Biopython, or the R report, run when those are
available and are reported as *skipped* otherwise.

How to run (from anywhere)
--------------------------
    python test/test_pipeline.py              # run everything
    python test/test_pipeline.py -v           # one line per test, with timings
    python test/test_pipeline.py -k assess    # only tests whose name contains "assess"
    python test/test_pipeline.py --list       # list tests without running them
    python test/test_pipeline.py --no-r       # skip the (slower) R report test

Exit status is 0 when nothing failed.  Tests marked ``known_bug`` document a
current defect in the workflow: they are reported but do not fail the run, and
they *do* fail once the bug is fixed so the marker can be removed.

The tests also run under pytest if you happen to have it:  pytest test/

Test modules
------------
  test_mutate_reference.py    mutation models, FASTA I/O, full script runs
  test_find_homologs.py       cluster parsing, Prodigal headers, BED output
  test_assess_variants.py     TP/FP/FN classification via the real script
  test_annotate_repeats.py    m8 → annotation TSV post-processing
  test_mmseqs_search.py       search orchestration with a fake `mmseqs`
  test_blend_reads.py         FASTQ subsampling / blending script
  test_rules.py               Python helpers and `run:` blocks inside the .smk files
  test_workflow.py            config validation, rule/file consistency, shell
                              placeholder checks and a dry-run DAG build
  test_integration.py         mutate → (synthetic) VCF → assess → aggregate
  test_report.py              knits genome_report.Rmd on synthetic data (R)

Shared plumbing lives in harness.py (runner, snakemake stub, fixtures) and
smk.py (dependency-free Snakefile/YAML reader and DAG resolver).
"""

import argparse
import importlib
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True          # keep __pycache__ out of the repo
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import harness  # noqa: E402

MODULES = [
    "test_mutate_reference",
    "test_find_homologs",
    "test_assess_variants",
    "test_annotate_repeats",
    "test_mmseqs_search",
    "test_blend_reads",
    "test_rules",
    "test_workflow",
    "test_integration",
    "test_report",
]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("-k", dest="pattern", default="",
                    help="only run tests whose name contains this substring")
    ap.add_argument("--list", action="store_true", help="list tests and exit")
    ap.add_argument("--no-r", action="store_true", help="skip tests that need R")
    args = ap.parse_args(argv)

    if args.no_r:
        os.environ["CALLCONF_SKIP_R"] = "1"

    tests = []
    for name in MODULES:
        tests.extend(harness.collect(importlib.import_module(name)))
    if args.pattern:
        tests = [(n, f) for n, f in tests if args.pattern in n]
    if args.list:
        for n, f in tests:
            flag = "  [known bug]" if getattr(f, "known_bug", None) else ""
            print(n + flag)
        return 0
    return harness.run_tests(tests, verbose=args.verbose)


if __name__ == "__main__":
    sys.exit(main())
