"""
tests/test_integration.py
==========================
End-to-end integration test for the CallConfidence Snakemake pipeline.

What this test does
-------------------
Rather than running on a real multi-megabase genome (which would take hours),
the test builds a minimal but *biologically valid* self-contained workspace:

  1.  Generates a synthetic 10 kbp reference FASTA in-process.
  2.  Creates a matching minimal mmseqs2 database (a FASTA containing one
      "mobile element" protein, pre-indexed with mmseqs createdb).
  3.  Writes a test config.yaml that points at those resources, sets
      coverage=5x and replicates=1 so the run is fast, and disables
      report generation (R/rmarkdown are not required for the core test).
  4.  Runs `snakemake --use-conda --cores 4` up to (but not including) the
      report step, targeting the aggregated assessment TSV.
  5.  Parses that TSV and asserts that:
        a. The file is non-empty and has the expected columns.
        b. Every row has a valid Truthiness value (TP / FP / FN).
        c. At least one mutation was introduced (the mutation engine ran).
        d. The overall recall (TP / (TP + FN)) on a 0.5 mutated-fraction
           scenario is above a minimal sanity floor (>10 %).  At 5× coverage
           and 50 % mutant fraction the effective depth on the mutant allele
           is ~2.5×, so perfect recall is not expected; we just confirm the
           caller produced non-trivial output.

Prerequisites (must already be in PATH or loadable via conda envs)
------------------------------------------------------------------
  snakemake  ≥ 7
  bwa
  samtools
  gatk
  art_illumina  (from the ART package)
  seqtk
  prodigal
  mmseqs

The test uses `--use-conda` so every rule picks up its own conda env.
If your envs are already activated and the tools are in PATH, you can
add --conda-not-use-env-modules or simply make sure the tools are available.

Usage
-----
    # Run from the repository root:
    pytest tests/test_integration.py -v -s

    # Increase timeout if running on a slow machine (default: 30 min):
    pytest tests/test_integration.py -v -s --timeout=7200

Configuration
-------------
Two environment variables let you override defaults without editing this file:

    CALLCONF_TEST_CORES   number of CPU cores for snakemake  (default: 4)
    CALLCONF_TEST_TIMEOUT timeout in seconds                  (default: 1800)
"""

import os
import random
import shutil
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pandas as pd
import pytest
import yaml


# ─────────────────────────────────────────────────────────────────────────────
# Tuneable constants (override with environment variables)
# ─────────────────────────────────────────────────────────────────────────────

CORES   = int(os.environ.get("CALLCONF_TEST_CORES",   4))
TIMEOUT = int(os.environ.get("CALLCONF_TEST_TIMEOUT", 1800))  # 30 min default

# The test sequence length.  Long enough to give ART enough room to simulate
# reads at the configured fragment length (mean 300 bp), but short enough that
# GATK + BWA finish in a reasonable time.
REF_LEN = 10_000

# Substitution rate for the test run.  High enough that we introduce a useful
# number of mutations into a 10 kbp sequence (expected ~50 SNPs).
SUBSTITUTION_RATE = 0.005

# Minimum acceptable recall for the sanity check (see module docstring).
MIN_RECALL_FLOOR = 0.05  # 5 % — deliberately lenient given low coverage


# ─────────────────────────────────────────────────────────────────────────────
# Helper: generate a synthetic reference FASTA
# ─────────────────────────────────────────────────────────────────────────────

def _make_reference_fasta(path: Path, length: int = REF_LEN, seed: int = 1) -> None:
    """
    Write a single-contig FASTA with a random nucleotide sequence.

    We use a fixed seed so the test is deterministic across runs.
    The sequence is named 'test_contig' to give predictable contig IDs in
    the VCF output.
    """
    rng = random.Random(seed)
    bases = "ACGT"
    # Slightly GC-biased (55 %) to mimic bacterial genomes
    weights = [0.225, 0.275, 0.275, 0.225]
    seq = "".join(rng.choices(bases, weights=weights, k=length))

    with open(path, "w") as fh:
        fh.write(">test_contig\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i : i + 70] + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# Helper: build a minimal mmseqs2 database
# ─────────────────────────────────────────────────────────────────────────────

def _make_mmseqs_db(db_dir: Path) -> Path:
    """
    Create a minimal mmseqs2 protein database from a single fake IS-element
    protein sequence and return the database prefix path.

    mmseqs createdb only needs a FASTA; the resulting set of index files is
    what the pipeline's mmseqs_search rule expects.

    Returns the path prefix (without extension) that the pipeline config
    should point at.
    """
    db_dir.mkdir(parents=True, exist_ok=True)
    fasta = db_dir / "te_sequences.faa"

    # A plausible 100-aa transposase-like protein sequence (entirely synthetic)
    fake_protein = (
        "MSKQELRAAAERPFKQRLIQNCLGQVVNHFGHQEALDEATMQELLEELQALNLEKQKME"
        "LAAQARLQGWLQSHRQAELERLKQLQEELKAQRQAELERLKQLQEELKAQR"
    )
    fasta.write_text(f">IS1_transposase#DNA/IS1\n{fake_protein}\n")

    db_prefix = db_dir / "test_te_db"
    result = subprocess.run(
        ["mmseqs", "createdb", str(fasta), str(db_prefix)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            f"mmseqs createdb failed — is mmseqs in PATH?\n{result.stderr}"
        )
    return db_prefix


# ─────────────────────────────────────────────────────────────────────────────
# Helper: write the test config.yaml
# ─────────────────────────────────────────────────────────────────────────────

def _write_config(config_path: Path, ref_fasta: Path, mmseqs_db: Path) -> None:
    """
    Write a config.yaml that is structurally identical to the production one
    but tuned for fast integration testing:

      - 1 replicate instead of 3  (keeps the run time short)
      - 5× coverage instead of 50×  (produces ~3 300 read pairs on a 10 kbp ref)
      - fragment length 300 / sdev 50  (safely above ART's minimum for 150 bp reads)
      - higher substitution rate so we get ~50 mutations to detect
      - report step is omitted from the test targets (R not required)
    """
    cfg = {
        "references": {
            "test_ref": str(ref_fasta),
        },
        "replicates": 1,
        "mutation": {
            "model": "tamura_nei",
            "substitution_rate": SUBSTITUTION_RATE,
            "kappa": 2.0,
            "gc_freq": 0.5,
            "base_seed": 42,
        },
        "annotation": {
            "mmseqs_databases": {
                "test_te_db": str(mmseqs_db),
            },
            "split_mem_limit": "4G",
            "max_memory": "8G",
            "threads": 2,
        },
        "simulation": {
            "platform": "MSv3",
            "read_length": 150,
            "mean_fragment_length": 300,
            "std_fragment_length": 50,
            "coverage": 5,
            "empirical_reads_R1": None,
            "empirical_reads_R2": None,
            "threads": 2,
        },
        "scenarios": {
            "test_scenario": [
                {
                    "ref_id": "test_ref",
                    "mutated_fraction": 0.5,
                    "abundance": 1.0,
                }
            ]
        },
        "variant_calling": {
            "threads": 2,
            "mem_mb": 4000,
            "min_base_quality_score": 20,
            "gatk_extra_flags": "",
            "ploidy": 2,
        },
        "assessment": {
            "min_quality": 20,
        },
        "report": {
            "genome_length": 0,
            "window_size": 1000,
        },
    }
    with open(config_path, "w") as fh:
        yaml.dump(cfg, fh, default_flow_style=False)


# ─────────────────────────────────────────────────────────────────────────────
# Helper: check required tools are available
# ─────────────────────────────────────────────────────────────────────────────

REQUIRED_TOOLS = [
    "snakemake",
    "bwa",
    "samtools",
    "gatk",
    "art_illumina",
    "seqtk",
    "prodigal",
    "mmseqs",
]


def _check_tools() -> list[str]:
    """Return a list of tools that are not found in PATH."""
    missing = []
    for tool in REQUIRED_TOOLS:
        if shutil.which(tool) is None:
            missing.append(tool)
    return missing


# ─────────────────────────────────────────────────────────────────────────────
# Pytest fixture: build the workspace once per test session
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def pipeline_workspace(tmp_path_factory):
    """
    Set up a complete, self-contained pipeline workspace and run Snakemake.

    The workspace is created under pytest's tmp directory and is preserved
    for the duration of the test session so that all test functions can
    inspect the same outputs.

    The fixture runs Snakemake targeting the aggregated assessment TSV — the
    final output of the core pipeline (before the R report step).  If
    Snakemake fails, the fixture calls pytest.fail() with the captured stderr
    so the error is visible in the test output.
    """
    # ── 0. Check tools ────────────────────────────────────────────────────────
    missing = _check_tools()
    if missing:
        pytest.skip(
            f"Integration test skipped — the following tools are not in PATH: "
            f"{', '.join(missing)}\n"
            f"Install them (or activate the appropriate conda environment) and "
            f"re-run."
        )

    # ── 1. Build workspace directory tree ────────────────────────────────────
    root = tmp_path_factory.mktemp("callconfidence_integration")
    (root / "resources").mkdir()
    (root / "config").mkdir()
    (root / "logs").mkdir()

    # Copy the workflow directory tree from the repository.
    # This test file lives at tests/test_integration.py; the workflow lives
    # at workflow/ relative to the repository root.
    repo_root = Path(__file__).parent.parent
    workflow_src = repo_root / "workflow"
    if not workflow_src.exists():
        pytest.fail(
            f"Could not find workflow/ directory at {workflow_src}.\n"
            f"Make sure you are running pytest from the repository root."
        )
    shutil.copytree(workflow_src, root / "workflow")

    # ── 2. Create test inputs ─────────────────────────────────────────────────
    ref_fasta = root / "resources" / "test_ref.fa"
    _make_reference_fasta(ref_fasta)

    mmseqs_db_dir = root / "resources" / "mmseqs"
    mmseqs_db_prefix = _make_mmseqs_db(mmseqs_db_dir)

    # ── 3. Write config ───────────────────────────────────────────────────────
    config_path = root / "config" / "config.yaml"
    _write_config(config_path, ref_fasta, mmseqs_db_prefix)

    # ── 4. Define the Snakemake targets ──────────────────────────────────────
    # We target the aggregated assessment TSV, which is the final output of
    # the core pipeline excluding the R report.  This exercises every rule
    # except generate_report.
    targets = [
        # Ground-truth mutation table
        "results/mutated/test_ref/rep1/test_ref.mutations.tsv",
        # Annotation (challenging regions)
        "results/annotation/test_ref/test_ref.challenging.tsv",
        # Blended reads
        "results/blended/test_scenario/rep1/test_scenario_R1.fastq.gz",
        "results/blended/test_scenario/rep1/test_scenario_R2.fastq.gz",
        # Variant calls
        "results/variant_calling/test_scenario/rep1/output.vcf",
        # Per-replicate assessment
        "results/assessment/test_scenario/rep1_assessment.tsv",
        # Aggregated assessment
        "results/assessment/test_scenario_all_replicates.tsv",
    ]

    # ── 5. Run Snakemake ──────────────────────────────────────────────────────
    cmd = [
        "snakemake",
        "--snakefile",   str(root / "workflow" / "Snakefile"),
        "--configfile",  str(config_path),
        "--use-conda",
        "--conda-prefix", str(root / ".snakemake" / "conda"),
        "--cores",       str(CORES),
        "--directory",   str(root),
        "--rerun-incomplete",
        # Print shell commands to stdout so failures are diagnosable
        "--printshellcmds",
        *targets,
    ]

    print(f"\n[integration] Running Snakemake in {root}")
    print(f"[integration] Command: {' '.join(cmd)}\n")

    start = time.monotonic()
    result = subprocess.run(
        cmd,
        capture_output=False,   # let stdout/stderr flow through to pytest -s
        timeout=TIMEOUT,
        cwd=str(root),
    )
    elapsed = time.monotonic() - start
    print(f"\n[integration] Snakemake finished in {elapsed:.0f}s "
          f"(exit code {result.returncode})")

    if result.returncode != 0:
        pytest.fail(
            f"Snakemake exited with code {result.returncode} after {elapsed:.0f}s.\n"
            f"Check the output above for the failing rule and log file path.\n"
            f"Workspace is preserved at: {root}"
        )

    # Expose the workspace root to all test functions
    return root


# ─────────────────────────────────────────────────────────────────────────────
# Test functions — each asserts one property of the pipeline output
# ─────────────────────────────────────────────────────────────────────────────

class TestPipelineOutputsExist:
    """Every expected output file is present and non-empty."""

    def test_mutations_tsv_exists(self, pipeline_workspace):
        p = pipeline_workspace / "results/mutated/test_ref/rep1/test_ref.mutations.tsv"
        assert p.exists(), "Ground-truth mutations TSV was not created"
        assert p.stat().st_size > 0, "Ground-truth mutations TSV is empty"

    def test_mutated_fasta_exists(self, pipeline_workspace):
        p = pipeline_workspace / "results/mutated/test_ref/rep1/test_ref.mutated.fasta"
        assert p.exists(), "Mutated reference FASTA was not created"
        assert p.stat().st_size > 0, "Mutated reference FASTA is empty"

    def test_blended_reads_exist(self, pipeline_workspace):
        for read in ("R1", "R2"):
            p = (pipeline_workspace /
                 f"results/blended/test_scenario/rep1/test_scenario_{read}.fastq.gz")
            assert p.exists(), f"Blended {read} FASTQ was not created"
            assert p.stat().st_size > 0, f"Blended {read} FASTQ is empty"

    def test_vcf_exists(self, pipeline_workspace):
        p = pipeline_workspace / "results/variant_calling/test_scenario/rep1/output.vcf"
        assert p.exists(), "VCF output was not created"
        assert p.stat().st_size > 0, "VCF output is empty"

    def test_per_replicate_assessment_tsv_exists(self, pipeline_workspace):
        p = pipeline_workspace / "results/assessment/test_scenario/rep1_assessment.tsv"
        assert p.exists(), "Per-replicate assessment TSV was not created"
        assert p.stat().st_size > 0, "Per-replicate assessment TSV is empty"

    def test_aggregated_assessment_tsv_exists(self, pipeline_workspace):
        p = pipeline_workspace / "results/assessment/test_scenario_all_replicates.tsv"
        assert p.exists(), "Aggregated assessment TSV was not created"
        assert p.stat().st_size > 0, "Aggregated assessment TSV is empty"

    def test_challenging_tsv_exists(self, pipeline_workspace):
        p = pipeline_workspace / "results/annotation/test_ref/test_ref.challenging.tsv"
        assert p.exists(), "Challenging-regions TSV was not created"

    def test_bam_exists(self, pipeline_workspace):
        p = pipeline_workspace / "results/variant_calling/test_scenario/rep1/aligned.bam"
        assert p.exists(), "Sorted BAM was not created"
        assert p.stat().st_size > 0, "Sorted BAM is empty"

    def test_bam_index_exists(self, pipeline_workspace):
        p = pipeline_workspace / "results/variant_calling/test_scenario/rep1/aligned.bam.bai"
        assert p.exists(), "BAM index (.bai) was not created"


class TestMutationsFile:
    """The ground-truth mutations TSV has the right structure and content."""

    @pytest.fixture(autouse=True)
    def load(self, pipeline_workspace):
        path = pipeline_workspace / "results/mutated/test_ref/rep1/test_ref.mutations.tsv"
        self.df = pd.read_csv(path, sep="\t")

    def test_expected_columns_present(self):
        # The mutations TSV written by mutate_reference.py has these columns
        # (mutation_type is currently commented out in the script, so we
        # don't require it here).
        for col in ("seq_id", "position", "ref_base", "alt_base"):
            assert col in self.df.columns, f"Column '{col}' missing from mutations TSV"

    def test_at_least_one_mutation_introduced(self):
        """
        At a substitution rate of 0.5 % over 10 kbp we expect ~50 mutations.
        Even with stochasticity the probability of getting zero is negligible.
        """
        assert len(self.df) > 0, \
            "No mutations were introduced — mutate_reference.py may be broken"

    def test_positions_are_one_based(self):
        assert self.df["position"].min() >= 1, \
            "Positions should be 1-based (minimum position must be ≥ 1)"

    def test_positions_within_reference_length(self):
        assert self.df["position"].max() <= REF_LEN, \
            f"A mutation position exceeds the reference length ({REF_LEN})"

    def test_ref_and_alt_are_valid_bases(self):
        valid = {"A", "C", "G", "T"}
        bad_ref = set(self.df["ref_base"].unique()) - valid
        bad_alt = set(self.df["alt_base"].unique()) - valid
        assert not bad_ref, f"Invalid ref_base values: {bad_ref}"
        assert not bad_alt, f"Invalid alt_base values: {bad_alt}"

    def test_ref_differs_from_alt(self):
        same = self.df[self.df["ref_base"] == self.df["alt_base"]]
        assert len(same) == 0, \
            f"{len(same)} rows have ref_base == alt_base (not a real mutation)"

    def test_contig_name_matches_reference(self):
        assert (self.df["seq_id"] == "test_contig").all(), \
            "seq_id in mutations TSV does not match the contig name in the reference FASTA"


class TestVcfFile:
    """The VCF produced by GATK is a valid VCF with at least a header."""

    @pytest.fixture(autouse=True)
    def load(self, pipeline_workspace):
        self.vcf_path = (
            pipeline_workspace /
            "results/variant_calling/test_scenario/rep1/output.vcf"
        )
        self.lines = self.vcf_path.read_text().splitlines()

    def test_vcf_has_header(self):
        meta_lines = [l for l in self.lines if l.startswith("##")]
        assert len(meta_lines) > 0, "VCF file has no ## meta-information header lines"

    def test_vcf_has_chrom_header(self):
        col_header = [l for l in self.lines if l.startswith("#CHROM")]
        assert len(col_header) == 1, "VCF file is missing the #CHROM column header line"

    def test_vcf_data_lines_are_tab_separated(self):
        data_lines = [l for l in self.lines if not l.startswith("#") and l.strip()]
        if not data_lines:
            pytest.skip("VCF contains no variant calls — cannot validate format")
        for line in data_lines[:5]:  # check the first few data lines
            fields = line.split("\t")
            assert len(fields) >= 8, \
                f"VCF data line has fewer than 8 tab-separated fields:\n{line}"

    def test_vcf_chrom_matches_reference(self):
        data_lines = [l for l in self.lines if not l.startswith("#") and l.strip()]
        if not data_lines:
            pytest.skip("VCF contains no variant calls")
        chroms = {l.split("\t")[0] for l in data_lines}
        assert chroms == {"test_contig"}, \
            f"Unexpected CHROM values in VCF: {chroms}"

    def test_vcf_positions_are_positive_integers(self):
        data_lines = [l for l in self.lines if not l.startswith("#") and l.strip()]
        if not data_lines:
            pytest.skip("VCF contains no variant calls")
        for line in data_lines:
            pos = int(line.split("\t")[1])
            assert 1 <= pos <= REF_LEN, \
                f"VCF position {pos} is outside the expected range [1, {REF_LEN}]"


class TestAssessmentFile:
    """The aggregated assessment TSV has correct structure and sensible values."""

    @pytest.fixture(autouse=True)
    def load(self, pipeline_workspace):
        path = (pipeline_workspace /
                "results/assessment/test_scenario_all_replicates.tsv")
        self.df = pd.read_csv(path, sep="\t")

    def test_expected_columns_present(self):
        required = {"CHROM", "POS", "REF", "ALT", "Truthiness", "replicate"}
        missing = required - set(self.df.columns)
        assert not missing, f"Assessment TSV is missing columns: {missing}"

    def test_truthiness_values_are_valid(self):
        valid = {"TP", "FP", "FN"}
        observed = set(self.df["Truthiness"].unique())
        invalid = observed - valid
        assert not invalid, \
            f"Unexpected Truthiness values in assessment: {invalid}"

    def test_replicate_column_is_populated(self):
        assert self.df["replicate"].notna().all(), \
            "Some rows are missing a replicate label"
        assert (self.df["replicate"] == "rep1").all(), \
            "Unexpected replicate label (expected 'rep1')"

    def test_at_least_one_expected_mutation_in_assessment(self):
        """
        The assessment must contain at least one ground-truth mutation entry
        (TP or FN).  A table with only FP rows would mean the ground-truth
        join is broken.
        """
        ground_truth_rows = self.df[self.df["Truthiness"].isin({"TP", "FN"})]
        assert len(ground_truth_rows) > 0, \
            "Assessment contains no TP or FN rows — the ground-truth join may be broken"

    def test_no_duplicate_ground_truth_positions(self):
        """
        Each (CHROM, POS, REF, ALT) tuple from the ground truth should appear
        exactly once in the assessment.  Duplicates would indicate a bug in
        the merge logic.
        """
        gt_rows = self.df[self.df["Truthiness"].isin({"TP", "FN"})]
        dupes = gt_rows.duplicated(subset=["CHROM", "POS", "REF", "ALT"])
        assert not dupes.any(), \
            f"{dupes.sum()} duplicate ground-truth positions in assessment"

    def test_recall_above_floor(self):
        """
        Sanity-check that the variant caller is actually finding *some*
        mutations.  At 5× coverage and a 0.5 mutated fraction, the effective
        allele depth is ~2.5×, so recall will be low — but it should not be
        zero.

        MIN_RECALL_FLOOR is set to 5 % to catch complete failures while
        not requiring high performance from an intentionally low-coverage run.
        """
        tp = (self.df["Truthiness"] == "TP").sum()
        fn = (self.df["Truthiness"] == "FN").sum()
        denominator = tp + fn
        if denominator == 0:
            pytest.fail("No TP or FN rows — cannot compute recall")
        recall = tp / denominator
        print(f"\n[integration] Recall: {tp}/{denominator} = {recall:.1%}")
        assert recall >= MIN_RECALL_FLOOR, (
            f"Recall {recall:.1%} is below the sanity floor {MIN_RECALL_FLOOR:.0%}. "
            f"The variant caller may not be running correctly."
        )

    def test_positions_match_between_mutations_and_assessment(self, pipeline_workspace):
        """
        Every position recorded in the ground-truth mutations TSV should
        appear in the assessment as either TP or FN.  A missing position
        means the join in assess_variants.py dropped a mutation silently.
        """
        mut_path = (pipeline_workspace /
                    "results/mutated/test_ref/rep1/test_ref.mutations.tsv")
        mutations = pd.read_csv(mut_path, sep="\t")
        gt_positions = set(zip(mutations["seq_id"], mutations["position"]))

        assessment_gt = self.df[self.df["Truthiness"].isin({"TP", "FN"})]
        assessed_positions = set(zip(assessment_gt["CHROM"], assessment_gt["POS"]))

        missing = gt_positions - assessed_positions
        assert not missing, (
            f"{len(missing)} mutation(s) from the ground-truth TSV are absent "
            f"from the assessment (showing up to 5):\n"
            + "\n".join(str(p) for p in list(missing)[:5])
        )


class TestReproducibility:
    """
    Re-running Snakemake on an already-complete workspace should be a no-op
    (all outputs up to date, no jobs re-run).
    """

    def test_snakemake_reports_nothing_to_do(self, pipeline_workspace):
        """
        A second dry-run should report 0 jobs to execute.
        This confirms that all outputs are correctly declared and that
        Snakemake's dependency tracking is consistent.
        """
        repo_root = Path(__file__).parent.parent
        cmd = [
            "snakemake",
            "--snakefile", str(pipeline_workspace / "workflow" / "Snakefile"),
            "--configfile", str(pipeline_workspace / "config" / "config.yaml"),
            "--cores", str(CORES),
            "--directory", str(pipeline_workspace),
            "--dry-run",
            "--quiet",
            "results/assessment/test_scenario_all_replicates.tsv",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True,
                                cwd=str(pipeline_workspace))

        assert result.returncode == 0, (
            f"Snakemake dry-run failed on an already-complete workspace.\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
        # Snakemake prints "Nothing to be done." when all targets are current.
        assert "Nothing to be done" in result.stdout or "nothing to be done" in result.stdout.lower(), (
            f"Snakemake dry-run suggests jobs need re-running after a complete "
            f"pipeline execution.  This may indicate missing output declarations "
            f"in one of the .smk rules.\nstdout:\n{result.stdout}"
        )


class TestLogFiles:
    """Log files are created alongside every rule's output."""

    def _check_log(self, pipeline_workspace, rel_path: str) -> None:
        p = pipeline_workspace / rel_path
        assert p.exists(), f"Log file not found: {rel_path}"
        assert p.stat().st_size > 0, f"Log file is empty: {rel_path}"

    def test_mutate_log_exists(self, pipeline_workspace):
        self._check_log(pipeline_workspace, "logs/mutate/test_ref/rep1.log")

    def test_simulation_log_exists(self, pipeline_workspace):
        # ART writes one log per (ref_id, replicate, mutated) combination
        self._check_log(
            pipeline_workspace,
            "logs/simulation/test_ref/rep1_mutated.log"
        )

    def test_variant_calling_log_exists(self, pipeline_workspace):
        self._check_log(
            pipeline_workspace,
            "logs/variant_calling/test_scenario/rep1.gatk.log"
        )

    def test_assessment_log_exists(self, pipeline_workspace):
        self._check_log(
            pipeline_workspace,
            "logs/assessment/test_scenario/rep1.log"
        )
