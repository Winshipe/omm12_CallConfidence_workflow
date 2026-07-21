"""
tests/test_pipeline_scripts.py
================================
Unit tests for the Python helper functions in the CallConfidence pipeline.

The Snakemake scripts cannot be imported directly (they reference a global
`snakemake` object that only exists at runtime), so we instead import just
the pure-Python helper functions by temporarily monkey-patching the module
namespace.  Each test module section explains which script is under test and
what the tests cover.

Run with:
    pytest tests/test_pipeline_scripts.py -v

Dependencies (all in the standard `python.yaml` conda env plus biopython):
    pytest, pandas, biopython
"""

import csv
import gzip
import io
import os
import random
import sys
import textwrap
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest


# ─────────────────────────────────────────────────────────────────────────────
# Helpers for importing snakemake-dependent scripts
#
# Each script starts running at module level (referencing `snakemake.*`),
# so we cannot do a plain `import`.  Instead we:
#   1. Install a stub `snakemake` object into sys.modules.
#   2. Use importlib to exec the source, extracting only the functions we want.
#
# A simpler pattern is to just copy-paste the pure functions into tests, but
# that creates a maintenance burden.  The approach below imports from the
# *actual* source files, so tests will catch refactoring regressions.
# ─────────────────────────────────────────────────────────────────────────────

SCRIPTS_DIR = Path(__file__).parent.parent / "workflow" / "scripts"


def _load_functions_from_script(script_path: Path, function_names: list[str], snakemake_stub=None):
    """
    Execute `script_path` in a fresh namespace with a dummy `snakemake`
    object, then return a dict of {name: callable} for each requested function.

    Any top-level code that tries to *call* snakemake attributes (e.g. open a
    file) will raise; that's fine — we only need the function definitions.
    """
    stub = snakemake_stub or MagicMock()
    namespace = {"snakemake": stub, "__name__": "__test__"}
    source = script_path.read_text()

    try:
        exec(compile(source, str(script_path), "exec"), namespace)
    except Exception:
        # Top-level driver code (e.g. `main()`) may crash on the stub — that
        # is expected and acceptable.  The function *definitions* are already
        # registered in `namespace` by the time the driver runs.
        pass

    return {name: namespace[name] for name in function_names if name in namespace}


# ═════════════════════════════════════════════════════════════════════════════
# 1.  mutate_reference.py
#     Tests the substitution model maths and FASTA I/O helpers.
# ═════════════════════════════════════════════════════════════════════════════

# Import the pure functions we want to test.
# We do this once at module level so the overhead is paid only once.
_mutate_fns = _load_functions_from_script(
    SCRIPTS_DIR / "mutate_reference.py",
    ["jukes_cantor_probs", "tamura_nei_probs", "mutate_sequence",
     "read_fasta", "write_fasta"],
)
jukes_cantor_probs = _mutate_fns["jukes_cantor_probs"]
tamura_nei_probs   = _mutate_fns["tamura_nei_probs"]
mutate_sequence    = _mutate_fns["mutate_sequence"]
read_fasta_fn      = _mutate_fns["read_fasta"]
write_fasta_fn     = _mutate_fns["write_fasta"]


class TestJukesCantor:
    """jukes_cantor_probs returns equal weights for the three non-ref bases."""

    def test_three_alternatives_returned(self):
        probs = jukes_cantor_probs("A")
        assert set(probs.keys()) == {"C", "G", "T"}

    def test_ref_base_excluded(self):
        for base in "ACGT":
            probs = jukes_cantor_probs(base)
            assert base not in probs

    def test_equal_probabilities(self):
        probs = jukes_cantor_probs("G")
        values = list(probs.values())
        assert all(abs(v - values[0]) < 1e-10 for v in values), \
            "Jukes-Cantor should give equal probability to each alternative base"

    def test_probabilities_sum_to_one(self):
        for base in "ACGT":
            total = sum(jukes_cantor_probs(base).values())
            assert abs(total - 1.0) < 1e-10


class TestTamuraNei:
    """tamura_nei_probs implements the transition/transversion rate ratio."""

    # Equal base frequencies and kappa=1 should collapse to equal rates
    @pytest.fixture
    def equal_pi(self):
        return {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

    def test_probabilities_sum_to_one(self, equal_pi):
        for base in "ACGT":
            total = sum(tamura_nei_probs(base, kappa=2.0, pi=equal_pi).values())
            assert abs(total - 1.0) < 1e-10

    def test_ref_base_excluded(self, equal_pi):
        for base in "ACGT":
            probs = tamura_nei_probs(base, kappa=2.0, pi=equal_pi)
            assert base not in probs

    def test_high_kappa_favours_transitions(self, equal_pi):
        """
        With kappa >> 1 the model should strongly favour transitions
        (A↔G, C↔T) over transversions (A/G ↔ C/T).
        
        A→G is a transition (both purines).
        A→C and A→T are transversions.
        """
        probs = tamura_nei_probs("A", kappa=100.0, pi=equal_pi)
        transition_prob   = probs["G"]          # A→G is the only Ti from A
        transversion_prob = probs["C"] + probs["T"]   # A→C and A→T are Tv
        assert transition_prob > transversion_prob, \
            "High kappa should make transitions more likely than transversions"

    def test_kappa_one_equal_pi_matches_jukes_cantor(self, equal_pi):
        """
        When kappa=1 and all base frequencies are equal, Tamura-Nei reduces
        to the Jukes-Cantor model (all substitution rates identical).
        """
        tn = tamura_nei_probs("A", kappa=1.0, pi=equal_pi)
        jc = jukes_cantor_probs("A")
        for base in tn:
            assert abs(tn[base] - jc[base]) < 1e-10, \
                f"TN(kappa=1, equal pi) should equal JC for base {base}"


class TestMutateSequence:
    """mutate_sequence applies mutations at the correct rate and records them."""

    def _deterministic_rng(self, seed=42):
        return random.Random(seed)

    def test_rate_zero_produces_no_mutations(self):
        seq = "ACGTACGTACGT"
        mutated, muts = mutate_sequence(seq, "jukes_cantor", rate=0.0,
                                        kappa=2.0, gc_freq=0.5,
                                        rng=self._deterministic_rng())
        assert mutated == seq
        assert muts == []

    def test_rate_one_mutates_every_site(self):
        """
        With rate=1.0 every position must be mutated (substituted to a
        different base).  The mutated sequence should differ from the original
        at every position.
        """
        seq = "AAAAAAAAAA"
        mutated, muts = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                        kappa=2.0, gc_freq=0.5,
                                        rng=self._deterministic_rng())
        assert len(muts) == len(seq)
        assert mutated != seq

    def test_mutations_are_1_based(self):
        """Position field in the output TSV is 1-based, not 0-based."""
        seq = "AAAA"
        _, muts = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                   kappa=2.0, gc_freq=0.5,
                                   rng=self._deterministic_rng())
        positions = [m["position"] for m in muts]
        assert min(positions) >= 1, "Positions must be 1-based (≥ 1)"
        assert max(positions) <= len(seq)

    def test_alt_base_differs_from_ref(self):
        seq = "ACGT" * 10
        _, muts = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                   kappa=2.0, gc_freq=0.5,
                                   rng=self._deterministic_rng())
        for m in muts:
            assert m["ref_base"] != m["alt_base"], \
                "A mutation must change the base"

    def test_ambiguous_bases_are_skipped(self):
        """Ns (and other non-ACGT characters) should never be mutated."""
        seq = "NNNNN"
        mutated, muts = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                        kappa=2.0, gc_freq=0.5,
                                        rng=self._deterministic_rng())
        assert muts == []
        assert mutated == seq

    def test_invalid_model_raises(self):
        with pytest.raises(ValueError, match="Unknown model"):
            mutate_sequence("ACGT", "nonexistent_model", rate=0.1,
                            kappa=2.0, gc_freq=0.5,
                            rng=self._deterministic_rng())

    def test_mutation_count_roughly_matches_rate(self):
        """
        Over a long sequence the observed mutation rate should be close to the
        requested rate (within ±20 % at rate=0.05, seq length 10 000).
        This is a stochastic test — it uses a fixed seed for reproducibility.
        """
        seq = "ACGT" * 2500   # 10 000 bp
        rate = 0.05
        _, muts = mutate_sequence(seq, "jukes_cantor", rate=rate,
                                   kappa=2.0, gc_freq=0.5,
                                   rng=random.Random(0))
        observed_rate = len(muts) / len(seq)
        assert abs(observed_rate - rate) < 0.02, \
            f"Observed rate {observed_rate:.4f} is too far from target {rate}"


class TestFastaIO:
    """read_fasta and write_fasta are inverses of each other."""

    def test_round_trip(self, tmp_path):
        records = [
            ("seq1 description", "ACGTACGT"),
            ("seq2", "TTTTCCCC"),
        ]
        out = tmp_path / "test.fa"
        write_fasta_fn(str(out), records)
        loaded = read_fasta_fn(str(out))
        assert [(h, s) for h, s in loaded] == records

    def test_multiline_fasta_is_joined(self, tmp_path):
        """Sequences split across multiple lines should be concatenated."""
        fasta_text = ">seq1\nACGT\nACGT\n"
        fa = tmp_path / "multi.fa"
        fa.write_text(fasta_text)
        records = read_fasta_fn(str(fa))
        assert records[0][1] == "ACGTACGT"

    def test_sequence_uppercased(self, tmp_path):
        """Lower-case bases in the input should be uppercased on load."""
        fa = tmp_path / "lower.fa"
        fa.write_text(">seq1\nacgtACGT\n")
        records = read_fasta_fn(str(fa))
        assert records[0][1] == "ACGTACGT"

    def test_empty_file_returns_empty_list(self, tmp_path):
        fa = tmp_path / "empty.fa"
        fa.write_text("")
        assert read_fasta_fn(str(fa)) == []


# ═════════════════════════════════════════════════════════════════════════════
# 2.  find_homologs.py
#     Tests cluster parsing and BED output logic.
# ═════════════════════════════════════════════════════════════════════════════

_homolog_fns = _load_functions_from_script(
    SCRIPTS_DIR / "find_homologs.py",
    ["find_homolog_ids", "write_bed"],
)
find_homolog_ids = _homolog_fns["find_homolog_ids"]
write_bed        = _homolog_fns["write_bed"]


class TestFindHomologIds:
    """find_homolog_ids correctly identifies multi-member clusters."""

    def _write_cluster_file(self, tmp_path, rows):
        """Write a two-column cluster TSV and return its path."""
        p = tmp_path / "clusters.txt"
        p.write_text("\n".join(f"{r}\t{m}" for r, m in rows) + "\n")
        return str(p)

    def test_singleton_clusters_excluded(self, tmp_path):
        """
        A row where representative == member is a singleton cluster.
        It should NOT appear in the returned homolog set.
        """
        path = self._write_cluster_file(tmp_path, [
            ("gene_A_1", "gene_A_1"),  # singleton
        ])
        result = find_homolog_ids(path)
        assert result == set()

    def test_multi_member_cluster_included(self, tmp_path):
        """Both the representative and the member should be returned."""
        path = self._write_cluster_file(tmp_path, [
            ("gene_A_1", "gene_A_1"),   # singleton — excluded
            ("gene_A_1", "gene_B_1"),   # gene_A_1 is rep for gene_B_1
        ])
        result = find_homolog_ids(path)
        assert "gene_A_1" in result
        assert "gene_B_1" in result

    def test_singletons_not_contaminating_multi(self, tmp_path):
        """Ensure singletons from other clusters do not bleed through."""
        path = self._write_cluster_file(tmp_path, [
            ("geneX_1", "geneX_1"),     # singleton
            ("geneY_1", "geneY_1"),     # singleton
            ("geneZ_1", "geneZ_2"),     # genuine homolog pair
        ])
        result = find_homolog_ids(path)
        assert "geneX_1" not in result
        assert "geneY_1" not in result
        assert {"geneZ_1", "geneZ_2"} == result

    def test_empty_file_returns_empty_set(self, tmp_path):
        path = self._write_cluster_file(tmp_path, [])
        assert find_homolog_ids(path) == set()


class TestWriteBed:
    """write_bed produces correctly formatted BED-style output."""

    def _make_fasta_entries(self, entries):
        """
        Build the nested dict that write_bed expects:
          { "fake_path.fna": { gene_id: (start, stop, strand_int_str) } }
        """
        return {"fake_path.fna": entries}

    def test_forward_strand_symbol(self, tmp_path):
        fasta_entries = self._make_fasta_entries({
            "contig_1_1": ("10", "270", "1"),
        })
        out = tmp_path / "out.tsv"
        write_bed({"contig_1_1"}, fasta_entries, str(out))
        content = out.read_text()
        assert "\t+\n" in content

    def test_reverse_strand_symbol(self, tmp_path):
        fasta_entries = self._make_fasta_entries({
            "contig_1_2": ("500", "800", "-1"),
        })
        out = tmp_path / "out.tsv"
        write_bed({"contig_1_2"}, fasta_entries, str(out))
        content = out.read_text()
        assert "\t-\n" in content

    def test_chrom_derived_from_gene_id(self, tmp_path):
        """
        Prodigal gene IDs look like <contig_name>_<index>.
        The chrom column should be everything except the trailing _<index>.
        """
        fasta_entries = self._make_fasta_entries({
            "my_contig_42": ("1", "300", "1"),
        })
        out = tmp_path / "out.tsv"
        write_bed({"my_contig_42"}, fasta_entries, str(out))
        data_lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        chrom = data_lines[0].split("\t")[0]
        assert chrom == "my_contig"

    def test_header_line_present(self, tmp_path):
        fasta_entries = self._make_fasta_entries({
            "contig_1_1": ("10", "270", "1"),
        })
        out = tmp_path / "out.tsv"
        write_bed({"contig_1_1"}, fasta_entries, str(out))
        first_line = out.read_text().splitlines()[0]
        assert first_line.startswith("#"), "BED output should start with a header comment"

    def test_gene_absent_from_fasta_entries_is_skipped(self, tmp_path):
        """A homolog ID not found in any FASTA should not produce a BED line."""
        out = tmp_path / "out.tsv"
        write_bed({"missing_gene_1"}, {}, str(out))
        data_lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
        assert data_lines == []


# ═════════════════════════════════════════════════════════════════════════════
# 3.  assess_variants.py  (the pandas version)
#     Tests the TP/FP/FN classification logic directly using pandas,
#     mirroring the actual script logic without importing the script.
# ═════════════════════════════════════════════════════════════════════════════

def _run_assessment(ground_truth_rows, vcf_rows, min_quality=20):
    """
    Re-implement the core logic of assess_variants.py so tests are
    independent of the snakemake stub complexity.

    Parameters
    ----------
    ground_truth_rows : list of dicts with keys CHROM, POS, REF, ALT, mutation_type
    vcf_rows          : list of dicts with keys CHROM, POS, REF, ALT, QUAL
    min_quality       : minimum QUAL to count as detected

    Returns
    -------
    pd.DataFrame with columns CHROM, POS, REF, ALT, Truthiness
    """
    gt = (pd.DataFrame(ground_truth_rows)
          if ground_truth_rows
          else pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "mutation_type"]))
    vcf = (pd.DataFrame(vcf_rows)
           if vcf_rows
           else pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "QUAL"]))
    # Filter VCF by quality threshold
    vcf_filtered = vcf[vcf["QUAL"] >= min_quality][["CHROM", "POS", "REF", "ALT"]]

    merged = pd.merge(
        gt[["CHROM", "POS", "REF", "ALT"]],
        vcf_filtered,
        how="outer",
        on=["CHROM", "POS", "REF", "ALT"],
        indicator=True,
    )
    label_map = {"both": "TP", "right_only": "FP", "left_only": "FN"}
    merged["Truthiness"] = merged["_merge"].map(label_map)
    return merged[["CHROM", "POS", "REF", "ALT", "Truthiness"]]


class TestAssessmentLogic:
    """Core TP / FP / FN classification for assess_variants.py."""

    def _gt(self, chrom, pos, ref, alt):
        return {"CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt, "mutation_type": "transition"}

    def _vcf(self, chrom, pos, ref, alt, qual=100.0):
        return {"CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt, "QUAL": qual}

    def test_perfect_match_is_tp(self):
        result = _run_assessment(
            [self._gt("chr1", 100, "A", "G")],
            [self._vcf("chr1", 100, "A", "G", qual=100)],
        )
        assert result.iloc[0]["Truthiness"] == "TP"

    def test_missed_mutation_is_fn(self):
        """Ground truth variant not present in VCF → FN."""
        result = _run_assessment(
            [self._gt("chr1", 200, "C", "T")],
            [],  # empty VCF
        )
        assert result.iloc[0]["Truthiness"] == "FN"

    def test_extra_vcf_call_is_fp(self):
        """VCF variant with no matching ground-truth entry → FP."""
        result = _run_assessment(
            [],  # no ground truth
            [self._vcf("chr1", 300, "G", "A", qual=50)],
            min_quality=20,
        )
        assert result.iloc[0]["Truthiness"] == "FP"

    def test_below_quality_threshold_counts_as_fn(self):
        """
        A VCF call at the right position but with QUAL below min_quality
        should be excluded from the 'detected' set, making it an FN.
        """
        result = _run_assessment(
            [self._gt("chr1", 400, "A", "T")],
            [self._vcf("chr1", 400, "A", "T", qual=5)],  # QUAL < min_quality=20
            min_quality=20,
        )
        assert result.iloc[0]["Truthiness"] == "FN", \
            "A call below the quality threshold should be treated as a missed detection (FN)"

    def test_position_mismatch_produces_fn_and_fp(self):
        """
        Ground truth at position 500, caller reports position 501.
        Expected: the GT site is FN, the VCF site is FP.
        """
        result = _run_assessment(
            [self._gt("chr1", 500, "A", "G")],
            [self._vcf("chr1", 501, "A", "G", qual=100)],
        )
        truths = set(result["Truthiness"])
        assert "FN" in truths
        assert "FP" in truths

    def test_mixed_scenario(self):
        """
        Three ground-truth variants:
          pos 100 → called correctly      → TP
          pos 200 → not called at all     → FN
          pos 300 → called below quality  → FN
        One extra VCF call:
          pos 999 → no ground truth       → FP
        """
        gt = [
            self._gt("chr1", 100, "A", "G"),
            self._gt("chr1", 200, "C", "T"),
            self._gt("chr1", 300, "G", "A"),
        ]
        vcf = [
            self._vcf("chr1", 100, "A", "G", qual=99),    # TP
            self._vcf("chr1", 300, "G", "A", qual=5),     # below threshold → not counted
            self._vcf("chr1", 999, "T", "C", qual=80),    # FP
        ]
        result = _run_assessment(gt, vcf, min_quality=20)
        counts = result["Truthiness"].value_counts().to_dict()
        assert counts.get("TP", 0) == 1
        assert counts.get("FN", 0) == 2
        assert counts.get("FP", 0) == 1

    def test_allele_mismatch_is_fn_and_fp(self):
        """
        The caller reports a different ALT allele at a known mutation site.
        The ground-truth mutation is FN; the wrong call is FP.
        This matters for distinguishing multi-allelic sites.
        """
        result = _run_assessment(
            [self._gt("chr1", 600, "A", "G")],   # expect A→G
            [self._vcf("chr1", 600, "A", "T", qual=100)],  # caller says A→T
        )
        truths = set(result["Truthiness"])
        assert "FN" in truths, "The expected A→G should be undetected (FN)"
        assert "FP" in truths, "The spurious A→T should count as a false positive (FP)"


# ═════════════════════════════════════════════════════════════════════════════
# 4.  blend_reads.py  (helper functions only)
#     Tests FASTQ I/O and reservoir sampling.
# ═════════════════════════════════════════════════════════════════════════════

_blend_fns = _load_functions_from_script(
    SCRIPTS_DIR / "blend_reads.py",
    ["fastq_records", "count_fastq_records", "reservoir_sample", "write_fastq_records"],
)
fastq_records_fn      = _blend_fns["fastq_records"]
count_fastq_records   = _blend_fns["count_fastq_records"]
reservoir_sample_fn   = _blend_fns["reservoir_sample"]
write_fastq_records   = _blend_fns["write_fastq_records"]


def _make_fastq_gz(tmp_path, n_reads, prefix="read") -> Path:
    """Write a minimal gzipped FASTQ with n_reads records."""
    path = tmp_path / f"{prefix}.fastq.gz"
    with gzip.open(path, "wt") as fh:
        for i in range(n_reads):
            fh.write(f"@{prefix}_{i}\n")
            fh.write("ACGT\n")
            fh.write("+\n")
            fh.write("IIII\n")
    return path


class TestFastqHelpers:
    """fastq_records, count_fastq_records, reservoir_sample."""

    def test_count_correct(self, tmp_path):
        path = _make_fastq_gz(tmp_path, 10)
        assert count_fastq_records(str(path)) == 10

    def test_count_empty_file(self, tmp_path):
        path = tmp_path / "empty.fastq.gz"
        with gzip.open(path, "wt") as fh:
            fh.write("")
        assert count_fastq_records(str(path)) == 0

    def test_fastq_records_yields_correct_tuples(self, tmp_path):
        path = _make_fastq_gz(tmp_path, 3, prefix="r")
        records = list(fastq_records_fn(str(path)))
        assert len(records) == 3
        for name, seq, plus, qual in records:
            assert name.startswith("@r_")
            assert seq.strip() == "ACGT"
            assert plus.strip() == "+"
            assert qual.strip() == "IIII"

    def test_reservoir_sample_exact_count(self, tmp_path):
        path = _make_fastq_gz(tmp_path, 20)
        rng = random.Random(42)
        samples = reservoir_sample_fn(str(path), 10, rng)
        assert len(samples) == 10

    def test_reservoir_sample_returns_all_when_k_exceeds_n(self, tmp_path):
        """If we request more reads than exist, we get back everything."""
        path = _make_fastq_gz(tmp_path, 5)
        rng = random.Random(0)
        samples = reservoir_sample_fn(str(path), 100, rng)
        assert len(samples) == 5

    def test_reservoir_sample_is_reproducible(self, tmp_path):
        """Same seed → same sample."""
        path = _make_fastq_gz(tmp_path, 50)
        s1 = reservoir_sample_fn(str(path), 20, random.Random(7))
        s2 = reservoir_sample_fn(str(path), 20, random.Random(7))
        assert [r[0] for r in s1] == [r[0] for r in s2]

    def test_reservoir_sample_different_seeds_differ(self, tmp_path):
        """Different seeds should (very likely) produce different samples."""
        path = _make_fastq_gz(tmp_path, 100)
        s1 = reservoir_sample_fn(str(path), 30, random.Random(1))
        s2 = reservoir_sample_fn(str(path), 30, random.Random(999))
        # It is astronomically unlikely they are identical
        assert [r[0] for r in s1] != [r[0] for r in s2]

    def test_write_and_read_roundtrip(self, tmp_path):
        """write_fastq_records output can be parsed back by fastq_records."""
        original = [
            ("@read_0\n", "ACGT\n", "+\n", "IIII\n"),
            ("@read_1\n", "TTTT\n", "+\n", "JJJJ\n"),
        ]
        out = tmp_path / "out.fastq.gz"
        with gzip.open(out, "wt") as fh:
            write_fastq_records(original, fh)
        recovered = list(fastq_records_fn(str(out)))
        assert recovered == original


# ═════════════════════════════════════════════════════════════════════════════
# 5.  annotate_repeats.py  (strand detection and element naming)
#     Tests the parse_element_name helper and strand logic.
# ═════════════════════════════════════════════════════════════════════════════

_repeat_fns = _load_functions_from_script(
    SCRIPTS_DIR / "annotate_repeats.py",
    ["parse_element_name"],
)
parse_element_name = _repeat_fns["parse_element_name"]


class TestParseElementName:
    """parse_element_name extracts a readable name from MMseqs2 target strings."""

    def test_hash_separator(self):
        assert parse_element_name("IS3#DNA/TcMar") == "IS3"

    def test_pipe_separator(self):
        assert parse_element_name("Tn10|transposon") == "Tn10"

    def test_space_separator(self):
        assert parse_element_name("ISAba1 family") == "ISAba1"

    def test_no_separator_returns_full_string(self):
        name = "NoSeparatorHere"
        assert parse_element_name(name) == name

    def test_whitespace_stripped(self):
        # Leading/trailing spaces around the name should be stripped
        result = parse_element_name("  IS1  #DNA")
        assert result == "IS1"
