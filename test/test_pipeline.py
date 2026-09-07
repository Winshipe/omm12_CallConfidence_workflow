#!/usr/bin/env python3
"""
test_pipeline.py
================
Plain-Python test suite for the CallConfidence variant-calling benchmark
pipeline.  No external test framework (pytest, unittest, etc.) is used —
every test is just a function; the runner at the bottom calls them all and
reports pass / fail counts.

How to run
----------
    python test_pipeline.py          # run everything
    python test_pipeline.py -v       # verbose: print each test name as it runs

Coverage
--------
  mutate_reference.py   – FASTA I/O, Jukes-Cantor model, Tamura-Nei model,
                          mutate_sequence(), biological properties of mutations
  find_homologs.py      – find_homolog_ids(), get_gene_coords(), write_bed()
  assess_variants.py    – the merge / classification logic (TP / FP / FN)
  Integration           – end-to-end: mutate → write TSV → assess

Structure
---------
Each test function:
  • Has a name starting with "test_"
  • Raises AssertionError (or any Exception) on failure
  • Returns normally on success

Helper utilities at the top create temporary files and directories and
clean them up at the end of every test, so tests are fully isolated.
"""

import io
import os
import random
import shutil
import sys
import tempfile
import textwrap
from pathlib import Path

# ── Make the project scripts importable ──────────────────────────────────────
# The scripts live in the project directory.  We import the pure-Python
# functions directly; Snakemake-specific top-level code is avoided by
# importing only the function definitions.

sys.path.insert(0, str(Path(__file__).parent))

# Import only the reusable functions from each script.
# (The scripts' top-level code references `snakemake`, so we cannot simply
#  `import mutate_reference`; instead we exec just the function definitions.)

class _SnakemakeStub:
    """
    Minimal stand-in for the `snakemake` object that Snakemake injects at
    runtime.  Scripts reference it at module level (e.g. for logging setup),
    so we need *something* in the namespace when we exec them in tests.
    None of the attribute values matter — we only care about extracting the
    function definitions that appear later in the same file.

    The one exception is `snakemake.log[0]`: logging.basicConfig writes to
    that path, so it must be a real (writable) file path rather than a
    directory.  We create a fresh temp file for each stub instance.
    """
    import tempfile as _tf

    class _Bunch:
        """Attribute bag: every attribute and index access returns '.'."""
        def __getattr__(self, _):
            return "."
        def __getitem__(self, _):
            return "."

    class _LogList:
        """Behaves like snakemake.log: indexing returns a temp-file path."""
        def __init__(self):
            import tempfile
            self._path = tempfile.mktemp(suffix=".log")
        def __getitem__(self, _):
            return self._path
        def __getattr__(self, _):
            return self._path

    def __init__(self):
        self.log = self._LogList()

    def __getattr__(self, name):
        if name == "log":          # already set in __init__
            raise AttributeError(name)
        return self._Bunch()

    def __getitem__(self, _):
        return "."


def _load_functions_from_script(script_path: str, names: list) -> dict:
    """
    Execute the given script in a namespace that already contains a stub
    `snakemake` object, then return the requested function objects.

    Scripts that call `snakemake.log[0]` or similar at module level will no
    longer crash; they'll get a harmless placeholder value from the stub.
    Any remaining errors (e.g. a top-level `main()` call that itself uses
    snakemake attributes in a way the stub doesn't cover) are silently
    swallowed — we only care about the function definitions.
    """
    with open(script_path) as fh:
        source = fh.read()

    ns = {"snakemake": _SnakemakeStub()}
    try:
        exec(compile(source, script_path, "exec"), ns)
    except Exception:
        pass  # tolerate top-level snakemake usage beyond what the stub covers

    missing = [n for n in names if n not in ns]
    if missing:
        raise ImportError(
            f"Could not find functions {missing} in {script_path}. "
            "Check that the function names match exactly."
        )
    return {n: ns[n] for n in names}


# Locate scripts relative to this test file
_PROJECT = Path(__file__).parent.parent.joinpath("workflow")

_mutate_fns = _load_functions_from_script(
    str(_PROJECT / "scripts" / "mutate_reference.py"),
    ["read_fasta", "write_fasta", "jukes_cantor_probs",
     "tamura_nei_probs", "mutate_sequence"],
)

_homolog_fns = _load_functions_from_script(
    str(_PROJECT / "scripts" / "find_homologs.py"),
    ["find_homolog_ids", "get_gene_coords", "write_bed"],
)

# assess_variants.py has no extractable pure functions — its logic is a
# short pandas pipeline.  We replicate / call that logic directly in tests.
import pandas as pd

# Unpack imported names for convenient use
read_fasta        = _mutate_fns["read_fasta"]
write_fasta       = _mutate_fns["write_fasta"]
jukes_cantor_probs = _mutate_fns["jukes_cantor_probs"]
tamura_nei_probs   = _mutate_fns["tamura_nei_probs"]
mutate_sequence    = _mutate_fns["mutate_sequence"]

find_homolog_ids = _homolog_fns["find_homolog_ids"]
get_gene_coords  = _homolog_fns["get_gene_coords"]
write_bed        = _homolog_fns["write_bed"]


# ═══════════════════════════════════════════════════════════════════════════════
# Shared helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _tmpdir():
    """Return a freshly created temporary directory path (caller must clean up)."""
    return tempfile.mkdtemp(prefix="callconf_test_")


def _write(path: str, text: str) -> str:
    """Write *text* to *path* and return the path."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(textwrap.dedent(text))
    return path


# ─────────────────────────────────────────────────────────────────────────────
# Replication of the assess_variants.py merge logic
# (extracted here so tests can call it without Snakemake)
# ─────────────────────────────────────────────────────────────────────────────

def _assess(ground_truth_df: pd.DataFrame,
            vcf_df: pd.DataFrame,
            min_quality: float,
            replicate: str = "rep1") -> pd.DataFrame:
    """
    Mirror the logic in assess_variants.py:
      1. Filter VCF rows below min_quality.
      2. Outer-join ground truth and VCF on (CHROM, POS, REF, ALT).
      3. Label each row TP / FP / FN.
      4. Add a replicate column.

    Parameters
    ----------
    ground_truth_df
        DataFrame with columns CHROM, POS, REF, ALT, mutation_type.
        (mutation_type is carried through but not used in the merge.)
    vcf_df
        DataFrame with columns CHROM, POS, ID, REF, ALT, QUAL, FILTER,
        INFO, FORMAT, scenario — matching the VCF column list in the script.
    min_quality
        Minimum QUAL value to keep a VCF row.
    replicate
        Label string added to every output row.

    Returns
    -------
    DataFrame with columns CHROM, POS, REF, ALT, Truthiness, replicate.
    """
    # Step 1 – quality filter
    passed = vcf_df[vcf_df["QUAL"] >= min_quality].copy()

    # Step 2 – outer join
    merged = pd.merge(
        ground_truth_df[["CHROM", "POS", "REF", "ALT"]],
        passed[["CHROM", "POS", "REF", "ALT"]],
        how="outer",
        on=["CHROM", "POS", "REF", "ALT"],
        indicator=True,
    )

    # Step 3 – label
    label_map = {"both": "TP", "right_only": "FP", "left_only": "FN"}
    merged["Truthiness"] = merged["_merge"].map(label_map)
    merged = merged.drop(columns=["_merge"])

    # Step 4 – replicate tag
    merged["replicate"] = replicate

    return merged


# ═══════════════════════════════════════════════════════════════════════════════
# ── mutate_reference.py tests ─────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

def test_read_fasta_single_record():
    """read_fasta correctly parses a single-record FASTA."""
    td = _tmpdir()
    try:
        fa = _write(f"{td}/seq.fa", """\
            >contig1 some description
            ACGTACGT
            NNNNACGT
        """)
        records = read_fasta(fa)
        assert len(records) == 1, f"Expected 1 record, got {len(records)}"
        header, seq = records[0]
        assert header == "contig1 some description", f"Unexpected header: {header!r}"
        # Sequence should be uppercased and concatenated
        assert seq == "ACGTACGTNNNNACGT", f"Unexpected sequence: {seq!r}"
    finally:
        shutil.rmtree(td)


def test_read_fasta_multi_record():
    """read_fasta handles multiple records in one file."""
    td = _tmpdir()
    try:
        fa = _write(f"{td}/multi.fa", """\
            >seq_A
            AAAA
            >seq_B
            CCCC
            >seq_C
            TTTT
        """)
        records = read_fasta(fa)
        assert len(records) == 3, f"Expected 3 records, got {len(records)}"
        assert records[0][0] == "seq_A"
        assert records[1][1] == "CCCC"
        assert records[2][1] == "TTTT"
    finally:
        shutil.rmtree(td)


def test_read_fasta_lowercases_are_uppercased():
    """read_fasta converts lowercase bases to uppercase (required for mutation model)."""
    td = _tmpdir()
    try:
        fa = _write(f"{td}/lower.fa", """\
            >contig
            acgtACGT
        """)
        _, seq = read_fasta(fa)[0]
        assert seq == "ACGTACGT", f"Lowercase not uppercased: {seq!r}"
    finally:
        shutil.rmtree(td)


def test_write_then_read_fasta_roundtrip():
    """write_fasta followed by read_fasta returns identical records."""
    td = _tmpdir()
    try:
        original = [("hdr1 extra", "ACGTACGT"), ("hdr2", "GGGGCCCC")]
        out_path = f"{td}/out.fa"
        write_fasta(out_path, original)
        recovered = read_fasta(out_path)
        assert len(recovered) == len(original)
        for (oh, os_), (rh, rs) in zip(original, recovered):
            assert rh == oh, f"Header mismatch: {rh!r} vs {oh!r}"
            assert rs == os_, f"Sequence mismatch for {oh!r}"
    finally:
        shutil.rmtree(td)


def test_write_fasta_line_wrapping():
    """write_fasta wraps long sequences at the specified line width."""
    td = _tmpdir()
    try:
        seq = "ACGT" * 30   # 120 bases
        out_path = f"{td}/wrapped.fa"
        write_fasta(out_path, [("hdr", seq)], line_width=10)
        with open(out_path) as fh:
            lines = [l.rstrip() for l in fh if not l.startswith(">")]
        assert all(len(l) <= 10 for l in lines), \
            f"Line too long: {max(len(l) for l in lines)}"
        assert "".join(lines) == seq
    finally:
        shutil.rmtree(td)


# ── Jukes-Cantor model ────────────────────────────────────────────────────────

def test_jukes_cantor_probs_sum_to_one():
    """Jukes-Cantor probabilities for alternatives must sum to 1.0."""
    for base in ("A", "C", "G", "T"):
        probs = jukes_cantor_probs(base)
        total = sum(probs.values())
        assert abs(total - 1.0) < 1e-9, \
            f"JC probs for {base} sum to {total}, not 1.0"


def test_jukes_cantor_excludes_ref_base():
    """Jukes-Cantor must never return the reference base as an alternative."""
    for base in ("A", "C", "G", "T"):
        probs = jukes_cantor_probs(base)
        assert base not in probs, \
            f"JC returned ref base {base!r} as an alternative"


def test_jukes_cantor_equal_rates():
    """Jukes-Cantor assigns equal probability to every alternative base."""
    for base in ("A", "C", "G", "T"):
        probs = jukes_cantor_probs(base)
        values = list(probs.values())
        assert all(abs(v - values[0]) < 1e-9 for v in values), \
            f"JC probabilities are not equal for ref={base}: {probs}"


# ── Tamura-Nei model ──────────────────────────────────────────────────────────

def _make_pi(gc: float) -> dict:
    """Convenience: build a base-frequency dict from GC content."""
    at = (1.0 - gc) / 2.0
    gc_ = gc / 2.0
    return {"A": at, "T": at, "G": gc_, "C": gc_}


def test_tamura_nei_probs_sum_to_one():
    """Tamura-Nei probabilities for alternatives must sum to 1.0."""
    pi = _make_pi(0.5)
    for base in ("A", "C", "G", "T"):
        probs = tamura_nei_probs(base, kappa=2.0, pi=pi)
        total = sum(probs.values())
        assert abs(total - 1.0) < 1e-9, \
            f"TN probs for {base} sum to {total}"


def test_tamura_nei_excludes_ref_base():
    """Tamura-Nei must never return the reference base as an alternative."""
    pi = _make_pi(0.5)
    for base in ("A", "C", "G", "T"):
        probs = tamura_nei_probs(base, kappa=2.0, pi=pi)
        assert base not in probs, \
            f"TN returned ref base {base!r} as an alternative"


def test_tamura_nei_kappa_favours_transitions():
    """
    With kappa > 1, the transition probability should exceed the transversion
    probability after accounting for base frequencies.

    A → G is a transition (both purines).
    A → C is a transversion.

    At equal base frequencies (pi=0.25 for all) and kappa=4:
        P(A→G) ∝ kappa * pi_G = 4 * 0.25 = 1.0
        P(A→C) ∝ pi_C          = 0.25
    So P(A→G) > P(A→C).
    """
    pi = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    probs = tamura_nei_probs("A", kappa=4.0, pi=pi)
    assert probs["G"] > probs["C"], \
        f"Expected P(A→G) > P(A→C) with kappa=4, but got {probs}"


def test_tamura_nei_kappa_1_equals_base_frequencies():
    """
    When kappa=1 (equal Ti/Tv rates) every alternative's weight equals its
    base frequency, so the normalised probabilities should be proportional to
    pi, regardless of whether the substitution is a transition or transversion.
    """
    pi = {"A": 0.1, "C": 0.4, "G": 0.1, "T": 0.4}
    probs = tamura_nei_probs("A", kappa=1.0, pi=pi)
    # With kappa=1 all weights equal pi[alt]; normalise manually
    expected_unnorm = {"C": 0.4, "G": 0.1, "T": 0.4}
    total = sum(expected_unnorm.values())
    expected = {b: v / total for b, v in expected_unnorm.items()}
    for base in ("C", "G", "T"):
        assert abs(probs[base] - expected[base]) < 1e-9, \
            f"TN kappa=1 gave {probs[base]:.4f} for {base}, expected {expected[base]:.4f}"


# ── mutate_sequence() ─────────────────────────────────────────────────────────

def test_mutate_sequence_rate_zero_no_mutations():
    """At rate=0 no mutations should be introduced."""
    rng = random.Random(42)
    seq = "ACGTACGTACGT"
    mutated, mutations = mutate_sequence(seq, "jukes_cantor", rate=0.0,
                                         kappa=2.0, gc_freq=0.5, rng=rng)
    assert mutated == seq, "Rate=0 should produce no changes"
    assert mutations == [], "Rate=0 should produce no mutation records"


def test_mutate_sequence_rate_one_all_sites_mutated():
    """
    At rate=1.0 every unambiguous site must be mutated.
    (The output base will differ from the input because the substitution
     model always picks an *alternative* base.)
    """
    rng = random.Random(42)
    seq = "ACGTACGT"   # 8 unambiguous bases
    mutated, mutations = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                         kappa=2.0, gc_freq=0.5, rng=rng)
    assert len(mutations) == 8, \
        f"Expected 8 mutations at rate=1, got {len(mutations)}"
    # Every position must differ from the original
    for ref, alt in zip(seq, mutated):
        assert ref != alt, f"Position with ref={ref!r} was not mutated"


def test_mutate_sequence_skips_ambiguous_bases():
    """
    The mutator must skip 'N' (and any character not in ACGT) and never
    record mutations at those positions.
    """
    rng = random.Random(0)
    seq = "ANNNACGT"
    _, mutations = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                   kappa=2.0, gc_freq=0.5, rng=rng)
    positions = {m["position"] for m in mutations}
    # Positions 2, 3, 4 are 'N' (1-based); they must not appear
    for bad_pos in (2, 3, 4):
        assert bad_pos not in positions, \
            f"Ambiguous position {bad_pos} was mutated"


def test_mutate_sequence_positions_are_1_based():
    """Mutation records must use 1-based coordinates, matching VCF convention."""
    rng = random.Random(42)
    seq = "ACGT"
    _, mutations = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                   kappa=2.0, gc_freq=0.5, rng=rng)
    positions = {m["position"] for m in mutations}
    assert min(positions) >= 1, "Positions should be ≥ 1 (1-based)"
    assert max(positions) <= len(seq), \
        f"Position {max(positions)} exceeds sequence length {len(seq)}"


def test_mutate_sequence_mutation_records_match_output_sequence():
    """
    Each mutation record should accurately reflect what ended up in the
    mutated sequence: ref_base matches the original and alt_base matches
    the changed position.
    """
    rng = random.Random(99)
    seq = "ACGTACGT"
    mutated, mutations = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                         kappa=2.0, gc_freq=0.5, rng=rng)
    for m in mutations:
        pos_0 = m["position"] - 1   # back to 0-based for string indexing
        assert m["ref_base"] == seq[pos_0], \
            f"ref_base mismatch at position {m['position']}: " \
            f"record says {m['ref_base']!r}, original has {seq[pos_0]!r}"
        assert m["alt_base"] == mutated[pos_0], \
            f"alt_base mismatch at position {m['position']}: " \
            f"record says {m['alt_base']!r}, mutated has {mutated[pos_0]!r}"


def test_mutate_sequence_reproducible_with_same_seed():
    """Two runs with the same RNG seed must produce identical output."""
    seq = "ACGT" * 50
    r1 = random.Random(7)
    r2 = random.Random(7)
    mut1, muts1 = mutate_sequence(seq, "tamura_nei", 0.01, 2.0, 0.5, r1)
    mut2, muts2 = mutate_sequence(seq, "tamura_nei", 0.01, 2.0, 0.5, r2)
    assert mut1 == mut2, "Same seed gave different mutated sequences"
    assert muts1 == muts2, "Same seed gave different mutation lists"


def test_mutate_sequence_different_seeds_differ():
    """
    Two runs with different seeds should (with overwhelming probability on a
    200-base sequence) produce different results.  The test is probabilistic
    but the probability of a false failure is astronomically small.
    """
    seq = "ACGT" * 50
    r1 = random.Random(1)
    r2 = random.Random(2)
    mut1, _ = mutate_sequence(seq, "jukes_cantor", 0.1, 2.0, 0.5, r1)
    mut2, _ = mutate_sequence(seq, "jukes_cantor", 0.1, 2.0, 0.5, r2)
    assert mut1 != mut2, "Different seeds produced identical sequences (very unlikely)"


def test_mutate_sequence_output_length_preserved():
    """Mutations are substitutions only; sequence length must not change."""
    seq = "ACGT" * 25
    rng = random.Random(42)
    mutated, _ = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                  kappa=2.0, gc_freq=0.5, rng=rng)
    assert len(mutated) == len(seq), \
        f"Sequence length changed: {len(seq)} → {len(mutated)}"


def test_mutate_sequence_tamura_nei_transition_bias():
    """
    Over many sites with kappa=10, transitions (Ti) should vastly outnumber
    transversions (Tv), as expected biologically and by the TN93 model.

    Ti: A↔G, C↔T   (same biochemical class)
    Tv: A↔C, A↔T, G↔C, G↔T  (different classes)
    """
    PURINES = {"A", "G"}
    seq = "ACGT" * 500   # 2 000 bases — large enough to detect the bias
    rng = random.Random(42)
    _, mutations = mutate_sequence(seq, "tamura_nei", rate=1.0,
                                   kappa=10.0, gc_freq=0.5, rng=rng)
    ti = sum(
        1 for m in mutations
        if (m["ref_base"] in PURINES) == (m["alt_base"] in PURINES)
    )
    tv = len(mutations) - ti
    assert ti > tv, f"Expected Ti > Tv with kappa=10, but Ti={ti}, Tv={tv}"


# ═══════════════════════════════════════════════════════════════════════════════
# ── find_homologs.py tests ────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

def test_find_homolog_ids_basic():
    """
    A cluster file with one multi-member cluster should return both IDs,
    while singleton entries (rep == member) are ignored.
    """
    td = _tmpdir()
    try:
        cluster_file = _write(f"{td}/clusters.tsv", """\
            geneA\tgeneA
            geneA\tgeneB
            geneC\tgeneC
        """)
        homologs = find_homolog_ids(cluster_file)
        assert "geneA" in homologs, "Representative of multi-member cluster missing"
        assert "geneB" in homologs, "Member of multi-member cluster missing"
        assert "geneC" not in homologs, "Singleton should be excluded"
    finally:
        shutil.rmtree(td)


def test_find_homolog_ids_empty_file():
    """An empty cluster file should return an empty set (not raise)."""
    td = _tmpdir()
    try:
        cluster_file = _write(f"{td}/empty.tsv", "")
        homologs = find_homolog_ids(cluster_file)
        assert homologs == set(), f"Expected empty set, got {homologs}"
    finally:
        shutil.rmtree(td)


def test_find_homolog_ids_all_singletons():
    """When every cluster has exactly one member, return an empty set."""
    td = _tmpdir()
    try:
        cluster_file = _write(f"{td}/singletons.tsv", """\
            gene1\tgene1
            gene2\tgene2
            gene3\tgene3
        """)
        homologs = find_homolog_ids(cluster_file)
        assert homologs == set()
    finally:
        shutil.rmtree(td)


def test_find_homolog_ids_many_members_same_cluster():
    """All members of a multi-member cluster must appear in the output set."""
    td = _tmpdir()
    try:
        # MMseqs2 lists each member on its own row with the same representative
        cluster_file = _write(f"{td}/multi.tsv", """\
            rep1\trep1
            rep1\tmemberA
            rep1\tmemberB
            rep1\tmemberC
        """)
        homologs = find_homolog_ids(cluster_file)
        for expected in ("rep1", "memberA", "memberB", "memberC"):
            assert expected in homologs, f"{expected!r} missing from homologs"
    finally:
        shutil.rmtree(td)


def test_get_gene_coords_forward_strand():
    """
    get_gene_coords correctly parses a Prodigal-style FASTA header and
    returns (start, stop, strand) for a forward-strand gene.

    Prodigal header format:
        >gene_name # start # stop # strand # attributes
    """
    from Bio import SeqIO
    td = _tmpdir()
    try:
        # Write a minimal Prodigal-style .fna with a single gene
        fa_path = _write(f"{td}/genes.fna", """\
            >contig_1_1 # 10 # 270 # 1 # ID=1_1;partial=00
            ATGATGATG
        """)
        coords = get_gene_coords(fa_path)
        assert "contig_1_1" in coords, "Gene ID not found in coords dict"
        start, stop, strand = coords["contig_1_1"]
        assert start == "10",   f"Expected start=10, got {start!r}"
        assert stop  == "270",  f"Expected stop=270, got {stop!r}"
        assert strand == "1",   f"Expected strand=1, got {strand!r}"
    finally:
        shutil.rmtree(td)


def test_get_gene_coords_reverse_strand():
    """Reverse-strand genes (strand = -1) are parsed correctly."""
    td = _tmpdir()
    try:
        fa_path = _write(f"{td}/rev.fna", """\
            >contig_99_3 # 500 # 800 # -1 # ID=99_3;partial=00
            ATGATGATG
        """)
        coords = get_gene_coords(fa_path)
        _, _, strand = coords["contig_99_3"]
        assert strand == "-1", f"Expected strand=-1, got {strand!r}"
    finally:
        shutil.rmtree(td)


def test_get_gene_coords_malformed_header_skipped():
    """
    A header that lacks the Prodigal ' # ' delimiter should be silently
    skipped rather than raising an exception.
    """
    td = _tmpdir()
    try:
        fa_path = _write(f"{td}/bad.fna", """\
            >gene_without_prodigal_format
            ACGT
        """)
        coords = get_gene_coords(fa_path)
        # Should return an empty dict (the bad record was skipped)
        assert coords == {}, f"Expected empty dict for malformed header, got {coords}"
    finally:
        shutil.rmtree(td)


def test_write_bed_content_and_strand_symbol():
    """
    write_bed should produce one BED-style line per homolog, converting
    the Prodigal strand integer (1 / -1) to BED convention (+ / -).
    """
    td = _tmpdir()
    try:
        homologs = {"contig_1_1", "contig_2_3"}
        fasta_entries = {
            "fake_path.fna": {
                "contig_1_1": ("100", "400", "1"),
                "contig_2_3": ("500", "900", "-1"),
            }
        }
        out_path = f"{td}/out.bed"
        write_bed(homologs, fasta_entries, out_path)

        with open(out_path) as fh:
            lines = [l.strip() for l in fh if not l.startswith("#")]

        assert len(lines) == 2, f"Expected 2 data lines, got {len(lines)}"

        # Build a lookup by gene name (column 4)
        by_name = {l.split("\t")[3]: l.split("\t") for l in lines}

        assert "contig_1_1" in by_name
        assert by_name["contig_1_1"][5] == "+", "Forward strand should be '+'"

        assert "contig_2_3" in by_name
        assert by_name["contig_2_3"][5] == "-", "Reverse strand should be '-'"
    finally:
        shutil.rmtree(td)


def test_write_bed_chrom_derived_from_gene_name():
    """
    The chrom column is derived by stripping the trailing gene-index suffix
    from the Prodigal gene name (e.g. 'contig_23611_5' → 'contig_23611').
    """
    td = _tmpdir()
    try:
        homologs = {"contig_23611_5"}
        fasta_entries = {
            "a.fna": {
                "contig_23611_5": ("1", "300", "1"),
            }
        }
        out_path = f"{td}/chrom_test.bed"
        write_bed(homologs, fasta_entries, out_path)

        with open(out_path) as fh:
            lines = [l for l in fh if not l.startswith("#")]

        assert len(lines) == 1
        chrom = lines[0].split("\t")[0]
        assert chrom == "contig_23611", \
            f"Expected chrom='contig_23611', got {chrom!r}"
    finally:
        shutil.rmtree(td)


def test_write_bed_header_line_present():
    """The output file must start with a comment header line."""
    td = _tmpdir()
    try:
        out_path = f"{td}/empty.bed"
        write_bed(set(), {}, out_path)
        with open(out_path) as fh:
            first_line = fh.readline()
        assert first_line.startswith("#"), \
            f"Expected a '#' comment header, got: {first_line!r}"
    finally:
        shutil.rmtree(td)


def test_write_bed_gene_missing_from_fasta_skipped():
    """
    If a homolog gene ID is not found in any FASTA entry dict,
    no line should be written for it (silently skipped).
    """
    td = _tmpdir()
    try:
        homologs = {"gene_exists", "gene_missing"}
        fasta_entries = {
            "a.fna": {
                "gene_exists": ("1", "100", "1"),
                # gene_missing is deliberately absent
            }
        }
        out_path = f"{td}/skip.bed"
        write_bed(homologs, fasta_entries, out_path)

        with open(out_path) as fh:
            lines = [l for l in fh if not l.startswith("#")]
        assert len(lines) == 1, \
            f"Expected 1 data line (for gene_exists only), got {len(lines)}"
    finally:
        shutil.rmtree(td)


# ═══════════════════════════════════════════════════════════════════════════════
# ── assess_variants.py logic tests ───────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

def _make_ground_truth(rows):
    """
    Convenience: build a ground-truth DataFrame from a list of dicts.
    Each dict must have keys: CHROM, POS, REF, ALT
    (mutation_type is optional and not needed for the merge).
    Always returns a DataFrame with the correct columns, even when rows=[].
    """
    cols = ["CHROM", "POS", "REF", "ALT"]
    if rows:
        return pd.DataFrame(rows, columns=cols)
    return pd.DataFrame(columns=cols)


def _make_vcf(rows):
    """
    Convenience: build a minimal VCF-like DataFrame.
    Each dict must have at minimum: CHROM, POS, REF, ALT, QUAL.
    Always returns a DataFrame with the full VCF column set, even when rows=[].
    """
    vcf_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                "FILTER", "INFO", "FORMAT", "scenario"]
    full_rows = []
    for r in rows:
        full_rows.append({
            "CHROM":    r["CHROM"],
            "POS":      r["POS"],
            "ID":       ".",
            "REF":      r["REF"],
            "ALT":      r["ALT"],
            "QUAL":     r["QUAL"],
            "FILTER":   ".",
            "INFO":     ".",
            "FORMAT":   "GT",
            "scenario": "test_scenario",
        })
    if full_rows:
        return pd.DataFrame(full_rows, columns=vcf_cols)
    return pd.DataFrame(columns=vcf_cols)


def test_assess_true_positive():
    """A variant present in both ground truth and VCF (above threshold) is a TP."""
    gt  = _make_ground_truth([{"CHROM": "chr1", "POS": 100, "REF": "A", "ALT": "T"}])
    vcf = _make_vcf([{"CHROM": "chr1", "POS": 100, "REF": "A", "ALT": "T", "QUAL": 50}])
    result = _assess(gt, vcf, min_quality=20)
    tp_rows = result[result["Truthiness"] == "TP"]
    assert len(tp_rows) == 1, f"Expected 1 TP, got {len(tp_rows)}"


def test_assess_false_negative_below_quality():
    """
    A variant in the ground truth that is also in the VCF but *below* the
    quality threshold should be classified as FN (same as not being called).
    The low-quality VCF call is filtered out before the merge, so the
    ground-truth row has no match and becomes FN.
    """
    gt  = _make_ground_truth([{"CHROM": "chr1", "POS": 200, "REF": "C", "ALT": "G"}])
    vcf = _make_vcf([{"CHROM": "chr1", "POS": 200, "REF": "C", "ALT": "G", "QUAL": 10}])
    result = _assess(gt, vcf, min_quality=20)
    fn_rows = result[result["Truthiness"] == "FN"]
    assert len(fn_rows) == 1, \
        f"Low-quality call should be FN, got: {result['Truthiness'].tolist()}"


def test_assess_false_negative_not_called():
    """A ground-truth mutation with no VCF match at all is a FN."""
    gt  = _make_ground_truth([{"CHROM": "chr1", "POS": 300, "REF": "G", "ALT": "A"}])
    vcf = _make_vcf([])   # empty VCF
    result = _assess(gt, vcf, min_quality=20)
    fn_rows = result[result["Truthiness"] == "FN"]
    assert len(fn_rows) == 1, f"Missing call should be FN, got {result['Truthiness'].tolist()}"


def test_assess_false_positive():
    """A VCF call (above threshold) with no ground-truth match is a FP."""
    gt  = _make_ground_truth([])   # nothing expected
    vcf = _make_vcf([{"CHROM": "chr1", "POS": 400, "REF": "T", "ALT": "C", "QUAL": 60}])
    result = _assess(gt, vcf, min_quality=20)
    fp_rows = result[result["Truthiness"] == "FP"]
    assert len(fp_rows) == 1, f"Novel call should be FP, got {result['Truthiness'].tolist()}"


def test_assess_mixed_tp_fn_fp():
    """
    A realistic mixture: one TP, one FN, one FP in the same run.
    """
    gt = _make_ground_truth([
        {"CHROM": "chr1", "POS": 100, "REF": "A", "ALT": "T"},  # will be TP
        {"CHROM": "chr1", "POS": 200, "REF": "C", "ALT": "G"},  # will be FN (not called)
    ])
    vcf = _make_vcf([
        {"CHROM": "chr1", "POS": 100, "REF": "A", "ALT": "T", "QUAL": 50},  # TP
        {"CHROM": "chr1", "POS": 999, "REF": "G", "ALT": "A", "QUAL": 40},  # FP
    ])
    result = _assess(gt, vcf, min_quality=20)
    counts = result["Truthiness"].value_counts().to_dict()
    assert counts.get("TP", 0) == 1, f"Expected 1 TP, got {counts}"
    assert counts.get("FN", 0) == 1, f"Expected 1 FN, got {counts}"
    assert counts.get("FP", 0) == 1, f"Expected 1 FP, got {counts}"


def test_assess_quality_threshold_boundary():
    """
    A call with QUAL exactly equal to min_quality should be retained (≥, not >).
    """
    gt  = _make_ground_truth([{"CHROM": "chr1", "POS": 500, "REF": "A", "ALT": "C"}])
    vcf = _make_vcf([{"CHROM": "chr1", "POS": 500, "REF": "A", "ALT": "C", "QUAL": 20}])
    result = _assess(gt, vcf, min_quality=20)   # exactly at the boundary
    assert result.iloc[0]["Truthiness"] == "TP", \
        "Call at exactly min_quality should be classified TP"


def test_assess_alt_allele_must_match():
    """
    A VCF call at the correct position but with the *wrong* ALT allele does
    NOT count as a true positive — the merge is on (CHROM, POS, REF, ALT).
    The expected variant becomes FN; the wrong call becomes FP.
    """
    gt  = _make_ground_truth([{"CHROM": "chr1", "POS": 600, "REF": "A", "ALT": "T"}])
    vcf = _make_vcf([{"CHROM": "chr1", "POS": 600, "REF": "A", "ALT": "G", "QUAL": 50}])
    result = _assess(gt, vcf, min_quality=20)
    counts = result["Truthiness"].value_counts().to_dict()
    assert counts.get("FN", 0) == 1, \
        f"Wrong ALT should give FN for expected mutation, got {counts}"
    assert counts.get("FP", 0) == 1, \
        f"Wrong ALT call should give FP, got {counts}"
    assert counts.get("TP", 0) == 0, \
        f"Should be no TP when ALT mismatches, got {counts}"


def test_assess_replicate_label_propagated():
    """The replicate label must appear on every row of the output."""
    gt  = _make_ground_truth([{"CHROM": "chr1", "POS": 100, "REF": "A", "ALT": "T"}])
    vcf = _make_vcf([{"CHROM": "chr1", "POS": 100, "REF": "A", "ALT": "T", "QUAL": 50}])
    result = _assess(gt, vcf, min_quality=20, replicate="rep3")
    assert (result["replicate"] == "rep3").all(), \
        "Replicate label should appear on all rows"


def test_assess_empty_ground_truth_all_fp():
    """With an empty ground truth, every VCF call should be a FP."""
    gt  = _make_ground_truth([])
    vcf = _make_vcf([
        {"CHROM": "chr1", "POS": 1, "REF": "A", "ALT": "T", "QUAL": 30},
        {"CHROM": "chr1", "POS": 2, "REF": "C", "ALT": "G", "QUAL": 40},
    ])
    result = _assess(gt, vcf, min_quality=20)
    assert (result["Truthiness"] == "FP").all(), \
        f"All rows should be FP, got {result['Truthiness'].tolist()}"


def test_assess_empty_vcf_all_fn():
    """With an empty (or entirely below-threshold) VCF, every mutation is FN."""
    gt = _make_ground_truth([
        {"CHROM": "chr1", "POS": 10, "REF": "A", "ALT": "G"},
        {"CHROM": "chr1", "POS": 20, "REF": "T", "ALT": "C"},
    ])
    vcf = _make_vcf([])
    result = _assess(gt, vcf, min_quality=20)
    assert (result["Truthiness"] == "FN").all(), \
        f"All rows should be FN with empty VCF, got {result['Truthiness'].tolist()}"


# ═══════════════════════════════════════════════════════════════════════════════
# ── Integration tests ─────────────────────────────────────────────────────────
# Tests that exercise multiple pipeline functions together.
# ═══════════════════════════════════════════════════════════════════════════════

def test_integration_mutate_then_assess_all_detected():
    """
    End-to-end smoke test:
      1. Generate mutations (rate=1.0 → every site mutated).
      2. Build a synthetic VCF containing all those mutations above threshold.
      3. Run assess; every expected mutation should be a TP.
    """
    seq = "ACGTACGT" * 10   # 80 bases
    rng = random.Random(42)
    mutated, mutations = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                          kappa=2.0, gc_freq=0.5, rng=rng)

    # Build ground truth DataFrame (mimic the TSV columns)
    gt_rows = [
        {"CHROM": "contig1", "POS": m["position"], "REF": m["ref_base"], "ALT": m["alt_base"]}
        for m in mutations
    ]
    gt = pd.DataFrame(gt_rows)

    # Build a perfect VCF (all mutations, all high-quality)
    vcf_rows = [
        {"CHROM": "contig1", "POS": m["position"], "REF": m["ref_base"],
         "ALT": m["alt_base"], "QUAL": 100}
        for m in mutations
    ]
    vcf = _make_vcf(vcf_rows)

    result = _assess(gt, vcf, min_quality=20)
    n_tp = (result["Truthiness"] == "TP").sum()
    n_fn = (result["Truthiness"] == "FN").sum()
    n_fp = (result["Truthiness"] == "FP").sum()

    assert n_tp == len(mutations), f"Expected {len(mutations)} TPs, got {n_tp}"
    assert n_fn == 0, f"Expected 0 FNs, got {n_fn}"
    assert n_fp == 0, f"Expected 0 FPs, got {n_fp}"


def test_integration_mutate_then_assess_none_detected():
    """
    The opposite extreme: all mutations exist in the ground truth but the
    VCF is empty.  Every expected mutation should be a FN.
    """
    seq = "ACGTACGT" * 10
    rng = random.Random(0)
    _, mutations = mutate_sequence(seq, "jukes_cantor", rate=1.0,
                                   kappa=2.0, gc_freq=0.5, rng=rng)

    gt_rows = [
        {"CHROM": "contig1", "POS": m["position"], "REF": m["ref_base"], "ALT": m["alt_base"]}
        for m in mutations
    ]
    gt  = pd.DataFrame(gt_rows)
    vcf = _make_vcf([])

    result = _assess(gt, vcf, min_quality=20)
    assert (result["Truthiness"] == "FN").all(), \
        "With empty VCF all mutations should be FN"


def test_integration_fasta_roundtrip_then_mutate():
    """
    Write a FASTA, read it back, mutate, and verify:
      • Sequence length is preserved.
      • Mutation positions are within the sequence bounds.
      • ref_base in each mutation record matches the original sequence.
    """
    td = _tmpdir()
    try:
        seq = "ACGTTTGCAAGTCCGGAATCGTACGTTTGCAAGTC"
        fa_in  = f"{td}/ref.fa"
        fa_out = f"{td}/mut.fa"
        write_fasta(fa_in, [("test_contig", seq)])

        records = read_fasta(fa_in)
        assert len(records) == 1
        _, recovered_seq = records[0]
        assert recovered_seq == seq

        rng = random.Random(1)
        mutated, mutations = mutate_sequence(recovered_seq, "tamura_nei",
                                             rate=0.3, kappa=2.0,
                                             gc_freq=0.5, rng=rng)
        assert len(mutated) == len(seq)

        for m in mutations:
            p = m["position"] - 1
            assert 0 <= p < len(seq), f"Position {p} out of range"
            assert m["ref_base"] == seq[p], \
                f"ref_base mismatch at {m['position']}"
    finally:
        shutil.rmtree(td)


# ═══════════════════════════════════════════════════════════════════════════════
# ── Test runner ───────────────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

def _collect_tests(module_globals: dict) -> list:
    """Return all callables whose name starts with 'test_', in definition order."""
    return [v for k, v in module_globals.items() if k.startswith("test_") and callable(v)]


def _run_all(verbose: bool = False) -> None:
    """
    Execute every test function.  Print a summary line and exit with code 1
    if any tests fail, so CI pipelines pick up failures automatically.
    """
    tests   = _collect_tests(globals())
    passed  = []
    failed  = []

    print(f"\nCallConfidence pipeline test suite — {len(tests)} tests\n" + "─" * 60)

    for fn in tests:
        name = fn.__name__
        if verbose:
            print(f"  RUNNING  {name} ... ", end="", flush=True)
        try:
            fn()
            passed.append(name)
            if verbose:
                print("PASS")
        except Exception as exc:
            failed.append((name, exc))
            if verbose:
                print(f"FAIL\n           → {type(exc).__name__}: {exc}")

    # Summary
    print(f"\n{'─' * 60}")
    print(f"Results: {len(passed)} passed, {len(failed)} failed")

    if failed:
        print("\nFailed tests:")
        for name, exc in failed:
            print(f"  ✗  {name}")
            # Indent the error message for readability
            for line in str(exc).splitlines():
                print(f"       {line}")
        print()
        sys.exit(1)
    else:
        print("\nAll tests passed ✓\n")


if __name__ == "__main__":
    verbose = "-v" in sys.argv or "--verbose" in sys.argv
    _run_all(verbose=verbose)
