"""
Tests for workflow/scripts/blend_reads.py

The script pools every read pair from every input (no subsampling), keeps
mates together, and shuffles the pairs with the configured seed.

Note: rule blend_reads in blend_reads.smk currently shells out to `seqtk`
instead of calling this script (that path is tested in test_rules.py).  These
tests keep the script honest in case it is wired back in.
"""

from collections import Counter

from harness import load_defs, read_fastq, run_script, tmpdir, write_fastq

_ns = load_defs("blend_reads.py")
fastq_records = _ns["fastq_records"]
read_pairs = _ns["read_pairs"]


def _inputs(td, scenario, n_reads=100):
    """Create ART-like paired FASTQs; read names encode their origin."""
    files = {}
    for i, c in enumerate(scenario):
        n = n_reads[i] if isinstance(n_reads, (list, tuple)) else n_reads
        for kind in ("unmutated", "mutated"):
            for mate in (1, 2):
                path = td / f"{c['ref_id']}_{kind}_R{mate}.fastq.gz"
                write_fastq(path, [(f"{c['ref_id']}.{kind}.{j}/{mate}", "ACGT" * 5)
                                   for j in range(n)])
                files[f"contrib_{i}_{kind}_R{mate}"] = path
    return files


def _blend(td, scenario, seed=1, n_reads=100, tag="out"):
    r1, r2 = td / f"{tag}_R1.fastq.gz", td / f"{tag}_R2.fastq.gz"
    run_script("blend_reads.py", input=_inputs(td, scenario, n_reads),
               output={"r1": r1, "r2": r2},
               params={"scenario_cfg": scenario, "seed": seed})
    return read_fastq(r1), read_fastq(r2)


def _origin(name):
    ref, kind, idx = name.lstrip("@").split("/")[0].split(".")
    return ref, kind, idx


def _fragment(rec):
    return rec[0].split("/")[0]


ONE_REF = [{"ref_id": "refA", "mutated_fraction": 0.5, "abundance": 1.0}]
TWO_REFS = [{"ref_id": "refA", "mutated_fraction": 0.5, "abundance": 3.0},
            {"ref_id": "refB", "mutated_fraction": 0.2, "abundance": 1.0}]


# ── helpers ──────────────────────────────────────────────────────────────────

def test_fastq_records_plain_and_gz():
    with tmpdir() as td:
        for name in ("r.fastq", "r.fastq.gz"):
            p = write_fastq(td / name, [(f"r{i}", "ACGT") for i in range(7)])
            recs = list(fastq_records(p))
            assert len(recs) == 7
            assert recs[0] == ("@r0\n", "ACGT\n", "+\n", "IIII\n")


def test_read_pairs_keeps_mates_in_lockstep():
    with tmpdir() as td:
        r1 = write_fastq(td / "a_R1.fq", [(f"frag{i}/1", "A") for i in range(20)])
        r2 = write_fastq(td / "a_R2.fq", [(f"frag{i}/2", "C") for i in range(20)])
        pairs = read_pairs(r1, r2)
        assert len(pairs) == 20
        assert all(_fragment(a) == _fragment(b) for a, b in pairs)


def test_read_pairs_rejects_unequal_mate_files():
    with tmpdir() as td:
        r1 = write_fastq(td / "a_R1.fq", [(f"f{i}/1", "A") for i in range(5)])
        r2 = write_fastq(td / "a_R2.fq", [(f"f{i}/2", "A") for i in range(4)])
        try:
            read_pairs(r1, r2)
        except ValueError as exc:
            assert "same number of reads" in str(exc)
        else:
            raise AssertionError("unequal R1/R2 files must raise")


# ── full script ──────────────────────────────────────────────────────────────

def test_every_input_read_is_used_exactly_once():
    with tmpdir() as td:
        r1, r2 = _blend(td, TWO_REFS, n_reads=(120, 45))
        expected = 2 * 120 + 2 * 45            # (unmutated + mutated) per ref
        assert len(r1) == len(r2) == expected
        names = Counter(_fragment(r) for r in r1)
        assert all(v == 1 for v in names.values()), "a read was written twice"
        by_source = Counter(_origin(r[0])[:2] for r in r1)
        assert by_source == {("refA", "unmutated"): 120, ("refA", "mutated"): 120,
                             ("refB", "unmutated"): 45, ("refB", "mutated"): 45}


def test_abundance_and_mutated_fraction_do_not_drop_reads():
    scenario = [{"ref_id": "refA", "mutated_fraction": 0.0, "abundance": 0.1},
                {"ref_id": "refB", "mutated_fraction": 1.0, "abundance": 10.0}]
    with tmpdir() as td:
        r1, _ = _blend(td, scenario, n_reads=30)
        assert len(r1) == 4 * 30


def test_mates_stay_together_after_shuffle():
    with tmpdir() as td:
        r1, r2 = _blend(td, TWO_REFS, n_reads=200)
        mism = sum(_fragment(a) != _fragment(b) for a, b in zip(r1, r2))
        assert mism == 0, f"{mism}/{len(r1)} pairs have mismatched mate names"
        assert all(a[0].endswith("/1") and b[0].endswith("/2") for a, b in zip(r1, r2))


def test_reads_are_shuffled_across_sources():
    with tmpdir() as td:
        r1, _ = _blend(td, TWO_REFS, n_reads=100)
        first_half_sources = {_origin(r[0])[:2] for r in r1[: len(r1) // 2]}
        assert len(first_half_sources) == 4, "output looks ordered by source file"


def test_records_are_written_intact():
    with tmpdir() as td:
        r1, r2 = _blend(td, ONE_REF)
        for rec in r1 + r2:
            assert rec[0].startswith("@") and rec[2] == "+" and len(rec[1]) == len(rec[3])
            assert rec[1] == "ACGT" * 5


def test_same_seed_is_reproducible_and_different_seed_differs():
    with tmpdir() as td:
        a = _blend(td, TWO_REFS, seed=5, tag="a")
        b = _blend(td, TWO_REFS, seed=5, tag="b")
        c = _blend(td, TWO_REFS, seed=6, tag="c")
        assert a == b
        assert a != c


def test_empty_inputs_give_empty_outputs():
    with tmpdir() as td:
        r1, r2 = _blend(td, ONE_REF, n_reads=0)
        assert r1 == [] and r2 == []
