"""
Tests for workflow/scripts/mutate_reference.py

Unit tests call the script's functions directly (loaded without running its
top-level code); end-to-end tests run the whole script with a stub
``snakemake`` object, exactly as the ``mutate_reference`` rule would.
"""

import math
import random

from harness import (load_defs, read_tsv, run_script, tmpdir,
                     write, write_fasta)

_ns = load_defs("mutate_reference.py")
read_fasta = _ns["read_fasta"]
write_fasta_script = _ns["write_fasta"]
jukes_cantor_probs = _ns["jukes_cantor_probs"]
tamura_nei_probs = _ns["tamura_nei_probs"]
mutate_sequence = _ns["mutate_sequence"]

PURINES = {"A", "G"}


def _is_transition(ref, alt):
    return (ref in PURINES) == (alt in PURINES)


def _pi(gc):
    return {"A": (1 - gc) / 2, "T": (1 - gc) / 2, "G": gc / 2, "C": gc / 2}


def _run(td, fasta_records, model="tamura_nei", rate=0.01, kappa=2.0,
         gc_freq=0.5, seed=42, tag="run"):
    """Run the whole script; return (mutated_records, tsv_header, tsv_rows)."""
    ref = write_fasta(td / "ref.fa", fasta_records)
    out_fa, out_tsv = td / f"{tag}.mutated.fasta", td / f"{tag}.mutations.tsv"
    run_script(
        "mutate_reference.py",
        input={"ref": ref},
        output={"mutated_fasta": out_fa, "mutations_tsv": out_tsv},
        params={"model": model, "rate": rate, "kappa": kappa,
                "gc_freq": gc_freq, "seed": seed},
        log=[td / f"{tag}.log"],
    )
    header, rows = read_tsv(out_tsv)
    return read_fasta(str(out_fa)), header, rows


def _random_seq(n, seed=0, alphabet="ACGT"):
    rng = random.Random(seed)
    return "".join(rng.choice(alphabet) for _ in range(n))


# ── FASTA I/O ────────────────────────────────────────────────────────────────

def test_read_fasta_single_record():
    with tmpdir() as td:
        fa = write(td / "seq.fa", """\
            >contig1 some description
            ACGTACGT
            NNNNACGT
        """)
        assert read_fasta(str(fa)) == [("contig1 some description", "ACGTACGTNNNNACGT")]


def test_read_fasta_multi_record():
    with tmpdir() as td:
        fa = write(td / "multi.fa", ">seq_A\nAAAA\n>seq_B\nCCCC\n>seq_C\nTTTT\n")
        assert read_fasta(str(fa)) == [("seq_A", "AAAA"), ("seq_B", "CCCC"), ("seq_C", "TTTT")]


def test_read_fasta_uppercases_soft_masked_bases():
    """Soft-masked (lowercase) bases must be uppercased or they'd never mutate."""
    with tmpdir() as td:
        fa = write(td / "lower.fa", ">contig\nacgtACGT\n")
        assert read_fasta(str(fa))[0][1] == "ACGTACGT"


def test_read_fasta_handles_crlf_and_blank_lines():
    """Windows line endings / trailing blank lines must not leak into sequences."""
    with tmpdir() as td:
        fa = td / "crlf.fa"
        fa.write_bytes(b">c1\r\nACGT\r\nAC\r\n\r\n>c2\r\nGG\r\n")
        assert read_fasta(str(fa)) == [("c1", "ACGTAC"), ("c2", "GG")]


def test_read_fasta_empty_file_returns_empty_list():
    with tmpdir() as td:
        assert read_fasta(str(write(td / "e.fa", ""))) == []


def test_write_then_read_fasta_roundtrip():
    with tmpdir() as td:
        original = [("hdr1 extra", "ACGTACGT" * 20), ("hdr2", "GGGGCCCC")]
        write_fasta_script(str(td / "out.fa"), original)
        assert read_fasta(str(td / "out.fa")) == original


def test_write_fasta_line_wrapping():
    with tmpdir() as td:
        seq = "ACGT" * 30
        write_fasta_script(str(td / "w.fa"), [("hdr", seq)], line_width=10)
        lines = [l.rstrip() for l in open(td / "w.fa") if not l.startswith(">")]
        assert all(len(l) <= 10 for l in lines)
        assert "".join(lines) == seq


def test_write_fasta_exact_multiple_of_width_has_no_blank_line():
    with tmpdir() as td:
        write_fasta_script(str(td / "w.fa"), [("h", "A" * 20)], line_width=10)
        assert open(td / "w.fa").read() == ">h\n" + "A" * 10 + "\n" + "A" * 10 + "\n"


# ── Substitution models ──────────────────────────────────────────────────────

def test_jukes_cantor_probs_are_uniform_over_the_three_alternatives():
    for base in "ACGT":
        probs = jukes_cantor_probs(base)
        assert set(probs) == set("ACGT") - {base}
        assert all(math.isclose(p, 1 / 3) for p in probs.values())


def test_tamura_nei_probs_sum_to_one_and_exclude_ref():
    for gc in (0.2, 0.5, 0.7):
        for kappa in (0.5, 1.0, 2.0, 10.0):
            for base in "ACGT":
                probs = tamura_nei_probs(base, kappa=kappa, pi=_pi(gc))
                assert base not in probs
                assert math.isclose(sum(probs.values()), 1.0), (gc, kappa, base, probs)
                assert all(p > 0 for p in probs.values())


def test_tamura_nei_exact_values():
    """P(ref→alt) ∝ kappa·pi[alt] for transitions, pi[alt] for transversions."""
    pi = {"A": 0.1, "C": 0.2, "G": 0.3, "T": 0.4}
    probs = tamura_nei_probs("A", kappa=3.0, pi=pi)
    weights = {"C": 0.2, "G": 3.0 * 0.3, "T": 0.4}
    total = sum(weights.values())
    for b, w in weights.items():
        assert math.isclose(probs[b], w / total), (b, probs)
    probs = tamura_nei_probs("C", kappa=3.0, pi=pi)
    weights = {"A": 0.1, "G": 0.3, "T": 3.0 * 0.4}
    total = sum(weights.values())
    for b, w in weights.items():
        assert math.isclose(probs[b], w / total), (b, probs)


def test_tamura_nei_kappa_1_equal_pi_matches_jukes_cantor():
    for base in "ACGT":
        tn = tamura_nei_probs(base, kappa=1.0, pi=_pi(0.5))
        jc = jukes_cantor_probs(base)
        assert all(math.isclose(tn[b], jc[b]) for b in jc)


# ── mutate_sequence() ────────────────────────────────────────────────────────

def test_mutate_sequence_rate_zero_changes_nothing():
    seq = _random_seq(500)
    mutated, muts = mutate_sequence(seq, "jukes_cantor", 0.0, 2.0, 0.5, random.Random(1))
    assert mutated == seq and muts == []


def test_mutate_sequence_rate_one_mutates_every_unambiguous_site():
    seq = "ACGTNACGT-RY"
    mutated, muts = mutate_sequence(seq, "jukes_cantor", 1.0, 2.0, 0.5, random.Random(1))
    assert len(muts) == 8
    for i, (a, b) in enumerate(zip(seq, mutated)):
        if a in "ACGT":
            assert a != b, f"site {i + 1} not mutated"
        else:
            assert a == b, f"ambiguous/gap character at {i + 1} was changed"


def test_mutate_sequence_records_match_output_and_are_1_based():
    seq = _random_seq(2000, seed=3)
    mutated, muts = mutate_sequence(seq, "tamura_nei", 0.05, 2.0, 0.5, random.Random(9))
    changed = {i + 1 for i, (a, b) in enumerate(zip(seq, mutated)) if a != b}
    assert changed == {m["position"] for m in muts}, "records disagree with sequence diff"
    for m in muts:
        assert seq[m["position"] - 1] == m["ref_base"]
        assert mutated[m["position"] - 1] == m["alt_base"]
        assert m["ref_base"] != m["alt_base"]
        assert set(m) >= {"position", "ref_base", "alt_base"}
    positions = [m["position"] for m in muts]
    assert positions == sorted(positions) and len(set(positions)) == len(positions)


def test_mutate_sequence_preserves_length():
    seq = _random_seq(1000)
    for model in ("jukes_cantor", "tamura_nei"):
        mutated, _ = mutate_sequence(seq, model, 0.5, 2.0, 0.5, random.Random(0))
        assert len(mutated) == len(seq)


def test_mutate_sequence_same_seed_is_reproducible():
    seq = _random_seq(1000)
    a = mutate_sequence(seq, "tamura_nei", 0.02, 2.0, 0.5, random.Random(7))
    b = mutate_sequence(seq, "tamura_nei", 0.02, 2.0, 0.5, random.Random(7))
    assert a == b


def test_mutate_sequence_different_seeds_differ():
    seq = _random_seq(1000)
    a, _ = mutate_sequence(seq, "jukes_cantor", 0.1, 2.0, 0.5, random.Random(1))
    b, _ = mutate_sequence(seq, "jukes_cantor", 0.1, 2.0, 0.5, random.Random(2))
    assert a != b


def test_mutate_sequence_unknown_model_raises():
    try:
        mutate_sequence("ACGT", "hky85", 1.0, 2.0, 0.5, random.Random(0))
    except ValueError as exc:
        assert "hky85" in str(exc)
    else:
        raise AssertionError("unknown model should raise ValueError")


def test_mutate_sequence_empty_and_all_n_sequences():
    assert mutate_sequence("", "tamura_nei", 1.0, 2.0, 0.5, random.Random(0)) == ("", [])
    assert mutate_sequence("NNNN", "tamura_nei", 1.0, 2.0, 0.5, random.Random(0)) == ("NNNN", [])


# ── Statistical properties (fixed seeds, generous tolerances) ────────────────

def test_mutation_count_matches_rate():
    """Observed count must lie within 5 SD of the binomial expectation."""
    n, rate = 50_000, 0.01
    _, muts = mutate_sequence(_random_seq(n, 11), "tamura_nei", rate, 2.0, 0.5, random.Random(5))
    expected, sd = n * rate, math.sqrt(n * rate * (1 - rate))
    assert abs(len(muts) - expected) < 5 * sd, (len(muts), expected)


def test_jukes_cantor_alt_bases_are_uniform():
    _, muts = mutate_sequence("A" * 30_000, "jukes_cantor", 1.0, 2.0, 0.5, random.Random(5))
    counts = {b: sum(m["alt_base"] == b for m in muts) for b in "CGT"}
    for b, c in counts.items():
        assert abs(c / 30_000 - 1 / 3) < 0.02, counts


def test_tamura_nei_ti_tv_ratio_matches_kappa():
    """With equal base frequencies each site has 1 Ti and 2 Tv options, so the
    expected Ti:Tv ratio is kappa/2."""
    kappa = 6.0
    _, muts = mutate_sequence(_random_seq(40_000, 2), "tamura_nei", 1.0, kappa, 0.5,
                              random.Random(8))
    ti = sum(_is_transition(m["ref_base"], m["alt_base"]) for m in muts)
    ratio = ti / (len(muts) - ti)
    assert abs(ratio - kappa / 2) < 0.15, ratio


def test_tamura_nei_gc_freq_biases_alt_bases():
    """kappa=1, gc_freq=0.8: from A, P(C)=P(G)=0.4/0.9, P(T)=0.1/0.9."""
    _, muts = mutate_sequence("A" * 30_000, "tamura_nei", 1.0, 1.0, 0.8, random.Random(4))
    frac_t = sum(m["alt_base"] == "T" for m in muts) / len(muts)
    assert abs(frac_t - 0.1 / 0.9) < 0.015, frac_t


# ── End-to-end script runs (as Snakemake would call it) ──────────────────────

def test_script_writes_fasta_and_tsv_consistently():
    with tmpdir() as td:
        seqs = [("contig_1 desc here", _random_seq(3000, 1)),
                ("contig_2", _random_seq(1500, 2))]
        mutated, header, rows = _run(td, seqs, rate=0.02)

        assert [h for h, _ in mutated] == ["contig_1 desc here [mutated]", "contig_2 [mutated]"]
        assert header[:4] == ["seq_id", "position", "ref_base", "alt_base"]
        assert rows, "expected some mutations at rate 0.02 on 4.5 kb"

        originals = {h.split()[0]: s for h, s in seqs}
        muts_by_seq = {h.split()[0]: list(s) for h, s in seqs}
        for row in rows:
            seq_id, pos, ref, alt = row[0], int(row[1]), row[2], row[3]
            assert seq_id in originals, f"seq_id {seq_id!r} should be the first header token"
            assert originals[seq_id][pos - 1] == ref
            muts_by_seq[seq_id][pos - 1] = alt
        # Replaying the TSV onto the originals must reproduce the mutated FASTA.
        for (h, s) in mutated:
            assert "".join(muts_by_seq[h.split()[0]]) == s, f"TSV does not explain {h}"


def test_script_rate_zero_outputs_identical_sequence_and_header_only_tsv():
    with tmpdir() as td:
        seq = _random_seq(500)
        mutated, header, rows = _run(td, [("c", seq)], rate=0.0)
        assert mutated == [("c [mutated]", seq)]
        assert rows == [] and header[0] == "seq_id"


def test_script_writes_log():
    with tmpdir() as td:
        _run(td, [("c", _random_seq(200))], rate=0.1, tag="logged")
        log = (td / "logged.log").read_text()
        assert "Model: tamura_nei" in log and "Done." in log


def test_script_accepts_string_params_from_config():
    """Snakemake params may arrive as strings; the script casts them."""
    with tmpdir() as td:
        _, _, rows = _run(td, [("c", _random_seq(500))], rate="1.0", kappa="2",
                          gc_freq="0.5")
        assert len(rows) == 500


def test_script_on_bundled_reference_smoke():
    """Run on the first 20 kb of the shipped A. muris genome."""
    from harness import ROOT
    ref = ROOT / "resources" / "CP065321.fa"
    with tmpdir() as td:
        header, seq = read_fasta(str(ref))[0]
        _, _, rows = _run(td, [(header, seq[:20_000])], rate=0.001)
        assert all(r[0] == "CP065321.1" for r in rows)


def test_script_same_seed_gives_identical_mutations():
    with tmpdir() as td:
        seq = [("c", _random_seq(20_000, 5))]
        _, _, rows_a = _run(td, seq, rate=0.01, seed=123, tag="a")
        _, _, rows_b = _run(td, seq, rate=0.01, seed=123, tag="b")
        assert rows_a == rows_b, "same seed produced different mutation sets"


def test_script_tsv_rows_have_same_width_as_header():
    with tmpdir() as td:
        _, header, rows = _run(td, [("c", _random_seq(2000))], rate=0.05)
        bad = [r for r in rows if len(r) != len(header)]
        assert not bad, f"{len(bad)} rows have {len(bad[0])} fields, header has {len(header)}"
