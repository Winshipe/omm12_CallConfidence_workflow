"""
Tests for workflow/scripts/find_homologs.py

The script imports Biopython at module level, so every test here is skipped
when Bio is not importable.
"""

import contextlib
import io

from harness import SkipTest, has_module, load_defs, read_tsv, run_script, tmpdir, write

if has_module("Bio"):
    _ns = load_defs("find_homologs.py")
    find_homolog_ids = _ns["find_homolog_ids"]
    get_gene_coords = _ns["get_gene_coords"]
    read_fasta_files = _ns["read_fasta_files"]
    write_bed = _ns["write_bed"]


def _need_bio():
    if not has_module("Bio"):
        raise SkipTest("needs Python module(s): Bio")


def _bed_rows(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    assert lines and lines[0] == "#chrom\tstart\tstop\tname\tscore\tstrand", lines[:1]
    return [l.split("\t") for l in lines[1:]]


# ── find_homolog_ids ─────────────────────────────────────────────────────────

def test_find_homolog_ids_multi_member_cluster_and_singletons():
    _need_bio()
    with tmpdir() as td:
        f = write(td / "c.tsv", "geneA\tgeneA\ngeneA\tgeneB\ngeneC\tgeneC\n")
        assert find_homolog_ids(str(f)) == {"geneA", "geneB"}


def test_find_homolog_ids_empty_and_all_singletons():
    _need_bio()
    with tmpdir() as td:
        assert find_homolog_ids(str(write(td / "e.tsv", ""))) == set()
        assert find_homolog_ids(str(write(td / "s.tsv", "g1\tg1\ng2\tg2\n"))) == set()


def test_find_homolog_ids_many_clusters_and_blank_lines():
    _need_bio()
    with tmpdir() as td:
        f = write(td / "c.tsv", "r1\tr1\nr1\tm1\nr1\tm2\n\nr2\tr2\nr2\tm3\nr3\tr3\n\n")
        assert find_homolog_ids(str(f)) == {"r1", "m1", "m2", "r2", "m3"}


def test_find_homolog_ids_tolerates_crlf():
    _need_bio()
    with tmpdir() as td:
        f = td / "c.tsv"
        f.write_bytes(b"a\ta\r\na\tb\r\n")
        assert find_homolog_ids(str(f)) == {"a", "b"}


# ── get_gene_coords ──────────────────────────────────────────────────────────

def test_get_gene_coords_parses_prodigal_headers():
    _need_bio()
    with tmpdir() as td:
        fa = write(td / "g.fna", """\
            >contig_1_1 # 10 # 270 # 1 # ID=1_1;partial=00
            ATGATGATG
            >contig_1_2 # 500 # 800 # -1 # ID=1_2;partial=00
            ATGATGATG
        """)
        assert get_gene_coords(str(fa)) == {
            "contig_1_1": ("10", "270", "1"),
            "contig_1_2": ("500", "800", "-1"),
        }


def test_get_gene_coords_real_ncbi_style_contig_names():
    """NCBI accessions contain dots and underscores (e.g. NZ_CP065321.1_12)."""
    _need_bio()
    with tmpdir() as td:
        fa = write(td / "g.fna", ">NZ_CP065321.1_12 # 1 # 99 # 1 # ID=1_12\nATG\n")
        assert get_gene_coords(str(fa)) == {"NZ_CP065321.1_12": ("1", "99", "1")}


def test_get_gene_coords_skips_malformed_headers():
    _need_bio()
    with tmpdir() as td:
        fa = write(td / "bad.fna", ">no_prodigal_fields\nACGT\n>ok_1 # 1 # 3 # 1 # x\nATG\n")
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            assert get_gene_coords(str(fa)) == {"ok_1": ("1", "3", "1")}
        assert "no_prodigal_fields" in stdout.getvalue(), "malformed header should be reported"


def test_read_fasta_files_keys_by_path():
    _need_bio()
    with tmpdir() as td:
        a = write(td / "a.fna", ">a_1 # 1 # 3 # 1 # x\nATG\n")
        b = write(td / "b.fna", ">b_1 # 4 # 9 # -1 # x\nATG\n")
        got = read_fasta_files([str(a), str(b)])
        assert got == {str(a): {"a_1": ("1", "3", "1")}, str(b): {"b_1": ("4", "9", "-1")}}


# ── write_bed ────────────────────────────────────────────────────────────────

def test_write_bed_strand_symbols_chrom_and_sorting():
    _need_bio()
    with tmpdir() as td:
        entries = {"x.fna": {"contig_2_3": ("500", "900", "-1"),
                             "contig_1_1": ("100", "400", "1"),
                             "NZ_CP1.1_7": ("7", "70", "1")}}
        write_bed({"contig_2_3", "contig_1_1", "NZ_CP1.1_7"}, entries, str(td / "o.bed"))
        rows = _bed_rows(td / "o.bed")
        assert [r[3] for r in rows] == sorted(["contig_2_3", "contig_1_1", "NZ_CP1.1_7"])
        by_name = {r[3]: r for r in rows}
        assert by_name["contig_1_1"] == ["contig_1", "100", "400", "contig_1_1", ".", "+"]
        assert by_name["contig_2_3"] == ["contig_2", "500", "900", "contig_2_3", ".", "-"]
        assert by_name["NZ_CP1.1_7"][0] == "NZ_CP1.1"
        assert all(len(r) == 6 for r in rows)


def test_write_bed_skips_homologs_not_in_any_fasta():
    """Homologs from *other* references are expected in the cluster file."""
    _need_bio()
    with tmpdir() as td:
        write_bed({"here_1", "other_ref_1"}, {"a.fna": {"here_1": ("1", "9", "1")}},
                  str(td / "o.bed"))
        assert [r[3] for r in _bed_rows(td / "o.bed")] == ["here_1"]


def test_write_bed_empty_input_writes_header_only():
    _need_bio()
    with tmpdir() as td:
        write_bed(set(), {}, str(td / "o.bed"))
        assert _bed_rows(td / "o.bed") == []


# ── End-to-end (rule find_homologs) ──────────────────────────────────────────

def test_script_end_to_end_with_named_inputs():
    """The rule passes input as clusters=..., fa_files=...; the script reads
    them positionally (input[0], input[1:]).  Verify that contract holds."""
    _need_bio()
    with tmpdir() as td:
        clusters = write(td / "clusters.txt",
                         "refA_1\trefA_1\nrefA_1\trefB_4\nrefA_2\trefA_2\nrefA_3\trefA_3\n"
                         "refA_3\trefA_5\n")
        fna = write(td / "refA.fna", """\
            >refA_1 # 1 # 300 # 1 # ID=1_1
            ATG
            >refA_2 # 400 # 600 # -1 # ID=1_2
            ATG
            >refA_3 # 700 # 900 # -1 # ID=1_3
            ATG
            >refA_5 # 1000 # 1300 # 1 # ID=1_5
            ATG
        """)
        out = td / "refA.challenging.tsv"
        run_script("find_homologs.py",
                   input={"clusters": clusters, "fa_files": fna},
                   output={"ref_hits": out})
        rows = _bed_rows(out)
        assert [r[3] for r in rows] == ["refA_1", "refA_3", "refA_5"]
        assert [r[5] for r in rows] == ["+", "-", "+"]


def test_script_end_to_end_multiple_fasta_inputs():
    _need_bio()
    with tmpdir() as td:
        clusters = write(td / "c.txt", "a_1\ta_1\na_1\tb_1\n")
        a = write(td / "a.fna", ">a_1 # 1 # 9 # 1 # x\nATG\n")
        b = write(td / "b.fna", ">b_1 # 11 # 19 # -1 # x\nATG\n")
        out = td / "o.tsv"
        run_script("find_homologs.py", input=[clusters, a, b], output=[out])
        assert [(r[0], r[3], r[5]) for r in _bed_rows(out)] == [("a", "a_1", "+"), ("b", "b_1", "-")]
