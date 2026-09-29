"""
Tests for workflow/scripts/annotate_repeats.py

Note: no rule in workflow/rules/ currently calls this script (annotate.smk
uses `mmseqs easy-search` directly).  It is tested so it stays correct if it
is wired back in.
"""

from harness import load_defs, read_tsv, run_script, tmpdir, write

parse_element_name = load_defs("annotate_repeats.py")["parse_element_name"]

M8_HEADER = ("query\ttarget\tpident\talign_len\tmismatches\tgap_opens\t"
             "q_start\tq_end\tt_start\tt_end\tevalue\tbitscore\tdatabase\n")
OUT_HEADER = ["seq_id", "start", "end", "strand", "element_name", "database",
              "evalue", "bitscore"]


def _run(td, hit_rows):
    hits = write(td / "hits.m8", M8_HEADER + "".join("\t".join(map(str, r)) + "\n"
                                                    for r in hit_rows), dedent=False)
    out = td / "annot.tsv"
    run_script("annotate_repeats.py", input={"hits": hits}, output={"tsv": out},
               wildcards={"ref_id": "refA"})
    return read_tsv(out)


def _hit(query, target, q_start, q_end, db="mobileOG", evalue="1e-30", bits="200"):
    return [query, target, "0.99", "300", "1", "0", q_start, q_end, "1", "300",
            evalue, bits, db]


def test_parse_element_name_separators():
    assert parse_element_name("IS256#ISfinder") == "IS256"
    assert parse_element_name("Tn916|TnCentral|x") == "Tn916"
    assert parse_element_name("IS3 family transposase") == "IS3"
    assert parse_element_name("plain_name") == "plain_name"


def test_parse_element_name_separator_priority():
    """'#' is checked before '|' and ' ', whatever their position."""
    assert parse_element_name("a|b#c") == "a|b"
    assert parse_element_name("a b|c") == "a b"


def test_forward_and_reverse_hits_normalise_coordinates():
    with tmpdir() as td:
        header, rows = _run(td, [_hit("c1", "IS1#x", 100, 400),
                                 _hit("c1", "IS2#x", 900, 600)])
        assert header == OUT_HEADER
        assert rows[0][:5] == ["c1", "100", "400", "+", "IS1"]
        assert rows[1][:5] == ["c1", "600", "900", "-", "IS2"]
        for r in rows:
            assert int(r[1]) <= int(r[2])


def test_output_sorted_by_seq_then_numeric_start():
    """Sorting must be numeric (1000 after 200), not lexicographic."""
    with tmpdir() as td:
        _, rows = _run(td, [_hit("c2", "a", 5, 50), _hit("c1", "b", 1000, 1100),
                            _hit("c1", "c", 200, 300)])
        assert [(r[0], r[1]) for r in rows] == [("c1", "200"), ("c1", "1000"), ("c2", "5")]


def test_database_evalue_and_bitscore_passed_through():
    with tmpdir() as td:
        _, rows = _run(td, [_hit("c1", "x", 1, 2, db="phrog", evalue="3.2e-45", bits="512")])
        assert rows[0][5:] == ["phrog", "3.2e-45", "512"]


def test_single_base_hit_is_forward():
    with tmpdir() as td:
        _, rows = _run(td, [_hit("c1", "x", 7, 7)])
        assert rows[0][1:4] == ["7", "7", "+"]


def test_no_hits_writes_header_only():
    with tmpdir() as td:
        header, rows = _run(td, [])
        assert header == OUT_HEADER and rows == []
