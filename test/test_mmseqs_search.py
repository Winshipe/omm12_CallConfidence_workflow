"""
Tests for workflow/scripts/mmseqs_search.py

MMseqs2 is not needed: a tiny fake `mmseqs` executable is put first on $PATH.
It records every call and, for `convertalis`, writes a deterministic m8 file,
so we can check the command lines the script builds and how it merges output.

Note: like annotate_repeats.py, no rule currently calls this script.
"""

import json
import os
import stat
import sys

from harness import prepend_path, read_tsv, run_script, tmpdir

FAKE_MMSEQS = r'''#!{python}
import json, os, sys
with open(os.environ["FAKE_MMSEQS_LOG"], "a") as fh:
    fh.write(json.dumps(sys.argv[1:]) + "\n")
if os.environ.get("FAKE_MMSEQS_FAIL") == sys.argv[1]:
    sys.stderr.write("simulated failure\n")
    sys.exit(3)
if sys.argv[1] == "convertalis":
    db, out = sys.argv[3], sys.argv[5]
    name = os.path.basename(db)
    with open(out, "w") as fh:
        fh.write(f"contig_1\t{{name}}_hitA\t0.9\t100\t1\t0\t1\t100\t1\t100\t1e-20\t150\n")
        fh.write(f"contig_2\t{{name}}_hitB\t0.8\t90\t2\t0\t300\t210\t5\t95\t1e-10\t90\n")
'''


def _setup(td):
    bindir = td / "bin"
    bindir.mkdir()
    exe = bindir / "mmseqs"
    exe.write_text(FAKE_MMSEQS.format(python=sys.executable))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return bindir, td / "calls.jsonl"


def _run(td, databases, fail=None):
    bindir, calls = _setup(td)
    out = td / "hits.m8"
    os.environ["FAKE_MMSEQS_LOG"] = str(calls)
    if fail:
        os.environ["FAKE_MMSEQS_FAIL"] = fail
    try:
        with prepend_path(bindir):
            run_script("mmseqs_search.py",
                       input={"query_db": td / "querydb"},
                       output={"hits": out},
                       params={"databases": [str(td / d) for d in databases],
                               "sensitivity": 7.5, "evalue": 1e-5,
                               "tmp": str(td / "tmp")},
                       threads=4, log=[td / "search.log"])
    finally:
        os.environ.pop("FAKE_MMSEQS_LOG", None)
        os.environ.pop("FAKE_MMSEQS_FAIL", None)
    invocations = [json.loads(l) for l in calls.read_text().splitlines()] if calls.exists() else []
    return out, invocations


def test_search_then_convert_per_database_with_expected_flags():
    with tmpdir() as td:
        _, calls = _run(td, ["mobileOG", "phrog"])
        assert [c[0] for c in calls] == ["search", "convertalis", "search", "convertalis"]
        search = calls[0]
        assert search[1] == str(td / "querydb" / "queryDB")
        assert search[2] == str(td / "mobileOG")
        opts = dict(zip(search[5::2], search[6::2]))
        assert opts == {"-s": "7.5", "-e": "1e-05", "--threads": "4", "--search-type": "3"}
        assert calls[2][2] == str(td / "phrog")


def test_merged_output_has_header_and_database_column():
    with tmpdir() as td:
        out, _ = _run(td, ["mobileOG", "phrog"])
        header, rows = read_tsv(out)
        assert header[-1] == "database" and len(header) == 13
        assert all(len(r) == 13 for r in rows)
        assert [r[-1] for r in rows] == ["mobileOG", "mobileOG", "phrog", "phrog"]
        assert rows[2][1] == "phrog_hitA"


def test_merged_output_feeds_annotate_repeats():
    """mmseqs_search.py's output is the documented input of annotate_repeats.py."""
    with tmpdir() as td:
        out, _ = _run(td, ["mobileOG"])
        annot = td / "annot.tsv"
        run_script("annotate_repeats.py", input={"hits": out}, output={"tsv": annot},
                   wildcards={"ref_id": "r"})
        _, rows = read_tsv(annot)
        assert [(r[0], r[1], r[2], r[3]) for r in rows] == [
            ("contig_1", "1", "100", "+"), ("contig_2", "210", "300", "-")]


def test_no_databases_gives_header_only():
    with tmpdir() as td:
        out, calls = _run(td, [])
        header, rows = read_tsv(out)
        assert calls == [] and rows == [] and header[0] == "query"


def test_failing_mmseqs_raises_and_logs_stderr():
    with tmpdir() as td:
        try:
            _run(td, ["mobileOG"], fail="search")
        except RuntimeError as exc:
            assert "mmseqs search" in str(exc)
        else:
            raise AssertionError("a failing mmseqs call must raise")
        assert "simulated failure" in (td / "search.log").read_text()
