"""
Tests for workflow/scripts/genome_report.Rmd

The report is knitted with knitr::knit() (which, unlike rmarkdown::render(),
needs no pandoc) on small synthetic inputs shaped like the real pipeline
outputs.  After knitting, the report's intermediate tables are exported from
the R session so their contents can be checked from Python.

Skipped when Rscript or the R packages the report loads are unavailable, or
when run with --no-r.
"""

import csv
import os
import shutil
import subprocess
import textwrap

from harness import SCRIPTS, SetupError, SkipTest, write

R_PACKAGES = ("knitr", "tidyverse", "cowplot", "kableExtra")
_R_OK = {}


def _need_r():
    if os.environ.get("CALLCONF_SKIP_R"):
        raise SkipTest("--no-r given")
    if "ok" not in _R_OK:
        rscript = shutil.which("Rscript")
        if not rscript:
            _R_OK["ok"] = "Rscript not found"
        else:
            check = "; ".join(f"if (!requireNamespace('{p}', quietly=TRUE)) cat('{p} ')"
                              for p in R_PACKAGES)
            out = subprocess.run([rscript, "-e", check], capture_output=True, text=True)
            _R_OK["ok"] = f"R packages missing: {out.stdout.strip()}" if out.stdout.strip() else ""
    if _R_OK["ok"]:
        raise SkipTest(_R_OK["ok"])


WINDOW = 1000


def _inputs(td):
    """Synthetic inputs for one reference 'refA' with two contigs."""
    # Prodigal GFF: contig refA_chr genes 1-3, plasmid refA_p1 gene 1
    gff = write(td / "refA.gff", textwrap.dedent("""\
        ##gff-version  3
        # Sequence Data: seqnum=1;seqlen=5000;seqhdr="refA_chr"
        refA_chr\tProdigal_v2.6.3\tCDS\t101\t700\t50.1\t+\t0\tID=1_1;partial=00;
        refA_chr\tProdigal_v2.6.3\tCDS\t1201\t1800\t40.2\t-\t0\tID=1_2;partial=00;
        refA_chr\tProdigal_v2.6.3\tCDS\t3001\t3900\t60.0\t+\t0\tID=1_3;partial=00;
        # Sequence Data: seqnum=2;seqlen=2000;seqhdr="refA_p1"
        refA_p1\tProdigal_v2.6.3\tCDS\t11\t610\t30.0\t+\t0\tID=2_1;partial=00;
        """))
    # Aggregated assessment (note: the duplicated replicate column is real)
    rows = []
    # window 0 (POS 0-999): 4 TP → recall 100 %
    rows += [("refA_chr", p, "A", "G", "TP") for p in (100, 200, 300, 400)]
    # window 1 (1000-1999): 1 TP, 3 FN → recall 25 %
    rows += [("refA_chr", 1100, "A", "G", "TP")]
    rows += [("refA_chr", p, "C", "T", "FN") for p in (1200, 1300, 1400)]
    # window 3 (3000-3999): 2 TP, 2 FP → recall 100 %, precision 50 %
    rows += [("refA_chr", p, "G", "A", "TP") for p in (3100, 3200)]
    rows += [("refA_chr", p, "T", "C", "FP") for p in (3300, 3400)]
    assess = td / "results" / "assessment" / "scenario_mix_all_replicates.tsv"
    write(assess, "CHROM\tPOS\tREF\tALT\tTruthiness\treplicate\treplicate\n" +
          "".join(f"{c}\t{p}\t{r}\t{a}\t{t}\trep1\trep1\n" for c, p, r, a, t in rows),
          dedent=False)
    challenging = write(td / "refA.challenging.tsv",
                        "#chrom\tstart\tstop\tname\tscore\tstrand\n"
                        "refA_chr\t1201\t1800\trefA_chr_2\t.\t-\n", dedent=False)
    # MMseqs easy-search (protein) hits: query, target, fident, alnlen, mismatch,
    # gapopen, qstart, qend, tstart, tend, evalue, bits — no header.
    hits = write(td / "mobileOG_hits",
                 "refA_chr_3\tmobileOG_000123|IS256\t0.95\t250\t5\t0\t11\t260\t1\t250\t1e-80\t400\n"
                 "refA_chr_3\tmobileOG_000999|IS3\t0.60\t250\t5\t0\t11\t260\t1\t250\t1e-20\t100\n"
                 "refA_chr_1\tmobileOG_000555|weak\t0.40\t120\t5\t0\t1\t120\t1\t120\t1e-5\t30\n"
                 # reverse-strand gene (1201-1800, '-') and a gene on the second contig
                 "refA_chr_2\tmobileOG_000777|Tn916\t0.90\t200\t5\t0\t11\t190\t1\t180\t1e-50\t300\n"
                 "refA_p1_1\tmobileOG_000888|IS3\t0.90\t200\t5\t0\t1\t100\t1\t100\t1e-50\t300\n",
                 dedent=False)
    sam_hdr = "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:refA_chr\tLN:5000\n"
    # source reads (ART .sam): read r{i} generated at pos 100*i, fragment 300
    src = "".join(f"r{i}\t99\trefA_chr\t{100 * i}\t99\t150M\t=\t{100 * i + 150}\t300\t"
                  f"{'A' * 150}\t{'I' * 150}\n" for i in range(1, 40))
    unmut = write(td / "unmut.sam", sam_hdr + src, dedent=False)
    # ART names reads identically in the mutated and unmutated simulations
    mut_src = "".join(f"r{i}\t99\trefA_chr\t{100 * i + 50}\t99\t150M\t=\t{100 * i + 200}\t300\t"
                      f"{'A' * 150}\t{'I' * 150}\n" for i in range(1, 40))
    mut = write(td / "mut.sam", sam_hdr + mut_src, dedent=False)
    # mapped reads (samtools view, no header, with tags): r1-r29 map where they
    # came from, r30-r39 map 2 kb away (cross-mapped within the same genome)
    mapped = "".join(
        f"r{i}\t99\trefA_chr\t{100 * i if i < 30 else 100 * i - 2000}\t60\t150M\t=\t0\t300\t"
        f"{'A' * 150}\t{'I' * 150}\tNM:i:0\tRG:Z:mix\n" for i in range(1, 40))
    mapped_sam = write(td / "results" / "variant_calling" / "mix" / "rep1" / "aligned.sam",
                       mapped, dedent=False)
    return dict(ref_id="refA", workingdir=str(td), scenario_list="",
                assessment_tsvs=str(assess.relative_to(td)), challenging_tsv=str(challenging),
                mapped_sam=str(mapped_sam), unmut_sam=str(unmut), mut_sam=str(mut),
                gff=str(gff), db_hit_files=str(hits), db_names="mobileOG",
                genome_length=5000, window_size=WINDOW)


def _knit(td, params):
    """Knit the report; export intermediate tables as CSV.  Returns (rc, log)."""
    r_params = ", ".join(
        f"{k} = {v}" if isinstance(v, int) else f"{k} = '{v}'" for k, v in params.items())
    script = textwrap.dedent(f"""\
        suppressMessages({{
          env <- new.env()
          env$params <- list({r_params})
          # knit() defaults to error=TRUE (errors are printed into the output);
          # rmarkdown::render() uses error=FALSE.  Match render().
          knitr::opts_chunk$set(error = FALSE)
          knitr::knit('{SCRIPTS / "genome_report.Rmd"}', output = 'report.md',
                      envir = env, quiet = TRUE)
          readr::write_csv(env$vcr_summ, 'vcr_summ.csv')
          readr::write_csv(env$problem_areas, 'problem_areas.csv')
          readr::write_csv(env$gene_start_stops, 'genes.csv')
          readr::write_csv(env$me_hits, 'me_hits.csv')
          readr::write_csv(env$mapped_reads_sum, 'xmap.csv')
        }})
        """)
    write(td / "run.R", script, dedent=False)
    proc = subprocess.run(["Rscript", "--vanilla", "run.R"], cwd=td,
                          capture_output=True, text=True, timeout=300)
    return proc.returncode, proc.stdout + proc.stderr


def _csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


_RESULT = {}


def _knitted():
    """Knit once and share the outputs between the tests in this module."""
    _need_r()
    if "td" not in _RESULT:
        import tempfile
        td = __import__("pathlib").Path(tempfile.mkdtemp(prefix="callconf_report_"))
        import atexit
        atexit.register(shutil.rmtree, td, True)
        rc, log = _knit(td, _inputs(td))
        _RESULT.update(td=td, rc=rc, log=log)
    if _RESULT["rc"] != 0:
        raise SetupError(f"knitting genome_report.Rmd failed:\n{_RESULT['log'][-3000:]}")
    return _RESULT["td"]


def test_report_knits_and_contains_all_sections():
    td = _knitted()
    md = (td / "report.md").read_text()
    for needle in ("Detection summary for reference: refA",
                   "Duplicated genes at 90pct ANI or higher: refA",
                   "Hits against supplied ME databases: refA"):
        assert needle in md, f"missing from report: {needle!r}"
    assert list((td / "figure").glob("*.png")), "no plots were produced"


def test_report_windowed_recall_and_precision():
    td = _knitted()
    got = {int(float(r["window"])): r for r in _csv(td / "vcr_summ.csv")}
    assert set(got) == {0, 1, 3}
    assert float(got[0]["recall_pct"]) == 100
    assert float(got[1]["recall_pct"]) == 25
    assert float(got[3]["recall_pct"]) == 100 and float(got[3]["precision_pct"]) == 50
    assert int(float(got[1]["n"])) == 4


def test_report_keeps_only_best_significant_db_hit_per_gene():
    """Filter is evalue < 1e-10 and alnlen >= 200, then best evalue per query."""
    td = _knitted()
    hits = _csv(td / "me_hits.csv")
    assert sorted((h["query"], h["target"]) for h in hits) == [
        ("refA_chr_2", "mobileOG_000777|Tn916"),
        ("refA_chr_3", "mobileOG_000123|IS256"),
        ("refA_p1_1", "mobileOG_000888|IS3")], hits


def test_report_flags_duplicated_gene_region():
    td = _knitted()
    areas = _csv(td / "problem_areas.csv")
    dup = [a for a in areas if a["label"] == "Has Close Homolog"]
    assert [(a["chrom"], a["start"], a["stop"]) for a in dup] == [("refA_chr", "1201", "1800")]


def test_report_cross_mapping_summary():
    """r30–r39 map 2 kb from their origin: windows 1–1.9 kb should show them."""
    td = _knitted()
    rows = _csv(td / "xmap.csv")
    same = {int(float(r["window"])): float(r["values"]) for r in rows if r["label"] == "Same genome"}
    assert same[0] == 0, same
    assert same[1] > 0, same


def test_report_gene_ids_match_prodigal_ids_on_every_contig():
    """query must be <contig>_<gene index from ID=contig_gene>, which is how
    Prodigal names genes in the .faa that MMseqs searches."""
    td = _knitted()
    genes = sorted(g["query"] for g in _csv(td / "genes.csv"))
    assert genes == ["refA_chr_1", "refA_chr_2", "refA_chr_3", "refA_p1_1"], genes


def _db_regions(td):
    return {a["target"].split("|")[1]: (a["chrom"], int(float(a["start"])), int(float(a["stop"])))
            for a in _csv(td / "me_hits.csv")}


def test_report_db_hit_forward_strand_coordinates():
    """Gene 3001-3900 (+), protein hit aa 11-260 → nt 3001+30 .. 3001+780-1."""
    assert _db_regions(_knitted())["IS256"] == ("refA_chr", 3031, 3780)


def test_report_db_hit_reverse_strand_coordinates():
    """Gene 1201-1800 (-), protein hit aa 11-190 → counted back from the gene
    end: 1800-570+1 .. 1800-30."""
    assert _db_regions(_knitted())["Tn916"] == ("refA_chr", 1231, 1770)


def test_report_db_hit_on_second_contig():
    """Gene refA_p1_1 is 11-610 (+); hit aa 1-100 → nt 11..310 on refA_p1."""
    assert _db_regions(_knitted())["IS3"] == ("refA_p1", 11, 310)


def test_report_db_hit_regions_lie_within_their_genes():
    td = _knitted()
    genes = {g["query"]: (int(float(g["gstart"])), int(float(g["gend"])))
             for g in _csv(td / "genes.csv")}
    for a in _csv(td / "me_hits.csv"):
        gstart, gend = genes[a["query"]]
        assert gstart <= int(float(a["start"])) <= int(float(a["stop"])) <= gend, \
            (a, (gstart, gend))
