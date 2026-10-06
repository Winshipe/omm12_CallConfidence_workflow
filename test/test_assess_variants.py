"""
Tests for workflow/scripts/assess_variants.py

The script is a straight pandas pipeline with no functions, so every test runs
the *real* script end-to-end through the snakemake stub (the old suite tested
a hand-copied replica of the logic, which could silently drift).
"""

from harness import VCF_HEADER, read_tsv, requires, run_script, tmpdir, write, write_vcf

GT_HEADER = "seq_id\tposition\tref_base\talt_base\n"


def _truth(path, rows):
    """rows: (seq_id, pos, ref, alt) — written exactly like mutate_reference.py."""
    return write(path, GT_HEADER + "".join(f"{s}\t{p}\t{r}\t{a}\n" for s, p, r, a in rows),
                 dedent=False)


def _assess(td, truth_rows, calls, min_quality=20, replicate="rep1", truth_files=None):
    """Run the script; return {(chrom, pos, ref, alt): truthiness} and the table."""
    gts = truth_files if truth_files is not None else [_truth(td / "gt.tsv", truth_rows)]
    vcf = write_vcf(td / "calls.vcf", calls)
    out = td / "assessment.tsv"
    run_script("assess_variants.py",
               input={"vcf": vcf, "ground_truth": gts},
               output={"tsv": out},
               params={"min_quality": min_quality, "replicate": replicate})
    header, rows = read_tsv(out)
    idx = {c: i for i, c in enumerate(header)}
    labels = {(r[idx["CHROM"]], int(r[idx["POS"]]), r[idx["REF"]], r[idx["ALT"]]):
              r[idx["Truthiness"]] for r in rows}
    return labels, header, rows


def _counts(labels):
    out = {"TP": 0, "FP": 0, "FN": 0}
    for v in labels.values():
        out[v] += 1
    return out


@requires("pandas")
def test_output_columns():
    with tmpdir() as td:
        _, header, _ = _assess(td, [("c1", 5, "A", "G")], [("c1", 5, "A", "G", 50)])
        assert header == ["CHROM", "POS", "REF", "ALT", "Truthiness", "replicate"]


@requires("pandas")
def test_true_positive():
    with tmpdir() as td:
        labels, _, _ = _assess(td, [("c1", 100, "A", "T")], [("c1", 100, "A", "T", 50)])
        assert labels == {("c1", 100, "A", "T"): "TP"}


@requires("pandas")
def test_false_negative_not_called():
    with tmpdir() as td:
        labels, _, _ = _assess(td, [("c1", 300, "G", "A")], [("c1", 7, "C", "T", 99)])
        assert labels[("c1", 300, "G", "A")] == "FN"
        assert labels[("c1", 7, "C", "T")] == "FP"


@requires("pandas")
def test_quality_threshold_is_inclusive_and_filters_below():
    with tmpdir() as td:
        labels, _, _ = _assess(
            td,
            [("c1", 10, "A", "C"), ("c1", 20, "A", "C"), ("c1", 30, "A", "C")],
            [("c1", 10, "A", "C", 20), ("c1", 20, "A", "C", 19.99), ("c1", 30, "A", "C", 20.01)],
            min_quality=20)
        assert labels == {("c1", 10, "A", "C"): "TP",
                          ("c1", 20, "A", "C"): "FN",
                          ("c1", 30, "A", "C"): "TP"}


@requires("pandas")
def test_low_quality_false_call_is_dropped_not_fp():
    with tmpdir() as td:
        labels, _, _ = _assess(td, [], [("c1", 5, "A", "G", 3)])
        assert labels == {}


@requires("pandas")
def test_alt_allele_must_match():
    with tmpdir() as td:
        labels, _, _ = _assess(td, [("c1", 600, "A", "T")], [("c1", 600, "A", "G", 50)])
        assert labels == {("c1", 600, "A", "T"): "FN", ("c1", 600, "A", "G"): "FP"}


@requires("pandas")
def test_chromosome_must_match():
    with tmpdir() as td:
        labels, _, _ = _assess(td, [("c1", 5, "A", "T")], [("c2", 5, "A", "T", 50)])
        assert _counts(labels) == {"TP": 0, "FP": 1, "FN": 1}


@requires("pandas")
def test_mixed_realistic_calls():
    truth = [("c1", p, "A", "G") for p in range(10, 200, 10)]          # 19 mutations
    calls = ([("c1", p, "A", "G", 60) for p in range(10, 100, 10)]     # 9 TP
             + [("c1", p, "A", "G", 5) for p in range(100, 150, 10)]   # 5 low-qual → FN
             + [("c1", 1001, "C", "T", 45), ("c1", 1002, "G", "A", 45)])  # 2 FP
    with tmpdir() as td:
        labels, _, rows = _assess(td, truth, calls)
        assert _counts(labels) == {"TP": 9, "FP": 2, "FN": 10}
        assert len(rows) == 21, "one output row per distinct variant"


@requires("pandas")
def test_replicate_label_on_every_row():
    with tmpdir() as td:
        _, header, rows = _assess(td, [("c1", 1, "A", "T"), ("c1", 2, "A", "T")],
                                  [("c1", 1, "A", "T", 50), ("c1", 9, "A", "T", 50)],
                                  replicate="rep7")
        col = header.index("replicate")
        assert {r[col] for r in rows} == {"rep7"} and len(rows) == 3


@requires("pandas")
def test_multiple_ground_truth_files_are_concatenated():
    """Scenarios mixing several references pass one TSV per reference."""
    with tmpdir() as td:
        gts = [_truth(td / "a.tsv", [("refA", 5, "A", "G")]),
               _truth(td / "b.tsv", [("refB", 5, "C", "T"), ("refB", 9, "G", "A")])]
        labels, _, _ = _assess(td, None, [("refA", 5, "A", "G", 50), ("refB", 9, "G", "A", 50)],
                               truth_files=gts)
        assert labels == {("refA", 5, "A", "G"): "TP", ("refB", 5, "C", "T"): "FN",
                          ("refB", 9, "G", "A"): "TP"}


@requires("pandas")
def test_header_only_vcf_gives_all_fn():
    with tmpdir() as td:
        labels, _, _ = _assess(td, [("c1", 10, "A", "G"), ("c1", 20, "T", "C")], [])
        assert _counts(labels) == {"TP": 0, "FP": 0, "FN": 2}


@requires("pandas")
def test_ground_truth_with_no_mutations_gives_all_fp():
    """A low substitution rate on a short contig can yield an empty truth TSV."""
    with tmpdir() as td:
        labels, _, _ = _assess(td, [], [("c1", 1, "A", "T", 30), ("c1", 2, "C", "G", 40)])
        assert _counts(labels) == {"TP": 0, "FP": 2, "FN": 0}


@requires("pandas")
def test_both_empty_gives_empty_table():
    with tmpdir() as td:
        labels, header, rows = _assess(td, [], [])
        assert rows == [] and "Truthiness" in header


@requires("pandas")
def test_vcf_info_and_sample_columns_are_ignored():
    """GATK writes a 10th sample column; extra INFO content must not matter."""
    with tmpdir() as td:
        vcf = write(td / "calls.vcf",
                    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                    "c1\t5\trs1\tA\tG\t812.64\tPASS\tAC=2;AF=1.00;AN=2;DP=19\tGT:AD:DP:GQ:PL\t1/1:0,19:19:57:826,57,0\n",
                    dedent=False)
        gt = _truth(td / "gt.tsv", [("c1", 5, "A", "G")])
        out = td / "o.tsv"
        run_script("assess_variants.py", input={"vcf": vcf, "ground_truth": [gt]},
                   output={"tsv": out}, params={"min_quality": 20, "replicate": "rep1"})
        header, rows = read_tsv(out)
        assert rows == [["c1", "5", "A", "G", "TP", "rep1"]]


@requires("pandas")
def test_log_records_counts():
    with tmpdir() as td:
        gt = _truth(td / "gt.tsv", [("c1", 5, "A", "G"), ("c1", 6, "A", "G")])
        vcf = write_vcf(td / "calls.vcf", [("c1", 5, "A", "G", 50)])
        log = td / "a.log"
        run_script("assess_variants.py", input={"vcf": vcf, "ground_truth": [gt]},
                   output={"tsv": td / "o.tsv"}, params={"min_quality": 20, "replicate": "r"},
                   log=[log])
        text = log.read_text()
        assert "Loaded 2 expected mutations" in text and "Parsed 1 SNP calls" in text


@requires("pandas")
def test_scenario_without_mutated_references():
    with tmpdir() as td:
        labels, _, _ = _assess(td, None, [("c1", 5, "A", "G", 50)], truth_files=[])
        assert _counts(labels) == {"TP": 0, "FP": 1, "FN": 0}


@requires("pandas")
def test_multiallelic_record_is_split_per_alt():
    """'G,T' is split into one call per ALT allele, so G matches the truth and
    T is a false positive."""
    with tmpdir() as td:
        labels, _, _ = _assess(td, [("c1", 5, "A", "G")], [("c1", 5, "A", "G,T", 50)])
        assert labels == {("c1", 5, "A", "G"): "TP", ("c1", 5, "A", "T"): "FP"}


@requires("pandas")
def test_gvcf_reference_blocks_and_non_ref_alleles_are_ignored():
    """HaplotypeCaller -ERC GVCF: <NON_REF>-only blocks have QUAL '.', and
    variant records carry an extra <NON_REF> ALT."""
    with tmpdir() as td:
        vcf = write(td / "calls.g.vcf", VCF_HEADER +
                    "c1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=4\tGT:DP\t0/0:30\n"
                    "c1\t5\t.\tA\tG,<NON_REF>\t300.6\t.\tDP=30\tGT:AD\t0/1:15,15,0\n"
                    "c1\t6\t.\tC\t<NON_REF>\t.\t.\tEND=99\tGT:DP\t0/0:30\n",
                    dedent=False)
        _, rows = _run_vcf(td, vcf, [("c1", 5, "A", "G"), ("c1", 50, "C", "T")])
        assert sorted(rows) == [["c1", "5", "A", "G", "TP", "rep1"],
                                ["c1", "50", "C", "T", "FN", "rep1"]]


# ── Mutect2 / FilterMutectCalls output ───────────────────────────────────────

MUTECT2_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##FILTER=<ID=weak_evidence,Description="Mutation does not meet likelihood threshold">\n'
    "##source=FilterMutectCalls\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmix_rep1\n"
)


def _mutect_vcf(td, records):
    """records: (chrom, pos, ref, alt, filter) — QUAL is always '.' in Mutect2."""
    return write(td / "mutect2.filtered.vcf", MUTECT2_HEADER + "".join(
        f"{c}\t{p}\t.\t{r}\t{a}\t.\t{f}\tDP=40;TLOD=25.3\tGT:AD:AF\t0/1:20,20:0.5\n"
        for c, p, r, a, f in records), dedent=False)


def _run_vcf(td, vcf, truth_rows, min_quality=20):
    gt = _truth(td / "gt.tsv", truth_rows)
    out = td / "o.tsv"
    run_script("assess_variants.py", input={"vcf": vcf, "ground_truth": [gt]},
               output={"tsv": out}, params={"min_quality": min_quality, "replicate": "rep1"})
    return read_tsv(out)


@requires("pandas")
def test_mutect2_pass_calls_count_despite_missing_qual():
    with tmpdir() as td:
        vcf = _mutect_vcf(td, [("c1", 5, "A", "G", "PASS"), ("c1", 9, "C", "T", "PASS")])
        _, rows = _run_vcf(td, vcf, [("c1", 5, "A", "G")])
        assert sorted((r[1], r[4]) for r in rows) == [("5", "TP"), ("9", "FP")]


@requires("pandas")
def test_mutect2_filtered_calls_are_ignored():
    """Calls that fail FilterMutectCalls (e.g. weak_evidence) are not calls:
    a filtered true mutation is a FN and a filtered artefact is not a FP."""
    with tmpdir() as td:
        vcf = _mutect_vcf(td, [("c1", 5, "A", "G", "weak_evidence"),
                               ("c1", 9, "C", "T", "strand_bias;weak_evidence"),
                               ("c1", 12, "G", "A", "PASS")])
        _, rows = _run_vcf(td, vcf, [("c1", 5, "A", "G"), ("c1", 12, "G", "A")])
        assert sorted((r[1], r[4]) for r in rows) == [("12", "TP"), ("5", "FN")]


@requires("pandas")
def test_mutect2_multiallelic_pass_record():
    with tmpdir() as td:
        vcf = _mutect_vcf(td, [("c1", 5, "A", "G,T", "PASS")])
        _, rows = _run_vcf(td, vcf, [("c1", 5, "A", "T")])
        assert sorted((r[3], r[4]) for r in rows) == [("G", "FP"), ("T", "TP")]


@requires("pandas")
def test_haplotypecaller_filter_column_is_respected_too():
    """A non-PASS FILTER from any caller (e.g. after VariantFiltration) drops
    the call even when its QUAL is high."""
    with tmpdir() as td:
        vcf = write(td / "c.vcf", VCF_HEADER +
                    "c1\t5\t.\tA\tG\t500\tLowQD\t.\tGT\t0/1\n"
                    "c1\t6\t.\tA\tG\t500\tPASS\t.\tGT\t0/1\n", dedent=False)
        _, rows = _run_vcf(td, vcf, [("c1", 5, "A", "G")])
        assert sorted((r[1], r[4]) for r in rows) == [("5", "FN"), ("6", "FP")]
