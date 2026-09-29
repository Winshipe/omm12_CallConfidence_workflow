"""
Tests for the Python that lives *inside* the Snakemake rule files:
input/params helper functions and the aggregate_replicates ``run:`` block.

The workflow is loaded with smk.py (no Snakemake needed) using a synthetic
config that exercises more branches than config/config.yaml: two references,
unequal abundances, a zero-mutation control, and an empirical ART profile.
"""

import copy
import csv
from pathlib import Path

import smk
from harness import ROOT, chdir, read_tsv, tmpdir, write

BASE_CONFIG = {
    "references": {"refA": "resources/refA.fa", "refB": "resources/refB.fa"},
    "replicates": 2,
    "mutation": {"model": "tamura_nei", "substitution_rate": 0.001, "kappa": 2.0,
                 "gc_freq": 0.5, "base_seed": 42},
    "annotation": {"mmseqs_databases": {"mobileOG": "db/mobileOG", "phrog": "db/phrog"},
                   "split_mem_limit": "164G", "max_memory": "196G", "threads": 8},
    "simulation": {"platform": "MSv3", "read_length": 150, "mean_fragment_length": 400,
                   "std_fragment_length": 100, "coverage": 50,
                   "empirical_reads_R1": None, "empirical_reads_R2": None, "threads": 4},
    "scenarios": {
        "equal_mix": [{"ref_id": "refA", "mutated_fraction": 0.5, "abundance": 1.0},
                      {"ref_id": "refB", "mutated_fraction": 0.2, "abundance": 1.0}],
        "control": [{"ref_id": "refA", "mutated_fraction": 0.0, "abundance": 1.0}],
        "skewed": [{"ref_id": "refA", "mutated_fraction": 0.5, "abundance": 3.0},
                   {"ref_id": "refB", "mutated_fraction": 0.0, "abundance": 1.0}],
        "very_skewed": [{"ref_id": "refA", "mutated_fraction": 0.2, "abundance": 3.0},
                        {"ref_id": "refB", "mutated_fraction": 0.0, "abundance": 1.0}],
    },
    "variant_calling": {"threads": 8, "mem_mb": 16000, "min_base_quality_score": 20,
                        "gatk_extra_flags": ""},
    "assessment": {"min_quality": 20},
    "report": {"genome_length": 0, "window_size": 2000},
}


def _wf(**overrides):
    cfg = copy.deepcopy(BASE_CONFIG)
    for key, value in overrides.items():
        cfg[key] = value
    return smk.load_workflow(ROOT, cfg)


def _param(wf, rule, name):
    return wf.rules[rule].kwargs["params"][name]


def W(**kw):
    return smk.Wildcards(kw)


# ── mutate.smk / simulate_reads.smk seeds ────────────────────────────────────

def test_mutate_seed_is_distinct_per_replicate():
    seed = _param(_wf(), "mutate_reference", "seed")
    seeds = [seed(W(ref_id="refA", replicate=f"rep{i}")) for i in range(1, 11)]
    assert len(set(seeds)) == 10 and seeds[0] == 42 + 10


def test_simulation_seed_is_distinct_per_replicate():
    seed = _param(_wf(), "simulate_reads", "seed")
    seeds = [seed(W(ref_id="refA", replicate=f"rep{i}", mutated="mutated")) for i in range(1, 11)]
    assert len(set(seeds)) == 10


def test_seed_helpers_are_bound_per_rule_despite_shared_name():
    """mutate.smk and simulate_reads.smk both define replicate_seed(); each rule
    must keep the version defined in its own file (bound at definition time)."""
    wf = _wf()
    wc = W(ref_id="refA", replicate="rep3", mutated="mutated")
    assert _param(wf, "mutate_reference", "seed")(wc) == 42 + 30
    assert _param(wf, "simulate_reads", "seed")(wc) == 42 + 3


# ── simulate_reads.smk ───────────────────────────────────────────────────────

def test_art_ref_fasta_picks_mutated_or_original():
    fn = _wf().ns["art_ref_fasta"]
    assert fn(W(ref_id="refA", replicate="rep2", mutated="mutated")) == \
        "results/mutated/refA/rep2/refA.mutated.fasta"
    assert fn(W(ref_id="refA", replicate="rep2", mutated="unmutated")) == "resources/refA.fa"


def test_art_profile_flags_builtin_and_empirical():
    assert _wf().ns["art_profile_flags"](W()) == "-ss MSv3"
    sim = dict(BASE_CONFIG["simulation"], empirical_reads_R1="r1.fq", empirical_reads_R2="r2.fq")
    flags = _wf(simulation=sim).ns["art_profile_flags"](W())
    assert "--qprof1 results/art_profile/empirical_R1.txt" in flags
    assert "--qprof2 results/art_profile/empirical_R2.txt" in flags


# ── blend_reads.smk ──────────────────────────────────────────────────────────

def test_scenario_input_reads_keys_and_paths():
    fn = _wf().ns["scenario_input_reads"]
    got = fn(W(scenario="equal_mix", replicate="rep1"))
    assert len(got) == 8
    assert got["contrib_1_mutated_R2"] == "results/simulated/refB/rep1/mutated/refB_R2.fastq.gz"
    assert got["contrib_0_unmutated_R1"] == "results/simulated/refA/rep1/unmutated/refA_R1.fastq.gz"


def _seqtk_commands(wf, scenario, replicate="rep1"):
    """Run the generator in a scratch cwd; return ({file: fraction}, script, cwd log)."""
    fn = wf.ns["prepare_seqtk_sampling_script"]
    with tmpdir() as td, chdir(td):
        script = fn(W(scenario=scenario, replicate=replicate))
        log = (td / "logs" / "blend" / scenario / f"{replicate}.log").read_text()
    fractions = {}
    for line in script.splitlines():
        if line.startswith("seqtk sample"):
            parts = line.split()
            assert parts[2] == "-s", line
            fractions[parts[4]] = float(parts[5])
    return fractions, script, log


def test_seqtk_script_samples_every_input_and_gzips_outputs():
    fractions, script, log = _seqtk_commands(_wf(), "equal_mix")
    assert len(fractions) == 8, "2 refs × {mutated, unmutated} × {R1, R2}"
    assert script == log, "the generated script is also written to the log"
    assert script.rstrip().endswith("equal_mix_R2.fastq;")
    assert "gzip results/blended/equal_mix/rep1/equal_mix_R1.fastq;" in script


def test_seqtk_mates_use_same_seed_and_fraction():
    """seqtk keeps pairs together only if R1 and R2 get identical -s and fraction."""
    fractions, script, _ = _seqtk_commands(_wf(), "equal_mix")
    seeds = {line.split()[3] for line in script.splitlines() if line.startswith("seqtk")}
    assert len(seeds) == 1
    for path, frac in fractions.items():
        if "_R1.fastq.gz" in path:
            assert fractions[path.replace("_R1.fastq.gz", "_R2.fastq.gz")] == frac


def test_seqtk_mutated_fraction_split():
    fractions, _, _ = _seqtk_commands(_wf(), "equal_mix")
    get = lambda ref, kind: fractions[f"results/simulated/{ref}/rep1/{kind}/{ref}_R1.fastq.gz"]
    assert abs(get("refA", "mutated") / (get("refA", "mutated") + get("refA", "unmutated")) - 0.5) < 1e-9
    assert abs(get("refB", "mutated") / (get("refB", "mutated") + get("refB", "unmutated")) - 0.2) < 1e-9


def test_seqtk_control_scenario_takes_no_mutated_reads():
    fractions, _, _ = _seqtk_commands(_wf(), "control")
    assert all(f == 0 for p, f in fractions.items() if "/mutated/" in p)
    assert all(f == 1 for p, f in fractions.items() if "/unmutated/" in p)


def test_seqtk_outputs_are_per_replicate():
    _, script, _ = _seqtk_commands(_wf(), "equal_mix", replicate="rep2")
    assert "/rep1/" not in script and "results/blended/equal_mix/rep2/" in script


def test_seqtk_fractions_are_valid_for_unequal_abundances():
    # very_skewed: refA abundance 3, mutated_fraction 0.2 (the old formula gave 1.2)
    fractions, _, _ = _seqtk_commands(_wf(), "very_skewed")
    bad = {p: f for p, f in fractions.items() if not 0 <= f <= 1}
    assert not bad, f"fractions outside [0, 1]: {bad}"


def test_seqtk_fractions_preserve_abundance_ratio():
    """Each reference is simulated at the same coverage, so the reads taken
    per reference (mutated + unmutated) must follow the abundance ratio, and
    the most abundant reference uses all of its reads."""
    for scenario, want in (("very_skewed", {"refA": 1.0, "refB": 1 / 3}),
                           ("skewed", {"refA": 1.0, "refB": 1 / 3}),
                           ("equal_mix", {"refA": 1.0, "refB": 1.0})):
        fractions, _, _ = _seqtk_commands(_wf(), scenario)
        total = {ref: sum(f for p, f in fractions.items()
                          if f"/{ref}/" in p and p.endswith("_R1.fastq.gz"))
                 for ref in want}
        for ref, w in want.items():
            assert abs(total[ref] - w) < 1e-9, (scenario, total)


# ── variant_calling.smk / assess.smk / report.smk helpers ────────────────────

def test_scenario_references_dedupes_in_order():
    scen = dict(BASE_CONFIG["scenarios"], dup=[
        {"ref_id": "refB", "mutated_fraction": 0.1, "abundance": 1},
        {"ref_id": "refA", "mutated_fraction": 0.1, "abundance": 1},
        {"ref_id": "refB", "mutated_fraction": 0.3, "abundance": 1}])
    fn = _wf(scenarios=scen).ns["scenario_references"]
    assert fn(W(scenario="dup")) == ["resources/refB.fa", "resources/refA.fa"]


def test_scenario_ground_truth_tsvs_only_mutated_refs():
    wf = _wf()
    fn = wf.rules["assess_variants"].kwargs["input"]["ground_truth"]
    assert fn(W(scenario="equal_mix", replicate="rep2")) == [
        "results/mutated/refA/rep2/refA.mutations.tsv",
        "results/mutated/refB/rep2/refB.mutations.tsv"]
    assert fn(W(scenario="skewed", replicate="rep1")) == [
        "results/mutated/refA/rep1/refA.mutations.tsv"]


def test_all_replicate_tsvs():
    fn = _wf().ns["all_replicate_tsvs"]
    assert fn(W(scenario="control")) == ["results/assessment/control/rep1_assessment.tsv",
                                         "results/assessment/control/rep2_assessment.tsv"]


def test_scenarios_for_ref_skips_unmutated_contributions():
    fn = _wf().ns["scenarios_for_ref"]
    assert fn("refA") == ["equal_mix", "skewed", "very_skewed"]
    assert fn("refB") == ["equal_mix"]


def test_report_annotation_hits_dict_and_list_configs():
    wf = _wf()
    assert wf.ns["report_annotation_hits"](W(ref_id="refA")) == [
        "results/annotation/refA/mobileOG_hits", "results/annotation/refA/phrog_hits"]
    # report.smk also accepts a *list* of database paths.  annotate.smk (which
    # runs first) needs a mapping, so patch the config after the rules are read.
    wf.config["annotation"]["mmseqs_databases"] = ["/dbs/mobileOG", "/dbs/phrog"]
    assert wf.ns["report_annotation_hits"](W(ref_id="refB")) == [
        "results/annotation/refB/mobileOG_hits", "results/annotation/refB/phrog_hits"]
    db_names = wf.rules["generate_report"].kwargs["params"]["db_names"]
    assert db_names(W(ref_id="refB")) == "mobileOG,phrog"


def test_report_mapped_reads_uses_rep1_of_each_scenario():
    fn = _wf().ns["report_mapped_reads"]
    assert fn(W(ref_id="refA")) == ["results/variant_calling/equal_mix/rep1/aligned.sam",
                                    "results/variant_calling/skewed/rep1/aligned.sam",
                                    "results/variant_calling/very_skewed/rep1/aligned.sam"]


def test_bwa_align_read_group_string():
    rg = _param(_wf(), "bwa_align", "rg")(W(scenario="mix", replicate="rep2"))
    assert rg == r"@RG\tID:mix_rep2\tSM:mix\tPL:ILLUMINA\tLB:mix_rep2"


def test_mmseqs_memory_resource_parsed_from_config():
    wf = _wf()
    assert wf.rules["mmseqs_search"].kwargs["resources"]["mem_mb"] == 196 * 1024
    assert wf.rules["mmseqs_search"].kwargs["resources"]["cpus_per_task"] == 8


# ── aggregate_replicates run: block ──────────────────────────────────────────

def _assessment(path, rep, rows):
    return write(path, "CHROM\tPOS\tREF\tALT\tTruthiness\treplicate\n" +
                 "".join(f"{c}\t{p}\t{r}\t{a}\t{t}\t{rep}\n" for c, p, r, a, t in rows),
                 dedent=False)


def _aggregate(td, files):
    out, log = td / "all.tsv", td / "agg.log"
    _wf().run_block("aggregate_replicates", input={"tsvs": [str(f) for f in files]},
                    output={"tsv": str(out)}, log=[str(log)])
    return out, log


def test_aggregate_concatenates_all_rows_with_replicate_from_filename():
    with tmpdir() as td:
        f1 = _assessment(td / "s" / "rep1_assessment.tsv", "rep1",
                         [("c", 1, "A", "G", "TP"), ("c", 2, "A", "G", "FN")])
        f2 = _assessment(td / "s" / "rep2_assessment.tsv", "rep2", [("c", 9, "C", "T", "FP")])
        out, log = _aggregate(td, [f2, f1])
        with open(out) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert [(r["POS"], r["replicate"]) for r in rows] == [("1", "rep1"), ("2", "rep1"),
                                                              ("9", "rep2")]
        assert "Aggregated 3 rows across 2 replicates" in log.read_text()


def test_aggregate_ten_replicates_sorted_order():
    """sorted() is lexicographic: rep10 sorts before rep2.  Rows must still all
    be present with the right labels."""
    with tmpdir() as td:
        files = [_assessment(td / f"rep{i}_assessment.tsv", f"rep{i}", [("c", i, "A", "G", "TP")])
                 for i in range(1, 11)]
        out, _ = _aggregate(td, files)
        with open(out) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert {(r["POS"], r["replicate"]) for r in rows} == {(str(i), f"rep{i}") for i in range(1, 11)}


def test_aggregate_header_has_unique_columns():
    with tmpdir() as td:
        f1 = _assessment(td / "rep1_assessment.tsv", "rep1", [("c", 1, "A", "G", "TP")])
        out, _ = _aggregate(td, [f1])
        header, _ = read_tsv(out)
        assert len(header) == len(set(header)), f"duplicate columns in {header}"


def test_aggregate_all_empty_inputs_still_has_header():
    with tmpdir() as td:
        f1 = _assessment(td / "rep1_assessment.tsv", "rep1", [])
        out, _ = _aggregate(td, [f1])
        assert out.read_text().startswith("CHROM"), repr(out.read_text())
