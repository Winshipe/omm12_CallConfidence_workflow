"""
Integration tests: execute the workflow's own Python steps in dependency
order, wired exactly as the rules wire them.

The DAG is built by smk.py from the real Snakefile, and every job whose rule
uses ``script:`` is run with that job's real input/output/params/log names.
This checks the rule ↔ script contract (e.g. that assess_variants.py's
``snakemake.input.vcf`` really is declared by the rule), which unit tests of
either side cannot.  Steps that need external tools (ART, BWA, GATK, Prodigal,
MMseqs2) are replaced by small synthetic outputs written to the paths the DAG
expects.
"""

import copy
import csv
import random
from collections import defaultdict
from pathlib import Path

import smk
from harness import (CONFIG, ROOT, chdir, requires, run_job, tmpdir, write,
                     write_fasta, write_vcf)


def _config(td, scenarios, replicates=2):
    rng = random.Random(0)
    seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    refa = write_fasta(td / "refs" / "refA.fa", [("refA_chr", seq(6000))])
    refb = write_fasta(td / "refs" / "refB.fa", [("refB_chr", seq(4000)), ("refB_p1", seq(1500))])
    write(td / "db" / "mobileOG", "")
    cfg = copy.deepcopy(smk.load_config(CONFIG))
    cfg["references"] = {"refA": str(refa), "refB": str(refb)}
    cfg["replicates"] = replicates
    cfg["annotation"]["mmseqs_databases"] = {"mobileOG": str(td / "db" / "mobileOG")}
    cfg["mutation"]["substitution_rate"] = 0.01
    cfg["scenarios"] = scenarios
    return cfg


def _truth(path):
    with open(path) as fh:
        return [(r["seq_id"], int(r["position"]), r["ref_base"], r["alt_base"])
                for r in csv.DictReader(fh, delimiter="\t")]


class FakeCaller:
    """Stands in for BWA + the selected GATK caller: calls ~80 % of true variants with high QUAL,
    ~10 % with QUAL below threshold, misses the rest, adds 2 false positives."""

    def __init__(self, seed=0):
        self.rng = random.Random(seed)
        self.expected = defaultdict(lambda: {"TP": 0, "FN": 0, "FP": 0})

    def write_vcf(self, path, truth_rows, key):
        calls = []
        for chrom, pos, ref, alt in truth_rows:
            u = self.rng.random()
            if u < 0.8:
                calls.append((chrom, pos, ref, alt, 60.0))
                self.expected[key]["TP"] += 1
            else:
                if u < 0.9:
                    calls.append((chrom, pos, ref, alt, 5.0))
                self.expected[key]["FN"] += 1
        chrom = truth_rows[0][0] if truth_rows else "refA_chr"
        for fp_pos in (999_991, 999_992):
            calls.append((chrom, fp_pos, "A", "C", 45.0))
            self.expected[key]["FP"] += 1
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        write_vcf(path, sorted(calls, key=lambda c: (c[0], c[1])))


def _run_python_pipeline(td, cfg):
    wf = smk.load_workflow(ROOT, cfg)
    jobs, _, errors = smk.resolve_dag(wf, cwd=td)
    assert not errors, errors
    by_rule = defaultdict(list)
    for j in jobs:
        by_rule[j.rule.name].append(j)
    caller = FakeCaller()
    with chdir(td):
        for job in by_rule["mutate_reference"]:
            run_job(job)
        # Stand-in for the variant caller: write the VCF that assess_variants reads
        for job in by_rule["select_variant_calls"]:
            truth = []
            scen = cfg["scenarios"][job.wildcards["scenario"]]
            for c in scen:
                if c["mutated_fraction"] > 0:
                    truth += _truth(f"results/mutated/{c['ref_id']}/{job.wildcards['replicate']}"
                                    f"/{c['ref_id']}.mutations.tsv")
            caller.write_vcf(job.output.vcf, truth, (job.wildcards["scenario"],
                                                     job.wildcards["replicate"]))
        for job in by_rule["assess_variants"]:
            run_job(job)
        for job in by_rule["aggregate_replicates"]:
            wf.run_block("aggregate_replicates", input=job.input, output=job.output, log=job.log)
    return wf, by_rule, caller


@requires("pandas")
def test_mutate_assess_aggregate_through_real_rule_wiring():
    scenarios = {
        "mix": [{"ref_id": "refA", "mutated_fraction": 0.5, "abundance": 1.0},
                {"ref_id": "refB", "mutated_fraction": 0.3, "abundance": 1.0}],
        "only_a": [{"ref_id": "refA", "mutated_fraction": 1.0, "abundance": 1.0}],
        "a_with_unmutated_b": [{"ref_id": "refA", "mutated_fraction": 0.5, "abundance": 1.0},
                               {"ref_id": "refB", "mutated_fraction": 0.0, "abundance": 1.0}],
    }
    with tmpdir() as td:
        cfg = _config(td, scenarios)
        _, by_rule, caller = _run_python_pipeline(td, cfg)

        # mutate: every output + log exists, mutated FASTA explained by the TSV
        for job in by_rule["mutate_reference"]:
            for p in list(job.output) + list(job.log):
                assert (td / p).exists(), f"{job.rule.name} did not create {p}"
            truth = _truth(td / job.output.mutations_tsv)
            assert truth, "rate 0.01 on ≥5 kb should yield mutations"

        # assess: per-replicate counts equal what the fake caller did
        for job in by_rule["assess_variants"]:
            key = (job.wildcards["scenario"], job.wildcards["replicate"])
            with open(td / job.output.tsv) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            got = {k: sum(r["Truthiness"] == k for r in rows) for k in ("TP", "FN", "FP")}
            assert got == caller.expected[key], (key, got, caller.expected[key])
            assert {r["replicate"] for r in rows} == {job.wildcards["replicate"]}

        # scenario without mutated refB: no refB truth may leak in
        job = next(j for j in by_rule["assess_variants"]
                   if j.wildcards["scenario"] == "a_with_unmutated_b")
        assert all("refB" not in p for p in job.input.ground_truth)

        # aggregate: one file per scenario holding every replicate's rows
        for job in by_rule["aggregate_replicates"]:
            scen = job.wildcards["scenario"]
            with open(td / job.output.tsv) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            want = sum(sum(caller.expected[(scen, f"rep{i}")].values())
                       for i in range(1, cfg["replicates"] + 1))
            assert len(rows) == want, (scen, len(rows), want)
            assert {r["replicate"] for r in rows} == {"rep1", "rep2"}


@requires("Bio")
def test_find_homologs_through_real_rule_wiring():
    """Feed fake Prodigal/MMseqs outputs to the find_homologs jobs."""
    scenarios = {"mix": [{"ref_id": "refA", "mutated_fraction": 0.5, "abundance": 1.0},
                         {"ref_id": "refB", "mutated_fraction": 0.5, "abundance": 1.0}]}
    with tmpdir() as td:
        cfg = _config(td, scenarios, replicates=1)
        jobs, _, errors = smk.resolve_dag(smk.load_workflow(ROOT, cfg), cwd=td)
        assert not errors, errors
        fh_jobs = [j for j in jobs if j.rule.name == "find_homologs"]
        assert len(fh_jobs) == 2
        with chdir(td):
            write("results/annotation/clusters.txt",
                  "refA_chr_1\trefA_chr_1\nrefA_chr_1\trefB_chr_2\nrefA_chr_2\trefA_chr_2\n"
                  "refB_p1_1\trefB_p1_1\nrefB_p1_1\trefB_chr_1\n")
            write("results/annotation/refA/refA.fna",
                  ">refA_chr_1 # 1 # 300 # 1 # ID=1_1\nATG\n>refA_chr_2 # 400 # 700 # -1 # ID=1_2\nATG\n")
            write("results/annotation/refB/refB.fna",
                  ">refB_chr_1 # 5 # 305 # -1 # ID=1_1\nATG\n>refB_chr_2 # 900 # 1200 # 1 # ID=1_2\nATG\n"
                  ">refB_p1_1 # 10 # 400 # 1 # ID=2_1\nATG\n")
            for job in fh_jobs:
                run_job(job)
            a = Path("results/annotation/refA/refA.challenging.tsv").read_text().splitlines()[1:]
            b = Path("results/annotation/refB/refB.challenging.tsv").read_text().splitlines()[1:]
        assert [l.split("\t")[3] for l in a] == ["refA_chr_1"]
        assert [l.split("\t")[3] for l in b] == ["refB_chr_1", "refB_chr_2", "refB_p1_1"]
        assert [l.split("\t")[0] for l in b] == ["refB_chr", "refB_chr", "refB_p1"]
