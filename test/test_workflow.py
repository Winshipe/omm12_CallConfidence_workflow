"""
Workflow-level checks that normally need Snakemake (or a real run) to surface:

  * config/config.yaml is well formed and internally consistent
  * every script / conda env / include referenced by a rule exists, and each
    conda env provides the tools and Python packages its rule uses
  * wildcards are consistent between output, log and input of every rule
  * a dry-run-style DAG can be built from `rule all` for the shipped config
    and for a richer synthetic config (multi-reference, empirical profile)
  * every job's shell command formats without placeholder errors
  * shell blocks have no orphaned continuation lines or undefined variables
  * the report rule only passes parameters that genome_report.Rmd declares
  * every config key is actually read by the workflow

All of it runs on the standard library via smk.py.
"""

import ast
import copy
import re
import sys
from collections import Counter
from pathlib import Path

import smk
from harness import (CONFIG, ENVS, ROOT, RULES, SCRIPTS, WORKFLOW,
                     tmpdir, write, write_fasta)

# ── shared fixtures ──────────────────────────────────────────────────────────

_CACHE = {}


def real_config():
    return copy.deepcopy(smk.load_config(CONFIG))


def real_workflow():
    if "wf" not in _CACHE:
        _CACHE["wf"] = smk.load_workflow(ROOT, real_config())
    return _CACHE["wf"]


def real_dag():
    if "dag" not in _CACHE:
        with tmpdir() as td:
            _CACHE["dag"] = smk.resolve_dag(real_workflow(), cwd=td)
    return _CACHE["dag"]


def _synthetic_config(td):
    """Two small references, three scenarios, empirical ART profile."""
    refa = write_fasta(td / "refs" / "refA.fa", [("refA_chr", "ACGT" * 50)])
    refb = write_fasta(td / "refs" / "refB.fa", [("refB_chr", "GGCC" * 50),
                                                 ("refB_plasmid", "ATAT" * 20)])
    write(td / "db" / "mobileOG", "")
    write(td / "db" / "phrog", "")
    write(td / "emp_R1.fq.gz", "")
    write(td / "emp_R2.fq.gz", "")
    cfg = real_config()
    cfg["references"] = {"refA": str(refa), "refB": str(refb)}
    cfg["replicates"] = 2
    cfg["annotation"]["mmseqs_databases"] = {"mobileOG": str(td / "db" / "mobileOG"),
                                             "phrog": str(td / "db" / "phrog")}
    cfg["simulation"]["empirical_reads_R1"] = str(td / "emp_R1.fq.gz")
    cfg["simulation"]["empirical_reads_R2"] = str(td / "emp_R2.fq.gz")
    cfg["scenarios"] = {
        "mix": [{"ref_id": "refA", "mutated_fraction": 0.5, "abundance": 1.0},
                {"ref_id": "refB", "mutated_fraction": 0.1, "abundance": 1.0}],
        "only_b": [{"ref_id": "refB", "mutated_fraction": 1.0, "abundance": 2.0}],
        "control": [{"ref_id": "refA", "mutated_fraction": 0.0, "abundance": 1.0}],
    }
    return cfg


# ═══════════════════════════════════════════════════════════════════════════════
# YAML subset parser (used when PyYAML is absent)
# ═══════════════════════════════════════════════════════════════════════════════

def test_yaml_subset_parser_features():
    text = '''
    # comment
    a: 1            # trailing comment
    b: "quoted # not a comment"
    c:
      nested: 0.001
      empty: null
      flag: true
      flow: [x, 2, 'y']
      next_line:
        "value on next line"
    items:
      - ref_id: r1
        frac: 0.5
      - ref_id: r2
        frac: 1
    plain_list:
    - one
    - 2
    '''
    import textwrap
    got = smk.load_yaml_subset(textwrap.dedent(text))
    assert got == {
        "a": 1, "b": "quoted # not a comment",
        "c": {"nested": 0.001, "empty": None, "flag": True, "flow": ["x", 2, "y"],
              "next_line": "value on next line"},
        "items": [{"ref_id": "r1", "frac": 0.5}, {"ref_id": "r2", "frac": 1}],
        "plain_list": ["one", 2],
    }, got


def test_yaml_subset_parser_rejects_duplicates_and_unsupported():
    for bad in ("a: 1\na: 2\n", "a: &anchor 1\nb: *anchor\n", "a: |\n  text\n"):
        try:
            smk.load_yaml_subset(bad)
        except smk.YamlSubsetError:
            continue
        raise AssertionError(f"should reject: {bad!r}")


def test_yaml_subset_parser_matches_pyyaml_on_config():
    """When PyYAML is present, the fallback parser must agree with it."""
    try:
        import yaml
    except ImportError:
        from harness import SkipTest
        raise SkipTest("PyYAML not installed (fallback parser is in use)")
    text = CONFIG.read_text()
    assert smk.load_yaml_subset(text) == yaml.safe_load(text)


# ═══════════════════════════════════════════════════════════════════════════════
# Config validation
# ═══════════════════════════════════════════════════════════════════════════════

IUPAC = set("ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.*")


def _num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def validate_config(cfg, root=ROOT):
    """Return a list of human-readable problems (empty when valid)."""
    p = []
    req = lambda cond, msg: None if cond else p.append(msg)

    refs = cfg.get("references")
    req(isinstance(refs, dict) and refs, "references: must be a non-empty mapping")
    refs = refs if isinstance(refs, dict) else {}
    for rid, path in refs.items():
        req(re.fullmatch(r"[A-Za-z0-9_.-]+", str(rid)),
            f"references.{rid}: id must be a single path-safe word")
        full = Path(path) if Path(str(path)).is_absolute() else root / str(path)
        req(full.is_file(), f"references.{rid}: file not found: {path}")

    reps = cfg.get("replicates")
    req(isinstance(reps, int) and not isinstance(reps, bool) and reps >= 1,
        "replicates: must be an integer >= 1")

    mut = cfg.get("mutation") or {}
    req(mut.get("model") in ("jukes_cantor", "tamura_nei"),
        "mutation.model: must be 'jukes_cantor' or 'tamura_nei'")
    rate = mut.get("substitution_rate")
    req(_num(rate) and 0 < rate <= 1, "mutation.substitution_rate: must be in (0, 1]")
    req(_num(mut.get("kappa", 2.0)) and mut.get("kappa", 2.0) > 0, "mutation.kappa: must be > 0")
    gc = mut.get("gc_freq", 0.5)
    req(_num(gc) and 0 < gc < 1, "mutation.gc_freq: must be in (0, 1)")
    req(isinstance(mut.get("base_seed"), int), "mutation.base_seed: must be an integer")

    ann = cfg.get("annotation") or {}
    dbs = ann.get("mmseqs_databases")
    req(isinstance(dbs, dict) and dbs, "annotation.mmseqs_databases: must be a non-empty mapping")
    for name, path in (dbs or {}).items() if isinstance(dbs, dict) else []:
        full = Path(str(path)) if Path(str(path)).is_absolute() else root / str(path)
        req(full.exists(), f"annotation.mmseqs_databases.{name}: not found: {path}")
    for key in ("split_mem_limit", "max_memory"):
        if key in ann:
            req(re.fullmatch(r"\d+G", str(ann[key])), f"annotation.{key}: must look like '164G'")

    sim = cfg.get("simulation") or {}
    for key in ("read_length", "mean_fragment_length", "std_fragment_length", "threads"):
        req(isinstance(sim.get(key), int) and sim.get(key) > 0,
            f"simulation.{key}: must be a positive integer")
    req(_num(sim.get("coverage")) and sim.get("coverage") > 0, "simulation.coverage: must be > 0")
    if isinstance(sim.get("read_length"), int) and isinstance(sim.get("mean_fragment_length"), int):
        req(sim["mean_fragment_length"] >= sim["read_length"],
            "simulation.mean_fragment_length: should be >= read_length")
    r1, r2 = sim.get("empirical_reads_R1"), sim.get("empirical_reads_R2")
    req((r1 is None) == (r2 is None),
        "simulation.empirical_reads_R1/R2: set both or neither")

    scen = cfg.get("scenarios")
    req(isinstance(scen, dict) and scen, "scenarios: must be a non-empty mapping")
    for name, contribs in (scen or {}).items() if isinstance(scen, dict) else []:
        req(re.fullmatch(r"[A-Za-z0-9_.-]+", str(name)),
            f"scenarios.{name}: name must be a single path-safe word")
        if not isinstance(contribs, list) or not contribs:
            p.append(f"scenarios.{name}: must be a non-empty list")
            continue
        ids = [c.get("ref_id") for c in contribs if isinstance(c, dict)]
        req(len(ids) == len(set(ids)), f"scenarios.{name}: a ref_id is listed twice")
        for i, c in enumerate(contribs):
            where = f"scenarios.{name}[{i}]"
            if not isinstance(c, dict):
                p.append(f"{where}: must be a mapping")
                continue
            req(set(c) == {"ref_id", "mutated_fraction", "abundance"},
                f"{where}: keys must be ref_id, mutated_fraction, abundance (got {sorted(c)})")
            req(c.get("ref_id") in refs, f"{where}: unknown ref_id {c.get('ref_id')!r}")
            mf = c.get("mutated_fraction")
            req(_num(mf) and 0 <= mf <= 1, f"{where}: mutated_fraction must be in [0, 1]")
            ab = c.get("abundance")
            req(_num(ab) and ab > 0, f"{where}: abundance must be > 0")

    vc = cfg.get("variant_calling") or {}
    req(isinstance(vc.get("threads"), int) and vc.get("threads") > 0,
        "variant_calling.threads: must be a positive integer")
    mq = (cfg.get("assessment") or {}).get("min_quality")
    req(_num(mq) and mq >= 0, "assessment.min_quality: must be a number >= 0")
    rep = cfg.get("report") or {}
    req(isinstance(rep.get("window_size", 2000), int) and rep.get("window_size", 2000) > 0,
        "report.window_size: must be a positive integer")
    req(isinstance(rep.get("genome_length", 0), int) and rep.get("genome_length", 0) >= 0,
        "report.genome_length: must be an integer >= 0")
    return p


def test_shipped_config_is_valid():
    problems = validate_config(real_config())
    assert not problems, "config/config.yaml problems:\n  " + "\n  ".join(problems)


def test_validator_catches_common_mistakes():
    """Guard the validator itself: each mutation must be reported."""
    mutations = {
        "unknown ref_id": lambda c: c["scenarios"]["scenario_single_ref"][0].update(ref_id="nope"),
        "mutated_fraction": lambda c: c["scenarios"]["scenario_single_ref"][0].update(mutated_fraction=1.5),
        "abundance": lambda c: c["scenarios"]["scenario_single_ref"][0].update(abundance=0),
        "keys must be": lambda c: c["scenarios"]["scenario_single_ref"][0].update(abundace=1),
        "mutation.model": lambda c: c["mutation"].update(model="HKY"),
        "substitution_rate": lambda c: c["mutation"].update(substitution_rate="0.001"),
        "replicates": lambda c: c.update(replicates=0),
        "file not found": lambda c: c["references"].update(Acutalibacter_muris="missing.fa"),
        "set both or neither": lambda c: c["simulation"].update(empirical_reads_R1="x.fq"),
        "not found: nowhere": lambda c: c["annotation"]["mmseqs_databases"].update(x="nowhere"),
        "window_size": lambda c: c["report"].update(window_size=-5),
    }
    for expected, mutate in mutations.items():
        cfg = real_config()
        mutate(cfg)
        problems = validate_config(cfg)
        assert any(expected in msg for msg in problems), (expected, problems)


def test_reference_fastas_are_valid_and_contig_names_unique():
    """bwa_index concatenates the references of a scenario, so contig names
    must be unique across references or GATK will reject the merged FASTA."""
    seen = {}
    for rid, path in real_config()["references"].items():
        full = ROOT / path
        names, seq_chars, n_seq = [], set(), 0
        with open(full) as fh:
            first = fh.readline()
            assert first.startswith(">"), f"{path}: does not start with '>'"
            names.append(first[1:].split()[0])
            for line in fh:
                if line.startswith(">"):
                    names.append(line[1:].split()[0])
                else:
                    s = line.strip()
                    seq_chars.update(s)
                    n_seq += len(s)
        assert n_seq > 0, f"{path}: no sequence"
        bad = seq_chars - IUPAC
        assert not bad, f"{path}: non-IUPAC characters {sorted(bad)}"
        assert len(names) == len(set(names)), f"{path}: duplicate contig names"
        for n in names:
            assert n not in seen, f"contig {n!r} is in both {seen.get(n)} and {rid}"
            seen[n] = rid


def test_every_config_key_is_used():
    cfg = smk.TrackingDict(real_config())
    wf = smk.load_workflow(ROOT, cfg)
    with tmpdir() as td:
        smk.resolve_dag(wf, cwd=td)
    # (leaf_paths walks a plain copy: iterating the TrackingDict would mark everything)
    unused = [".".join(map(str, p)) for p in smk.leaf_paths(real_config()) if p not in cfg.seen]
    assert not unused, f"config keys never read by the workflow: {unused}"


# ═══════════════════════════════════════════════════════════════════════════════
# Files referenced by the workflow
# ═══════════════════════════════════════════════════════════════════════════════

def test_all_python_scripts_compile():
    """Syntax-check every script (compile only; nothing is written to disk)."""
    # assess_variants_debugging.py is a scratch copy used by no rule; it does
    # not currently parse, so it is excluded here.
    scripts = sorted(p for p in SCRIPTS.glob("*.py") if p.name != "assess_variants_debugging.py")
    assert scripts
    for path in scripts:
        compile(path.read_text(), str(path), "exec")


def test_snakefile_includes_exist_and_rules_are_unique():
    wf = real_workflow()
    assert not wf.duplicate_rules, f"duplicate rule names: {wf.duplicate_rules}"
    included = {p.name for p in wf.files}
    for f in RULES.glob("*.smk"):
        if f.name == "variant_calling_breseq.smk":       # documented alternative caller
            continue
        assert f.name in included, f"{f.name} is not included by the Snakefile"


def test_default_target_is_rule_all_without_wildcards():
    wf = real_workflow()
    first = wf.first_rule
    assert first.name == "all", f"default target is {first.name!r}"
    for pat in first.patterns("input"):
        assert not smk.wildcard_names(pat), f"rule all input has wildcards: {pat}"


def test_script_and_conda_paths_exist():
    for rule in real_workflow().rules.values():
        base = rule.file.parent
        for directive in ("script", "conda"):
            value = rule.single(directive)
            if value:
                assert (base / value).is_file(), f"{rule.location}: {directive} {value!r} missing"


def _env(path):
    data = smk.load_yaml_subset(Path(path).read_text())
    deps = [d for d in data.get("dependencies", []) if isinstance(d, str)]
    return {re.split(r"[=<>! ]", d, 1)[0].lower() for d in deps}


def test_conda_env_files_are_well_formed():
    for env in sorted(ENVS.glob("*.yaml")):
        data = smk.load_yaml_subset(env.read_text())
        assert data.get("name"), f"{env.name}: no name"
        assert data.get("channels"), f"{env.name}: no channels"
        assert data.get("dependencies"), f"{env.name}: no dependencies"


PY_PACKAGES = {"pandas": "pandas", "Bio": "biopython", "numpy": "numpy", "yaml": "pyyaml"}


def test_python_script_imports_are_provided_by_rule_env():
    for rule in real_workflow().rules.values():
        script, conda = rule.single("script"), rule.single("conda")
        if not script or not script.endswith(".py"):
            continue
        assert conda, f"{rule.location}: python script rule has no conda env"
        pkgs = _env(rule.file.parent / conda)
        assert "python" in pkgs, f"{rule.location}: env {conda} lacks python"
        tree = ast.parse((rule.file.parent / script).read_text())
        mods = {n.names[0].name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import)}
        mods |= {n.module.split(".")[0] for n in ast.walk(tree)
                 if isinstance(n, ast.ImportFrom) and n.module}
        third_party = {m for m in mods if m not in sys.stdlib_module_names}
        for m in third_party:
            assert m in PY_PACKAGES, f"{script}: unknown third-party module {m!r}; add it to PY_PACKAGES"
            assert PY_PACKAGES[m] in pkgs, f"{rule.location}: {script} imports {m} but {conda} lacks {PY_PACKAGES[m]}"


def test_report_env_provides_everything_the_rmd_loads():
    rmd = (SCRIPTS / "genome_report.Rmd").read_text()
    libs = set(re.findall(r"library\((\w+)\)", rmd)) | {"rmarkdown"}
    pkgs = _env(ENVS / "r_report.yaml")
    missing = {l for l in libs if f"r-{l.lower()}" not in pkgs}
    assert not missing, f"r_report.yaml lacks: {sorted(missing)}"
    assert "pandoc" in pkgs, "rmarkdown::render needs pandoc"


# ═══════════════════════════════════════════════════════════════════════════════
# Wildcard consistency (static, per rule)
# ═══════════════════════════════════════════════════════════════════════════════

def test_outputs_of_a_rule_share_wildcards():
    for rule in real_workflow().rules.values():
        sets = {frozenset(smk.wildcard_names(p)) for p in rule.patterns("output")}
        assert len(sets) <= 1, f"{rule.location}: outputs use different wildcards {sets}"


def test_log_and_input_wildcards_are_determined_by_output():
    for rule in real_workflow().rules.values():
        outs = rule.patterns("output")
        if not outs:
            continue
        out_wc = set(smk.wildcard_names(outs[0]))
        for directive in ("log", "input", "benchmark"):
            for pat in rule.patterns(directive):
                extra = set(smk.wildcard_names(pat)) - out_wc
                assert not extra, f"{rule.location}: {directive} {pat!r} uses {extra} not in output"


# ═══════════════════════════════════════════════════════════════════════════════
# Dry-run DAG
# ═══════════════════════════════════════════════════════════════════════════════

def test_dag_builds_for_shipped_config():
    jobs, leaves, errors = real_dag()
    assert not errors, "DAG errors:\n  " + "\n  ".join(errors)
    counts = Counter(j.rule.name for j in jobs)
    reps = real_config()["replicates"]
    n_scen = len(real_config()["scenarios"])
    n_refs = len(real_config()["references"])
    assert counts["mutate_reference"] == n_refs * reps
    assert counts["simulate_reads"] == n_refs * reps * 2
    for rule in ("blend_reads", "bwa_index", "bwa_align", "gatk_haplotype_caller",
                 "assess_variants"):
        assert counts[rule] == n_scen * reps, (rule, counts)
    assert counts["aggregate_replicates"] == n_scen
    assert counts["generate_report"] == n_refs
    assert counts["cluster_genes_by_ani"] == 1
    assert "resources/CP065321.fa" in leaves


def test_dag_builds_for_multi_reference_config_with_empirical_profile():
    with tmpdir() as td:
        cfg = _synthetic_config(td)
        problems = validate_config(cfg, root=td)
        assert not problems, problems
        wf = smk.load_workflow(ROOT, cfg)
        jobs, leaves, errors = smk.resolve_dag(wf, cwd=td)
    assert not errors, "DAG errors:\n  " + "\n  ".join(errors)
    counts = Counter(j.rule.name for j in jobs)
    assert counts["art_learn_profile"] == 1
    assert counts["mutate_reference"] == 4 and counts["simulate_reads"] == 8
    assert counts["mmseqs_search"] == 4, "2 refs × 2 databases"
    assert counts["generate_report"] == 2
    assert counts["assess_variants"] == 6
    sim = next(j for j in jobs if j.rule.name == "simulate_reads")
    assert "results/art_profile/empirical_R1.txt" in sim.input


def test_every_job_output_is_unique():
    jobs, _, _ = real_dag()
    produced = Counter(p for j in jobs for p in j.output)
    dup = [p for p, n in produced.items() if n > 1]
    assert not dup, f"files produced by more than one job: {dup}"


def test_every_job_shell_command_renders():
    """Catches {input.typo}/{params.missing} errors that Snakemake only raises
    when the job is about to run."""
    with tmpdir() as td:
        cfg = _synthetic_config(td)
        jobs, _, errors = smk.resolve_dag(smk.load_workflow(ROOT, cfg), cwd=td)
    assert not errors
    for job in jobs + real_dag()[0]:
        try:
            smk.render_shell(job)
        except Exception as exc:
            raise AssertionError(f"{job.rule.location} {job.wildcards}: {type(exc).__name__}: {exc}")


TOOL_PACKAGES = {
    "bwa": "bwa", "samtools": "samtools", "gatk": "gatk4",
    "art_illumina": "art", "art_profiler_illumina": "art",
    "mmseqs": "mmseqs2", "prodigal": "prodigal", "seqtk": "seqtk",
    "Rscript": "r-base", "breseq": "breseq",
}


def _commands(shell):
    shell = shell.replace("\\\n", " ")
    for line in shell.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        for part in re.split(r"\|\||&&|[;|]", line):
            words = [w for w in part.split() if not re.match(r"^\w+=", w)]
            if words:
                yield words[0]


def test_shell_tools_are_provided_by_rule_env():
    for job in {j.rule.name: j for j in real_dag()[0]}.values():
        shell = smk.render_shell(job)
        if not shell:
            continue
        tools = {t for t in _commands(shell) if t in TOOL_PACKAGES}
        if not tools:
            continue
        conda = job.rule.single("conda")
        assert conda, f"{job.rule.location}: uses {tools} but has no conda env"
        pkgs = _env(job.rule.file.parent / conda)
        missing = {t for t in tools if TOOL_PACKAGES[t] not in pkgs}
        assert not missing, f"{job.rule.location}: {conda} lacks packages for {missing}"


# ═══════════════════════════════════════════════════════════════════════════════
# Shell-block lint
# ═══════════════════════════════════════════════════════════════════════════════

CONTINUATION_START = ("&>", ">", "2>", "|", "&&", "--")
ENV_VARS = {"SLURM_JOB_ID", "TMPDIR", "HOME", "PATH", "USER", "PWD"}


def _orphaned_lines(shell):
    prev, bad = None, []
    for line in shell.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if prev is not None and not prev.endswith("\\") and s.startswith(CONTINUATION_START):
            bad.append(s)
        prev = s
    return bad


def _undefined_vars(shell):
    text = shell.replace("{{", "{").replace("}}", "}")
    assigned = set(re.findall(r"(?m)^\s*(\w+)=", text)) | set(re.findall(r"\bfor\s+(\w+)\s+in\b", text))
    used = set(re.findall(r"\$\{?([A-Za-z_]\w*)", text))
    return used - assigned - ENV_VARS


def _shell_problems(check):
    found = {}
    for rule in real_workflow().rules.values():
        shell = rule.single("shell")
        if isinstance(shell, str):
            bad = check(shell)
            if bad:
                found[rule.name] = bad
    return found


def test_no_orphaned_shell_continuations():
    found = _shell_problems(_orphaned_lines)
    assert not found, f"lines that look like continuations of a line without '\\': {found}"


def test_no_undefined_shell_variables():
    found = _shell_problems(_undefined_vars)
    assert not found, f"undefined shell variables: {found}"


def test_mmseqs_search_shell_is_well_formed():
    shell = real_workflow().rules["mmseqs_search"].single("shell")
    assert not _orphaned_lines(shell), _orphaned_lines(shell)
    assert not _undefined_vars(shell), _undefined_vars(shell)
    assert "[[ dbpath" not in shell


FAKE_MMSEQS = """#!/bin/sh
echo "$@" >> "$FAKE_MMSEQS_CALLS"
echo "mmseqs $1 ran"
"""


def _run_mmseqs_search_shell(td, db_path):
    """Render the mmseqs_search shell for one job and run it in bash with a
    fake mmseqs on PATH; return (list of mmseqs argument strings, log text)."""
    import os
    import shutil
    import stat
    import subprocess
    if not shutil.which("bash"):
        from harness import SkipTest
        raise SkipTest("bash not available")
    cfg = real_config()
    cfg["annotation"]["mmseqs_databases"] = {"mobileOG": str(db_path)}
    wf = smk.load_workflow(ROOT, cfg)
    ref = cfg["references"]["Acutalibacter_muris"]
    jobs, _, errors = smk.resolve_dag(
        wf, targets=["results/annotation/Acutalibacter_muris/mobileOG_hits"], cwd=td)
    assert not errors, errors
    job = next(j for j in jobs if j.rule.name == "mmseqs_search")
    shell = smk.render_shell(job)
    for p in list(job.output) + list(job.log):
        (td / p).parent.mkdir(parents=True, exist_ok=True)
    fake = td / "bin" / "mmseqs"
    fake.parent.mkdir()
    fake.write_text(FAKE_MMSEQS)
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    env = dict(os.environ, PATH=f"{fake.parent}{os.pathsep}{os.environ['PATH']}",
               FAKE_MMSEQS_CALLS=str(td / "calls.txt"))
    # Snakemake runs shell blocks with bash strict mode
    proc = subprocess.run(["bash", "-c", "set -euo pipefail; " + shell], cwd=td, env=env,
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    calls = (td / "calls.txt").read_text().splitlines()
    return calls, (td / job.log[0]).read_text()


def test_mmseqs_search_shell_with_fasta_database_builds_db_first():
    with tmpdir() as td:
        db = write(td / "dbs" / "mobileOG.faa", ">p\nMKV\n")
        calls, log = _run_mmseqs_search_shell(td, db)
        assert calls[0] == f"createdb {db} {td / 'dbs' / 'mobileOG'}", calls
        assert calls[1].startswith("easy-search ") and f" {td / 'dbs' / 'mobileOG'} " in calls[1]
        assert "mmseqs createdb ran" in log and "mmseqs easy-search ran" in log


def test_mmseqs_search_shell_with_prebuilt_database_searches_directly():
    with tmpdir() as td:
        db = write(td / "dbs" / "mobileOG", "")
        calls, log = _run_mmseqs_search_shell(td, db)
        assert len(calls) == 1 and calls[0].startswith("easy-search "), calls
        assert f" {db} " in calls[0]
        assert "mmseqs easy-search ran" in log


# ═══════════════════════════════════════════════════════════════════════════════
# Report rule ↔ R Markdown contract
# ═══════════════════════════════════════════════════════════════════════════════

def test_report_passes_only_declared_rmd_params():
    """rmarkdown::render() errors on params not declared in the YAML header."""
    rmd = (SCRIPTS / "genome_report.Rmd").read_text()
    front = rmd.split("---")[1]
    declared = set(smk.load_yaml_subset(front)["params"])
    job = next(j for j in real_dag()[0] if j.rule.name == "generate_report")
    shell = smk.render_shell(job)
    block = shell[shell.index("params"):]
    block = block[block.index("list(") + 5:]
    passed = set(re.findall(r"^\s*(\w+)\s*=", block, flags=re.M))
    assert passed, "could not find the params list in the report shell"
    assert passed <= declared, f"undeclared Rmd params: {passed - declared}"
    for name, value in re.findall(r"^\s*(\w+)\s*=\s*'([^']*)'", block, flags=re.M):
        assert "{" not in value, f"unrendered placeholder in {name}: {value}"
