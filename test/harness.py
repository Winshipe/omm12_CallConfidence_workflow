"""
harness.py
==========
Shared, dependency-free plumbing for the CallConfidence test suite.

Nothing in here is a test.  It provides:

  * ``SkipTest`` / ``requires`` / ``known_bug``  – test markers
  * ``tmpdir``                                  – self-cleaning temp directory
  * ``NamedList`` / ``make_snakemake``          – a faithful stand-in for the
                                                  ``snakemake`` object that
                                                  Snakemake injects into scripts
  * ``run_script``                              – execute a workflow script
                                                  end-to-end with that stub
  * ``load_defs``                               – import only the ``def``/``class``
                                                  statements of a script (no
                                                  top-level side effects)
  * small fixture writers (FASTA / FASTQ / VCF / TSV)
  * ``run_tests``                               – the plain-Python test runner

Only the Python standard library is used here, so the suite can always start.
Individual tests that need pandas / Biopython / R declare that with
``@requires(...)`` and are reported as *skipped* when it is missing.
"""

import ast
import contextlib
import functools
import gzip
import importlib
import importlib.util
import logging
import os
import shutil
import sys
import tempfile
import textwrap
import time
import traceback
import unittest
from pathlib import Path

# ── Locations ────────────────────────────────────────────────────────────────

TEST_DIR = Path(__file__).resolve().parent
ROOT = TEST_DIR.parent                      # the workflow's working directory
WORKFLOW = ROOT / "workflow"
SCRIPTS = WORKFLOW / "scripts"
RULES = WORKFLOW / "rules"
ENVS = WORKFLOW / "envs"
CONFIG = ROOT / "config" / "config.yaml"


# ═══════════════════════════════════════════════════════════════════════════════
# Markers
# ═══════════════════════════════════════════════════════════════════════════════

# unittest.SkipTest is honoured by pytest and unittest too, so these test files
# still work if someone does run them under a framework.
SkipTest = unittest.SkipTest


def has_module(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def requires(*modules):
    """Skip the decorated test unless every named Python module is importable."""
    def deco(fn):
        @functools.wraps(fn)
        def wrapper(*a, **kw):
            missing = [m for m in modules if not has_module(m)]
            if missing:
                raise SkipTest(f"needs Python module(s): {', '.join(missing)}")
            return fn(*a, **kw)
        return wrapper
    return deco


class SetupError(Exception):
    """A shared fixture failed; never counted as a known bug."""


class KnownBug(Exception):
    """Raised by a ``known_bug`` test that failed exactly as expected."""


class UnexpectedPass(AssertionError):
    """Raised by a ``known_bug`` test that now passes (bug fixed?)."""


def known_bug(reason: str):
    """
    Mark a test that documents a *current* defect in the workflow.

    The test body is written to assert the *correct* behaviour.  While the bug
    exists the assertion fails and the runner reports the test as KNOWN BUG
    (it does not fail the run).  Once the bug is fixed the test starts to pass
    and the runner reports it as a failure so the marker gets removed and the
    test becomes an ordinary regression test.
    """
    def deco(fn):
        @functools.wraps(fn)
        def wrapper(*a, **kw):
            try:
                fn(*a, **kw)
            except (SkipTest, SetupError):
                raise
            except Exception as exc:
                if "pytest" in sys.modules:           # play nicely with pytest
                    import pytest
                    pytest.xfail(reason)
                raise KnownBug(f"{reason}\n  ({type(exc).__name__}: {exc})") from exc
            raise UnexpectedPass(
                f"marked known_bug but now passes — remove the marker: {reason}")
        wrapper.known_bug = reason
        return wrapper
    return deco


# ═══════════════════════════════════════════════════════════════════════════════
# Temporary files
# ═══════════════════════════════════════════════════════════════════════════════

@contextlib.contextmanager
def tmpdir():
    """Yield a fresh temporary directory (as a Path) and always delete it."""
    d = Path(tempfile.mkdtemp(prefix="callconf_test_"))
    try:
        yield d
    finally:
        shutil.rmtree(d, ignore_errors=True)


@contextlib.contextmanager
def chdir(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


@contextlib.contextmanager
def prepend_path(directory):
    """Temporarily put *directory* first on $PATH (for fake executables)."""
    old = os.environ.get("PATH", "")
    os.environ["PATH"] = f"{directory}{os.pathsep}{old}"
    try:
        yield
    finally:
        os.environ["PATH"] = old


def write(path, text: str, dedent: bool = True) -> Path:
    """Write *text* to *path* (creating parents) and return the Path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text) if dedent else text)
    return path


def write_fasta(path, records, width=60) -> Path:
    """records: iterable of (header, sequence)."""
    lines = []
    for header, seq in records:
        lines.append(f">{header}")
        lines.extend(seq[i:i + width] for i in range(0, len(seq), width))
    return write(path, "\n".join(lines) + "\n", dedent=False)


def write_fastq(path, records) -> Path:
    """records: iterable of (name, seq).  Gzipped when *path* ends in .gz."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
    return path


def read_fastq(path):
    """Return a list of (name, seq, plus, qual) with newlines stripped."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = [l.rstrip("\n") for l in fh]
    assert len(lines) % 4 == 0, f"{path}: line count {len(lines)} is not a multiple of 4"
    return [tuple(lines[i:i + 4]) for i in range(0, len(lines), 4)]


VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##source=HaplotypeCaller\n"
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
)


def write_vcf(path, calls) -> Path:
    """calls: iterable of (chrom, pos, ref, alt, qual) → GATK-shaped VCF."""
    body = "".join(
        f"{c}\t{p}\t.\t{r}\t{a}\t{q}\t.\tAF=0.5\tGT:AD\t0/1:10,10\n"
        for c, p, r, a, q in calls
    )
    return write(path, VCF_HEADER + body, dedent=False)


def read_tsv(path):
    """Return (header, rows) for a tab-separated file; '#' lines are data."""
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    header = lines[0].split("\t")
    return header, [l.split("\t") for l in lines[1:]]


# ═══════════════════════════════════════════════════════════════════════════════
# A faithful stand-in for the injected `snakemake` object
# ═══════════════════════════════════════════════════════════════════════════════

class NamedList(list):
    """
    Mimics ``snakemake.io.Namedlist``: a list whose items can also be reached by
    name.  Named items are part of the positional list too, in declaration
    order, and a named *list* value is flattened into the positional list —
    exactly like Snakemake does for ``input: a="x", b=["y", "z"]``.
    """

    def __init__(self, positional=(), named=None):
        named = dict(named or {})
        flat = []
        for item in list(positional) + list(named.values()):
            if isinstance(item, (list, tuple)):
                flat.extend(item)
            else:
                flat.append(item)
        super().__init__(flat)
        object.__setattr__(self, "_names", named)

    @classmethod
    def exact(cls, positional, names):
        """Use *positional* as-is and *names* as the name → value mapping."""
        obj = cls()
        list.extend(obj, positional)
        object.__setattr__(obj, "_names", dict(names))
        return obj

    def __getattr__(self, name):
        try:
            return self._names[name]
        except KeyError:
            raise AttributeError(
                f"NamedList has no item {name!r} (have: {sorted(self._names)})"
            ) from None

    def __getitem__(self, key):
        if isinstance(key, str):
            return self._names[key]
        return super().__getitem__(key)

    def keys(self):
        return self._names.keys()

    def items(self):
        return self._names.items()

    def get(self, key, default=None):
        return self._names.get(key, default)


class SnakemakeStub:
    """Plain attribute bag holding NamedLists, like snakemake.script.Snakemake."""

    def __init__(self, input, output, params, log, wildcards, threads, resources):
        self.input = input
        self.output = output
        self.params = params
        self.log = log
        self.wildcards = wildcards
        self.threads = threads
        self.resources = resources
        self.config = {}
        self.rule = "test_rule"


def _as_namedlist(value):
    if value is None:
        return NamedList()
    if isinstance(value, NamedList):
        return value
    if isinstance(value, dict):
        return NamedList(named=value)
    if isinstance(value, (list, tuple)):
        return NamedList(positional=value)
    return NamedList(positional=[value])


def make_snakemake(input=None, output=None, params=None, log=None,
                   wildcards=None, threads=1, resources=None):
    """Build a stub; dicts become named items, lists positional items."""
    to_str = lambda v: [str(x) for x in v] if isinstance(v, (list, tuple)) else (
        str(v) if isinstance(v, Path) else v)
    stringify = lambda d: {k: to_str(v) for k, v in d.items()} if isinstance(d, dict) \
        else to_str(d)
    return SnakemakeStub(
        input=_as_namedlist(stringify(input)),
        output=_as_namedlist(stringify(output)),
        params=_as_namedlist(params),
        log=_as_namedlist(stringify(log)),
        wildcards=_as_namedlist(wildcards),
        threads=threads,
        resources=_as_namedlist(resources),
    )


def _reset_root_logging():
    """Scripts call logging.basicConfig(); that is a no-op once configured, so
    clear root handlers between runs to give each script its own log file."""
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)
        with contextlib.suppress(Exception):
            h.close()


def run_script(script, **snakemake_kwargs):
    """
    Execute workflow/scripts/<script> exactly as Snakemake's ``script:``
    directive would: top-level code runs with a global ``snakemake`` object.
    A log file is created automatically if none is given.

    Returns the script's global namespace (handy for inspecting helpers).
    Exceptions raised by the script propagate unchanged.
    """
    path = Path(script) if Path(script).is_absolute() else SCRIPTS / script
    if not snakemake_kwargs.get("log"):
        fd, logpath = tempfile.mkstemp(prefix="callconf_log_", suffix=".log")
        os.close(fd)
        snakemake_kwargs["log"] = [logpath]
    sm = make_snakemake(**snakemake_kwargs)
    ns = {"__name__": "__main__", "__file__": str(path), "snakemake": sm}
    _reset_root_logging()
    try:
        exec(compile(path.read_text(), str(path), "exec"), ns)
    finally:
        _reset_root_logging()
    return ns


def run_job(job):
    """
    Run a DAG job (from smk.resolve_dag) whose rule uses ``script:``, passing
    the rule's real input/output/params/log names.  Paths are relative to the
    current working directory, as under Snakemake; output/log parent
    directories are created first, as Snakemake does.
    """
    script = (job.rule.file.parent / job.rule.single("script")).resolve()
    for p in list(job.output) + list(job.log):
        Path(p).parent.mkdir(parents=True, exist_ok=True)
    nl = lambda attrlist: NamedList.exact(list(attrlist), attrlist._names)
    sm = SnakemakeStub(
        input=nl(job.input), output=nl(job.output), log=nl(job.log),
        params=NamedList.exact(list(job.params.values()), job.params),
        wildcards=NamedList.exact(list(job.wildcards.values()), job.wildcards),
        threads=job.threads,
        resources=NamedList.exact(list(job.resources.values()), job.resources))
    ns = {"__name__": "__main__", "__file__": str(script), "snakemake": sm}
    _reset_root_logging()
    try:
        exec(compile(script.read_text(), str(script), "exec"), ns)
    finally:
        _reset_root_logging()
    return ns


def load_defs(script, extra_globals=None):
    """
    Import only imports, function/class definitions and constant assignments
    from a script, skipping any top-level statement that touches ``snakemake``
    (logging setup, ``main()`` calls, pipelines).  This is safer than the old
    exec-and-swallow-errors approach: nothing with side effects runs, and a
    genuinely broken definition still raises.
    """
    path = SCRIPTS / script
    source = path.read_text()
    tree = ast.parse(source, filename=str(path))
    keep = []
    for node in tree.body:
        seg = ast.get_source_segment(source, node) or ""
        if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                             ast.ClassDef, ast.AsyncFunctionDef)):
            keep.append(node)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)) and node.value is not None:
            # constants (BASES = ("A", ...)) and loggers only; anything computed
            # at top level may depend on snakemake inputs.
            try:
                ast.literal_eval(node.value)
                keep.append(node)
            except ValueError:
                if seg.split("=", 1)[-1].strip().startswith("logging.getLogger"):
                    keep.append(node)
    module = ast.Module(body=keep, type_ignores=[])
    ns = {"__name__": f"callconf_{path.stem}", "__file__": str(path)}
    ns.update(extra_globals or {})
    exec(compile(module, str(path), "exec"), ns)
    return ns


# ═══════════════════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════════════════

def collect(module):
    """All test_* callables of *module*, in definition order."""
    return [(f"{module.__name__}.{name}", obj)
            for name, obj in vars(module).items()
            if name.startswith("test_") and callable(obj)]


def run_tests(tests, verbose=False):
    """
    Run (name, fn) pairs.  Returns a process exit code: 0 when nothing failed
    (skips and known bugs are fine), 1 otherwise.
    """
    results = {"pass": [], "fail": [], "skip": [], "bug": []}
    width = max((len(n) for n, _ in tests), default=0)
    t0 = time.time()
    print(f"\nCallConfidence workflow test suite — {len(tests)} tests\n" + "─" * 70)

    for name, fn in tests:
        started = time.time()
        try:
            fn()
            status, detail = "pass", ""
        except SkipTest as exc:
            status, detail = "skip", str(exc)
        except KnownBug as exc:
            status, detail = "bug", str(exc)
        except Exception as exc:
            status = "fail"
            detail = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        results[status].append((name, detail))
        label = {"pass": "PASS", "fail": "FAIL", "skip": "SKIP", "bug": "KNOWN BUG"}[status]
        if verbose or status == "fail":
            took = time.time() - started
            print(f"  {label:<9} {name:<{width}}  {took:5.2f}s")
            if verbose and status in ("skip",):
                print(f"            ↳ {detail}")
        elif not verbose:
            print({"pass": ".", "skip": "s", "bug": "b"}[status], end="", flush=True)

    print(f"\n{'─' * 70}")
    if results["bug"]:
        print("\nKnown bugs (documented defects; these do not fail the run):")
        for name, detail in results["bug"]:
            first, *rest = detail.splitlines()
            print(f"  ⚑  {name}\n       {first}")
            for line in rest:
                print(f"       {line}")
    if results["skip"]:
        print("\nSkipped:")
        for name, detail in results["skip"]:
            print(f"  –  {name}: {detail}")
    if results["fail"]:
        print("\nFailures:")
        for name, detail in results["fail"]:
            print(f"\n  ✗  {name}")
            print(textwrap.indent(detail.rstrip(), "       "))
    print(
        f"\nResults: {len(results['pass'])} passed, {len(results['fail'])} failed, "
        f"{len(results['bug'])} known bugs, {len(results['skip'])} skipped "
        f"in {time.time() - t0:.1f}s\n"
    )
    return 1 if results["fail"] else 0
