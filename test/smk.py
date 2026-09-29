"""
smk.py
======
A small, dependency-free reader for *this* workflow's Snakefile / .smk files,
plus a YAML-subset loader for config/config.yaml.

Snakemake and PyYAML are not needed to run the test suite, but we still want
to catch the errors a ``snakemake --dry-run`` would catch: missing rules,
missing input files, bad config keys, wildcard mismatches, ambiguous rules,
typos in shell placeholders, and so on.  This module gets most of the way by
emulating the relevant parts of Snakemake's semantics:

  * Top-level Python in each file (helpers, imports, constants) is executed in
    one shared namespace in include order, just as Snakemake does.
  * Each rule's ``input/output/params/log/threads/resources/conda/script/shell``
    directive is evaluated *at definition time* in that namespace, so
    ``config[...]`` look-ups and ``rules.<name>.output.<key>`` references
    behave as they would under Snakemake.
  * ``resolve_dag`` then walks backwards from ``rule all`` like the real DAG
    builder: it matches each requested file against rule outputs, extracts
    wildcards, calls input functions, and recurses until it reaches files that
    exist on disk.

It is deliberately limited to the syntax this workflow uses; if the workflow
grows new syntax the parser raises a clear error instead of silently passing.
"""

import inspect
import itertools
import os
import re
import textwrap
from dataclasses import dataclass, field
from pathlib import Path

# ═══════════════════════════════════════════════════════════════════════════════
# YAML subset
# ═══════════════════════════════════════════════════════════════════════════════

class YamlSubsetError(ValueError):
    pass


def _strip_comment(line):
    out, quote = [], None
    for i, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = None
        elif ch in "\"'":
            quote = ch
        elif ch == "#" and (i == 0 or line[i - 1] in " \t"):
            break
        out.append(ch)
    return "".join(out).rstrip()


def _scalar(text):
    text = text.strip()
    if text == "" or text in ("~", "null", "Null", "NULL"):
        return None
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
        return text[1:-1]
    if text in ("true", "True", "TRUE"):
        return True
    if text in ("false", "False", "FALSE"):
        return False
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        return [_scalar(x) for x in inner.split(",")] if inner else []
    if text.startswith("{") and text.endswith("}"):
        inner = text[1:-1].strip()
        return {k.strip(): _scalar(v) for k, v in
                (kv.split(":", 1) for kv in inner.split(","))} if inner else {}
    for cast in (int, float):
        try:
            return cast(text)
        except ValueError:
            pass
    return text


_KEY_RE = re.compile(r"""^(?P<key>"[^"]*"|'[^']*'|[^:#"'\s][^:#]*?)\s*:(?:\s+(?P<val>.*))?$""")


def load_yaml_subset(text):
    """
    Parse block-style YAML made of mappings, sequences (including sequences of
    mappings), flow lists/maps of scalars and plain/quoted scalars.  Anchors,
    multi-line strings and documents are rejected loudly.
    """
    lines = []
    for n, raw in enumerate(text.splitlines(), 1):
        if "\t" in raw[: len(raw) - len(raw.lstrip())]:
            raise YamlSubsetError(f"line {n}: tab used for indentation")
        line = _strip_comment(raw)
        if not line.strip() or line.strip() == "---":
            continue
        if any(tok in line for tok in (" &", " *", ": |", ": >")):
            raise YamlSubsetError(f"line {n}: unsupported YAML feature: {raw!r}")
        lines.append([len(line) - len(line.lstrip(" ")), line.strip(), n])
    if not lines:
        return None
    value, i = _parse_block(lines, 0, lines[0][0])
    if i != len(lines):
        raise YamlSubsetError(f"line {lines[i][2]}: unexpected indentation")
    return value


def _is_item(content):
    return content == "-" or content.startswith("- ")


def _parse_block(lines, i, indent):
    if _is_item(lines[i][1]):
        return _parse_list(lines, i, indent)
    return _parse_mapping(lines, i, indent)


def _parse_mapping(lines, i, indent):
    result = {}
    while i < len(lines) and lines[i][0] == indent and not _is_item(lines[i][1]):
        m = _KEY_RE.match(lines[i][1])
        if not m:
            raise YamlSubsetError(f"line {lines[i][2]}: expected 'key: value', got {lines[i][1]!r}")
        key = _scalar(m.group("key"))
        if key in result:
            raise YamlSubsetError(f"line {lines[i][2]}: duplicate key {key!r}")
        val = m.group("val")
        i += 1
        if val is not None and val.strip():
            result[key] = _scalar(val)
        elif (i < len(lines) and lines[i][0] > indent
              and not _is_item(lines[i][1]) and not _KEY_RE.match(lines[i][1])):
            # value written on the following, more-indented line(s)
            parts = []
            while i < len(lines) and lines[i][0] > indent:
                parts.append(lines[i][1])
                i += 1
            result[key] = _scalar(" ".join(parts))
        elif i < len(lines) and lines[i][0] > indent:
            result[key], i = _parse_block(lines, i, lines[i][0])
        elif i < len(lines) and lines[i][0] == indent and _is_item(lines[i][1]):
            result[key], i = _parse_list(lines, i, indent)
        else:
            result[key] = None
    return result, i


def _parse_list(lines, i, indent):
    items = []
    while i < len(lines) and lines[i][0] == indent and _is_item(lines[i][1]):
        content = lines[i][1][1:].strip()
        if not content:
            i += 1
            value, i = _parse_block(lines, i, lines[i][0])
            items.append(value)
        elif _KEY_RE.match(content) and not content.startswith(("'", '"', "[", "{")):
            # "- key: value" starts a mapping whose keys sit at indent + 2
            lines[i] = [indent + 2, content, lines[i][2]]
            value, i = _parse_mapping(lines, i, indent + 2)
            items.append(value)
        else:
            items.append(_scalar(content))
            i += 1
    return items, i


def load_config(path):
    """Load a config file with PyYAML if available, else the subset parser."""
    text = Path(path).read_text()
    try:
        import yaml
    except ImportError:
        yaml = None
    return yaml.safe_load(text) if yaml else load_yaml_subset(text)


# ═══════════════════════════════════════════════════════════════════════════════
# Snakemake built-ins
# ═══════════════════════════════════════════════════════════════════════════════

# {name} or {name,regex}; ignores {{escaped}} braces.
WILDCARD_RE = re.compile(r"(?<!\{)\{\s*(\w+)\s*(?:,\s*((?:[^{}]|\{[^{}]*\})+))?\}(?!\})")


def wildcard_names(pattern):
    return [m.group(1) for m in WILDCARD_RE.finditer(pattern)]


def format_wildcards(pattern, values, strict=True):
    """Substitute {name} with values[name]; {{ }} become literal braces."""
    def repl(m):
        name = m.group(1)
        if name not in values:
            if strict:
                raise WildcardError(f"no value for wildcard {{{name}}} in {pattern!r}")
            return m.group(0)
        return str(values[name])
    out = WILDCARD_RE.sub(repl, pattern)
    return out.replace("{{", "{").replace("}}", "}") if strict else out


def pattern_regex(pattern):
    """Compile an output pattern into a regex with named groups (default '.+')."""
    parts, last, seen = [], 0, set()
    for m in WILDCARD_RE.finditer(pattern):
        parts.append(re.escape(pattern[last:m.start()]))
        name, constraint = m.group(1), m.group(2)
        if name in seen:
            parts.append(f"(?P={name})")
        else:
            parts.append(f"(?P<{name}>{constraint or '.+'})")
            seen.add(name)
        last = m.end()
    parts.append(re.escape(pattern[last:]))
    return re.compile("".join(parts) + r"\Z")


class WildcardError(Exception):
    pass


class DagError(Exception):
    pass


def expand(patterns, **wildcards):
    if isinstance(patterns, str):
        patterns = [patterns]
    keys = list(wildcards)
    values = [[v] if isinstance(v, (str, int, float)) else list(v)
              for v in wildcards.values()]
    out = []
    for p in patterns:
        for combo in itertools.product(*values):
            out.append(format_wildcards(p, dict(zip(keys, combo))))
    return out


class Unpack:
    def __init__(self, fn):
        self.fn = fn


class Wildcards:
    def __init__(self, values):
        self.__dict__.update(values)
        self._values = dict(values)

    def __getitem__(self, k):
        return self._values[k]

    def get(self, k, default=None):
        return self._values.get(k, default)

    def __repr__(self):
        return f"Wildcards({self._values})"


class AttrList(list):
    """List with named access; used for evaluated input/output of one job."""
    def __init__(self, items=(), names=None):
        super().__init__(items)
        self._names = dict(names or {})

    def __getattr__(self, name):
        try:
            return self._names[name]
        except KeyError:
            raise AttributeError(name) from None

    def __getitem__(self, k):
        return self._names[k] if isinstance(k, str) else super().__getitem__(k)

    def keys(self):
        return self._names.keys()

    def __str__(self):
        # Snakemake renders a file list in shell commands space-separated
        return " ".join(str(x) for x in self)


class RulesProxy:
    def __init__(self):
        self._rules = {}

    def __getattr__(self, name):
        try:
            return self._rules[name]
        except KeyError:
            raise AttributeError(f"rules.{name} referenced before definition") from None


class _RuleRef:
    def __init__(self, output):
        self.output = output


def _identity(x, *a, **k):
    return x


def snakemake_globals(config):
    return {
        "config": config,
        "expand": expand,
        "temp": _identity, "protected": _identity, "ancient": _identity,
        "directory": _identity, "touch": _identity, "pipe": _identity,
        "unpack": Unpack,
        "rules": RulesProxy(),
        "__builtins__": __builtins__,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Parsing
# ═══════════════════════════════════════════════════════════════════════════════

RULE_RE = re.compile(r"^(rule|checkpoint)\s+(\w+)\s*:\s*$")
TOP_DIRECTIVE_RE = re.compile(
    r"^(configfile|include|ruleorder|localrules|wildcard_constraints|workdir|report)\s*:\s*(.*)$")
RULE_DIRECTIVE_RE = re.compile(r"^    (\w+)\s*:\s*(.*)$")
EVALUATED = ("input", "output", "log", "params", "threads", "resources",
             "conda", "script", "shell", "benchmark", "message")
KNOWN_DIRECTIVES = set(EVALUATED) | {
    "run", "shadow", "priority", "group", "wildcard_constraints", "retries",
    "container", "envmodules", "notebook", "wrapper", "cache", "localrule",
    "handover", "default_target", "__doc__",
}


@dataclass
class Rule:
    name: str
    file: Path
    lineno: int
    raw: dict = field(default_factory=dict)          # directive -> (lineno, text)
    args: dict = field(default_factory=dict)         # directive -> list
    kwargs: dict = field(default_factory=dict)       # directive -> dict

    def items(self, directive):
        """All (key|None, value) pairs of an evaluated directive."""
        return ([(None, v) for v in self.args.get(directive, [])]
                + list(self.kwargs.get(directive, {}).items()))

    def patterns(self, directive):
        """Flattened list of string patterns of a directive (functions skipped)."""
        out = []
        for _, v in self.items(directive):
            for x in (v if isinstance(v, (list, tuple)) else [v]):
                if isinstance(x, str):
                    out.append(x)
        return out

    def single(self, directive):
        vals = self.args.get(directive, [])
        return vals[0] if vals else None

    @property
    def location(self):
        return f"{self.file.name}:{self.lineno} (rule {self.name})"


class Workflow:
    def __init__(self, snakefile, config, root):
        self.root = Path(root).resolve()
        self.config = config
        self.ns = snakemake_globals(config)
        self.rules = {}                 # name -> Rule (definition order)
        self.files = []                 # files read, in include order
        self.duplicate_rules = []
        self._read(Path(snakefile))

    # ── file walking ──────────────────────────────────────────────────────────
    def _read(self, path):
        path = path.resolve()
        self.files.append(path)
        lines = path.read_text().splitlines()
        i, code, code_start, in_triple = 0, [], 1, False

        def flush():
            nonlocal code
            if code:
                src = "\n" * (code_start - 1) + "\n".join(code)
                exec(compile(src, str(path), "exec"), self.ns)
            code = []

        while i < len(lines):
            line = lines[i]
            if not in_triple and (m := RULE_RE.match(line)):
                flush()
                start = i + 1
                i += 1
                body = []
                while i < len(lines) and (not lines[i].strip()
                                          or lines[i][0] in " \t#"):
                    body.append(lines[i])
                    i += 1
                self._define_rule(m.group(2), path, start, body)
                code_start = i + 1
                continue
            if not in_triple and (m := TOP_DIRECTIVE_RE.match(line)):
                flush()
                self._top_directive(m.group(1), m.group(2), path)
                i += 1
                code_start = i + 1
                continue
            if not code:
                code_start = i + 1
            code.append(line)
            if line.count('"""') % 2:
                in_triple = not in_triple
            i += 1
        flush()

    def _top_directive(self, name, value, path):
        if name == "include":
            target = eval(value, self.ns)
            self._read(path.parent / target)
        elif name == "configfile":
            pass    # config is supplied by the caller (like --configfile)
        else:
            raise NotImplementedError(f"{path.name}: top-level '{name}:' not supported by smk.py")

    def _define_rule(self, name, path, lineno, body):
        rule = Rule(name, path, lineno)
        current, buf, buf_line, in_triple = "__doc__", [], lineno + 1, False
        for offset, line in enumerate(body, start=lineno + 1):
            m = None if in_triple else RULE_DIRECTIVE_RE.match(line)
            if m:
                rule.raw[current] = (buf_line, "\n".join(buf))
                current, buf, buf_line = m.group(1), [m.group(2)], offset
            else:
                buf.append(line)
            if line.count('"""') % 2:
                in_triple = not in_triple
        rule.raw[current] = (buf_line, "\n".join(buf))

        unknown = set(rule.raw) - KNOWN_DIRECTIVES
        if unknown:
            raise NotImplementedError(f"{rule.location}: unknown directive(s) {sorted(unknown)}")

        for directive in EVALUATED:
            if directive not in rule.raw:
                continue
            line, text = rule.raw[directive]
            src = "\n" * (line - 1) + "__collect__(" + text + "\n)"
            self.ns["__collect__"] = lambda *a, **k: (list(a), dict(k))
            try:
                rule.args[directive], rule.kwargs[directive] = eval(
                    compile(src, str(path), "eval"), self.ns)
            except Exception as exc:
                raise type(exc)(f"{rule.location}: evaluating '{directive}': {exc}") from exc

        if name in self.rules:
            self.duplicate_rules.append(name)
        self.rules[name] = rule
        out = AttrList(rule.patterns("output"), {
            k: v for k, v in rule.kwargs.get("output", {}).items()})
        self.ns["rules"]._rules[name] = _RuleRef(out)

    def run_block(self, rule_name, input=(), output=(), log=(), params=None,
                  wildcards=None, threads=1):
        """Execute a rule's ``run:`` block with the given concrete files."""
        rule = self.rules[rule_name]
        line, text = rule.raw["run"]
        body = textwrap.dedent("\n".join(text.splitlines()[1:]))
        src = "\n" * line + body

        def named(v):
            return v if isinstance(v, AttrList) else AttrList(
                list(v.values()) if isinstance(v, dict) else list(v),
                v if isinstance(v, dict) else {})
        ns = dict(self.ns)
        ns.update(input=named(input), output=named(output), log=named(log),
                  params=named(params or {}), wildcards=Wildcards(wildcards or {}),
                  threads=threads)
        exec(compile(src, str(rule.file), "exec"), ns)

    # ── queries ───────────────────────────────────────────────────────────────
    @property
    def first_rule(self):
        """The default target: first rule of the *main* Snakefile.  Rules from
        include:d files never become the default target (Snakemake docs,
        "Includes")."""
        main = self.files[0]
        for rule in self.rules.values():
            if rule.file == main:
                return rule
        return next(iter(self.rules.values()))

    def producers(self, path):
        """[(rule, wildcard dict)] for every rule whose output matches *path*."""
        found = []
        for rule in self.rules.values():
            for pat in rule.patterns("output"):
                m = pattern_regex(pat).match(path)
                if m:
                    found.append((rule, m.groupdict()))
                    break
        return found


# ═══════════════════════════════════════════════════════════════════════════════
# DAG resolution (the core of `snakemake --dry-run`)
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class Job:
    rule: Rule
    wildcards: dict
    input: AttrList
    output: AttrList
    log: AttrList
    params: dict
    threads: object
    resources: dict


def _call(fn, wildcards, **available):
    """Call an input/params function the way Snakemake does (by arg name)."""
    sig = inspect.signature(fn)
    args, kwargs = [], {}
    for i, (pname, p) in enumerate(sig.parameters.items()):
        if pname in available:
            kwargs[pname] = available[pname]
        elif i == 0:
            args.append(wildcards)
        elif p.default is inspect.Parameter.empty:
            raise TypeError(f"cannot supply argument {pname!r} to {fn}")
    return fn(*args, **kwargs)


def _flatten(value):
    if isinstance(value, (list, tuple)):
        return [y for x in value for y in _flatten(x)]
    return [value]


def _expand_items(rule, directive, wc, strict=True, **available):
    """Evaluate one directive for concrete wildcards → AttrList of paths."""
    wcobj = Wildcards(wc)
    positional, named = [], {}
    for key, value in rule.items(directive):
        if isinstance(value, Unpack):
            result = _call(value.fn, wcobj, **available)
            if not isinstance(result, dict):
                raise DagError(f"{rule.location}: unpack() function returned {type(result).__name__}, not dict")
            for k, v in result.items():
                vals = [format_wildcards(x, wc, strict) for x in _flatten(v)]
                named[k] = vals if isinstance(v, (list, tuple)) else vals[0]
                positional.extend(vals)
            continue
        if callable(value):
            value = _call(value, wcobj, **available)
        if not isinstance(value, (str, list, tuple)):
            raise DagError(f"{rule.location}: {directive} item {key or ''} evaluated to "
                           f"{type(value).__name__} ({value!r}); expected str or list")
        vals = [format_wildcards(x, wc, strict) for x in _flatten(value)]
        for v in vals:
            if not isinstance(v, str):
                raise DagError(f"{rule.location}: {directive} contains non-path {v!r}")
        positional.extend(vals)
        if key is not None:
            named[key] = vals if isinstance(value, (list, tuple)) else vals[0]
    return AttrList(positional, named)


def resolve_dag(workflow, targets=None, cwd=None):
    """
    Build the job graph for *targets* (default: inputs of the first rule).

    Returns (jobs, leaves, errors): jobs is a list of Job, leaves are existing
    source files, errors is a list of human-readable problems.  Input and
    params functions are executed with cwd set to *cwd* (they may write files).
    """
    jobs, leaves, errors = {}, set(), []
    cwd = Path(cwd) if cwd else None

    def exists(path):
        return (Path(path) if os.path.isabs(path) else workflow.root / path).exists()

    if targets is None:
        first = workflow.first_rule
        try:
            targets = list(_expand_items(first, "input", {}))
        except Exception as exc:
            return [], set(), [f"{first.location}: cannot evaluate inputs: {exc}"]

    stack = [(t, "<target>") for t in targets]
    visiting = set()

    def job_for(path, requester):
        producers = workflow.producers(path)
        if not producers:
            if exists(path):
                leaves.add(path)
            else:
                errors.append(f"MissingInput: {path!r} (needed by {requester}) is "
                              f"not produced by any rule and does not exist")
            return
        if len(producers) > 1:
            names = ", ".join(r.name for r, _ in producers)
            errors.append(f"AmbiguousRule: {path!r} can be produced by: {names}")
            return
        rule, wc = producers[0]
        key = (rule.name, tuple(sorted(wc.items())))
        if key in jobs:
            return jobs[key]
        if key in visiting:
            errors.append(f"Cycle detected at {path!r} ({rule.name})")
            return
        visiting.add(key)
        try:
            output = _expand_items(rule, "output", wc)
            log = _expand_items(rule, "log", wc)
            inp = _expand_items(rule, "input", wc)
            threads = rule.single("threads")
            if callable(threads):
                threads = _call(threads, Wildcards(wc), input=inp)
            threads = 1 if threads is None else threads
            params = {}
            for k, v in rule.kwargs.get("params", {}).items():
                params[k] = _call(v, Wildcards(wc), input=inp, output=output,
                                  threads=threads) if callable(v) else (
                    format_wildcards(v, wc) if isinstance(v, str) else v)
            resources = {}
            for k, v in rule.kwargs.get("resources", {}).items():
                resources[k] = _call(v, Wildcards(wc), input=inp,
                                     threads=threads) if callable(v) else v
        except Exception as exc:
            errors.append(f"{rule.location} with wildcards {wc}: "
                          f"{type(exc).__name__}: {exc}")
            visiting.discard(key)
            return
        job = Job(rule, wc, inp, output, log, params, threads, resources)
        jobs[key] = job
        for p in inp:
            job_for(p, f"{rule.name}{wc}")
        visiting.discard(key)
        return job

    old = os.getcwd()
    if cwd:
        os.chdir(cwd)
    try:
        for target, requester in stack:
            job_for(target, requester)
    finally:
        os.chdir(old)
    return list(jobs.values()), leaves, errors


def render_shell(job):
    """Format a job's shell command exactly like Snakemake's {input.x} etc.
    Raises KeyError/AttributeError/IndexError on a bad placeholder."""
    shell = job.rule.single("shell")
    if shell is None:
        return None
    wrap = lambda d: AttrList(
        [AttrList(v) if isinstance(v, list) else v for v in d.values()],
        {k: AttrList(v) if isinstance(v, list) else v for k, v in d.items()})
    return shell.format(
        input=job.input, output=job.output, log=job.log,
        params=wrap(job.params), resources=wrap(job.resources),
        wildcards=Wildcards(job.wildcards), threads=job.threads, rule=job.rule.name)


# ═══════════════════════════════════════════════════════════════════════════════
# Config access tracking (to find config keys the workflow never reads)
# ═══════════════════════════════════════════════════════════════════════════════

class TrackingDict(dict):
    """dict that records every key path read through it."""

    def __init__(self, data, path=(), seen=None):
        self._path = path
        self._seen = seen if seen is not None else set()
        super().__init__({k: _track(v, path + (k,), self._seen) for k, v in data.items()})

    def _mark(self, key):
        self._seen.add(self._path + (key,))

    def _mark_all(self):
        for k in dict.keys(self):
            self._mark(k)

    def __getitem__(self, key):
        self._mark(key)
        return super().__getitem__(key)

    def get(self, key, default=None):
        self._mark(key)
        return super().get(key, default)

    def __contains__(self, key):
        self._mark(key)
        return super().__contains__(key)

    def keys(self):
        self._mark_all()
        return super().keys()

    def values(self):
        self._mark_all()
        return super().values()

    def items(self):
        self._mark_all()
        return super().items()

    def __iter__(self):
        self._mark_all()
        return super().__iter__()

    @property
    def seen(self):
        return self._seen


def _track(value, path, seen):
    if isinstance(value, dict):
        return TrackingDict(value, path, seen)
    if isinstance(value, list):
        return [_track(v, path + (i,), seen) if isinstance(v, dict) else v
                for i, v in enumerate(value)]
    return value


def leaf_paths(data, path=()):
    """Every key path to a non-dict value (list items are not descended)."""
    out = []
    for k, v in data.items():
        if isinstance(v, dict) and v:
            out.extend(leaf_paths(v, path + (k,)))
        else:
            out.append(path + (k,))
    return out


def load_workflow(root, config, snakefile=None):
    root = Path(root)
    return Workflow(snakefile or root / "workflow" / "Snakefile", config, root)
