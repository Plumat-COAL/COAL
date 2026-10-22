"""Shared helpers for the app tests: paths, a private runs folder, temp files, the Gradio stand-in (so the GUI's
layout and event wiring are tested without Gradio), and process checks."""
import gc
import os
import random
import re
import shutil
import sys
import tempfile
import unittest

from coal_app import core

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)                         # apps/coal-app: src/coal_app, desktop/, space/, pyproject.toml
ROOT = os.path.dirname(os.path.dirname(APP))        # the repository root: src/coal, examples/, tests/, demo/, site/
FX = os.path.join(ROOT, "tests", "fixtures")


def _tmp(test):
    d = tempfile.mkdtemp(prefix="coal_test_")
    test.addCleanup(shutil.rmtree, d, True)
    return d


def _write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def _load(name, path):
    """A module from a file path (scripts/, demo/, desktop/)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class RunsDirCase(unittest.TestCase):
    """Point core.RUNS_DIR at a private folder for the test (and clean it up)."""

    def setUp(self):
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)

    def tmp(self):
        d = tempfile.mkdtemp(prefix="coal_test_")
        self.addCleanup(shutil.rmtree, d, True)
        return d


def _procs_mentioning(text):
    """PIDs of live, non-zombie processes whose command line contains `text` (Linux /proc)."""
    found = []
    for pid in os.listdir("/proc"):
        if not pid.isdigit() or int(pid) == os.getpid():
            continue
        try:
            with open(f"/proc/{pid}/cmdline", "rb") as fh:
                cmd = fh.read().decode("utf-8", "replace")
            with open(f"/proc/{pid}/stat", "rb") as fh:
                state = fh.read().decode().rsplit(")", 1)[1].split()[0]
        except OSError:
            continue
        if text in cmd and state != "Z":
            found.append(int(pid))
    return found


def _wait_gone(text, limit=10.0):
    import time
    end = time.monotonic() + limit
    left = _procs_mentioning(text)
    while left and time.monotonic() < end:
        time.sleep(0.2)
        left = _procs_mentioning(text)
    return left


# ------------------------------------------------------------------ Gradio stand-in
_EVENTS = {"click", "change", "upload", "load", "then", "submit", "input"}  # input: settings


class _Comp:
    """Stand-in for any Gradio object: records constructor arguments and event wiring."""
    REG: list = []

    def __init__(self, name="", parent=None, args=(), kw=None):
        self._name, self._parent, self._args, self._kw = name, parent, args, kw or {}
        self._origin = None  # for a chain, the (component, event) it follows (".then" wiring)

    def __getattr__(self, n):
        if n.startswith("__"):
            raise AttributeError(n)
        return _Comp(n, self)

    def __call__(self, *a, **kw):
        if self._name in _EVENTS and self._parent is not None:
            _Comp.REG.append((self._parent, self._name, a, kw))
            chain = _Comp("chain")
            chain._origin = (self._parent, self._name)
            return chain
        return _Comp(self._name, None, a, kw)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __iter__(self):
        return iter(())


def _stub_gradio():
    import types
    mod = types.ModuleType("gradio")
    mod.__getattr__ = lambda n: _Comp(n)
    mod.update = lambda **kw: _Comp("update", None, (), kw)
    return mod


class StubUiCase(unittest.TestCase):
    """coal_app.ui imported against the Gradio stand-in; self.src is ui.py's source."""

    def setUp(self):
        saved = {k: sys.modules.get(k) for k in ("gradio", "coal_app.ui")}

        def restore():
            for k, v in saved.items():
                if v is None:
                    sys.modules.pop(k, None)
                else:
                    sys.modules[k] = v
        self.addCleanup(restore)
        sys.modules["gradio"] = _stub_gradio()
        sys.modules.pop("coal_app.ui", None)
        import importlib
        self.ui = importlib.import_module("coal_app.ui")
        _Comp.REG = []
        with open(self.ui.__file__, encoding="utf-8") as fh:
            self.src = fh.read()

    def comps(self, name, **kw):
        return [o for o in gc.get_objects() if isinstance(o, _Comp) and o._name == name
                and all(o._kw.get(k) == v for k, v in kw.items())]


class SettingsDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="coal_test_settings_")
        self.addCleanup(shutil.rmtree, self.dir, True)
        for k, v in (("COAL_SETTINGS_DIR", self.dir), ("COAL_HOSTED", None)):
            old = os.environ.get(k)
            self.addCleanup(lambda k=k, old=old: os.environ.__setitem__(k, old) if old is not None
                            else os.environ.pop(k, None))
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


class UiCase(SettingsDirCase):
    def setUp(self):
        super().setUp()
        saved = {k: sys.modules.get(k) for k in ("gradio", "coal_app.ui")}

        def restore():
            for k, v in saved.items():
                if v is None:
                    sys.modules.pop(k, None)
                else:
                    sys.modules[k] = v
        self.addCleanup(restore)
        sys.modules["gradio"] = _stub_gradio()
        sys.modules.pop("coal_app.ui", None)
        import importlib
        self.ui = importlib.import_module("coal_app.ui")
        _Comp.REG = []
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)


class Page:
    """The built page under the stub: component values and the click -> .then chains, as the browser runs them."""

    def __init__(self, ui):
        _Comp.REG = []
        ui.build()
        self.events = list(_Comp.REG)
        self.vals: dict[int, object] = {}
        self.comps = []
        seen = set()

        def add(c):
            if isinstance(c, _Comp) and id(c) not in seen and c._name != "chain":
                seen.add(id(c))
                self.comps.append(c)
        for c, _ev, a, kw in self.events:
            add(c)
            for io in (a[1:3] if len(a) > 1 else ()):
                for x in (io if isinstance(io, (list, tuple)) else [io]):
                    add(x)

    def find(self, **kw):
        for c in self.comps:
            if all((c._kw.get(k) if k != "text" else (c._args[0] if c._args else None)) == v
                   for k, v in kw.items()):
                return c
        raise LookupError(kw)

    def get(self, c):
        if id(c) in self.vals:
            return self.vals[id(c)]
        if "value" in c._kw:
            return c._kw["value"]
        if c._name == "State":
            return c._args[0] if c._args else None
        return "" if c._name == "Textbox" else None

    def _set(self, c, v):
        if isinstance(v, _Comp) and v._name == "update":
            if "value" in v._kw:
                self.vals[id(c)] = v._kw["value"]
        else:
            self.vals[id(c)] = v

    def _run(self, a):
        fn, inputs, outputs = a[0], a[1], a[2]
        args = [] if inputs is None else ([self.get(x) for x in inputs] if isinstance(inputs, (list, tuple))
                                          else [self.get(inputs)])
        res = fn(*args)
        if isinstance(outputs, (list, tuple)):
            for c, v in zip(outputs, res):
                self._set(c, v)
        elif outputs is not None:
            self._set(outputs, res)

    def fire(self, comp, event="click"):
        for c, ev, a, _kw in self.events:
            if c is comp and ev == event:
                self._run(a)
                for ch, ev2, a2, _ in self.events:  # the .then steps of this event, in order
                    if ch._name == "chain" and ch._origin == (comp, event) and ev2 == "then":
                        self._run(a2)


# ------------------------------------------------------------------ examples/clean variants
CLEAN_P, CLEAN_D = core.example_paths("clean")


def _read_fasta(path):
    recs, name = [], None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
                recs.append([name, ""])
            elif name is not None:
                recs[-1][1] += line.strip()
    return recs


def _write_fasta(path, recs):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.writelines(f">{n}\n{s}\n" for n, s in recs)


class ExampleVariantCase(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="coal_fin_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def variant(self, prot=None, dna=None, name="v"):
        """examples/clean with the protein / DNA records changed by the given functions (list -> list)."""
        p = _read_fasta(CLEAN_P)
        d = _read_fasta(CLEAN_D)
        if prot:
            p = prot(p)
        if dna:
            d = dna(d)
        pp, dd = os.path.join(self.d, f"{name}.pep.fa"), os.path.join(self.d, f"{name}.dna.fa")
        _write_fasta(pp, p)
        _write_fasta(dd, d)
        return pp, dd

    def run_app(self, p, d, probe=False, **kw):
        return core.run_coal(p, d, out_dir=os.path.join(self.d, "run_" + str(random.random())[2:]),
                             probe_codes=probe, **kw)


def _randomise(i):
    def f(recs):
        rng = random.Random(7)
        recs[i][1] = "ATG" + "".join(rng.choice("ACGT") for _ in recs[i][1][3:])
        return recs
    return f


def _recode(fn, genes=None):
    """DNA records with codon i of each chosen gene replaced by fn(gene index, codon index, codon)."""
    def f(recs):
        for g, rec in enumerate(recs):
            if genes is None or g in genes:
                cods = [rec[1][i:i + 3] for i in range(0, len(rec[1]), 3)]
                rec[1] = "".join(fn(g, i, c) for i, c in enumerate(cods))
        return recs
    return f


def _anchor_ids(page):
    return set(re.findall(r'id="(coalv-[^"]+)"', page))
