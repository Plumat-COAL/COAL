"""Round 3 review fixes (overnight/integration/CHANGES.md, R3-A*, R3-D*, R3-E*)."""
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

from coal_app import core, tools

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FX = os.path.join(ROOT, "tests", "fixtures")


class RunsDirCase(unittest.TestCase):
    def setUp(self):
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)

    def tmp(self):
        d = tempfile.mkdtemp(prefix="coal_test_r3_")
        self.addCleanup(shutil.rmtree, d, True)
        return d


class TestR3D1ExamplesAreLF(unittest.TestCase):  # R3-D1
    def test_gitattributes_forces_lf_for_examples(self):
        with open(os.path.join(ROOT, ".gitattributes"), encoding="utf-8") as fh:
            rules = [ln.split() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
        self.assertIn(["examples/**", "text", "eol=lf"], rules)

    def test_git_reports_lf_for_example_files(self):
        if not os.path.isdir(os.path.join(ROOT, ".git")) and not os.path.isfile(os.path.join(ROOT, ".git")):
            self.skipTest("not a git checkout")
        if shutil.which("git") is None:
            self.skipTest("git not installed")
        files = sorted(glob.glob(os.path.join(ROOT, "examples", "*", "*.fasta")))
        self.assertTrue(files)
        rel = [os.path.relpath(f, ROOT).replace(os.sep, "/") for f in files]
        out = subprocess.run(["git", "check-attr", "eol", "--"] + rel, cwd=ROOT, capture_output=True,
                             text=True, timeout=60).stdout
        for line in out.splitlines():
            self.assertTrue(line.endswith(": eol: lf"), line)

    def test_example_files_have_no_carriage_returns(self):
        files = glob.glob(os.path.join(ROOT, "examples", "*", "*.fasta"))
        self.assertTrue(files)
        for f in files:
            with open(f, "rb") as fh:
                self.assertNotIn(b"\r", fh.read(), f)


class TestR3D2CiConditions(unittest.TestCase):  # R3-D2
    def setUp(self):
        import re
        with open(os.path.join(ROOT, ".github", "workflows", "ci.yml"), encoding="utf-8") as fh:
            self.ci = fh.read()
        self.re = re

    def _job(self, name):
        m = self.re.search(rf"^  {name}:\n(.*?)(?=^  [\w-]+:\n|\Z)", self.ci, self.re.S | self.re.M)
        self.assertIsNotNone(m, name)
        return m.group(1)

    def _step(self, job, name):
        m = self.re.search(rf"^(\s*)- name: {self.re.escape(name)}\n(.*?)(?=^\1- name:|\Z)", self._job(job),
                           self.re.S | self.re.M)
        self.assertIsNotNone(m, name)
        return m.group(2)

    def test_desktop_step_runs_even_after_app_tests_fail(self):
        step = self._step("python", "Desktop wrapper tests")
        cond = self.re.search(r"^\s*if: (.*)$", step, self.re.M)
        self.assertIsNotNone(cond)
        self.assertIn("!cancelled()", cond.group(1))
        self.assertIn("matrix.os != 'macos-latest'", cond.group(1))
        self.assertNotIn("always()", cond.group(1))

    def test_explainer_job_is_required(self):
        self.assertIsNone(self.re.search(r"^\s+continue-on-error:", self._job("explainer"), self.re.M))

    # CHANGED(FINAL-1): was test_windows_gui_job_still_informational (asserted continue-on-error on Windows);
    # Windows passed in CI, so the GUI job and the desktop step must now be required on every OS.
    def test_windows_gui_job_is_required(self):
        self.assertIsNone(self.re.search(r"^\s+continue-on-error:", self._job("gui"), self.re.M))

    def test_desktop_step_is_required_on_windows(self):
        self.assertIsNone(self.re.search(r"^\s+continue-on-error:", self._step("python", "Desktop wrapper tests"), self.re.M))


BOM = "\ufeff"


def _variants(text: str, dna: bool):
    """Edge-case spellings of one file that the engine accepts (name -> text)."""
    v = {"plain": text, "bom": BOM + text, "crlf": text.replace("\n", "\r\n"), "cr": text.replace("\n", "\r"),
         "blank_lines": "\n\n" + text.replace("\n>", "\n\n>")}
    if dna:
        v["hash_preamble"] = "# exported by tool\n# second line\n" + text
        v["bom_crlf_hash"] = BOM + ("# exported by tool\n\n" + text).replace("\n", "\r\n")
        v["hash_inside"] = text.replace("\n>", "\n# note\n>")
    return v


class TestR3A1ReadersMirrorEngine(RunsDirCase):  # R3-A1
    def write(self, d, name, text):
        p = os.path.join(d, name)
        with open(p, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        return p

    def example(self, name="clean"):
        p, d = core.example_paths(name)
        with open(p, encoding="utf-8") as fh:
            pt = fh.read()
        with open(d, encoding="utf-8") as fh:
            dt = fh.read()
        return pt, dt

    def test_app_readers_give_the_engine_ids_and_sequences(self):
        from coal import seqio
        from coal_app import viewer
        pt, dt = self.example()
        d = self.tmp()
        clustal = "CLUSTAL W\n\ns1  MKV-AL\ns2  MKVSAL\n        ** **\n\ns1  WW\ns2  WW\n"
        cases = [(k, t, False) for k, t in _variants(pt, False).items()]
        cases += [("clustal_" + k, t, False) for k, t in _variants(clustal, False).items() if k != "blank_lines"]
        cases += [(k, t, True) for k, t in _variants(dt, True).items()]
        for k, text, dna in cases:
            with self.subTest(k, dna=dna):
                path = self.write(d, f"{k}_{int(dna)}.fa", text)
                if dna:
                    idx = seqio.DnaIndex([path])
                    want = [(i, idx.get(n)) for n, i in enumerate(idx.ids)]
                    idx.close()
                else:
                    want = seqio.read_alignment(path)
                self.assertTrue(want)
                self.assertEqual(core.read_ids(path), [i for i, _ in want])
                self.assertEqual(viewer.read_alignment(path, dna=dna), want)

    def run_pair(self, ptext, dtext, name="clean", **kw):
        d = self.tmp()
        p = self.write(d, "p.fa", ptext)
        n = self.write(d, "d.fa", dtext)
        return tools.codon_align(p, n, **kw)

    def test_bom_dna_with_mito_protein_runs(self):
        pt, dt = self.example("mito")
        r = self.run_pair(pt, BOM + dt, codontable=1)
        self.assertTrue(r["completed"], r.get("next_step"))
        self.assertTrue(any(h.get("suggested_codontable") == 2 for h in r["hints"]), r["hints"])
        self.assertEqual(r["summary"]["total"], 24, r["summary"])

    def test_hash_preamble_dna_with_mito_protein_runs(self):
        pt, dt = self.example("mito")
        r = self.run_pair(pt, "# exported by tool\n" + dt, codontable=1)
        self.assertTrue(r["completed"], r.get("next_step"))
        self.assertTrue(any(h.get("suggested_codontable") == 2 for h in r["hints"]), r["hints"])
        self.assertTrue(r["problems"])  # per-gene warnings, not only run-level tips

    def test_clean_example_variants_pair_by_id(self):
        pt, dt = self.example()
        for k, (ptext, dtext) in {"bom_dna": (pt, BOM + dt), "hash_dna": (pt, "# x\n" + dt),
                                  "bom_both": (BOM + pt, BOM + dt),
                                  "bom_crlf_both": (BOM + pt.replace("\n", "\r\n"),
                                                    BOM + ("# x\n" + dt).replace("\n", "\r\n"))}.items():
            with self.subTest(k):
                r = self.run_pair(ptext, dtext)
                self.assertTrue(r["completed"], r.get("next_step"))
                s = r["summary"]
                self.assertEqual((s["total"], s["ok"], s["pairing"]), (8, 8, "id"), s)
                self.assertFalse([h for h in r["hints"] if h.get("code") == "id_mismatch"], r["hints"])
                self.assertNotIn("Of 0 sequences", r["methods_text"])


def _masked_files(d, frac=0.55, n=100, seed=7):
    """One gene with `frac` of its codons NNN (correctly placed) plus one clean gene."""
    import random
    from coal_app import gcode
    sense = sorted(c for c in (a + b + e for a in "TCAG" for b in "TCAG" for e in "TCAG")
                   if gcode.translate(c, 1) not in ("*", "M"))
    rng = random.Random(seed)
    cod = ["ATG"] + [rng.choice(sense) for _ in range(n - 1)]
    prot = "".join(gcode.translate(c, 1) for c in cod)
    for i in rng.sample(range(1, n), int(frac * n)):
        cod[i] = "NNN"
    p, q = os.path.join(d, "p.fa"), os.path.join(d, "n.fa")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(f">masked\n{prot}\n>clean\nMKVAL{'-' * (n - 5)}\n")
    with open(q, "w", encoding="utf-8") as fh:
        fh.write(f">masked\n{''.join(cod)}\n>clean\nATGAAAGTTGCTCTG\n")
    return p, q, int(frac * n)


class TestR3A2AmbiguousCodons(RunsDirCase):  # R3-A2
    def test_nnn_codons_are_not_called_a_genetic_code_problem(self):
        p, q, k = _masked_files(self.tmp())
        r = tools.codon_align(p, q)
        self.assertTrue(r["completed"], r.get("next_step"))
        prob = next(x for x in r["problems"] if x["id"] == "masked")
        self.assertIn("unknown or ambiguous bases", prob["explanation"])
        self.assertIn(f"{k} codons", prob["explanation"])
        self.assertNotIn("genetic code", prob["explanation"])
        self.assertNotIn("do not match their amino acid", prob["explanation"])
        self.assertFalse([h for h in r["hints"] if h.get("code") == "genetic_code"], r["hints"])

    def test_nomismatch_says_why_columns_were_removed(self):
        p, q, k = _masked_files(self.tmp())
        r = tools.codon_align(p, q, nomismatch=True)
        self.assertTrue(r["completed"], r.get("next_step"))
        self.assertIn(f"{k} codon columns were removed", r["summary_text"])
        self.assertIn("unknown or ambiguous bases", r["summary_text"])
        self.assertNotIn("because they didn't match", r["summary_text"])
        prob = next(x for x in r["problems"] if x["id"] == "masked")
        self.assertIn("removed", prob["explanation"])
        self.assertNotIn("genetic code", prob["explanation"])

    def test_mixed_definite_and_ambiguous(self):
        msgs = [{"code": "mismatch", "column": 3, "residue": "G", "codon": "NNN"},
                {"code": "mismatch", "column": 7, "residue": "W", "codon": "TGA"}]
        text = core.sequence_sentence({"status": "warning", "messages": msgs})
        self.assertIn("1 codon does not match its amino acid (first at column 7", text)
        self.assertIn("1 codon contains unknown or ambiguous bases", text)
        self.assertIn("genetic code", text)

    def test_real_mismatches_keep_the_genetic_code_wording(self):  # mito example under code 1
        p, d = core.example_paths("mito")
        r = tools.codon_align(p, d, codontable=1)
        self.assertTrue(any(h.get("suggested_codontable") == 2 for h in r["hints"]), r["hints"])
        self.assertTrue(all("genetic code" in x["explanation"] for x in r["problems"]))


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


def _load_demo(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(f"r3_{name}", os.path.join(ROOT, "demo", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@unittest.skipUnless(os.path.isdir("/proc/self") and os.name != "nt", "R3-A4: needs Linux /proc and POSIX signals")
class TestR3A4DemoTimeoutsStopWorkers(unittest.TestCase):  # R3-A4
    @classmethod
    def setUpClass(cls):
        cls.data = tempfile.mkdtemp(prefix="coal_test_r3a4_")
        r = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks", "gen_data.py"), "40000", "600",
                            cls.data], capture_output=True, text=True, timeout=300)
        assert r.returncode == 0, r.stderr
        cls.pep, cls.nuc = os.path.join(cls.data, "pep.fasta"), os.path.join(cls.data, "nuc.fasta")

    @classmethod
    def tearDownClass(cls):
        for pid in _procs_mentioning(cls.data):
            try:
                os.kill(pid, 9)
            except OSError:
                pass
        shutil.rmtree(cls.data, True)

    def out_dir(self):
        d = tempfile.mkdtemp(prefix="coal_test_r3a4_out_", dir=self.data)
        return d

    def test_race_batched_timeout(self):
        import time
        race = _load_demo("race")
        r = race.Racer("COAL", 1)
        out = self.out_dir()
        race.run_batched(r, [sys.executable, "-m", "coal", self.pep, self.nuc], os.path.join(out, "o.fasta"),
                         time.perf_counter() + 0.8)
        self.assertIsNotNone(r.error)
        self.assertEqual(_wait_gone(self.data), [])
        self.assertEqual([n for n in os.listdir(out) if n.startswith(".coal_")], [])

    def test_diamond_run_batch_timeout(self):
        import argparse
        dp = _load_demo("diamond_pipeline")
        recs = dp.read_fasta(self.pep)
        dna = {rid: seq for rid, _, seq in dp.read_fasta(self.nuc)}
        items, cds = [], {"a": {}, "b": {}}
        for (x, _, rx), (y, _, ry) in zip(recs[0::2], recs[1::2]):
            items.append((x, y, rx, ry))
            cds["a"][x], cds["b"][y] = dna[x], dna[y]
        out = self.out_dir()
        a = argparse.Namespace(a_name="a", b_name="b", out=out, codontable=1, timeout=0.8)
        res = dp.run_batch(items, a, "", "batch", [], cds)
        self.assertIsInstance(res, tuple)
        self.assertIn("did not finish", res[1])
        self.assertEqual(_wait_gone(self.data), [])
        self.assertEqual([n for n in os.listdir(out) if n.startswith(".coal_")], [])


_EVENTS = {"click", "change", "upload", "load", "then", "submit", "input"}  # CHANGED(GUI-9): + input (settings)


class _Comp:
    """Stand-in for any Gradio object: records constructor arguments and event wiring."""
    REG: list = []

    def __init__(self, name="", parent=None, args=(), kw=None):
        self._name, self._parent, self._args, self._kw = name, parent, args, kw or {}
        self._origin = None  # CHANGED(PRE-1): for a chain, the (component, event) it follows (".then" wiring)

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


class TestR3A5UploadAfterPath(RunsDirCase):  # R3-A5
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

    def _label(self, comp):
        return str(comp._kw.get("label", ""))

    def test_upload_clears_its_path_box(self):
        self.ui.build()
        found = {}
        for comp, ev, a, kw in _Comp.REG:
            if ev != "upload":
                continue
            outs = kw.get("outputs", a[2] if len(a) > 2 else None)
            outs = outs if isinstance(outs, (list, tuple)) else [outs]
            fn = kw.get("fn", a[0] if a else None)
            found[self._label(comp)] = ([self._label(o) for o in outs], fn)
        # CHANGED(GUI-13): the labels were "1. Protein alignment (...)" and "2. Coding DNA (...)"; now "Protein
        # alignment" and "Coding DNA" (the formats are a hint inside the drop zone)
        prot = next((k for k in found if k.startswith("Protein")), None)
        dna = next((k for k in found if k.startswith("Coding DNA")), None)
        self.assertIsNotNone(prot, found)
        self.assertIsNotNone(dna, found)
        for k in (prot, dna):
            labels, fn = found[k]
            self.assertEqual(len(labels), 1, labels)
            # CHANGED(GUI-1): was labels "or type/paste a path ..."; the path boxes now sit under "Use a file path
            # instead (large files)" and are labelled "<file> path"
            self.assertTrue(labels[0].endswith(" path"), labels)
            self.assertEqual(labels[0].split()[0], {prot: "Protein", dna: "Coding"}[k], labels)
            self.assertEqual(fn("/tmp/x.fa"), "")
            self.assertEqual(fn(None), "")

    def test_banner_names_the_files_used(self):
        p, d = core.example_paths("clean")
        out = self.ui.run_ui(None, None, "fasta", 1, False, False, p, d)
        banner = out[0]._kw["value"]
        self.assertIn("protein.aln.fasta", banner)
        self.assertIn("dna.fasta", banner)
        self.assertIn("read in place", banner)


class TestR3E1AppParsesNewNotice(unittest.TestCase):  # R3-E1
    def test_new_notice_with_ambiguous_count(self):
        line = ("coal: warning: g: no 10-residue segment of the protein was found in the DNA; codons were assigned "
                "by position (1 of 12 checkable codons do not match; 30 contain N or other ambiguous bases)")
        sid, msg = core._parse_warning(line)
        self.assertEqual((sid, msg["code"]), ("g", "no_anchor"))
        self.assertIn("1 of 12 checkable codons do not match; 30 contain N", msg["text"])


class TestR3D3SmokeTestMatchesWrapper(unittest.TestCase):  # R3-D3
    def test_open_with_and_outside_drop_load_by_path(self):
        with open(os.path.join(ROOT, "desktop", "SMOKE_TEST.md"), encoding="utf-8") as fh:
            row = next(ln for ln in fh if ln.startswith("| 6 |"))
        self.assertNotIn("doesn't load the file yet", row)
        self.assertNotIn("box to drag it onto", row)
        self.assertIn("puts the file into the right box", row)
        self.assertIn("Only one DNA file can be loaded at a time".lower(), row.lower())


if __name__ == "__main__":
    unittest.main()
