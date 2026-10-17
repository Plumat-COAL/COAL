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

    def test_windows_gui_job_still_informational(self):
        self.assertIn("continue-on-error: ${{ matrix.os == 'windows-latest' }}", self._job("gui"))


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
        self.assertTrue(r["ok"], r.get("error"))
        self.assertTrue(any(h.get("suggested_codontable") == 2 for h in r["hints"]), r["hints"])
        self.assertEqual(r["summary"]["total"], 24, r["summary"])

    def test_hash_preamble_dna_with_mito_protein_runs(self):
        pt, dt = self.example("mito")
        r = self.run_pair(pt, "# exported by tool\n" + dt, codontable=1)
        self.assertTrue(r["ok"], r.get("error"))
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
                self.assertTrue(r["ok"], r.get("error"))
                s = r["summary"]
                self.assertEqual((s["total"], s["ok"], s["pairing"]), (8, 8, "id"), s)
                self.assertFalse([h for h in r["hints"] if h.get("code") == "id_mismatch"], r["hints"])
                self.assertNotIn("Of 0 sequences", r["methods_text"])


class TestR3E1AppParsesNewNotice(unittest.TestCase):  # R3-E1
    def test_new_notice_with_ambiguous_count(self):
        line = ("coal: warning: g: no 10-residue segment of the protein was found in the DNA; codons were assigned "
                "by position (1 of 12 checkable codons do not match; 30 contain N or other ambiguous bases)")
        sid, msg = core._parse_warning(line)
        self.assertEqual((sid, msg["code"]), ("g", "no_anchor"))
        self.assertIn("1 of 12 checkable codons do not match; 30 contain N", msg["text"])


if __name__ == "__main__":
    unittest.main()
