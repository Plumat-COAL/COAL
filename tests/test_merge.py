"""Merge of origin/overnight/core into the reviewed engine: one test class per reconciliation decision (MRG-1 ...
MRG-8 in overnight/integration/CHANGES.md). Inputs are built in temp dirs; every subprocess has a timeout."""
from __future__ import annotations

import io
import json
import os
import stat
import subprocess
import sys
import tempfile
import threading
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import __version__, pipeline  # noqa: E402
from coal.cli import main as coal_main  # noqa: E402
from test_report import write_genes  # noqa: E402

FORMATS = ("fasta", "clustal", "paml", "codon")
ENV = {**os.environ, "PYTHONPATH": str(ROOT / "src")}


def run(d, args):
    out, err = io.StringIO(), io.StringIO()
    cwd = os.getcwd()
    os.chdir(d)
    try:
        with redirect_stdout(out), redirect_stderr(err):
            try:
                rc = coal_main(list(args))
            except SystemExit as e:
                rc = e.code
    finally:
        os.chdir(cwd)
    return rc, out.getvalue(), err.getvalue()


def sub(args, cwd, env=None, timeout=120):
    return subprocess.run([sys.executable, "-m", "coal", *args], cwd=cwd, env={**ENV, **(env or {})},
                          capture_output=True, timeout=timeout)


class _Tmp(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = Path(self._t.name)

    def tearDown(self):
        self._t.cleanup()


# ------------------------------------------------------------------ MRG-1 safe ID pairing
class TestMrg1Pairing(_Tmp):
    """Proteins g0..g4; the DNA renamed and reordered."""

    def setUp(self):
        super().setUp()
        write_genes(self.d, 5, seed=6)
        self.recs = (self.d / "nuc.fasta").read_text().split(">")[1:]
        rc, self.ref, _ = run(self.d, ["pep.fasta", "nuc.fasta"])

    def dna(self, names, order):
        (self.d / "nuc.fasta").write_text("".join(">" + names[i] + "\n" + self.recs[i].split("\n", 1)[1]
                                                  for i in order), newline="\n")

    def test_mispairing_order_uses_normalised_ids(self):
        self.dna([f"G{i}.1" for i in range(5)], (4, 2, 0, 1, 3))
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta"])
        self.assertEqual((rc, out), (0, self.ref))  # = PAL2NAL's output on the DNA put in protein order
        self.assertIn("coal: warning: 5 of 5 protein IDs match a DNA ID only after normalising (version suffix, "
                      "case or database fields), so sequences were paired by those IDs, not by file order as "
                      "PAL2NAL pairs them: by order, 5 of 5 proteins would get the DNA of another gene (e.g. by file "
                      "order, protein 'g0' would get DNA 'G4.1'; COAL paired it with 'G0.1')", err)  # CHANGED(ENG-UX-1)
        self.assertEqual(run(self.d, ["pep.fasta", "nuc.fasta", "-q"])[2], "")  # DEC-2

    def test_same_order_stays_pal2nal(self):
        # file order already gives every protein its own DNA: PAL2NAL's pairing, output and (silent) stderr
        self.dna([f"G{i}.1" for i in range(5)], range(5))
        self.assertEqual(run(self.d, ["pep.fasta", "nuc.fasta"]), (0, self.ref, ""))
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
        r = json.loads((self.d / "r.json").read_text())
        self.assertEqual(r["summary"]["pairing"], "order")
        self.assertEqual([h["code"] for h in r["hints"]], ["id_mismatch"])

    def test_partial_exact_match_same_order_keeps_eng_d_notice(self):
        self.dna(["g0", "G1.1", "g2", "g3", "g4"], range(5))
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta"])
        self.assertEqual((rc, out), (0, self.ref))
        self.assertEqual(err, "coal: warning: only 4 of 5 protein IDs are DNA IDs; pairing sequences by order instead "
                              "(e.g. protein 'g1' gets DNA 'G1.1') (--strict makes this an error)\n")

    def test_strict_partial_match_is_an_error_even_if_normalised_ids_would_pair(self):
        self.dna(["g0", "G1.1", "g2", "g3", "g4"], (1, 0, 2, 3, 4))
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "--strict"])
        self.assertEqual((rc, out), (1, ""))
        self.assertIn("coal: error: only 4 of 5 protein IDs are DNA IDs, so sequences would be paired by order", err)
        rc, out, _ = run(self.d, ["pep.fasta", "nuc.fasta"])  # without --strict: paired by normalised ID
        self.assertEqual((rc, out), (0, self.ref))

    def test_no_id_matches_is_silent(self):  # MRG-1 decision 2 (ENG-d, test_no_ids_match_stays_silent)
        self.dna([f"x{i}" for i in range(5)], range(5))
        self.assertEqual(run(self.d, ["pep.fasta", "nuc.fasta"]), (0, self.ref, ""))

    def test_ambiguous_normalisation_falls_back_to_order(self):
        self.dna(["g0.1", "G0.2", "g2.1", "g3.1", "g4.1"], (1, 0, 2, 3, 4))
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
        self.assertEqual(json.loads((self.d / "r.json").read_text())["summary"]["pairing"], "order")


# ------------------------------------------------------------------ MRG-2 report / summary / hint
class TestMrg2Report(_Tmp):
    def test_notices_are_report_messages_and_options_complete(self):
        write_genes(self.d, 4, seed=8)
        lines = (self.d / "nuc.fasta").read_text().split("\n")
        import random
        rng = random.Random(1)
        lines[3] = "ATG" + "".join(rng.choice("ACGT") for _ in range(len(lines[3]) - 3))  # g1: unrelated DNA
        (self.d / "nuc.fasta").write_text("\n".join(lines), newline="\n")
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
        self.assertEqual(rc, 0)
        r = json.loads((self.d / "r.json").read_text())
        g1 = r["sequences"][1]
        self.assertEqual(g1["status"], "failed")  # DEC-1-A: written but probably wrong
        self.assertEqual(g1["messages"][0]["code"], "suspect_unrelated")
        self.assertEqual(r["options"]["strict"], False)
        self.assertIn("blockonly", r["options"])
        self.assertIn("1 was written to the output but flagged as probably wrong", r["methods_text"])
        self.assertEqual((self.d / "o.fa").read_text().count(">"), 4)

    def test_nomismatch_report_marks_removed_columns_and_stderr_stays_silent(self):
        write_genes(self.d, 3, wrong={1})
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "-nomismatch", "-o", "o.fa", "--report", "r.json"])
        self.assertEqual((rc, err), (0, ""))  # ENG-i2: as PAL2NAL; no summary line about mismatches either
        r = json.loads((self.d / "r.json").read_text())
        self.assertEqual(r["summary"]["columns_removed"], 1)
        m = r["sequences"][1]["messages"][0]
        self.assertEqual((m["code"], m["column"], m.get("removed")), ("mismatch", 21, True))

    def test_start_codon_only_ties_prefer_code_11(self):
        write_genes(self.d, 6, seed=3)
        lines = (self.d / "nuc.fasta").read_text().split("\n")
        for k in range(1, len(lines), 2):
            if lines[k]:
                lines[k] = "ATT" + lines[k][3:]  # ATT start: code 11 (and others) accept it, code 1 does not
        (self.d / "nuc.fasta").write_text("\n".join(lines), newline="\n")
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
        h = [x for x in json.loads((self.d / "r.json").read_text())["hints"] if x["code"] == "genetic_code"]
        self.assertEqual(h[0]["suggested_codontable"], 11, h)

    def test_hint_tests_every_gene(self):
        # code 2 fixes the one mito-like gene but breaks none of the others here; with 1 problem gene out of 10, the
        # hint needs the code to fit more genes in all, not just the problem one
        write_genes(self.d, 10, seed=5)
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
        self.assertEqual(json.loads((self.d / "r.json").read_text())["hints"], [])

    @unittest.skipUnless((ROOT / "src" / "coal_app").is_dir(), "needs src/coal_app")
    def test_engine_report_equals_app_fallback_on_examples(self):
        import report_vs_app as rva
        saved = rva.OPTION_SETS
        rva.OPTION_SETS = [{}, {"nomismatch": True}, {"codontable": 2}]
        try:
            total, bad = rva.compare_all(verbose=False, only={"examples/mito", "examples/one_broken_gene",
                                                               "examples/messy_ids", "examples/clean"})
        finally:
            rva.OPTION_SETS = saved
        self.assertEqual(total, 12)
        self.assertEqual(bad, {})


# ------------------------------------------------------------------ MRG-3 --keep-going composes
class TestMrg3KeepGoing(_Tmp):
    def unrelated(self, n=6, k=2, seed=4):
        write_genes(self.d, n, seed=seed)
        lines = (self.d / "nuc.fasta").read_text().split("\n")
        import random
        rng = random.Random(seed)
        lines[2 * k + 1] = "ATG" + "".join(rng.choice("ACGT") for _ in range(len(lines[2 * k + 1]) - 3))
        (self.d / "nuc.fasta").write_text("\n".join(lines), newline="\n")

    def test_strict_violation_is_skipped(self):
        self.unrelated()
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "--strict", "--keep-going", "-o", "o.fa",
                                    "--report", "r.json"])
        self.assertEqual(rc, 2)
        self.assertIn("coal: warning: skipped g2: protein and DNA are inconsistent (", err)
        self.assertEqual((self.d / "o.fa").read_text().count(">"), 5)
        r = json.loads((self.d / "r.json").read_text())
        self.assertEqual(r["sequences"][2]["status"], "failed")
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "--strict"])
        self.assertEqual((rc, out), (1, ""))  # without --keep-going: DEC-1 error, as before

    def test_default_dec1_notice_still_writes_the_row(self):
        self.unrelated()
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "--keep-going"])
        self.assertEqual(rc, 0)  # nothing skipped: DEC-1 writes it as PAL2NAL does
        self.assertEqual(out.count(">"), 6)
        self.assertIn("g2: protein and DNA are probably unrelated", err)

    def test_partial_id_strict_error_is_not_skippable(self):
        write_genes(self.d, 3)
        t = (self.d / "nuc.fasta").read_text().replace(">g1\n", ">other\n")
        (self.d / "nuc.fasta").write_text(t)
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "--strict", "--keep-going"])
        self.assertEqual((rc, out), (1, ""))

    def test_every_format_jobs_and_chunked_blocks(self):
        write_genes(self.d, 30, codons=25, broken={4, 17}, seed=7)
        old, oldm = pipeline.BATCH_NT, pipeline.BLOCK_MEMORY
        try:
            for fmt in FORMATS:
                for extra in ([], ["-nogap"], ["-nomismatch"], ["--strict"]):
                    if fmt == "codon" and extra and extra != ["--strict"]:
                        continue
                    outs = set()
                    for jobs, batch, mem in ((1, 2_000_000, 16_000_000), (2, 300, 16_000_000), (4, 200, 1)):
                        pipeline.BATCH_NT, pipeline.BLOCK_MEMORY = batch, mem
                        args = ["pep.fasta", "nuc.fasta", "--keep-going", "-output", fmt, "-j", str(jobs), *extra]
                        cwd = os.getcwd()
                        os.chdir(self.d)
                        try:
                            buf = io.StringIO()
                            with redirect_stdout(buf), redirect_stderr(io.StringIO()):
                                n = pipeline.run("pep.fasta", ["nuc.fasta"], jobs=jobs, min_nt_per_job=0,
                                                 out_format=fmt, keep_going=True, nogap="-nogap" in extra,
                                                 nomismatch="-nomismatch" in extra, strict="--strict" in extra)
                        finally:
                            os.chdir(cwd)
                        self.assertEqual(n, 2)
                        outs.add(buf.getvalue())
                    with self.subTest(fmt=fmt, extra=extra):
                        self.assertEqual(len(outs), 1)
        finally:
            pipeline.BATCH_NT, pipeline.BLOCK_MEMORY = old, oldm


# ------------------------------------------------------------------ MRG-4 --resume / --progress compose
def interrupt(cmd, cwd, after=2, extra_env=None):
    env = {**ENV, "COAL_BATCH_NT": "20000", "COAL_CHECKPOINT_EVERY": "0",
           "COAL_TEST_HOOKS": "1",  # CHANGED(IMP-E-I): the hooks need this too
           "COAL_TEST_PAUSE_AFTER_CHECKPOINT": str(after), **(extra_env or {})}
    p = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    dog = threading.Timer(60, p.kill)
    dog.start()
    try:
        for line in p.stderr:
            if line.startswith("coal: test: paused"):
                break
        p.terminate()
        p.wait(timeout=30)
    finally:
        dog.cancel()
        p.stderr.close()
    return p.returncode


@unittest.skipIf(os.name == "nt", "POSIX signals, symlinks and modes")
class TestMrg4ResumeComposes(_Tmp):
    def setUp(self):
        super().setUp()
        sys.path.insert(0, str(ROOT / "benchmarks"))
        from gen_data import generate
        generate(3000, 15, self.d, seed=5)
        self.ref = sub(["pep.fasta", "nuc.fasta"], self.d).stdout

    def test_resume_through_a_symlink_keeps_link_and_mode(self):  # ENG-a, R2-E5
        target = self.d / "real.fa"
        target.write_text("old\n")
        os.chmod(target, 0o640)
        os.symlink("real.fa", self.d / "link.fa")
        cmd = [sys.executable, "-m", "coal", "pep.fasta", "nuc.fasta", "-o", "link.fa"]
        self.assertEqual(interrupt(cmd, self.d), 143)
        self.assertTrue((self.d / "real.fa.coal-partial").exists())  # next to the real target, not the link
        self.assertTrue((self.d / "real.fa.coal-checkpoint.json").exists())
        r = sub(["pep.fasta", "nuc.fasta", "-o", "link.fa", "--resume"], self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn(b"coal: notice: resuming link.fa", r.stderr)
        self.assertTrue(os.path.islink(self.d / "link.fa"))
        self.assertEqual(target.read_bytes(), self.ref)
        self.assertEqual(stat.S_IMODE(os.stat(target).st_mode), 0o640)
        self.assertEqual(sorted(p.name for p in self.d.iterdir() if "coal" in p.name), [])

    def test_sigterm_before_any_checkpoint_cleans_up(self):  # ENG-f: nothing to resume -> nothing left
        cmd = [sys.executable, "-m", "coal", "pep.fasta", "nuc.fasta", "-o", "o.fa"]
        rc = interrupt(cmd, self.d, after=99, extra_env={"COAL_CHECKPOINT_EVERY": "3600",
                                                         "COAL_BATCH_NT": "2000000"})
        self.assertIn(rc, (0, 143))
        left = sorted(p.name for p in self.d.iterdir() if "coal" in p.name and p.name != "o.fa")
        self.assertEqual(left, [])

    def test_piped_input_cannot_resume_and_says_so(self):  # ENG-b
        r = subprocess.run([sys.executable, "-m", "coal", "pep.fasta", "/dev/stdin", "-o", "o.fa", "--resume"],
                           cwd=self.d, env=ENV, input=(self.d / "nuc.fasta").read_bytes(), capture_output=True,
                           timeout=120)  # a real pipe (a file on stdin is a regular file)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn(b"coal: notice: o.fa cannot be resumed (an input is not a regular file); starting a normal run",
                      r.stderr)
        self.assertEqual((self.d / "o.fa").read_bytes(), self.ref)
        self.assertFalse((self.d / "o.fa.coal-checkpoint.json").exists())

    def test_resume_blockonly_and_strict_keep_going(self):
        for opts in (["-output", "codon", "--strict", "--keep-going"], ["-output", "clustal", "-nomismatch"]):
            with self.subTest(opts=opts):
                ref = sub(["pep.fasta", "nuc.fasta", *opts], self.d).stdout
                cmd = [sys.executable, "-m", "coal", "pep.fasta", "nuc.fasta", "-o", "r.out", *opts]
                interrupt(cmd, self.d)
                r = sub(["pep.fasta", "nuc.fasta", "-o", "r.out", "--resume", *opts], self.d)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertEqual((self.d / "r.out").read_bytes(), ref)

    def test_progress_does_not_touch_stdout_and_q(self):
        r = sub(["pep.fasta", "nuc.fasta", "--progress", "-q"], self.d,
                env={"COAL_PROGRESS_MIN_SEQS": "100", "COAL_PROGRESS_EVERY": "0", "COAL_BATCH_NT": "5000"})
        self.assertEqual(r.stdout, self.ref)
        self.assertIn(b"coal: progress: 100% (3,000 of 3,000 sequences), done in", r.stderr)


# ------------------------------------------------------------------ MRG-8 version
class TestMrg8Version(unittest.TestCase):
    def test_version_is_0_4_0_everywhere(self):
        self.assertEqual(__version__, "0.4.0")
        self.assertIn('version = "0.4.0"', (ROOT / "pyproject.toml").read_text())
        self.assertRegex((ROOT / "CMakeLists.txt").read_text(), r"VERSION 0\.4\.0")
        site = (ROOT / "site" / "index.html").read_text(encoding="utf-8")
        self.assertRegex(site, r'id="version"[^>]*>0\.4\.0<')
        self.assertIn("## 0.4.0", (ROOT / "CHANGELOG.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
