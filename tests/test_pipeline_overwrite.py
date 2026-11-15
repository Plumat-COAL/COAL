"""`coal pipeline`: files in the output folder that no COAL run wrote, and the ARGS options' values
(review of 2026-11-15, L1 and L3; docs/pipeline.md "Outputs" and "Options").

  * A planned output name (genes.tsv, dnds.tsv, report.json, orthologs.tsv, clusters.tsv, ...) that exists in -o DIR
    and that no COAL run recorded there (work/outputs.json, or an older report.json) is refused with exit 2 before
    anything in the folder changes; --overwrite replaces it. Re-runs into COAL's own folder, --resume, and
    `--pairs DIR/orthologs.tsv` followed by a run that searches again, go on as before.
  * `--diamond-args '-c 1 -b 4'` (a value that starts with one of COAL's own short options) is the string of
    DIAMOND options; a missing value names the `=` form.
Needs build/coal-cpp (COAL_CPP_BIN); the DeepClust part needs DIAMOND.
"""
from __future__ import annotations

import io
import json
import os
import shutil
import subprocess
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import dnds_pipeline as P  # noqa: E402
from test_dnds_pipeline import (FAM, PLAIN, Q, S, PipelineCase, needs_cpp, needs_diamond,  # noqa: E402
                                read_fasta)



def env() -> dict:
    return {**os.environ, "PYTHONPATH": str(ROOT / "src")}   # (at call time: PipelineCase sets COAL_CPP)


def run(fn, *args, **kw):
    """(result, stderr text) of a pipeline call."""
    with redirect_stderr(io.StringIO()) as err:
        res = fn(*args, **kw)
    return res, err.getvalue()


def cli(*args, timeout=300):
    r = subprocess.run([sys.executable, "-m", "coal", "pipeline", *map(str, args)], capture_output=True, text=True,
                       encoding="utf-8", timeout=timeout, stdin=subprocess.DEVNULL, env=env())
    return r.returncode, r.stderr.replace("\r\n", "\n")


def snapshot(folder: Path) -> dict:
    return {str(p.relative_to(folder)): (p.read_bytes() if p.is_file() else None) for p in sorted(folder.rglob("*"))}


@needs_cpp
class ForeignFiles(PipelineCase):
    def small(self):
        q, s = self.tmp / "q.fa", self.tmp / "s.fa"
        allq, alls = read_fasta(Q), read_fasta(S)
        q.write_text("".join(f">{a}\n{allq[a]}\n" for a, _ in PLAIN), encoding="utf-8")
        s.write_text("".join(f">{b}\n{alls[b]}\n" for _, b in PLAIN), encoding="utf-8")
        pf = self.tmp / "pairs.tsv"
        pf.write_text("".join(f"{a}\t{b}\n" for a, b in PLAIN), encoding="utf-8")
        return q, s, pf

    def mine(self, folder: Path, *names: str) -> dict:
        folder.mkdir(parents=True, exist_ok=True)
        for n in names:
            (folder / n).write_text(f"the user's own {n}\n", encoding="utf-8")
        return snapshot(folder)

    def test_refused_before_anything_changes(self):
        q, s, pf = self.small()
        out = self.tmp / "mine"
        before = self.mine(out, "genes.tsv", "dnds.tsv", "report.json", "notes.txt")
        res, err = run(P.run_pairs, q, s, out, pairs=pf)
        self.assertEqual(res.exit_code, 2, err)
        self.assertEqual(res.report["status"], "failed")
        self.assertEqual(res.report["exit_code"], 2)
        self.assertIn("already holds genes.tsv, dnds.tsv and report.json", err)
        self.assertIn("no COAL run wrote there", err)
        self.assertIn("Use a new -o folder, or add --overwrite", err)
        self.assertIn("Nothing was changed", err)
        self.assertEqual(snapshot(out), before)            # no work/, no report, the files as they were
        self.assertFalse(Path(str(out) + ".coal-lock").exists())
        # the command line: exit 2, one error line
        code, err = cli("pairs", "--query", q, "--subject", s, "--pairs", pf, "-o", out)
        self.assertEqual(code, 2, err)
        self.assertEqual(err.count("error:"), 1, err)
        self.assertEqual(snapshot(out), before)

    def test_a_link_with_an_output_name_counts(self):
        if not hasattr(os, "symlink"):
            self.skipTest("no symbolic links")
        q, s, pf = self.small()
        out = self.tmp / "out"
        out.mkdir()
        target = self.tmp / "elsewhere.tsv"
        target.write_text("kept\n", encoding="utf-8")
        try:
            (out / "dnds.tsv").symlink_to(target)
        except OSError:
            self.skipTest("cannot create symbolic links here")
        res, err = run(P.run_pairs, q, s, out, pairs=pf)
        self.assertEqual(res.exit_code, 2, err)
        self.assertIn("already holds dnds.tsv,", err)

    def test_overwrite_replaces_them(self):
        q, s, pf = self.small()
        ref, _ = run(P.run_pairs, q, s, self.tmp / "ref", pairs=pf, quiet=True)
        out = self.tmp / "mine"
        self.mine(out, "genes.tsv", "dnds.tsv", "report.json", "failed.tsv", "notes.txt")
        res, err = run(P.run_pairs, q, s, out, pairs=pf, overwrite=True)
        self.assertEqual(res.exit_code, 0, err)
        self.assertEqual((out / "dnds.tsv").read_bytes(), (self.tmp / "ref" / "dnds.tsv").read_bytes())
        self.assertEqual((out / "notes.txt").read_text(encoding="utf-8"), "the user's own notes.txt\n")
        self.assertIs(res.report["parameters"].get("overwrite", True), True)
        code, err = cli("pairs", "--query", q, "--subject", s, "--pairs", pf, "-o", self.tmp / "mine2", "-q")
        self.assertEqual(code, 0, err)
        self.mine(self.tmp / "mine3", "orthologs.tsv")
        code, err = cli("pairs", "--query", q, "--subject", s, "--pairs", pf, "-o", self.tmp / "mine3", "-q",
                        "--overwrite")
        self.assertEqual(code, 0, err)
        self.assertEqual((self.tmp / "mine3" / "dnds.tsv").read_bytes(), (self.tmp / "ref" / "dnds.tsv").read_bytes())

    def test_only_the_names_this_run_writes(self):
        """A file with an output name that this run does not write (dnds.kaks without --kaks-format, clusters.tsv
        in pair mode) is no reason to refuse: it is not touched."""
        q, s, pf = self.small()
        out = self.tmp / "out"
        self.mine(out, "dnds.kaks", "clusters.tsv", "cluster_summary.tsv", "unpaired.tsv", "dnds.tsv.tmp")
        res, err = run(P.run_pairs, q, s, out, pairs=pf)
        self.assertEqual(res.exit_code, 0, err)
        for n in ("dnds.kaks", "clusters.tsv", "cluster_summary.tsv", "unpaired.tsv"):
            self.assertEqual((out / n).read_text(encoding="utf-8"), f"the user's own {n}\n", n)
        res, err = run(P.run_pairs, q, s, out, pairs=pf, kaks_format=True)   # now dnds.kaks is planned
        self.assertEqual(res.exit_code, 2, err)
        self.assertIn("already holds dnds.kaks,", err)

    def test_the_repository_example_folder(self):
        """`-o examples/dnds` (here a copy of it) replaced the repository's orthologs.tsv."""
        ex = self.tmp / "dnds"
        shutil.copytree(ROOT / "examples" / "dnds", ex)
        before = snapshot(ex)
        code, err = cli("pairs", "--query", ex / "species_a_cds.fasta", "--subject", ex / "species_b_cds.fasta",
                        "-o", ex)
        self.assertEqual(code, 2, err)
        self.assertIn("orthologs.tsv", err)
        self.assertEqual(snapshot(ex), before)

    def test_reruns_into_coals_own_folder(self):
        """Unchanged: the same folder again, with other options, with --resume, from an older report.json."""
        q, s, pf = self.small()
        out = self.tmp / "out"
        first, _ = run(P.run_pairs, q, s, out, pairs=pf, kaks_format=True, keep_alignments=True)
        self.assertEqual(first.exit_code, 0)
        for kw in ({}, {"kaks_format": True}, {"resume": True}, {"keep_protein_alignments": True},
                   {"report": str(self.tmp / "elsewhere.json")}, {}):
            res, err = run(P.run_pairs, q, s, out, pairs=pf, **kw)
            self.assertEqual(res.exit_code, 0, (kw, err))
            self.assertNotIn("error:", err, kw)
        self.assertEqual((out / "dnds.tsv").read_bytes(), Path(first.outputs["dnds"]).read_bytes())
        # a run killed after it recorded its outputs: what it left is its own
        (out / "dnds.tsv").unlink()
        (out / "dnds.tsv.tmp").write_text("half\n", encoding="utf-8")
        res, err = run(P.run_pairs, q, s, out, pairs=pf)
        self.assertEqual(res.exit_code, 0, err)
        # a folder of a version without work/outputs.json: the outputs its report.json lists
        (out / "work" / "outputs.json").unlink()
        res, err = run(P.run_pairs, q, s, out, pairs=pf)
        self.assertEqual(res.exit_code, 0, err)

    def test_pairs_of_an_earlier_run_read_in_place_then_searched_again(self):
        """`--pairs DIR/orthologs.tsv -o DIR` keeps that file out of the clean-up list; it stays COAL's, so a
        later run that writes orthologs.tsv again is not refused."""
        q, s, pf = self.small()
        out = self.tmp / "out"
        self.assertEqual(run(P.run_pairs, q, s, out, pairs=pf)[0].exit_code, 0)
        orth = out / "orthologs.tsv"
        res, err = run(P.run_pairs, q, s, out, pairs=orth)
        self.assertEqual(res.exit_code, 0, err)
        rec = json.loads((out / "work" / "outputs.json").read_text(encoding="utf-8"))
        self.assertNotIn("orthologs.tsv", rec["files"])
        self.assertEqual(rec["read_in_place"], ["orthologs.tsv"])
        res, err = run(P.run_pairs, q, s, out, pairs=orth)             # twice in place: still recorded
        self.assertEqual(res.exit_code, 0, err)
        res, err = run(P.run_pairs, q, s, out, pairs=pf)               # writes orthologs.tsv again
        self.assertEqual(res.exit_code, 0, err)
        rec = json.loads((out / "work" / "outputs.json").read_text(encoding="utf-8"))
        self.assertIn("orthologs.tsv", rec["files"])
        self.assertEqual(rec["read_in_place"], [])
        # the user's own orthologs.tsv put into a new folder and read in place: allowed, but it stays theirs
        own = self.tmp / "own"
        own.mkdir()
        shutil.copy(pf, own / "orthologs.tsv")
        res, err = run(P.run_pairs, q, s, own, pairs=own / "orthologs.tsv")
        self.assertEqual(res.exit_code, 0, err)
        res, err = run(P.run_pairs, q, s, own, pairs=pf)
        self.assertEqual(res.exit_code, 2, err)
        self.assertIn("already holds orthologs.tsv,", err)
        self.assertEqual((own / "orthologs.tsv").read_bytes(), pf.read_bytes())

    def test_cluster_mode_clusters_tsv(self):
        out = self.tmp / "c"
        before = self.mine(out, "clusters.tsv")
        res, err = run(P.run_clusters, FAM, out)
        self.assertEqual(res.exit_code, 2, err)
        self.assertIn("already holds clusters.tsv,", err)
        self.assertEqual(snapshot(out), before)
        code, err = cli("clusters", "--cds", FAM, "-o", out)
        self.assertEqual(code, 2, err)
        self.assertEqual(snapshot(out), before)
        # the user's clusters.tsv given as --clusters: read, never written; no refusal
        res, err = run(P.run_clusters, FAM, out, clusters=out / "clusters.tsv")
        self.assertNotIn("already holds", err)

    def test_option_in_help_and_docs(self):
        code = subprocess.run([sys.executable, "-m", "coal", "pipeline", "pairs", "--help"], capture_output=True,
                              text=True, encoding="utf-8", env=env())
        self.assertIn("--overwrite", code.stdout)
        doc = (ROOT / "docs" / "pipeline.md").read_text(encoding="utf-8")
        self.assertTrue("`--overwrite`" in doc)


@needs_diamond
class ForeignFilesDeepClust(PipelineCase):
    def test_clusters_of_an_earlier_run_read_in_place_then_clustered_again(self):
        out = self.tmp / "c"
        kw = {"max_cluster_size": 2, "quiet": True}   # no cluster needs an MSA program
        first, err = run(P.run_clusters, FAM, out, **kw)
        self.assertEqual(first.exit_code in (0, 2), True, err)
        self.assertTrue((out / "clusters.tsv").is_file(), err)
        res, err = run(P.run_clusters, FAM, out, clusters=out / "clusters.tsv", **kw)
        self.assertNotIn("already holds", err)
        rec = json.loads((out / "work" / "outputs.json").read_text(encoding="utf-8"))
        self.assertEqual(rec["read_in_place"], ["clusters.tsv"])
        for extra in ({}, {"resume": True}):
            res, err = run(P.run_clusters, FAM, out, **kw, **extra)
            self.assertNotIn("already holds", err, extra)
            self.assertIn(res.exit_code, (0, 2), err)


@needs_cpp
class ArgsStrings(PipelineCase):
    def test_bind(self):
        ap = P.build_parser()
        b = P._bind_args_values
        self.assertEqual(b(["pairs", "--diamond-args", "-c 1 -b 4", "-q"], ap),
                         (["pairs", "--diamond-args=-c 1 -b 4", "-q"], None))
        self.assertEqual(b(["clusters", "--deepclust-args", "-o x", "--realign-args", "--block-size 2"], ap),
                         (["clusters", "--deepclust-args=-o x", "--realign-args=--block-size 2"], None))
        self.assertEqual(b(["pairs", "--diamond-args", "--ultra-fast-thing"], ap)[0],
                         ["pairs", "--diamond-args=--ultra-fast-thing"])
        self.assertEqual(b(["pairs", "--diamond-args", "block-size 4"], ap)[0],
                         ["pairs", "--diamond-args", "block-size 4"])          # argparse takes it as it is
        self.assertEqual(b(["pairs", "--diamond-args=-c 1"], ap)[0], ["pairs", "--diamond-args=-c 1"])
        self.assertEqual(b(["pairs", "--", "--diamond-args", "-c 1"], ap)[0], ["pairs", "--", "--diamond-args", "-c 1"])
        for argv in (["pairs", "--diamond-args"], ["pairs", "--diamond-args", "--resume"],
                     ["pairs", "--diamond-args", "-c"], ["clusters", "--realign-args", "-o", "x"]):
            self.assertIn("needs a value", b(argv, ap)[1], argv)
            self.assertIn("='-c 1 -b 4'", b(argv, ap)[1], argv)

    def test_space_form_on_the_command_line(self):
        """A value starting with one of COAL's short options reaches the ARGS check (before: "expected one
        argument"): -o is DIAMOND's here, and COAL's own, so it is refused by name."""
        for mode, opt in (("pairs", "--diamond-args"), ("clusters", "--deepclust-args"),
                          ("clusters", "--realign-args")):
            req = (["--query", "a", "--subject", "b"] if mode == "pairs" else ["--cds", "a"])
            extra = ["--realign"] if opt == "--realign-args" else []
            for form in ([opt, "-o x"], [f"{opt}=-o x"]):
                code, err = cli(mode, *req, "-o", self.tmp / "never", *extra, *form, timeout=60)
                self.assertEqual(code, 2, (form, err))
                self.assertNotIn("expected one argument", err, form)
                self.assertIn(f"{opt}: -o is an option COAL sets itself", err, form)
            code, err = cli(mode, *req, "-o", self.tmp / "never", *extra, opt, "-c", timeout=60)
            self.assertEqual(code, 2, err)
            self.assertIn(f"{opt} needs a value", err)
            self.assertIn(f"{opt}='-c 1 -b 4'", err)
        self.assertFalse((self.tmp / "never").exists())

    def test_help_and_docs(self):
        r = subprocess.run([sys.executable, "-m", "coal", "pipeline", "pairs", "--help"], capture_output=True,
                           text=True, encoding="utf-8", env=env())
        self.assertIn("'-c 1 -b 4'", " ".join(r.stdout.split()))
        doc = (ROOT / "docs" / "pipeline.md").read_text(encoding="utf-8")
        self.assertTrue("--diamond-args '-c 1 -b 4'" in doc)


@needs_diamond
class ArgsStringsWithDiamond(PipelineCase):
    def test_short_options_reach_diamond(self):
        q, s = self.tmp / "q.fa", self.tmp / "s.fa"
        allq, alls = read_fasta(Q), read_fasta(S)
        q.write_text("".join(f">{a}\n{allq[a]}\n" for a, _ in PLAIN), encoding="utf-8")
        s.write_text("".join(f">{b}\n{alls[b]}\n" for _, b in PLAIN), encoding="utf-8")
        ref = self.tmp / "ref"
        code, err = cli("pairs", "--query", q, "--subject", s, "-o", ref, "-q")
        self.assertEqual(code, 0, err)
        for form in (["--diamond-args", "-c 1 -b 4"], ["--diamond-args=-c 1 -b 4"]):
            out = self.tmp / ("o" + str(len(form)))
            code, err = cli("pairs", "--query", q, "--subject", s, "-o", out, "-q", *form)
            self.assertEqual(code, 0, (form, err))
            rep = json.loads((out / "report.json").read_text(encoding="utf-8"))
            self.assertEqual(rep["command"][-len(form):], form)          # the command as typed
            blastp = [c for st in rep["stages"] for c in st.get("commands", []) if c[1:2] == ["blastp"]]
            self.assertTrue(blastp)
            for c in blastp:
                self.assertEqual(c[-4:], ["-c", "1", "-b", "4"])
            self.assertEqual((out / "dnds.tsv").read_bytes(), (ref / "dnds.tsv").read_bytes())


if __name__ == "__main__":
    unittest.main()
