"""`coal pipeline`: every table and alignment file is complete or absent (docs/pipeline.md "Output").

Each file is written as NAME.tmp, created anew (a leftover or a symbolic link of that name is removed, never written
through; O_EXCL | O_NOFOLLOW), and renamed when it is complete: the engine's files (dnds.tsv, dnds.kaks,
codon_alignments.fasta, protein_alignments.fasta, the per-cluster files) when the engine has ended with success, the
tables Python writes (genes.tsv, orthologs.tsv, unpaired.tsv, cluster_summary.tsv, failed.tsv, report.json) when
their last line is written, clusters.tsv when DIAMOND has ended.

  * A symbolic link in the output folder under an output's name (or its .tmp name) is replaced: its target is
    untouched (POSIX; Windows needs a privilege to make links, the test skips there without it).
  * An engine stage that fails in the middle (test hook COAL_TEST_DNDS_FAIL_AFTER_GROUPS, a test build of
    coal-cpp): exit 1, no dnds.tsv / dnds.kaks / alignment file under a final name, no .tmp left, and the report
    lists only the files that exist.
  * A Python-side failure while a table is written (injected): the table is absent, not cut short.
  * What a killed run left (NAME.tmp) is removed by the next run.
Needs build/coal-cpp (COAL_CPP_BIN).
"""
from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
import cpp_hooks  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402
from test_dnds_pipeline import CPP, FAM, FAM_CLUSTERS, PLAIN, Q, S, PipelineCase, needs_cpp  # noqa: E402

FINAL = ("dnds.tsv", "dnds.kaks", "codon_alignments.fasta", "protein_alignments.fasta")
KEEP = dict(kaks_format=True, keep_alignments=True, keep_protein_alignments=True, jobs=1)


def run(fn, *args, **kw):
    with redirect_stderr(io.StringIO()) as err:
        res = fn(*args, **kw)
    return res, err.getvalue()


def leftovers(folder: Path) -> list[str]:
    """Every .tmp / .body file under the output folder, work/ excluded."""
    return sorted(str(p.relative_to(folder)) for p in folder.rglob("*")
                  if p.name.endswith((".tmp", ".body")) and "work" not in p.relative_to(folder).parts)


def symlink(target: Path, link: Path, test: unittest.TestCase) -> None:
    try:
        os.symlink(target, link)
    except (OSError, NotImplementedError, AttributeError) as e:
        test.skipTest(f"cannot make a symbolic link here: {e}")


class Units(unittest.TestCase):
    def test_atomic_writes_through_a_temporary_file(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            out = d / "t.tsv"
            out.write_text("old\n", encoding="utf-8")
            (d / "t.tsv.tmp").write_text("a leftover of a killed run\n", encoding="utf-8")
            with P._atomic(str(out)) as fh:
                fh.write("new\n")
                self.assertEqual(out.read_text(encoding="utf-8"), "old\n")      # not yet
            self.assertEqual(out.read_bytes(), b"new\n")                         # \n on every system
            self.assertFalse((d / "t.tsv.tmp").exists())
            with self.assertRaises(RuntimeError):
                with P._atomic(str(out)) as fh:
                    fh.write("half a ta")
                    raise RuntimeError("the stage failed")
            self.assertEqual(out.read_bytes(), b"new\n")                         # the complete file stays
            self.assertEqual(sorted(p.name for p in d.iterdir()), ["t.tsv"])

    def test_a_symbolic_link_is_replaced_not_followed(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            victim = d / "victim.txt"
            victim.write_text("mine\n", encoding="utf-8")
            symlink(victim, d / "t.tsv", self)
            symlink(victim, d / "t.tsv.tmp", self)
            with P._atomic(str(d / "t.tsv")) as fh:
                fh.write("table\n")
            self.assertEqual(victim.read_text(encoding="utf-8"), "mine\n")
            self.assertFalse((d / "t.tsv").is_symlink())
            self.assertEqual((d / "t.tsv").read_text(encoding="utf-8"), "table\n")

    def test_every_table_has_its_temporary_name_on_record(self):
        """work/outputs.json names NAME.tmp for every table, so the next run removes what a killed run left."""
        for name in P._RUN_TABLES:
            self.assertIn(name + ".tmp", P._RUN_OUTPUTS)
        for name in ("codon_alignments.fasta", "protein_alignments.fasta", "dnds.kaks", "genes.tsv", "failed.tsv",
                     "orthologs.tsv", "unpaired.tsv", "clusters.tsv", "cluster_summary.tsv"):
            self.assertIn(name, P._RUN_TABLES)


@needs_cpp
class Engine(unittest.TestCase):
    """coal-cpp codon-dnds --new-files, and its per-cluster files."""

    def groups(self) -> str:
        seqs = {}
        cur = None
        for line in Q.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = ""
            else:
                seqs[cur] += line.strip()
        ids = [i for i in seqs if len(seqs[i]) % 3 == 0][:4]
        text = ""
        for g, (a, b) in enumerate(((ids[0], ids[1]), (ids[2], ids[3]))):
            text += f"G\tc{g}\tall\n" + "".join(f"M\t{i}\t{seqs[i]}\n" for i in (a, b)) + "E\n"
        return text

    def test_new_files_replaces_a_link_and_cluster_files_are_renamed(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "in.txt").write_text(self.groups(), encoding="utf-8", newline="\n")
            r = subprocess.run([str(CPP), "codon-dnds", "-i", str(d / "in.txt"), "-o", str(d / "probe.tsv"), "-q"],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            victim = d / "victim.txt"
            victim.write_text("mine\n", encoding="utf-8")
            symlink(victim, d / "out.tsv", self)
            symlink(victim, d / "k.kaks", self)
            (d / "cd").mkdir()
            cmd = [str(CPP), "codon-dnds", "-i", str(d / "in.txt"), "-o", str(d / "out.tsv"), "--kaks-out",
                   str(d / "k.kaks"), "--codon-dir", str(d / "cd"), "--new-files", "-q"]
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(victim.read_text(encoding="utf-8"), "mine\n")
            self.assertFalse((d / "out.tsv").is_symlink())
            self.assertEqual((d / "out.tsv").read_bytes(), (d / "probe.tsv").read_bytes())
            self.assertTrue((d / "k.kaks").read_text(encoding="utf-8").startswith("Sequence\t"))
            self.assertEqual(sorted(p.name for p in (d / "cd").iterdir()), ["c0.fasta", "c1.fasta"])
            # without the option a file is written in place, as before (-o through a link is the user's choice)
            symlink(victim, d / "direct.tsv", self)
            r = subprocess.run([str(CPP), "codon-dnds", "-i", str(d / "in.txt"), "-o", str(d / "direct.tsv"), "-q"],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(victim.read_bytes(), (d / "probe.tsv").read_bytes())


@needs_cpp
class Pipeline(PipelineCase):
    def test_outputs_are_complete_and_no_temporary_file_stays(self):
        res, err = run(P.run_pairs, Q, S, self.tmp / "o", pairs=self.pairs_file(PLAIN), **KEEP)
        self.assertIn(res.exit_code, (0, 2), err)
        out = self.tmp / "o"
        for name in FINAL + ("genes.tsv", "orthologs.tsv", "failed.tsv", "report.json"):
            self.assertTrue((out / name).is_file(), name)
        self.assertEqual(leftovers(out), [])
        # the engine was asked for NAME.tmp, created anew; the report shows the command as it ran
        cmd = next(c for c in res.report["tools"]["coal-cpp"]["commands"] if "codon-dnds" in c)
        self.assertIn("--new-files", cmd)
        for opt, name in (("-o", "dnds.tsv"), ("--kaks-out", "dnds.kaks"), ("--codon-out", "codon_alignments.fasta"),
                          ("--protein-out", "protein_alignments.fasta")):
            self.assertEqual(cmd[cmd.index(opt) + 1], str(out / (name + ".tmp")))
        rec = json.loads((out / "work" / "outputs.json").read_text(encoding="utf-8"))
        for name in FINAL + ("genes.tsv", "failed.tsv", "orthologs.tsv"):
            self.assertIn(name + ".tmp", rec["files"])

    def test_symlinks_in_the_output_folder_are_replaced(self):
        out = self.tmp / "o"
        out.mkdir()
        victim = self.tmp / "victim.txt"
        victim.write_text("mine\n", encoding="utf-8")
        names = ("dnds.kaks", "dnds.kaks.tmp", "dnds.tsv.tmp", "genes.tsv", "orthologs.tsv.tmp", "failed.tsv",
                 "codon_alignments.fasta", "protein_alignments.fasta.tmp")
        for name in names:
            symlink(victim, out / name, self)
        res, err = run(P.run_pairs, Q, S, out, pairs=self.pairs_file(PLAIN), **KEEP)
        self.assertIn(res.exit_code, (0, 2), err)
        self.assertEqual(victim.read_text(encoding="utf-8"), "mine\n")     # the target is untouched
        for name in FINAL + ("genes.tsv", "orthologs.tsv", "failed.tsv"):
            self.assertTrue((out / name).is_file() and not (out / name).is_symlink(), name)
        self.assertTrue((out / "dnds.kaks").read_text(encoding="utf-8").startswith("Sequence\t"))
        self.assertEqual(leftovers(out), [])
        base, _ = run(P.run_pairs, Q, S, self.tmp / "base", pairs=self.pairs_file(PLAIN), **KEEP)
        for name in FINAL:
            self.assertEqual((out / name).read_bytes(), (self.tmp / "base" / name).read_bytes(), name)

    def test_engine_failure_in_the_middle_leaves_no_truncated_file(self):
        cpp_hooks.need_test_hooks(self, CPP)
        out = self.tmp / "o"
        with mock.patch.dict(os.environ, {"COAL_TEST_HOOKS": "1", "COAL_TEST_DNDS_FAIL_AFTER_GROUPS": "1"}):
            res, err = run(P.run_pairs, Q, S, out, pairs=self.pairs_file(PLAIN), **KEEP)
        self.assertEqual(res.exit_code, 1, err)
        self.assertIn("coal-cpp codon-dnds failed", err)
        for name in FINAL:
            self.assertFalse(os.path.lexists(out / name), f"{name} is there after a failed stage")
        self.assertEqual(leftovers(out), [])
        # the tables that were complete before the stage are there, and the report lists what exists
        self.assertTrue((out / "genes.tsv").is_file() and (out / "failed.tsv").is_file())
        self.assertEqual(res.report["status"], "failed")
        for key, rel in res.report["outputs"].items():
            self.assertTrue((out / rel).exists(), f"the report lists {key}: {rel}, which is not there")
        self.assertNotIn("dnds", res.report["outputs"])
        self.assertNotIn("kaks", res.report["outputs"])
        on_disk = json.loads((out / "report.json").read_text(encoding="utf-8"))
        self.assertEqual(on_disk["status"], "failed")
        # the same folder again, without the failure: a complete run
        res, err = run(P.run_pairs, Q, S, out, pairs=self.pairs_file(PLAIN), **KEEP)
        self.assertIn(res.exit_code, (0, 2), err)
        self.assertEqual(leftovers(out), [])
        self.assertTrue(all((out / n).is_file() for n in FINAL))

    def test_engine_failure_in_cluster_mode(self):
        cpp_hooks.need_test_hooks(self, CPP)
        out = self.tmp / "o"
        with mock.patch.dict(os.environ, {"COAL_TEST_HOOKS": "1", "COAL_TEST_DNDS_FAIL_AFTER_GROUPS": "1"}):
            res, err = run(P.run_clusters, FAM, out, clusters=FAM_CLUSTERS, max_cluster_size=2, **KEEP)
        self.assertEqual(res.exit_code, 1, err)
        for name in ("dnds.tsv", "dnds.kaks", "cluster_summary.tsv"):
            self.assertFalse(os.path.lexists(out / name), name)
        self.assertEqual(leftovers(out), [])
        # per-cluster files: those of the clusters that were done are complete (each was renamed when written)
        for d in ("codon_alignments", "protein_alignments"):
            for f in (out / d).glob("*"):
                self.assertTrue(f.name.endswith(".fasta"), f.name)
                text = f.read_text(encoding="utf-8")
                self.assertTrue(text.startswith(">") and text.endswith("\n"), f.name)
        good, err = run(P.run_clusters, FAM, out, clusters=FAM_CLUSTERS, max_cluster_size=2, **KEEP)
        self.assertIn(good.exit_code, (0, 2), err)
        self.assertEqual(leftovers(out), [])
        fresh, _ = run(P.run_clusters, FAM, self.tmp / "fresh", clusters=FAM_CLUSTERS, max_cluster_size=2, **KEEP)
        for d in ("codon_alignments", "protein_alignments"):
            self.assertEqual(sorted(p.name for p in (out / d).iterdir()),
                             sorted(p.name for p in (self.tmp / "fresh" / d).iterdir()))

    def test_python_side_failure_while_a_table_is_written(self):
        """genes.tsv: the disk fills while the second set is copied (injected). No genes.tsv, cut short or not."""
        out = self.tmp / "o"
        real = P._Run.read_genes
        calls = []

        def failing(self_, genes_path, set_name, out_fh, *a, **kw):
            calls.append(set_name)
            if len(calls) == 2:
                out_fh.write("subject\thalf a li")
                raise OSError(28, "No space left on device", str(out / "genes.tsv.tmp"))
            return real(self_, genes_path, set_name, out_fh, *a, **kw)
        with mock.patch.object(P._Run, "read_genes", failing):
            res, err = run(P.run_pairs, Q, S, out, pairs=self.pairs_file(PLAIN), **KEEP)
        self.assertEqual(res.exit_code, 1, err)
        self.assertNotIn("Traceback", err)
        self.assertFalse(os.path.lexists(out / "genes.tsv"))
        self.assertEqual(leftovers(out), [])
        self.assertTrue((out / "failed.tsv").is_file())      # closed and named by the end of the run

    def test_the_next_run_removes_what_a_killed_run_left(self):
        out = self.tmp / "o"
        res, err = run(P.run_pairs, Q, S, out, pairs=self.pairs_file(PLAIN), **KEEP)
        self.assertIn(res.exit_code, (0, 2), err)
        for name in ("dnds.tsv.tmp", "codon_alignments.fasta.tmp", "dnds.kaks.tmp", "genes.tsv.tmp"):
            (out / name).write_text("cut sho", encoding="utf-8")       # as a kill -9 leaves them
        (out / "notes.tmp").write_text("the user's own\n", encoding="utf-8")
        # a run with fewer outputs: the earlier run's files and leftovers go, the user's file stays
        res, err = run(P.run_pairs, Q, S, out, pairs=self.pairs_file(PLAIN), jobs=1)
        self.assertIn(res.exit_code, (0, 2), err)
        self.assertEqual(leftovers(out), ["notes.tmp"])
        self.assertFalse((out / "dnds.kaks").exists())


if __name__ == "__main__":
    unittest.main()
