"""`coal pipeline`: inputs and machines that do not cooperate (docs/pipeline.md "Limits", "Failures and exit codes").

  * gzip-compressed inputs (--query, --subject, --cds, --pairs, --clusters), found by their first bytes: decompressed
    into work/inputs/ and read from there; the tables are those of the plain files; --resume reuses the stages.
  * An operating-system error anywhere (a full disk): one error line, exit 1, a report that says "failed" or no
    report, never a traceback and never a report cut short. By injection here; on a small tmpfs by hand
    (the commit message has the transcript).
  * --pairs / --clusters whose IDs match (almost) nothing: one line with example IDs of both sides.
  * A DIAMOND older than 2.1.0: cluster mode stops with the reason, pair mode warns.
  * -j above the CPU count is clamped for DIAMOND and the aligners; a tool ended by a signal is reported as that.
  * The same file as --query and --subject: a warning.
Needs build/coal-cpp (COAL_CPP_BIN); DIAMOND is replaced by shell scripts where one is needed (POSIX).
"""
from __future__ import annotations

import errno
import gzip
import io
import json
import os
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import dnds_pipeline as P  # noqa: E402
from coal import msa  # noqa: E402
from test_dnds_pipeline import (FAM, FAM_CLUSTERS, PLAIN, Q, S, PipelineCase, needs_cpp,  # noqa: E402
                                needs_diamond)

POSIX = os.name != "nt"


def run(fn, *args, **kw):
    """(result, stderr text) of a pipeline call."""
    with redirect_stderr(io.StringIO()) as err:
        res = fn(*args, **kw)
    return res, err.getvalue()


def gz(src: Path, dst: Path) -> Path:
    dst.write_bytes(gzip.compress(src.read_bytes()))
    return dst


@needs_cpp
class GzipInputs(PipelineCase):
    def test_pairs_mode_reads_gz_inputs(self):
        pairs = self.pairs_file(PLAIN)
        plain, _ = run(P.run_pairs, Q, S, self.tmp / "plain", pairs=pairs, jobs=1)
        q, s = gz(Q, self.tmp / "q.fasta.gz"), gz(S, self.tmp / "no_suffix")
        pz = gz(pairs, self.tmp / "pairs.tsv.gz")
        res, err = run(P.run_pairs, q, s, self.tmp / "gz", pairs=pz, jobs=1)
        self.assertEqual(res.exit_code, plain.exit_code, err)
        self.check_report(res)
        for name in ("dnds", "genes", "orthologs"):
            self.assertEqual(Path(res.outputs[name]).read_bytes(), Path(plain.outputs[name]).read_bytes(), name)
        self.assertIn("3 gzip-compressed input file(s) decompressed", err)
        made = sorted(p.name for p in (self.tmp / "gz" / "work" / "inputs").iterdir())
        self.assertEqual(made, ["0_q.fasta", "1_no_suffix.txt", "2_pairs.tsv", "sources.json"])
        self.assertEqual(res.report["stages"][0]["name"], "decompress_inputs")
        self.assertEqual(res.report["parameters"]["query"], str(q))   # the report names what was given
        # --resume: the copies are kept (same compressed file), so translation is reused
        again, _ = run(P.run_pairs, q, s, self.tmp / "gz", pairs=pz, jobs=1, resume=True)
        cached = {st["name"]: st["cached"] for st in again.report["stages"]}
        self.assertTrue(cached["translate_query"] and cached["translate_subject"])
        self.assertEqual(Path(again.outputs["dnds"]).read_bytes(), Path(plain.outputs["dnds"]).read_bytes())

    def test_clusters_mode_reads_gz_inputs(self):
        plain, _ = run(P.run_clusters, FAM, self.tmp / "plain", clusters=FAM_CLUSTERS, jobs=1, max_cluster_size=2)
        res, err = run(P.run_clusters, gz(FAM, self.tmp / "fam.fa.gz"), self.tmp / "gz",
                       clusters=gz(FAM_CLUSTERS, self.tmp / "cl.tsv.gz"), jobs=1, max_cluster_size=2)
        self.assertEqual(res.exit_code, plain.exit_code, err)
        for name in ("dnds", "cluster_summary"):
            self.assertEqual(Path(res.outputs[name]).read_bytes(), Path(plain.outputs[name]).read_bytes(), name)

    def test_a_damaged_gz_is_one_error_line(self):
        bad = self.tmp / "cut.fa.gz"
        bad.write_bytes(gzip.compress(Q.read_bytes())[:200])
        res, err = run(P.run_pairs, bad, S, self.tmp / "o", pairs=self.pairs_file(PLAIN))
        self.assertEqual(res.exit_code, 1)
        self.assertEqual(len(err.strip().splitlines()), 1, err)
        self.assertIn("--query: cannot decompress", err)
        self.assertIn("damaged or cut short", err)
        self.assertEqual(res.report["status"], "failed")

    def test_command_line(self):
        out = self.tmp / "o"
        with redirect_stderr(io.StringIO()):
            rc = P.main(["pairs", "--query", str(gz(Q, self.tmp / "q.gz")), "--subject", str(S), "--pairs",
                         str(self.pairs_file(PLAIN)), "-o", str(out), "-q"])
        self.assertEqual(rc, 2)   # (the fixture has genes that cannot be used: exit 2, as for the plain file)
        self.assertGreater(len((out / "dnds.tsv").read_text(encoding="utf-8").splitlines()), 10)


@needs_cpp
class OsErrors(PipelineCase):
    """A full disk at different moments of a run, injected."""

    def full(self, name="x"):
        return OSError(errno.ENOSPC, os.strerror(errno.ENOSPC), str(self.tmp / "o" / name))

    def check_failed(self, res, err, expect_report=True):
        self.assertEqual(res.exit_code, 1, err)
        self.assertNotIn("Traceback", err)
        lines = [ln for ln in err.strip().splitlines() if ln.startswith("coal pipeline: error:")]
        self.assertEqual(len(lines), 1, err)
        self.assertRegex(lines[0], "(?i)no space left on device")
        self.assertIn("free space or use another -o folder", lines[0])
        rp = self.tmp / "o" / "report.json"
        if expect_report:
            rep = json.loads(rp.read_text(encoding="utf-8"))   # complete JSON, and it says the truth
            self.assertEqual((rep["exit_code"], rep["status"]), (1, "failed"))
            self.assertRegex(rep["error"], "(?i)no space left")
        self.assertFalse((self.tmp / "o.coal-lock").exists(), "the lock was left behind")
        self.assertFalse((self.tmp / "o" / "report.json.tmp").exists())

    def test_while_saving_the_stage_state(self):
        with mock.patch.object(P._Run, "_save_state", side_effect=self.full("work/state.json.tmp")):
            res, err = run(P.run_pairs, Q, S, self.tmp / "o", pairs=self.pairs_file(PLAIN))
        self.check_failed(res, err)
        self.assertIn("state.json.tmp", err)

    def test_while_writing_the_genes_table(self):
        with mock.patch.object(P._Run, "read_genes", side_effect=self.full("genes.tsv")):
            res, err = run(P.run_pairs, Q, S, self.tmp / "o", pairs=self.pairs_file(PLAIN))
        self.check_failed(res, err)

    def test_when_the_failures_table_cannot_be_flushed(self):
        """The scout's traceback: failed.tsv is closed in finish(), and a full disk shows there."""
        real_open = open

        class Full(io.TextIOWrapper):
            def close(self_inner):
                if not self_inner.closed:
                    try:
                        super().close()
                    finally:
                        pass
                    raise OSError(errno.ENOSPC, os.strerror(errno.ENOSPC), self_inner.name)

        def fake_open(path, *a, **k):
            if str(path).endswith("failed.tsv.tmp") and a and "w" in a[0]:   # (written as .tmp, named at the end)
                return Full(io.FileIO(path, "w"), encoding="utf-8", newline="\n")
            return real_open(path, *a, **k)
        with mock.patch("builtins.open", fake_open):
            res, err = run(P.run_pairs, Q, S, self.tmp / "o", pairs=self.pairs_file(PLAIN))
        self.assertEqual(res.exit_code, 1, err)
        self.assertNotIn("Traceback", err)
        self.assertRegex(err, r"coal pipeline: error: cannot write .*failed\.tsv: (?i:no space left)")
        rep = json.loads((self.tmp / "o" / "report.json").read_text(encoding="utf-8"))
        self.assertEqual((rep["exit_code"], rep["status"]), (1, "failed"))
        # a failures table that could not be completed is not there under its name (nor as .tmp)
        self.assertFalse((self.tmp / "o" / "failed.tsv").exists())
        self.assertFalse((self.tmp / "o" / "failed.tsv.tmp").exists())

    def test_when_the_report_cannot_be_written(self):
        real_dump = json.dump
        calls = []

        def dump(obj, fh, **k):
            if getattr(fh, "name", "").endswith("report.json.tmp"):
                calls.append(len(obj))
                if len(calls) == 1:   # the whole report does not fit; the short one does
                    raise OSError(errno.ENOSPC, os.strerror(errno.ENOSPC), fh.name)
            return real_dump(obj, fh, **k)
        with mock.patch.object(P.json, "dump", dump):
            res, err = run(P.run_pairs, Q, S, self.tmp / "o", pairs=self.pairs_file(PLAIN), quiet=True)
        # the fixture's unusable genes make this run exit 2, and it stays 2: the tables are complete
        self.assertEqual(res.exit_code, 2)
        self.assertRegex(err, r"could not write the report .*report\.json: (?i:no space left)")
        rep = json.loads((self.tmp / "o" / "report.json").read_text(encoding="utf-8"))
        self.assertEqual((rep["exit_code"], rep["status"], rep["report_truncated"]),
                         (2, "completed_with_failures", True))
        self.assertEqual(len(calls), 2)
        self.assertNotIn("counts", rep)

    def test_a_clean_run_whose_report_cannot_be_written_is_exit_1(self):
        q, s = self.tmp / "q.fa", self.tmp / "s.fa"
        q.write_text(">a\nATGAAAGTTCTGGCAGGTCCGTTA\n", encoding="utf-8")
        s.write_text(">b\nATGAAGGTTTTGGCAGGACCGCTA\n", encoding="utf-8")
        real_dump = json.dump

        def dump(obj, fh, **k):
            if getattr(fh, "name", "").endswith("report.json.tmp"):
                raise OSError(errno.ENOSPC, os.strerror(errno.ENOSPC), fh.name)
            return real_dump(obj, fh, **k)
        ok, _ = run(P.run_pairs, q, s, self.tmp / "ok", pairs=self.pairs_file([("a", "b")]))
        self.assertEqual(ok.exit_code, 0)
        with mock.patch.object(P.json, "dump", dump):
            res, err = run(P.run_pairs, q, s, self.tmp / "o", pairs=self.pairs_file([("a", "b")]), quiet=True)
        self.assertEqual(res.exit_code, 1)
        self.assertEqual((res.report["status"], res.report["exit_code"]), ("failed", 1))
        self.assertFalse((self.tmp / "o" / "report.json").exists())     # absent, not cut short
        self.assertFalse((self.tmp / "o" / "report.json.tmp").exists())

    def test_os_error_text(self):
        e = OSError(errno.ENOSPC, "No space left on device", "/out/work/x")
        self.assertEqual(P._os_error_text(e, "/out"),
                         "cannot write /out/work/x: No space left on device. The disk (or quota) of the output "
                         "folder is full: free space or use another -o folder; --resume reuses the finished stages")
        self.assertEqual(P._os_error_text(OSError(errno.EIO, "Input/output error"), "/out"),
                         "cannot write into the output folder /out: Input/output error")
        self.assertEqual(P._os_error_text(PermissionError(errno.EACCES, "Permission denied", "/in/q.fa"), "/out"),
                         "cannot read or write /in/q.fa: Permission denied")


@needs_cpp
class UnmatchedIds(PipelineCase):
    def test_pairs_with_other_ids_get_one_summary_line(self):
        pairs = self.pairs_file([(q + ".1", s) for q, s in PLAIN])   # a version suffix the file does not have
        res, err = run(P.run_pairs, Q, S, self.tmp / "o", pairs=pairs)
        self.assertEqual(res.exit_code, 2)
        lines = [ln for ln in err.splitlines() if "are not in the" in ln]
        self.assertEqual(len(lines), 1, err)
        self.assertIn("12 of the 12 query IDs of --pairs", lines[0])
        self.assertIn("'g01.1' in pairs.tsv", lines[0])
        self.assertIn("while the sequences are named like 'g01'", lines[0])
        self.assertIn("version or transcript suffix", lines[0])
        self.assertEqual(res.report["warnings"], [lines[0].split("warning: ", 1)[1]])
        self.assertEqual(res.report["counts"]["pairs"]["error"], 12)   # every pair still has its row

    def test_suffix_in_the_file_and_not_in_the_pairs(self):
        q = self.tmp / "q.fa"
        q.write_text("".join(f">{r.split()[0][1:]}.1\n" + r.split("\n", 1)[1] for r in
                             (">" + x for x in Q.read_text(encoding="utf-8").split(">")[1:])), encoding="utf-8")
        res, err = run(P.run_pairs, q, S, self.tmp / "o", pairs=self.pairs_file(PLAIN))
        self.assertIn("'g01' in pairs.tsv, while the sequences are named like 'g01.1'", err)
        self.assertIn("version or transcript suffix", err)

    def test_swapped_columns(self):
        res, err = run(P.run_pairs, Q, S, self.tmp / "o", pairs=self.pairs_file([(s, q) for q, s in PLAIN]))
        lines = [ln for ln in err.splitlines() if "are not in the" in ln]
        self.assertEqual(len(lines), 2, err)   # one per side
        self.assertIn("the columns are in the right order", lines[0])

    def test_a_few_missing_ids_get_no_summary_line(self):
        res, err = run(P.run_pairs, Q, S, self.tmp / "o", pairs=self.pairs_file(PLAIN + [("nobody", "h01")]))
        self.assertNotIn("are not in the", err)
        self.assertNotIn("warnings", res.report)

    def test_clusters(self):
        cl = self.tmp / "cl.tsv"
        cl.write_text("".join(f"{a}.t1\t{b}.t1\n" for a, b in
                              (ln.split("\t") for ln in FAM_CLUSTERS.read_text(encoding="utf-8").splitlines()
                               if ln and not ln.startswith(("#", "centroid", "representative")))), encoding="utf-8")
        res, err = run(P.run_clusters, FAM, self.tmp / "o", clusters=cl)
        lines = [ln for ln in err.splitlines() if "are not in the" in ln]
        self.assertEqual(len(lines), 1, err)
        self.assertIn("CDS IDs of --clusters", lines[0])


@needs_cpp
@unittest.skipUnless(POSIX, "the stand-in DIAMOND is a POSIX shell script")
class Tools(PipelineCase):
    def fake(self, body: str, version: str = "2.2.8") -> str:
        f = self.tmp / "diamond"
        f.write_text(f'#!/bin/sh\nif [ "$1" = version ]; then echo "diamond version {version}"; exit 0; fi\n'
                     f'echo "$@" >> "{self.tmp}/calls.txt"\n{body}\n', encoding="utf-8")
        f.chmod(0o755)
        return str(f)

    def test_old_diamond_stops_cluster_mode(self):
        res, err = run(P.run_clusters, FAM, self.tmp / "o", diamond=self.fake("exit 0", "2.0.14"))
        self.assertEqual(res.exit_code, 1)
        self.assertIn("DIAMOND 2.0.14", err)
        self.assertIn("`diamond deepclust` needs DIAMOND 2.1.0 or later (COAL is tested with 2.2.8)", err)
        self.assertIn("coal pipeline install diamond", err)
        self.assertFalse((self.tmp / "calls.txt").exists(), "the old DIAMOND was started")

    def test_old_diamond_with_given_clusters_is_not_needed(self):
        res, err = run(P.run_clusters, FAM, self.tmp / "o", clusters=FAM_CLUSTERS, max_cluster_size=2,
                       diamond=self.fake("exit 0", "2.0.14"))
        self.assertIn(res.exit_code, (0, 2), err)
        self.assertNotIn("too old", err)

    def test_old_diamond_warns_in_pair_mode(self):
        res, err = run(P.run_pairs, Q, S, self.tmp / "o", diamond=self.fake("exit 3", "v2.0.9"))
        self.assertIn("warning: DIAMOND v2.0.9", err)
        self.assertIn("older than 2.1.0", err)
        self.assertEqual(res.exit_code, 1)   # (the stand-in fails; a real 2.0 would go on)
        self.assertIn("diamond makedb failed (exit 3)", err)

    def test_version_tuple(self):
        self.assertEqual([P._version_tuple(v) for v in ("2.1.0", "v2.0.14", "2.1.10.160", "2.1.9-beta", "x")],
                         [(2, 1, 0), (2, 0, 14), (2, 1, 10, 160), (2, 1, 9), None])
        self.assertLess(P._version_tuple("2.0.15"), P.DIAMOND_MIN_VERSION)
        self.assertGreaterEqual(P._version_tuple("2.1"), (2, 1))

    def test_jobs_above_the_cpu_count_are_clamped_for_the_tools(self):
        cpus = os.cpu_count() or 1
        res, err = run(P.run_pairs, Q, S, self.tmp / "o", diamond=self.fake("exit 3"), jobs=cpus + 4999)
        self.assertIn(f"-j {cpus + 4999} is more than the {cpus} CPU(s) of this machine", err)
        calls = (self.tmp / "calls.txt").read_text(encoding="utf-8")
        self.assertIn(f"--threads {cpus} ", calls)
        self.assertNotIn(str(cpus + 4999), calls)
        # within the CPU count nothing is said and nothing changes
        res, err = run(P.run_pairs, Q, S, self.tmp / "o2", diamond=self.fake("exit 3"), jobs=1)
        self.assertNotIn("is more than", err)

    def test_a_tool_killed_by_a_signal_is_reported_as_that(self):
        res, err = run(P.run_pairs, Q, S, self.tmp / "o", diamond=self.fake("kill -9 $$"))
        self.assertEqual(res.exit_code, 1)
        self.assertIn("diamond makedb failed (killed by signal 9 (SIGKILL; out of memory?))", err)
        self.assertNotIn("exit -9", err)

    def test_exit_text(self):
        self.assertEqual(msa.exit_text(3), "exit 3")
        self.assertEqual(msa.exit_text(0), "exit 0")
        self.assertEqual(msa.exit_text(-9), "killed by signal 9 (SIGKILL; out of memory?)")
        self.assertEqual(msa.exit_text(-11), "killed by signal 11 (SIGSEGV; the program crashed)")
        self.assertEqual(msa.exit_text(-250), "killed by signal 250")


@needs_cpp
class SameFile(PipelineCase):
    def test_query_and_subject_the_same_file_warns(self):
        link = self.tmp / "also_q.fasta"
        try:
            os.symlink(Q, link)
        except OSError:
            link = Q
        pairs = self.pairs_file([(q, q) for q, _ in PLAIN])
        res, err = run(P.run_pairs, Q, link, self.tmp / "o", pairs=pairs)
        self.assertEqual(res.exit_code, 2, err)   # (2: the fixture's unusable genes)
        self.assertIn("warning: --query and --subject are the same file", err)
        self.assertIn("coal pipeline clusters", err)
        self.assertEqual(len(res.report["warnings"]), 1)
        res, err = run(P.run_pairs, Q, S, self.tmp / "o2", pairs=self.pairs_file(PLAIN))
        self.assertNotIn("same file", err)


@needs_diamond
class RealDiamond(PipelineCase):
    def test_the_installed_diamond_is_new_enough_and_gz_inputs_search(self):
        v = P._version_tuple(P._diamond_version(P.find_diamond()))
        self.assertGreaterEqual(v, P.DIAMOND_MIN_VERSION)
        plain, _ = run(P.run_pairs, Q, S, self.tmp / "plain", jobs=2)
        res, err = run(P.run_pairs, gz(Q, self.tmp / "q.fa.gz"), gz(S, self.tmp / "s.fa.gz"), self.tmp / "gz", jobs=2)
        self.assertEqual(res.exit_code, plain.exit_code, err)
        self.assertEqual(Path(res.outputs["dnds"]).read_bytes(), Path(plain.outputs["dnds"]).read_bytes())


if __name__ == "__main__":
    unittest.main()
