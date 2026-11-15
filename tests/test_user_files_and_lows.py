"""The low findings of the final review of 2026-11-15 (REVIEW-COAL-2 L1, L2, L3, L5, L6; AUDIT-91's Windows MAFFT
notes) and the small texts that came with them.

  * `coal pipeline -o DIR`: a `.fasta` of the user's in DIR/codon_alignments/ (or protein_alignments/) is kept when
    the next run clears the earlier run's files; a run that writes that folder again refuses (exit 1) instead of
    possibly replacing it. A DIR/work/ that no COAL run made is not written into (exit 1, unless --overwrite), and
    `--report FILE` onto an existing file that is not a COAL report is refused the same way.
  * `coal-cpp kaks`: the gzip refusal names `coal kaks` / `coal dnds`, the Python commands that decompress.
  * `--digits N` / digits=N reaches dnds.kaks (codon-dnds --kaks-digits); without it dnds.kaks keeps
    KaKs_Calculator's 6 digits, as before. digits out of 1..17 is refused from Python as on the command line.
  * coal.msa on Windows (mocked): a user's batch file found on PATH gives way to the pinned package in COAL's cache;
    COAL's own current folder is no reason to refuse a batch file of the user's; MAFFT's package in a folder
    with ';' is refused, one with non-ASCII letters is run as any other.
  * CITATION.cff, docs/python-api.md, coalr's NEWS.md, scripts/release_notes.py, CHANGELOG.md.
Needs build/coal-cpp (COAL_CPP_BIN) for the pipeline and engine parts.
"""
from __future__ import annotations

import contextlib
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import __version__, msa  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402
from test_dnds_pipeline import CPP, FAM, PLAIN, Q, S, PipelineCase, needs_cpp, read_fasta  # noqa: E402


# a cluster run on the fixture ends with exit 2: one of its records is not usable (listed in failed.tsv)
OK = 2


def run(fn, *args, **kw):
    """(result, stderr text) of a pipeline call."""
    with redirect_stderr(io.StringIO()) as err:
        res = fn(*args, **kw)
    return res, err.getvalue()


def env() -> dict:
    return {**os.environ, "PYTHONPATH": str(ROOT / "src")}


def cli(*args, timeout=300):
    r = subprocess.run([sys.executable, "-m", "coal", "pipeline", *map(str, args)], capture_output=True, text=True,
                       encoding="utf-8", timeout=timeout, stdin=subprocess.DEVNULL, env=env(), check=False)
    return r.returncode, r.stderr.replace("\r\n", "\n")


def snapshot(folder: Path) -> dict:
    return {str(p.relative_to(folder)): (p.read_bytes() if p.is_file() else None) for p in sorted(folder.rglob("*"))}


@needs_cpp
class UserFilesInTheOutputFolder(PipelineCase):
    def pairs_inputs(self):
        q, s = self.tmp / "q.fa", self.tmp / "s.fa"
        allq, alls = read_fasta(Q), read_fasta(S)
        q.write_text("".join(f">{a}\n{allq[a]}\n" for a, _ in PLAIN[:4]), encoding="utf-8")
        s.write_text("".join(f">{b}\n{alls[b]}\n" for _, b in PLAIN[:4]), encoding="utf-8")
        pf = self.tmp / "pairs.tsv"
        pf.write_text("".join(f"{a}\t{b}\n" for a, b in PLAIN[:4]), encoding="utf-8")
        return q, s, pf

    def two_member_clusters(self) -> Path:
        """Clusters of two (aligned by NW: no Clustal Omega needed)."""
        ids = list(read_fasta(FAM))
        tsv = self.tmp / "clusters.tsv"
        tsv.write_text("".join(f"{ids[k]}\t{ids[k]}\n{ids[k]}\t{ids[k + 1]}\n" for k in range(0, 8, 2)),
                       encoding="utf-8")
        return tsv

    # -- L3: codon_alignments/ and protein_alignments/
    def test_a_users_fasta_in_codon_alignments_is_kept(self):
        tsv = self.two_member_clusters()
        out = self.tmp / "c"
        first, err = run(P.run_clusters, FAM, out, clusters=tsv, keep_alignments=True, keep_protein_alignments=True)
        self.assertEqual(first.exit_code, OK, err)
        written = sorted(p.name for p in (out / "codon_alignments").iterdir())
        self.assertEqual(len(written), 4)
        rec = json.loads((out / "work" / "outputs.json").read_text(encoding="utf-8"))
        self.assertEqual(rec["codon_alignments_files"], written)
        self.assertEqual(rec["protein_alignments_files"], written)
        mine = out / "codon_alignments" / "mine.fasta"
        mine.write_text(">mine\nATGAAA\n", encoding="utf-8")
        (out / "protein_alignments" / "mine.fasta").write_text(">mine\nMK\n", encoding="utf-8")
        # a run that does not write the folders: the earlier run's files go, the user's stay
        res, err = run(P.run_clusters, FAM, out, clusters=tsv)
        self.assertEqual(res.exit_code, OK, err)
        self.assertEqual(sorted(p.name for p in (out / "codon_alignments").iterdir()), ["mine.fasta"])
        self.assertEqual(mine.read_text(encoding="utf-8"), ">mine\nATGAAA\n")
        self.assertEqual(sorted(p.name for p in (out / "protein_alignments").iterdir()), ["mine.fasta"])

    def test_a_run_that_writes_the_folder_again_refuses(self):
        tsv = self.two_member_clusters()
        out = self.tmp / "c"
        res, err = run(P.run_clusters, FAM, out, clusters=tsv, keep_alignments=True)
        self.assertEqual(res.exit_code, OK, err)
        (out / "codon_alignments" / "mine.fasta").write_text(">mine\nATGAAA\n", encoding="utf-8")
        before = snapshot(out)
        for kw in ({}, {"overwrite": True}):
            res, err = run(P.run_clusters, FAM, out, clusters=tsv, keep_alignments=True, **kw)
            self.assertEqual(res.exit_code, 1, (kw, err))
            self.assertIn("holds mine.fasta, which no COAL run wrote there", err)
            self.assertIn("Nothing was changed", err)
            self.assertEqual(snapshot(out), before, kw)
        code, err = cli("clusters", "--cds", FAM, "--clusters", tsv, "-o", out, "--keep-alignments")
        self.assertEqual(code, 1, err)
        self.assertEqual(snapshot(out), before)
        # once the user's file is moved away the run goes on as before
        (out / "codon_alignments" / "mine.fasta").unlink()
        res, err = run(P.run_clusters, FAM, out, clusters=tsv, keep_alignments=True)
        self.assertEqual(res.exit_code, OK, err)

    def test_a_folder_recorded_without_a_list_is_cleared_as_before(self):
        """0.6.0's work/outputs.json (or a run killed while it wrote the folder) has no list: every .fasta there is
        the earlier run's, as before."""
        tsv = self.two_member_clusters()
        out = self.tmp / "c"
        res, err = run(P.run_clusters, FAM, out, clusters=tsv, keep_alignments=True)
        self.assertEqual(res.exit_code, OK, err)
        rp = out / "work" / "outputs.json"
        rec = json.loads(rp.read_text(encoding="utf-8"))
        del rec["codon_alignments_files"]
        rp.write_text(json.dumps(rec), encoding="utf-8")
        (out / "codon_alignments" / "other.fasta").write_text(">x\nATG\n", encoding="utf-8")
        res, err = run(P.run_clusters, FAM, out, clusters=tsv)
        self.assertEqual(res.exit_code, OK, err)
        self.assertFalse((out / "codon_alignments").exists())

    def test_normal_reruns_are_unchanged(self):
        tsv = self.two_member_clusters()
        ref = self.tmp / "ref"
        res, err = run(P.run_clusters, FAM, ref, clusters=tsv, keep_alignments=True)
        self.assertEqual(res.exit_code, OK, err)
        out = self.tmp / "c"
        for kw in ({"keep_alignments": True}, {"keep_alignments": True}, {"keep_alignments": True, "resume": True},
                   {}, {"keep_alignments": True}):
            res, err = run(P.run_clusters, FAM, out, clusters=tsv, **kw)
            self.assertEqual(res.exit_code, OK, (kw, err))
            self.assertNotIn("error:", err)
        for name in ("dnds.tsv", "cluster_summary.tsv"):
            self.assertEqual((out / name).read_bytes(), (ref / name).read_bytes(), name)
        self.assertEqual(snapshot(out / "codon_alignments"), snapshot(ref / "codon_alignments"))

    # -- L2: work/
    def test_a_users_work_folder_is_not_written_into(self):
        q, s, pf = self.pairs_inputs()
        out = self.tmp / "o"
        (out / "work").mkdir(parents=True)
        (out / "work" / "notes.txt").write_text("mine\n", encoding="utf-8")
        before = snapshot(out)
        res, err = run(P.run_pairs, q, s, out, pairs=pf)
        self.assertEqual(res.exit_code, 1, err)
        self.assertIn("work exists and no COAL run made it", err)
        self.assertIn("--overwrite", err)
        self.assertEqual(snapshot(out), before)
        self.assertFalse(Path(str(out) + ".coal-lock").exists())
        code, err = cli("pairs", "--query", q, "--subject", s, "--pairs", pf, "-o", out)
        self.assertEqual((code, err.count("error:")), (1, 1), err)
        self.assertEqual(snapshot(out), before)
        res, err = run(P.run_pairs, q, s, out, pairs=pf, overwrite=True)
        self.assertEqual(res.exit_code, 0, err)
        self.assertEqual((out / "work" / "notes.txt").read_text(encoding="utf-8"), "mine\n")
        res, err = run(P.run_pairs, q, s, out, pairs=pf)      # now COAL's (work/outputs.json): as before
        self.assertEqual(res.exit_code, 0, err)

    def test_work_cases_that_go_on(self):
        q, s, pf = self.pairs_inputs()
        out = self.tmp / "empty"
        (out / "work").mkdir(parents=True)                       # an empty work/: nothing to mix with
        res, err = run(P.run_pairs, q, s, out, pairs=pf)
        self.assertEqual(res.exit_code, 0, err)
        # a COAL work/ of a version without work/outputs.json: the folder's COAL report.json says so
        (out / "work" / "outputs.json").unlink()
        res, err = run(P.run_pairs, q, s, out, pairs=pf)
        self.assertEqual(res.exit_code, 0, err)
        # a file named work: never removed, not even with --overwrite
        out2 = self.tmp / "file"
        out2.mkdir()
        (out2 / "work").write_text("x\n", encoding="utf-8")
        for kw in ({}, {"overwrite": True}):
            res, err = run(P.run_pairs, q, s, out2, pairs=pf, **kw)
            self.assertEqual(res.exit_code, 1, err)
            self.assertIn("is not a folder", err)
        self.assertEqual(sorted(p.name for p in out2.iterdir()), ["work"])

    # -- L1: --report FILE
    def test_report_onto_a_users_file_is_refused(self):
        q, s, pf = self.pairs_inputs()
        for out, where in ((self.tmp / "o0", self.tmp / "elsewhere" / "notes.json"),
                           (self.tmp / "o1", self.tmp / "o1" / "mine.json")):
            where.parent.mkdir(parents=True, exist_ok=True)
            where.write_text('{"my": "data"}\n', encoding="utf-8")
            res, err = run(P.run_pairs, q, s, out, pairs=pf, report=str(where))
            self.assertEqual(res.exit_code, 1, err)
            self.assertIn(f"--report {where} already exists and is not a COAL report", err)
            self.assertEqual(where.read_text(encoding="utf-8"), '{"my": "data"}\n')
            self.assertEqual(sorted(p.name for p in out.iterdir()) if out.exists() else [],
                             [where.name] if where.parent == out else [])
            code, err = cli("pairs", "--query", q, "--subject", s, "--pairs", pf, "-o", out, "--report", where)
            self.assertEqual(code, 1, err)
            self.assertEqual(where.read_text(encoding="utf-8"), '{"my": "data"}\n')
            res, err = run(P.run_pairs, q, s, out, pairs=pf, report=str(where), overwrite=True)
            self.assertEqual(res.exit_code, 0, err)
            self.assertEqual(json.loads(where.read_text(encoding="utf-8"))["schema"], P.SCHEMA)

    def test_report_onto_a_coal_report_goes_on(self):
        q, s, pf = self.pairs_inputs()
        rp = self.tmp / "reports" / "latest.json"
        rp.parent.mkdir()
        for out in (self.tmp / "a", self.tmp / "b", self.tmp / "a"):   # the same report file for several folders
            res, err = run(P.run_pairs, q, s, out, pairs=pf, report=str(rp))
            self.assertEqual(res.exit_code, 0, err)
            self.assertEqual(json.loads(rp.read_text(encoding="utf-8"))["out_dir"], str(out))
        # a `coal --report` report (coal.report/1) is a COAL report too
        rp.write_text(json.dumps({"schema": "coal.report/1", "exit_code": 0}), encoding="utf-8")
        res, err = run(P.run_pairs, q, s, self.tmp / "c", pairs=pf, report=str(rp))
        self.assertEqual(res.exit_code, 0, err)
        self.assertTrue(P._is_coal_report(str(rp)))
        for text in ("", "[]", '{"schema": "other/1"}', '{"schema": "coal.pipeline/1 x"}', "not json"):
            rp.write_text(text, encoding="utf-8")
            self.assertFalse(P._is_coal_report(str(rp)), text)

    # -- L6: --digits and dnds.kaks
    def test_digits_reach_dnds_kaks(self):
        q, s, pf = self.pairs_inputs()

        def numbers(path: Path) -> list[str]:
            rows = path.read_text(encoding="utf-8").splitlines()[1:]
            return [f for row in rows for f in row.split("\t")[2:5] if re.fullmatch(r"[0-9.e+-]+", f)]

        def digits(x: str) -> int:
            return len(x.split("e")[0].replace(".", "").replace("-", "").lstrip("0"))

        base, err = run(P.run_pairs, q, s, self.tmp / "d0", pairs=pf, kaks_format=True)
        self.assertEqual(base.exit_code, 0, err)
        default = numbers(self.tmp / "d0" / "dnds.kaks")
        self.assertTrue(default and max(map(digits, default)) == 6, default)       # KaKs_Calculator's, as before
        engine = [c for st in base.report["stages"] for c in st.get("commands", ()) if "codon-dnds" in c]
        self.assertTrue(engine and "--kaks-out" in engine[0] and "--kaks-digits" not in engine[0], engine)
        res, err = run(P.run_pairs, q, s, self.tmp / "d4", pairs=pf, kaks_format=True, digits=4)
        self.assertEqual(res.exit_code, 0, err)
        four = numbers(self.tmp / "d4" / "dnds.kaks")
        self.assertTrue(four and max(map(digits, four)) <= 4 and len(four) == len(default), four)
        self.assertEqual(res.report["parameters"]["digits"], 4)
        self.assertIn("--digits 4", " ".join(res.report["command"]))
        code, err = cli("pairs", "--query", q, "--subject", s, "--pairs", pf, "-o", self.tmp / "c4", "--kaks-format",
                        "--digits", "4", "-q")
        self.assertEqual(code, 0, err)
        self.assertEqual((self.tmp / "c4" / "dnds.kaks").read_bytes(), (self.tmp / "d4" / "dnds.kaks").read_bytes())
        self.assertEqual((self.tmp / "c4" / "dnds.tsv").read_bytes(), (self.tmp / "d4" / "dnds.tsv").read_bytes())
        # --digits 8 given: 8 in dnds.kaks too; dnds.tsv as without it
        res, err = run(P.run_pairs, q, s, self.tmp / "d8", pairs=pf, kaks_format=True, digits=8)
        self.assertEqual((self.tmp / "d8" / "dnds.tsv").read_bytes(), (self.tmp / "d0" / "dnds.tsv").read_bytes())
        self.assertEqual(base.report["parameters"]["digits"], 8)
        self.assertNotIn("_digits_given", base.report["parameters"])
        for bad in (0, 18, True, 2.5):
            with self.subTest(bad=bad), self.assertRaises(ValueError) as cm:
                P.run_pairs(q, s, self.tmp / "bad", pairs=pf, digits=bad)
            self.assertIn("--digits: expected 1 to 17", str(cm.exception))
        self.assertFalse((self.tmp / "bad").exists())

    def test_codon_dnds_kaks_digits(self):
        groups = "G\tg1\tall\nM\ta\tATGAAAGTTCTGGCAGGTCCGTTA\nM\tb\tATGAAGGTTTTGGCAGGACCGCTA\nE\n"
        outs = {}
        for extra in ([], ["--kaks-digits", "3"]):
            k = self.tmp / f"k{len(extra)}.kaks"
            r = subprocess.run([str(CPP), "codon-dnds", "-i", "-", "-o", str(self.tmp / "o.tsv"), "--kaks-out", str(k),
                                *extra], input=groups, capture_output=True, text=True, check=False)
            self.assertEqual(r.returncode, 0, r.stderr)
            outs[len(extra)] = k.read_text(encoding="utf-8").splitlines()[1].split("\t")
        self.assertNotEqual(outs[0], outs[2])
        for f in outs[2][2:5]:        # Ka, Ks, Ka/Ks: at most 3 significant digits
            self.assertLessEqual(len(f.split("e")[0].replace(".", "").lstrip("0")), 3, outs[2])
        self.assertEqual(max(len(f.split("e")[0].replace(".", "").lstrip("0")) for f in outs[0][2:5]), 6, outs[0])
        r = subprocess.run([str(CPP), "codon-dnds", "-o", str(self.tmp / "o.tsv"), "--kaks-digits", "0"],
                           input="", capture_output=True, text=True, check=False)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("--kaks-digits: expected 1 to 17", r.stderr)


@needs_cpp
class KaksGzipMessage(unittest.TestCase):
    def test_names_the_commands_that_decompress(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            gz = Path(td, "aln.fa.gz")
            gz.write_bytes(gzip.compress(b">a\nATGAAA\n>b\nATGAAG\n"))
            for sub, args in (("kaks", ["-i", gz, "-o", Path(td, "o")]), ("dnds", [gz])):
                r = subprocess.run([str(CPP), sub, *map(str, args)], capture_output=True, text=True, check=False)
                self.assertEqual(r.returncode, 1, r.stderr)
                self.assertIn("is gzip-compressed: decompress it first (gunzip -k", r.stderr)
                self.assertIn("`coal kaks` and `coal dnds` decompress .gz inputs themselves", r.stderr)
                self.assertNotIn("`coal pipeline`", r.stderr)
                self.assertNotIn("the `coal` command", r.stderr)


class WindowsMafft(unittest.TestCase):
    """coal.msa on Windows, mocked (runs on every system)."""

    def setUp(self):
        import tempfile
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        msa._path_bat_told.clear()

    def package(self, where: str) -> str:
        for f in ("mafft.bat", "usr/bin/bash.exe", "usr/bin/mafft"):
            p = Path(where, "mafft-win", f)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(b"")
        return str(Path(where, "mafft-win", "mafft.bat"))

    @contextlib.contextmanager
    def windows(self, on_path: str | None, pinned_ok: bool, coal_mafft: str | None = None):
        pinned = os.path.join(self.tmp, "cache", "tools", "mafft", "mafft-win", "mafft.bat")
        env = {"COAL_CACHE_DIR": os.path.join(self.tmp, "cache")}
        with mock.patch.dict(os.environ, env), mock.patch.object(msa, "_platform", return_value=("Windows", "AMD64")), \
                mock.patch.object(msa.shutil, "which", return_value=on_path), \
                mock.patch.object(msa, "_trusted_cached", lambda tool, exe: pinned_ok and exe == pinned):
            if coal_mafft is None:
                os.environ.pop("COAL_MAFFT", None)
            else:
                os.environ["COAL_MAFFT"] = coal_mafft
            yield pinned

    def test_a_users_batch_file_on_path_gives_way_to_the_pinned_package(self):
        mine = os.path.join(self.tmp, "bin", "mafft.bat")
        Path(mine).parent.mkdir()
        Path(mine).write_text("@echo off\r\n", encoding="utf-8")
        with self.windows(mine, True) as pinned, redirect_stderr(io.StringIO()) as err:
            self.assertEqual(msa.find_mafft(), pinned)
            self.assertEqual(msa.find_mafft(), pinned)
        self.assertEqual(err.getvalue().count("coal: note:"), 1, err.getvalue())   # once
        self.assertIn(f"{mine} (mafft on PATH) is a batch file of your own", err.getvalue())
        self.assertIn("set COAL_MAFFT to it", err.getvalue())

    def test_when_the_path_one_is_kept(self):
        mine = os.path.join(self.tmp, "bin", "mafft.bat")
        Path(mine).parent.mkdir()
        Path(mine).write_text("@echo off\r\n", encoding="utf-8")
        pkg = self.package(os.path.join(self.tmp, "own"))
        for on_path, ok, named, want in ((mine, False, None, mine),        # no pinned install: the PATH one
                                         (pkg, True, None, pkg),           # MAFFT's own package (no cmd.exe)
                                         (None, True, mine, mine),         # COAL_MAFFT: always the user's
                                         (os.path.join(self.tmp, "mafft.exe"), True, None,
                                          os.path.join(self.tmp, "mafft.exe"))):
            Path(os.path.join(self.tmp, "mafft.exe")).write_bytes(b"")
            with self.subTest(on_path=on_path, named=named), self.windows(on_path, ok, named), \
                    redirect_stderr(io.StringIO()) as err:
                self.assertEqual(msa.find_mafft(), want)
            self.assertNotIn("note:", err.getvalue())
        # not Windows: PATH first, as for every tool
        with mock.patch.object(msa.shutil, "which", return_value=mine), mock.patch.dict(os.environ, {}), \
                mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            os.environ.pop("COAL_MAFFT", None)
            self.assertEqual(msa.find_mafft(), mine)

    def test_a_package_folder_with_a_semicolon_is_refused(self):
        bat = self.package(os.path.join(self.tmp, "a;b"))
        work = os.path.join(self.tmp, "run")
        os.makedirs(work)
        with self.assertRaises(msa.MsaError) as cm:
            msa._launch([bat, "--auto", os.path.join(work, "in.fasta")], work)
        self.assertIn("has a ';' in its path", str(cm.exception))
        self.assertIn("COAL_CACHE_DIR", str(cm.exception))
        # a ';' in the run's folder is only a current folder and a temporary folder (not in PATH): fine
        odd = os.path.join(self.tmp, "r;un")
        os.makedirs(odd)
        _, extra = msa._launch([self.package(os.path.join(self.tmp, "ok")), "--auto"], odd)
        self.assertEqual(extra["cwd"], odd)

    def test_a_package_folder_with_non_ascii_letters_runs(self):
        where = os.path.join(self.tmp, "Zo\u00eb \u6570\u636e caf\u00e9")
        bat = self.package(where)
        work = os.path.join(self.tmp, "r\u00fcn")
        os.makedirs(work)
        args, extra = msa._launch([bat, "--auto", "--quiet", os.path.join(work, "in.fasta")], work)
        usr_bin = os.path.join(where, "mafft-win", "usr", "bin")
        self.assertEqual(args, [os.path.join(usr_bin, "bash.exe"), "/usr/bin/mafft", "--auto", "--quiet", "in.fasta"])
        self.assertTrue(extra["env"]["PATH"].startswith(usr_bin + os.pathsep))
        self.assertEqual(extra["cwd"], work)

    def test_the_current_folder_does_not_matter_for_a_users_batch_file(self):
        odd = os.path.join(self.tmp, "R&D ^x")
        os.makedirs(odd)
        here = os.getcwd()
        os.chdir(odd)
        self.addCleanup(os.chdir, here)
        hidden = {**os.environ, "=C:": "C:\\R&D ^x", "=D:": "D:\\"}
        with mock.patch.object(msa.os, "environ", hidden):
            _, extra = msa._launch([os.path.join(self.tmp, "m.bat"), "--auto"], os.path.join(self.tmp, "run"))
        self.assertEqual(extra["cwd"], os.path.join(self.tmp, "run"))
        self.assertNotIn("=C:", extra["env"])
        self.assertNotIn("=D:", extra["env"])
        # the batch file's own folder, the run's folder and TMP / TEMP are still checked
        with self.assertRaises(msa.MsaError):
            msa._launch([os.path.join(odd, "m.bat")], os.path.join(self.tmp, "run"))

    def test_the_real_download_skip_is_visible(self):
        import test_windows_mafft_real as real
        summary = Path(self.tmp, "summary.md")
        real._SKIP_SHOWN.clear()
        self.addCleanup(real._SKIP_SHOWN.clear)
        out, err = io.StringIO(), io.StringIO()
        with mock.patch.dict(os.environ, {"GITHUB_ACTIONS": "true", "GITHUB_STEP_SUMMARY": str(summary)}), \
                contextlib.redirect_stdout(out), redirect_stderr(err):
            why = real.show_skip("MAFFT's site could not be reached: download failed (timed out)\n100%")
            real.show_skip("again")
        self.assertEqual(why, "MAFFT's site could not be reached: download failed (timed out) 100%")
        self.assertEqual(out.getvalue(), "::warning title=MAFFT real-download tests skipped::MAFFT's site could not "
                                         "be reached: download failed (timed out) 100%25\n")
        self.assertIn("SKIPPED, not run: MAFFT's site could not be reached", err.getvalue())
        self.assertEqual(summary.read_text(encoding="utf-8").count("MAFFT real-download tests skipped"), 1)


class Texts(unittest.TestCase):
    def read(self, rel: str) -> str:
        return (ROOT / rel).read_text(encoding="utf-8")

    def test_citation_has_the_release_date(self):
        cff = self.read("CITATION.cff")
        dev = re.match(r"^\d+(?:\.\d+)+\.dev\d+$", __version__)
        if not dev:
            self.assertRegex(cff, r"(?m)^date-released: \d{4}-\d{2}-\d{2}$")
        if __version__ == "0.6.0":
            self.assertIn("\ndate-released: 2026-11-15\n", cff)

    def test_no_stale_version_notes(self):
        self.assertNotRegex(self.read("docs/python-api.md"), r"run with COAL \d+\.\d+\.\d+\.dev\d")
        self.assertNotIn("notes: the development version number", self.read("packaging/r/coalr/NEWS.md"))
        doc = self.read("scripts/release_notes.py")
        self.assertNotIn("about 127,000", doc)
        sys.path.insert(0, str(ROOT / "scripts"))
        import release_notes
        n = len(release_notes.changelog_section("0.6.0") or "")
        self.assertGreater(n, release_notes.MAX_BODY)      # the docstring's "0.6.0's does not [fit]"

    def test_changelog_unreleased_names_the_fixes(self):
        log = self.read("CHANGELOG.md")
        self.assertLess(log.index("## Unreleased"), log.index("## 0.6.0"))
        head = log[log.index("## Unreleased"):log.index("## 0.6.0")]
        for words in ("codon_alignments/", "work/", "--report", "--digits", "gzip", "mafft.bat", "date-released"):
            self.assertIn(words, head, words)


if __name__ == "__main__":
    unittest.main()
