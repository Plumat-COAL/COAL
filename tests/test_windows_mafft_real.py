"""MAFFT's real Windows package, downloaded (COAL_TEST_REAL_MAFFT=1; 27 MB from mafft.cbrc.jp).

CI's C++ job sets COAL_TEST_REAL_MAFFT=1 on its Windows runner. There `coal pipeline install mafft` installs the
pinned package into a new cache folder, and every run here has the cache folder, the current folder and TMP / TEMP
in paths with spaces, `&`, `^` and `( )` (cmd.exe's special characters: COAL runs the package's bash.exe, not
mafft.bat, so they are only text). The tests check the record of the install, that a second install downloads
nothing, that MAFFT runs (--version, an alignment, the paths MSYS2 reads), and one `coal pipeline clusters
--aligner mafft` run on examples/dnds, end to end with coal-cpp (COAL_CPP_BIN), with the numbers of before. On another system the download is checked and unpacked the same way, with the
platform taken to be x86-64 Windows (nothing of it is run there).

A download that cannot be made (the site cannot be reached) is a skip that says so; a file that is not the pinned
one (another size or SHA-256: MAFFT replaced the file at that address) fails. Without COAL_TEST_REAL_MAFFT=1
everything here skips: tests/test_windows_tools.py tests the same code with the download mocked."""
from __future__ import annotations

import csv
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from coal import msa  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402

REAL = os.environ.get("COAL_TEST_REAL_MAFFT") == "1"
WINDOWS = os.name == "nt"
EX = ROOT / "examples" / "dnds"


def cpp_bin() -> Path | None:
    """coal-cpp as the C++ job builds it (COAL_CPP_BIN), as tests/test_dnds_pipeline.py finds it."""
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


def unreachable(text: str) -> bool:
    """The install failed because MAFFT's site could not be reached (not because the file was another one)."""
    return "download failed (" in text


@unittest.skipUnless(REAL, "COAL_TEST_REAL_MAFFT=1 is not set (downloads MAFFT's Windows package, 27 MB)")
class TestRealMafftWindowsPackage(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="coal real & ^ mafft ")
        cls.cache = os.path.join(cls.tmp, "cache & ^ (x) dir")
        cls.cwd = os.path.join(cls.tmp, "work & ^ (x)")          # the current folder of every run
        cls.temp = os.path.join(cls.tmp, "temp & ^ (x)")         # TMP / TEMP of every run
        for d in (cls.cwd, cls.temp):
            os.makedirs(d)
        cls.env = {**os.environ, "COAL_CACHE_DIR": cls.cache, "TMP": cls.temp, "TEMP": cls.temp,
                   "PYTHONPATH": str(ROOT / "src") + os.pathsep + os.environ.get("PYTHONPATH", "")}
        cls.env.pop("COAL_MAFFT", None)
        cls.first = cls.install()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    @classmethod
    def install(cls) -> tuple[int, str, str]:
        """`coal pipeline install mafft` (rc, stdout, stderr): a process of its own on Windows; elsewhere in this
        process, with the platform taken to be x86-64 Windows and the version not read (nothing runs)."""
        if WINDOWS:
            r = subprocess.run([sys.executable, "-m", "coal", "pipeline", "install", "mafft"], capture_output=True,
                               text=True, env=cls.env, timeout=900, stdin=subprocess.DEVNULL, encoding="utf-8",
                               errors="backslashreplace", cwd=cls.cwd)
            return r.returncode, r.stdout, r.stderr
        out, err = io.StringIO(), io.StringIO()
        with mock.patch.dict(os.environ, {"COAL_CACHE_DIR": cls.cache}), \
                mock.patch.object(msa, "_platform", return_value=("Windows", "AMD64")), \
                mock.patch.object(msa.platform, "system", return_value="Windows"), \
                mock.patch.object(msa, "tool_version", return_value=f"v{msa.MAFFT_VERSION} (2024/Apr/26)"), \
                redirect_stdout(out), redirect_stderr(err):
            rc = P.main(["install", "mafft"])
        return rc, out.getvalue(), err.getvalue()

    def setUp(self):
        rc, out, err = self.first
        if rc != 0 and unreachable(err):
            self.skipTest(f"MAFFT's site could not be reached: {err.strip()[-300:]}")
        self.assertEqual(rc, 0, err)
        self.entry = os.path.join(self.cache, "tools", "mafft", "mafft-win", "mafft.bat")

    def test_installed_and_recorded(self):
        rc, out, err = self.first
        self.assertEqual(out.strip(), self.entry)
        w = msa.MAFFT_WINDOWS
        self.assertIn(f"coal: MAFFT package: {w.url} (win-64, {w.filename}, SHA-256 {w.sha256}, the pinned build)",
                      err)
        rec = json.loads(Path(self.entry + msa.MANIFEST_SUFFIX).read_text(encoding="utf-8"))
        self.assertEqual(rec["package"]["sha256"], w.sha256)
        self.assertEqual(rec["package"]["size"], w.size)
        files = rec["files"]
        self.assertIn("mafft.bat", files)
        self.assertIn("usr/bin/mafft", files)
        self.assertTrue(any(n.startswith("usr/lib/mafft/") and n.endswith(".exe") for n in files), sorted(files))
        on_disk = sorted(str(p.relative_to(Path(self.entry).parent)).replace(os.sep, "/")
                         for p in Path(self.entry).parent.rglob("*") if p.is_file())
        self.assertEqual(on_disk, sorted([*files, "mafft.bat" + msa.MANIFEST_SUFFIX]))
        dirs = [p for p in Path(self.entry).parent.rglob("*") if p.is_dir()]
        # the zip's 114 entries: every file and (at most) every folder, under the one folder mafft-win/
        self.assertLessEqual(len(files), 114)
        self.assertLessEqual(len(files) + len(dirs) + 1, 114 + 1)
        print(f"\n[real MAFFT] {len(files)} files, {len(dirs)} folders under mafft-win/", file=sys.stderr)
        self.assertEqual(sorted(os.listdir(os.path.join(self.cache, "tools", "mafft"))), ["mafft-win"])

    def test_a_second_install_downloads_nothing(self):
        rc, out, err = self.install()
        self.assertEqual(rc, 0, err)
        self.assertEqual(out.strip(), self.entry)
        self.assertNotIn("downloading", err)

    @unittest.skipUnless(WINDOWS, "mafft.bat runs on Windows")
    def test_mafft_bat_runs(self):
        print("\n[real MAFFT] mafft.bat:\n" + Path(self.entry).read_text(encoding="utf-8", errors="replace"),
              file=sys.stderr)
        here = os.getcwd()
        os.chdir(self.cwd)
        self.addCleanup(os.chdir, here)
        recs = [("a", "MKVLAAGIVGLLLAQ*"), ("b", "MKVLSAGIVGLLAQ"), ("c", "MRVLAAGIIGLLLAQ"), ("d", "MKVAAGIVGLLLAQ")]
        with mock.patch.dict(os.environ, {"TMP": self.temp, "TEMP": self.temp}), \
                mock.patch.object(msa.tempfile, "tempdir", self.temp):
            self.assertEqual(msa.tool_version(self.entry), f"v{msa.MAFFT_VERSION} (2024/Apr/26)")
            rows = msa.run_msa(recs, "mafft", exe=self.entry, threads=2)
        self.assertEqual([r[0] for r in rows], ["a", "b", "c", "d"])
        self.assertEqual(len({len(r[1]) for r in rows}), 1)
        for (_, s), (_, row) in zip(recs, rows):
            self.assertEqual(row.replace("-", ""), s)
        print(f"\n[real MAFFT] {msa.tool_version(self.entry)}: " + " ".join(r for _, r in rows), file=sys.stderr)

    @unittest.skipUnless(WINDOWS, "MAFFT's package runs on Windows")
    def test_paths_msys2_reads(self):
        """Audit note 6: a drive path with '/' as an argument is read; a UNC path (\\\\localhost\\C$\\...) is not
        (MSYS2 turns its two leading backslashes into one: seen in CI, 2026-11-15), so COAL never passes a path of
        the user's (the input is given relative to the run's folder) and refuses a run folder on a network path."""
        work = tempfile.mkdtemp(prefix="paths & ^ ", dir=self.temp)
        inp = os.path.join(work, "in.fasta")
        Path(inp).write_text(">s0\nMKVLAAGIV\n>s1\nMKVLSAGIV\n>s2\nMRVAAGIV\n", encoding="utf-8", newline="\n")
        args, extra = msa._launch([self.entry, "--auto", "--amino", "--anysymbol", "--quiet",
                                   os.path.abspath(inp).replace("\\", "/")], work)
        extra.pop("_tree")
        r = subprocess.run(args, capture_output=True, text=True, timeout=300, **extra)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual([x[1:].strip() for x in r.stdout.splitlines() if x.startswith(">")], ["s0", "s1", "s2"],
                         r.stdout + r.stderr)
        drive, rest = os.path.splitdrive(os.path.abspath(self.temp))
        unc = "\\\\localhost\\" + drive[0] + "$" + rest
        if not os.path.isdir(unc):
            self.skipTest(f"the administrative share of {drive} is not reachable here ({unc})")
        with mock.patch.object(msa.tempfile, "tempdir", unc), self.assertRaises(msa.MsaError) as cm:
            msa.run_msa([("a", "MKV"), ("b", "MKI"), ("c", "MRV")], "mafft", exe=self.entry)
        self.assertIn("on a network path", str(cm.exception))
        print(f"\n[real MAFFT] drive path with '/': read; run folder {unc}: refused", file=sys.stderr)

    @unittest.skipUnless(WINDOWS, "mafft.bat runs on Windows")
    def test_pipeline_clusters_with_mafft(self):
        """`coal pipeline clusters --aligner mafft` on examples/dnds: three clusters of four genes (two ortholog
        pairs each), so every cluster goes to MAFFT, found in COAL's cache."""
        cpp = cpp_bin()
        if cpp is None:
            if os.environ.get("COAL_REQUIRE_CPP"):
                self.fail("build/coal-cpp (COAL_CPP_BIN) is not built (COAL_REQUIRE_CPP is set)")
            self.skipTest("build/coal-cpp (COAL_CPP_BIN) is not built")
        work = Path(self.tmp, "pipeline run")
        work.mkdir(exist_ok=True)
        tsv = work / "clusters.tsv"
        groups = [("01", "02"), ("03", "04"), ("05", "06")]
        with open(tsv, "w", encoding="utf-8", newline="\n") as fh:
            for g, h in groups:
                for m in (f"spA_g{g}", f"spB_g{g}", f"spA_g{h}", f"spB_g{h}"):
                    fh.write(f"spA_g{g}\t{m}\n")
        out = Path(self.cwd) / "out & ^"
        env = {**self.env, "COAL_CPP": str(cpp)}
        env.pop("COAL_CLUSTALO", None)
        r = subprocess.run([sys.executable, "-m", "coal", "pipeline", "clusters", "--cds", str(EX / "species_a_cds.fasta"),
                            str(EX / "species_b_cds.fasta"), "--clusters", str(tsv), "--aligner", "mafft", "-o",
                            str(out)], capture_output=True, text=True, env=env, timeout=900,
                           stdin=subprocess.DEVNULL, encoding="utf-8", errors="backslashreplace", cwd=self.cwd)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out / "cluster_summary.tsv", encoding="utf-8") as fh:
            summary = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(sorted(row["cluster_id"] for row in summary), sorted(f"spA_g{g}" for g, _ in groups))
        for row in summary:
            self.assertEqual(row["aligner"], "mafft", row)
            self.assertEqual(row["n_members"], "4", row)
        with open(out / "dnds.tsv", encoding="utf-8") as fh:
            pairs = list(csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"))
        got = {tuple(sorted((p["seq1"], p["seq2"]))) for p in pairs}
        self.assertEqual(len(got), 3 * 6)        # every pair of each cluster of four
        for g, h in groups:                      # the ortholog pairs have estimates (the others may be saturated)
            for k in (g, h):
                rows = [p for p in pairs if {p["seq1"], p["seq2"]} == {f"spA_g{k}", f"spB_g{k}"}]
                self.assertTrue(rows and all(p["status"] == "ok" for p in rows), rows)
        g01 = [p for p in pairs if {p["seq1"], p["seq2"]} == {"spA_g01", "spB_g01"}]
        self.assertEqual([p["dNdS"] for p in g01], ["0.025483147"])     # as in the run through mafft.bat
        report = json.loads((out / "report.json").read_text(encoding="utf-8"))
        tool = json.dumps(report.get("tools", report))
        self.assertIn("mafft.bat", tool)
        self.assertIn(msa.MAFFT_VERSION, tool)
        print(f"\n[real MAFFT] pipeline clusters --aligner mafft: {len(pairs)} pairs, "
              + ", ".join(f"{p['seq1']}-{p['seq2']} dN/dS {p.get('dNdS')}" for p in pairs[:3]), file=sys.stderr)


if __name__ == "__main__":
    unittest.main()
