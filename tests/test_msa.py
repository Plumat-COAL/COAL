"""Cluster MSAs for the codon step: src/coal/msa.py (docs/alignment.md).

  * Stop handling without any tool: trailing '*' becomes a last column of its own, internal '*' is aligned as X
    and put back; residues and order are kept.
  * With Clustal Omega / MAFFT found (PATH, COAL_CLUSTALO / COAL_MAFFT, or COAL's cache): run_msa and
    align_clusters on real tools (skipped otherwise).
  * align_clusters isolates a failing cluster, aligns single sequences without a tool and 2-sequence clusters with
    the pairwise aligner (pairs_with="nw").
  * install_clustalo: the cached copy is reused; the official-binary route (download faked, nothing fetched)
    extracts the Windows zip layout safely and checks the version; a binary that does not run is refused.
    The real download is not run here (network); docs/alignment.md records a manual run.
"""
from __future__ import annotations

import io
import os
import shutil
import stat
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from coal import msa  # noqa: E402

RECS = [("a", "MKVLAAGIVALLLW*"), ("b", "MKVLGIVAL*LLW"), ("c", "MKVAAGIVLLLW")]


def fake_tool(folder: str, name: str, script: str) -> str:
    """A POSIX shell stand-in for a tool."""
    p = os.path.join(folder, name)
    with open(p, "w", newline="\n") as fh:
        fh.write("#!/bin/sh\n" + script + "\n")
    os.chmod(p, os.stat(p).st_mode | stat.S_IXUSR)
    return p


class Stops(unittest.TestCase):
    def test_prepare_restore(self):
        seqs, tails, internal = msa._prepare(RECS)
        self.assertEqual(seqs, ["MKVLAAGIVALLLW", "MKVLGIVALXLLW", "MKVAAGIVLLLW"])
        self.assertEqual(tails, [True, False, False])
        self.assertEqual(internal, [[], [9], []])
        rows = ["MKVLAAGIVALLLW", "MKVL-GIVALXLLW", "MKV-AAGIVLLLW-"]
        out = msa._restore(rows, seqs, tails, internal)
        self.assertEqual(out, ["MKVLAAGIVALLLW*", "MKVL-GIVAL*LLW-", "MKV-AAGIVLLLW--"])
        for (_, s), r in zip(RECS, out):
            self.assertEqual(r.replace("-", ""), s)

    def test_no_stops_no_extra_column(self):
        seqs, tails, internal = msa._prepare([("a", "MK"), ("b", "MR")])
        self.assertEqual(msa._restore(["MK", "MR"], seqs, tails, internal), ["MK", "MR"])

    def test_single_sequence_needs_no_tool(self):
        with mock.patch.dict(os.environ, {"COAL_CLUSTALO": "/nonexistent"}):
            self.assertEqual(msa.run_msa([("x", "MKV*")]), [("x", "MKV*")])
            self.assertEqual(msa.run_msa([]), [])

    def test_missing_tool(self):
        with mock.patch.dict(os.environ, {"COAL_CLUSTALO": "/nonexistent"}):
            with self.assertRaises(msa.MsaError) as cm:
                msa.run_msa(RECS)
            self.assertIn("install_clustalo", str(cm.exception))


class Clusters(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    @unittest.skipIf(os.name == "nt", "POSIX shell stand-in")
    def test_failing_cluster_is_isolated_and_pairs_use_nw(self):
        # stand-in "clustalo": copies its input (all sequences of equal length align trivially), fails on 'FAIL'
        tool = fake_tool(self.tmp, "clustalo", 'in=""; out=""\nwhile [ $# -gt 0 ]; do case "$1" in -i) in="$2"; '
                         'shift;; -o) out="$2"; shift;; esac; shift; done\n'
                         'if grep -q FAIL "$in"; then echo "boom" >&2; exit 3; fi\ncp "$in" "$out"')
        clusters = {"ok": [("x", "MKVL"), ("y", "MRVL"), ("z", "MKIL")], "bad": [("x", "FAIL"), ("y", "MKVL"),
                                                                                 ("z", "MKVL")],
                    "pair": [("p", "MKVLW*"), ("q", "MKW*")], "one": [("s", "MK")],
                    "empty": [("e", ""), ("f", "MK")]}
        res = msa.align_clusters(clusters, exe=tool, jobs=2)
        self.assertEqual(list(res), list(clusters))
        self.assertEqual(res["ok"].rows, [("x", "MKVL"), ("y", "MRVL"), ("z", "MKIL")])
        self.assertIn("boom", res["bad"].error)
        self.assertEqual(res["one"].rows, [("s", "MK")])
        self.assertIn("is empty", res["empty"].error)
        # the pair went to the pairwise aligner (BLOSUM62, needle defaults), stops restored
        rows = dict(res["pair"].rows)
        self.assertEqual(rows["p"].replace("-", ""), "MKVLW*")
        self.assertEqual(rows["q"].replace("-", ""), "MKW*")
        self.assertTrue(rows["p"].endswith("*") and rows["q"].endswith("*"))

    def test_real_tools(self):
        for tool, find in (("clustalo", msa.find_clustalo), ("mafft", msa.find_mafft)):
            exe = find()
            if exe is None:
                continue
            with self.subTest(tool=tool):
                rows = msa.run_msa(RECS, tool=tool, threads=2, exe=exe)
                self.assertEqual([r for r, _ in rows], ["a", "b", "c"])
                self.assertEqual(len({len(r) for _, r in rows}), 1)
                for (_, s), (_, r) in zip(RECS, rows):
                    self.assertEqual(r.replace("-", ""), s)
                self.assertTrue(rows[0][1].endswith("*") and rows[1][1].endswith("-"))
                res = msa.align_clusters([RECS, RECS[:2], [("x", "M" * 3000), ("y", "M" * 2900), ("z", "MK")]],
                                         tool=tool, jobs=2, pairs_with="msa", exe=exe)
                self.assertTrue(all(r.error is None for r in res.values()), res)
        if msa.find_clustalo() is None and msa.find_mafft() is None:
            self.skipTest("neither clustalo nor mafft found")


@unittest.skipIf(os.name == "nt", "POSIX shell stand-ins")
class Install(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        p = mock.patch.dict(os.environ, {"COAL_CACHE_DIR": os.path.join(self.tmp, "cache")})
        p.start()
        self.addCleanup(p.stop)

    def test_cache_dir_and_find(self):
        self.assertEqual(msa.cache_dir(), os.path.join(self.tmp, "cache"))
        with mock.patch.dict(os.environ, {"PATH": self.tmp}):
            os.environ.pop("COAL_CLUSTALO", None)
            self.assertIsNone(msa.find_clustalo())
            d = os.path.join(self.tmp, "cache", "tools", "clustalo")
            os.makedirs(d)
            exe = fake_tool(d, "clustalo", "echo 1.2.4")
            self.assertEqual(msa.find_clustalo(), exe)
            self.assertEqual(msa.install_clustalo(quiet=True), exe)  # cached: nothing installed

    def test_official_zip_route(self):
        payload = io.BytesIO()
        with zipfile.ZipFile(payload, "w") as zf:
            zf.writestr("clustal-omega-1.2.2-win64/clustalo.exe", "#!/bin/sh\necho 1.2.2\n")
            zf.writestr("clustal-omega-1.2.2-win64/libargtable2-0.dll", "dll")
            zf.writestr("clustal-omega-1.2.2-win64/../../evil.txt", "no")

        def fake_download(url, dest, quiet):
            self.assertTrue(url.endswith(".zip"))
            with open(dest, "wb") as fh:
                fh.write(payload.getvalue())
        with mock.patch.object(msa, "_download", fake_download), mock.patch.object(msa, "_EXE", ".exe"):
            exe = msa._install_official("Windows", "AMD64", True)
        d = os.path.join(self.tmp, "cache", "tools", "clustalo")
        self.assertEqual(exe, os.path.join(d, "clustalo.exe"))
        self.assertTrue(os.path.isfile(os.path.join(d, "libargtable2-0.dll")))
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "evil.txt")))
        self.assertEqual(msa.tool_version(exe), "1.2.2")

    def test_official_binary_that_does_not_run_is_refused(self):
        def fake_download(url, dest, quiet):
            with open(dest, "w") as fh:
                fh.write("#!/bin/sh\nexit 1\n")
        with mock.patch.object(msa, "_download", fake_download), \
                mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            with self.assertRaises(msa.MsaError) as cm:
                msa.install_clustalo(method="official", quiet=True)
        self.assertIn("is not Clustal Omega 1.2", str(cm.exception))
        self.assertIn("conda install -c bioconda clustalo", str(cm.exception))

    def test_auto_order(self):
        tried = []

        def fake_conda_install(tool, spec, conda, quiet, timeout):
            tried.append(("conda", conda, spec))
            raise msa.MsaError("offline")

        def fake_official(sysname, machine, quiet):
            tried.append(("official", sysname))
            raise msa.MsaError("offline")
        with mock.patch.object(msa, "_conda_install", fake_conda_install), \
                mock.patch.object(msa, "_install_official", fake_official), \
                mock.patch.object(msa, "_conda_exe", return_value="/opt/conda/bin/mamba"), \
                mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            with self.assertRaises(msa.MsaError):
                msa.install_clustalo(quiet=True)
        self.assertEqual(tried, [("conda", "/opt/conda/bin/mamba", "clustalo=1.2.4"), ("official", "Linux")])
        tried.clear()
        with mock.patch.object(msa, "_install_official", fake_official), \
                mock.patch.object(msa, "_platform", return_value=("Windows", "AMD64")):
            with self.assertRaises(msa.MsaError):
                msa.install_clustalo(quiet=True)
        self.assertEqual(tried, [("official", "Windows")])


if __name__ == "__main__":
    unittest.main()
