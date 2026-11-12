"""Cluster MSAs for the codon step: src/coal/msa.py (docs/alignment.md).

  * Stop handling without any tool: trailing '*' becomes a last column of its own, internal '*' is aligned as X
    and put back; residues and order are kept.
  * With Clustal Omega / MAFFT found (PATH, COAL_CLUSTALO / COAL_MAFFT, or COAL's cache): run_msa and
    align_clusters on real tools (skipped otherwise).
  * align_clusters isolates a failing cluster, aligns single sequences without a tool and 2-sequence clusters with
    the pairwise aligner (pairs_with="nw").
  * install_clustalo: the cached copy is reused; the official-binary route (download faked, nothing fetched)
    extracts the Windows zip layout safely and checks the version; a binary that does not run is refused.
    No downloaded file is run before its SHA-256 was checked: the official binary (no published checksum, none
    pinned) is refused without sha256= or allow_unverified=True, a wrong hash deletes the download, and
    micromamba is one pinned release with pinned hashes. The real download is not run here (network);
    docs/alignment.md records a manual run.
"""
from __future__ import annotations

import hashlib
import io
import os
import shutil
import stat
import sys
import tempfile
import unittest
import urllib.error
import zipfile
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from coal import msa  # noqa: E402

RECS = [("a", "MKVLAAGIVALLLW*"), ("b", "MKVLGIVAL*LLW"), ("c", "MKVAAGIVLLLW")]


def fake_tool(folder: str, name: str, script: str) -> str:
    """A POSIX shell stand-in for a tool."""
    p = os.path.join(folder, name)
    with open(p, "w", newline="\n", encoding="utf-8") as fh:
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
        digest = hashlib.sha256(payload.getvalue()).hexdigest()
        seen = []

        def fake_download(url, dest, quiet):
            seen.append(url)
            self.assertTrue(url.endswith(".zip"))
            with open(dest + ".part", "wb") as fh:
                fh.write(payload.getvalue())
            return digest
        with mock.patch.object(msa, "_download", fake_download), mock.patch.object(msa, "_EXE", ".exe"):
            exe = msa._install_official("Windows", "AMD64", True, sha256=digest)
        self.assertTrue(seen[0].startswith("https://"), "https is tried before http")
        d = os.path.join(self.tmp, "cache", "tools", "clustalo")
        self.assertEqual(exe, os.path.join(d, "clustalo.exe"))
        self.assertTrue(os.path.isfile(os.path.join(d, "libargtable2-0.dll")))
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "evil.txt")))
        self.assertEqual(msa.tool_version(exe), "1.2.2")

    def _serve(self, body: bytes, fail_https: bool = False):
        """A stand-in for urllib.request.urlopen that serves `body` and records the URLs asked for."""
        calls = []

        def urlopen(url, timeout=None):
            calls.append(url)
            if fail_https and url.startswith("https://"):
                raise urllib.error.URLError("connection refused")
            return io.BytesIO(body)
        return calls, mock.patch.object(msa.urllib.request, "urlopen", urlopen)

    def test_official_binary_is_refused_without_a_hash(self):
        """clustal.org publishes no checksums and none is pinned: nothing is downloaded, nothing is run."""
        marker = os.path.join(self.tmp, "ran")
        calls, served = self._serve(f"#!/bin/sh\ntouch {marker}\necho 1.2.4\n".encode())
        with served, mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            with self.assertRaises(msa.MsaError) as cm:
                msa.install_clustalo(method="official", quiet=True)
        text = str(cm.exception)
        self.assertIn("refusing to download and run the official Clustal Omega binary", text)
        self.assertIn("no SHA-256 is known", text)
        self.assertIn("plain HTTP", text)
        self.assertIn("--sha256 HEX", text)
        self.assertIn("--allow-unverified", text)
        self.assertIn("conda install -c bioconda clustalo", text)
        self.assertEqual(calls, [], "refused before any download")
        self.assertFalse(os.path.exists(marker))
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "cache", "tools", "clustalo", "clustalo")))
        for (key, (path, member, pinned)) in msa._OFFICIAL.items():   # no hash was guessed for these files
            self.assertIsNone(pinned, key)

    def test_official_binary_with_the_wrong_hash_is_deleted_not_run(self):
        marker = os.path.join(self.tmp, "ran")
        body = f"#!/bin/sh\ntouch {marker}\necho 1.2.4\n".encode()
        calls, served = self._serve(body)
        with served, mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            with self.assertRaises(msa.MsaError) as cm:
                msa.install_clustalo(method="official", quiet=True, sha256="0" * 64)
        self.assertIn("failed its SHA-256 check and was deleted", str(cm.exception))
        self.assertIn(hashlib.sha256(body).hexdigest(), str(cm.exception))
        d = os.path.join(self.tmp, "cache", "tools", "clustalo")
        self.assertEqual(os.listdir(d), [], "nothing is left in the cache")
        self.assertFalse(os.path.exists(marker), "the file was never run")
        with self.assertRaises(ValueError):
            msa.install_clustalo(method="official", quiet=True, sha256="abc")

    def test_official_binary_with_the_right_hash_or_explicit_permission(self):
        body = b"#!/bin/sh\necho 1.2.4\n"
        digest = hashlib.sha256(body).hexdigest()
        calls, served = self._serve(body, fail_https=True)
        with served, mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            exe = msa.install_clustalo(method="official", quiet=True, sha256=digest.upper())
            self.assertEqual(msa.tool_version(exe), "1.2.4")
            self.assertEqual([u.split(":")[0] for u in calls], ["https", "http"])   # https first, then http
            err = io.StringIO()
            with redirect_stderr(err):
                exe2 = msa.install_clustalo(method="official", quiet=True, force=True, allow_unverified=True)
            self.assertEqual(exe2, exe)
            self.assertIn("was NOT verified", err.getvalue())
            self.assertIn(digest, err.getvalue())   # the hash to pin next time

    def test_official_binary_that_does_not_run_is_refused(self):
        body = b"#!/bin/sh\nexit 1\n"
        calls, served = self._serve(body)
        with served, mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            with self.assertRaises(msa.MsaError) as cm:
                msa.install_clustalo(method="official", quiet=True, sha256=hashlib.sha256(body).hexdigest())
        self.assertIn("is not Clustal Omega 1.2", str(cm.exception))
        self.assertIn("conda install -c bioconda clustalo", str(cm.exception))

    def test_micromamba_is_a_pinned_release_with_a_pinned_hash(self):
        self.assertRegex(msa.MICROMAMBA_VERSION, r"^\d+\.\d+\.\d+-\d+$")
        self.assertNotIn("latest", msa._MICROMAMBA_URL)
        self.assertTrue(msa._MICROMAMBA_URL.startswith("https://github.com/mamba-org/micromamba-releases/"))
        for key, (asset, digest) in msa._MICROMAMBA.items():
            self.assertRegex(digest, r"^[0-9a-f]{64}$", key)
        self.assertEqual(len({d for _, d in msa._MICROMAMBA.values()}), 4)
        marker = os.path.join(self.tmp, "ran")
        calls, served = self._serve(f"#!/bin/sh\ntouch {marker}\n".encode())   # not the pinned file
        with served, mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            with self.assertRaises(msa.MsaError) as cm:
                msa._get_micromamba(True)
        self.assertIn("failed its SHA-256 check", str(cm.exception))
        self.assertEqual(calls, [f"https://github.com/mamba-org/micromamba-releases/releases/download/"
                                 f"{msa.MICROMAMBA_VERSION}/micromamba-linux-64"])
        self.assertEqual(os.listdir(os.path.join(self.tmp, "cache", "tools", "micromamba")), [])
        self.assertFalse(os.path.exists(marker))

    def test_fetch_verified(self):
        body = b"payload"
        digest = hashlib.sha256(body).hexdigest()
        dest = os.path.join(self.tmp, "f")
        calls, served = self._serve(body)
        with served:
            self.assertEqual(msa._fetch_verified(["https://example.org/f"], dest, True, digest), digest)
            with open(dest, "rb") as fh:
                self.assertEqual(fh.read(), body)
            os.remove(dest)
            with self.assertRaises(msa.MsaError):   # a caller's hash never overrides the pinned one
                msa._fetch_verified(["https://example.org/f"], dest, True, digest, sha256="1" * 64)
            with self.assertRaises(msa.MsaError):
                msa._fetch_verified(["https://example.org/f"], dest, True, None)
            self.assertFalse(os.path.exists(dest) or os.path.exists(dest + ".part"))

    def test_auto_order(self):
        tried = []

        def fake_conda_install(tool, spec, conda, quiet, timeout):
            tried.append(("conda", conda, spec))
            raise msa.MsaError("offline")

        def fake_official(sysname, machine, quiet, sha256=None, allow_unverified=False):
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
