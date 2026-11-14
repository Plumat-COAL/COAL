"""Cluster MSAs for the codon step: src/coal/msa.py (docs/alignment.md).

  * Stop handling without any tool: trailing '*' becomes a last column of its own, internal '*' is aligned as X
    and put back; residues and order are kept.
  * With Clustal Omega / MAFFT found (PATH, COAL_CLUSTALO / COAL_MAFFT, or COAL's cache): run_msa and
    align_clusters on real tools (skipped otherwise).
  * align_clusters isolates a failing cluster, aligns single sequences without a tool and 2-sequence clusters with
    the pairwise aligner (pairs_with="nw").
  * install_clustalo: the cached copy is reused; Bioconda (conda / mamba / micromamba) is the only route, the
    former "official" download of the Clustal site is refused whatever is passed and Windows is told what works
    (ClustalSiteGone); a clustalo that does not run is refused.
    No downloaded file is run before its SHA-256 was checked: a download with no known hash is refused without
    sha256= or allow_unverified=True, a wrong hash deletes the download (faked downloads, nothing fetched), and
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
import threading
import time
import unittest
import urllib.error
from contextlib import redirect_stderr, redirect_stdout
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
@unittest.skipIf(os.name == "nt", "POSIX shell stand-ins for the aligner")
class Interrupts(unittest.TestCase):
    """Independent review 2026-11-13, finding 2: the aligner was started first and registered afterwards, so an
    interrupt in between left it running while kill_running() reported nothing left. The window is widened here
    with the module's test hook (_after_start), not with luck."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = self._td.name
        self.started = os.path.join(self.tmp, "started")
        os.mkdir(self.started)
        self.tool = fake_tool(self.tmp, "clustalo", f'touch "{self.started}/$$"\nexec sleep 120')
        self.recs = [("a", "MKV"), ("b", "MKL"), ("c", "MRV")]

    def tearDown(self):
        msa._after_start = None
        msa.kill_running()
        self._td.cleanup()

    def alive(self, pid: int) -> bool:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        # (a killed child that was not yet waited for still answers: read its state)
        try:
            with open(f"/proc/{pid}/stat", encoding="ascii", errors="replace") as fh:
                return fh.read().rsplit(")", 1)[1].split()[0] != "Z"
        except OSError:
            return True

    def test_interrupt_between_start_and_registration(self):
        in_window, interrupted, procs, errors = threading.Event(), threading.Event(), [], []

        def hook(proc):
            procs.append(proc)
            in_window.set()
            interrupted.wait(20)   # the aligner runs, and the interrupt arrives while the start is still here
                                   # (the old code had not registered the process yet)

        def work():
            try:
                msa.run_msa(self.recs, "clustalo", exe=self.tool)
            except msa.MsaError as e:
                errors.append(str(e))
        msa._after_start = hook
        t = threading.Thread(target=work)
        t.start()
        self.assertTrue(in_window.wait(20))
        killed = msa.kill_running()        # the interrupt arrives inside the window
        interrupted.set()
        t.join(20)
        self.assertFalse(t.is_alive(), "run_msa did not return: its aligner was not ended")
        self.assertEqual(killed, 1, "kill_running() did not see the aligner that had just been started")
        self.assertEqual(len(errors), 1, errors)
        self.assertIn("clustalo failed", errors[0])
        self.assertFalse(self.alive(procs[0].pid), "the aligner was left running")

    def test_nothing_is_started_after_the_interrupt(self):
        """A batch that was interrupted starts no further aligner (its worker threads may be about to); the next
        batch and the next run_msa() run normally."""
        begun = msa._interrupt_count()
        self.assertEqual(msa.kill_running(), 0)
        with self.assertRaisesRegex(msa.MsaError, "^interrupted$"):
            msa.run_msa(self.recs, "clustalo", exe=self.tool, _begun=begun)
        self.assertEqual(os.listdir(self.started), [], "an aligner was started after kill_running()")
        ok = fake_tool(self.tmp, "clustalo_ok", 'while [ $# -gt 0 ]; do case "$1" in -i) in="$2"; shift;; '
                                               '-o) out="$2"; shift;; esac; shift; done\ncp "$in" "$out"')
        self.assertEqual(msa.run_msa(self.recs, "clustalo", exe=ok), self.recs)
        res = msa.align_clusters({"x": self.recs}, exe=ok, jobs=1)
        self.assertEqual((res["x"].error, res["x"].rows), (None, self.recs))

    def test_batch_interrupted_while_aligners_start(self):
        """Many clusters, two workers, a delay in every start, and kill_running() from outside while the batch
        runs: every aligner that was started is dead afterwards and every cluster has a result or an error."""
        started_procs = []

        def hook(proc):
            started_procs.append(proc)
            time.sleep(0.05)
        msa._after_start = hook
        out = {}
        t = threading.Thread(target=lambda: out.update(
            msa.align_clusters({k: self.recs for k in range(30)}, exe=self.tool, jobs=2, pairs_with="msa")))
        t.start()
        while len(started_procs) < 2:
            time.sleep(0.01)
        msa.kill_running()
        t.join(30)
        self.assertFalse(t.is_alive())
        self.assertEqual(len(out), 30)
        self.assertLessEqual(len(started_procs), 4, "aligners were started after the interrupt")
        self.assertTrue(all(r.rows is None for r in out.values()))
        self.assertGreaterEqual(sum(r.error == "interrupted" for r in out.values()), 26)
        time.sleep(0.2)
        self.assertEqual([p.pid for p in started_procs if self.alive(p.pid)], [])


class Install(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        p = mock.patch.dict(os.environ, {"COAL_CACHE_DIR": os.path.join(self.tmp, "cache")})
        p.start()
        self.addCleanup(p.stop)

    def _cache_binary_needs_its_record(self, name: str):
        """A binary in the cache without COAL's checksum record is not run and is named once; with the record it
        is found, until the file changes. `name`: the file name the cache lookup uses on this platform."""
        with mock.patch.dict(os.environ, {"PATH": self.tmp}):
            os.environ.pop("COAL_CLUSTALO", None)
            self.assertIsNone(msa.find_clustalo())
            d = os.path.join(self.tmp, "cache", "tools", "clustalo")
            os.makedirs(d)
            exe = fake_tool(d, name, "echo 1.2.4")
            self.assertEqual(msa._cached("clustalo")[0], exe)
            # Finding 5b: a binary in the cache without COAL's checksum record (0.5.0 and earlier fetched Clustal
            # Omega unverified; anything can have been put there) is not run, and that is said once
            err = io.StringIO()
            with redirect_stderr(err):
                self.assertIsNone(msa.find_clustalo())
                self.assertIsNone(msa.find_clustalo())
            self.assertEqual(err.getvalue().count("has no checksum record"), 1, err.getvalue())
            self.assertIn("`coal pipeline install clustalo --force`", err.getvalue())
            self.assertIn(exe, err.getvalue())
            with mock.patch.object(msa, "tool_version", side_effect=AssertionError("the cache binary was run")), \
                    mock.patch.object(msa, "_conda_exe", return_value=None), \
                    mock.patch.object(msa, "_get_micromamba", side_effect=msa.MsaError("no network in the test")), \
                    self.assertRaises(msa.MsaError), redirect_stderr(io.StringIO()):
                msa.install_clustalo(quiet=True)      # it installs afresh instead of trusting what is there
            # with the record written at install time it is found, until the file changes
            msa._write_manifest(exe, "test", "test")
            self.assertEqual(msa.find_clustalo(), exe)
            # (the stand-in is a shell script: its version is given here, so the same lines hold where it cannot run)
            with mock.patch.object(msa, "tool_version", return_value="1.2.4"):
                self.assertEqual(msa.install_clustalo(quiet=True), exe)  # cached: nothing installed
            with open(exe, "a", encoding="utf-8") as fh:
                fh.write("echo changed\n")
            with redirect_stderr(io.StringIO()) as err:
                self.assertIsNone(msa.find_clustalo())
            self.assertIn("is not the file COAL installed", err.getvalue())

    def test_cache_dir_and_find(self):
        self.assertEqual(msa.cache_dir(), os.path.join(self.tmp, "cache"))
        # the cache holds <tool>.exe on Windows, <tool> elsewhere (msa._EXE). The test used to plant `clustalo`
        # everywhere: on Windows the lookup never saw that file, so nothing was said about it (CI, 2026-11-13)
        self._cache_binary_needs_its_record("clustalo" + msa._EXE)

    def test_cache_lookup_with_windows_names(self):
        """The Windows side of the cache lookup, run everywhere: the binary is <tool>.exe, there is no executable
        bit to ask for, and a file named <tool> without .exe is not a cache binary (not used, not reported)."""
        with mock.patch.object(msa, "_EXE", ".exe"), mock.patch.object(msa.os, "name", "nt"):
            self.assertTrue(msa._cached("clustalo")[0].endswith(os.path.join("tools", "clustalo", "clustalo.exe")))
            self.assertTrue(msa._cached("clustalo")[1].endswith(os.path.join("Library", "bin", "clustalo.exe")))
            self._cache_binary_needs_its_record("clustalo.exe")
            exe = msa._cached("clustalo")[0]
            # no executable bit on Windows: a recorded file without one is still the cache binary
            os.remove(exe)
            with open(exe, "w", encoding="utf-8") as fh:
                fh.write("MZ")
            os.chmod(exe, 0o644)
            msa._write_manifest(exe, "test", "test")
            with mock.patch.dict(os.environ, {"PATH": self.tmp}):
                os.environ.pop("COAL_CLUSTALO", None)
                self.assertEqual(msa.find_clustalo(), exe)
                # a file without .exe next to it is not looked at
                os.remove(exe)
                os.remove(exe + msa.MANIFEST_SUFFIX)
                fake_tool(os.path.dirname(exe), "clustalo", "echo 1.2.4")
                with redirect_stderr(io.StringIO()) as err:
                    self.assertIsNone(msa.find_clustalo())
                self.assertEqual(err.getvalue(), "")

    @unittest.skipIf(os.name == "nt", "symbolic links need a privilege on Windows")
    def test_download_is_not_written_through_a_planted_part_file(self):
        """Finding 5a: the download went to the fixed name `<dest>.part`; a link there was written through."""
        body = b"#!/bin/sh\necho 1.2.4\n"
        victim = os.path.join(self.tmp, "victim.txt")
        with open(victim, "w", encoding="utf-8") as fh:
            fh.write("the user's file\n")
        dest = os.path.join(self.tmp, "dl", "tool")
        os.mkdir(os.path.dirname(dest))
        os.symlink(victim, dest + ".part")
        calls, served = self._serve(body)
        with served:
            with self.assertRaises(msa.MsaError):
                msa._fetch_verified(["https://example.org/t"], dest, True, None, sha256="0" * 64)
            with open(victim, encoding="utf-8") as fh:
                self.assertEqual(fh.read(), "the user's file\n")
            msa._fetch_verified(["https://example.org/t"], dest, True, None, sha256=hashlib.sha256(body).hexdigest())
        with open(victim, encoding="utf-8") as fh:
            self.assertEqual(fh.read(), "the user's file\n")
        self.assertEqual(sorted(os.listdir(os.path.dirname(dest))), ["tool", "tool.part"])
        with open(dest, "rb") as fh:
            self.assertEqual(fh.read(), body)

    def _serve(self, body: bytes, fail_https: bool = False):
        """A stand-in for urllib.request.urlopen that serves `body` and records the URLs asked for."""
        calls = []

        def urlopen(url, timeout=None):
            calls.append(url)
            if fail_https and url.startswith("https://"):
                raise urllib.error.URLError("connection refused")
            return io.BytesIO(body)
        return calls, mock.patch.object(msa.urllib.request, "urlopen", urlopen)

    def test_a_download_without_a_known_hash_is_refused_before_it_is_fetched(self):
        calls, served = self._serve(b"anything")
        dest = os.path.join(self.tmp, "tool")
        with served, self.assertRaises(msa.MsaError) as cm:
            msa._fetch_verified(["https://example.org/tool"], dest, True, None, what="the tool")
        text = str(cm.exception)
        self.assertIn("refusing to download and run the tool", text)
        self.assertIn("no SHA-256 is known", text)
        self.assertIn("--sha256 HEX", text)
        self.assertIn("--allow-unverified", text)
        self.assertEqual(calls, [], "refused before any download")
        self.assertFalse(os.path.exists(dest) or os.path.exists(dest + ".part"))

    def test_a_download_with_the_wrong_hash_is_deleted(self):
        body = b"#!/bin/sh\necho 1.2.4\n"
        calls, served = self._serve(body)
        dest = os.path.join(self.tmp, "tool")
        with served, self.assertRaises(msa.MsaError) as cm:
            msa._fetch_verified(["https://example.org/tool"], dest, True, None, sha256="0" * 64, what="the tool")
        self.assertIn("failed its SHA-256 check and was deleted", str(cm.exception))
        self.assertIn(hashlib.sha256(body).hexdigest(), str(cm.exception))
        self.assertEqual(os.listdir(self.tmp), [], "nothing is left")
        with self.assertRaises(ValueError):
            msa._fetch_verified(["https://example.org/tool"], dest, True, None, sha256="abc")

    def test_a_download_with_the_right_hash_or_explicit_permission(self):
        body = b"#!/bin/sh\necho 1.2.4\n"
        digest = hashlib.sha256(body).hexdigest()
        dest = os.path.join(self.tmp, "tool")
        calls, served = self._serve(body, fail_https=True)
        with served:
            got = msa._fetch_verified(["https://example.org/t", "http://example.org/t"], dest, True, None,
                                      sha256=digest.upper())
            self.assertEqual(got, digest)
            self.assertEqual([u.split(":")[0] for u in calls], ["https", "http"])   # the URLs in the order given
            err = io.StringIO()
            with redirect_stderr(err):
                msa._fetch_verified(["https://example.org/t", "http://example.org/t"], dest, True, None,
                                    allow_unverified=True)
            self.assertIn("was NOT verified", err.getvalue())
            self.assertIn(digest, err.getvalue())   # the hash to pin next time

    def test_a_clustalo_that_does_not_run_is_refused(self):
        def fake_conda_install(tool, spec, conda, quiet, timeout):
            return fake_tool(self.tmp, "clustalo", "exit 1")
        with mock.patch.object(msa, "_conda_install", fake_conda_install), \
                mock.patch.object(msa, "_conda_exe", return_value="/opt/conda/bin/mamba"), \
                mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            with self.assertRaises(msa.MsaError) as cm:
                msa.install_clustalo(quiet=True)
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
        # Bioconda only: an installed conda / mamba / micromamba, else the pinned micromamba; nothing after it
        tried = []

        def fake_conda_install(tool, spec, conda, quiet, timeout):
            tried.append((conda, spec))
            raise msa.MsaError("offline")
        with mock.patch.object(msa, "_conda_install", fake_conda_install), \
                mock.patch.object(msa, "_get_micromamba", return_value="/cache/micromamba"), \
                mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            with mock.patch.object(msa, "_conda_exe", return_value="/opt/conda/bin/mamba"):
                with self.assertRaises(msa.MsaError) as cm:
                    msa.install_clustalo(quiet=True)
                self.assertIn("conda: offline", str(cm.exception))
            with mock.patch.object(msa, "_conda_exe", return_value=None):
                with self.assertRaises(msa.MsaError):
                    msa.install_clustalo(quiet=True)
                with self.assertRaises(msa.MsaError):
                    msa.install_clustalo(method="micromamba", quiet=True)
        spec = "bioconda::clustalo=1.2.4=h503566f_10"      # the pinned build for linux-64 (D6b)
        self.assertEqual(tried, [("/opt/conda/bin/mamba", spec), ("/cache/micromamba", spec),
                                 ("/cache/micromamba", spec)])


@unittest.skipIf(os.name == "nt", "shell scripts stand in for conda and the tools")
class PinnedBuilds(unittest.TestCase):
    """D6b (2026-11-14): MAFFT and Clustal Omega are installed as one exact conda package per platform (channel,
    build string and the package's SHA-256), and the manifest, `coal pipeline install` and report.json say which."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        p = mock.patch.dict(os.environ, {"COAL_CACHE_DIR": os.path.join(self.tmp, "cache")})
        p.start()
        self.addCleanup(p.stop)

    def fake_conda(self, sha: str, version_line: str, fail: str = "") -> str:
        """A 'conda' that makes the environment `-p PREFIX` with bin/<tool> and the conda-meta record of the last
        argument's package, with the SHA-256 `sha`; or prints `fail` and exits 1."""
        if fail:
            return fake_tool(self.tmp, "conda", f"echo '{fail}' >&2; exit 1")
        return fake_tool(self.tmp, "conda", f"""
while [ $# -gt 0 ]; do case "$1" in -p) prefix="$2"; shift 2;; *) spec="$1"; shift;; esac; done
pkg="${{spec#*::}}"; name="${{pkg%%=*}}"; rest="${{pkg#*=}}"; ver="${{rest%%=*}}"; build="${{rest#*=}}"
[ "$build" = "$rest" ] && build=hfake_0
mkdir -p "$prefix/bin" "$prefix/conda-meta"
printf '#!/bin/sh\necho "{version_line}"\n' > "$prefix/bin/$name"; chmod +x "$prefix/bin/$name"
printf '{{"name": "%s", "version": "%s", "build": "%s", "channel": "c", "subdir": "linux-64", "url": "https://x/%s-%s-%s.conda", "sha256": "{sha}"}}' "$name" "$ver" "$build" "$name" "$ver" "$build" > "$prefix/conda-meta/$name-$ver-$build.json"
""")

    def test_pins_are_complete(self):
        for tool, channel, version in (("mafft", "conda-forge", msa.MAFFT_VERSION),
                                       ("clustalo", "bioconda", msa.CLUSTALO_VERSION)):
            pins = msa._PACKAGE_PINS[tool]
            self.assertEqual(set(pins), {("Linux", "x86_64"), ("Linux", "aarch64"), ("Darwin", "x86_64"),
                                         ("Darwin", "arm64")})
            for p in pins.values():
                self.assertEqual((p.channel, p.name, p.version), (channel, tool, version))
                self.assertRegex(p.sha256, r"^[0-9a-f]{64}$")
                self.assertRegex(p.build, r"^h[0-9a-f]{7}_\d+$")
                self.assertEqual(p.spec, f"{channel}::{tool}={version}={p.build}")
            self.assertEqual(len({p.sha256 for p in pins.values()}), 4)
        with mock.patch.object(msa, "_platform", return_value=("Linux", "ppc64le")):
            self.assertIsNone(msa.package_pin("mafft"))

    def test_pinned_build_installed_and_recorded(self):
        with mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            pin = msa.package_pin("mafft")
            conda = self.fake_conda(pin.sha256, "v7.526 (2024/Apr/26)")
            with redirect_stderr(io.StringIO()) as err:
                exe = msa._install_pinned("mafft", "MAFFT", conda, False, 60)
            self.assertIn(f"coal: MAFFT build: {pin.spec} (linux-64, {pin.filename}, SHA-256 {pin.sha256}), "
                          "verified", err.getvalue())
            pkg = msa.installed_package("mafft", exe)
            self.assertEqual(pkg, {"spec": pin.spec, "channel": "conda-forge", "subdir": "linux-64",
                                   "filename": pin.filename, "sha256": pin.sha256, "pinned": True})
            self.assertIsNone(msa.pin_note("mafft", exe, "v7.526 (2024/Apr/26)"))
            self.assertEqual(msa._cache_state(exe), "ok")
            # `coal pipeline install mafft` shows the package
            from coal import dnds_pipeline as P
            with mock.patch.object(msa, "install_mafft", return_value=exe), redirect_stderr(io.StringIO()) as err, \
                    redirect_stdout(io.StringIO()) as out:
                self.assertEqual(P.main(["install", "mafft"]), 0)
            self.assertEqual(out.getvalue().strip(), exe)
            self.assertIn(f"coal: MAFFT package: {pin.spec} (linux-64, {pin.filename}, SHA-256 {pin.sha256}, the "
                          "pinned build)", err.getvalue())

    def test_another_package_is_refused(self):
        with mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            pin = msa.package_pin("clustalo")
            conda = self.fake_conda("0" * 64, "1.2.4")
            with self.assertRaises(msa.MsaError) as cm:
                msa._install_pinned("clustalo", "Clustal Omega", conda, True, 60)
            self.assertIn(f"is not the pinned package: expected SHA-256 {pin.sha256}", str(cm.exception))
            self.assertFalse(os.path.exists(os.path.join(msa.cache_dir(), "envs", "clustalo")))

    def test_fallback_to_the_version_pin(self):
        # the channel no longer serves the build: the version pin, with a warning, recorded as not pinned
        specs = []
        real = msa._conda_install

        def install(tool, spec, conda, quiet, timeout):
            specs.append(spec)
            if spec.count("=") == 2:
                raise msa.MsaError(f"installing {spec} with micromamba failed: nothing provides requested {spec}")
            return real(tool, spec, self.fake_conda("ab" * 32, "v7.526 (2024/Apr/26)"), quiet, timeout)
        with mock.patch.object(msa, "_platform", return_value=("Darwin", "arm64")), \
                mock.patch.object(msa, "_conda_install", install), redirect_stderr(io.StringIO()) as err:
            exe = msa._install_pinned("mafft", "MAFFT", "/x/conda", True, 60)
            pkg = msa.installed_package("mafft", exe)
            note = msa.pin_note("mafft", exe, "v7.526 (2024/Apr/26)")
        self.assertEqual(specs, ["conda-forge::mafft=7.526=h99b78c6_0", "mafft=7.526"])
        self.assertIn("warning: the pinned MAFFT build conda-forge::mafft=7.526=h99b78c6_0 is not available",
                      err.getvalue())
        self.assertFalse(pkg["pinned"])
        self.assertEqual(pkg["sha256"], "ab" * 32)
        self.assertIn("not the pinned build mafft-7.526-h99b78c6_0.conda", note)
        # a network failure is not a missing build: no second attempt
        specs.clear()
        with mock.patch.object(msa, "_platform", return_value=("Darwin", "arm64")), \
                mock.patch.object(msa, "_conda_install", side_effect=msa.MsaError("connection reset")):
            with self.assertRaises(msa.MsaError):
                msa._install_pinned("mafft", "MAFFT", "/x/conda", True, 60)

    def test_no_pin_for_the_platform(self):
        calls = []
        with mock.patch.object(msa, "_platform", return_value=("Linux", "ppc64le")), \
                mock.patch.object(msa, "_conda_install", lambda *a: calls.append(a[1]) or "/x/mafft"), \
                redirect_stderr(io.StringIO()) as err:
            self.assertEqual(msa._install_pinned("mafft", "MAFFT", "/x/conda", False, 60), "/x/mafft")
        self.assertEqual(calls, ["mafft=7.526"])
        self.assertIn("COAL pins no MAFFT build for Linux/ppc64le", err.getvalue())


class ClustalSiteGone(unittest.TestCase):
    """The pre-compiled Clustal Omega binaries of the Clustal site no longer exist (2026-11-12: the address serves
    an unrelated site, the files are "not found" pages). COAL holds no address of that site, the "official" route
    is refused whatever is passed, and Windows, which had no other route, is told what works. Nothing here needs
    a shell, so it runs on Windows too."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        p = mock.patch.dict(os.environ, {"COAL_CACHE_DIR": os.path.join(self.tmp, "cache")})
        p.start()
        self.addCleanup(p.stop)
        self.calls = []

        def urlopen(url, timeout=None):
            self.calls.append(url)
            raise AssertionError(f"a download was attempted: {url}")
        p = mock.patch.object(msa.urllib.request, "urlopen", urlopen)
        p.start()
        self.addCleanup(p.stop)

    def test_the_official_route_cannot_be_brought_back(self):
        for system in (("Linux", "x86_64"), ("Darwin", "arm64"), ("Windows", "AMD64")):
            for kw in ({}, {"sha256": "a" * 64}, {"allow_unverified": True},
                       {"sha256": "a" * 64, "allow_unverified": True, "force": True}):
                with self.subTest(system=system, kw=kw), mock.patch.object(msa, "_platform", return_value=system):
                    with self.assertRaises(msa.MsaError) as cm:
                        msa.install_clustalo(method="official", quiet=True, **kw)
                    self.assertIn("binaries are gone", str(cm.exception))
                    self.assertIn("Bioconda", str(cm.exception))
        # the options of that route apply to nothing else: they are an error, not a way round
        for method in ("auto", "conda", "micromamba"):
            for kw in ({"sha256": "a" * 64}, {"allow_unverified": True}):
                with self.subTest(method=method, kw=kw), \
                        mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")), \
                        mock.patch.object(msa, "_conda_install", side_effect=AssertionError("installed")):
                    with self.assertRaises(msa.MsaError) as cm:
                        msa.install_clustalo(method=method, quiet=True, **kw)
                    self.assertIn("do not apply to Clustal Omega", str(cm.exception))
        with self.assertRaises(ValueError):
            msa.install_clustalo(method="official", sha256="abc")     # (a malformed hash is still a ValueError)
        with self.assertRaises(ValueError):
            msa.install_clustalo(method="ftp")
        self.assertEqual(self.calls, [])
        self.assertFalse(hasattr(msa, "_OFFICIAL") or hasattr(msa, "_OFFICIAL_HOST") or hasattr(msa, "_install_official"))

    def test_windows_is_told_what_works(self):
        with mock.patch.object(msa, "_platform", return_value=("Windows", "AMD64")), \
                mock.patch.object(msa, "_conda_install", side_effect=AssertionError("installed")), \
                mock.patch.object(msa, "_get_micromamba", side_effect=AssertionError("fetched")):
            for method in ("auto", "conda", "micromamba"):
                with self.subTest(method=method), self.assertRaises(msa.MsaError) as cm:
                    msa.install_clustalo(method=method, quiet=True)
                text = str(cm.exception)
                self.assertIn("no automatic Clustal Omega install on Windows", text)
                for option in ("WSL", "COAL_CLUSTALO", "--aligner mafft", "COAL_MAFFT", "not tested"):
                    self.assertIn(option, text)
        with mock.patch.object(msa.platform, "system", return_value="Windows"):
            with self.assertRaises(msa.MsaError) as cm:
                msa.install_mafft(quiet=True)
        self.assertIn("no automatic MAFFT install on Windows", str(cm.exception))
        self.assertIn("COAL_MAFFT", str(cm.exception))
        self.assertEqual(self.calls, [])

    def test_no_address_of_that_site_is_left_in_the_repository(self):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        site = "clustal" + ".org"
        hits = []
        for top in ("src", "apps", "cpp", "packaging", "docs", "site", "scripts", "tests", "benchmarks", "examples",
                    ".github", "README.md", "CHANGELOG.md", "NOTICE", "CONTRIBUTING.md"):
            path = os.path.join(root, top)
            files = [path] if os.path.isfile(path) else [os.path.join(d, f) for d, _, fs in os.walk(path) for f in fs]
            for f in files:
                if not f.endswith((".py", ".md", ".html", ".R", ".Rd", ".cpp", ".hpp", ".yml", ".yaml", ".toml", ".txt",
                                   ".json", ".cff", ".js", ".css", ".sh", ".ps1", ".bat")) and os.path.basename(f) != "NOTICE":
                    continue
                with open(f, encoding="utf-8", errors="backslashreplace") as fh:
                    text = fh.read().lower()
                for prefix in ("http://", "https://", "www.", "//"):
                    if prefix + site in text or prefix + "www." + site in text:
                        hits.append(os.path.relpath(f, root))
                        break
        self.assertEqual(hits, [], "a link to the former Clustal site (it now serves an unrelated site)")


if __name__ == "__main__":
    unittest.main()
