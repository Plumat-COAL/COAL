"""Clustal Omega and MAFFT on native Windows (WIN-TOOLS 2026-11-14, WIN-MAFFT 2026-11-15), and the small texts that
came with them.

Clustal Omega: no automatic install on Windows (neither Bioconda nor conda-forge has a win-64 package, and COAL
adds a download only with a SHA-256 computed from the file itself): the installer downloads and runs nothing there,
and its messages (and the pipeline's) name what works. MAFFT: on x86-64 Windows `coal pipeline install mafft`
installs MAFFT's own Windows package, one zip pinned by size and SHA-256 (coal.msa.MAFFT_WINDOWS), checked before
anything is unpacked, unpacked member by member under checked names, recorded like every other install; on any
other Windows machine nothing is downloaded. Batch files (mafft.bat) are started with a command line COAL builds
for cmd.exe (coal.msa.batch_command_line). Also: the Windows-path hint of `coal`'s missing-file message (U17), the
app's "needs coal-cpp" text in a development version, and the D9 wording of the help.

Everything here is mocked (a zip made by the test stands in for MAFFT's, with its own size and hash pinned in its
place), so it runs on every OS; the tests marked for Windows start real batch files. The real download is
tests/test_windows_mafft_real.py (COAL_TEST_REAL_MAFFT=1)."""
from __future__ import annotations

import argparse
import dataclasses
import errno
import hashlib
import io
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time
import unittest
import zipfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from coal import msa, seqio  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402

VERSION_LINE = f"v{msa.MAFFT_VERSION} (2024/Apr/26)"


class WindowsCase(unittest.TestCase):
    machine = "AMD64"

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.calls = []
        for target, name in ((msa.urllib.request, "urlopen"), (msa, "_download"), (msa, "_conda_install"),
                             (msa, "_get_micromamba"), (msa.subprocess, "run"), (msa.subprocess, "Popen")):
            p = mock.patch.object(target, name, self._refuse(name))
            p.start()
            self.addCleanup(p.stop)
        for p in (mock.patch.dict(os.environ, {"COAL_CACHE_DIR": os.path.join(self.tmp, "cache")}),
                  mock.patch.object(msa, "_platform", return_value=("Windows", self.machine)),
                  mock.patch.object(msa.platform, "system", return_value="Windows")):
            p.start()
            self.addCleanup(p.stop)

    def _refuse(self, name):
        def f(*a, **k):
            self.calls.append((name, a))
            raise AssertionError(f"{name} was called on Windows: {a}")
        return f


class TestNothingIsDownloadedOnWindows(WindowsCase):
    """Clustal Omega on Windows, and MAFFT on a Windows machine that is not x86-64: nothing is downloaded or run."""
    machine = "ARM64"

    def test_clustalo_every_method(self):
        for machine in ("AMD64", "ARM64"):
            for method in ("auto", "conda", "micromamba"):
                for force in (False, True):
                    with self.subTest(machine=machine, method=method, force=force), \
                            mock.patch.object(msa, "_platform", return_value=("Windows", machine)), \
                            self.assertRaises(msa.MsaError) as cm:
                        msa.install_clustalo(method=method, force=force, quiet=True)
                    self.assertEqual(str(cm.exception), msa.NO_WINDOWS_CLUSTALO)
        self.assertEqual(self.calls, [])

    def test_mafft_every_method_on_arm64(self):
        for method in ("auto", "official", "conda", "micromamba"):
            for force in (False, True):
                with self.subTest(method=method, force=force), self.assertRaises(msa.MsaError) as cm:
                    msa.install_mafft(method=method, force=force, quiet=True)
                if method in ("auto", "official"):
                    self.assertEqual(str(cm.exception), msa.NO_WINDOWS_MAFFT)
                else:
                    self.assertIn("no Windows package of MAFFT", str(cm.exception))
        self.assertEqual(self.calls, [])

    def test_mafft_conda_methods_on_amd64(self):
        with mock.patch.object(msa, "_platform", return_value=("Windows", "AMD64")):
            for method in ("conda", "micromamba"):
                with self.subTest(method=method), self.assertRaises(msa.MsaError) as cm:
                    msa.install_mafft(method=method, quiet=True)
                self.assertIn("(method='auto' or 'official')", str(cm.exception))
            with self.assertRaises(ValueError):
                msa.install_mafft(method="ftp")
        self.assertEqual(self.calls, [])

    def test_the_install_command(self):
        for tool, text, machine in (("clustalo", msa.NO_WINDOWS_CLUSTALO, "AMD64"),
                                    ("mafft", msa.NO_WINDOWS_MAFFT, "ARM64")):
            err = io.StringIO()
            with self.subTest(tool=tool), redirect_stderr(err), \
                    mock.patch.object(msa, "_platform", return_value=("Windows", machine)):
                rc = P.main(["install", tool, "-q"])
            self.assertEqual(rc, 1)
            self.assertIn(text, err.getvalue())
        self.assertEqual(self.calls, [])

    def test_the_one_windows_download_is_pinned(self):
        # MAFFT's Windows package, by size and SHA-256 (computed by the organiser from two identical downloads of
        # MAFFT's page's file, 2026-11-15); nothing else for Windows: no micromamba, no conda build
        w = msa.MAFFT_WINDOWS
        self.assertEqual(w.url, "https://mafft.cbrc.jp/alignment/software/mafft-7.526-win64-signed.zip")
        self.assertEqual(w.filename, f"mafft-{msa.MAFFT_VERSION}-win64-signed.zip")
        self.assertEqual(w.size, 27_476_513)
        self.assertEqual(w.sha256, "7c3a508c8512af04007d68019c6fbd4ff00d891244b58f87ec9aa0e53fbae8a7")
        self.assertEqual((w.top, w.entry), ("mafft-win", "mafft.bat"))
        self.assertNotIn(("Windows", "AMD64"), msa._MICROMAMBA)
        for machine in ("AMD64", "ARM64"):
            with mock.patch.object(msa, "_platform", return_value=("Windows", machine)):
                for tool in ("clustalo", "mafft"):
                    self.assertIsNone(msa.package_pin(tool))
        src = (ROOT / "src" / "coal" / "msa.py").read_text(encoding="utf-8")
        # one Windows archive address in the download tables: MAFFT's, built from the version pin
        self.assertEqual(re.findall(r"https?://\S*\.zip", src),
                         ["https://mafft.cbrc.jp/alignment/software/mafft-{MAFFT_VERSION}-win64-signed.zip"])


class TestWhatTheMessagesSay(unittest.TestCase):
    def test_clustal_omega(self):
        t = msa.NO_WINDOWS_CLUSTALO
        for must in ("no automatic Clustal Omega install on Windows", "Bioconda and conda-forge", "WSL",
                     "coal pipeline install clustalo", "MSYS2", "https://www.msys2.org",
                     "pacman -S mingw-w64-ucrt-x86_64-clustal-omega", "COAL_CLUSTALO",
                     "C:\\msys64\\ucrt64\\bin\\clustalo.exe", "PATH", "--aligner mafft", "COAL_MAFFT", "mafft.bat",
                     "not tested", "`coal pipeline install mafft` installs MAFFT's own Windows package"):
            self.assertIn(must, t)
        self.assertIn(msa.CLUSTALO_VERSION, t)       # MSYS2 packages the version COAL pins elsewhere
        self.assertNotIn("clustal" + ".org", t)      # (no address of the former Clustal site anywhere)

    def test_mafft(self):
        t = msa.NO_WINDOWS_MAFFT
        for must in ("no automatic MAFFT install for this Windows machine", "Bioconda and conda-forge",
                     "x86-64 (AMD64) Windows only",
                     "https://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html", "LF line ends",
                     "COAL_MAFFT", "mafft-win\\mafft.bat", "not checked by COAL", "WSL",
                     "coal pipeline install mafft"):
            self.assertIn(must, t)
        self.assertEqual(msa.MAFFT_VERSION, "7.526")   # the zip pinned is that version's (MAFFT's page)

    def test_exported(self):
        for name in ("NO_WINDOWS_MAFFT", "NO_WINDOWS_CLUSTALO", "MAFFT_WINDOWS", "batch_command_line",
                     "is_batch_file"):
            self.assertIn(name, msa.__all__)

    def test_install_help_says_what_windows_gets(self):
        out = io.StringIO()
        with mock.patch("sys.stdout", out), self.assertRaises(SystemExit):
            P.main(["install", "--help"])
        text = " ".join(out.getvalue().split())
        self.assertIn("Clustal Omega and MAFFT: Linux and macOS (conda packages)", text)
        self.assertIn("MAFFT on x86-64 Windows: MAFFT's own Windows package (mafft-7.526-win64-signed.zip from "
                      "mafft.cbrc.jp, pinned by size and SHA-256)", text)
        self.assertIn("Clustal Omega on Windows: the command says what works there instead", text)
        self.assertIn("(conda packages, and MAFFT's pinned Windows package)", text)
        self.assertNotIn("Bioconda only", text)       # MAFFT comes from conda-forge


def make_zip(members) -> bytes:
    """A zip of (name, data) members; data None: a folder entry; ("link", target): a symbolic link."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, data in members:
            if data is None:
                zf.writestr(zipfile.ZipInfo(name), b"")
            elif isinstance(data, tuple):
                info = zipfile.ZipInfo(name)
                info.create_system = 3
                info.external_attr = (stat.S_IFLNK | 0o777) << 16
                zf.writestr(info, data[1])
            else:
                zf.writestr(name, data)
    return buf.getvalue()


GOOD = [("mafft-win/", None), ("mafft-win/mafft.bat", b"@echo off\n\"%~dp0usr\\bin\\bash\" /usr/bin/mafft %*\n"),
        ("mafft-win/usr/", None), ("mafft-win/usr/bin/mafft", b"#!/bin/bash\necho MAFFT\n"),
        ("mafft-win/usr/bin/msys-2.0.dll", b"MZ dll"), ("mafft-win/usr/lib/mafft/disttbfast.exe", b"MZ exe"),
        ("mafft-win/license.txt", b"BSD; the package: GPL\n")]


ORIG_DOWNLOAD = msa._download


class Served:
    """urlopen stand-in: serves `data` for the pinned address and counts the requests."""
    def __init__(self):
        self.data, self.urls = b"", []

    def __call__(self, url, timeout=None):
        self.urls.append(url)
        if url != msa.MAFFT_WINDOWS.url:
            raise AssertionError(f"unexpected download: {url}")
        return io.BytesIO(self.data)


class TestMafftWindowsInstall(WindowsCase):
    """`coal pipeline install mafft` on x86-64 Windows with the download mocked: the file served has to be the
    pinned one (size, SHA-256) before anything is unpacked."""

    def setUp(self):
        super().setUp()
        self.real_pin = msa.MAFFT_WINDOWS
        self.served = Served()
        self.versions = []
        self.unpacked = []
        real_unpack = msa._unpack_windows_package

        def unpack(*a, **k):
            self.unpacked.append(a[0])
            return real_unpack(*a, **k)

        def version(exe, timeout=60):
            self.versions.append(exe)
            return self.version_line
        self.version_line = VERSION_LINE
        for p in (mock.patch.object(msa.urllib.request, "urlopen", self.served),
                  mock.patch.object(msa, "_download", ORIG_DOWNLOAD),   # (WindowsCase refuses it)
                  mock.patch.object(msa, "tool_version", version),
                  mock.patch.object(msa, "_unpack_windows_package", unpack),
                  mock.patch.object(msa.shutil, "which", return_value=None),
                  mock.patch.dict(os.environ, {"COAL_MAFFT": ""})):
            p.start()
            self.addCleanup(p.stop)
        os.environ.pop("COAL_MAFFT")
        self.serve(GOOD)

    def serve(self, members, pin=True, size=None, sha=None):
        """Serve a zip of `members`; pin=True: pin its own size and hash in place of MAFFT's."""
        self.served.data = make_zip(members)
        if pin:
            p = mock.patch.object(msa, "MAFFT_WINDOWS", dataclasses.replace(
                self.real_pin, size=len(self.served.data) if size is None else size,
                sha256=sha or hashlib.sha256(self.served.data).hexdigest()))
            p.start()
            self.addCleanup(p.stop)

    @property
    def folder(self) -> Path:
        return Path(self.tmp, "cache", "tools", "mafft")

    def entry(self) -> str:
        return str(self.folder / "mafft-win" / "mafft.bat")

    def install(self, **kw) -> str:
        kw.setdefault("quiet", True)
        return msa.install_mafft(**kw)

    def assert_nothing_installed(self):
        left = sorted(p.name for p in self.folder.iterdir()) if self.folder.is_dir() else []
        self.assertEqual(left, [])            # no package, no archive, no half-unpacked folder
        self.assertIsNone(msa.find_mafft())

    def test_good_hash_installs(self):
        err = io.StringIO()
        with redirect_stderr(err):
            path = self.install(quiet=False)
        self.assertEqual(path, self.entry())
        self.assertEqual(self.served.urls, [msa.MAFFT_WINDOWS.url])
        for name, data in GOOD:
            f = self.folder / name
            self.assertTrue(f.is_dir() if data is None else f.read_bytes() == data, name)
        self.assertEqual(sorted(p.name for p in self.folder.iterdir()), ["mafft-win"])   # the archive is gone
        rec = json.loads(Path(path + msa.MANIFEST_SUFFIX).read_text(encoding="utf-8"))
        self.assertEqual(sorted(rec["files"]), ["license.txt", "mafft.bat", "usr/bin/mafft", "usr/bin/msys-2.0.dll",
                                                "usr/lib/mafft/disttbfast.exe"])
        self.assertEqual(rec["source"], msa.MAFFT_WINDOWS.url)
        self.assertEqual(rec["verified"], "pinned size and SHA-256 of the archive")
        pkg = msa.installed_package("mafft", path)
        self.assertEqual(pkg, rec["package"])
        self.assertEqual(pkg, {"spec": msa.MAFFT_WINDOWS.url, "channel": "mafft.cbrc.jp", "subdir": "win-64",
                               "filename": msa.MAFFT_WINDOWS.filename, "sha256": msa.MAFFT_WINDOWS.sha256,
                               "size": msa.MAFFT_WINDOWS.size, "pinned": True})
        self.assertEqual(msa._cache_state(path), "ok")
        self.assertEqual(msa.find_mafft(), path)
        self.assertIsNone(msa.pin_note("mafft", path, VERSION_LINE))
        self.assertIn(f"coal: downloading {msa.MAFFT_WINDOWS.url}", err.getvalue())
        self.assertIn(f"coal: MAFFT {msa.MAFFT_VERSION} installed: {path}", err.getvalue())
        self.assertEqual(self.versions, [path])     # the version is read from the batch file once installed

    def test_bad_hash_refuses_before_extraction(self):
        self.serve(GOOD, sha="0" * 64)
        with self.assertRaises(msa.MsaError) as cm:
            self.install()
        self.assertIn("failed its SHA-256 check and was deleted", str(cm.exception))
        self.assertIn("expected " + "0" * 64, str(cm.exception))
        self.assertEqual(self.unpacked, [])
        self.assertEqual(self.versions, [])
        self.assert_nothing_installed()

    def test_the_real_pin_refuses_another_file(self):
        self.serve(GOOD, pin=False)       # MAFFT's pin: this file is neither its size nor its hash
        with mock.patch.object(msa, "MAFFT_WINDOWS", self.real_pin), self.assertRaises(msa.MsaError) as cm:
            self.install()
        self.assertIn(f"not the {self.real_pin.size:,} bytes COAL pins for it", str(cm.exception))
        self.assertEqual(self.unpacked, [])
        self.assert_nothing_installed()

    def test_wrong_size(self):
        n = len(make_zip(GOOD))
        for size, said in ((n + 1, f"is {n:,} bytes, not the {n + 1:,} bytes"),
                           (n - 1, f"larger than the {n - 1:,} bytes COAL pins for it; the download was stopped")):
            with self.subTest(size=size):
                self.serve(GOOD, size=size)      # the right hash, the wrong size: refused all the same
                with self.assertRaises(msa.MsaError) as cm:
                    self.install()
                self.assertIn(said, str(cm.exception))
                self.assertEqual(self.unpacked, [])
                self.assert_nothing_installed()

    def test_malicious_member_names_rejected(self):
        outside = Path(self.tmp, "evil")
        bad = ["/abs/evil", "../evil", "mafft-win/../../evil", "mafft-win/./x", "mafft-win/a\\b", "C:/evil",
               "mafft-win/c:evil", "mafft-win/x:stream", "mafft-win/CON", "mafft-win/nul.txt", "mafft-win/.hidden",
               "mafft-win/trailing.", "mafft-win/sp ace", "mafft-win//double", "other/evil", "evil",
               "mafft-win/caf\u00e9"]
        for name in bad:
            with self.subTest(name=name):
                self.serve(GOOD + [(name, b"x")])
                with self.assertRaises(msa.MsaError) as cm:
                    self.install()
                self.assertIn("unsafe name", str(cm.exception))
                self.assertIn(repr(name), str(cm.exception))
                self.assertFalse(outside.exists())
                self.assert_nothing_installed()
        for members, said in (([("mafft-win/usr/bin/sh", ("link", b"/bin/sh"))], "unsafe name or type"),
                              ([("mafft-win/MAFFT.BAT", b"x")], "twice (letter case aside)"),
                              ([("mafft-win/license.txt/x", b"x")], "a file that is also a folder")):
            with self.subTest(members=members):
                self.serve(GOOD + members)
                with self.assertRaises(msa.MsaError) as cm:
                    self.install()
                self.assertIn(said, str(cm.exception))
                self.assert_nothing_installed()
        self.serve([m for m in GOOD if m[0] != "mafft-win/mafft.bat"])
        with self.assertRaises(msa.MsaError) as cm:
            self.install()
        self.assertIn("holds no mafft-win/mafft.bat", str(cm.exception))
        self.assert_nothing_installed()

    def test_rerun_is_a_no_op(self):
        path = self.install()
        mtime = os.stat(path).st_mtime_ns
        err = io.StringIO()
        with redirect_stderr(err):
            self.assertEqual(self.install(quiet=False), path)
        self.assertEqual(len(self.served.urls), 1)              # not downloaded again
        self.assertEqual(os.stat(path).st_mtime_ns, mtime)
        self.assertNotIn("downloading", err.getvalue())
        # an earlier MAFFT in the cache: kept, and named (as on Linux)
        self.version_line = "v7.525 (2024/Mar/13)"
        with redirect_stderr(err):
            self.assertEqual(self.install(quiet=False), path)
        self.assertIn("is not the pinned 7.526", err.getvalue())
        self.assertEqual(len(self.served.urls), 1)

    def test_force_reinstalls(self):
        path = self.install()
        Path(self.folder, "mafft-win", "stale.txt").write_text("from before", encoding="utf-8")
        self.assertEqual(self.install(force=True), path)
        self.assertEqual(len(self.served.urls), 2)
        self.assertFalse(Path(self.folder, "mafft-win", "stale.txt").exists())    # a clean folder, not a merge
        self.assertEqual(sorted(p.name for p in self.folder.iterdir()), ["mafft-win"])
        self.assertEqual(msa._cache_state(path), "ok")

    def test_a_changed_file_is_not_used_and_is_reinstalled(self):
        path = self.install()
        Path(self.folder, "mafft-win", "usr", "lib", "mafft", "disttbfast.exe").write_bytes(b"MZ changed")
        msa._trust_seen.clear()
        with redirect_stderr(io.StringIO()) as err:
            self.assertIsNone(msa.find_mafft())
        self.assertIn("is not the file COAL installed", err.getvalue())
        self.assertEqual(self.install(), path)                  # not trusted: installed again, without --force
        self.assertEqual(len(self.served.urls), 2)
        self.assertEqual(msa.find_mafft(), path)

    def test_a_version_that_does_not_run_fails(self):
        self.version_line = ""
        with self.assertRaises(msa.MsaError) as cm:
            self.install()
        self.assertIn(f"does not run or is not MAFFT {msa.MAFFT_VERSION} (no output)", str(cm.exception))

    def test_the_install_command(self):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            self.assertEqual(P.main(["install", "mafft"]), 0)
        w = msa.MAFFT_WINDOWS
        self.assertEqual(out.getvalue().strip(), self.entry())
        self.assertIn(f"coal: MAFFT package: {w.url} (win-64, {w.filename}, SHA-256 {w.sha256}, the pinned build)",
                      err.getvalue())



class TestPipelineOnWindows(WindowsCase):
    """`coal pipeline clusters` with no aligner found on Windows: Clustal Omega (and MAFFT on a machine that is not
    x86-64) get the installer's text; MAFFT on x86-64 is installed with --install-tools, or the message says how."""

    def run_tool(self, tool, install_tools, machine="AMD64"):
        r = P._Run.__new__(P._Run)
        r.tools, r.quiet = {}, True
        r.a = argparse.Namespace(clustalo=None, mafft=None, install_tools=install_tools)
        with mock.patch.object(msa, "find_clustalo", return_value=None), \
                mock.patch.object(msa, "find_mafft", return_value=None), \
                mock.patch.object(msa, "_platform", return_value=("Windows", machine)), \
                mock.patch.object(P.platform, "system", return_value="Windows"), \
                mock.patch.dict(os.environ, {"COAL_CLUSTALO": "", "COAL_MAFFT": ""}), \
                self.assertRaises(P.PipelineError) as cm:
            r.msa_tool(tool)
        return str(cm.exception)

    def test_not_found(self):
        for install_tools in (False, True):
            with self.subTest(install_tools=install_tools):
                self.assertIn(msa.NO_WINDOWS_CLUSTALO, self.run_tool("clustalo", install_tools))
                self.assertIn(msa.NO_WINDOWS_MAFFT, self.run_tool("mafft", install_tools, "ARM64"))
        self.assertIn("Run again with --install-tools, or run `coal pipeline install mafft`",
                      self.run_tool("mafft", False))
        self.assertEqual(self.calls, [])

    def test_install_tools_installs_mafft(self):
        with mock.patch.object(msa, "install_mafft", side_effect=msa.MsaError("download failed (test)")) as inst:
            self.assertEqual(self.run_tool("mafft", True), "download failed (test)")
        inst.assert_called_once_with(quiet=True)


class TestBatchFiles(unittest.TestCase):
    """mafft.bat: Windows runs a batch file through cmd.exe, which reads its command line by its own rules."""

    def test_quoting(self):
        line = msa.batch_command_line(["C:\\Users\\A B\\mafft-win\\mafft.bat", "--auto", "--thread", "2",
                                       "C:\\T\\R&D (x)\\in.fasta", "a^b", "a,b;c=d", "caf\u00e9"])
        self.assertEqual(line, '"C:\\Users\\A B\\mafft-win\\mafft.bat" --auto --thread 2 "C:\\T\\R&D (x)\\in.fasta" '
                               '"a^b" "a,b;c=d" "caf\u00e9"')
        self.assertEqual(msa.batch_command_line(["C:\\m\\mafft.bat", "C:/t/in.fasta", "-x+1@y"]),
                         "C:\\m\\mafft.bat C:/t/in.fasta -x+1@y")

    def test_refused(self):
        for arg, said in (("C:\\T\\%USERNAME%\\in.fasta", "'%'"), ("C:\\T\\Hi!\\in.fasta", "'!'"),
                          ('C:\\T\\a"b\\in.fasta', "'\"'"), ("C:\\T\\a\nb", "'\\n'"),
                          ("C:\\T\\dir\\", "a backslash at its end"), ("", "an empty argument")):
            with self.subTest(arg=arg), self.assertRaises(msa.MsaError) as cm:
                msa.batch_command_line(["C:\\m\\mafft.bat", "--auto", arg])
            self.assertIn(said, str(cm.exception))
            self.assertIn("cmd.exe", str(cm.exception))
            self.assertIn("TMP / TEMP", str(cm.exception))

    def test_launch(self):
        self.assertTrue(msa.is_batch_file("C:\\x\\MAFFT.BAT") and msa.is_batch_file("a.cmd"))
        self.assertFalse(msa.is_batch_file("/usr/bin/mafft") or msa.is_batch_file("mafft.exe"))
        self.assertEqual(msa._launch(["/usr/bin/mafft", "--auto"]), (["/usr/bin/mafft", "--auto"], {}))
        args, extra = msa._launch(["C:\\m\\mafft.bat", "--auto", "C:\\a b\\in.fasta"])
        self.assertEqual(extra, {"stdin": subprocess.DEVNULL})
        if os.name == "nt":
            self.assertEqual(args, 'C:\\m\\mafft.bat --auto "C:\\a b\\in.fasta"')

    def test_run_msa_refuses_before_starting(self):
        with mock.patch.object(msa.subprocess, "Popen", side_effect=AssertionError("started")), \
                mock.patch.object(msa.tempfile, "tempdir", None), \
                mock.patch.dict(os.environ, {"TMPDIR": "", "TMP": "", "TEMP": ""}):
            d = tempfile.mkdtemp(suffix="%x%")
            self.addCleanup(shutil.rmtree, d, True)
            with mock.patch.object(msa.tempfile, "tempdir", d), self.assertRaises(msa.MsaError) as cm:
                msa.run_msa([("a", "MKV"), ("b", "MKI"), ("c", "MRV")], "mafft", exe="C:\\m\\mafft.bat")
        self.assertIn("cannot pass", str(cm.exception))
        self.assertEqual(msa.tool_version("C:\\m\\%x%\\mafft.bat"), "")

    def test_stopping_a_batch_file_ends_its_tree(self):
        proc = mock.Mock(pid=4242, _coal_batch=True)
        proc.poll.return_value = None
        with mock.patch.object(msa.os, "name", "nt"), mock.patch.object(msa.subprocess, "run") as run:
            msa._end(proc)
        self.assertEqual(run.call_args[0][0], ["taskkill", "/F", "/T", "/PID", "4242"])
        proc.kill.assert_called_once_with()
        plain = mock.Mock(pid=1, _coal_batch=False)
        with mock.patch.object(msa.subprocess, "run") as run:
            msa._end(plain)
        run.assert_not_called()
        plain.kill.assert_called_once_with()


ECHO = """import json, sys
out = sys.argv[1]
with open(out, "w", encoding="utf-8") as fh:
    json.dump(sys.argv[2:], fh)
"""
FAKE_MAFFT = """import sys, time
args = sys.argv[1:]
if args == ["--version"]:
    sys.stderr.write("v7.526 (2024/Apr/26)\\n")
    sys.exit(0)
if "--sleep" in args:
    time.sleep(120)
with open(args[-1], encoding="utf-8") as fh:   # the input FASTA: equal lengths here, so it is its own alignment
    sys.stdout.write(fh.read())
"""


@unittest.skipUnless(os.name == "nt", "a real batch file needs Windows (cmd.exe)")
class TestRealBatchFile(unittest.TestCase):
    """On Windows: a batch file started as COAL starts mafft.bat gets every argument as it was, in a folder whose
    path has spaces and cmd.exe's special characters, and stopping it ends what it started."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="coal bat & (x) ^ ")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def bat(self, script: str, name: str) -> str:
        py = os.path.join(self.tmp, name + ".py")
        Path(py).write_text(script, encoding="utf-8")
        bat = os.path.join(self.tmp, name + ".bat")
        Path(bat).write_text(f'@echo off\r\n"{sys.executable}" "{py}" %*\r\n', encoding="utf-8")
        return bat

    def test_arguments_arrive_as_they_were(self):
        bat = self.bat(ECHO, "echo")
        out = os.path.join(self.tmp, "args.json")
        args = [out, "--auto", "--thread", "2", os.path.join(self.tmp, "in.fasta"), "a&b", "x(1)", "p^q", "s;c,e=1",
                "C:/fwd/slash", "caf\u00e9"]
        cmd, extra = msa._launch([bat, *args])
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120, **extra)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(json.loads(Path(out).read_text(encoding="utf-8")), args[1:])

    def test_run_msa_and_version_through_a_batch_file(self):
        bat = self.bat(FAKE_MAFFT, "mafft")
        self.assertEqual(msa.tool_version(bat), "v7.526 (2024/Apr/26)")
        with mock.patch.object(msa.tempfile, "tempdir", self.tmp):     # the input FASTA's folder: spaces, & ( ) ^
            rows = msa.run_msa([("a", "MKV"), ("b", "MKI"), ("c", "MRV")], "mafft", exe=bat)
        self.assertEqual(rows, [("a", "MKV"), ("b", "MKI"), ("c", "MRV")])

    def test_a_timeout_ends_the_whole_tree(self):
        bat = self.bat(FAKE_MAFFT, "mafft")
        t = time.monotonic()
        with self.assertRaises(subprocess.TimeoutExpired):
            msa.run_msa([("a", "MKV"), ("b", "MKI"), ("c", "MRV")], "mafft", exe=bat, extra_args=["--sleep"],
                        timeout=3)
        self.assertLess(time.monotonic() - t, 60)     # cmd.exe alone killed: Python held the pipes for 120 s


class TestWindowsPathHint(unittest.TestCase):
    """U17 for `coal PEP NUC` (the pipeline had it already): a Windows path missing on Linux / macOS names its WSL
    form. Python engine here; coal-cpp says the same (cpp/src/seqio.cpp open_error), compared byte for byte by
    tests/test_engine_ux.py (cases windows_path, windows_path_slash)."""

    def enoent(self, path):
        return OSError(errno.ENOENT, "No such file or directory", path)

    @unittest.skipIf(os.name == "nt", "on Windows a Windows path is just a path")
    def test_hint(self):
        for path, wsl in (("C:\\Users\\me\\p.fa", "/mnt/c/Users/me/p.fa"), ("D:/data/x y/n.fa", "/mnt/d/data/x y/n.fa")):
            self.assertEqual(seqio.open_error_message(path, self.enoent(path)),
                             f"cannot open '{path}': no such file. This is a Windows path; under WSL (Windows "
                             f"Subsystem for Linux) the same file is {wsl}")
        # the same words as the pipeline's
        self.assertEqual("cannot open 'C:\\a.fa': no such file" + seqio.windows_path_hint("C:\\a.fa"),
                         P._cannot_open("C:\\a.fa"))

    def test_no_hint(self):
        for path in ("p.fa", "/data/p.fa", "C:p.fa", "CD:\\p.fa", "\\\\server\\share\\p.fa"):
            self.assertEqual(seqio.open_error_message(path, self.enoent(path)), f"cannot open '{path}': no such file")
        self.assertEqual(seqio.open_error_message("C:\\x.fa", OSError(errno.EACCES, "denied", "C:\\x.fa")),
                         "cannot open 'C:\\x.fa': permission denied")
        with mock.patch.object(seqio.os, "name", "nt"):
            self.assertEqual(seqio.windows_path_hint("C:\\x.fa"), "")

    @unittest.skipIf(os.name == "nt", "on Windows a Windows path is just a path")
    def test_coal_command(self):
        import subprocess
        prot = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, prot, True)
        (prot / "p.fa").write_text(">a\nMK\n", encoding="utf-8")
        env = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "COAL_ENGINE": "python"}
        r = subprocess.run([sys.executable, "-m", "coal", "C:\\Users\\me\\n.fa", str(prot / "p.fa")],
                           capture_output=True, text=True, env=env, timeout=120, stdin=subprocess.DEVNULL)
        self.assertEqual(r.returncode, 1, r.stderr)
        self.assertIn("coal: error: cannot open 'C:\\Users\\me\\n.fa': no such file. This is a Windows path; under "
                      "WSL (Windows Subsystem for Linux) the same file is /mnt/c/Users/me/n.fa", r.stderr)


class TestHelpD9Wording(unittest.TestCase):
    """--comeron-min-sites: the D9 refinement (a class that would leave the distance at 0 stays NA) in both
    engines' help (the copy in _cpp_help.py is compared with coal-cpp's own by tests/test_cli_commands.py)."""

    def test_both_copies(self):
        new = "dS NA unless it is then 0 (its logarithm is 0, status warning; default 3; 0: any can)"
        cpp = (ROOT / "cpp" / "src" / "kaks_cli.cpp").read_text(encoding="utf-8")
        from coal._cpp_help import HELP
        for name in ("kaks", "dnds"):
            self.assertIn(new, HELP[name], name)
        self.assertEqual(cpp.count(new), 2)
        self.assertNotIn("any class can)", cpp)


class TestAppNeedsCpp(unittest.TestCase):
    """The app's "needs coal-cpp" text is development-version aware, as coal.cli.needs_cpp_text is (U01)."""

    def setUp(self):
        app_src = ROOT / "apps" / "coal-app" / "src"
        if not app_src.is_dir():
            self.skipTest("no apps/coal-app")
        sys.path.insert(0, str(app_src))
        self.addCleanup(sys.path.remove, str(app_src))
        try:
            from coal_app import dnds
        except ImportError as e:   # the app's own dependencies
            self.skipTest(f"coal_app does not import here: {e}")
        self.dnds = dnds

    def test_texts(self):
        dev = self.dnds.needs_cpp_text("0.6.0.dev0")
        self.assertIn("development version 0.6.0.dev0", dev)
        self.assertIn("COAL_CPP", dev)
        self.assertNotIn("releases/latest", dev)
        rel = self.dnds.needs_cpp_text("0.6.0")
        self.assertIn("https://github.com/plumat-coal/COAL/releases/latest", rel)
        self.assertNotIn("development", rel)
        from coal import __version__
        self.assertEqual(self.dnds.NEEDS_CPP, self.dnds.needs_cpp_text(__version__))
        self.assertEqual(".dev" in __version__, "releases/latest" not in self.dnds.NEEDS_CPP)


class TestDocs(unittest.TestCase):
    def test_install_md_has_the_windows_routes(self):
        text = (ROOT / "docs" / "install.md").read_text(encoding="utf-8")
        self.assertIn("## Clustal Omega and MAFFT on Windows", text)
        self.assertIn("(#clustal-omega-and-mafft-on-windows)", text)
        for must in ("pacman -S mingw-w64-ucrt-x86_64-clustal-omega", f"mafft-{msa.MAFFT_VERSION}-win64-signed.zip",
                     "COAL_CLUSTALO", "COAL_MAFFT", "WSL"):
            self.assertIn(must, text)

    def test_usage_link_resolves(self):
        usage = (ROOT / "docs" / "usage.md").read_text(encoding="utf-8")
        install = (ROOT / "docs" / "install.md").read_text(encoding="utf-8")
        for anchor in re.findall(r"\(install\.md#([\w-]+)\)", usage):
            heads = {re.sub(r"[^\w\- ]", "", h).strip().lower().replace(" ", "-")
                     for h in re.findall(r"(?m)^#+ (.+)$", install)}
            self.assertIn(anchor, heads)


if __name__ == "__main__":
    unittest.main()
