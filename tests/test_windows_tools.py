"""Clustal Omega and MAFFT on native Windows (WIN-TOOLS, 2026-11-14), and the small texts that came with it.

No automatic install was added for Windows: neither Bioconda nor conda-forge has a win-64 package of clustalo or
mafft, and COAL adds a download only with a SHA-256 computed from the file itself, which could not be done for
MSYS2's Clustal Omega or MAFFT's own Windows zip (src/coal/msa.py, the comment above NO_WINDOWS_CLUSTALO). What is
checked here: on Windows the installers download and run nothing, whatever is asked, and their messages (and the
pipeline's) name what works there, with the version COAL pins. Also: the Windows-path hint of `coal`'s missing-file
message (U17), the app's "needs coal-cpp" text in a development version, and the D9 wording of the help. Everything
is mocked, so it runs on every OS."""
from __future__ import annotations

import argparse
import errno
import io
import os
import re
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from coal import msa, seqio  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402


class WindowsCase(unittest.TestCase):
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
                  mock.patch.object(msa, "_platform", return_value=("Windows", "AMD64")),
                  mock.patch.object(msa.platform, "system", return_value="Windows")):
            p.start()
            self.addCleanup(p.stop)

    def _refuse(self, name):
        def f(*a, **k):
            self.calls.append((name, a))
            raise AssertionError(f"{name} was called on Windows: {a}")
        return f


class TestNothingIsDownloadedOnWindows(WindowsCase):
    def test_clustalo_every_method(self):
        for method in ("auto", "conda", "micromamba"):
            for force in (False, True):
                with self.subTest(method=method, force=force), self.assertRaises(msa.MsaError) as cm:
                    msa.install_clustalo(method=method, force=force, quiet=True)
                self.assertEqual(str(cm.exception), msa.NO_WINDOWS_CLUSTALO)
        self.assertEqual(self.calls, [])

    def test_mafft_every_method(self):
        for method in ("auto", "conda", "micromamba"):
            for force in (False, True):
                with self.subTest(method=method, force=force), self.assertRaises(msa.MsaError) as cm:
                    msa.install_mafft(method=method, force=force, quiet=True)
                self.assertEqual(str(cm.exception), msa.NO_WINDOWS_MAFFT)
        self.assertEqual(self.calls, [])

    def test_the_install_command(self):
        for tool, text in (("clustalo", msa.NO_WINDOWS_CLUSTALO), ("mafft", msa.NO_WINDOWS_MAFFT)):
            err = io.StringIO()
            with self.subTest(tool=tool), redirect_stderr(err):
                rc = P.main(["install", tool, "-q"])
            self.assertEqual(rc, 1)
            self.assertIn(text, err.getvalue())
        self.assertEqual(self.calls, [])

    def test_no_windows_download_is_pinned(self):
        # what would have to change to add one: an entry with its SHA-256 (never without); there is none today
        self.assertNotIn(("Windows", "AMD64"), msa._MICROMAMBA)
        for tool in ("clustalo", "mafft"):
            self.assertIsNone(msa.package_pin(tool))
        src = (ROOT / "src" / "coal" / "msa.py").read_text(encoding="utf-8")
        self.assertNotRegex(src, r"https?://\S*\.zip")       # no Windows archive address in the download tables


class TestWhatTheMessagesSay(unittest.TestCase):
    def test_clustal_omega(self):
        t = msa.NO_WINDOWS_CLUSTALO
        for must in ("no automatic Clustal Omega install on Windows", "Bioconda and conda-forge", "WSL",
                     "coal pipeline install clustalo", "MSYS2", "https://www.msys2.org",
                     "pacman -S mingw-w64-ucrt-x86_64-clustal-omega", "COAL_CLUSTALO",
                     "C:\\msys64\\ucrt64\\bin\\clustalo.exe", "PATH", "--aligner mafft", "COAL_MAFFT", "mafft.bat",
                     "not tested"):
            self.assertIn(must, t)
        self.assertIn(msa.CLUSTALO_VERSION, t)       # MSYS2 packages the version COAL pins elsewhere
        self.assertNotIn("clustal" + ".org", t)      # (no address of the former Clustal site anywhere)

    def test_mafft(self):
        t = msa.NO_WINDOWS_MAFFT
        for must in ("no automatic MAFFT install on Windows", "Bioconda and conda-forge",
                     f"mafft-{msa.MAFFT_VERSION}-win64-signed.zip",
                     "https://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html", "LF line ends",
                     "COAL_MAFFT", "mafft-win\\mafft.bat", "not tested with COAL", "WSL",
                     "coal pipeline install mafft", "publishes none"):
            self.assertIn(must, t)
        self.assertEqual(msa.MAFFT_VERSION, "7.526")   # the zip named above is that version's (MAFFT's page)

    def test_exported(self):
        self.assertIn("NO_WINDOWS_MAFFT", msa.__all__)
        self.assertIn("NO_WINDOWS_CLUSTALO", msa.__all__)

    def test_install_help_says_what_windows_gets(self):
        out = io.StringIO()
        with mock.patch("sys.stdout", out), self.assertRaises(SystemExit):
            P.main(["install", "--help"])
        text = " ".join(out.getvalue().split())
        self.assertIn("Clustal Omega and MAFFT: Linux and macOS", text)
        self.assertIn("on Windows the command says what works there instead", text)
        self.assertIn("(conda packages only)", text)
        self.assertNotIn("Bioconda only", text)       # MAFFT comes from conda-forge


class TestPipelineOnWindows(WindowsCase):
    """`coal pipeline clusters` with no aligner found on Windows: the same text as the installer."""

    def run_tool(self, tool, install_tools):
        r = P._Run.__new__(P._Run)
        r.tools, r.quiet = {}, True
        r.a = argparse.Namespace(clustalo=None, mafft=None, install_tools=install_tools)
        with mock.patch.object(msa, "find_clustalo", return_value=None), \
                mock.patch.object(msa, "find_mafft", return_value=None), \
                mock.patch.object(P.platform, "system", return_value="Windows"), \
                mock.patch.dict(os.environ, {"COAL_CLUSTALO": "", "COAL_MAFFT": ""}), \
                self.assertRaises(P.PipelineError) as cm:
            r.msa_tool(tool)
        return str(cm.exception)

    def test_not_found(self):
        for install_tools in (False, True):
            with self.subTest(install_tools=install_tools):
                self.assertIn(msa.NO_WINDOWS_CLUSTALO, self.run_tool("clustalo", install_tools))
                self.assertIn(msa.NO_WINDOWS_MAFFT, self.run_tool("mafft", install_tools))
        self.assertEqual(self.calls, [])


class TestWindowsPathHint(unittest.TestCase):
    """U17 for `coal PEP NUC` (the pipeline had it already): a Windows path missing on Linux / macOS names its WSL
    form. Python engine; coal-cpp's own message is unchanged (a proposal for its owner)."""

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
