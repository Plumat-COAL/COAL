"""The R package (packaging/r/coalr): its vendored C++ sources are in sync with cpp/src, its version is COAL's, and
its Makevars compile every vendored source. (The R tests themselves run with R: docs/r.md.)"""
import os
import re
import subprocess
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, "packaging", "r", "coalr")


def read(*parts):
    with open(os.path.join(*parts), encoding="utf-8") as fh:
        return fh.read()


class RBindings(unittest.TestCase):
    def test_vendored_sources_are_in_sync(self):
        r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "sync_r_sources.py"), "--check"],
                           capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_sync_check_and_line_ends(self):
        # A Windows checkout with core.autocrlf gives cpp/src CRLF line ends; without the .gitattributes rule the
        # vendored copies had them too, and --check called every file "out of date" (CI, Windows). cpp/src is read
        # with either line end; the copies are LF (an R source package), and a CRLF copy is named as that.
        import contextlib
        import importlib.util
        import io
        import shutil
        import tempfile
        from pathlib import Path
        spec = importlib.util.spec_from_file_location("sync_r_sources", os.path.join(ROOT, "scripts", "sync_r_sources.py"))
        sync = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(sync)

        def run(*argv):
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                return sync.main(list(argv)), out.getvalue() + err.getvalue()

        with tempfile.TemporaryDirectory() as td:
            src, dst = Path(td, "src"), Path(td, "dst")
            shutil.copytree(sync.SRC, src)
            shutil.copytree(sync.DST, dst)
            for f in list(src.iterdir()) + list(dst.iterdir()):   # LF, whatever this checkout has
                f.write_bytes(f.read_bytes().replace(b"\r\n", b"\n"))
            sync.SRC, sync.DST = src, dst
            self.assertEqual(run("--check")[0], 0)
            for f in src.iterdir():                                # cpp/src as autocrlf checks it out
                f.write_bytes(f.read_bytes().replace(b"\n", b"\r\n"))
            self.assertEqual(run("--check")[0], 0, "CRLF in cpp/src is not a difference")
            for f in dst.iterdir():                                # the copies too: not an R source package
                f.write_bytes(f.read_bytes().replace(b"\n", b"\r\n"))
            rc, text = run("--check")
            self.assertEqual(rc, 1)
            self.assertIn("CRLF line ends", text)
            self.assertNotIn("out of date", text)
            self.assertEqual(run()[0], 0)                          # the sync writes LF
            self.assertEqual(run("--check")[0], 0)
            self.assertFalse(any(b"\r" in f.read_bytes() for f in dst.iterdir()))
            (dst / "core.cpp").write_bytes((dst / "core.cpp").read_bytes() + b"// edited\n")
            rc, text = run("--check")
            self.assertEqual(rc, 1)
            self.assertIn("out of date: core.cpp", text)

    def test_gitattributes_keep_byte_exact_files_lf(self):
        # every folder whose files are compared byte for byte (or hashed) is LF in every checkout
        text = read(ROOT, ".gitattributes")
        for rule in ("tests/fixtures/** text eol=lf", "tests/golden/** text eol=lf", "tests/kaks/** eol=lf",
                     "tests/realdata/** eol=lf", "benchmarks/results/** eol=lf", "packaging/r/** eol=lf",
                     "examples/** text eol=lf", "*.sh eol=lf"):
            self.assertRegex(text, rf"(?m)^{re.escape(rule)}$")
        if os.path.isdir(os.path.join(ROOT, ".git")) or os.path.isfile(os.path.join(ROOT, ".git")):
            for sub in (("packaging", "r", "coalr", "src", "coal"), ("tests", "kaks", "golden"),
                        ("tests", "fixtures", "pipeline"), ("tests", "fixtures", "align")):
                for f in os.listdir(os.path.join(ROOT, *sub)):
                    if not os.path.isfile(os.path.join(ROOT, *sub, f)):
                        continue
                    with open(os.path.join(ROOT, *sub, f), "rb") as fh:
                        self.assertNotIn(b"\r\n", fh.read(), f"{'/'.join(sub)}/{f} has CRLF line ends (a checkout from before "
                                         f".gitattributes kept it LF: delete the folder and `git checkout -- "
                                         f"{'/'.join(sub)}`)")

    def test_version_is_coals(self):
        version = re.search(r'^version\s*=\s*"([^"]+)"', read(ROOT, "pyproject.toml"), re.M).group(1)
        self.assertRegex(read(PKG, "DESCRIPTION"), rf"(?m)^Version: {re.escape(version)}$")
        self.assertIn(f'"{version}"', read(PKG, "src", "coal", "coal_version.h"))

    def test_makevars_compile_every_vendored_source(self):
        cpps = sorted(f[:-4] for f in os.listdir(os.path.join(PKG, "src", "coal")) if f.endswith(".cpp"))
        self.assertEqual(cpps, ["align", "core", "dnds", "kaks", "tables"])
        for name in ("Makevars", "Makevars.win"):
            text = read(PKG, "src", name)
            for c in cpps:
                self.assertIn(f"coal/{c}.o", text, name)
            self.assertIn("CXX_STD = CXX17", text)

    def test_no_output_or_exit_in_the_vendored_core(self):
        # R packages must not write to stdout/stderr or end the process from compiled code
        for f in os.listdir(os.path.join(PKG, "src", "coal")):
            if f.endswith((".cpp", ".hpp")):
                text = read(PKG, "src", "coal", f)
                for bad in ("std::cout", "std::cerr", "exit(", "abort(", "printf("):
                    if bad == "printf(":
                        text = text.replace("snprintf(", "")
                    self.assertNotIn(bad, text, f"{f}: {bad}")

    def test_licence_and_exports(self):
        self.assertRegex(read(PKG, "DESCRIPTION"), r"(?m)^License: GPL-3$")  # version 3 only, as COAL
        self.assertFalse(os.path.exists(os.path.join(PKG, "LICENSE")), "GPL-3 needs no LICENSE file in an R package")
        # one string that looks like a path but names no file is "file not found" (documented heuristic)
        self.assertIn("looks_like_path", read(PKG, "R", "utils.R"))
        self.assertIn("file not found", read(PKG, "R", "utils.R"))
        self.assertIn("Sequence inputs", read(PKG, "man", "coalr-package.Rd"))
        self.assertIn("file not found", read(ROOT, "docs", "r.md"))
        ns = read(PKG, "NAMESPACE")
        for fn in ("coal_dnds_pairs", "coal_dnds_clusters", "coal_kaks", "coal_codon_align", "coal_align_pairs",
                   "coal_translate", "coal_pipeline_pairs", "coal_pipeline_clusters"):
            self.assertIn(f"export({fn})", ns)


if __name__ == "__main__":
    unittest.main()
