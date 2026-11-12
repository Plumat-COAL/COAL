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
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)

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
        ns = read(PKG, "NAMESPACE")
        for fn in ("coal_dnds_pairs", "coal_dnds_clusters", "coal_kaks", "coal_codon_align", "coal_align_pairs",
                   "coal_translate", "coal_pipeline_pairs", "coal_pipeline_clusters"):
            self.assertIn(f"export({fn})", ns)


if __name__ == "__main__":
    unittest.main()
