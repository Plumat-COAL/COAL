"""packaging/bioconda: the recipe as Bioconda will read it (what can be checked without conda).

  * meta.yaml is valid YAML once its Jinja is resolved; selectors are well formed.
  * A compiled, non-noarch recipe: `python` unpinned in host and run, Pythons below 3.10 skipped by a selector,
    ARM Linux and Apple silicon listed as additional platforms.
  * The test command that failed (scout K1): `pal2nal` without arguments exits 1 with the usage on stderr and
    nothing on stdout. The recipe's three lines for it are run here, in a shell under `set -e` as conda-build does.
  * build.sh installs the two programs only: no static library, headers or CMake package.
A real `conda build` / `bioconda-utils lint` is not run here (packaging/bioconda/README.md says what was).
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RECIPE = ROOT / "packaging" / "bioconda"
SELECTOR = re.compile(r"#\s*\[([^\]]+)\]\s*$")


def rendered(version: str = "0.6.0") -> str:
    """meta.yaml with its Jinja resolved the way conda-build does for this recipe's four constructs."""
    text = (RECIPE / "meta.yaml").read_text(encoding="utf-8")
    sets = dict(re.findall(r'\{%\s*set\s+(\w+)\s*=\s*"([^"]*)"\s*%\}', text))
    sets["version"] = version
    text = re.sub(r"\{%.*?%\}", "", text)
    text = re.sub(r"\{\{\s*compiler\('cxx'\)\s*\}\}", "cxx-compiler", text)
    text = re.sub(r"\{\{\s*stdlib\('c'\)\s*\}\}", "c-stdlib", text)
    text = re.sub(r"\{\{\s*pin_subpackage\([^)]*\)\s*\}\}", sets["name"], text)
    text = re.sub(r"\{\{\s*(\w+)\s*\}\}", lambda m: sets[m.group(1)], text)
    assert "{{" not in text and "{%" not in text
    return text


class Recipe(unittest.TestCase):
    def setUp(self):
        self.raw = (RECIPE / "meta.yaml").read_text(encoding="utf-8")

    def load(self):
        try:
            import yaml
        except ImportError:
            self.skipTest("PyYAML is not installed")
        return yaml.safe_load(rendered())

    def test_valid_yaml_and_selectors(self):
        d = self.load()
        self.assertEqual((d["package"]["name"], str(d["package"]["version"])), ("coal-align", "0.6.0"))
        for line in self.raw.splitlines():
            m = SELECTOR.search(line)
            if m and not line.lstrip().startswith("#"):
                # a selector is a Python expression over conda-build's names
                for word in re.findall(r"[A-Za-z_]\w*", m.group(1)):
                    self.assertIn(word, {"win", "linux", "osx", "unix", "or", "and", "not", "py", "aarch64", "arm64"}
                                  | {w for w in [word] if re.fullmatch(r"py\d*", w)}, line)
                compile(m.group(1), "selector", "eval")

    def test_non_noarch_python_handling(self):
        d = self.load()
        self.assertNotIn("noarch", d["build"])
        self.assertIn("python", d["requirements"]["host"])
        self.assertIn("pip", d["requirements"]["host"])
        self.assertEqual(d["requirements"]["run"], ["python"])
        for sec in ("host", "run"):
            for dep in d["requirements"][sec]:
                self.assertNotRegex(dep, r"^python\s*[<>=]", f"{sec}: a compiled recipe leaves python unpinned")
        self.assertRegex(self.raw, r"(?m)^  skip: True  # \[win or py<310\]$")
        self.assertEqual(d["build"]["skip"], True)
        self.assertEqual(d["extra"]["additional-platforms"], ["linux-aarch64", "osx-arm64"])
        self.assertIn("cxx-compiler", d["requirements"]["build"])
        # the floor of the recipe is the floor of the package
        need = re.search(r'requires-python = ">=3\.(\d+)"', (ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        self.assertEqual(need.group(1), "10")
        st = re.search(r'requires = \["setuptools>=([\d.]+)"\]', (ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        self.assertIn(f"setuptools >={st.group(1)}", d["requirements"]["host"])

    def test_only_the_programs_are_installed(self):
        sh = (RECIPE / "build.sh").read_text(encoding="utf-8")
        self.assertNotIn("\r", sh)
        code = "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))
        self.assertNotIn("cmake --install", code)
        self.assertIn("-DCOAL_INSTALL=OFF", code)
        for prog in ("coal-cpp", "coal-cpp-align"):
            self.assertIn(f'install -m 0755 build-conda/{prog} "${{PREFIX}}/bin/{prog}"', code)
        self.assertNotIn("libcoal_core", code)
        cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("option(COAL_INSTALL", cmake)
        self.assertIn("if(COAL_INSTALL)", cmake)

    def test_every_test_file_is_there_and_lf(self):
        d = self.load()
        for name in d["test"]["files"]:
            data = (RECIPE / name).read_bytes()
            self.assertNotIn(b"\r", data, name)

    @unittest.skipIf(os.name == "nt" or shutil.which("bash") is None, "the recipe's tests are POSIX shell lines")
    def test_the_pal2nal_usage_lines_pass_as_written(self):
        """The three lines of the recipe, with `pal2nal` the drop-in of this tree, under `set -e`."""
        d = self.load()
        cmds = d["test"]["commands"]
        k = next(i for i, c in enumerate(cmds) if c.startswith("pal2nal > usage_stdout.txt"))
        lines = cmds[k:k + 3]
        self.assertEqual(lines[1:], ["test ! -s usage_stdout.txt", 'grep -q "pal2nal pep.aln nuc.fasta" usage.txt'])
        self.assertNotIn('pal2nal > usage.txt && grep -q "pal2nal pep.aln nuc.fasta" usage.txt', cmds)
        with tempfile.TemporaryDirectory() as tmp:
            shim = Path(tmp) / "pal2nal"
            shim.write_text(f'#!/bin/sh\nexec "{sys.executable}" -m coal_compat "$@"\n', encoding="utf-8", newline="\n")
            shim.chmod(0o755)
            env = {**os.environ, "PATH": tmp + os.pathsep + os.environ.get("PATH", ""),
                   "PYTHONPATH": str(ROOT / "src")}
            script = "set -ex\n" + "\n".join(lines) + "\necho passed\n"
            r = subprocess.run(["bash", "-c", script], cwd=tmp, env=env, capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("passed", r.stdout)
            # the line as it was fails, which is why it was replaced
            old = 'set -ex\npal2nal > usage.txt && grep -q "pal2nal pep.aln nuc.fasta" usage.txt\n'
            r = subprocess.run(["bash", "-c", old], cwd=tmp, env=env, capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 1)
            # and the next line of the recipe: -h is exit 0 with the help on stdout
            r = subprocess.run(["bash", "-c", "set -ex\n" + cmds[k + 3]], cwd=tmp, env=env, capture_output=True,
                               text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)


if __name__ == "__main__":
    unittest.main()
