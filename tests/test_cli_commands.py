"""`coal --help` names every command, `coal help [COMMAND]` works, and the Ka/Ks commands print their help without
coal-cpp.

`coal kaks`, `coal dnds`, `coal translate` and `coal codon-dnds` are run by coal-cpp, which an install from the source
package does not have. Their --help must still work there: src/coal/_cpp_help.py holds a copy of coal-cpp's help texts.
The copy is generated from the C++ sources (`python tests/test_cli_commands.py --write`), and this module checks it
against the sources and, when coal-cpp is built, against what the binary prints.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import cli  # noqa: E402
from test_pair_table import _c_string, cpp_bin  # noqa: E402

SOURCES = {"kaks": ("kaks_cli.cpp", "const char* KAKS_HELP ="), "dnds": ("kaks_cli.cpp", "const char* DNDS_HELP ="),
           "translate": ("dnds_cli.cpp", "const char* TRANSLATE_HELP ="),
           "codon-dnds": ("dnds_cli.cpp", "const char* CODON_DNDS_HELP =")}
MODULE = ROOT / "src" / "coal" / "_cpp_help.py"


def source_help() -> dict[str, str]:
    out = {}
    for name, (file, start) in SOURCES.items():
        text = (ROOT / "cpp" / "src" / file).read_text(encoding="utf-8")
        out[name] = _c_string(text, start).replace('\\"', '"').replace("\\\\", "\\")
    return out


def module_text() -> str:
    lines = ['"""The --help texts of the commands coal-cpp runs (kaks, dnds, translate, codon-dnds), for an install',
             "without coal-cpp. Generated from cpp/src/kaks_cli.cpp and cpp/src/dnds_cli.cpp by",
             "`python tests/test_cli_commands.py --write`; do not edit by hand (tests/test_cli_commands.py compares",
             'this copy with the sources and with the built binary)."""', "", "HELP = {"]
    for name, text in source_help().items():
        lines.append(f"    {name!r}: (")
        lines += [f"        {ln!r}" for ln in text.splitlines(keepends=True)]
        lines.append("    ),")
    lines.append("}")
    return "\n".join(lines) + "\n"


def coal(*args, env=None):
    e = {**os.environ, "PYTHONPATH": str(ROOT / "src"), **(env or {})}
    return subprocess.run([sys.executable, "-m", "coal", *map(str, args)], capture_output=True, text=True,
                          timeout=120, stdin=subprocess.DEVNULL, env=e, encoding="utf-8")


def no_cpp_env(tmp: str) -> dict:
    """An environment in which `coal` finds no coal-cpp: no COAL_CPP, nothing on PATH. (A binary inside the
    package or next to the interpreter would still be found; the callers skip then.)"""
    return {"COAL_CPP": "", "PATH": tmp}


class HelpCopy(unittest.TestCase):
    def test_the_copy_is_the_sources_text(self):
        self.assertEqual(MODULE.read_text(encoding="utf-8"), module_text(),
                         "src/coal/_cpp_help.py is stale: python tests/test_cli_commands.py --write")
        from coal._cpp_help import HELP
        self.assertEqual(HELP, source_help())
        for name, text in HELP.items():
            self.assertTrue(text.startswith(f"usage: coal {name} "), name)
            self.assertGreater(len(text.splitlines()), 10, name)

    def test_the_binary_prints_the_same(self):
        cpp = cpp_bin()
        if cpp is None:
            if os.environ.get("COAL_REQUIRE_CPP"):
                self.fail("build/coal-cpp (COAL_CPP_BIN) is not built (COAL_REQUIRE_CPP is set)")
            self.skipTest("build/coal-cpp (COAL_CPP_BIN) is not built")
        from coal._cpp_help import HELP
        for name, text in HELP.items():
            r = subprocess.run([str(cpp), name, "--help"], capture_output=True, timeout=60)
            self.assertEqual(r.returncode, 0, name)
            self.assertEqual(r.stdout.decode("utf-8").replace("\r\n", "\n"), text, name)


class Commands(unittest.TestCase):
    def test_help_lists_every_command(self):
        r = coal("--help")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("Other commands (`coal COMMAND --help`):", r.stdout)
        for name, text in cli.COMMANDS:
            self.assertRegex(r.stdout, rf"(?m)^  coal {re.escape(name)} +{re.escape(text)}$")
        self.assertIn("coal --cpp-path", r.stdout)
        self.assertIn("Convert a protein alignment + unaligned coding DNA into a codon alignment.", r.stdout)
        # every word main() dispatches on is in the list
        src = (ROOT / "src" / "coal" / "cli.py").read_text(encoding="utf-8")
        body = src[src.index("def main("):]
        dispatched = set(re.findall(r'argv\[0\] == "([a-z-]+)"', body))
        dispatched |= set(re.findall(r'"([a-z-]+)"', re.search(r"argv\[0\] in \(([^)]*)\)", body).group(1)))
        self.assertEqual(dispatched, {name for name, _ in cli.COMMANDS})

    def test_coal_help(self):
        self.assertEqual(coal("help").stdout, coal("--help").stdout)
        r = coal("help", "pipeline")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("usage: coal pipeline", r.stdout)
        r = coal("help", "pipeline", "pairs")
        self.assertIn("--query", r.stdout)
        self.assertIn("--longest-isoform", r.stdout)
        self.assertIn("--min-codons", r.stdout)
        r = coal("help", "cite")
        self.assertIn("usage: coal cite", r.stdout)
        r = coal("help", "nonsense")
        self.assertEqual(r.returncode, 2)
        self.assertIn("is not a command of coal: kaks, dnds, pipeline", r.stderr)
        self.assertEqual(coal("pipeline", "--help").returncode, 0)

    def test_kaks_and_dnds_help_without_coal_cpp(self):
        from coal import dispatch
        from coal._cpp_help import HELP
        with tempfile.TemporaryDirectory() as t:
            env = no_cpp_env(t)
            probe = subprocess.run([sys.executable, "-c", "from coal import dispatch; print(dispatch.find_cpp())"],
                                   capture_output=True, text=True, env={**os.environ, "PYTHONPATH": str(ROOT / "src"),
                                                                        **env})
            if probe.stdout.strip() != "None":
                self.skipTest(f"a coal-cpp is found without COAL_CPP and PATH ({probe.stdout.strip()})")
            for name in HELP:
                for flag in ("--help", "-h"):
                    r = coal(name, flag, env=env)
                    self.assertEqual(r.returncode, 0, (name, r.stderr))
                    self.assertEqual(r.stdout, HELP[name], name)
                    self.assertIn("note: this command cannot run here: coal-cpp was not found", r.stderr)
                    self.assertIn("platform wheels", r.stderr)
            # running it says what this install is and how to get the other one
            aln = Path(t) / "a.fasta"
            aln.write_text(">a\nATGAAA\n>b\nATGAAG\n", encoding="utf-8")
            r = coal("dnds", aln, env=env)
            self.assertEqual(r.returncode, 1)
            for must in ("coal dnds: needs coal-cpp, which was not found", "this install has the Python engine only",
                         "codon alignments", "coal_align-<version>-py3-none-<platform>.whl",
                         "https://github.com/plumat-coal/COAL/releases/latest", "docs/install.md"):
                self.assertIn(must, r.stderr)
            self.assertNotIn("pip install coal-align", r.stderr.replace("python -m pip install` that file", ""))
            r = coal("pipeline", "pairs", "--query", aln, "--subject", aln, "-o", Path(t) / "o", env=env)
            self.assertEqual(r.returncode, 1)
            self.assertIn("coal pipeline needs coal-cpp", r.stderr)
            self.assertIn("this install has the Python engine only", r.stderr)
            # the codon alignment itself works in such an install
            ex = ROOT / "examples" / "clean"
            r = coal(ex / "protein.aln.fasta", ex / "dna.fasta", env=env)
            self.assertEqual(r.returncode, 0, r.stderr)

    def test_no_document_tells_users_to_pip_install_from_pypi_as_if_it_were_there(self):
        """coal-align is not on PyPI yet: a bare `pip install coal-align` in the documentation or in a message
        sends users to a 404. It may appear only where the text says that this is for later ("Once it is on
        PyPI"), next to "PyPI" in the same paragraph."""
        files = [ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md")), *sorted((ROOT / "site").glob("*.html")),
                 *sorted((ROOT / "packaging" / "r" / "coalr" / "R").glob("*.R")),
                 *sorted((ROOT / "packaging" / "r" / "coalr" / "man").glob("*.Rd")),
                 ROOT / "packaging" / "r" / "coalr" / "README.md", *sorted((ROOT / "src" / "coal").glob("*.py"))]
        bare = re.compile(r"pip install (?:--upgrade )?coal-align(?![\w\[-]| @)")
        for f in files:
            if not f.is_file() or f.name.startswith("stress-test"):
                continue
            text = f.read_text(encoding="utf-8")
            for para in re.split(r"\n\s*\n", text):
                if bare.search(para):
                    self.assertRegex(para, r"(?i)\bPyPI\b", f"{f.relative_to(ROOT)}: `pip install coal-align` "
                                     "without saying that it is not on PyPI yet")


if __name__ == "__main__":
    if "--write" in sys.argv:
        MODULE.write_text(module_text(), encoding="utf-8", newline="\n")
        print(f"wrote {MODULE}")
    else:
        unittest.main()
