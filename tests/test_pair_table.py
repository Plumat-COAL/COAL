"""One pair table, written in four places: its header must be the same everywhere.

The columns of dnds.tsv are spelled out in cpp/src/dnds_cli.cpp (`coal codon-dnds`, the pipeline's engine),
cpp/src/kaks_cli.cpp (`coal dnds`: the same table without the two protein-alignment columns, which only a run that
aligned proteins has), src/coal/dnds_pipeline.py (PAIR_COLUMNS) and packaging/r/coalr/R/utils.R (PAIR_COLUMNS of the R
package, whose C++ side fills the data frame by these names). A column added, renamed or reordered in one of them and
not in the others would give tables that parse into the wrong fields; this module reads the four sources, and asks the
built coal-cpp too when there is one (COAL_CPP_BIN).
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
from coal import dnds_pipeline as P  # noqa: E402

PROTEIN_ONLY = ("prot_aln_length", "prot_identity")   # not in `coal dnds` (its input is already a codon alignment)


def _c_string(text: str, start: str) -> str:
    """The C string literal (adjacent pieces joined) that follows `start` in a C++ source."""
    at = text.index(start) + len(start)
    pieces = []
    while True:
        m = re.compile(r'\s*"((?:[^"\\]|\\.)*)"').match(text, at)
        if not m:
            break
        pieces.append(m.group(1))
        at = m.end()
    return "".join(pieces).replace("\\t", "\t").replace("\\n", "\n")


def cpp_bin() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


class TestPairTableHeader(unittest.TestCase):
    def setUp(self):
        self.want = list(P.PAIR_COLUMNS)
        self.tidy = [c for c in self.want if c not in PROTEIN_ONLY]

    def test_python_columns_are_the_documented_ones(self):
        self.assertEqual(self.want[:3], ["seq1", "seq2", "method"])
        self.assertEqual(self.want[-2:], ["status", "note"])
        for c in PROTEIN_ONLY:
            self.assertIn(c, self.want)
        doc = (ROOT / "docs" / "pipeline.md").read_text(encoding="utf-8")
        for c in self.want:
            self.assertIn(f"`{c}`", doc, f"docs/pipeline.md does not describe the column {c}")

    def test_codon_dnds_source(self):
        src = (ROOT / "cpp" / "src" / "dnds_cli.cpp").read_text(encoding="utf-8")
        self.assertEqual(_c_string(src, "        fout <<").rstrip("\n").split("\t"), self.want)
        # ... and its --help names the same columns in the same order
        help_text = _c_string(src, "const char* CODON_DNDS_HELP =")
        listed = re.search(r"the pair table: \[NAME\] (.*?)\n\s+--group-column", help_text, re.S).group(1).split()
        self.assertEqual(listed, self.want)

    def test_coal_dnds_source(self):
        src = (ROOT / "cpp" / "src" / "kaks_cli.cpp").read_text(encoding="utf-8")
        self.assertEqual(_c_string(src, "const char* TIDY_HEADER =").rstrip("\n").split("\t"), self.tidy)
        help_text = _c_string(src, "const char* DNDS_HELP =")
        listed = re.search(r"one row per pair and method:\n(.*?)\n", help_text).group(1).split()
        self.assertEqual(listed, self.tidy)

    def test_coalr_source(self):
        src = (ROOT / "packaging" / "r" / "coalr" / "R" / "utils.R").read_text(encoding="utf-8")
        block = re.search(r"(?m)^PAIR_COLUMNS <- c\((.*?)\)$", src, re.S).group(1)
        self.assertEqual(re.findall(r'"([^"]+)"', block), self.want)
        # the vendored engine is cpp/src (scripts/sync_r_sources.py --check, tests/test_r_bindings.py), and the
        # R glue names every column it fills
        glue = (ROOT / "packaging" / "r" / "coalr" / "src" / "coalr.cpp").read_text(encoding="utf-8")
        for c in self.want[2:]:   # (seq1 and seq2 are member numbers there: m1, m2; R puts the IDs in)
            self.assertTrue(f'_["{c}"]' in glue, f"packaging/r/coalr/src/coalr.cpp never fills the column {c}")

    def test_the_built_programs_write_these_headers(self):
        cpp = cpp_bin()
        if cpp is None:
            if os.environ.get("COAL_REQUIRE_CPP"):
                self.fail("build/coal-cpp (COAL_CPP_BIN) is not built (COAL_REQUIRE_CPP is set)")
            self.skipTest("build/coal-cpp (COAL_CPP_BIN) is not built")
        with tempfile.TemporaryDirectory() as t:
            aln = Path(t) / "pair.fasta"
            aln.write_text(">a\nATGAAAGTTCTGGCA\n>b\nATGAAGGTTTTGGCA\n", encoding="utf-8")
            r = subprocess.run([str(cpp), "dnds", str(aln), "-q"], capture_output=True, timeout=60)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(r.stdout.decode().splitlines()[0].split("\t"), self.tidy)
            out = Path(t) / "pairs.tsv"
            groups = b"G\tg1\tall\nM\ta\tATGAAAGTTCTGGCA\nM\tb\tATGAAGGTTTTGGCA\nE\n"
            for extra, lead in (([], []), (["--group-column", "cluster_id"], ["cluster_id"])):
                r = subprocess.run([str(cpp), "codon-dnds", "-o", str(out), "-q", *extra], input=groups,
                                   capture_output=True, timeout=60)
                self.assertEqual(r.returncode, 0, r.stderr)
                lines = out.read_text(encoding="utf-8").splitlines()
                self.assertEqual(lines[0].split("\t"), lead + self.want)
                self.assertEqual(len(lines[1].split("\t")), len(lead) + len(self.want))


if __name__ == "__main__":
    unittest.main()
