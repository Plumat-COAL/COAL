"""Protein files with two defects: coal-cpp must report the same one as the Python engine.

Found by a differential fuzz of the two engines. The Python reader (seqio.iter_fasta) yields a record before it
checks the '>' line that ended it, so pass 1's row-length check sees the record first; and the pipeline parses a
CLUSTAL file (read_clustal) before it reads the '#' mark line (clustal_mask). Each case runs `python -m coal` and
coal-cpp on the same files and asserts identical exit code, stdout and stderr.

Skipped if the binary is not built (COAL_CPP_BIN, default ./build/coal-cpp).
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
from test_differential import cpp_binary  # noqa: E402

SRC = ROOT / "src"
TIMEOUT = 60
DNA = ">a\nATGAAA\n>b\nATGAAA\n>a\nATGAAA\n"


def _run(cmd, cwd):
    env = dict(os.environ, PYTHONPATH=str(SRC))
    r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, timeout=TIMEOUT, stdin=subprocess.DEVNULL)
    return r.returncode, r.stdout, r.stderr


class TestCppErrorPrecedence(unittest.TestCase):
    def setUp(self):
        self.binary = cpp_binary()
        if self.binary is None:
            self.skipTest("coal-cpp not built")
        self._td = tempfile.TemporaryDirectory(prefix="coal_prec_")
        self.d = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def both(self, pep: str, data: bytes, expect: str):
        (self.d / pep).write_bytes(data)
        (self.d / "n.fasta").write_text(DNA)
        py = _run([sys.executable, "-m", "coal", pep, "n.fasta"], self.d)
        cpp = _run([str(self.binary), pep, "n.fasta"], self.d)
        self.assertEqual(py[0], 1, py[2])
        self.assertIn(expect, py[2].decode("utf-8", "replace"))
        self.assertEqual(cpp, py)

    # finding 1: the duplicate / empty header is checked only after the pending row is handed back
    def test_duplicate_header_after_short_row(self):
        self.both("p3.fasta", b">a\nMKV\n>b\nMK\n>a\nMK\n", "rows have different lengths (3 and 2, at b)")

    def test_empty_header_after_short_row(self):
        self.both("pe.fasta", b">a\nMKV\n>b\nMK\n>\nMK\n", "rows have different lengths (3 and 2, at b)")

    def test_bad_header_ending_the_file(self):
        self.both("pt.fasta", b">a\nMKV\n>b\nMK\n>a\n", "rows have different lengths (3 and 2, at b)")
        self.both("pz.fasta", b">a\nMKV\n>\n", "pz.fasta:3: empty FASTA header")

    def test_duplicate_header_alone_still_reported(self):
        self.both("pd.fasta", b">a\nMKV\n>b\nMKV\n>a\nMKV\n", "duplicate sequence ID 'a'")

    def test_non_ascii_header_is_reported_at_read_time(self):
        # both engines check the bytes of a line as it is read, before the previous record is handed back
        self.both("pn.fasta", b">a\nMKV\n>b\nMK\n>g\xc3\xa9\nMK\n", "pn.fasta:5: non-ASCII character")

    # finding 2: CLUSTAL rows are parsed before the '#' mark line is read (whole file, with the non-ASCII check)
    def test_malformed_clustal_line_before_non_ascii(self):
        self.both("c1.aln", b"CLUSTAL W (1.83) multiple sequence alignment\n\na\ng\xc3\xa9 MK\n",
                  "c1.aln:3: malformed CLUSTAL line 'a'")

    def test_malformed_clustal_id_only_lines_then_non_ascii_id(self):
        self.both("c0.aln", b"CLUSTAL W (1.83) multiple sequence alignment\n\na\nb\ng\xc3\xa9\n",
                  "c0.aln:3: malformed CLUSTAL line 'a'")

    def test_clustal_non_ascii_alone_still_reported(self):
        self.both("cn.aln", b"CLUSTAL W (1.83) multiple sequence alignment\n\na MKV\ng\xc3\xa9 MKV\n",
                  "cn.aln:4: non-ASCII character")


if __name__ == "__main__":
    unittest.main()
