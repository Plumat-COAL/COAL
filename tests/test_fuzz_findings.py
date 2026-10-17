"""FZ-1 .. FZ-6: findings of the independent differential fuzz against the real PAL2NAL v14
(overnight/integration/CHANGES.md). Expected outputs marked "PAL2NAL v14:" were produced by running the real PAL2NAL.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
TIMEOUT = 60


def run(module, args, cwd):
    env = dict(os.environ, PYTHONPATH=str(SRC))
    return subprocess.run([sys.executable, "-m", module, *args], cwd=cwd, env=env, capture_output=True,
                          timeout=TIMEOUT, stdin=subprocess.DEVNULL)


class Tmp(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory(prefix="coal_fz_")
        self.d = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def write(self, name, data):
        (self.d / name).write_bytes(data if isinstance(data, bytes) else data.encode())

    def coal(self, *args, module="coal"):
        r = run(module, list(args), self.d)
        return r.returncode, r.stdout.decode("latin-1"), r.stderr.decode("utf-8", "replace")


NUC = ">a\nATGAAAGTT\n>b\nATGCGT\n"
OUT = ">a\nATGAAAGTT\n>b\nATGCGT---\n"  # PAL2NAL v14 for every input below


class TestFZ1ClustalHashAndBlankLines(Tmp):  # FZ-1
    def test_hash_line_after_block_is_a_mark_line_not_a_sequence(self):
        self.write("p.aln", "CLUSTAL W\n\na    MKV\nb    MR-\n#    ###\n")
        self.write("n.fa", NUC)
        self.assertEqual(self.coal("p.aln", "n.fa")[:2], (0, OUT))
        # PAL2NAL v14: the '#' line is the mark line under the block
        rc, out, err = self.coal("p.aln", "n.fa", "-output", "clustal")
        self.assertEqual((rc, out), (0, "CLUSTAL W multiple sequence alignment\n\na             ATGAAAGTT\n"
                                        "b             ATGCGT---\n              #########\n\n"), err)

    def test_hash_line_inside_a_block_is_skipped(self):
        self.write("p.aln", "CLUSTAL W\n\na    MKV\n# c  ###\nb    MR-\n")
        self.write("n.fa", NUC)
        rc, out, err = self.coal("p.aln", "n.fa", "-output", "clustal")
        self.assertEqual((rc, out), (0, "CLUSTAL W multiple sequence alignment\n\na             ATGAAAGTT\n"
                                        "b             ATGCGT---\n              #########\n\n"), err)

    def test_clustal_after_leading_blank_lines(self):
        self.write("p.aln", "\n\nCLUSTAL W\n\na    MKV\nb    MR-\n")
        self.write("n.fa", NUC)
        self.assertEqual(self.coal("p.aln", "n.fa")[:2], (0, OUT))

    def test_app_reads_the_same_ids(self):
        self.write("p.aln", "\nCLUSTAL W\n\na    MKV\nb    MR-\n#    ###\n")
        sys.path.insert(0, str(SRC))
        from coal_app.seqformat import read_ids
        self.assertEqual(read_ids(str(self.d / "p.aln")), ["a", "b"])


class TestFZ2NonAsciiHeaderText(Tmp):  # FZ-2
    def test_non_ascii_description_in_protein_header(self):
        self.write("p.fa", b">a Kr\xc3\xbcppel-like factor 4\nMKV\n>b\nMR-\n")
        self.write("n.fa", NUC)
        rc, out, err = self.coal("p.fa", "n.fa")
        self.assertEqual((rc, out), (0, OUT), err)

    def test_non_ascii_id_is_a_clear_error(self):
        self.write("p.fa", b">\xc3\xbc1 x\nMKV\n>b\nMR-\n")
        self.write("n.fa", NUC)
        rc, out, err = self.coal("p.fa", "n.fa")
        self.assertEqual((rc, out), (1, ""))
        self.assertIn("p.fa:1: non-ASCII character in the sequence ID", err)

    def test_non_ascii_in_sequence_still_an_error(self):
        self.write("p.fa", b">a\nMK\xc3\xbcV\n>b\nMR-\n")
        self.write("n.fa", NUC)
        rc, out, err = self.coal("p.fa", "n.fa")
        self.assertEqual(rc, 1)
        self.assertIn("p.fa:2: non-ASCII character", err)


if __name__ == "__main__":
    unittest.main()
