"""`coal dnds` / `coal kaks` / `coal translate` on inputs that cannot work, and on outputs that cannot be written.

  * A run in which no pair has an estimate because of the input (every result is `error` or `no_data`) ends with
    exit status 2 and one line that names the usual cause: a protein file, or sequences that are not aligned. A run
    with at least one estimate, saturated or not-converged result ends with 0 as before (docs/kaks.md "Exit status").
  * A gzip-compressed input is named as such by coal-cpp ("decompress it first"), not read as text; the `coal`
    command decompresses it into a temporary file and runs coal-cpp on that (the same table as for the plain file).
  * A write error names the file and the cause.
Needs build/coal-cpp (COAL_CPP_BIN).
"""
from __future__ import annotations

import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
from test_dnds_pipeline import CPP, needs_cpp  # noqa: E402

EX = ROOT / "examples" / "clean"
ALN = ">a\nATGAAAGTTCTGGCAGGTCCGTTA\n>b\nATGAAGGTTTTGGCAGGACCGCTA\n"
PROT = ">a\nMKVLAGPLEEFQ\n>b\nMKVLSGPLEEFQ\n>c\nMRVLSGPIEEFQ\n"


def cpp(*args, **kw):
    return subprocess.run([str(CPP), *map(str, args)], capture_output=True, text=True, timeout=120,
                          stdin=subprocess.DEVNULL, **kw)


def coal(*args, env=None):
    e = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "COAL_CPP": str(CPP), **(env or {})}
    return subprocess.run([sys.executable, "-m", "coal", *map(str, args)], capture_output=True, text=True,
                          timeout=120, stdin=subprocess.DEVNULL, env=e)


@needs_cpp
class NoEstimate(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def test_protein_input_is_exit_2_with_the_hint(self):
        p = self.tmp / "prot.fasta"
        p.write_text(PROT, encoding="utf-8")
        for sub, args in (("dnds", [p]), ("kaks", ["-i", p]), ("dnds", [p, "-q"])):
            r = cpp(sub, *args)
            self.assertEqual(r.returncode, 2, r.stderr)
            self.assertIn(f"coal {sub}: error: none of the 3 pair(s) has an estimate", r.stderr)
            self.assertIn("looks like protein sequences", r.stderr)
            self.assertIn("coal PROTEIN_ALN CDS.fasta -o codon.fasta", r.stderr)
            self.assertIn("coal pipeline", r.stderr)
        # the table is still written, with every row (dnds) / the header only (kaks leaves such pairs in: no_data)
        r = cpp("dnds", p)
        rows = r.stdout.splitlines()
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(ln.split("\t")[13] == "no_data" for ln in rows[1:]))
        r = cpp("dnds", EX / "protein.aln.fasta", "-q")   # the scout's case: 28 pairs, all no_data
        self.assertEqual(r.returncode, 2)

    def test_unaligned_input_all_pairs(self):
        p = self.tmp / "unaligned.fasta"
        p.write_text(">a\nATGAAAGTTCTGGCA\n>b\nATGAAGGTTTTG\n>c\nATGAAGGTT\n", encoding="utf-8")
        r = cpp("dnds", p)
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("none of the 3 pair(s) has an estimate (3 result(s) error, 0 no_data)", r.stderr)
        self.assertIn("differ in length, so the input is not an alignment", r.stderr)
        self.assertEqual(len(r.stdout.splitlines()), 4)   # header + 3 error rows: nothing is dropped
        r = cpp("kaks", "-i", p)
        self.assertEqual(r.returncode, 2)
        self.assertEqual(len(r.stdout.splitlines()), 1)   # KaKs_Calculator's format leaves invalid pairs out

    def test_some_estimates_is_exit_0_with_a_hint_for_unaligned_pairs(self):
        r = cpp("dnds", EX / "dna.fasta")   # unaligned CDS: 25 of 28 pairs differ in length, 3 happen to match
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("warning: the two sequences of 25 pair(s) differ in length: is the input aligned?", r.stderr)
        r = cpp("dnds", EX / "dna.fasta", "-q")
        self.assertEqual((r.returncode, r.stderr), (0, ""))

    def test_good_and_saturated_inputs_are_exit_0_and_unchanged(self):
        p = self.tmp / "aln.fasta"
        p.write_text(ALN, encoding="utf-8")
        r = cpp("dnds", p)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotIn("differ in length", r.stderr)
        self.assertNotIn("error", r.stderr.replace("0 error", ""))
        sat = self.tmp / "sat.fasta"   # every codon differs at all three positions: no correction possible
        sat.write_text(">a\n" + "AAACCCGGGTTT" * 5 + "\n>b\n" + "CCCGGGTTTAAA" * 5 + "\n", encoding="utf-8")
        r = cpp("dnds", sat, "-m", "NG")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn(r.stdout.splitlines()[1].split("\t")[13], ("saturated", "ok", "warning"))

    def test_empty_table_is_not_an_input_failure(self):
        p = self.tmp / "empty.axt"
        p.write_text("\n", encoding="utf-8")
        r = cpp("dnds", p)
        self.assertEqual(r.returncode, 0, r.stderr)


@needs_cpp
class Gzip(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)
        self.plain = self.tmp / "aln.fasta"
        self.plain.write_text(ALN, encoding="utf-8")
        self.gz = self.tmp / "aln.fasta.gz"
        self.gz.write_bytes(gzip.compress(ALN.encode()))
        self.odd = self.tmp / "compressed_without_suffix"     # found by its first bytes, not by its name
        self.odd.write_bytes(gzip.compress(ALN.encode()))

    def tearDown(self):
        self._td.cleanup()

    def test_coal_cpp_names_the_problem(self):
        for sub, args in (("dnds", [self.gz]), ("kaks", ["-i", self.odd]), ("translate", [self.gz]),
                          ("codon-dnds", ["-i", self.gz, "-o", self.tmp / "o.tsv"])):
            r = cpp(sub, *args)
            self.assertEqual(r.returncode, 1, (sub, r.stderr))
            self.assertIn("is gzip-compressed: decompress it first (gunzip -k", r.stderr)
            self.assertNotIn("control character", r.stderr)

    def test_the_coal_command_reads_it(self):
        want = coal("dnds", self.plain, "-q")
        self.assertEqual(want.returncode, 0, want.stderr)
        for args in (["dnds", self.gz, "-q"], ["dnds", "-i", self.odd, "-q"], ["dnds", f"--in={self.gz}", "-q"]):
            got = coal(*args)
            self.assertEqual((got.returncode, got.stdout), (0, want.stdout), (args, got.stderr))
        out = self.tmp / "o.kaks"
        r = coal("kaks", "-i", self.gz, "-o", out, "-q")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(out.read_text(encoding="utf-8"),
                         coal("kaks", "-i", self.plain, "-o", "-", "-q").stdout)
        cds = self.tmp / "cds.fa.gz"
        cds.write_bytes(gzip.compress(b">g1\nATGAAAGTTTAA\n"))
        r = coal("translate", cds, "-q")
        self.assertEqual((r.returncode, r.stdout), (0, ">g1\nMKV\n"), r.stderr)
        # nothing is left behind next to the input or in the temporary folder of this test
        self.assertEqual(sorted(p.name for p in self.tmp.iterdir()),
                         ["aln.fasta", "aln.fasta.gz", "cds.fa.gz", "compressed_without_suffix", "o.kaks"])

    def test_a_broken_gzip_file_is_one_error_line(self):
        bad = self.tmp / "cut.fa.gz"
        bad.write_bytes(gzip.compress(ALN.encode() * 50)[:60])
        r = coal("dnds", bad)
        self.assertEqual(r.returncode, 1)
        self.assertIn("coal dnds: error:", r.stderr)
        self.assertIn("cut.fa.gz", r.stderr)
        self.assertNotIn("Traceback", r.stderr)

    def test_an_output_is_never_decompressed(self):
        # a file named after -o that happens to be gzip data is an output: replaced, not read
        r = coal("dnds", self.plain, "-o", self.odd, "-q")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(self.odd.read_text(encoding="utf-8").startswith("seq1\tseq2"))


@needs_cpp
class WriteErrors(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)
        self.aln = self.tmp / "aln.fasta"
        self.aln.write_text(ALN, encoding="utf-8")

    def tearDown(self):
        self._td.cleanup()

    def test_missing_folder_and_folder_as_output(self):
        r = cpp("dnds", self.aln, "-o", self.tmp / "no_such_folder" / "o.tsv")
        self.assertEqual(r.returncode, 1)
        self.assertRegex(r.stderr, r"coal dnds: error: cannot write '.*no_such_folder.o\.tsv': \S")
        r = cpp("kaks", "-i", self.aln, "-o", self.tmp)
        self.assertEqual(r.returncode, 1)
        self.assertIn("it is a folder, not a file", r.stderr)

    @unittest.skipUnless(os.path.exists("/dev/full"), "needs /dev/full (Linux)")
    def test_full_device_names_the_file_and_the_cause(self):
        for sub, args in (("dnds", [self.aln, "-o", "/dev/full"]), ("kaks", ["-i", self.aln, "-o", "/dev/full"])):
            r = cpp(sub, *args)
            self.assertEqual(r.returncode, 1, r.stderr)
            self.assertIn(f"coal {sub}: error: write error on '/dev/full': ", r.stderr)
            self.assertRegex(r.stderr, "(?i)no space left")


if __name__ == "__main__":
    unittest.main()
