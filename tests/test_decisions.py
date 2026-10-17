"""DEC-1 .. DEC-4: the four decisions of overnight/integration/CHANGES.md ("Decisions").

Expected outputs marked "PAL2NAL v14:" were produced by running the real PAL2NAL v14 on the same input.
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


def run(module, args, cwd, timeout=TIMEOUT):
    env = dict(os.environ, PYTHONPATH=str(SRC))
    return subprocess.run([sys.executable, "-m", module, *args], cwd=cwd, env=env, capture_output=True, text=True,
                          timeout=timeout, stdin=subprocess.DEVNULL)


def coal(args, cwd):
    return run("coal", args, cwd)


# 40-residue protein; DNA = 2 nt of UTR + a ~90% NNN-masked copy (reviewer's R3-E1 input)
MASKED_PROT = "MTVDCLAHRFISDLVKSPSFGSVLACYVISPLDRNVCGQH"
MASKED_DNA = ("GCATGNNNNNNNNNNNNTTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATCAGTNNNNNNNNNNNN"
              "NNNNNNTGTNNNNNNNNN")
# PAL2NAL v14: the first 120 nt of the DNA (shifted by the 2-nt UTR), exit 0, a warning per codon
MASKED_PAL2NAL = (">g\nGCATGNNNNNNNNNNNNTTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
                  "NNNNNNNNNNNNNNNNNNNNNNNNNNATCAGTNNNNNNNNNNNNNNNNNNTGTNNNNNNN\n")
# Unrelated protein and DNA of the same length: PAL2NAL v14 writes the DNA unchanged (codons by position)
UNREL_PROT = "MKVALWYFHIKLDEPQRSTGNCMKVALWYFHI"
UNREL_DNA = "ATGCCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT"
# Partial ID match: PAL2NAL v14 pairs by order
PART_PEP = ">geneA\nMKVAL\n>geneB\nMKVAL\n"
PART_NUC = ">geneB\nATGAAAGTTGCTCTG\n>geneX\nATGAAGGTCGCACTA\n"
PART_PAL2NAL = ">geneA\nATGAAAGTTGCTCTG\n>geneB\nATGAAGGTCGCACTA\n"


class Tmp(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory(prefix="coal_dec_")
        self.d = Path(self._td.name)
        (self.d / "mp.fa").write_text(f">g\n{MASKED_PROT}\n")
        (self.d / "mn.fa").write_text(f">g\n{MASKED_DNA}\n")
        (self.d / "up.fa").write_text(f">u\n{UNREL_PROT}\n")
        (self.d / "un.fa").write_text(f">u\n{UNREL_DNA}\n")
        (self.d / "ip.fa").write_text(PART_PEP)
        (self.d / "in.fa").write_text(PART_NUC)

    def tearDown(self):
        self._td.cleanup()


class TestDec1DefaultIsPal2nal(Tmp):  # DEC-1
    def test_unrelated_pair_gives_pal2nal_output_with_a_notice(self):
        r = coal(["up.fa", "un.fa"], self.d)
        self.assertEqual((r.returncode, r.stdout), (0, f">u\n{UNREL_DNA[:60]}\n{UNREL_DNA[60:]}\n"), r.stderr)
        self.assertIn("u: protein and DNA are probably unrelated", r.stderr)
        self.assertIn("--strict makes this an error", r.stderr)
        self.assertIn("pepAlnPos", r.stderr)  # per-codon warnings as PAL2NAL

    def test_masked_gene_with_utr_gives_pal2nal_output_with_a_notice(self):
        r = coal(["mp.fa", "mn.fa"], self.d)
        self.assertEqual((r.returncode, r.stdout), (0, MASKED_PAL2NAL), r.stderr)
        self.assertIn("g: position of the gene in the DNA is a guess", r.stderr)
        self.assertIn("may be shifted", r.stderr)

    def test_partial_id_match_pairs_by_order_and_says_strict(self):
        r = coal(["ip.fa", "in.fa"], self.d)
        self.assertEqual((r.returncode, r.stdout), (0, PART_PAL2NAL), r.stderr)
        self.assertIn("only 1 of 2 protein IDs are DNA IDs", r.stderr)
        self.assertIn("--strict makes this an error", r.stderr)

    def test_library_default_is_not_strict(self):
        sys.path.insert(0, str(SRC))
        from coal.core import codon_align
        rows, warns = codon_align([("u", UNREL_PROT)], [("u", UNREL_DNA)])
        self.assertEqual(rows, [("u", UNREL_DNA)])
        self.assertTrue(any("probably unrelated" in w for w in warns))


class TestDec1Strict(Tmp):  # DEC-1
    def test_strict_spellings_turn_each_rule_into_an_error(self):
        for flag in ("--strict", "-strict"):
            for pep, nuc, text in (("up.fa", "un.fa", "u: protein and DNA are inconsistent (30 of 32 checkable"),
                                   ("mp.fa", "mn.fa", "too few to tell where the gene starts"),
                                   ("ip.fa", "in.fa", "only 1 of 2 protein IDs are DNA IDs")):
                with self.subTest(flag=flag, pep=pep):
                    r = coal([pep, nuc, flag], self.d)
                    self.assertEqual((r.returncode, r.stdout), (1, ""), r.stderr)
                    self.assertIn("coal: error:", r.stderr)
                    self.assertIn(text, r.stderr)

    def test_strict_keeps_good_input_unchanged(self):
        (self.d / "p.fa").write_text(">s1\nMKV-AL\n>s2\nMKVSAL\n")
        (self.d / "n.fa").write_text(">s1\nATGAAAGTTGCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n")
        r = coal(["p.fa", "n.fa", "--strict"], self.d)
        self.assertEqual((r.returncode, r.stdout, r.stderr),
                         (0, ">s1\nATGAAAGTT---GCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n", ""))

    def test_help_names_the_rules(self):
        r = coal(["--help"], self.d)
        h = " ".join(r.stdout.split())
        for phrase in ("--strict", "more than half", "checkable codons", "10-residue", "protein IDs"):
            self.assertIn(phrase, h)


class TestDec2Quiet(Tmp):  # DEC-2
    def test_quiet_silences_warnings_and_notices_but_not_stdout(self):
        for pep, nuc in (("up.fa", "un.fa"), ("mp.fa", "mn.fa"), ("ip.fa", "in.fa")):
            loud = coal([pep, nuc], self.d)
            for q in ("-q", "--quiet"):
                with self.subTest(pep=pep, q=q):
                    r = coal([pep, nuc, q], self.d)
                    self.assertEqual((r.returncode, r.stdout, r.stderr), (0, loud.stdout, ""))
                    self.assertNotEqual(loud.stderr, "")

    def test_errors_are_always_printed(self):
        r = coal(["up.fa", "un.fa", "-q", "--strict"], self.d)
        self.assertEqual((r.returncode, r.stdout), (1, ""))
        self.assertIn("coal: error: u: protein and DNA are inconsistent", r.stderr)
        r = coal(["up.fa", "missing.fa", "-q"], self.d)
        self.assertEqual(r.returncode, 1)
        self.assertIn("coal: error:", r.stderr)

    def test_help_says_errors_are_still_printed(self):
        h = " ".join(coal(["--help"], self.d).stdout.split())
        self.assertIn("-q, --quiet do not print warnings or notices (errors are always printed)", h)

    def test_nostderr_is_quiet_in_the_drop_in(self):
        r = run("coal_compat", ["up.fa", "un.fa", "-output", "fasta", "-nostderr"], self.d)
        self.assertEqual((r.returncode, r.stderr), (0, ""))
        self.assertTrue(r.stdout.startswith(">u\n"))


class TestDec3ClosedPipe(Tmp):  # DEC-3
    class _Closed(__import__("io").StringIO):
        def __init__(self, err):
            super().__init__()
            self.err = err

        def write(self, s):
            raise self.err

        def fileno(self):
            raise __import__("io").UnsupportedOperation("no fileno")

    def _main(self, module, windows, err):
        import contextlib
        import errno
        import io
        from unittest import mock
        sys.path.insert(0, str(SRC))
        from coal import pipeline
        main = __import__(module + ".cli", fromlist=["main"]).main
        e = io.StringIO()
        cwd = os.getcwd()
        os.chdir(self.d)
        try:
            with mock.patch.object(pipeline, "_WINDOWS", windows), \
                    mock.patch.object(sys, "stdout", self._Closed(err)), contextlib.redirect_stderr(e):
                args = ["up.fa", "un.fa"] + (["-output", "fasta"] if module == "coal_compat" else [])
                rc = main(args)
        finally:
            os.chdir(cwd)
        return rc, e.getvalue()

    def test_posix_broken_pipe_is_141_and_silent_in_both_commands(self):
        import errno
        for module in ("coal", "coal_compat"):
            with self.subTest(module=module):
                self.assertEqual(self._main(module, False, BrokenPipeError(errno.EPIPE, "Broken pipe"))[0], 141)

    def test_windows_einval_is_141_too(self):
        import errno
        for module in ("coal", "coal_compat"):
            with self.subTest(module=module):
                rc, err = self._main(module, True, OSError(errno.EINVAL, "Invalid argument"))
                self.assertEqual(rc, 141)
                self.assertNotIn("Traceback", err)
                self.assertNotIn("coal: error", err)


if __name__ == "__main__":
    unittest.main()
