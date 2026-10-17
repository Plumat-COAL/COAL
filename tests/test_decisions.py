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


if __name__ == "__main__":
    unittest.main()
