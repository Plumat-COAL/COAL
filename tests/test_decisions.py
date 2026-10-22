"""Design decisions: PAL2NAL-compatible default vs --strict, --quiet, closed-pipe exit codes, and the output formats,
drop-in, parallel/block writing and PAL2NAL fetching.

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


# 40-residue protein; DNA = 2 nt of UTR + a ~90% NNN-masked copy (a reviewer's input)
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


class TestDec1DefaultIsPal2nal(Tmp):
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


class TestDec1Strict(Tmp):
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


class TestDec2Quiet(Tmp):
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


class TestDec3ClosedPipe(Tmp):
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


# CLUSTAL input with PAL2NAL '#' block marks under each block; expected outputs from the real PAL2NAL v14
MARKED_ALN = ("CLUSTAL W (1.83) multiple sequence alignment\n\n\ns1      MKV-ALWB\ns2      MKVSALWY\n        ### ## #\n\n"
              "s1      KDE\ns2      KDE\n         ##\n")
MARKED_NUC = ">s1\nATGAAAGTTGCTCTGTGGAACAAAGATGAA\n>s2\nATGAAGGTCTCAGCACTATGGTATAAGGACGAG\n"
MARKED_PAL2NAL = {
    ("-output", "clustal"): ("CLUSTAL W multiple sequence alignment\n\n"
                             "s1            ATGAAAGTT---GCTCTGTGGAACAAAGATGAA\n"
                             "s2            ATGAAGGTCTCAGCACTATGGTATAAGGACGAG\n"
                             "              #########   ######   ###   ######\n\n"),
    ("-output", "codon"): ("              M   K   V   -   A   L   W   B   K   D   E\n"
                           "s1            ATG AAA GTT --- GCT CTG TGG AAC AAA GAT GAA\n"
                           "              M   K   V   S   A   L   W   Y   K   D   E\n"
                           "s2            ATG AAG GTC TCA GCA CTA TGG TAT AAG GAC GAG\n"
                           "              ### ### ###     ### ###     ###     ### ###\n\n"),
    ("-output", "fasta", "-blockonly"): ">s1\nATGAAAGTTGCTCTGAACGATGAA\n>s2\nATGAAGGTCGCACTATATGACGAG\n",
    ("-output", "clustal", "-nogap"): ("CLUSTAL W multiple sequence alignment\n\n"
                                       "s1            ATGAAAGTTGCTCTGTGGAACAAAGATGAA\n"
                                       "s2            ATGAAGGTCGCACTATGGTATAAGGACGAG\n"
                                       "              ###############   ###   ######\n\n"),
    ("-output", "clustal", "-blockonly"): ("CLUSTAL W multiple sequence alignment\n\n"
                                           "s1            ATGAAAGTTGCTCTGAACGATGAA\n"
                                           "s2            ATGAAGGTCGCACTATATGACGAG\n\n"),
    ("-output", "paml", "-nomismatch"): ("   2     30\ns1\nATGAAAGTT---GCTCTGTGGAAAGATGAA\n"
                                         "s2\nATGAAGGTCTCAGCACTATGGAAGGACGAG\n"),
}


class TestDec4FormatsAndDropIn(Tmp):
    def setUp(self):
        super().setUp()
        (self.d / "pep.aln").write_text(MARKED_ALN)
        (self.d / "nuc.fa").write_text(MARKED_NUC)

    def test_marked_clustal_input_in_both_commands(self):
        for opts, want in MARKED_PAL2NAL.items():
            for module in ("coal", "coal_compat"):
                with self.subTest(opts=opts, module=module):
                    r = run(module, ["pep.aln", "nuc.fa", *opts], self.d)
                    self.assertEqual((r.returncode, r.stdout), (0, want), r.stderr)

    def test_drop_in_defaults_to_clustal_like_pal2nal(self):
        fx = ROOT / "tests" / "fixtures" / "clean"
        gold = (ROOT / "tests" / "golden" / "clean" / "clustal.stdout").read_text()
        r = run("coal_compat", [str(fx / "pep.fasta"), str(fx / "nuc.fasta")], self.d)
        self.assertEqual((r.returncode, r.stdout), (0, gold), r.stderr)
        r = coal([str(fx / "pep.fasta"), str(fx / "nuc.fasta")], self.d)  # coal's own default stays fasta
        self.assertEqual(r.stdout, (ROOT / "tests" / "golden" / "clean" / "fasta.stdout").read_text())

    def test_codon_with_column_removal_is_an_error_as_in_pal2nal(self):
        for opt in ("-nogap", "-nomismatch", "-blockonly"):
            for module in ("coal", "coal_compat"):
                with self.subTest(opt=opt, module=module):
                    r = run(module, ["pep.aln", "nuc.fa", "-output", "codon", opt], self.d)
                    self.assertEqual((r.returncode, r.stdout), (1, ""))
                    self.assertIn("-output codon cannot be combined", r.stderr)

    def test_html_is_a_clear_error_and_nostderr_is_quiet(self):
        r = run("coal_compat", ["pep.aln", "nuc.fa", "-html"], self.d)
        self.assertEqual((r.returncode, r.stdout), (1, ""))
        self.assertIn("-html", r.stderr)
        self.assertIn("not supported", r.stderr)
        r = run("coal_compat", ["pep.aln", "nuc.fa", "-nostderr"], self.d)
        self.assertEqual((r.returncode, r.stdout, r.stderr), (0, MARKED_PAL2NAL[("-output", "clustal")], ""))

    def test_unknown_residue_column_removed_by_nomismatch(self):  # PAL2NAL v14: ATGAAAGTTCTG
        (self.d / "b.fa").write_text(">a\nMKBVL\n>b\nMKAVL\n")
        (self.d / "bn.fa").write_text(">a\nATGAAAGATGTTCTG\n>b\nATGAAAGCTGTTCTG\n")
        r = coal(["b.fa", "bn.fa", "-nomismatch"], self.d)
        self.assertEqual((r.returncode, r.stdout), (0, ">a\nATGAAAGTTCTG\n>b\nATGAAAGTTCTG\n"), r.stderr)

    def test_help_lists_the_formats_for_the_app(self):
        h = coal(["--help"], self.d).stdout
        self.assertIn("{fasta,clustal,paml,codon}", h)
        self.assertIn("-blockonly", h)


class TestDec4ParallelAndBlocks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import random
        sys.path.insert(0, str(SRC))
        from coal.tables import get_code
        cls._td = tempfile.TemporaryDirectory(prefix="coal_dec_j_")
        d = Path(cls._td.name)
        rng = random.Random(4)
        code = get_code(1).code
        sense = sorted(c for c, a in code.items() if a != "*")
        n, L = 40, 150
        cods = [["ATG"] + [rng.choice(sense) for _ in range(L - 1)] for _ in range(n)]
        rows = []
        for k, cs in enumerate(cods):
            row = [code[c] for c in cs]
            for i in rng.sample(range(1, L), 10):
                row[i] = "-"
            if k % 7 == 3:
                row[5] = "4"
            rows.append("".join(row))
        dna = []
        for k, (cs, row) in enumerate(zip(cods, rows)):
            out = []
            for c, ch in zip(cs, row):
                out.append("" if ch == "-" else ("ACGT" if ch == "4" else c))
            if k % 9 == 2:
                i = next(i for i in range(20, L) if row[i].isalpha())
                out[i] = "TGG" if row[i] != "W" else "TTT"  # a mismatch
            dna.append("CC" + "".join(out) + "TAA")
        mask = "".join(rng.choice("# ") for _ in range(L))
        text = "CLUSTAL W (1.83) multiple sequence alignment\n\n\n"
        for b in range(0, L, 60):
            for k, row in enumerate(rows):
                text += f"seq{k:<6}{row[b:b + 60]}\n"
            text += " " * 9 + mask[b:b + 60] + "\n\n"
        (d / "p.aln").write_text(text)
        (d / "n.fa").write_text("".join(f">seq{k}\n{s}\n" for k, s in enumerate(dna)))
        cls.d = d

    @classmethod
    def tearDownClass(cls):
        cls._td.cleanup()

    def _run(self, fmt, opts, jobs, memory=None):
        import io
        from contextlib import redirect_stdout
        from unittest import mock
        from coal import pipeline
        out = self.d / f"o_{fmt}_{jobs}_{'_'.join(opts)}_{memory}.txt"
        kw = dict(nogap="-nogap" in opts, nomismatch="-nomismatch" in opts, blockonly="-blockonly" in opts)
        with mock.patch.object(pipeline, "BATCH_NT", 3000), \
                mock.patch.object(pipeline, "BLOCK_MEMORY", memory or pipeline.BLOCK_MEMORY):
            pipeline.run(str(self.d / "p.aln"), [str(self.d / "n.fa")], out_path=str(out), jobs=jobs,
                         min_nt_per_job=1, out_format=fmt, **kw)
        return out.read_text()

    def test_jobs_1_2_4_identical_for_every_format_and_option(self):
        for fmt in ("fasta", "clustal", "paml", "codon"):
            for opts in ((), ("-nogap",), ("-nomismatch",), ("-blockonly",), ("-nogap", "-nomismatch", "-blockonly")):
                if fmt == "codon" and opts:
                    continue
                with self.subTest(fmt=fmt, opts=opts):
                    one = self._run(fmt, opts, 1)
                    self.assertTrue(one)
                    self.assertEqual(self._run(fmt, opts, 2), one)
                    self.assertEqual(self._run(fmt, opts, 4), one)

    def test_block_writer_memory_limit_does_not_change_output(self):
        for fmt in ("clustal", "codon"):
            with self.subTest(fmt=fmt):
                self.assertEqual(self._run(fmt, (), 1, memory=1), self._run(fmt, (), 1))


class TestDec4FetchPal2nal(unittest.TestCase):
    def setUp(self):
        import hashlib
        import importlib.util
        spec = importlib.util.spec_from_file_location("fetch_pal2nal", ROOT / "scripts" / "fetch_pal2nal.py")
        self.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.m)
        self._td = tempfile.TemporaryDirectory(prefix="coal_fetch_")
        self.m.DEST = Path(self._td.name) / "PAL2NAL"
        self.fake = b"#!/usr/bin/perl\n# stand-in for pal2nal.pl\n"
        self.m.SCRIPT_SHA256 = hashlib.sha256(self.fake).hexdigest()
        self.script = self.m.DEST / "pal2nal.v14" / "pal2nal.pl"

    def tearDown(self):
        self._td.cleanup()

    def _main(self, download):
        import contextlib
        import io
        from unittest import mock
        with mock.patch.object(self.m, "_download", download), contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()) as err:
            try:
                return self.m.main([]), err.getvalue()
            except SystemExit as e:
                return e.code, err.getvalue()

    def test_falls_back_to_github_when_the_official_site_fails(self):
        def download(url):
            if url == self.m.URL:
                raise OSError("403 Forbidden")
            self.assertEqual(url, self.m.GITHUB_URL)
            return self.fake
        rc, err = self._main(download)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.script.read_bytes(), self.fake)
        self.assertIn("403", err)

    def test_wrong_checksum_from_github_is_refused(self):
        def download(url):
            if url == self.m.URL:
                raise OSError("unreachable")
            return self.fake + b"tampered\n"
        rc, err = self._main(download)
        self.assertNotEqual(rc, 0)
        self.assertIn("checksum mismatch", str(rc))
        self.assertFalse(self.script.exists())

    def test_official_tarball_still_used_first(self):
        import hashlib
        import io
        import tarfile
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tar:
            info = tarfile.TarInfo("pal2nal.v14/pal2nal.pl")
            info.size = len(self.fake)
            tar.addfile(info, io.BytesIO(self.fake))
        data = buf.getvalue()
        self.m.ARCHIVE_SHA256 = hashlib.sha256(data).hexdigest()
        seen = []

        def download(url):
            seen.append(url)
            return data
        rc, err = self._main(download)
        self.assertEqual((rc, seen), (0, [self.m.URL]), err)
        self.assertEqual(self.script.read_bytes(), self.fake)


if __name__ == "__main__":
    unittest.main()
