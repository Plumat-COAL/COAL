"""The public Python API as docs/python-api.md documents it: every example on that page, run.

`import coal` exports codon_align, align_files, CoalError and __version__; the examples' outputs are the documented
ones. A change here is a change to the documented API: update docs/python-api.md with it."""
from __future__ import annotations

import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import coal  # noqa: E402
from coal import CoalError, align_files, codon_align  # noqa: E402

DOC = ROOT / "docs" / "python-api.md"
EXAMPLE = ROOT / "examples" / "one_broken_gene"

# The gene of the "one wrong codon" example: column 11 is H, its codon GAT (D)
PROT = "MFWDKESRSPHESAPQYARKIWEMAAAVAP"
DNA = "ATGTTCTGGGACAAAGAAAGCAGAAGCCCAGATGAAAGCGCACCACAATACGCAAGAAAAATATGGGAAATGGCAGCAGCAGTAGCACCA"


class TestExports(unittest.TestCase):
    def test_all(self):
        self.assertEqual(coal.__all__, ["CoalError", "codon_align", "align_files", "__version__"])
        for name in coal.__all__:
            self.assertTrue(hasattr(coal, name), name)
        from coal.pipeline import run
        self.assertIs(align_files, run)
        self.assertTrue(issubclass(CoalError, Exception))

    def test_doc_names_every_export(self):
        text = DOC.read_text(encoding="utf-8")
        for name in coal.__all__:
            self.assertIn(f"coal.{name}", text, f"docs/python-api.md does not document coal.{name}")
        self.assertIn(f"COAL {coal.__version__}", text, "docs/python-api.md names another version than the package")


class TestCodonAlignExamples(unittest.TestCase):
    def test_pairing_by_id_any_order(self):
        rows, warnings = codon_align([("g1", "MK-VL"), ("g2", "MKQVL")],
                                     [("g2", "ATGAAACAGGTTCTG"), ("g1", "ATGAAAGTTCTGTAA")])
        self.assertEqual(rows, [("g1", "ATGAAA---GTTCTG"), ("g2", "ATGAAACAGGTTCTG")])
        self.assertEqual(warnings, [])

    def test_one_wrong_codon(self):
        rows, warnings = codon_align([("x", PROT)], [("x", DNA)])
        self.assertEqual(rows, [("x", DNA)])
        self.assertEqual(warnings, ["x pepAlnPos 11: H does not correspond to GAT"])
        rows, warnings = codon_align([("x", PROT)], [("x", DNA)], nomismatch=True)
        self.assertEqual((len(rows[0][1]), warnings), (87, []))
        self.assertEqual(rows[0][1], DNA[:30] + DNA[33:])

    def test_warning_formats(self):
        # B/Z/J/O/U
        rows, warnings = codon_align([("a", "MUK")], [("a", "ATGTGAAAA")])
        self.assertEqual(rows, [("a", "ATGTGAAAA")])
        self.assertEqual(warnings, ["a pepAlnPos 2: U (selenocysteine) is not checked; any codon is accepted"])
        # pairing by order with a partial ID match: a pairing warning first
        rows, warnings = codon_align([("a", "MK"), ("b", "MK")], [("a", "ATGAAA"), ("c", "ATGAAG")])
        self.assertEqual(rows, [("a", "ATGAAA"), ("b", "ATGAAG")])
        self.assertEqual(len(warnings), 1)
        self.assertTrue(warnings[0].startswith("only 1 of 2 protein IDs are DNA IDs; pairing sequences by order instead"))
        # no ID matches: silent, by order
        self.assertEqual(codon_align([("a", "MK"), ("b", "MK")], [("x", "ATGAAA"), ("y", "ATGAAG")]),
                         ([("a", "ATGAAA"), ("b", "ATGAAG")], []))
        # normalised IDs where order would mispair
        rows, warnings = codon_align([("AT1.1", "MK"), ("AT2.1", "MV")], [("AT2", "ATGGTT"), ("AT1", "ATGAAA")])
        self.assertEqual(rows, [("AT1.1", "ATGAAA"), ("AT2.1", "ATGGTT")])
        self.assertTrue(warnings[0].startswith("2 of 2 protein IDs match a DNA ID only after normalising"))

    def test_options(self):
        self.assertEqual(codon_align([("a", "mk")], [("a", "atgaaa")]), ([("a", "atgaaa")], []))  # case kept
        self.assertEqual(codon_align([("a", "M2K")], [("a", "ATGCCAAA")]), ([("a", "ATGCC-AAA")], []))
        self.assertEqual(codon_align([("a", "MK*"), ("b", "M-*")], [("a", "ATGAAATAA"), ("b", "ATGTAG")], nogap=True),
                         ([("a", "ATG"), ("b", "ATG")], []))
        self.assertEqual(codon_align([("a", "MK")], [("a", "ATGAAA")], table=11), ([("a", "ATGAAA")], []))

    def test_errors(self):
        cases = [
            ("row_lengths", None, [("a", "MK"), ("b", "M")], [("a", "ATGAAA"), ("b", "ATG")], {}),
            ("count_mismatch", None, [("a", "MK")], [("b", "ATGAAA"), ("c", "ATGAAA")], {}),
            ("no_sequences", None, [], [], {}),
            ("partial_id_match", None, [("a", "MK"), ("b", "MK")], [("a", "ATGAAA"), ("c", "ATGAAG")],
             {"strict": True}),
            ("invalid_protein", "a", [("a", "M@")], [("a", "ATGAAA")], {}),
            ("empty", "a", [("a", "--")], [("a", "ATG")], {}),
            ("inconsistent", "a", [("a", "MK")], [("a", "CCCCCC")], {}),
        ]
        for code, seq_id, p, d, kw in cases:
            with self.assertRaises(CoalError) as cm:
                codon_align(p, d, **kw)
            e = cm.exception
            self.assertEqual((e.code, e.seq_id, e.index), (code, seq_id, None), str(e))
        with self.assertRaises(CoalError) as cm:
            codon_align([("a", "MK")], [("a", "CCCCCC")])
        self.assertEqual(str(cm.exception), "a: protein and DNA are inconsistent (no match found)")

    def test_unknown_table_is_a_coal_error(self):  # issue #18: not tables.get_code's ValueError
        from coal.tables import SUPPORTED
        for bad in (7, 0, -1, 99, "1"):
            with self.subTest(table=bad):
                with self.assertRaises(CoalError) as cm:
                    codon_align([("a", "MK")], [("a", "ATGAAA")], table=bad)
                e = cm.exception
                self.assertEqual((e.code, e.seq_id, e.index), ("unknown_table", None, None))
                self.assertNotIsInstance(e, ValueError)
                for n in SUPPORTED:  # the message lists the supported codes
                    self.assertIn(str(n), str(e))
        self.assertEqual(codon_align([("a", "MK")], [("a", "ATGAAA")], table=SUPPORTED[-1]), ([("a", "ATGAAA")], []))

    def test_duplicate_ids_are_refused(self):  # as the command line does
        with self.assertRaises(CoalError) as cm:
            codon_align([("a", "MK"), ("a", "MK")], [("a", "ATGAAA"), ("b", "ATGAAA")])
        self.assertEqual((cm.exception.code, str(cm.exception)),
                         ("duplicate_id", "the protein alignment: duplicate sequence ID 'a': every sequence needs its own ID"))
        with self.assertRaises(CoalError) as cm:
            codon_align([("a", "MK")], [("a", "ATGAAA"), ("a", "ATGAAA")])
        self.assertEqual((cm.exception.code, str(cm.exception)),
                         ("duplicate_id", "the DNA: duplicate sequence ID 'a': every sequence needs its own ID"))

    def test_wrong_types_are_type_errors(self):  # a clear TypeError, not a KeyError / AttributeError from inside
        ok_p, ok_d = [("a", "MK")], [("a", "ATGAAA")]
        cases = [
            ({"a": "MK"}, ok_d, "proteins must be a list of (id, sequence) pairs, not dict"),
            ((x for x in ok_p), ok_d, "proteins must be a list of (id, sequence) pairs, not generator"),
            ("MK", ok_d, "proteins must be a list of (id, sequence) pairs, not str"),
            (ok_p, None, "dnas must be a list of (id, sequence) pairs, not NoneType"),
            ([(1, "MK")], ok_d, "proteins[0]: the ID must be a str, not int"),
            ([("a", None)], ok_d, "proteins[0] ('a'): the sequence must be a str, not NoneType"),
            ([("a", "MK", "x")], ok_d, "proteins[0] must be an (id, sequence) pair, not a tuple of 3"),
            (["MK"], ok_d, "proteins[0] must be an (id, sequence) pair, not str"),
            (ok_p, [("a", b"ATG")], "dnas[0] ('a'): the sequence must be a str, not bytes"),
        ]
        for p, d, msg in cases:
            with self.subTest(msg=msg):
                with self.assertRaises(TypeError) as cm:
                    codon_align(p, d)
                self.assertEqual(str(cm.exception), msg)
        self.assertEqual(codon_align(tuple(ok_p), [["a", "ATGAAA"]]), ([("a", "ATGAAA")], []))  # tuples / lists


class TestAlignFilesArguments(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        (self.tmp / "p.fa").write_text(">s1\nMKV-AL\n>s2\nMKVSAL\n")
        (self.tmp / "n.fa").write_text(">s1\nATGAAAGTTGCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_one_dna_path_as_a_string_or_pathlike(self):
        out = self.tmp / "o.fa"
        for nuc in (str(self.tmp / "n.fa"), self.tmp / "n.fa", [self.tmp / "n.fa"], (str(self.tmp / "n.fa"),)):
            with self.subTest(nuc=type(nuc).__name__):
                out.unlink(missing_ok=True)
                self.assertEqual(align_files(self.tmp / "p.fa", nuc, out_path=str(out)), 0)
                self.assertEqual(out.read_text(), ">s1\nATGAAAGTT---GCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n")

    def test_wrong_types_are_type_errors(self):
        p, n = str(self.tmp / "p.fa"), str(self.tmp / "n.fa")
        cases = [
            ((p, 5), "nuc_paths must be a path or a list of paths, not int"),
            ((p, {n: 1}), "nuc_paths must be a path or a list of paths, not dict"),
            ((p, (x for x in [n])), "nuc_paths must be a path or a list of paths, not generator"),
            ((None, n), "pep_path must be a path (str or os.PathLike), not NoneType"),
            ((p, [n, 3]), "nuc_paths[1] must be a path (str or os.PathLike), not int"),
            ((p, []), "nuc_paths must name at least one DNA file"),
        ]
        for args, msg in cases:
            with self.subTest(msg=msg):
                with self.assertRaises(TypeError) as cm:
                    align_files(*args, out_path=str(self.tmp / "o.fa"))
                self.assertEqual(str(cm.exception), msg)
        self.assertFalse((self.tmp / "o.fa").exists())


class TestAlignFilesExample(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        for f in ("protein.aln.fasta", "dna.fasta"):
            shutil.copy(EXAMPLE / f, self.tmp)
        self.cwd = os.getcwd()
        os.chdir(self.tmp)

    def tearDown(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_keep_going_with_report(self):
        from coal.report import RunReport
        pep, dna, out = "protein.aln.fasta", ["dna.fasta"], "codon.fasta"
        messages = []
        rec = RunReport([], pep, dna, out, "fasta", 1, False, False, keep_going=True)
        skipped = align_files(pep, dna, out_path=out, keep_going=True, warn=messages.append,
                              notice=messages.append, rec=rec, jobs=1)
        rec.exit_code = 2 if skipped else 0
        report = rec.to_dict()
        self.assertEqual(skipped, 1)
        self.assertEqual(messages, ["skipped OG0001_sp14: protein and DNA are inconsistent (no match found)"])
        self.assertEqual(report["summary"], {"total": 20, "ok": 19, "warning": 0, "failed": 1, "pairing": "id"})
        self.assertEqual(report["sequences"][13]["messages"][0]["code"], "frameshift")
        self.assertEqual(report["schema"], "coal.report/1")
        self.assertEqual(report["exit_code"], 2)
        self.assertEqual(rec.summary_lines(), [
            "1 sequence of 20 could not be converted and was left out (first: OG0001_sp14: protein and DNA are "
            "inconsistent (no match found). Frameshift at column 102: one base is missing in the DNA, so from here on "
            "the codons are read out of step.)"])
        with open(out, encoding="ascii") as fh:
            self.assertEqual(sum(line.startswith(">") for line in fh), 19)

    def test_without_keep_going_raises_and_writes_nothing(self):
        with self.assertRaises(CoalError) as cm:
            align_files("protein.aln.fasta", ["dna.fasta"], out_path="codon.fasta", jobs=1)
        e = cm.exception
        self.assertEqual((e.code, e.seq_id, e.index), ("inconsistent", "OG0001_sp14", 13))
        self.assertFalse(os.path.exists("codon.fasta"))
        self.assertEqual([f for f in os.listdir(".") if f.startswith(".coal")], [])

    def test_unknown_table_is_a_coal_error(self):  # issue #18
        with self.assertRaises(CoalError) as cm:
            align_files("protein.aln.fasta", ["dna.fasta"], out_path="codon.fasta", table=7, jobs=1)
        self.assertEqual(cm.exception.code, "unknown_table")
        self.assertFalse(os.path.exists("codon.fasta"))

    def test_input_format_is_a_coal_error(self):  # issue #18: InputFormatError is a CoalError (and a ValueError)
        from coal.seqio import InputFormatError
        self.assertTrue(issubclass(InputFormatError, CoalError))
        self.assertTrue(issubclass(InputFormatError, ValueError))
        with open("notes.txt", "w", encoding="ascii") as fh:
            fh.write("hello world\nthis is not FASTA\n")
        for pep, nuc in (("notes.txt", ["dna.fasta"]), ("protein.aln.fasta", ["notes.txt"])):
            with self.subTest(pep=pep, nuc=nuc):
                with self.assertRaises(CoalError) as cm:
                    align_files(pep, nuc, out_path="codon.fasta", jobs=1)
                e = cm.exception
                self.assertIsInstance(e, InputFormatError)
                self.assertEqual((e.code, e.seq_id, e.index), ("input_format", None, None))
                self.assertTrue(str(e).startswith("notes.txt:1: sequence data before the first '>' header"), str(e))
        self.assertFalse(os.path.exists("codon.fasta"))

    def test_fail_records_the_error(self):
        from coal.report import RunReport
        rec = RunReport([], "protein.aln.fasta", ["dna.fasta"], "codon.fasta", "fasta", 1, False, False)
        try:
            align_files("protein.aln.fasta", ["dna.fasta"], out_path="codon.fasta", rec=rec, jobs=1)
        except CoalError as e:
            rec.fail(e.code, str(e), e.seq_id)
            rec.exit_code = 1
        report = rec.to_dict()
        self.assertEqual(report["error"]["code"], "inconsistent")
        self.assertIsNone(report["output"])
        self.assertEqual(report["summary"]["not_checked"], 19)


if __name__ == "__main__":
    unittest.main()
