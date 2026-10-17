"""ENG-UX: every error path says what is wrong, which file / gene / position, and what to do; coal-cpp says it
byte for byte the same (stderr and exit code).

CHANGED(ENG-UX-1..ENG-UX-17): new. See overnight/integration/CHANGES.md. Python runs in-process; the C++ half is
skipped if coal-cpp is not built (COAL_CPP_BIN, default ./build/coal-cpp).
"""
from __future__ import annotations

import errno
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal.cli import main as coal_main  # noqa: E402
from coal.seqio import open_error_message  # noqa: E402
from test_differential import cpp_binary  # noqa: E402

EX = ROOT / "examples"
P = ">a\nMKV-L\n>b\nMKVAL\n"
D = ">a\nATGAAAGTTCTG\n>b\nATGAAAGTTGCTCTG\n"
LONG_DNA = ">a\nATGAAAGTTCTGATGAAAGTTCTG\n>b\nATGAAAGTTGCTATGAAAGTTGCT\n"  # equal lengths: read as an "alignment"
LONG_PROT = ">a\nMKVLMKVLEEFFIIPPQQLLEEFF\n>b\nMKVAMKVAEEFFIIPPQQLLEEFF\n"
BROKEN = ">a\nMKVAL\n>b\nMKVAL\n", ">a\nATGAAAGTTGCTCTG\n>b\nCCCCCC\n"

# name -> (files, arguments, fragments Python's stderr must contain)
CASES = {
    "missing": ({"d.fa": D}, ["nope.fa", "d.fa"], ["coal: error: cannot open 'nope.fa': no such file\n"]),
    "folder": ({"d.fa": D, "sub/x": ""}, ["sub", "d.fa"], ["coal: error: cannot open 'sub': it is a folder, not a file"]),
    "swapped": ({"p.fa": P, "d.fa": D}, ["d.fa", "p.fa"],
                ["coal: error: d.fa looks like DNA, not a protein alignment: give the protein alignment first and the "
                 "DNA second (coal PROTEIN_ALIGNMENT DNA_FASTA)"]),
    "unaligned": ({"p.fa": ">a\nMKVL\n>b\nMKVAL\n", "d.fa": D}, ["p.fa", "d.fa"],
                  ["protein alignment rows have different lengths (4 and 5, at b): is p.fa an aligned protein file, "
                   "e.g. from MAFFT or MUSCLE?"]),
    "swapped_equal_rows": ({"d.fa": LONG_DNA}, ["d.fa", "d.fa"],
                           ["a: protein and DNA are inconsistent (no match found); its protein row looks like DNA: "
                            "are the two files the other way round?"]),
    "protein_as_dna": ({"p.fa": P, "q.fa": LONG_PROT}, ["p.fa", "q.fa"],
                       ["its DNA looks like protein: the second file must hold the coding DNA"]),
    "more_dna": ({"p.fa": P, "d.fa": D + ">c\nATG\n"}, ["p.fa", "d.fa"],
                 ["number of sequences differs (protein: 2, DNA: 3): each protein needs exactly one DNA record; "
                  "remove the 1 extra DNA record (e.g. 'c', which is not in the protein alignment) or check you used "
                  "the right files"]),
    "more_protein": ({"p.fa": P + ">c\nMKVAL\n", "d.fa": D}, ["p.fa", "d.fa"],
                     ["1 protein has no DNA record (e.g. 'c'): add the missing DNA"]),
    "empty_protein": ({"p.fa": "", "d.fa": D}, ["p.fa", "d.fa"], ["coal: error: p.fa: no sequences found"]),
    "empty_dna": ({"p.fa": P, "d.fa": "\n"}, ["p.fa", "d.fa"],
                  ["coal: error: d.fa: no sequences found: the file is empty or has no '>' header lines"]),
    "not_fasta": ({"p.fa": P, "d.fa": "LOCUS x\nORIGIN\n"}, ["p.fa", "d.fa"],
                  ["d.fa:1: sequence data before the first '>' header: not a FASTA file?"]),
    "binary": ({"p.fa": b"\x89PNG\r\n\x1a\n\xff\x00", "d.fa": D}, ["p.fa", "d.fa"],
               ["p.fa:1: non-ASCII character: sequence files must be plain text"]),
    "empty_header": ({"p.fa": P, "d.fa": ">\nATG\n"}, ["p.fa", "d.fa"],
                     ["d.fa:1: empty FASTA header: a '>' line without a sequence ID"]),
    "duplicate": ({"p.fa": ">a\nMKVAL\n>a\nMKVAL\n", "d.fa": D}, ["p.fa", "d.fa"],
                  ["p.fa: duplicate sequence ID 'a': every sequence needs its own ID"]),
    "clustal_line": ({"p.fa": "CLUSTAL W\n\nfoo\n", "d.fa": D}, ["p.fa", "d.fa"],
                     ["p.fa:3: malformed CLUSTAL line 'foo': each line of a block must hold a sequence name"]),
    "invalid_char": ({"p.fa": ">a\nMK@AL\n>b\nMKVAL\n", "d.fa": D}, ["p.fa", "d.fa"],
                     ["['@'] (first at column 3): a protein alignment may hold only amino-acid letters"]),
    "all_gaps": ({"p.fa": ">a\n-----\n>b\nMKVAL\n", "d.fa": D}, ["p.fa", "d.fa"],
                 ["a: protein sequence is empty (its alignment row has only gaps): remove it from both files"]),
    "gene_stops_run": ({"p.fa": BROKEN[0], "d.fa": BROKEN[1]}, ["p.fa", "d.fa"],
                       ["coal: summary: The run stopped at b, so no output was written. Add --keep-going to write the "
                        "other sequences and leave this one out.\n"]),
    "none_converted": ({"p.fa": BROKEN[0], "d.fa": ">a\nCCCCCC\n>b\nCCCCCC\n"}, ["p.fa", "d.fa", "--keep-going"],
                       ["none of the 2 sequences could be converted: check that the DNA file holds the coding "
                        "sequences of these proteins"]),
    "out_folder": ({"p.fa": P, "d.fa": D, "od/x": ""}, ["p.fa", "d.fa", "-o", "od"],
                   ["coal: error: cannot write od: it is a folder, not a file (give a file name with -o)\n"]),
    "out_is_input": ({"p.fa": P, "d.fa": D}, ["p.fa", "d.fa", "-o", "d.fa"],
                     ["cannot write d.fa: it is one of the input files; choose another name for the output"]),
    "out_missing_folder": ({"p.fa": P, "d.fa": D}, ["p.fa", "d.fa", "-o", "nofolder/x.fa"],
                           ["nofolder does not exist\n"]),
    "report_folder": ({"p.fa": P, "d.fa": D, "rd/x": ""}, ["p.fa", "d.fa", "--report", "rd"],
                      ["coal: error: could not write the report rd: it is a folder, not a file\n"]),
    "resume_notice": ({"p.fa": P, "d.fa": D}, ["p.fa", "d.fa", "-o", "o.fa", "--resume"],
                      ["coal: notice: no checkpoint to resume for o.fa; starting a normal run\n"]),
}
# ENG-UX-18 / ENG-UX-19: real-data findings (REAL-1 review). A camelid VH gene (HyPhy camelid set, LGL237335) whose
# partial codons sit at alignment columns 21 and 24 was "Frameshift at column 94".
FS_PEP = ">s\n---LQQSGGGLVQAGGSLRL-CV-SGIIASDSAMAWYRQAPGKGREYVA-HITAGGLSNYEDFVNGRFTISRDNARNTMYLQMNSLKPEDTAVYYC\n"
FS_NUC = (">s\nCTGCAGCAGTCTGGGGGAGGCTTGGTGCAGGCTGGGGGGTCTCTGAGACTCTCTGTGTAGCTCTGGAATAATCGCCAGCGACAGTGCCATGGCCTGGTACCGC"
          "CAGGCTCCAGGGAAGGGGCGCGAGTATGTCGCACATATCACAGCTGGTGGTCTTTCAAATTATGAAGACTTCGTGAACGGCCGATTCACCATCTCCAGAGACAAC"
          "GCCAGGAACACCATGTATCTGCAAATGAACAGCCTGAAACCTGAGGACACGGCCGTCTATTATTGT\n")


def _dup_case():
    """A gene whose first 30 residues (and codons) recur later, its DNA missing the first 2 nt: was "The DNA ends
    before the protein does (at column 181)"."""
    import random
    from coal.tables import get_code
    rng, code = random.Random(5), get_code(1)
    aas = sorted(a for a in code.codons if a not in "*M")
    head = ["M"] + [rng.choice(aas) for _ in range(29)]
    prot = head + [rng.choice(aas) for _ in range(120)] + head[1:] + [rng.choice(aas) for _ in range(150)]
    codons = [rng.choice(sorted(code.codons[a])) for a in prot]
    for i in range(1, 30):
        codons[149 + i] = codons[i]
    return "".join(prot), "".join(codons)[2:]


DUP_ROW, DUP_DNA = _dup_case()
CASES["order_failure"] = ({"p.fa": ">NP_1.1\nMKVLDEWAGH\n>NP_2.1\nMPPQRSTV--\n",
                           "d.fa": ">NM_2.3 gene two mRNA\nATGCCGCCGCAGCGTTCTACCGTA\n>NM_1.3 gene one mRNA\n"
                                   "ATGAAAGTTCTGGACGAATGGGCAGGTCAT\n"}, ["p.fa", "d.fa"],
                          ["coal: error: NP_1.1: protein and DNA are inconsistent (no match found) (paired by file "
                           "order with DNA 'NM_2.3', because no protein ID matches a DNA ID)\n"])
CASES["frameshift_near"] = ({"p.fa": FS_PEP, "d.fa": FS_NUC}, ["p.fa", "d.fa"],
                            ["coal: summary: s: Frameshift near column 22: there are four extra bases in the DNA"])
CASES["no_false_location"] = ({"p.fa": f">g\n{DUP_ROW}\n", "d.fa": f">g\n{DUP_DNA}\n"}, ["p.fa", "d.fa"],
                              ["coal: summary: The run stopped at g, so no output was written.\n"])
if sys.platform.startswith("linux") and os.path.isdir("/proc"):
    CASES["out_proc"] = ({"p.fa": P, "d.fa": D}, ["p.fa", "d.fa", "-o", "/proc/coal-x.fa"],
                         ["coal: error: cannot write /proc/coal-x.fa: No such file or directory\n"])  # folder exists


def make(d: Path, files: dict) -> None:
    for name, content in files.items():
        p = d / name
        p.parent.mkdir(parents=True, exist_ok=True)
        if name.endswith("/x"):
            continue  # just the folder
        if isinstance(content, bytes):
            p.write_bytes(content)
        else:
            p.write_text(content, newline="\n")


def run_py(d: Path, args) -> tuple[int, str, str]:
    out, err = io.StringIO(), io.StringIO()
    cwd = os.getcwd()
    os.chdir(d)
    try:
        with redirect_stdout(out), redirect_stderr(err):
            try:
                rc = coal_main(list(args))
            except SystemExit as e:
                rc = e.code
    finally:
        os.chdir(cwd)
    return rc, out.getvalue(), err.getvalue()


class TestMessages(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.d = Path(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_each_message(self):
        for name, (files, args, frags) in CASES.items():
            with self.subTest(case=name):
                d = self.d / name
                d.mkdir()
                make(d, files)
                rc, _, err = run_py(d, args)
                self.assertNotEqual(rc, 0, err) if name != "resume_notice" else self.assertEqual(rc, 0, err)
                for f in frags:
                    self.assertIn(f, err)
                # ENG-UX-5: an error is said once (no "The run failed: ..." / "[Errno" repeats)
                self.assertNotIn("The run failed", err)
                self.assertNotIn("[Errno", err)

    def test_run_level_error_is_one_line(self):  # ENG-UX-5
        make(self.d, {"p.fa": P, "d.fa": D + ">c\nATG\n"})
        rc, out, err = run_py(self.d, ["p.fa", "d.fa"])
        self.assertEqual((rc, out, err.count("\n")), (1, "", 1), err)

    def test_open_error_reasons(self):  # ENG-UX-4 (root can read anything, so the reasons are tested directly)
        self.assertEqual(open_error_message("x.fa", PermissionError(errno.EACCES, "Permission denied", "x.fa")),
                         "cannot open 'x.fa': permission denied")
        self.assertEqual(open_error_message("x.fa", OSError(errno.ENAMETOOLONG, "File name too long", "x.fa")),
                         "cannot open 'x.fa': File name too long")

    def test_normalised_pairing_example_reads_the_right_way(self):  # ENG-UX-1
        rc, _, err = run_py(ROOT, [str(EX / "messy_ids" / "protein.aln.fasta"), str(EX / "messy_ids" / "dna.fasta")])
        self.assertEqual(rc, 0)
        self.assertIn("(e.g. by file order, protein 'AT1G79010.1' would get DNA 'AT2G68010'; COAL paired it with "
                      "'AT1G79010')", err)
        self.assertNotIn("gets DNA", err)

    def test_broken_gene_says_where_and_what_to_do(self):  # ENG-UX-2, ENG-UX-5
        rc, _, err = run_py(ROOT, [str(EX / "one_broken_gene" / "protein.aln.fasta"),
                                   str(EX / "one_broken_gene" / "dna.fasta")])
        self.assertEqual(rc, 1)
        self.assertEqual(err, "coal: error: OG0001_sp14: protein and DNA are inconsistent (no match found)\n"
                              "coal: summary: OG0001_sp14: Frameshift at column 102: one base is missing in the DNA, so "
                              "from here on the codons are read out of step.\n"
                              "coal: summary: The run stopped at OG0001_sp14: 6 sequences were not checked and no output "
                              "was written. Add --keep-going to write the other sequences and leave this one out.\n")

    def test_genetic_code_hint_names_the_option_report_keeps_app_text(self):  # ENG-UX-7, ENG-UX-9
        out = self.d / "o.fa"
        rep = self.d / "r.json"
        rc, _, err = run_py(ROOT, [str(EX / "mito" / "protein.aln.fasta"), str(EX / "mito" / "dna.fasta"),
                                   "-o", str(out), "--report", str(rep)])
        self.assertEqual(rc, 0)
        self.assertTrue(err.rstrip("\n").endswith("Re-run with -codontable 2."), err[-300:])
        r = json.loads(rep.read_text())
        hint = next(h for h in r["hints"] if h["code"] == "genetic_code")
        self.assertEqual(hint["suggested_codontable"], 2)
        self.assertTrue(hint["text"].endswith("Re-run with code 2?"))  # the app's contract
        self.assertTrue(r["methods_text"].startswith("Fix the warnings first: only 0 of 24 sequences matched their "
                                                     "protein exactly"), r["methods_text"][:200])
        self.assertIn("Methods text for this run as it stands: Codon alignments were built with COAL",
                      r["methods_text"])
        # a clean run has no such note
        rc, _, _ = run_py(ROOT, [str(EX / "clean" / "protein.aln.fasta"), str(EX / "clean" / "dna.fasta"),
                                 "-o", str(out), "--report", str(rep)])
        self.assertTrue(json.loads(rep.read_text())["methods_text"].startswith("Codon alignments were built"))

    def test_standard_code_is_suggested_when_it_fits(self):  # PRE-4
        # one_broken_gene run with code 2: codes 1, 4, 6, 10, 11, 15 and 16 all fit the 19 healthy genes; the
        # hint said "code 4 ... differs least from code 2"; now code 1, the standard code
        rep = self.d / "r.json"
        args = [str(EX / "one_broken_gene" / "protein.aln.fasta"), str(EX / "one_broken_gene" / "dna.fasta"),
                "-codontable", "2", "--keep-going", "-o", str(self.d / "o.fa"), "--report", str(rep)]
        rc, _, err = run_py(ROOT, args)
        self.assertIn("Under genetic code 1 (Standard (nuclear genes of most organisms)), 19 of 20 sequences match "
                      "their protein; under code 2, none do. Codes 4, 10, 11, 15, 16 and 6 fit as well, but code 1 "
                      "is the standard code, used by most nuclear genes. Re-run with -codontable 1.\n", err)
        self.assertNotIn("differs least", err)
        hint = next(h for h in json.loads(rep.read_text())["hints"] if h["code"] == "genetic_code")
        self.assertEqual(hint["suggested_codontable"], 1)
        self.assertTrue(hint["text"].endswith("but code 1 is the standard code, used by most nuclear genes. "
                                              "Re-run with code 1?"), hint["text"])
        # the mito case is unchanged: run with 1, code 2 suggested, 5 and 13 named
        args = [str(EX / "mito" / "protein.aln.fasta"), str(EX / "mito" / "dna.fasta"), "-o", str(self.d / "o.fa"),
                "--report", str(rep)]
        rc, _, err = run_py(ROOT, args)
        self.assertIn("Under genetic code 2 (Vertebrate mitochondrial), all 24 sequences match their protein; under "
                      "code 1, none do. Codes 5 and 13 fit these sequences too: choose the code of the organism and "
                      "organelle the genes come from. Re-run with -codontable 2.\n", err)

    def test_break_locations_are_not_wrong(self):  # ENG-UX-19
        from coal.report import break_event
        self.assertIsNone(break_event(DUP_ROW, DUP_DNA, 1))  # no evidence for a column: no column
        ev = break_event(FS_PEP.split("\n")[1], FS_NUC.split("\n")[1], 1)
        self.assertEqual((ev["code"], ev["column"], ev["shift"]), ("frameshift", 22, 4))
        self.assertTrue(ev["text"].startswith("Frameshift near column 22"))

    def test_order_pairing_named_only_on_failure(self):  # ENG-UX-18 (ENG-d: silent when every gene fits)
        make(self.d, {"p.fa": ">x\nMKV-L\n>y\nMKVAL\n", "d.fa": D})
        rc, _, err = run_py(self.d, ["p.fa", "d.fa"])
        self.assertEqual((rc, err), (0, ""))
        make(self.d, {"p.fa": ">x\nMKV-L\n>y\nMKWAL\n>z\nMKVAL\n", "d.fa": D + ">z\nCCCCCC\n"})
        rc, _, err = run_py(self.d, ["p.fa", "d.fa", "--keep-going"])
        self.assertIn("coal: warning: skipped y: protein and DNA are inconsistent (no match found) (paired by file "
                      "order with DNA 'b', because only 1 of 3 protein IDs are DNA IDs)", err)
        self.assertNotIn("skipped z: protein and DNA are inconsistent (no match found) (paired", err)  # same ID

    def test_clean_example_readme_counts(self):  # ENG-UX-8
        text = (EX / "clean" / "README.md").read_text()
        self.assertNotIn("150 codons", text)
        self.assertIn("140-144 codons each, aligned in 150 columns", text)


@unittest.skipUnless(cpp_binary(), "coal-cpp not built (COAL_CPP_BIN)")
class TestCppSaysTheSame(unittest.TestCase):
    """coal-cpp: the same stderr (byte for byte) and exit code on every case above and the examples."""

    def test_codon_output_when_the_last_sequence_is_skipped(self):  # ENG-UX-20 (TestCppMatchesPython seed 31337)
        pep = '>seq0\nmsr-vet-yhsrneldpyykkalmafwkqyvm-gpwlklhccv\n>seq1\nMSRSVECDY-SRNELD-YYKKD-2GR-KQ-VMPGPTLKWHCCV\n>seq2\nMS-SVETD-HNRNELDPY4PKALTM!WKWYVMPFPTFKLRCCD\n'
        nuc = '>seq0 desc\nATGTCACGCGTTGAGACATACCACTCACGTAATGAGCTCGACCCTTATTATAAGAAGGCCCTGATGGCCTTCTGGAAATAATACGTTATGGGGCCGTGGCTGAAACTGCACTGCTGTGWGTGAA\n>seq1 desc\nCGATACGATGTCGAGATCGGTCGAGTGCGATTACTCTCGAAACGAGTTGGACTACTATAAGAAGGACACGGACGTAAACAAGTAATGCCCGGRCCGACTCTCAAATGGCATTGCTGTGTG\n>seq2 desc\nATGTCTTCGGTTGAGACGGACCATAATCGGAACGAGCTAGATCCTTATTATACCGAAAGCGCTTACGATGACCTGGAAGTGGTATGTTATGCCTTTCCCAACGTTTAAGCTAAGATGCTGTGAC\n'
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            make(d, {"pep.fasta": pep, "nuc.fasta": nuc})
            args = ["pep.fasta", "nuc.fasta", "-output", "codon", "-codontable", "6", "--keep-going"]
            rc, out, err = run_py(d, args)
            r = subprocess.run([str(cpp_binary()), *args], cwd=d, capture_output=True, timeout=60)
            self.assertEqual(rc, 2)
            self.assertEqual((r.returncode, r.stdout.decode()), (rc, out))

    def test_same_messages(self):
        binary = str(cpp_binary())
        with tempfile.TemporaryDirectory() as tmp:
            runs = [(n, f, a) for n, (f, a, _) in CASES.items()]
            for ex in ("messy_ids", "one_broken_gene", "mito"):
                runs.append((ex, {}, [str(EX / ex / "protein.aln.fasta"), str(EX / ex / "dna.fasta"), "-o", "o.fa"]))
            # PRE-4: the standard-code preference (one_broken_gene and clean run with code 2)
            for ex in ("one_broken_gene", "clean"):
                runs.append((ex + "-code2", {}, [str(EX / ex / "protein.aln.fasta"), str(EX / ex / "dna.fasta"),
                                                 "-o", "o.fa", "-codontable", "2", "--keep-going"]))
            for name, files, args in runs:
                with self.subTest(case=name):
                    d = Path(tmp) / name
                    d.mkdir()
                    make(d, files)
                    rc, out, err = run_py(d, args)
                    for f in ("o.fa",):
                        if (d / f).exists():
                            (d / f).unlink()
                    r = subprocess.run([binary, *args], cwd=d, capture_output=True, timeout=60)
                    cerr = r.stderr.decode("utf-8", "replace").replace("\r\n", "\n")
                    self.assertEqual((r.returncode, cerr), (rc, err))


    def test_same_genetic_code_hint_in_the_report(self):  # PRE-4
        binary = str(cpp_binary())
        with tempfile.TemporaryDirectory() as tmp:
            for ex, code in (("one_broken_gene", "2"), ("clean", "2"), ("mito", "1")):
                with self.subTest(example=ex, code=code):
                    d = Path(tmp)
                    args = [str(EX / ex / "protein.aln.fasta"), str(EX / ex / "dna.fasta"), "-o", "o.fa",
                            "-codontable", code, "--keep-going", "--report", "r.json"]
                    run_py(d, args)
                    py = json.loads((d / "r.json").read_text())
                    subprocess.run([binary, *args], cwd=d, capture_output=True, timeout=60)
                    cpp = json.loads((d / "r.json").read_text())
                    self.assertEqual(cpp["hints"], py["hints"])
                    hint = next(h for h in py["hints"] if h["code"] == "genetic_code")
                    self.assertEqual(hint["suggested_codontable"], {"1": 2, "2": 1}[code])

if __name__ == "__main__":
    unittest.main()
