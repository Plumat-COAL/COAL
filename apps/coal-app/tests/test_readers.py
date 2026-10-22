"""The app's own FASTA/CLUSTAL readers give what the engine reads (BOM, CR/CRLF, '#' lines), the whole-file reader matches the line reader, and pasted paths are cleaned up."""
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from coal_app import core, seqformat, tools  # noqa: E402
from helpers import RunsDirCase  # noqa: E402


BOM = "\ufeff"


def _variants(text: str, dna: bool):
    """Edge-case spellings of one file that the engine accepts (name -> text)."""
    v = {"plain": text, "bom": BOM + text, "crlf": text.replace("\n", "\r\n"), "cr": text.replace("\n", "\r"),
         "blank_lines": "\n\n" + text.replace("\n>", "\n\n>")}
    if dna:
        v["hash_preamble"] = "# exported by tool\n# second line\n" + text
        v["bom_crlf_hash"] = BOM + ("# exported by tool\n\n" + text).replace("\n", "\r\n")
        v["hash_inside"] = text.replace("\n>", "\n# note\n>")
    return v


class TestReadersMirrorEngine(RunsDirCase):
    def write(self, d, name, text):
        p = os.path.join(d, name)
        with open(p, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        return p

    def example(self, name="clean"):
        p, d = core.example_paths(name)
        with open(p, encoding="utf-8") as fh:
            pt = fh.read()
        with open(d, encoding="utf-8") as fh:
            dt = fh.read()
        return pt, dt

    def test_app_readers_give_the_engine_ids_and_sequences(self):
        from coal import seqio
        from coal_app import viewer
        pt, dt = self.example()
        d = self.tmp()
        clustal = "CLUSTAL W\n\ns1  MKV-AL\ns2  MKVSAL\n        ** **\n\ns1  WW\ns2  WW\n"
        cases = [(k, t, False) for k, t in _variants(pt, False).items()]
        cases += [("clustal_" + k, t, False) for k, t in _variants(clustal, False).items() if k != "blank_lines"]
        cases += [(k, t, True) for k, t in _variants(dt, True).items()]
        for k, text, dna in cases:
            with self.subTest(k, dna=dna):
                path = self.write(d, f"{k}_{int(dna)}.fa", text)
                if dna:
                    idx = seqio.DnaIndex([path])
                    want = [(i, idx.get(n)) for n, i in enumerate(idx.ids)]
                    idx.close()
                else:
                    want = seqio.read_alignment(path)
                self.assertTrue(want)
                self.assertEqual(core.read_ids(path), [i for i, _ in want])
                self.assertEqual(viewer.read_alignment(path, dna=dna), want)

    def run_pair(self, ptext, dtext, name="clean", **kw):
        d = self.tmp()
        p = self.write(d, "p.fa", ptext)
        n = self.write(d, "d.fa", dtext)
        return tools.codon_align(p, n, **kw)

    def test_bom_dna_with_mito_protein_runs(self):
        pt, dt = self.example("mito")
        r = self.run_pair(pt, BOM + dt, codontable=1)
        self.assertTrue(r["completed"], r.get("next_step"))
        self.assertTrue(any(h.get("suggested_codontable") == 2 for h in r["hints"]), r["hints"])
        self.assertEqual(r["summary"]["total"], 24, r["summary"])

    def test_hash_preamble_dna_with_mito_protein_runs(self):
        pt, dt = self.example("mito")
        r = self.run_pair(pt, "# exported by tool\n" + dt, codontable=1)
        self.assertTrue(r["completed"], r.get("next_step"))
        self.assertTrue(any(h.get("suggested_codontable") == 2 for h in r["hints"]), r["hints"])
        self.assertTrue(r["problems"])  # per-gene warnings, not only run-level tips

    def test_clean_example_variants_pair_by_id(self):
        pt, dt = self.example()
        for k, (ptext, dtext) in {"bom_dna": (pt, BOM + dt), "hash_dna": (pt, "# x\n" + dt),
                                  "bom_both": (BOM + pt, BOM + dt),
                                  "bom_crlf_both": (BOM + pt.replace("\n", "\r\n"),
                                                    BOM + ("# x\n" + dt).replace("\n", "\r\n"))}.items():
            with self.subTest(k):
                r = self.run_pair(ptext, dtext)
                self.assertTrue(r["completed"], r.get("next_step"))
                s = r["summary"]
                self.assertEqual((s["total"], s["ok"], s["pairing"]), (8, 8, "id"), s)
                self.assertFalse([h for h in r["hints"] if h.get("code") == "id_mismatch"], r["hints"])
                self.assertNotIn("Of 0 sequences", r["methods_text"])


class TestFastFastaReader(unittest.TestCase):
    """The whole-file FASTA reader gives what the line-by-line reader gives (the engine's rules, seqformat)."""
    CASES = {
        "plain": ">a x\nMKV\n>b\nMK-\n",
        "wrapped_crlf": ">a\r\nMK\r\nV-\r\n>b desc\r\nMKV\r\n",
        "bom_comments": "﻿# made by me\n\n>a\nMK V\n#note\nL\n>b\n\nM\n",
        "lone_cr": ">a\rMKV\r>b\rM\r",
        "empty_header": ">\nMK\n>b\nK",
        "dna": ">a\nATG NNN\n>b\natg-ccc\n",
    }

    def test_same_records_as_the_line_reader(self):
        tmp = tempfile.mkdtemp(prefix="coal_test_fasta_")
        self.addCleanup(shutil.rmtree, tmp, True)
        old = seqformat.FAST_READ_MAX
        self.addCleanup(setattr, seqformat, "FAST_READ_MAX", old)
        for name, text in self.CASES.items():
            path = os.path.join(tmp, name + ".fa")
            with open(path, "w", encoding="utf-8", newline="") as fh:
                fh.write(text)
            for dna in (False, True):
                for want in (None, {"b"}):
                    seqformat.FAST_READ_MAX = 1 << 30
                    fast = seqformat.read_records(path, dna=dna, want=want)
                    seqformat.FAST_READ_MAX = -1
                    slow = seqformat.read_records(path, dna=dna, want=want)
                    self.assertEqual(fast, slow, (name, dna, want))
        bad = os.path.join(tmp, "bad.fa")
        with open(bad, "w") as fh:
            fh.write("MKV\n>a\nMK\n")
        seqformat.FAST_READ_MAX = 1 << 30
        with self.assertRaises(ValueError):
            seqformat.read_records(bad)


class TestPastedPaths(unittest.TestCase):
    def setUp(self):
        d = tempfile.mkdtemp(prefix="coal test r2 ")
        self.addCleanup(shutil.rmtree, d, True)
        self.f = os.path.join(d, "my protein.fa")
        with open(self.f, "w", encoding="utf-8") as fh:
            fh.write(">a\nMK\n")

    def test_file_url_with_percent_escapes(self):
        from pathlib import Path
        url = Path(self.f).as_uri()
        self.assertIn("%20", url)
        self.assertEqual(core.resolve_input(None, url)[0], os.path.abspath(self.f))
        self.assertEqual(core.resolve_input(None, "  '" + url + "'  ")[0], os.path.abspath(self.f))

    def test_invisible_marks_and_quotes(self):
        for junk in ("\u202a", "\u200e", "\u200f", "\ufeff", " \u202a"):
            got = core.resolve_input(None, junk + '"' + self.f + '"' + junk + " \n")[0]
            self.assertEqual(got, os.path.abspath(self.f), repr(junk))

    def test_bad_path_still_rejected(self):
        with self.assertRaises(ValueError):
            core.resolve_input(None, "file:///no/such%20file.fa")


if __name__ == "__main__":
    unittest.main()
