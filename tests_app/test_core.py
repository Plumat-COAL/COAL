import json
import os
import tempfile
import shutil
import unittest

from coal_app import core


def _tmpdir(test):
    d = tempfile.mkdtemp(prefix="coal_test_")
    test.addCleanup(shutil.rmtree, d, True)
    return d


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestSampleReport(unittest.TestCase):
    def setUp(self):
        self.rep = core.load_sample_report()

    def test_schema_fields(self):
        r = self.rep
        self.assertEqual(r["schema"], "coal.report/1")
        for k in ("coal_version", "command", "pal2nal_equivalent", "options", "inputs", "output",
                  "summary", "hints", "sequences", "runtime_seconds", "exit_code"):
            self.assertIn(k, r)
        self.assertEqual(set(r["options"]), {"output_format", "codontable", "nogap", "nomismatch", "keep_going",
                                             "strict", "blockonly"})  # FIN-A6: strict, blockonly
        self.assertEqual({i["role"] for i in r["inputs"]}, {"protein", "dna"})

    def test_counts(self):  # APP-15: the sample is now a real run of examples/mito with code 1
        r = self.rep
        self.assertEqual(len(r["sequences"]), 24)
        st = [s["status"] for s in r["sequences"]]
        self.assertEqual((st.count("ok"), st.count("warning"), st.count("failed")), (0, 24, 0))
        self.assertEqual(r["summary"], {"total": 24, "ok": 0, "warning": 24, "failed": 0, "pairing": "id"})
        self.assertEqual([h["code"] for h in r["hints"]], ["genetic_code"])
        self.assertEqual(r["hints"][0]["suggested_codontable"], 2)

    def test_explain_and_rows(self):
        text = core.explain(self.rep)
        self.assertIn("24 with warnings", text)
        self.assertIn("Homo_COX1", text)
        self.assertIn("genetic code 2", text)
        rows = core.health_rows(self.rep)
        self.assertEqual(len(rows), 24)
        for icon, word, sid, sentence in rows:
            self.assertTrue(icon and word and sid and sentence)  # never colour alone
        html = core.health_table_html(self.rep)
        self.assertIn("<caption>", html)
        self.assertIn("Homo_COX1", html)

    def test_load_report_rejects_other_json(self):
        with tempfile.TemporaryDirectory() as t:
            p = os.path.join(t, "x.json")
            with open(p, "w", encoding="utf-8") as fh:
                json.dump({"hello": 1}, fh)
            with self.assertRaises(ValueError):
                core.load_report(p)


class TestHelpers(unittest.TestCase):
    def test_codes_and_formats(self):
        self.assertEqual(sorted(core.GENETIC_CODES), [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23])
        f = core.detect_features()
        self.assertIn("fasta", f.output_formats)

    def test_read_ids_fasta_and_clustal(self):
        p, d = core.example_paths("messy_ids")
        self.assertTrue(all(i.endswith(".1") for i in core.read_ids(p)))
        with tempfile.TemporaryDirectory() as t:
            c = os.path.join(t, "a.aln")
            with open(c, "w") as fh:
                fh.write("CLUSTAL W\n\nseqA  MK-L\nseqB  MKQL\n     ** *\n\nseqA  AA\nseqB  AA\n")
            self.assertEqual(core.read_ids(c), ["seqA", "seqB"])

    def test_commands(self):
        c = core.coal_command("pep.aln", "nuc.fa", "paml", 2, True, False, "out.paml")
        self.assertEqual(c, "coal pep.aln nuc.fa -output paml -codontable 2 -nogap -o out.paml")
        p = core.pal2nal_command("pep.aln", "nuc.fa", "fasta", 1, False, True, "out.fa")
        self.assertEqual(p, "pal2nal.pl pep.aln nuc.fa -output fasta -nomismatch > out.fa")

    def test_bad_format_rejected(self):
        p, d = core.example_paths("clean")
        with self.assertRaises(ValueError):
            core.run_coal(p, d, output_format="nonsense")


class TestRunsOnExamples(unittest.TestCase):
    def run_example(self, name, **kw):
        p, d = core.example_paths(name)
        t = _tmpdir(self)
        return core.run_coal(p, d, out_dir=t, **kw)

    def test_clean(self):
        r = self.run_example("clean")
        self.assertTrue(r.ok)
        self.assertEqual(r.report["summary"]["ok"], 8)
        self.assertEqual(r.hints, [])
        self.assertTrue(os.path.isfile(r.output_path))
        self.assertTrue(os.path.isfile(r.report_path))
        self.assertIn("Nothing went wrong", core.explain(r.report))

    def test_mito_gets_code2_hint(self):
        r = self.run_example("mito")
        codes = [h.get("suggested_codontable") for h in r.hints if h["code"] == "genetic_code"]
        self.assertEqual(codes, [2])
        r2 = self.run_example("mito", codontable=2)
        self.assertEqual(r2.report["summary"]["ok"], 24)
        self.assertEqual(r2.hints, [])

    def test_messy_ids_hint(self):
        r = self.run_example("messy_ids")
        if core.detect_features().report:  # newer engine: its own report and ID normalisation
            self.assertEqual(r.report["summary"]["pairing"], "normalised")
        else:  # this engine: no normalisation, pairs by order, the app's hint explains why
            self.assertEqual(r.report["summary"]["pairing"], "order")
            self.assertIn("id_mismatch", [h["code"] for h in r.hints])
            self.assertIn("version", r.hints[0]["text"])

    def test_broken_gene_is_named(self):
        r = self.run_example("one_broken_gene")
        failed = [s["id"] for s in r.report["sequences"] if s["status"] == "failed"]
        self.assertEqual(failed, ["OG0001_sp14"])
        self.assertIn("OG0001_sp14", core.explain(r.report))
        if core.detect_features().keep_going:  # engine with --keep-going: the 19 good genes are written
            self.assertEqual(r.exit_code, 2)
            self.assertIsNotNone(r.output_path)
        else:  # this engine stops at the broken gene
            self.assertEqual(r.exit_code, 1)
            self.assertIsNone(r.output_path)



class TestMethodsText(unittest.TestCase):
    def test_from_sample_report(self):
        m = core.methods_text(core.load_sample_report())
        self.assertIn("COAL (demo)", m)  # APP-15
        self.assertIn("standard genetic code (NCBI translation table 1)", m)
        self.assertIn("Of 24 sequences, 0 matched", m)
        self.assertIn("24 contained codons that did not match", m)
        # CHANGED(MRG-10): was "... dna.fasta -o ..."; the app now passes --keep-going to the merged engine (MRG-6/2)
        self.assertIn("Command: coal protein.aln.fasta dna.fasta --keep-going -o protein.codon.fasta", m)

    def test_engine_text_wins(self):
        rep = dict(core.load_sample_report(), methods_text="Engine sentence.")
        self.assertEqual(core.methods_text(rep), "Engine sentence.")

    def test_real_run_uses_basenames_and_options(self):
        p, d = core.example_paths("mito")
        r = core.run_coal(p, d, codontable=2, nogap=True, out_dir=_tmpdir(self))
        m = core.methods_text(r.report)
        self.assertIn("vertebrate mitochondrial genetic code (NCBI translation table 2)", m)
        self.assertIn("gaps or stop codons were removed", m)
        self.assertIn("Command: coal protein.aln.fasta dna.fasta -codontable 2 -nogap", m)
        self.assertNotIn(os.path.dirname(p), m)
        self.assertIn("Methods text (ready to paste):", core.explain(r.report))

    def test_no_text_without_output(self):
        p, d = core.example_paths("one_broken_gene")
        r = core.run_coal(p, d, out_dir=_tmpdir(self))
        if core.detect_features().keep_going:  # output written: the methods text counts the failed gene
            self.assertIn("1 could not be aligned", core.methods_text(r.report))
        else:
            self.assertIsNone(r.output_path)
            self.assertEqual(core.methods_text(r.report), "")


if __name__ == "__main__":
    unittest.main()
