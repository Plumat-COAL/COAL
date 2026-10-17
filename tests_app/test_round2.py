"""Round 2 review fixes (overnight/integration/CHANGES.md, R2-A*, R2-D*, R2-E*)."""
import os
import shutil
import tempfile
import unittest

from coal_app import core, tools

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FX = os.path.join(ROOT, "tests", "fixtures")


class RunsDirCase(unittest.TestCase):
    def setUp(self):
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)

    def tmp(self):
        d = tempfile.mkdtemp(prefix="coal_test_r2_")
        self.addCleanup(shutil.rmtree, d, True)
        return d


class TestR2A1Nomismatch(RunsDirCase):  # R2-A1 / R2-D5 / R2-E1
    def test_mito_code1_nomismatch_is_not_all_ok(self):
        p, d = core.example_paths("mito")
        r = tools.codon_align(p, d, codontable=1, nomismatch=True)
        self.assertTrue(r["ok"])
        s = r["summary"]
        self.assertEqual(s["ok"], 0, s)
        self.assertEqual(s["warning"], 24, s)
        self.assertEqual(s["columns_removed"], 46, s)
        self.assertIn("46 codon columns were removed", r["summary_text"])
        self.assertTrue(any(h.get("suggested_codontable") == 2 for h in r["hints"]), r["hints"])
        self.assertNotIn("24 matched", r["methods_text"])
        self.assertIn("46 codon columns", r["methods_text"])
        for prob in r["problems"]:
            self.assertIn("removed", prob["explanation"])
        text = tools.explain_report(r["report_path"])
        self.assertNotIn("Nothing went wrong", text)
        self.assertNotEqual(core.banner_head(r["ok"], s), core.banner_head(True, {"ok": 24, "total": 24}))

    def test_probe_leaves_no_files(self):
        p, d = core.example_paths("mito")
        out = self.tmp()
        core.run_coal(p, d, codontable=1, nomismatch=True, out_dir=out)
        self.assertEqual(sorted(os.listdir(out)), ["protein.codon.fasta", "protein.report.json"])

    def test_output_unchanged_by_probe(self):
        p, d = core.example_paths("mito")
        a = core.run_coal(p, d, codontable=1, nomismatch=True, out_dir=self.tmp(), probe_codes=False)
        with open(a.output_path, encoding="utf-8") as fh:
            got = fh.read()
        import subprocess, sys
        ref = subprocess.run([sys.executable, "-m", "coal", p, d, "-nomismatch"], capture_output=True,
                             text=True, timeout=120).stdout
        self.assertEqual(got, ref)

    def test_table3_nomismatch_not_nothing_went_wrong(self):
        p, d = os.path.join(FX, "table3", "pep.fasta"), os.path.join(FX, "table3", "nuc.fasta")
        r = core.run_coal(p, d, codontable=1, nomismatch=True, out_dir=self.tmp())
        self.assertNotIn("Nothing went wrong", core.explain(r.report))
        self.assertGreater(r.report["summary"]["columns_removed"], 0)

    def test_long_mismatch_one_warning(self):  # R2-E1 repro
        p, d = os.path.join(FX, "long_mismatch", "pep.fasta"), os.path.join(FX, "long_mismatch", "nuc.fasta")
        r = core.run_coal(p, d, nomismatch=True, out_dir=self.tmp(), probe_codes=False)
        s = r.report["summary"]
        self.assertEqual((s["ok"], s["warning"], s["columns_removed"]), (1, 1, 1))
        bad = [x for x in r.report["sequences"] if x["status"] == "warning"]
        self.assertEqual(bad[0]["id"], "s2")
        self.assertEqual(bad[0]["messages"][0]["column"], 13)

    def test_clean_nomismatch_is_ok(self):
        p, d = core.example_paths("clean")
        r = core.run_coal(p, d, nomismatch=True, out_dir=self.tmp())
        s = r.report["summary"]
        self.assertEqual((s["ok"], s["warning"], s.get("columns_removed")), (8, 0, 0))
        self.assertIn("Nothing went wrong", core.explain(r.report))


class TestR2A2Notices(RunsDirCase):  # R2-A2 / R2-E2
    def write(self, pep, nuc):
        d = self.tmp()
        p, n = os.path.join(d, "p.fa"), os.path.join(d, "n.fa")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(pep)
        with open(n, "w", encoding="utf-8") as fh:
            fh.write(nuc)
        return p, n

    def test_parse_no_anchor_notice_strips_colon(self):
        line = ("coal: warning: s0: no 10-residue segment of the protein was found in the DNA; codons were "
                "assigned by position (7 of 25 do not match)")
        sid, msg = core._parse_warning(line)
        self.assertEqual(sid, "s0")
        self.assertEqual(msg["code"], "no_anchor")

    def test_table3_notice_marks_gene(self):
        p, d = os.path.join(FX, "table3", "pep.fasta"), os.path.join(FX, "table3", "nuc.fasta")
        r = core.run_coal(p, d, codontable=1, out_dir=self.tmp(), probe_codes=False)
        s0 = next(x for x in r.report["sequences"] if x["id"] == "s0")
        self.assertEqual(s0["status"], "warning")
        self.assertTrue(any(m["code"] == "no_anchor" for m in s0["messages"]))
        self.assertIn("assigned by position", core.sequence_sentence(s0))
        self.assertNotIn("Nothing went wrong", core.explain(r.report))

    def test_partial_id_match_is_a_prominent_hint(self):
        p, n = self.write(">geneA\nMKVAL\n>geneB\nMKVAL\n", ">geneB\nATGAAAGTTGCTCTG\n>geneX\nATGAAGGTCGCACTA\n")
        r = core.run_coal(p, n, out_dir=self.tmp(), probe_codes=False)
        self.assertIn("1 of 2 protein IDs", r.stderr)
        codes = [h["code"] for h in r.hints]
        self.assertIn("partial_id_match", codes)
        h = next(h for h in r.hints if h["code"] == "partial_id_match")
        self.assertIn("geneA", h["text"])
        self.assertIn("geneB", h["text"])
        self.assertNotIn("looks right", " ".join(x["text"] for x in r.hints))
        text = core.explain(r.report)
        self.assertNotIn("Nothing went wrong", text)
        self.assertIn("1 of 2", text)
        self.assertIn("problems", core.banner_head(r.ok, r.report["summary"], r.hints))

    def test_unknown_warning_line_is_never_dropped(self):
        rep = core.build_fallback_report(*self.write(">a\nMK\n", ">a\nATGAAA\n"), "out.fa", "fasta", 1, False,
                                         False, False, 0, "coal: warning: something new happened here\n", 0.1)
        self.assertTrue(any("something new happened here" in h["text"] for h in rep["hints"]), rep["hints"])
        self.assertNotIn("Nothing went wrong", core.explain(rep))


if __name__ == "__main__":
    unittest.main()
