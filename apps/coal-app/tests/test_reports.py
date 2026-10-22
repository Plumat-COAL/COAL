"""The run report the app builds from the engine: warnings and notices, stopped runs and usage errors, ambiguous codons, ID pairing, malformed reports, the sample report and the methods sentence."""
import os
import shutil
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from coal_app import core, tools, viewer  # noqa: E402
from helpers import ExampleVariantCase, FX, ROOT, RunsDirCase, _load, _randomise, _read_fasta, _tmp, _write, _write_fasta  # noqa: E402


class TestNomismatchResults(RunsDirCase):
    def test_mito_code1_nomismatch_is_not_all_ok(self):
        p, d = core.example_paths("mito")
        r = tools.codon_align(p, d, codontable=1, nomismatch=True)
        self.assertTrue(r["completed"])
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
        self.assertNotEqual(core.banner_head(r["completed"], s), core.banner_head(True, {"ok": 24, "total": 24}))

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
                             text=True, timeout=120, stdin=subprocess.DEVNULL).stdout
        self.assertEqual(got, ref)

    def test_table3_nomismatch_not_nothing_went_wrong(self):
        p, d = os.path.join(FX, "table3", "pep.fasta"), os.path.join(FX, "table3", "nuc.fasta")
        r = core.run_coal(p, d, codontable=1, nomismatch=True, out_dir=self.tmp())
        self.assertNotIn("Nothing went wrong", core.explain(r.report))
        self.assertGreater(r.report["summary"]["columns_removed"], 0)

    def test_long_mismatch_one_warning(self):
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


class TestEngineNotices(RunsDirCase):
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


class TestNoAnchorNoticeWithAmbiguousCount(unittest.TestCase):
    def test_new_notice_with_ambiguous_count(self):
        line = ("coal: warning: g: no 10-residue segment of the protein was found in the DNA; codons were assigned "
                "by position (1 of 12 checkable codons do not match; 30 contain N or other ambiguous bases)")
        sid, msg = core._parse_warning(line)
        self.assertEqual((sid, msg["code"]), ("g", "no_anchor"))
        self.assertIn("1 of 12 checkable codons do not match; 30 contain N", msg["text"])


class TestStoppedRuns(RunsDirCase):
    def report(self, code, stderr):
        p, d = core.example_paths("clean")
        return core.build_fallback_report(p, d, os.path.join(self.tmp(), "out.fa"), "fasta", 1, False, False, False,
                                          code, stderr, 0.1)

    def test_killed_or_interrupted_is_not_nothing_went_wrong(self):
        for code, err in ((143, ""), (-9, ""), (130, "coal: interrupted\n"), (1, "coal: interrupted\n"),
                          (137, "")):
            rep = self.report(code, err)
            text = core.explain(rep)
            self.assertNotIn("Nothing went wrong", text, code)
            self.assertIn("stopped before it finished", text, code)
            self.assertEqual(rep["summary"]["ok"], 0, code)
            self.assertIsNone(rep["output"])

    def test_timeout_is_not_nothing_went_wrong(self):
        rep = self.report(1, "coal: error: COAL did not finish within 5 seconds\n")
        self.assertNotIn("Nothing went wrong", core.explain(rep))

    def test_success_still_says_nothing_went_wrong(self):
        self.assertIn("Nothing went wrong", core.explain(self.report(0, "")))


class TestUsageErrors(RunsDirCase):
    def test_exit_2_usage_is_not_success(self):
        p, d = core.example_paths("clean")
        rep = core.build_fallback_report(p, d, os.path.join(self.tmp(), "o.fa"), "fasta", 1, False, False, False,
                                         2, "usage: coal [-h] ...\ncoal: error: the following arguments are "
                                            "required: nuc_fasta\n", 0.1)
        self.assertEqual(rep["summary"]["ok"], 0)
        self.assertNotIn("Nothing went wrong", core.explain(rep))
        self.assertEqual(core.methods_text(rep), "")
        rr = core.RunResult(2, p, "", rep, "usage: coal", "", "")
        self.assertFalse(rr.ok)

    def test_file_names_starting_with_dash(self):
        d = self.tmp()
        p0, d0 = core.example_paths("clean")
        shutil.copy(p0, os.path.join(d, "-p.fa"))
        shutil.copy(d0, os.path.join(d, "-d.fa"))
        old = os.getcwd()
        os.chdir(d)
        self.addCleanup(os.chdir, old)
        r = tools.codon_align("-p.fa", "-d.fa")
        self.assertTrue(r["completed"], r)
        self.assertEqual(r["summary"]["ok"], 8, r)


def _masked_files(d, frac=0.55, n=100, seed=7):
    """One gene with `frac` of its codons NNN (correctly placed) plus one clean gene."""
    import random
    from coal_app import gcode
    sense = sorted(c for c in (a + b + e for a in "TCAG" for b in "TCAG" for e in "TCAG")
                   if gcode.translate(c, 1) not in ("*", "M"))
    rng = random.Random(seed)
    cod = ["ATG"] + [rng.choice(sense) for _ in range(n - 1)]
    prot = "".join(gcode.translate(c, 1) for c in cod)
    for i in rng.sample(range(1, n), int(frac * n)):
        cod[i] = "NNN"
    p, q = os.path.join(d, "p.fa"), os.path.join(d, "n.fa")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(f">masked\n{prot}\n>clean\nMKVAL{'-' * (n - 5)}\n")
    with open(q, "w", encoding="utf-8") as fh:
        fh.write(f">masked\n{''.join(cod)}\n>clean\nATGAAAGTTGCTCTG\n")
    return p, q, int(frac * n)


class TestAmbiguousCodons(RunsDirCase):
    def test_nnn_codons_are_not_called_a_genetic_code_problem(self):
        p, q, k = _masked_files(self.tmp())
        r = tools.codon_align(p, q)
        self.assertTrue(r["completed"], r.get("next_step"))
        prob = next(x for x in r["problems"] if x["id"] == "masked")
        self.assertIn("unknown or ambiguous bases", prob["explanation"])
        self.assertIn(f"{k} codons", prob["explanation"])
        self.assertNotIn("genetic code", prob["explanation"])
        self.assertNotIn("do not match their amino acid", prob["explanation"])
        self.assertFalse([h for h in r["hints"] if h.get("code") == "genetic_code"], r["hints"])

    def test_nomismatch_says_why_columns_were_removed(self):
        p, q, k = _masked_files(self.tmp())
        r = tools.codon_align(p, q, nomismatch=True)
        self.assertTrue(r["completed"], r.get("next_step"))
        self.assertIn(f"{k} codon columns were removed", r["summary_text"])
        self.assertIn("unknown or ambiguous bases", r["summary_text"])
        self.assertNotIn("because they didn't match", r["summary_text"])
        prob = next(x for x in r["problems"] if x["id"] == "masked")
        self.assertIn("removed", prob["explanation"])
        self.assertNotIn("genetic code", prob["explanation"])

    def test_mixed_definite_and_ambiguous(self):
        msgs = [{"code": "mismatch", "column": 3, "residue": "G", "codon": "NNN"},
                {"code": "mismatch", "column": 7, "residue": "W", "codon": "TGA"}]
        text = core.sequence_sentence({"status": "warning", "messages": msgs})
        self.assertIn("1 codon does not match its amino acid (first at column 7", text)
        self.assertIn("1 codon contains unknown or ambiguous bases", text)
        self.assertIn("genetic code", text)

    def test_real_mismatches_keep_the_genetic_code_wording(self):  # mito example under code 1
        p, d = core.example_paths("mito")
        r = tools.codon_align(p, d, codontable=1)
        self.assertTrue(any(h.get("suggested_codontable") == 2 for h in r["hints"]), r["hints"])
        self.assertTrue(all("genetic code" in x["explanation"] for x in r["problems"]))


class TestEveryWarningCounts(unittest.TestCase):
    def test_unknown_aa_warning_is_a_warning(self):
        t = _tmp(self)
        p = _write(os.path.join(t, "p.fasta"), ">s1\nMBKL\n>s2\nMAKL\n")
        d = _write(os.path.join(t, "d.fasta"), ">s1\nATGGATAAACTG\n>s2\nATGGCTAAACTG\n")
        r = core.run_coal(p, d, out_dir=t, probe_codes=False)
        self.assertIn("B (D or N) is not checked", r.stderr)
        self.assertEqual((r.report["summary"]["ok"], r.report["summary"]["warning"]), (1, 1))
        s1 = next(s for s in r.report["sequences"] if s["id"] == "s1")
        self.assertEqual(s1["status"], "warning")
        sentence = core.sequence_sentence(s1)
        self.assertNotIn("match codon for codon", sentence)
        self.assertIn("B", sentence)

    def test_unrecognised_warning_format_still_counts(self):
        t = _tmp(self)
        p = _write(os.path.join(t, "p.fasta"), ">s1\nMAKL\n")
        d = _write(os.path.join(t, "d.fasta"), ">s1\nATGGCTAAACTG\n")
        rep = core.build_fallback_report(p, d, os.path.join(t, "o.fasta"), "fasta", 1, False, False, False, 0,
                                         "coal: warning: s1 something new happened here\n", 0.1)
        s1 = next(s for s in rep["sequences"] if s["id"] == "s1")
        self.assertEqual(s1["status"], "warning")
        self.assertIn("something new happened here", core.sequence_sentence(s1))


class TestNoPairingSentenceWhenNothingPaired(unittest.TestCase):
    def test_empty_protein_file(self):
        t = _tmp(self)
        p = _write(os.path.join(t, "empty.fasta"), "")
        d = core.example_paths("clean")[1]
        r = core.run_coal(p, d, out_dir=t, probe_codes=False)
        self.assertEqual(r.exit_code, 1)
        text = core.explain(r.report)
        self.assertNotIn("paired", text)
        self.assertIsNone(r.report["summary"]["pairing"])

    def test_dna_given_as_protein(self):
        d = core.example_paths("clean")[1]
        r = core.run_coal(d, d, out_dir=_tmp(self), probe_codes=False)
        self.assertEqual(r.exit_code, 1)
        self.assertNotIn("paired", core.explain(r.report))

    def test_successful_run_still_states_pairing(self):
        p, d = core.example_paths("clean")
        r = core.run_coal(p, d, out_dir=_tmp(self), probe_codes=False)
        self.assertIn("paired by matching their IDs exactly", core.explain(r.report))


class TestIdHints(unittest.TestCase):
    def _files(self, t, pid_suffix, did_suffix):
        p, d = core.example_paths("clean")
        prot = viewer.read_alignment(p)
        dna = dict(viewer.read_alignment(d))
        pp = _write(os.path.join(t, "p.fasta"), "".join(f">{i}{pid_suffix}\n{s}\n" for i, s in prot))
        dd = _write(os.path.join(t, "d.fasta"), "".join(f">{i}{did_suffix}\n{dna[i]}\n" for i, _ in prot))
        return pp, dd

    def test_example_shows_the_real_dna_id(self):
        t = _tmp(self)
        pp, dd = self._files(t, "", ".1")
        h = core.id_hints(pp, dd)
        self.assertEqual(len(h), 1)
        self.assertIn("protein 'gene01' vs DNA 'gene01.1'", h[0]["text"])

    def test_successful_order_pairing_is_not_called_wrong(self):
        t = _tmp(self)
        pp, dd = self._files(t, ".1", "")
        r = core.run_coal(pp, dd, out_dir=t, probe_codes=False)
        self.assertEqual((r.exit_code, r.report["summary"]["ok"]), (0, 8))
        texts = [h["text"] for h in r.hints if h["code"] == "id_mismatch"]
        self.assertEqual(len(texts), 1)
        self.assertNotIn("probably wrong", texts[0])
        self.assertIn("every pair matched", texts[0])

    def test_failed_run_keeps_the_warning(self):
        # the engine pairs messy_ids by their normalised form (order would mispair), so every gene aligns and the
        # report says how they were paired
        p, d = core.example_paths("messy_ids")
        r = core.run_coal(p, d, out_dir=_tmp(self), probe_codes=False)
        self.assertEqual(r.exit_code, 0)
        self.assertEqual(r.report["summary"]["pairing"], "normalised")
        self.assertEqual(r.report["summary"]["ok"], 10)
        self.assertEqual([h["code"] for h in r.hints], ["normalised_ids"])


class TestMessyIdsPairingCause(ExampleVariantCase):
    def test_order_pairing_is_named_as_the_cause(self):
        p, d = core.example_paths("messy_ids")
        # messy_ids aligns (its IDs match after normalising, so COAL pairs by them); the same
        # genes with DNA IDs that match nothing, in the same (shuffled) order, are paired by order and fail
        recs = _read_fasta(d)
        for k, rec in enumerate(recs):
            rec[0] = f"x{k}"
        d = os.path.join(self.d, "unrelated_ids.dna.fa")
        _write_fasta(d, recs)
        r = self.run_app(p, d)
        g = next(x for x in r.report["sequences"] if x["status"] == "failed")
        text = core.sequence_sentence(g)
        self.assertNotIn("truncated", text)
        # the engine's own failure message names the order pairing; no break location is claimed (one would be false
        # here, e.g. "The DNA ends before the protein does (at column 120)" for another gene's DNA)
        self.assertIn("paired by file order with DNA 'x", " ".join(m.get("text", "") for m in g["messages"]))
        hint = next(h for h in r.hints if h["code"] == "id_mismatch")["text"]
        self.assertNotIn("use a COAL version", hint)
        self.assertIn("AT1G79010", hint)


class TestNotCheckedWording(ExampleVariantCase):
    def test_genes_before_and_after_the_error(self):
        p, d = core.example_paths("one_broken_gene")
        r = self.run_app(p, d, keep_going=False)  # a stopped run (the app keeps going by default)
        seqs = {x["id"]: x for x in r.report["sequences"]}
        before = core.sequence_sentence(seqs["OG0001_sp01"])
        after = core.sequence_sentence(seqs["OG0001_sp20"])
        self.assertIn("Not written: COAL stopped at an error in OG0001_sp14", before)
        self.assertIn("this gene itself had no problem", before)
        self.assertNotIn("earlier error", before)
        self.assertIn("Not checked", after)
        self.assertIn("OG0001_sp14", after)
        self.assertEqual(r.report["summary"]["not_checked"], 19)


class TestMalformedReports(unittest.TestCase):
    def _report(self, obj=None, text=None, bom=False):
        import json
        t = _tmp(self)
        p = os.path.join(t, "r.json")
        data = text if text is not None else json.dumps(obj)
        with open(p, "w", encoding="utf-8-sig" if bom else "utf-8") as fh:
            fh.write(data)
        return p

    def test_wrong_types_give_a_clear_message(self):
        from coal_app import tools
        for obj in ({"schema": "coal.report/1", "summary": []},
                    {"schema": "coal.report/1", "output": "x"},
                    {"schema": "coal.report/1", "sequences": [1, 2]},
                    {"schema": "coal.report/1", "sequences": [{"id": "a", "messages": "oops"}]},
                    {"schema": "coal.report/1", "hints": "x"},
                    {"schema": "coal.report/1", "errors": [3]},
                    {"schema": "coal.report/1", "options": []},
                    {"schema": "coal.report/1", "inputs": {}},
                    [1, 2]):
            with self.assertRaises(tools.CoalToolError) as cm:  # raised, not returned
                tools.explain_report(self._report(obj))
            text = str(cm.exception)
            self.assertTrue(text.startswith("Could not read the report:"), (obj, text))
            self.assertIn("not a valid COAL report", text)

    def test_bom_is_accepted(self):
        import json
        from coal_app import tools
        with open(core.SAMPLE_REPORT, encoding="utf-8") as fh:
            obj = json.load(fh)
        text = tools.explain_report(self._report(obj, bom=True))
        self.assertFalse(text.startswith("Could not read"), text)

    def test_not_json(self):
        from coal_app import tools
        with self.assertRaisesRegex(tools.CoalToolError, "Could not read the report"):
            tools.explain_report(self._report(text="{nope"))


class TestSampleReportIsARealRun(unittest.TestCase):
    def test_sample_is_a_real_run_labelled_demo(self):
        rep = core.load_sample_report()
        self.assertEqual(rep["coal_version"], "demo")
        self.assertTrue(rep["demo"])
        self.assertIn("DEMO SAMPLE", rep["source"])
        self.assertEqual(rep["options"]["codontable"], 1)
        self.assertEqual([h["suggested_codontable"] for h in rep["hints"]], [2])
        # every mismatch is one that code 2 explains (TGA = W, ATA = M)
        from coal_app import gcode
        for s in rep["sequences"]:
            for m in s["messages"]:
                if m["code"] == "mismatch":
                    self.assertEqual(gcode.translate(m["codon"], 2), m["residue"], m)
        self.assertIn("Demo sample", core.sample_banner(rep))
        self.assertIn("not your data", core.sample_banner(rep))

    def test_sample_matches_a_fresh_run(self):
        mk = _load("make_sample_report_r", os.path.join(ROOT, "demo", "make_sample_report.py"))
        fresh, stored = mk.build(), core.load_sample_report()
        # every field, summary_text included
        self.assertEqual(fresh, stored)


class TestMethodsTextCounts(ExampleVariantCase):
    def test_suspect_gene_is_written_not_excluded(self):
        r = self.run_app(*self.variant(dna=_randomise(2)))
        self.assertEqual(r.report["summary"]["failed"], 1)
        with open(r.output_path) as fh:
            self.assertEqual(fh.read().count(">"), 8)  # the suspect gene IS in the output
        m = core.methods_text(r.report)
        self.assertNotIn("excluded", m)
        self.assertIn("1 was written to the output but flagged as probably wrong", m)
        self.assertIn("remove it before analysis", m)

    def test_non_standard_residue_is_not_a_mismatch(self):
        def ubz(recs):
            s = recs[0][1]
            i = s.index("W")  # gene01 has a W: replace by U (selenocysteine), B and Z further on
            s = s[:i] + "U" + s[i + 1:]
            recs[0][1] = s
            return recs
        r = self.run_app(*self.variant(prot=ubz))
        self.assertEqual(r.report["summary"]["warning"], 1, r.stderr)
        m = core.methods_text(r.report)
        self.assertNotIn("did not match", m)
        self.assertIn("non-standard amino acids", m)

    def test_pairing_only_warning_is_not_a_mismatch(self):
        def ren(recs):
            recs[2][0], recs[4][0] = "geneX3", "geneX5"
            return recs
        r = self.run_app(*self.variant(dna=ren))
        m = core.methods_text(r.report)
        self.assertNotIn("did not match", m)
        self.assertIn("2 were paired with their DNA only by their order in the files", m)

    def test_real_mismatch_still_counted(self):
        def mm(recs):
            s = recs[0][1]
            recs[0][1] = s[:9] + ("TTT" if s[9:12] != "TTT" else "GGG") + s[12:]  # 4th codon changed
            return recs
        r = self.run_app(*self.variant(dna=mm))
        m = core.methods_text(r.report)
        self.assertIn("1 contained codons that did not match the aligned amino acid", m)

    def test_proper_names_keep_their_capitals(self):
        rep = core.load_sample_report()
        for n, frag in ((22, "Scenedesmus obliquus mitochondrial genetic code"),
                        (23, "Thraustochytrium mitochondrial genetic code"),
                        (15, "Blepharisma nuclear genetic code"),
                        (6, "ciliate, dasycladacean and Hexamita nuclear genetic code"),
                        (4, "Mycoplasma/Spiroplasma"),
                        (1, "standard genetic code (NCBI translation table 1)")):
            r = dict(rep, options=dict(rep["options"], codontable=n))
            self.assertIn(frag, core.methods_text(r))


if __name__ == "__main__":
    unittest.main()
