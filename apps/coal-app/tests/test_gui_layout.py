"""The GUI's header, summary banner, health table grouping, viewer theming and layout.
Runs without Gradio: the layout is checked against the Gradio stand-in in helpers.py; also long warning groups."""
import os
import re
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from coal_app import core, gui_html, viewer  # noqa: E402
from helpers import _Comp, _stub_gradio  # noqa: E402


def _seq(i, status="ok", **kw):
    return {"id": f"g{i:02d}", "status": status, "messages": [], **kw}


class TestGui1Banner(unittest.TestCase):
    def test_headline_counts_aligned_genes(self):
        rep = {"summary": {"total": 20, "ok": 19, "failed": 1}, "sequences": [_seq(i) for i in range(19)] +
               [_seq(19, "failed")]}
        html = gui_html.run_banner("### ⚠ Done, with some problems (see the table)", rep, "s.", "Files used: x")
        self.assertIn("19 of 20 genes aligned", html)
        self.assertIn("1 failed", html)
        self.assertIn('role="status"', html)
        self.assertEqual(html.count('class="coal-tile '), 20)
        self.assertIn('title="g19: Failed"', html)  # tiles carry the gene and a word, not colour alone

    def test_all_ok_and_stopped(self):
        ok = {"summary": {"total": 3, "ok": 3}, "sequences": [_seq(i) for i in range(3)]}
        self.assertIn("All 3 genes match codon for codon",
                      gui_html.run_banner("### ✔ Done: every sequence matched", ok, "s."))
        stop = {"summary": {"total": 3, "failed": 1, "not_checked": 2}, "sequences": []}
        html = gui_html.run_banner("### ✖ COAL stopped: no alignment was written", stop, "s.")
        self.assertEqual(html.count("COAL stopped: no alignment was written"), 1)  # not twice

    def test_stopped_run_banner_says_why(self):
        # a run that stopped before pairing (an unreadable file): the engine's reason is in the banner, file names only
        rep = {"summary": {"total": 1, "ok": 0}, "sequences": [],
               "errors": [{"code": "error", "text": "/srv/up/abc/garbage.fasta:1: control character: sequence files "
                                                   "must be plain text (is this a binary or word-processor file?)"}]}
        why = core.stop_reason(rep, ["/srv/up/abc/garbage.fasta", "/srv/up/abc/dna.fasta"])
        self.assertTrue(why.startswith("garbage.fasta:1: control character"), why)
        self.assertNotIn("/srv/", why)
        html = gui_html.run_banner("### ✖ COAL stopped: no alignment was written", rep,
                                   f"{core.summary_sentence(rep)} {why}", "Files used: x")
        self.assertIn("control character", html)
        self.assertIn("word-processor file?", html)
        self.assertEqual(core.stop_reason({"summary": {}}), "")

    def test_code_fix_banner_says_when_the_engine_tested_a_sample(self):
        names = {1: "Standard (nuclear)", 2: "Vertebrate mitochondrial (x)"}
        fix = {"current": 2, "code": 1, "checked": 50, "fit": 50, "fit_current": 0, "total": 100000, "others": ""}
        html = gui_html.code_fix_banner(fix, {"summary": {"total": 100000, "warning": 50}}, names)
        self.assertIn("50 of the 50 genes COAL tested (a sample of the 100,000)", html)
        self.assertIn("all 50 of the sample match", html)
        small = dict(fix, total=24, checked=24, fit=24)
        self.assertIn("24 of 24 genes have mismatched codons", gui_html.code_fix_banner(small, {"summary": {"total": 24}}, names))
        rep = {"summary": {"total": 100000}, "options": {"codontable": 2},
               "hints": [{"code": "genetic_code", "suggested_codontable": 1, "checked": 50, "fit": 50, "fit_current": 0,
                          "text": "..."}]}
        self.assertEqual(core.code_fix(rep)["total"], 100000)

    def test_strip_is_capped_and_escaped(self):
        big = {"summary": {"total": 401, "ok": 401}, "sequences": [_seq(i) for i in range(401)]}
        self.assertNotIn("coal-strip", gui_html.run_banner("### ✔ Done", big, "s."))
        odd = {"summary": {"total": 1, "ok": 1}, "sequences": [{"id": '<b>"x', "status": "ok"}]}
        self.assertNotIn('<b>"x', gui_html.run_banner("### ✔ Done", odd, "s."))

    def test_header_has_no_external_resources(self):
        html = gui_html.header_html(core.intro_text())
        self.assertNotRegex(html, r"https?://(?!www\.w3\.org/2000/svg)")
        self.assertIn("<svg", html)
        # the header is the mark and "COAL" only (owner: "just COAL is fine"); the where-the-data-goes note is the
        # page footer. Checked with the hosted note (the local app has none, and "" is in every string)
        intro = core.INTRO_HOSTED
        self.assertNotIn(intro, gui_html.header_html(intro))
        self.assertNotIn("codon alignment", html)
        self.assertIn(intro, gui_html.footer_html(intro))
        self.assertEqual(gui_html.footer_html(""), "")


class TestGui4HealthTable(unittest.TestCase):
    def test_ok_genes_are_one_collapsed_row(self):
        rep = {"sequences": [_seq(i) for i in range(19)] + [
            _seq(19, "failed", messages=[{"code": "frameshift", "column": 5, "text": "Frameshift at column 5."}])]}
        html = core.health_table_html(rep, links=True)
        self.assertEqual(html.count("<tr"), 3)  # header, failed gene, one OK row
        self.assertIn("19 genes match codon for codon", html)
        self.assertIn("<details>", html)
        self.assertIn("<code>g00</code>", html)
        self.assertIn('href="#coalv-g19-c5"', html)

    def test_identical_problem_rows_are_grouped(self):
        nc = [{"code": "not_reached", "failed_id": "g00"}]
        rep = {"sequences": [_seq(0, "failed")] + [_seq(i, "not_checked", messages=nc) for i in range(1, 6)]}
        html = core.health_table_html(rep)
        self.assertIn("5 genes", html)
        self.assertEqual(html.count("before reaching this gene"), 1)
        two = {"sequences": [_seq(i, "not_checked", messages=nc) for i in range(2)]}
        self.assertEqual(core.health_table_html(two).count("before reaching this gene"), 2)  # 2: not grouped

    def test_empty_report_is_no_table(self):
        self.assertEqual(core.health_table_html({}), "")


class TestGui2ViewerTheme(unittest.TestCase):
    def test_dark_variables_for_gradio_and_for_the_page(self):
        self.assertIn(".dark div.coalv {", viewer.CSS)
        self.assertIn("@media (prefers-color-scheme: dark) { body.coalv-page div.coalv", viewer.CSS)
        page = viewer.render_page("<div class='coalv'></div>")
        self.assertIn('<body class="coalv-page">', page)
        self.assertNotRegex(page, r"https?://")

    def test_text_colours_are_pinned(self):
        for sel in ("div.coalv h2", "div.coalv a", "div.coalv td {", "div.coalv th.name"):
            rule = viewer.CSS[viewer.CSS.index(sel):]
            rule = rule[:rule.index("}")]
            self.assertRegex(rule, r"(?<![-\w])color:var\(--\w+\) !important", sel)

    def test_dark_tints_keep_contrast(self):
        def lum(h):
            r, g, b = (int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))
            f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4  # noqa: E731
            return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)

        def ratio(a, b):
            x, y = sorted((lum(a), lum(b)))
            return (y + 0.05) / (x + 0.05)
        for pal in (viewer._LIGHT, viewer._DARK):
            for k in ("hy", "ar", "po", "ne", "pl", "sp", "st", "hl", "bg"):
                self.assertGreaterEqual(ratio(pal["fg"], pal[k]), 4.5, k)
            for k in ("muted", "gap", "link"):
                self.assertGreaterEqual(ratio(pal[k], pal["bg"]), 4.5, k)

    def test_residue_classes_stay_apart_for_colour_blind_readers(self):
        # verification 2026-10-24, V-17: the pastel tints of equal lightness merged under red-green deficiency and were
        # one grey in greyscale. The palette now varies lightness too: every pair of the seven classes keeps a CIE76
        # distance under the Machado et al. (2009) simulations of protanopia, deuteranopia and tritanopia, and a
        # lightness step in greyscale (the letters, outlines and legend still carry the meaning on their own).
        import itertools
        import math
        sims = {
            "protan": [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
            "deutan": [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
            "tritan": [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]],
        }

        def lin(h):
            out = []
            for i in (1, 3, 5):
                c = int(h[i:i + 2], 16) / 255
                out.append(c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4)
            return out

        def lab(rgb):
            r, g, b = rgb
            x = (0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.95047
            y = 0.2126 * r + 0.7152 * g + 0.0722 * b
            z = (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.08883
            f = lambda t: t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116  # noqa: E731
            return 116 * f(y) - 16, 500 * (f(x) - f(y)), 200 * (f(y) - f(z))

        def sim(rgb, m):
            return [max(0.0, min(1.0, sum(m[i][j] * rgb[j] for j in range(3)))) for i in range(3)]

        classes = ("hy", "ar", "po", "ne", "pl", "sp", "st")
        for pal, floor in ((viewer._LIGHT, 12.0), (viewer._DARK, 10.0)):
            for name, m in sims.items():
                labs = {k: lab(sim(lin(pal[k]), m)) for k in classes}
                for a, b in itertools.combinations(classes, 2):
                    d = math.sqrt(sum((p - q) ** 2 for p, q in zip(labs[a], labs[b])))
                    self.assertGreaterEqual(d, floor, f"{name}: {a} vs {b} ({d:.1f})")
            lights = sorted(lab(lin(pal[k]))[0] for k in classes)
            for lo, hi in zip(lights, lights[1:]):
                self.assertGreaterEqual(hi - lo, 4.0, f"greyscale step {lo:.1f} -> {hi:.1f}")


class TestGui3Layout(unittest.TestCase):
    def setUp(self):
        saved = {k: sys.modules.get(k) for k in ("gradio", "coal_app.ui")}

        def restore():
            for k, v in saved.items():
                if v is None:
                    sys.modules.pop(k, None)
                else:
                    sys.modules[k] = v
        self.addCleanup(restore)
        sys.modules["gradio"] = _stub_gradio()
        sys.modules.pop("coal_app.ui", None)
        import importlib
        self.ui = importlib.import_module("coal_app.ui")
        _Comp.REG = []
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)

    def test_theme_uses_no_web_fonts(self):
        with open(self.ui.__file__, encoding="utf-8") as fh:
            self.assertNotIn("GoogleFont", fh.read().replace("never gr.themes.GoogleFont", ""))
        self.assertNotIn("@import", self.ui.CSS)
        self.assertNotRegex(self.ui.CSS, r"url\(")

    def test_strict_help_says_what_happens(self):
        self.assertIn("left out", self.ui.STRICT_INFO)
        self.assertIn("still aligned", self.ui.STRICT_INFO)
        self.assertNotIn("Stop with an error", self.ui.STRICT_INFO)

    def test_all_events_still_wired(self):
        self.ui.build(demo_mode=True)
        events = [ev for _, ev, _, _ in _Comp.REG]
        self.assertGreaterEqual(events.count("click"), 3 + len(core.available_examples()))  # run, sample, hint, ex
        self.assertEqual(events.count("upload"), 2)
        self.assertEqual(events.count("change"), 2)
        loads = [(a, kw) for c, ev, a, kw in _Comp.REG if ev == "load"]
        self.assertTrue(loads)
        # the URL hook fills both path boxes and opens their disclosure
        outs = [a[2] if len(a) > 2 else kw.get("outputs") for a, kw in loads]
        self.assertTrue(any(isinstance(o, list) and len(o) == 3 for o in outs), outs)

    def test_run_ui_returns_html_banner(self):
        p, d = core.example_paths("one_broken_gene")
        out = self.ui.run_ui(None, None, "fasta", 1, False, False, p, d)
        banner = out[0]._kw["value"]
        self.assertIn('class="coal-banner', banner)
        self.assertRegex(banner, r"of 20 genes aligned|COAL stopped")
        self.assertIn("read in place", banner)
        out = self.ui.run_ui(None, None, "fasta", 1, False, False)
        self.assertIn("Please add both files", out[0]._kw["value"])
        self.assertIn("coal-banner", self.ui.show_sample()[0]._kw["value"])
        self.assertIn("not your data", self.ui.show_sample()[0]._kw["value"])


class TestGui13LessText(unittest.TestCase):
    def setUp(self):
        saved = {k: sys.modules.get(k) for k in ("gradio", "coal_app.ui")}

        def restore():
            for k, v in saved.items():
                if v is None:
                    sys.modules.pop(k, None)
                else:
                    sys.modules[k] = v
        self.addCleanup(restore)
        sys.modules["gradio"] = _stub_gradio()
        sys.modules.pop("coal_app.ui", None)
        import importlib
        self.ui = importlib.import_module("coal_app.ui")
        _Comp.REG = []

    def test_banner_does_not_repeat_the_counts(self):
        rep = {"summary": {"total": 24, "warning": 24}, "sequences": [_seq(i, "warning") for i in range(24)]}
        sentence = "24 sequences: 0 OK, 24 with warnings, 0 failed."
        html = gui_html.run_banner("### ⚠ Done, with some problems (see the table)", rep, sentence, "Files used: x")
        self.assertNotIn("0 OK", html)
        self.assertNotIn("Done, with some problems", html)
        self.assertEqual(html.count("24 with warnings"), 1)
        removed = sentence + " 3 codon columns were removed because they didn't match (-nomismatch)."
        self.assertIn("3 codon columns were removed", gui_html.run_banner("### ⚠ Done", rep, removed))
        self.assertEqual(gui_html.extra_sentence("1 sequence: 1 OK, 0 with warnings, 0 failed."), "")

    def test_shared_closing_sentence_is_said_once(self):
        msg = [{"code": "mismatch", "column": 3, "text": "x"}]
        rep = {"sequences": [_seq(i, "warning", messages=msg) for i in range(5)]}
        tail = core._closing_sentence(core.sequence_sentence(rep["sequences"][0]))
        self.assertTrue(tail)
        html = core.health_table_html(rep, links=True)
        self.assertEqual(html.count(tail.replace("'", "&#x27;")), 1)
        self.assertIn("5 genes:", html)
        self.assertIn('href="#coalv-g00-c3"', html)
        self.assertIn("Column ", html)
        self.assertNotIn("Show in alignment", html)

    def test_form_has_short_labels_and_no_help_lines(self):
        self.ui.build()
        import gc
        comps = [o for o in gc.get_objects() if isinstance(o, _Comp) and o._name in ("Dropdown", "File", "Tab")]
        labels = {str(c._kw.get("label", c._args[0] if c._args else "")) for c in comps}
        self.assertTrue({"Protein alignment", "Coding DNA", "Genetic code"} <= labels, labels)
        with open(self.ui.__file__, encoding="utf-8") as fh:
            src = fh.read()
        self.assertIn('gr.Tab("Align"', src)
        for c in comps:
            # the main form's genetic-code box has one hint line (the Settings tab's has none)
            if c._name == "Dropdown" and c._kw.get("label") == "Genetic code":
                self.assertEqual(c._kw.get("info"), self.ui.CODE_HINT if c._kw.get("elem_id") == "coal-code"
                                 else None)
        self.assertIsNone(self.ui._format_info(core.Features(output_formats=("fasta", "paml"))))
        self.assertIn("Re-run", self.ui._strip_rerun("Use code 2. Re-run with code 2?", 2) + "Re-run")
        self.assertEqual(self.ui._strip_rerun("Use code 2. Re-run with code 2?", 2), "Use code 2.")


class TestLongWarningGroups(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import tempfile
        cls.tmp = tempfile.mkdtemp(prefix="coal_pre3_")
        p, d = core.example_paths("mito")
        cls.report = core.run_coal(p, d, "fasta", 1, False, False, out_dir=cls.tmp).report
        # with its genetic-code hint this run's warnings are one closed row (TestAmb4CodeFix); the grouping tested
        # here is for runs without such a hint, so the hint is left out
        cls.report = dict(cls.report, hints=[])

    @classmethod
    def tearDownClass(cls):
        import shutil
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_three_rows_then_show_more(self):
        h = core.health_table_html(self.report, links=True)
        self.assertIn("24 genes:</strong> Kept in the alignment; check the genetic code.", h)
        head, sep, tail = h.partition("<details><summary>Show 21 more</summary>")
        self.assertTrue(sep, h[:500])
        self.assertEqual(head.count('<tr class="coal-warning">'), 3)
        self.assertEqual(tail.count('<tr class="coal-warning">'), 21)
        self.assertEqual(h.count("Show "), 1)
        # every gene once, in the same order as before (problems first, then file order), with its links
        ids = re.findall(r"<td><code>([^<]+)</code>", h)
        self.assertEqual(ids, [s["id"] for s in self.report["sequences"]])
        for s in self.report["sequences"]:
            self.assertIn(f'href="#coalv-{s["id"]}-c', h)
        self.assertIn('<thead class="coal-sr">', tail)

    def test_three_or_fewer_rows_have_no_disclosure(self):
        import copy
        rep = copy.deepcopy(self.report)
        rep["sequences"] = rep["sequences"][:3]
        h = core.health_table_html(rep, links=True)
        self.assertNotIn("Show ", h)
        self.assertIn("3 genes:</strong>", h)
        rep["sequences"] = copy.deepcopy(self.report["sequences"][:4])
        self.assertIn("<summary>Show 1 more</summary>", core.health_table_html(rep, links=True))


if __name__ == "__main__":
    unittest.main()
