"""GUI-1..GUI-6: the redesigned GUI (header, summary banner, health table grouping, viewer theming, layout).
Runs without Gradio: the layout is checked against the stub from test_round3. See overnight/integration/CHANGES.md"""
import os
import re
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from coal_app import core, gui_html, viewer  # noqa: E402
from test_round3 import _Comp, _stub_gradio  # noqa: E402


def _seq(i, status="ok", **kw):
    return {"id": f"g{i:02d}", "status": status, "messages": [], **kw}


class TestGui1Banner(unittest.TestCase):  # GUI-1
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

    def test_strip_is_capped_and_escaped(self):
        big = {"summary": {"total": 401, "ok": 401}, "sequences": [_seq(i) for i in range(401)]}
        self.assertNotIn("coal-strip", gui_html.run_banner("### ✔ Done", big, "s."))
        odd = {"summary": {"total": 1, "ok": 1}, "sequences": [{"id": '<b>"x', "status": "ok"}]}
        self.assertNotIn('<b>"x', gui_html.run_banner("### ✔ Done", odd, "s."))

    def test_header_has_no_external_resources(self):
        html = gui_html.header_html(core.intro_text())
        self.assertNotRegex(html, r"https?://(?!www\.w3\.org/2000/svg)")
        self.assertIn("<svg", html)
        self.assertIn(core.intro_text(), html)


class TestGui4HealthTable(unittest.TestCase):  # GUI-4
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


class TestGui2ViewerTheme(unittest.TestCase):  # GUI-2
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


class TestGui3Layout(unittest.TestCase):  # GUI-1, GUI-3, GUI-5
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


if __name__ == "__main__":
    unittest.main()
