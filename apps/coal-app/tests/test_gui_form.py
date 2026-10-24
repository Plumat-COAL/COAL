"""The GUI form, built against the Gradio stand-in (no Gradio needed): controls and sections, upload vs path boxes, examples that set their own settings, the genetic-code fix banner, and no needless text. The click sequences are replayed through the recorded event wiring, so what is tested is what the buttons are wired to do. See docs/ui-principles.md."""
import os
import re
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from coal_app import core, gui_html  # noqa: E402
from helpers import APP, Page, RunsDirCase, StubUiCase, UiCase, _Comp, _stub_gradio  # noqa: E402


class TestFormControls(UiCase):
    def test_dropdown_labels_are_short(self):
        for label, n in self.ui._code_choices():
            self.assertLessEqual(len(label), 38, label)
            self.assertTrue(label.startswith(f"{n}: "), label)
        self.assertEqual(len(self.ui._code_choices()), len(core.GENETIC_CODES))
        f = core.detect_features()
        for label, _ in self.ui._format_choices(f):
            self.assertLessEqual(len(label), 20, label)
        self.assertNotIn("(nuclear genes of most organisms)", str(self.ui._code_choices()))

    def test_drop_zones_have_no_fixed_height(self):
        with open(self.ui.__file__, encoding="utf-8") as fh:
            src = fh.read()
        self.assertNotRegex(src, r"gr\.File\([^)]*height=")
        self.assertIn("button[aria-dropeffect] > div", self.ui.CSS)


class TestFormSections(StubUiCase):
    def test_one_heading_per_group_in_order(self):
        heads = re.findall(r'gr\.HTML\(_sec\("(\w+)"', self.src)
        self.assertEqual(self.ui._sec("Files"), '<h2 class="coal-sec">Files</h2>')
        for h in ("Files", "Output", "Results", "Defaults", "Preview", "About"):
            self.assertIn(h, heads)
        order = [self.src.index(s) for s in ('_sec("Files"', 'elem_id="coal-paths"', 'label="Genetic code"',
                                             '_sec("Output"', 'elem_id="coal-more"', 'elem_id="coal-run"',
                                             '_sec("Results"')]
        self.assertEqual(order, sorted(order))
        self.assertIn(".coal-sec {", self.ui.CSS)
        self.assertIn('block_title_text_color="#4A4744"', self.src)  # labels lighter than headings


class TestGeneticCodeFirst(StubUiCase):
    def test_code_is_its_own_section_and_format_secondary(self):
        self.ui.build()
        code = self.comps("Dropdown", elem_id="coal-code")[0]
        fmt = self.comps("Dropdown", elem_id="coal-format")[0]
        self.assertEqual(code._kw["label"], "Genetic code")
        self.assertEqual(code._kw["info"], self.ui.CODE_HINT)
        self.assertLessEqual(len(self.ui.CODE_HINT), 70)  # one line
        self.assertIs(fmt._kw.get("show_label"), False)
        self.assertEqual(fmt._kw["label"], "Output format")  # still its accessible name
        self.assertIn('#coal-code span[data-testid="block-info"]', self.ui.CSS)


class TestNoEmptyResults(StubUiCase):
    def test_results_hidden_until_a_result(self):
        self.ui.build()
        box = self.comps("Column", elem_id="coal-results")
        self.assertEqual(len(box), 1)
        self.assertFalse(box[0]._kw.get("visible"))
        self.assertEqual(gui_html.empty_banner(), "")
        p, d = core.example_paths("clean")
        for out in (self.ui.run_ui(None, None, "fasta", 1, False, False, p, d),
                    self.ui.run_ui(None, None, "fasta", 1, False, False),
                    self.ui.show_sample()):
            self.assertIs(out[-1]._kw.get("visible"), True)


class TestMethodsBoxIsLabelled(StubUiCase):
    def test_label_is_not_taken_by_a_copy_button_inside_it(self):
        # Gradio puts show_copy_button's button inside the <label>, which then labels the button, not the textarea
        # (axe: an unnamed text field): the box keeps a plain label and a separate button copies it
        self.ui.build()
        box = self.comps("Textbox", label="Methods paragraph")
        self.assertEqual(len(box), 1)
        self.assertFalse(box[0]._kw.get("show_copy_button"))
        self.assertEqual(box[0]._kw.get("elem_id"), "coal-methods")
        btn = self.comps("Button", elem_id="coal-methods-copy")
        self.assertEqual(len(btn), 1)
        self.assertEqual(btn[0]._args[0], "Copy methods paragraph")
        self.assertIn("navigator.clipboard.writeText", self.src)

    def test_hidden_duplicate_tab_bar_is_out_of_the_tab_order(self):
        # Gradio's second, visually hidden tab bar kept two focusable buttons (axe aria-hidden-focus)
        self.assertIn(".tab-container.visually-hidden, .tab-container.visually-hidden * { display:none !important; }",
                      self.ui.CSS)


class TestNoReassuranceLines(StubUiCase):
    def test_splash_and_footer_lines_are_gone(self):
        with open(os.path.join(APP, "desktop", "coal_desktop.py"), encoding="utf-8") as fh:
            desk = fh.read()
        desk = desk[desk.index('SPLASH = """'):]
        desk = desk[:desk.index('"""', 12) + 3]  # the splash page only (a comment around it may quote the old line)
        self.assertNotIn("No internet needed", desk)
        self.assertNotIn("This computer only", desk)
        self.assertIn("Starting COAL&hellip;", desk)
        self.assertEqual(core.INTRO_LOCAL, "")
        self.assertNotIn("Nothing leaves your computer", core.INTRO_HOSTED)
        self.assertIn("Hugging Face", core.INTRO_HOSTED)  # the hosted note is a fact people must know: kept

    def test_local_app_has_no_footer_note(self):
        if core.hosted():
            self.skipTest("hosted")
        self.assertEqual(core.intro_text(), "")
        self.assertIn("if core.intro_text():\n            gr.HTML(gui_html.footer_html", self.src)


class TestOneSettingsPlace(StubUiCase):
    def test_gradio_settings_link_hidden_api_link_kept(self):
        css = self.ui.CSS
        self.assertIn("footer > button.settings", css)
        self.assertIn("footer > div.divider:nth-last-child(2)", css)
        self.assertNotIn("show-api", css)
        self.assertNotIn("built-with", css)

    def test_no_language_choice(self):
        for m in ("ui.py", "gui_html.py", "core.py", "settings.py"):
            with open(os.path.join(APP, "src", "coal_app", m), encoding="utf-8") as fh:
                text = fh.read()
            self.assertIsNone(re.search(r"\b(i18n|I18n|locale|Language\b)", text), m)
            self.assertEqual(re.findall(r"language=\"(\w+)\"", text), ["shell"] if m == "ui.py" else [])


class TestUploadClearsPathBox(RunsDirCase):
    def setUp(self):
        super().setUp()
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

    def _label(self, comp):
        return str(comp._kw.get("label", ""))

    def test_upload_clears_its_path_box(self):
        self.ui.build()
        found = {}
        for comp, ev, a, kw in _Comp.REG:
            if ev != "upload":
                continue
            outs = kw.get("outputs", a[2] if len(a) > 2 else None)
            outs = outs if isinstance(outs, (list, tuple)) else [outs]
            fn = kw.get("fn", a[0] if a else None)
            found[self._label(comp)] = ([self._label(o) for o in outs], fn)
        # the labels are "Protein alignment" and "Coding DNA" (the formats are a hint inside the drop zone)
        prot = next((k for k in found if k.startswith("Protein")), None)
        dna = next((k for k in found if k.startswith("Coding DNA")), None)
        self.assertIsNotNone(prot, found)
        self.assertIsNotNone(dna, found)
        for k in (prot, dna):
            labels, fn = found[k]
            self.assertEqual(len(labels), 1, labels)
            # the path boxes sit under "Use a file path instead (large files)" and are labelled "<file> path"
            self.assertTrue(labels[0].endswith(" path"), labels)
            self.assertEqual(labels[0].split()[0], {prot: "Protein", dna: "Coding"}[k], labels)
            self.assertEqual(fn("/tmp/x.fa"), "")
            self.assertEqual(fn(None), "")

    def test_unusable_runs_dir_is_a_banner_not_an_error(self):
        core.RUNS_DIR = "/proc/nonexistent/runs" if os.name != "nt" else "?:\\nonexistent\\runs"
        p, d = core.example_paths("clean")
        out = self.ui.run_ui(None, None, "fasta", 1, False, False, p, d)
        banner = out[0]._kw["value"]
        self.assertIn("Could not create the results folder", banner)
        self.assertIn("nonexistent", banner)

    def test_launch_with_an_unusable_runs_dir_exits_with_a_message(self):
        import io
        from contextlib import redirect_stderr
        core.RUNS_DIR = "/proc/nonexistent/runs" if os.name != "nt" else "?:\\nonexistent\\runs"
        err = io.StringIO()
        with redirect_stderr(err), self.assertRaises(SystemExit) as cm:
            self.ui.launch()
        self.assertEqual(cm.exception.code, 1)
        self.assertIn("Could not create the results folder", err.getvalue())
        self.assertNotIn("Traceback", err.getvalue())

    def test_banner_names_the_files_used(self):
        p, d = core.example_paths("clean")
        out = self.ui.run_ui(None, None, "fasta", 1, False, False, p, d)
        banner = out[0]._kw["value"]
        self.assertIn("protein.aln.fasta", banner)
        self.assertIn("dna.fasta", banner)
        self.assertIn("read in place", banner)


class TestExamplesSetTheirSettings(UiCase):
    def setUp(self):
        super().setUp()
        self.page = Page(self.ui)
        pg = self.page
        self.code = pg.find(elem_id="coal-code")
        self.banner = pg.find(elem_id="coal-banner")
        self.hint_btn = pg.find(elem_id="coal-hint-rerun")
        self.table = pg.find(elem_id="coal-health-box")
        self.ex = {n: pg.find(text=core.EXAMPLES[n][0]) for n in core.available_examples()}
        self.opts = [pg.find(label=lbl) for lbl in ("Remove gap columns", "Remove mismatched columns",
                                                     "Refuse probably-wrong genes", "Keep only '#'-marked columns")]

    def test_owner_click_sequence(self):
        pg = self.page
        pg.fire(self.ex["mito"])  # deliberately code 1: 24 warnings and the tip
        self.assertEqual(pg.get(self.code), 1)
        self.assertIn("Re-run with genetic code 2", str(self.hint_btn and pg.vals.get(id(self.hint_btn))))
        self.assertIn("24 genes", pg.get(self.table))
        pg.fire(self.hint_btn)  # "Re-run with genetic code 2"
        self.assertEqual(pg.get(self.code), 2)
        self.assertIn("All 24 genes match codon for codon", pg.get(self.banner))
        pg.fire(self.ex["one_broken_gene"])  # must not keep code 2 from the previous run (19 warnings, tip for code 4)
        self.assertEqual(pg.get(self.code), 1)
        banner, table = pg.get(self.banner), pg.get(self.table)
        self.assertIn("19 of 20 genes aligned", banner)
        self.assertIn("19 genes match codon for codon", table)
        self.assertNotIn("check the genetic code", table)
        self.assertNotIn("coal-tile coal-warning", banner)
        self.assertEqual(pg.get(pg.find(elem_id="coal-hint")), "")  # no genetic-code tip

    def test_options_are_reset_and_shown(self):
        pg = self.page
        for c in self.opts:
            pg.vals[id(c)] = True
        pg.vals[id(self.code)] = 11
        pg.fire(self.ex["clean"])
        self.assertEqual([pg.get(c) for c in self.opts], [False] * 4)
        self.assertEqual(pg.get(self.code), 1)
        self.assertIn("All 8 genes match codon for codon", pg.get(self.banner))

    def test_own_files_keep_the_form(self):
        pg = self.page
        p, d = core.example_paths("mito")
        pg.vals[id(pg.find(label="Protein alignment"))] = p
        pg.vals[id(pg.find(label="Coding DNA"))] = d
        pg.vals[id(self.code)] = 2
        pg.fire(pg.find(elem_id="coal-run"))
        self.assertEqual(pg.get(self.code), 2)
        self.assertIn("All 24 genes match codon for codon", pg.get(self.banner))

    def test_every_example_has_settings(self):
        self.assertEqual(set(core.EXAMPLE_SETTINGS), set(core.EXAMPLES))
        self.assertEqual(core.EXAMPLE_SETTINGS["mito"]["codontable"], 1)  # the demo shows the tip


class TestGeneticCodeFixBanner(UiCase):
    """AMB-4: when the genetic-code hint explains the run, the result reads as "wrong code, here is the fix"."""

    def setUp(self):
        super().setUp()
        self.page = Page(self.ui)
        self.banner = self.page.find(elem_id="coal-banner")
        self.table = self.page.find(elem_id="coal-health-box")

    def test_mito_banner_names_cause_and_fix(self):
        pg = self.page
        pg.fire(pg.find(text=core.EXAMPLES["mito"][0]))
        b = pg.get(self.banner)
        self.assertIn('class="coal-banner coal-info coal-codefix" role="status"', b)  # informational, not a warning
        self.assertIn('aria-hidden="true">ℹ</span> Wrong genetic code?', b)  # icon + words
        self.assertIn("With code 1 (Standard) 24 of 24 genes have mismatched codons; with code 2 (Vertebrate "
                      "mitochondrial) all 24 match.", b)
        self.assertIn('<p class="coal-small">Codes 5 and 13 fit these sequences too', b)  # secondary line
        self.assertNotIn("coal-tile", b)
        self.assertNotIn("with warnings", b)
        self.assertEqual(pg.get(pg.find(elem_id="coal-hint")), "")  # the tip is in the banner now
        # the 24 warnings: one closed row, every gene inside it with its jump links
        t = pg.get(self.table)
        self.assertIn("<details><summary>24 genes with mismatched codons under code 1: expected with the wrong "
                      "code</summary>", t)
        self.assertEqual(t.count("coal-codefix-row"), 1)
        self.assertNotIn("Show ", t)
        self.assertEqual(t.count('<tr class="coal-warning">'), 24)

    def test_rerun_is_the_primary_action_in_the_banner(self):
        import ast
        import inspect
        p, d = core.example_paths("mito")
        rep = core.run_coal(p, d, "fasta", 1, False, False, out_dir=core.new_run_dir()).report
        md, btn, code = self.ui._hint_parts(rep)
        self.assertEqual(code, 2)
        self.assertEqual(btn._kw, {"value": "Re-run with genetic code 2", "visible": True, "variant": "primary"})
        self.assertEqual(md._kw.get("visible"), False)
        # the button and the banner are in one column, drawn as one box when the banner is the code fix
        import textwrap
        tree = ast.parse(textwrap.dedent(inspect.getsource(self.ui.build)))
        box = next(w for w in ast.walk(tree) if isinstance(w, ast.With)
                   and "coal-fixbox" in ast.unparse(w.items[0].context_expr))
        inside = ast.unparse(box)
        self.assertIn("elem_id='coal-banner'", inside)
        self.assertIn("elem_id='coal-hint-rerun'", inside)
        self.assertIn("#coal-fixbox:has(.coal-codefix) {", self.ui.CSS)

    def test_other_hints_and_runs_unchanged(self):
        pg = self.page
        pg.fire(pg.find(text=core.EXAMPLES["one_broken_gene"][0]))  # no hint: as before
        b = pg.get(self.banner)
        self.assertNotIn("coal-codefix", b)
        self.assertIn("19 of 20 genes aligned", b)
        rep = {"options": {"codontable": 1}, "summary": {"total": 10, "warning": 10},
               "hints": [{"code": "genetic_code", "text": "t", "suggested_codontable": 2, "checked": 10,
                          "fit": 6, "fit_current": 0}]}
        self.assertIsNone(core.code_fix(rep))  # 6 of 10 fit: a tip as before, not the code-fix banner
        md, btn, code = self.ui._hint_parts(rep)
        self.assertEqual(btn._kw.get("variant"), "secondary")
        self.assertEqual(md._kw.get("visible"), True)

    def test_one_broken_gene_with_code_2(self):
        # 19 of 20 fit code 1, one gene has a frameshift: the code-fix banner, with the failure still counted
        p, d = core.example_paths("one_broken_gene")
        rep = core.run_coal(p, d, "fasta", 2, False, False, out_dir=core.new_run_dir()).report
        fix = core.code_fix(rep)
        self.assertEqual((fix["current"], fix["code"], fix["fit"], fix["checked"]), (2, 1, 19, 20))
        from coal_app import gui_html
        b = gui_html.code_fix_banner(fix, rep, self.ui.CODE_SHORT)
        self.assertIn("with code 1 (Standard) 19 of 20 match.", b)
        self.assertIn("✖</span> 1 failed", b)
        t = core.health_table_html(rep)
        self.assertLess(t.index("<strong>Failed</strong>"), t.index("coal-codefix-row"))  # real problem first


class TestAccordionIcon(unittest.TestCase):
    def test_turned_icon_is_square(self):
        # Gradio 5.49.1 turns the closed accordion's icon by 90 degrees inside an overflow:auto block; a square
        # icon keeps its turned box inside it (no scrollbar under a collapsed accordion)
        import importlib
        import sys as _sys
        from helpers import _stub_gradio
        saved = _sys.modules.get("gradio"), _sys.modules.get("coal_app.ui")
        _sys.modules["gradio"] = _stub_gradio()
        _sys.modules.pop("coal_app.ui", None)
        try:
            css = importlib.import_module("coal_app.ui").CSS
        finally:
            for k, v in zip(("gradio", "coal_app.ui"), saved):
                if v is None:
                    _sys.modules.pop(k, None)
                else:
                    _sys.modules[k] = v
        rule = re.search(r"button\.label-wrap > span\.icon \{([^}]*)\}", css)
        self.assertIsNotNone(rule)
        w = re.search(r"width:([^;]+);", rule.group(1)).group(1)
        h = re.search(r"height:([^;]+);", rule.group(1)).group(1)
        self.assertEqual(w, h)
        self.assertIn("flex:none", rule.group(1))


if __name__ == "__main__":
    unittest.main()
