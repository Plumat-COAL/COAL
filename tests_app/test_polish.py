"""GUI polish: no reassurance lines, section headings,
the genetic code first and emphasised, no empty results area, Gradio's own settings link hidden."""
import gc
import os
import re
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from coal_app import core, gui_html  # noqa: E402
from test_round3 import _Comp, _stub_gradio  # noqa: E402


class _UiCase(unittest.TestCase):
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
        with open(self.ui.__file__, encoding="utf-8") as fh:
            self.src = fh.read()

    def comps(self, name, **kw):
        return [o for o in gc.get_objects() if isinstance(o, _Comp) and o._name == name
                and all(o._kw.get(k) == v for k, v in kw.items())]


class TestPol1NoReassuranceLines(_UiCase):
    def test_splash_and_footer_lines_are_gone(self):
        with open(os.path.join(ROOT, "desktop", "coal_desktop.py"), encoding="utf-8") as fh:
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


class TestPol2Sections(_UiCase):
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


class TestPol3GeneticCodeFirst(_UiCase):
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


class TestPol4NoEmptyResults(_UiCase):
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


class TestPol6OneSettingsPlace(_UiCase):
    def test_gradio_settings_link_hidden_api_link_kept(self):
        css = self.ui.CSS
        self.assertIn("footer > button.settings", css)
        self.assertIn("footer > div.divider:nth-last-child(2)", css)
        self.assertNotIn("show-api", css)
        self.assertNotIn("built-with", css)

    def test_no_language_choice(self):
        for m in ("ui.py", "gui_html.py", "core.py", "settings.py"):
            with open(os.path.join(ROOT, "src", "coal_app", m), encoding="utf-8") as fh:
                text = fh.read()
            self.assertIsNone(re.search(r"\b(i18n|I18n|locale|Language\b)", text), m)
            self.assertEqual(re.findall(r"language=\"(\w+)\"", text), ["shell"] if m == "ui.py" else [])


if __name__ == "__main__":
    unittest.main()
