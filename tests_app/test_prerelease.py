"""Pre-release fixes found by running the app in real Gradio 5.49.1 (PRE-1..PRE-3).
See overnight/integration/CHANGES.md. Runs without Gradio (the stub of test_round3): the click sequences are
replayed through the recorded event wiring, so what is tested is what the buttons are wired to do."""
import os
import re
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from coal_app import core  # noqa: E402
from test_gui2 import UiCase  # noqa: E402
from test_round3 import _Comp  # noqa: E402


class Page:
    """The built page under the stub: component values and the click -> .then chains, as the browser runs them."""

    def __init__(self, ui):
        _Comp.REG = []
        ui.build()
        self.events = list(_Comp.REG)
        self.vals: dict[int, object] = {}
        self.comps = []
        seen = set()

        def add(c):
            if isinstance(c, _Comp) and id(c) not in seen and c._name != "chain":
                seen.add(id(c))
                self.comps.append(c)
        for c, _ev, a, kw in self.events:
            add(c)
            for io in (a[1:3] if len(a) > 1 else ()):
                for x in (io if isinstance(io, (list, tuple)) else [io]):
                    add(x)

    def find(self, **kw):
        for c in self.comps:
            if all((c._kw.get(k) if k != "text" else (c._args[0] if c._args else None)) == v
                   for k, v in kw.items()):
                return c
        raise LookupError(kw)

    def get(self, c):
        if id(c) in self.vals:
            return self.vals[id(c)]
        if "value" in c._kw:
            return c._kw["value"]
        if c._name == "State":
            return c._args[0] if c._args else None
        return "" if c._name == "Textbox" else None

    def _set(self, c, v):
        if isinstance(v, _Comp) and v._name == "update":
            if "value" in v._kw:
                self.vals[id(c)] = v._kw["value"]
        else:
            self.vals[id(c)] = v

    def _run(self, a):
        fn, inputs, outputs = a[0], a[1], a[2]
        args = [] if inputs is None else ([self.get(x) for x in inputs] if isinstance(inputs, (list, tuple))
                                          else [self.get(inputs)])
        res = fn(*args)
        if isinstance(outputs, (list, tuple)):
            for c, v in zip(outputs, res):
                self._set(c, v)
        elif outputs is not None:
            self._set(outputs, res)

    def fire(self, comp, event="click"):
        for c, ev, a, _kw in self.events:
            if c is comp and ev == event:
                self._run(a)
                for ch, ev2, a2, _ in self.events:  # the .then steps of this event, in order
                    if ch._name == "chain" and ch._origin == (comp, event) and ev2 == "then":
                        self._run(a2)


class TestPre1ExamplesSetTheirSettings(UiCase):
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
        pg.fire(self.ex["one_broken_gene"])  # ran with code 2 before PRE-1 (19 warnings, tip for code 4)
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



class TestPre2AccordionIcon(unittest.TestCase):
    def test_turned_icon_is_square(self):
        # Gradio 5.49.1 turns the closed accordion's icon by 90 degrees inside an overflow:auto block; a square
        # icon keeps its turned box inside it (no scrollbar under a collapsed accordion)
        import importlib
        import sys as _sys
        from test_round3 import _stub_gradio
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


class TestPre3LongGroups(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import tempfile
        cls.tmp = tempfile.mkdtemp(prefix="coal_pre3_")
        p, d = core.example_paths("mito")
        cls.report = core.run_coal(p, d, "fasta", 1, False, False, out_dir=cls.tmp).report

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
