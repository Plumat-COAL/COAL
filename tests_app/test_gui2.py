"""GUI: settings, the alignment preview at scale, dropdown labels, drop zones.
See docs/ui-principles.md. Runs without Gradio (the stub of test_round3)."""
import json
import os
import re
import shutil
import sys
import tempfile
import time
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from coal_app import core, overview, seqformat, settings, viewer  # noqa: E402
from test_round3 import _Comp, _stub_gradio  # noqa: E402


class SettingsDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="coal_test_settings_")
        self.addCleanup(shutil.rmtree, self.dir, True)
        for k, v in (("COAL_SETTINGS_DIR", self.dir), ("COAL_HOSTED", None)):
            old = os.environ.get(k)
            self.addCleanup(lambda k=k, old=old: os.environ.__setitem__(k, old) if old is not None
                            else os.environ.pop(k, None))
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


class TestGui9Settings(SettingsDirCase):
    def test_defaults_without_a_file(self):
        data, note = settings.load()
        self.assertEqual(data, settings.DEFAULTS)
        self.assertEqual(note, "")

    def test_save_is_atomic_and_round_trips(self):
        data, err = settings.save({"codontable": 2, "output_format": "paml", "nogap": True, "preview": "small"})
        self.assertEqual(err, "")
        self.assertTrue(os.path.isfile(settings.settings_path()))
        self.assertEqual(sorted(os.listdir(self.dir)), ["settings.json"])  # no temporary file left behind
        back, note = settings.load()
        self.assertEqual(back, data)
        self.assertEqual((back["codontable"], back["output_format"], back["nogap"], back["preview"]),
                         (2, "paml", True, "small"))

    def test_damaged_or_odd_files_give_defaults(self):
        with open(settings.settings_path(), "w", encoding="utf-8") as fh:
            fh.write('{"codontable": 2, "nog')
        data, note = settings.load()
        self.assertEqual(data, settings.DEFAULTS)
        self.assertIn("damaged", note)
        with open(settings.settings_path(), "w", encoding="utf-8") as fh:
            json.dump({"codontable": 7, "output_format": "exe", "nogap": "yes", "preview": "huge", "x": 1,
                       "strict": True}, fh)
        data, note = settings.load()
        self.assertEqual(note, "")
        self.assertEqual(data, dict(settings.DEFAULTS, strict=True))  # only the valid value kept
        self.assertEqual(settings.clean([1, 2]), settings.DEFAULTS)

    def test_format_the_engine_cannot_write_is_not_kept(self):
        self.assertEqual(settings.clean({"output_format": "paml"}, formats=["fasta"])["output_format"], "fasta")

    def test_hosted_space_writes_no_file(self):
        os.environ["COAL_HOSTED"] = "1"
        data, err = settings.save({"codontable": 2})
        self.assertEqual((data["codontable"], err), (2, ""))
        self.assertEqual(os.listdir(self.dir), [])
        self.assertEqual(settings.load()[0], settings.DEFAULTS)

    def test_reset(self):
        settings.save({"codontable": 2})
        data, _ = settings.reset()
        self.assertEqual(data, settings.DEFAULTS)
        self.assertEqual(settings.load()[0], settings.DEFAULTS)

    def test_config_folder_per_platform(self):
        os.environ.pop("COAL_SETTINGS_DIR")
        d = settings.settings_dir()
        if sys.platform == "win32":
            self.assertTrue(d.endswith(os.path.join("", "COAL")))
        elif sys.platform == "darwin":
            self.assertIn("Application Support", d)
        else:
            self.assertTrue(d.endswith(os.path.join("", "coal")))

    def test_preview_limits(self):
        self.assertEqual(settings.preview_limits("large"), (viewer.MAX_SEQS, viewer.MAX_COLS, viewer.MAX_CELLS))
        self.assertEqual(settings.preview_limits("nonsense"), settings.PREVIEW_SIZES["standard"])


class UiCase(SettingsDirCase):
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
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)


class TestGui8Controls(UiCase):
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


class TestGui9SettingsUi(UiCase):
    def test_settings_tab_wired(self):
        self.ui.build()
        events = [(c, ev, a, kw) for c, ev, a, kw in _Comp.REG]
        inputs = [e for e in events if e[1] == "input"]
        self.assertEqual(len(inputs), 6)  # genetic code, format, three options, preview size
        outs = inputs[0][2][2]
        self.assertEqual(len(outs), 6 + 5 + 2)  # settings, main form, status, browser store
        loads = [e for e in events if e[1] == "load"]
        self.assertTrue(any(len(a[2]) == 13 for _, _, a, _ in loads if len(a) > 2 and isinstance(a[2], list)))
        clicks = [c._kw.get("value", c._args[0] if c._args else "") for c, ev, _, _ in events if ev == "click"]
        self.assertIn("Reset to defaults", clicks)
        self.assertIn("Show these columns", clicks)
        self.assertEqual([ev for _, ev, _, _ in events].count("change"), 2)  # unchanged (path boxes)

    def test_save_applies_to_the_form_and_persists(self):
        out = self.ui.save_settings_ui(2, "fasta", True, False, False, "small")
        self.assertEqual(list(out[:6]), [2, "fasta", True, False, False, "small"])
        self.assertEqual(list(out[6:11]), [2, "fasta", True, False, False])
        self.assertIn("Saved", out[11])
        self.assertEqual(out[12]["codontable"], 2)
        self.assertEqual(settings.load()[0]["preview"], "small")
        again = self.ui.load_settings_ui(None)
        self.assertEqual(again[0], 2)
        reset = self.ui.reset_settings_ui()
        self.assertEqual(reset[0], 1)
        self.assertIn("defaults", reset[11])

    def test_hosted_reads_the_browser_store(self):
        os.environ["COAL_HOSTED"] = "1"
        out = self.ui.load_settings_ui({"codontable": 11, "preview": "large"})
        self.assertEqual((out[0], out[5], out[6]), (11, "large", 11))
        self.assertIn("browser", self.ui._settings_where())
        self.assertEqual(os.listdir(self.dir), [])

    def test_run_uses_preview_size_and_offers_centring(self):
        p, d = core.example_paths("one_broken_gene")
        out = self.ui.run_ui(None, None, "fasta", 1, False, False, p, d, False, False, "small")
        self.assertEqual(len(out), 15)  # including the results area (shown by every result)
        ctx = out[12]
        self.assertEqual(ctx["preview"], "small")
        self.assertFalse(out[13]._kw.get("visible"))  # 20 genes, all columns drawn: no centring needed
        view, dl, table = self.ui.centre_ui(ctx, 30)
        self.assertIn("coalv", view._kw["value"])
        self.assertIn("coalv-OG0001_sp14", table._kw["value"])
        self.assertIn("Type a column number", self.ui.centre_ui(ctx, None)[0]._kw["value"])


def _synthetic(tmp, n, ncodons, problem_every=100):
    """Protein, codon output and report of n genes without running the engine (all rows alike)."""
    aa = ("MKV" + "ALG" * ncodons)[:ncodons]
    nt = ("ATGAAAGTT" + "GCTCTGGGT" * ncodons)[:3 * ncodons]
    pep, cod = os.path.join(tmp, "pep.fasta"), os.path.join(tmp, "codon.fasta")
    with open(pep, "w") as fp, open(cod, "w") as fc:
        for i in range(n):
            fp.write(f">g{i}\n{aa}\n")
            fc.write(f">g{i}\n{nt}\n")
    seqs = []
    for i in range(n):
        if i % problem_every == 7:
            seqs.append({"id": f"g{i}", "status": "warning", "messages": [
                {"code": "mismatch", "column": 1 + (i * 37) % ncodons, "text": "a mismatch"}]})
        else:
            seqs.append({"id": f"g{i}", "status": "ok", "messages": []})
    nw = sum(1 for s in seqs if s["status"] == "warning")
    report = {"sequences": seqs, "options": {"codontable": 1, "output_format": "fasta"},
              "summary": {"total": n, "ok": n - nw, "warning": nw, "failed": 0}}
    return pep, cod, report


class TestGui10Preview(unittest.TestCase):
    def test_select_rows_problems_first_then_context(self):
        ids = [f"g{i}" for i in range(1000)]
        rep = {"sequences": [{"id": f"g{i}", "status": "failed" if i in (500, 900) else
                              ("warning" if i % 50 == 3 else "ok")} for i in range(1000)]}
        chosen, sel = viewer.select_rows(ids, rep, 40)
        self.assertEqual(len(chosen), 40)
        self.assertEqual(chosen[:2], ["g500", "g900"])  # failed first
        self.assertTrue(all(i in ("g500", "g900") or int(i[1:]) % 50 == 3 for i in chosen[:22]))
        self.assertTrue(sel["problems_first"])
        self.assertEqual((sel["problems"], sel["problems_shown"], sel["context"]), (22, 22, 18))
        self.assertIn("g0", chosen)
        self.assertIn("g999", chosen)  # context spread through the whole file
        few, sel = viewer.select_rows(ids[:30], rep, 40)
        self.assertEqual(few, ids[:30])  # everything fits: file order
        many = {"sequences": [{"id": i, "status": "warning"} for i in ids]}
        chosen, sel = viewer.select_rows(ids, many, 40)
        self.assertEqual(sel["context"], 0)
        self.assertEqual(len(chosen), 40)

    def test_column_window(self):
        self.assertEqual(viewer.column_window(12000, 600), (1, 600))
        self.assertEqual(viewer.column_window(12000, 600, 3456), (3156, 3755))
        self.assertEqual(viewer.column_window(12000, 600, 11990), (11401, 12000))
        self.assertEqual(viewer.column_window(500, 600, 400), (1, 500))

    def test_preview_sentence(self):
        s = viewer.preview_sentence(40, 21800, (1, 600), 12000,
                                    {"problems_first": True, "problems": 30, "problems_shown": 30, "context": 10})
        self.assertEqual(s, "Preview: 40 of 21,800 sequences, problems first (all 30 with problems, then 10 OK ones "
                            "from across the file); columns 1–600 of 12,000.")

    def test_small_run_is_not_a_preview(self):
        p, d = core.example_paths("one_broken_gene")
        info = {}
        frag = viewer.view(p, None, None, d, page=False, info=info)
        self.assertFalse(info["preview"])
        self.assertNotIn("Preview:", frag)
        self.assertNotIn('class="overview"', frag)

    def test_large_run_preview_overview_and_links(self):
        tmp = tempfile.mkdtemp(prefix="coal_test_prev_")
        self.addCleanup(shutil.rmtree, tmp, True)
        pep, cod, rep = _synthetic(tmp, 3000, 1500, problem_every=60)
        info = {}
        frag = viewer.view(pep, cod, rep, None, page=False, max_seqs=40, max_cols=400, max_cells=16000, info=info)
        self.assertTrue(info["preview"])
        self.assertIn("Preview: 40 of 3,000 sequences, problems first", frag)
        self.assertIn('class="overview"', frag)
        self.assertLessEqual(frag.count("<path"), 8)  # one path per series, whatever the size
        self.assertNotIn("<rect", frag.replace('<rect class="ov-box"', ""))
        self.assertIn("Jalview", frag)
        c0, c1 = info["first"], info["last"]
        # every link in the fragment points at a cell that exists
        cells = set(re.findall(r'<td[^>]* id="([^"]+)"', frag))
        for target in re.findall(r'href="#(coalv-[^"]+-c\d+)"', frag):
            self.assertIn(target, cells)
        self.assertIn("(outside the drawn columns", frag)
        table = core.health_table_html(rep, links=True, preview=info)
        for target in re.findall(r'href="#(coalv-[^"]+-c\d+)"', table):
            self.assertIn(target, cells)
        self.assertIn("Centre the preview on column", table)
        # centred elsewhere
        info2 = {}
        frag2 = viewer.view(pep, cod, rep, None, page=False, max_seqs=40, max_cols=400, max_cells=16000,
                            centre=1200, info=info2)
        self.assertEqual((info2["first"], info2["last"]), (1000, 1399))
        self.assertIn("columns 1,000–1,399 of 1,500", frag2)
        self.assertNotEqual((c0, c1), (info2["first"], info2["last"]))
        # centred: the sequences with a problem in those columns are drawn first
        first_row = re.search(r'<tr id="coalv-([^"]+)"', frag2).group(1)
        col = next(m["column"] for sq in rep["sequences"] if sq["id"] == first_row for m in sq["messages"])
        self.assertTrue(1000 <= col <= 1399, (first_row, col))

    def test_rendering_is_fast_and_small_at_100000(self):
        tmp = tempfile.mkdtemp(prefix="coal_test_big_")
        self.addCleanup(shutil.rmtree, tmp, True)
        for n in (20_000, 100_000):
            pep, cod, rep = _synthetic(tmp, n, 300)
            info = {}
            t = time.perf_counter()
            frag, page = viewer.render_run(pep, cod, rep, None, 1, tmp, limits=settings.preview_limits("standard"),
                                           info=info)
            table = core.health_table_html(rep, links=True, preview=info)
            elapsed = time.perf_counter() - t
            self.assertLess(elapsed, 15.0, n)  # ~1 s here; generous for slow CI machines
            self.assertLess(len(frag.encode("utf-8")), 3_000_000, n)
            self.assertLess(os.path.getsize(page), 3_000_000, n)
            self.assertLess(len(table), 500_000, n)
            self.assertIn(f"of {n:,} sequences", frag)


class TestGui10Overview(unittest.TestCase):
    def test_gap_profile(self):
        rows = ["A-C-", "AC--", "----", "ACGT"]
        gaps, n = overview.gap_profile(rows, 4)
        self.assertEqual(n, 4)
        self.assertEqual(gaps, [0.25, 0.5, 0.5, 0.75])
        gaps, n = overview.gap_profile(["A-"] * 600, 3)  # several 255-row batches; the short end counts as gap
        self.assertEqual(gaps, [0.0, 1.0, 1.0])

    def test_problem_columns(self):
        rep = {"sequences": [{"messages": [{"column": 2}, {"column": 2}, {"column": 9}]},
                             {"messages": [{"column": 2}, {"column": "x"}]}]}
        self.assertEqual(overview.problem_columns(rep, 3), [0, 2, 0])

    def test_strips_are_bounded(self):
        st = ["ok"] * 100_000
        st[5] = "failed"
        svg = overview.sequences_strip(st, {0, 50_000})
        self.assertLessEqual(svg.count("<path"), 3)
        self.assertIn('viewBox="0 0 400 ', svg)
        self.assertIn("1 failed", svg)


class TestGui10FastReader(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
