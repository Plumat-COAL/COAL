"""Settings: the JSON file (atomic, damaged files, hosted Space writes none) and the Settings tab."""
import json
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from coal_app import core, settings, viewer  # noqa: E402
from helpers import SettingsDirCase, UiCase, _Comp  # noqa: E402


class TestSettingsFile(SettingsDirCase):
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

    def test_non_utf8_file_is_damaged_not_a_crash(self):
        with open(settings.settings_path(), "wb") as fh:
            fh.write(b"\xff\xfe\x00garbage")
        data, note = settings.load()
        self.assertEqual(data, settings.DEFAULTS)
        self.assertIn("damaged", note)

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


class TestSettingsTab(UiCase):
    def test_settings_tab_wired(self):
        self.ui.build()
        events = [(c, ev, a, kw) for c, ev, a, kw in _Comp.REG]
        inputs = [e for e in events if e[1] == "input"]
        self.assertEqual(len(inputs), 7)  # genetic code, format, three options, preview size, residue colours
        outs = inputs[0][2][2]
        self.assertEqual(len(outs), 7 + 5 + 2)  # settings, main form, status, browser store
        loads = [e for e in events if e[1] == "load"]
        self.assertTrue(any(len(a[2]) == 14 for _, _, a, _ in loads if len(a) > 2 and isinstance(a[2], list)))
        clicks = [c._kw.get("value", c._args[0] if c._args else "") for c, ev, _, _ in events if ev == "click"]
        self.assertIn("Reset to defaults", clicks)
        self.assertIn("Show these columns", clicks)
        self.assertEqual([ev for _, ev, _, _ in events].count("change"), 2)  # unchanged (path boxes)

    def test_save_applies_to_the_form_and_persists(self):
        out = self.ui.save_settings_ui(2, "fasta", True, False, False, "small", False)
        self.assertEqual(list(out[:7]), [2, "fasta", True, False, False, "small", False])
        self.assertEqual(list(out[7:12]), [2, "fasta", True, False, False])
        self.assertIn("Saved", out[12])
        self.assertEqual(out[13]["codontable"], 2)
        self.assertEqual(settings.load()[0]["preview"], "small")
        self.assertIs(settings.load()[0]["residue_colours"], False)  # the new key round-trips
        again = self.ui.load_settings_ui(None)
        self.assertEqual(again[0], 2)
        reset = self.ui.reset_settings_ui()
        self.assertEqual(reset[0], 1)
        self.assertIs(reset[6], True)
        self.assertIn("defaults", reset[12])

    def test_hosted_reads_the_browser_store(self):
        os.environ["COAL_HOSTED"] = "1"
        out = self.ui.load_settings_ui({"codontable": 11, "preview": "large", "residue_colours": "no"})
        self.assertEqual((out[0], out[5], out[6], out[7]), (11, "large", True, 11))  # a non-bool is the default
        self.assertIn("browser", self.ui._settings_where())
        self.assertEqual(os.listdir(self.dir), [])

    def test_residue_colours_off_draws_letters_only(self):
        # the Settings tab's "Colour residues by class" off: the view and its download page carry the coalv-plain
        # class (no residue tints; the problem marks and the key's mark entries stay), also after re-centring
        p, d = core.example_paths("one_broken_gene")
        out = self.ui.run_ui(None, None, "fasta", 1, False, False, p, d, False, False, "small", False)
        view = out[9]._kw["value"]
        self.assertIn('<div class="coalv coalv-plain">', view)
        self.assertIn("coalv-plain td.hy", view)  # the rule that blanks the tints is in the fragment's own style
        self.assertIn('aria-label="Key"', view)
        self.assertIn("frameshift (double outline, hatched)", view)
        self.assertIn("coalv-OG0001_sp14-c102", view)  # the problem marks are still drawn
        with open(out[10]._kw["value"], encoding="utf-8") as fh:
            self.assertIn('class="coalv coalv-plain"', fh.read())
        ctx = out[13]
        self.assertFalse(ctx["colours"])
        view2 = self.ui.centre_ui(ctx, 30)[0]._kw["value"]
        self.assertIn('<div class="coalv coalv-plain">', view2)
        on = self.ui.run_ui(None, None, "fasta", 1, False, False, p, d, False, False, "small", True)
        self.assertNotIn("coalv-plain\">", on[9]._kw["value"])
        self.assertIn('aria-label="Colour key"', on[9]._kw["value"])

    def test_run_uses_preview_size_and_offers_centring(self):
        p, d = core.example_paths("one_broken_gene")
        out = self.ui.run_ui(None, None, "fasta", 1, False, False, p, d, False, False, "small")
        self.assertEqual(len(out), 16)  # including the results area (shown by every result) and the methods copy button
        ctx = out[13]
        self.assertEqual(ctx["preview"], "small")
        self.assertFalse(out[14]._kw.get("visible"))  # 20 genes, all columns drawn: no centring needed
        self.assertTrue(ctx["colours"])
        view, dl, table = self.ui.centre_ui(ctx, 30)
        self.assertIn("coalv", view._kw["value"])
        self.assertNotIn('<div class="coalv coalv-plain">', view._kw["value"])  # (the style sheet names the class)
        self.assertIn("coalv-OG0001_sp14", table._kw["value"])
        self.assertIn("Type a column number", self.ui.centre_ui(ctx, None)[0]._kw["value"])


if __name__ == "__main__":
    unittest.main()
