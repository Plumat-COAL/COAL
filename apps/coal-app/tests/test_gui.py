"""Gradio GUI tests. Skipped when Gradio is not installed (pip install -e "apps/coal-app[gui]")."""
import importlib.util
import os
import threading
import unittest

os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")
HAS_GRADIO = importlib.util.find_spec("gradio") is not None
# APP-16: the CI job that installs apps/coal-app[gui] sets COAL_REQUIRE_EXTRAS=gui, so a missing Gradio fails instead of skipping
REQUIRE_GUI = "gui" in os.environ.get("COAL_REQUIRE_EXTRAS", "").split(",")


@unittest.skipUnless(HAS_GRADIO or REQUIRE_GUI, "gradio not installed (pip install -e 'apps/coal-app[gui]')")
class TestGui(unittest.TestCase):
    def setUp(self):
        import shutil
        import tempfile
        from coal_app import core
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)

    def test_run_ui_function(self):
        from coal_app import core, ui
        p, d = core.example_paths("mito")
        out = ui.run_ui(p, d, "fasta", 1, False, False)
        self.assertEqual(out[3], 2)  # suggested genetic code from the hint

    def test_download_buttons_hand_out_the_alignment_and_the_report(self):
        # mutation testing (verify 2026-10-24, V-16): no test looked at the file behind "Download codon alignment"
        import json
        from coal_app import core, ui
        p, d = core.example_paths("clean")
        for fmt, first in (("fasta", ">"), ("paml", " ")):
            with self.subTest(fmt=fmt):
                out = ui.run_ui(p, d, fmt, 1, False, False)
                aln, rep = out[4]["value"], out[5]["value"]
                self.assertTrue(out[4]["visible"] and out[5]["visible"])
                self.assertTrue(aln.endswith(".codon." + fmt), aln)
                with open(aln, encoding="utf-8") as fh:
                    text = fh.read()
                self.assertTrue(text.startswith(first), text[:40])
                self.assertEqual(text.count("gene0") + text.count("gene1"), 8)  # the 8 genes of the example
                with open(rep, encoding="utf-8") as fh:
                    report = json.load(fh)
                self.assertEqual(report["schema"], "coal.report/1")
                self.assertEqual(report["summary"]["ok"], 8)
                self.assertNotEqual(aln, rep)

    def test_path_box_wins_and_is_not_copied(self):  # APP-L1
        from coal_app import core, ui
        p, d = core.example_paths("mito")
        out = ui.run_ui(None, None, "fasta", 1, False, False, p, d)
        self.assertEqual(out[3], 2)
        out = ui.run_ui(None, None, "fasta", 1, False, False, "/no/such/file.fasta", d)
        self.assertIn("does not exist", str(out[0]))

    def test_build_and_serve_briefly(self):
        from coal_app import ui
        timer = None
        demo = ui.build(demo_mode=True)
        try:
            # watchdog: close the server even if something below hangs
            timer = threading.Timer(60, demo.close)
            timer.start()
            demo.launch(server_name="127.0.0.1", server_port=None, prevent_thread_lock=True,
                        inbrowser=False, mcp_server=False)
            import urllib.request
            with urllib.request.urlopen(demo.local_url, timeout=20) as r:
                self.assertEqual(r.status, 200)
        finally:
            if timer:
                timer.cancel()
            demo.close()


if __name__ == "__main__":
    unittest.main()
