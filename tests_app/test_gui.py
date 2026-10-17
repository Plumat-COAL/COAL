"""Gradio GUI tests. Skipped when Gradio is not installed (pip install -e .[gui])."""
import importlib.util
import os
import threading
import unittest

os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")
HAS_GRADIO = importlib.util.find_spec("gradio") is not None


@unittest.skipUnless(HAS_GRADIO, "gradio not installed")
class TestGui(unittest.TestCase):
    def test_run_ui_function(self):
        from coal_app import core, ui
        p, d = core.example_paths("mito")
        out = ui.run_ui(p, d, "fasta", 1, False, False)
        self.assertEqual(out[3], 2)  # suggested genetic code from the hint

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
