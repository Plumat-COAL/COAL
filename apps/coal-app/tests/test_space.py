import os
import re
import subprocess
import sys
import tempfile
import unittest

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # apps/coal-app


class TestSpace(unittest.TestCase):
    def test_front_matter_matches_gui_pin(self):
        with open(os.path.join(APP, "space", "README.md"), encoding="utf-8") as fh:
            readme = fh.read()
        with open(os.path.join(APP, "pyproject.toml"), encoding="utf-8") as fh:
            pin = re.search(r'gradio\[mcp\]==([\d.]+)', fh.read()).group(1)
        self.assertTrue(readme.startswith("---\n"))
        self.assertIn("sdk: gradio\n", readme)
        self.assertIn(f"sdk_version: {pin}\n", readme)
        self.assertIn("app_file: app.py\n", readme)
        with open(os.path.join(APP, "space", "requirements.txt"), encoding="utf-8") as fh:
            self.assertIn(f"=={pin}", fh.read())

    def test_build_space_folder_runs_engine(self):
        with tempfile.TemporaryDirectory() as t:
            r = subprocess.run([sys.executable, os.path.join(APP, "space", "build_space.py"), t],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            for f in ("app.py", "README.md", "requirements.txt", "LICENSE",
                      "src/coal/cli.py", "src/coal_app/space_app.py", "examples/mito/dna.fasta"):
                self.assertTrue(os.path.isfile(os.path.join(t, f)), f)
            # the copied engine works on its own (as the Space's subprocess will run it)
            env = dict(os.environ, PYTHONPATH=os.path.join(t, "src"))
            r = subprocess.run([sys.executable, "-m", "coal", "--version"], capture_output=True, text=True,
                               timeout=60, env=env, cwd=t, stdin=subprocess.DEVNULL)
            self.assertEqual(r.returncode, 0, r.stderr)


if __name__ == "__main__":
    unittest.main()
