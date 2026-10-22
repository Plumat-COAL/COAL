"""Texts in desktop/ that must match what the wrapper does: the smoke-test checklist and the installer's "Open with" task."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from helpers import APP  # noqa: E402


class TestSmokeTestMatchesWrapper(unittest.TestCase):
    def test_open_with_and_outside_drop_load_by_path(self):
        with open(os.path.join(APP, "desktop", "SMOKE_TEST.md"), encoding="utf-8") as fh:
            row = next(ln for ln in fh if ln.startswith("| 6 |"))
        self.assertNotIn("doesn't load the file yet", row)
        self.assertNotIn("box to drag it onto", row)
        self.assertIn("puts the file into the right box", row)
        self.assertIn("Only one DNA file can be loaded at a time".lower(), row.lower())


class TestInstallerText(unittest.TestCase):
    def test_open_with_text_is_current(self):
        with open(os.path.join(APP, "desktop", "installer.iss"), encoding="utf-8") as fh:
            iss = fh.read()
        line = next(ln for ln in iss.splitlines() if ln.startswith('Name: "openwith"'))
        self.assertNotIn("can't load the file by itself", line)
        self.assertIn("filled in", line)


if __name__ == "__main__":
    unittest.main()
