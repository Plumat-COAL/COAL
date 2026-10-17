import os
import subprocess
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestRaceDemo(unittest.TestCase):
    def test_missing_pal2nal_explains_and_exits_cleanly(self):
        r = subprocess.run([sys.executable, os.path.join(ROOT, "demo", "race.py"),
                            "--pal2nal", os.path.join(ROOT, "no", "such", "pal2nal.pl")],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("fetch_pal2nal.py", r.stdout)
        self.assertIn("not found", r.stdout)


if __name__ == "__main__":
    unittest.main()
