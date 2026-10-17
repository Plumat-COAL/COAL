import os
import shutil
import subprocess
import sys
import tempfile
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


STAND_IN = """import os, subprocess, sys
# stand-in for pal2nal.pl (NOT PAL2NAL): logs each call, then runs COAL on the same files
with open(os.environ["RACE_LOG"], "a", encoding="utf-8") as fh:
    fh.write(sys.argv[1] + "\\n")
r = subprocess.run([sys.executable, "-m", "coal", sys.argv[1], sys.argv[2]], capture_output=True, text=True)
sys.stdout.write(r.stdout); sys.stderr.write(r.stderr); sys.exit(r.returncode)
"""


class TestRaceIsPerFamily(unittest.TestCase):
    """Review finding 8: PAL2NAL is called once per gene family (not once per sequence)."""

    def test_stand_in_called_once_per_family_and_outputs_match(self):
        t = tempfile.mkdtemp(prefix="coal_test_")
        self.addCleanup(shutil.rmtree, t, True)
        stub = os.path.join(t, "stand_in_pal2nal.py")
        with open(stub, "w", encoding="utf-8") as fh:
            fh.write(STAND_IN)
        log = os.path.join(t, "calls.txt")
        env = dict(os.environ, RACE_LOG=log)
        r = subprocess.run([sys.executable, os.path.join(ROOT, "demo", "race.py"), "--pal2nal", stub,
                            "--perl", sys.executable, "--families", "6", "--seqs", "3", "--codons", "50"],
                           capture_output=True, text=True, timeout=300, env=env)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(log, encoding="utf-8") as fh:
            calls = fh.read().split()
        self.assertEqual(len(calls), 6)
        self.assertEqual(len(set(calls)), 6)
        self.assertIn("6 gene families (18 sequences)", r.stdout)
        self.assertIn("once per gene family", r.stdout)
        self.assertIn("output MATCHES PAL2NAL's", r.stdout)


if __name__ == "__main__":
    unittest.main()
