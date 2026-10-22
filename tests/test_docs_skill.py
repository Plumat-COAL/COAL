"""The Claude skill (.claude/skills/coal/SKILL.md) only uses what the engine has."""
import os
import re
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestSkillFile(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(ROOT, ".claude", "skills", "coal", "SKILL.md"), encoding="utf-8") as fh:
            self.text = fh.read()

    def test_current_section_only_uses_flags_the_engine_has(self):
        import subprocess
        import sys
        h = subprocess.run([sys.executable, "-m", "coal", "--help"], capture_output=True, text=True,
                           encoding="utf-8", timeout=60).stdout
        marker = "## After the engine update"
        self.assertIn(marker, self.text)
        current = self.text.split(marker)[0]
        flags = set(re.findall(r"(?<![\w-])(--?[a-z][\w-]*)", "\n".join(
            line for line in current.splitlines() if line.startswith("coal "))))
        self.assertTrue(flags)
        for fl in flags:
            self.assertIn(fl, h, f"{fl} is in SKILL.md's current commands but not in `coal --help`")
        for later in ("--keep-going", "--report", "normalised"):  # -output paml is current
            self.assertNotIn(later, current, later)
        self.assertIn("usage error", current)
        self.assertIn("coal --help", current)

    def test_extras_in_separate_environments(self):
        self.assertIn("separate environments", self.text)


if __name__ == "__main__":
    unittest.main()
