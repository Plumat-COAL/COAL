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
        # every `coal ...` command line in the skill (the common commands and the pipelines section: --keep-going and
        # --report are in the engine since 0.5.0, so there is no "after the engine update" section any more)
        self.assertNotIn("After the engine update", self.text)
        flags = set(re.findall(r"(?<![\w-])(--?[a-z][\w-]*)", "\n".join(
            line for line in self.text.splitlines() if line.startswith("coal "))))
        self.assertTrue(flags)
        self.assertIn("--keep-going", flags)
        for fl in flags:
            self.assertIn(fl, h, f"{fl} is in SKILL.md's commands but not in `coal --help`")
        self.assertIn("usage error", self.text)
        self.assertIn("coal --help", self.text)
        self.assertNotIn("not yet verified", self.text)  # the Windows exit codes were verified 2026-10-24

    def test_extras_in_separate_environments(self):
        self.assertIn("separate environments", self.text)


if __name__ == "__main__":
    unittest.main()
