"""R2-D2: numbers quoted in the docs must equal what the harness and the repo say now."""
import glob
import importlib.util
import os
import re
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _harness():
    spec = importlib.util.spec_from_file_location("coal_test_differential", os.path.join(ROOT, "tests",
                                                                                        "test_differential.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _read(rel):
    with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
        return re.sub(r"\s+", " ", fh.read())  # line breaks don't matter


def _lines(pattern):
    total = 0
    for p in glob.glob(os.path.join(ROOT, pattern)):
        with open(p, encoding="utf-8") as fh:
            total += sum(1 for _ in fh)
    return total


class TestDocsNumbers(unittest.TestCase):
    longMessage = False  # the message says what is missing; do not dump whole documents
    @classmethod
    def setUpClass(cls):
        h = _harness()
        cls.ncases, cls.n, cls.c = h.compare_all(unittest.TestCase())
        cls.line = h.summary_line(cls.ncases, cls.n, cls.c)

    def expect(self, rel, *phrases):
        text = _read(rel)
        for ph in phrases:
            self.assertIn(ph, text, f"{rel} should say: {ph!r} (harness: {self.line})")

    def test_pitch_and_qa(self):
        c, n = self.c, self.n
        byte, warn, both = c["stdout byte-identical"], c["warnings identical"], c["both fail"]
        other = c["deliberate error"] + c["deliberate differs"]
        for rel in ("pitch/PITCH.md", "pitch/QA.md"):
            self.expect(rel, f"{n} checks on {self.ncases} hand-made inputs", f"{byte} compare",
                        f"{warn} of those also compare the warnings", f"{both} check that COAL rejects",
                        f"{other} are documented, deliberate differences")

    def test_skill(self):
        self.expect(".claude/skills/coal/SKILL.md", f"all {self.c['stdout byte-identical']} byte-for-byte comparisons")

    def test_readme(self):
        c = self.c
        self.expect("README.md", f"= {self.n} checks", f"byte-identical in {c['stdout byte-identical']} checks",
                    f"identical in {c['warnings identical']} of them ({c['with -nomismatch']} with `-nomismatch`)",
                    f"in {c['both fail']} both tools fail")

    def test_tests_readme_quotes_the_printed_line(self):
        self.expect("tests/README.md", self.line)

    def test_qa_line_and_test_counts(self):
        text = _read("pitch/QA.md")
        m = re.search(r"about ([\d,]+) lines of dependency-free Python in the engine and ([\d,]+) in the app", text)
        self.assertIsNotNone(m, "QA.md should state the engine and app line counts")
        for said, actual in ((m.group(1), _lines("src/coal/*.py")), (m.group(2), _lines("src/coal_app/*.py"))):
            said = int(said.replace(",", ""))
            self.assertLess(abs(said - actual), 0.1 * actual, f"QA.md says {said} lines, the code has {actual}")
        m = re.search(r"more than ([\d,]+) app tests", text)
        self.assertIsNotNone(m, "QA.md should say 'more than N app tests'")
        here = os.path.join(ROOT, "tests_app")
        actual = unittest.TestLoader().discover(here, top_level_dir=here).countTestCases()
        said = int(m.group(1).replace(",", ""))
        self.assertTrue(0.8 * actual <= said < actual, f"QA.md says more than {said} app tests; there are {actual}")


if __name__ == "__main__":
    unittest.main()
