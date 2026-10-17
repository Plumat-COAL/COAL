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
        other = c["deliberate error"] + c["deliberate differs"] + c.get("deliberate nogap_code", 0)  # FZ-5
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

    def test_readme_exit_status_count(self):  # R3-D5
        n_exit = self.c["COAL exit status"]
        self.expect("README.md", f"Exit status in the other {n_exit}")
        self.assertNotIn(f"Exit status in all {self.n}", _read("README.md"))

    def test_readme_scaling_extrapolation(self):  # R3-D4
        import csv
        with open(os.path.join(ROOT, "benchmarks", "results", "bench_windows_20261016_2356.csv"),
                  encoding="utf-8") as fh:
            row = next(r for r in csv.DictReader(fh) if (r["n_seqs"], r["n_codons"]) == ("10000", "500"))
        t = float(row["pal2nal_median_s"])
        m100k, m500k = round(t * 10 / 60), round(t * 50 / 60)
        self.expect("README.md", f"10,000-sequence time ({t:.1f} s)", f"roughly {m100k} and {m500k} minutes")

    def test_compatibility_claim_is_qualified_everywhere(self):  # FZ-claim
        c = self.c
        claim = (f"byte-identical to PAL2NAL v14 on well-formed input: all {self.n} stored reference checks pass "
                 f"({c['stdout byte-identical']} byte-identical, {c['both fail']} where both tools fail, "
                 f"{c['deliberate error'] + c['deliberate differs'] + c['deliberate nogap_code']} documented deliberate "
                 f"differences)")
        for rel in ("README.md", "pitch/PITCH.md", "pitch/QA.md", ".claude/skills/coal/SKILL.md"):
            self.expect(rel, claim, "0 undocumented differences")
        self.expect("README.md", "5,000 cases: 20,000 comparisons, 14,368 byte-identical, 4,723 where both fail, 909")
        self.expect("tests/README.md", "| fasta | 5,000 | 4,096 | 612 | 292 | 0 |")
        for rel in ("README.md", "pitch/PITCH.md", "pitch/QA.md", ".claude/skills/coal/SKILL.md"):
            text = _read(rel)
            for m in re.finditer(r"[^.]*byte-identical to PAL2NAL[^.]*", text):
                self.assertTrue(any(q in m.group(0) for q in ("well-formed", "benchmark")), f"{rel}: {m.group(0)}")

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
