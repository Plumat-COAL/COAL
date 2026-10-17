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

    # CHANGED(FIN-C2/FIN-C3/FIN-C10): was the "byte-identical to PAL2NAL v14 on well-formed input ... 23 documented
    # deliberate differences ... every difference is one where PAL2NAL silently produces wrong output" claim (false:
    # well-formed mito -codontable 2 -nogap differs) and the fuzz numbers of the run before FIN-C10. Now one claim built
    # from the harness's counts and the FIN-C10 fuzz run, and every breakdown must add up to the total.
    FUZZ = {"cases": 5000, "comparisons": 20000, "identical": 14269, "both": 4692, "deliberate": 1039,
            "rows": {"clustal": (4067, 600, 333), "paml": (4067, 600, 333), "fasta": (4067, 600, 333),
                     "codon": (2068, 2892, 40)}}

    def claim(self):
        c, f = self.c, self.FUZZ
        k = c["deliberate error"] + c["deliberate differs"] + c["deliberate nogap_code"]
        return ("COAL's output is identical to PAL2NAL v14's except for a short list of documented, deliberate "
                "differences (newer NCBI genetic-code tables for codes 3, 10, 11 and 13; `-nogap` removing the selected "
                "code's stop codons; the selected code's start codons for the first Met; and inputs where PAL2NAL v14 "
                f"writes malformed output or fails). Tested: all {self.n} stored reference checks pass "
                f"({c['stdout byte-identical']} byte-identical, {c['deliberate same_output']} of them with deliberately "
                f"different warnings; {c['both fail']} where both tools fail; {k} deliberate output or exit-status "
                f"differences: {c['deliberate error']} deliberate errors, {c['deliberate differs']} different outputs, "
                f"{c['deliberate nogap_code']} `-nogap` checks with a non-standard code), and in "
                f"{f['comparisons']:,} random comparisons with the real PAL2NAL v14 ({f['cases']:,} cases x 4 output "
                f"formats) {f['identical']:,} were byte-identical, {f['both']:,} both failed and {f['deliberate']:,} "
                f"were deliberate differences, each checked to be exactly the documented change; 0 unexplained "
                f"differences.")

    def test_compatibility_claim_everywhere(self):  # FZ-claim, FIN-C3
        c, f = self.c, self.FUZZ
        self.assertEqual(c["stdout byte-identical"] + c["both fail"] + c["deliberate error"] + c["deliberate differs"]
                         + c["deliberate nogap_code"], self.n, "the harness's own breakdown must add up")
        self.assertEqual(f["identical"] + f["both"] + f["deliberate"], f["comparisons"])
        self.assertEqual(tuple(map(sum, zip(*f["rows"].values()))), (f["identical"], f["both"], f["deliberate"]))
        for rel in ("README.md", "pitch/PITCH.md", "pitch/QA.md", ".claude/skills/coal/SKILL.md"):
            self.expect(rel, self.claim())
        self.expect("README.md", f"seed 20261017, {f['cases']:,} cases: {f['comparisons']:,} comparisons, "
                                 f"{f['identical']:,} byte-identical, {f['both']:,} where both fail, "
                                 f"{f['deliberate']:,} deliberate differences")
        for fmt, (i, b, d) in f["rows"].items():
            self.expect("tests/README.md", f"| {fmt} | {f['cases']:,} | {i:,} | {b:,} | {d:,} | 0 |")
        for rel in ("README.md", "tests/README.md", "pitch/PITCH.md", "pitch/QA.md", ".claude/skills/coal/SKILL.md",
                    "site/index.html", "src/coal_compat/cli.py"):
            text = _read(rel)
            self.assertNotIn("on well-formed input", text, f"{rel}: the well-formed claim is false (FIN-C3)")
            self.assertNotIn("stratified", text, f"{rel}: the stratified run is not reproducible (FIN-C20)")
            for m in re.finditer(r"[^.]*\bidentical to PAL2NAL\b[^.]*", text):
                self.assertRegex(m.group(0), r"PAL2NAL v14", f"{rel}: say PAL2NAL v14: {m.group(0)}")
                self.assertTrue(any(q in m.group(0) for q in ("except", "benchmark")), f"{rel}: {m.group(0)}")

    def test_breakdowns_add_up(self):  # FIN-C2 / FIN-A3
        c = self.c
        k = c["deliberate error"] + c["deliberate differs"] + c["deliberate nogap_code"]
        text = _read("README.md")
        m = re.search(r"stdout is byte-identical in (\d+) checks.*?in (\d+) both tools fail\. (\d+) checks are "
                      r"deliberate output or exit-status differences: (\d+) deliberate errors.*?, (\d+) with different "
                      r"output and (\d+) `-nogap` checks", text)
        self.assertIsNotNone(m, "README.md Tests: the breakdown sentence is missing")
        byte, both, kk, e, d, g = map(int, m.groups())
        self.assertEqual((byte, both, kk, e, d, g), (c["stdout byte-identical"], c["both fail"], k,
                                                     c["deliberate error"], c["deliberate differs"],
                                                     c["deliberate nogap_code"]))
        self.assertEqual(e + d + g, kk)
        self.assertEqual(byte + both + kk, self.n)
        site = _read("site/index.html")
        spans = dict((v, v) for v in re.findall(r'data-check="([^"]+)"', site))
        nums = {}
        for key, pat in (("total", r"^(\d+) checks$"), ("byte", r"^byte-identical in (\d+) checks$"),
                         ("both", r"^in (\d+) both tools fail$"),
                         ("delib", r"^(\d+) checks are deliberate output or exit-status differences$")):
            hits = [int(re.match(pat, v).group(1)) for v in spans if re.match(pat, v)]
            self.assertEqual(len(hits), 1, f"site/index.html: one data-check for {key}")
            nums[key] = hits[0]
        self.assertEqual(nums["byte"] + nums["both"] + nums["delib"], nums["total"], f"site numbers: {nums}")
        self.assertEqual(nums["total"], self.n)
        for rel in ("README.md", "pitch/PITCH.md", "pitch/QA.md", ".claude/skills/coal/SKILL.md", "site/index.html"):
            self.assertNotIn("29 checks are deliberate", _read(rel))

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
