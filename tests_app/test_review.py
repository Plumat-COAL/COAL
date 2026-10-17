"""Regression tests for the independent review findings (see overnight/app/LOG.md, "Review fixes")."""
import os
import re
import shutil
import tempfile
import unittest

from coal_app import core, viewer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _tmp(test):
    d = tempfile.mkdtemp(prefix="coal_test_")
    test.addCleanup(shutil.rmtree, d, True)
    return d


def _write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def _anchor_ids(page):
    return set(re.findall(r'id="(coalv-[^"]+)"', page))


class TestFinding1CodonOnlyMode(unittest.TestCase):
    """-nogap / -nomismatch: report columns are protein columns, the output has fewer columns."""

    def setUp(self):
        p, d = core.example_paths("mito")
        self.p, self.d = p, d
        self.r = core.run_coal(p, d, nogap=True, out_dir=_tmp(self), probe_codes=False)
        self.page = viewer.view(p, self.r.output_path, self.r.report, d)

    def test_no_cell_marked_by_comparing_codon_to_itself(self):
        # every outlined cell's title would say "X is not coded by <codon>"; translated from that very codon
        titles = re.findall(r'title="Column \d+: (\S) is not coded by (\w+)"', self.page)
        self.assertEqual(titles, [])
        self.assertNotRegex(self.page, r'<td class="[^"]*\bmm\b')

    def test_no_dangling_jump_links(self):
        ids = _anchor_ids(self.page)
        hrefs = re.findall(r'href="#(coalv-[^"]+-c\d+)"', self.page)
        self.assertEqual([h for h in hrefs if h not in ids], [])
        table = core.health_table_html(self.r.report, links=True)
        self.assertEqual(re.findall(r'href="#coalv-[^"]+-c\d+"', table), [])

    def test_explains_why_no_marks(self):
        self.assertIn("protein-alignment columns", self.page)


class TestFinding2TolerantStart(unittest.TestCase):
    def test_mismatched_first_codon_still_finds_frameshift(self):
        p, d = core.example_paths("one_broken_gene")
        recs = viewer.read_alignment(p)
        t = _tmp(self)
        lines = []
        for sid, row in recs:
            if sid == "OG0001_sp14":
                i = next(k for k, a in enumerate(row) if a.isalpha())
                self.assertEqual(row[i], "M")
                row = row[:i] + "L" + row[i + 1:]
            lines += [f">{sid}", row]
        p2 = _write(os.path.join(t, "protein.aln.fasta"), "\n".join(lines) + "\n")
        prot = dict(viewer.read_alignment(p2))
        dna = dict(viewer.read_alignment(d))
        ev = viewer.diagnose(prot["OG0001_sp14"], dna["OG0001_sp14"])["events"]
        self.assertEqual([(e["code"], e["column"]) for e in ev], [("frameshift", 102)])
        r = core.run_coal(p2, d, out_dir=t, probe_codes=False)
        text = core.explain(r.report)
        self.assertIn("Frameshift at column 102", text)
        self.assertNotIn("does not spell the start", text)

    def test_mito_ata_start_code1(self):
        # ATA for the first M (code 1 has no ATA start): start still found, codon marked as one mismatch
        prot = "MKLAVGHWSTRE"
        dna = "ATA" + "AAACTGGCTGTTGGTCATTGGTCTACTCGTGAA"
        d = viewer.diagnose(prot, dna, 1)
        self.assertEqual(d["events"], [])
        self.assertEqual(d["cells"][1], ("ATA", "mismatch"))
        self.assertEqual(d["matched"], 11)


class TestFinding3AnyWarningCounts(unittest.TestCase):
    def test_unknown_aa_warning_is_a_warning(self):
        t = _tmp(self)
        p = _write(os.path.join(t, "p.fasta"), ">s1\nMBKL\n>s2\nMAKL\n")
        d = _write(os.path.join(t, "d.fasta"), ">s1\nATGGATAAACTG\n>s2\nATGGCTAAACTG\n")
        r = core.run_coal(p, d, out_dir=t, probe_codes=False)
        self.assertIn("unknown AA type", r.stderr)
        self.assertEqual((r.report["summary"]["ok"], r.report["summary"]["warning"]), (1, 1))
        s1 = next(s for s in r.report["sequences"] if s["id"] == "s1")
        self.assertEqual(s1["status"], "warning")
        sentence = core.sequence_sentence(s1)
        self.assertNotIn("match codon for codon", sentence)
        self.assertIn("B", sentence)

    def test_unrecognised_warning_format_still_counts(self):
        t = _tmp(self)
        p = _write(os.path.join(t, "p.fasta"), ">s1\nMAKL\n")
        d = _write(os.path.join(t, "d.fasta"), ">s1\nATGGCTAAACTG\n")
        rep = core.build_fallback_report(p, d, os.path.join(t, "o.fasta"), "fasta", 1, False, False, False, 0,
                                         "coal: warning: s1 something new happened here\n", 0.1)
        s1 = next(s for s in rep["sequences"] if s["id"] == "s1")
        self.assertEqual(s1["status"], "warning")
        self.assertIn("something new happened here", core.sequence_sentence(s1))


if __name__ == "__main__":
    unittest.main()
