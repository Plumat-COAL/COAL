import os
import re
import subprocess
import sys
import tempfile
import shutil
import unittest

from coal_app import core, viewer


def _tmpdir(test):
    d = tempfile.mkdtemp(prefix="coal_test_")
    test.addCleanup(shutil.rmtree, d, True)
    return d


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def coal(*args):
    return subprocess.run([sys.executable, "-m", "coal", *args], capture_output=True, text=True, timeout=120)


class TestViewer(unittest.TestCase):
    def test_mito_code1_marks_tga_mismatches(self):
        p, d = core.example_paths("mito")
        with tempfile.TemporaryDirectory() as t:
            out = os.path.join(t, "codon.fa")
            coal(p, d, "-o", out)
            page = viewer.view(p, out)
        # Homo_COX1 column 36 is W coded by TGA (a stop in code 1): outlined, with the symbol and an anchor
        m = re.search(r'<td class="([^"]*)" id="coalv-Homo_COX1-c36"[^>]*><b>W</b><i>TGA</i>', page)
        self.assertIsNotNone(m)
        self.assertIn("mm", m.group(1).split())
        self.assertIn('href="#coalv-Homo_COX1-c36"', page)
        self.assertIn("\\2260", page)  # the mismatch symbol in CSS

    def test_mito_code2_has_no_marks(self):
        p, d = core.example_paths("mito")
        with tempfile.TemporaryDirectory() as t:
            out = os.path.join(t, "codon.fa")
            coal(p, d, "-codontable", "2", "-o", out)
            page = viewer.view(p, out, codontable=2)
        self.assertNotIn('class="hy mm"', page)
        self.assertNotRegex(page, r'<td class="[^"]*\bmm\b')
        self.assertIn("Every codon matches", page)

    def test_broken_gene_shows_frameshift(self):
        p, d = core.example_paths("one_broken_gene")
        r = core.run_coal(p, d, out_dir=_tmpdir(self))
        page = viewer.view(p, r.output_path, r.report, d)
        m = re.search(r'<td class="([^"]*)" id="coalv-OG0001_sp14-c(\d+)"', page)
        self.assertIsNotNone(m)
        self.assertIn("fs", m.group(1).split())
        col = int(m.group(2))
        # the report and the health table point at the same column
        msg = [x for s in r.report["sequences"] if s["id"] == "OG0001_sp14" for x in s["messages"]][0]
        self.assertEqual((msg["code"], msg["column"]), ("frameshift", col))
        self.assertIn(f'href="#coalv-OG0001_sp14-c{col}"', core.health_table_html(r.report, links=True))
        self.assertIn("Frameshift at column", core.explain(r.report))
        self.assertEqual(len(re.findall(r'\bfs\b"', page)), 1)  # only the broken gene

    def test_diagnose_good_and_shifted(self):
        rows = dict(viewer.read_alignment(core.example_paths("one_broken_gene")[0]))
        dna = dict(viewer.read_alignment(core.example_paths("one_broken_gene")[1]))
        self.assertEqual(viewer.diagnose(rows["OG0001_sp01"], dna["OG0001_sp01"])["events"], [])
        ev = viewer.diagnose(rows["OG0001_sp14"], dna["OG0001_sp14"])["events"]
        self.assertEqual([(e["code"], e["shift"]) for e in ev], [("frameshift", -1)])

    def test_self_contained(self):
        p, d = core.example_paths("clean")
        page = viewer.view(p, None, None, d)
        self.assertNotRegex(page, r'(src|href)="(https?:)?//')
        self.assertNotIn("@import", page)
        self.assertTrue(page.startswith("<!doctype html>"))

    def test_large_input_is_capped(self):
        rows = [{"id": f"s{i}", "status": "ok", "note": "", "events": [],
                 "cells": [("A", "GCT", "")] * 2500} for i in range(300)]
        frag = viewer.render_fragment(rows, 2500)
        self.assertIn("Large alignment: showing the first 200 of 300 sequences", frag)
        self.assertLessEqual(frag.count("<td"), viewer.MAX_CELLS)
        self.assertNotIn('id="coalv-s200"', frag)

    def test_cli(self):
        p, d = core.example_paths("one_broken_gene")
        with tempfile.TemporaryDirectory() as t:
            out = os.path.join(t, "view.html")
            r = subprocess.run([sys.executable, "-m", "coal_app.viewer", p, "--dna", d, "-o", out],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(out, encoding="utf-8") as fh:
                self.assertIn("coalv-OG0001_sp14", fh.read())

    def test_nogap_output_still_renders(self):
        p, d = core.example_paths("clean")
        with tempfile.TemporaryDirectory() as t:
            out = os.path.join(t, "codon.fa")
            coal(p, d, "-nogap", "-o", out)
            page = viewer.view(p, out)
        self.assertIn("Columns were removed", page)


if __name__ == "__main__":
    unittest.main()
