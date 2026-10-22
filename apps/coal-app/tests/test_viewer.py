import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from coal_app import core, viewer  # noqa: E402
from helpers import CLEAN_D, CLEAN_P, ExampleVariantCase, RunsDirCase, _anchor_ids, _read_fasta, _tmp, _write, _write_fasta  # noqa: E402


def _tmpdir(test):
    d = tempfile.mkdtemp(prefix="coal_test_")
    test.addCleanup(shutil.rmtree, d, True)
    return d


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
        # the preview line
        self.assertIn("Preview: 200 of 300 sequences", frag)
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


class TestViewerNoFalseFailures(RunsDirCase):
    def test_mito_code1_diagnose_calls_no_break(self):
        from coal_app import viewer
        p, d = core.example_paths("mito")
        dna = dict(viewer.read_alignment(d))
        for pid, row in viewer.read_alignment(p):
            ev = viewer.diagnose(row, dna[pid], 1)["events"]
            self.assertEqual(ev, [], pid)

    def test_standalone_viewer_runs_the_engine(self):
        import subprocess, sys
        p, d = core.example_paths("mito")
        out = os.path.join(self.tmp(), "v.html")
        r = subprocess.run([sys.executable, "-m", "coal_app.viewer", p, "--dna", d, "--codontable", "1", "-o", out],
                           capture_output=True, text=True, timeout=120,
                           env=dict(os.environ, COAL_RUNS_DIR=core.RUNS_DIR))
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out, encoding="utf-8") as fh:
            page = fh.read()
        self.assertNotIn("Failed", page)
        self.assertNotIn("no longer spells", page)
        self.assertNotIn("COAL wrote no codon alignment", page)

    def test_engine_aligned_gene_is_never_failed(self):
        from coal_app import viewer
        p, d = core.example_paths("mito")
        rep = {"sequences": [{"id": pid, "status": "warning", "messages": []}
                             for pid, _ in viewer.read_alignment(p)]}
        rows, _, _ = viewer.build_rows(viewer.read_alignment(p), None, rep, viewer.read_alignment(d), 1)
        self.assertNotIn("failed", {r["status"] for r in rows})

    def test_real_frameshift_and_wrong_dna_still_found(self):
        from coal_app import viewer
        p, d = core.example_paths("one_broken_gene")
        prot, dna = dict(viewer.read_alignment(p)), dict(viewer.read_alignment(d))
        ev = viewer.diagnose(prot["OG0001_sp14"], dna["OG0001_sp14"])["events"]
        self.assertEqual(ev[0]["code"], "frameshift")
        import random
        rng = random.Random(1)
        row = prot["OG0001_sp01"]
        junk = dna["OG0001_sp01"][:60] + "".join(rng.choice("ACGT") for _ in range(len(dna["OG0001_sp01"]) - 60))
        ev = viewer.diagnose(row, junk)["events"]
        self.assertTrue(ev and ev[0]["code"] == "inconsistent", ev)


class TestViewerCodonOnlyMode(unittest.TestCase):
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


class TestViewerTolerantStart(unittest.TestCase):
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


class TestViewerWithoutCodonFile(unittest.TestCase):
    def test_mismatches_without_events_are_warnings(self):
        p, d = core.example_paths("mito")
        rows, _, _ = viewer.build_rows(viewer.read_alignment(p), None, None, viewer.read_alignment(d), 1)
        danio = next(r for r in rows if r["id"] == "Danio_COX1")
        self.assertEqual(danio["events"], [])
        self.assertTrue(any(f == "mismatch" for _, _, f in danio["cells"]))
        self.assertEqual(danio["status"], "warning")
        page = viewer.view(p, None, None, d)
        self.assertNotRegex(page, r"✔ OK</strong> <a[^>]*><code>Danio_COX1</code></a>: [^<]*Jump to")


class TestIupacCodons(unittest.TestCase):
    def test_ambiguous_codon_accepted_like_the_engine(self):
        from coal_app import gcode
        self.assertEqual(gcode.translate("GCN"), "A")
        self.assertEqual(gcode.translate("GAN"), "X")  # GAT/GAC = D, GAA/GAG = E
        self.assertEqual(gcode.translate("TAR"), "*")
        t = _tmp(self)
        p = _write(os.path.join(t, "p.fasta"), ">s1\nMAKL\n")
        d = _write(os.path.join(t, "d.fasta"), ">s1\nATGGCNAAACTG\n")
        r = core.run_coal(p, d, out_dir=t, probe_codes=False)
        self.assertEqual(r.exit_code, 0)
        self.assertEqual(r.report["summary"]["ok"], 1)
        page = viewer.view(p, r.output_path, r.report, d)
        self.assertNotIn("is not coded by", page)
        self.assertNotRegex(page, r'<td class="[^"]*\bmm\b')
        # and a genuinely ambiguous codon is still flagged
        self.assertTrue(viewer.diagnose("MDKL", "ATGGANAAACTG")["cells"][2][1] == "mismatch")


class TestViewerOtherFormats(ExampleVariantCase):
    def test_paml_and_codon_output_are_drawn(self):
        from coal_app import viewer
        for fmt in ("paml", "codon", "clustal"):
            with self.subTest(fmt=fmt):
                od = os.path.join(self.d, "v_" + fmt)
                r = core.run_coal(CLEAN_P, CLEAN_D, fmt, out_dir=od, probe_codes=False)
                frag, path = viewer.render_run(CLEAN_P, r.output_path, r.report, CLEAN_D, 1, od)
                self.assertIn("gene01", frag)
                self.assertNotIn("not FASTA or CLUSTAL", frag)
                self.assertNotIn(od, frag)  # no server path
                self.assertIn("ATG", frag)
                self.assertEqual([f for f in os.listdir(od) if f.startswith(".view")], [])

    def test_error_text_has_no_server_path(self):
        self.assertEqual(core.public_error(ValueError("/srv/run_x/a.paml: not FASTA"), ["/srv/run_x/a.paml"]),
                         "a.paml: not FASTA")


class TestViewerFrameshiftNote(ExampleVariantCase):
    def test_wide_frameshift_column_is_not_called_nogap(self):
        from coal_app import viewer

        # protein M4..., DNA with 4 extra bases after the first codon: a widened column (digit >= 4)
        p = _read_fasta(CLEAN_P)
        d = _read_fasta(CLEAN_D)
        p = [[n, s[0] + ("4" if i == 0 else "-") + s[1:]] for i, (n, s) in enumerate(p)]
        d[0][1] = d[0][1][:3] + "ACGT" + d[0][1][3:]
        pp, dd = os.path.join(self.d, "fs.pep"), os.path.join(self.d, "fs.dna")
        _write_fasta(pp, p)
        _write_fasta(dd, d)
        od = os.path.join(self.d, "fs")
        r = core.run_coal(pp, dd, out_dir=od, probe_codes=False)
        self.assertEqual(r.exit_code, 0, r.stderr)
        frag, _ = viewer.render_run(pp, r.output_path, r.report, dd, 1, od)
        self.assertNotIn("-nogap or -nomismatch", frag)


if __name__ == "__main__":
    unittest.main()
