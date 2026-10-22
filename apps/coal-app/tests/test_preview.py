"""The alignment preview at scale: problem rows first, column windows, the overview strips, and rendering 100,000 sequences fast and small."""
import os
import re
import shutil
import sys
import tempfile
import time
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from coal_app import core, overview, settings, viewer  # noqa: E402


def _synthetic(tmp, n, ncodons, problem_every=100):
    """Protein, codon output and report of n genes without running the engine (all rows alike)."""
    aa = ("MKV" + "ALG" * ncodons)[:ncodons]
    nt = ("ATGAAAGTT" + "GCTCTGGGT" * ncodons)[:3 * ncodons]
    pep, cod = os.path.join(tmp, "pep.fasta"), os.path.join(tmp, "codon.fasta")
    with open(pep, "w") as fp, open(cod, "w") as fc:
        for i in range(n):
            fp.write(f">g{i}\n{aa}\n")
            fc.write(f">g{i}\n{nt}\n")
    seqs = []
    for i in range(n):
        if i % problem_every == 7:
            seqs.append({"id": f"g{i}", "status": "warning", "messages": [
                {"code": "mismatch", "column": 1 + (i * 37) % ncodons, "text": "a mismatch"}]})
        else:
            seqs.append({"id": f"g{i}", "status": "ok", "messages": []})
    nw = sum(1 for s in seqs if s["status"] == "warning")
    report = {"sequences": seqs, "options": {"codontable": 1, "output_format": "fasta"},
              "summary": {"total": n, "ok": n - nw, "warning": nw, "failed": 0}}
    return pep, cod, report


class TestAlignmentPreview(unittest.TestCase):
    def test_select_rows_problems_first_then_context(self):
        ids = [f"g{i}" for i in range(1000)]
        rep = {"sequences": [{"id": f"g{i}", "status": "failed" if i in (500, 900) else
                              ("warning" if i % 50 == 3 else "ok")} for i in range(1000)]}
        chosen, sel = viewer.select_rows(ids, rep, 40)
        self.assertEqual(len(chosen), 40)
        self.assertEqual(chosen[:2], ["g500", "g900"])  # failed first
        self.assertTrue(all(i in ("g500", "g900") or int(i[1:]) % 50 == 3 for i in chosen[:22]))
        self.assertTrue(sel["problems_first"])
        self.assertEqual((sel["problems"], sel["problems_shown"], sel["context"]), (22, 22, 18))
        self.assertIn("g0", chosen)
        self.assertIn("g999", chosen)  # context spread through the whole file
        few, sel = viewer.select_rows(ids[:30], rep, 40)
        self.assertEqual(few, ids[:30])  # everything fits: file order
        many = {"sequences": [{"id": i, "status": "warning"} for i in ids]}
        chosen, sel = viewer.select_rows(ids, many, 40)
        self.assertEqual(sel["context"], 0)
        self.assertEqual(len(chosen), 40)

    def test_column_window(self):
        self.assertEqual(viewer.column_window(12000, 600), (1, 600))
        self.assertEqual(viewer.column_window(12000, 600, 3456), (3156, 3755))
        self.assertEqual(viewer.column_window(12000, 600, 11990), (11401, 12000))
        self.assertEqual(viewer.column_window(500, 600, 400), (1, 500))

    def test_preview_sentence(self):
        s = viewer.preview_sentence(40, 21800, (1, 600), 12000,
                                    {"problems_first": True, "problems": 30, "problems_shown": 30, "context": 10})
        self.assertEqual(s, "Preview: 40 of 21,800 sequences, problems first (all 30 with problems, then 10 OK ones "
                            "from across the file); columns 1–600 of 12,000.")

    def test_small_run_is_not_a_preview(self):
        p, d = core.example_paths("one_broken_gene")
        info = {}
        frag = viewer.view(p, None, None, d, page=False, info=info)
        self.assertFalse(info["preview"])
        self.assertNotIn("Preview:", frag)
        self.assertNotIn('class="overview"', frag)

    def test_large_run_preview_overview_and_links(self):
        tmp = tempfile.mkdtemp(prefix="coal_test_prev_")
        self.addCleanup(shutil.rmtree, tmp, True)
        pep, cod, rep = _synthetic(tmp, 3000, 1500, problem_every=60)
        info = {}
        frag = viewer.view(pep, cod, rep, None, page=False, max_seqs=40, max_cols=400, max_cells=16000, info=info)
        self.assertTrue(info["preview"])
        self.assertIn("Preview: 40 of 3,000 sequences, problems first", frag)
        self.assertIn('class="overview"', frag)
        self.assertLessEqual(frag.count("<path"), 8)  # one path per series, whatever the size
        self.assertNotIn("<rect", frag.replace('<rect class="ov-box"', ""))
        self.assertIn("Jalview", frag)
        c0, c1 = info["first"], info["last"]
        # every link in the fragment points at a cell that exists
        cells = set(re.findall(r'<td[^>]* id="([^"]+)"', frag))
        for target in re.findall(r'href="#(coalv-[^"]+-c\d+)"', frag):
            self.assertIn(target, cells)
        self.assertIn("(outside the drawn columns", frag)
        table = core.health_table_html(rep, links=True, preview=info)
        for target in re.findall(r'href="#(coalv-[^"]+-c\d+)"', table):
            self.assertIn(target, cells)
        self.assertIn("Centre the preview on column", table)
        # centred elsewhere
        info2 = {}
        frag2 = viewer.view(pep, cod, rep, None, page=False, max_seqs=40, max_cols=400, max_cells=16000,
                            centre=1200, info=info2)
        self.assertEqual((info2["first"], info2["last"]), (1000, 1399))
        self.assertIn("columns 1,000–1,399 of 1,500", frag2)
        self.assertNotEqual((c0, c1), (info2["first"], info2["last"]))
        # centred: the sequences with a problem in those columns are drawn first
        first_row = re.search(r'<tr id="coalv-([^"]+)"', frag2).group(1)
        col = next(m["column"] for sq in rep["sequences"] if sq["id"] == first_row for m in sq["messages"])
        self.assertTrue(1000 <= col <= 1399, (first_row, col))

    def test_rendering_is_fast_and_small_at_100000(self):
        tmp = tempfile.mkdtemp(prefix="coal_test_big_")
        self.addCleanup(shutil.rmtree, tmp, True)
        for n in (20_000, 100_000):
            pep, cod, rep = _synthetic(tmp, n, 300)
            info = {}
            t = time.perf_counter()
            frag, page = viewer.render_run(pep, cod, rep, None, 1, tmp, limits=settings.preview_limits("standard"),
                                           info=info)
            table = core.health_table_html(rep, links=True, preview=info)
            elapsed = time.perf_counter() - t
            self.assertLess(elapsed, 15.0, n)  # ~1 s here; generous for slow CI machines
            self.assertLess(len(frag.encode("utf-8")), 3_000_000, n)
            self.assertLess(os.path.getsize(page), 3_000_000, n)
            self.assertLess(len(table), 500_000, n)
            self.assertIn(f"of {n:,} sequences", frag)


class TestOverviewStrips(unittest.TestCase):
    def test_gap_profile(self):
        rows = ["A-C-", "AC--", "----", "ACGT"]
        gaps, n = overview.gap_profile(rows, 4)
        self.assertEqual(n, 4)
        self.assertEqual(gaps, [0.25, 0.5, 0.5, 0.75])
        gaps, n = overview.gap_profile(["A-"] * 600, 3)  # several 255-row batches; the short end counts as gap
        self.assertEqual(gaps, [0.0, 1.0, 1.0])

    def test_problem_columns(self):
        rep = {"sequences": [{"messages": [{"column": 2}, {"column": 2}, {"column": 9}]},
                             {"messages": [{"column": 2}, {"column": "x"}]}]}
        self.assertEqual(overview.problem_columns(rep, 3), [0, 2, 0])

    def test_strips_are_bounded(self):
        st = ["ok"] * 100_000
        st[5] = "failed"
        svg = overview.sequences_strip(st, {0, 50_000})
        self.assertLessEqual(svg.count("<path"), 3)
        self.assertIn('viewBox="0 0 400 ', svg)
        self.assertIn("1 failed", svg)


if __name__ == "__main__":
    unittest.main()
