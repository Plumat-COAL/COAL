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


class TestFinding4NoCodonFileStatus(unittest.TestCase):
    def test_mismatches_without_events_are_warnings(self):
        p, d = core.example_paths("mito")
        rows, _, _ = viewer.build_rows(viewer.read_alignment(p), None, None, viewer.read_alignment(d), 1)
        danio = next(r for r in rows if r["id"] == "Danio_COX1")
        self.assertEqual(danio["events"], [])
        self.assertTrue(any(f == "mismatch" for _, _, f in danio["cells"]))
        self.assertEqual(danio["status"], "warning")
        page = viewer.view(p, None, None, d)
        self.assertNotRegex(page, r"✔ OK</strong> <a[^>]*><code>Danio_COX1</code></a>: [^<]*Jump to")


class TestFinding5Iupac(unittest.TestCase):
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


class TestFinding6NoPairingSentenceWhenNothingPaired(unittest.TestCase):
    def test_empty_protein_file(self):
        t = _tmp(self)
        p = _write(os.path.join(t, "empty.fasta"), "")
        d = core.example_paths("clean")[1]
        r = core.run_coal(p, d, out_dir=t, probe_codes=False)
        self.assertEqual(r.exit_code, 1)
        text = core.explain(r.report)
        self.assertNotIn("paired", text)
        self.assertIsNone(r.report["summary"]["pairing"])

    def test_dna_given_as_protein(self):
        d = core.example_paths("clean")[1]
        r = core.run_coal(d, d, out_dir=_tmp(self), probe_codes=False)
        self.assertEqual(r.exit_code, 1)
        self.assertNotIn("paired", core.explain(r.report))

    def test_successful_run_still_states_pairing(self):
        p, d = core.example_paths("clean")
        r = core.run_coal(p, d, out_dir=_tmp(self), probe_codes=False)
        self.assertIn("paired by matching their IDs exactly", core.explain(r.report))


class TestFinding7IdHint(unittest.TestCase):
    def _files(self, t, pid_suffix, did_suffix):
        p, d = core.example_paths("clean")
        prot = viewer.read_alignment(p)
        dna = dict(viewer.read_alignment(d))
        pp = _write(os.path.join(t, "p.fasta"), "".join(f">{i}{pid_suffix}\n{s}\n" for i, s in prot))
        dd = _write(os.path.join(t, "d.fasta"), "".join(f">{i}{did_suffix}\n{dna[i]}\n" for i, _ in prot))
        return pp, dd

    def test_example_shows_the_real_dna_id(self):
        t = _tmp(self)
        pp, dd = self._files(t, "", ".1")
        h = core.id_hints(pp, dd)
        self.assertEqual(len(h), 1)
        self.assertIn("protein 'gene01' vs DNA 'gene01.1'", h[0]["text"])

    def test_successful_order_pairing_is_not_called_wrong(self):
        t = _tmp(self)
        pp, dd = self._files(t, ".1", "")
        r = core.run_coal(pp, dd, out_dir=t, probe_codes=False)
        self.assertEqual((r.exit_code, r.report["summary"]["ok"]), (0, 8))
        texts = [h["text"] for h in r.hints if h["code"] == "id_mismatch"]
        self.assertEqual(len(texts), 1)
        self.assertNotIn("probably wrong", texts[0])
        self.assertIn("every pair matched", texts[0])

    def test_failed_run_keeps_the_warning(self):
        p, d = core.example_paths("messy_ids")
        r = core.run_coal(p, d, out_dir=_tmp(self), probe_codes=False)
        self.assertNotEqual(r.exit_code, 0)
        texts = [h["text"] for h in r.hints if h["code"] == "id_mismatch"]
        self.assertEqual(len(texts), 1)
        self.assertIn("probably wrong", texts[0])


def _load(name, path):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestFinding9And10DiamondPipeline(unittest.TestCase):
    def setUp(self):
        self.dp = _load("diamond_pipeline_r", os.path.join(ROOT, "demo", "diamond_pipeline.py"))

    def test_every_pair_skipped_is_a_clear_error(self):
        import subprocess
        import sys
        sys.path.insert(0, os.path.join(ROOT, "tests_app"))
        try:
            from test_diamond_pipeline import make_input
        finally:
            sys.path.pop(0)
        t = _tmp(self)
        p = make_input(t)
        r = subprocess.run([sys.executable, os.path.join(ROOT, "demo", "diamond_pipeline.py"),
                            "--a-protein", p["ap"], "--a-cds", p["ac"], "--b-protein", p["bp"], "--b-cds", p["bc"],
                            "--out", os.path.join(t, "o"), "--diamond",
                            os.path.join(ROOT, "tests_app", "diamond_stub.py"), "--threads", "1",
                            "--max-length", "10"], capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertNotIn("Traceback", r.stderr)
        self.assertIn("--max-length", r.stderr)

    def test_alignment_memory_is_bytes_not_python_ints(self):
        import random
        import tracemalloc
        rng = random.Random(1)
        a = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(600))
        b = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(600))
        tracemalloc.start()
        ra, rb = self.dp.needleman_wunsch(a, b)
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        self.assertEqual((ra.replace("-", ""), rb.replace("-", "")), (a, b))
        self.assertLess(peak, 3_000_000)  # 3 x 601 x 601 bytes of traceback ~ 1.1 MB

    def test_same_alignment_as_full_matrix_reference(self):
        import random
        rng = random.Random(5)
        for _ in range(20):
            a = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(rng.randint(1, 40)))
            b = list(a)
            for _ in range(rng.randint(0, 8)):
                k = rng.randrange(len(b) + 1)
                if rng.random() < 0.5 and b:
                    del b[min(k, len(b) - 1)]
                else:
                    b.insert(k, rng.choice("ACDEFGHIKLMNPQRSTVWY"))
            b = "".join(b) or "M"
            got = self.dp.needleman_wunsch(a, b)
            self.assertEqual(_nw_score(self.dp, *got), _nw_best(self.dp, a, b))

    def test_default_max_length_is_2000(self):
        import subprocess
        import sys
        r = subprocess.run([sys.executable, os.path.join(ROOT, "demo", "diamond_pipeline.py"), "--help"],
                           capture_output=True, text=True, timeout=60)
        self.assertIn("default 2000", " ".join(r.stdout.split()))


def _nw_score(dp, ra, rb, go=11, ge=1):
    """Score of a gapped pair under BLOSUM62 + affine gaps."""
    score, state = 0, None
    for x, y in zip(ra, rb):
        if x != "-" and y != "-":
            score += dp._score(x, y)
            state = None
        else:
            st = "a" if y == "-" else "b"
            score -= ge if state == st else go
            state = st
    return score


def _nw_best(dp, s, t, go=11, ge=1):
    """Optimal score by a plain full-matrix Gotoh recursion (reference)."""
    NEG = -10 ** 9
    n, m = len(s), len(t)
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = -go - (i - 1) * ge
    for j in range(1, m + 1):
        Y[0][j] = -go - (j - 1) * ge
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            M[i][j] = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]) + dp._score(s[i - 1], t[j - 1])
            X[i][j] = max(M[i - 1][j] - go, X[i - 1][j] - ge, Y[i - 1][j] - go)
            Y[i][j] = max(M[i][j - 1] - go, X[i][j - 1] - go, Y[i][j - 1] - ge)
    return max(M[n][m], X[n][m], Y[n][m])


class TestFinding18Encodings(unittest.TestCase):
    FILES = ["demo/diamond_pipeline.py", "demo/race.py", "demo/make_mcp_demo.py", "tests_app/diamond_stub.py",
             "src/coal_app/core.py", "src/coal_app/viewer.py", "src/coal_app/tools.py", "src/coal_app/ui.py",
             "src/coal_compat/cli.py", "space/build_space.py", "examples/make_examples.py"]

    def test_text_files_and_pipes_name_utf8(self):
        import ast
        bad = []
        for rel in self.FILES:
            with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
                tree = ast.parse(fh.read())
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                fn = node.func
                name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
                kw = {k.arg for k in node.keywords}
                if name == "open":
                    mode = node.args[1].value if len(node.args) > 1 and isinstance(node.args[1], ast.Constant) \
                        else next((k.value.value for k in node.keywords if k.arg == "mode"), "r")
                    if "b" not in str(mode) and "encoding" not in kw:
                        bad.append(f"{rel}:{node.lineno} open")
                elif name == "run" and "text" in kw and "encoding" not in kw:
                    bad.append(f"{rel}:{node.lineno} subprocess.run(text=True)")
        self.assertEqual(bad, [])


if __name__ == "__main__":
    unittest.main()
