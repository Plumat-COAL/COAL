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


class _RunsDir(unittest.TestCase):
    """Point core.RUNS_DIR at a private folder for the test (and clean it up)."""

    def setUp(self):
        self._old_runs = core.RUNS_DIR
        core.RUNS_DIR = _tmp(self)
        self.addCleanup(setattr, core, "RUNS_DIR", self._old_runs)


class TestFinding11HostedSandbox(_RunsDir):
    def _hosted(self, on):
        old = os.environ.get("COAL_HOSTED")
        if on:
            os.environ["COAL_HOSTED"] = "1"
        else:
            os.environ.pop("COAL_HOSTED", None)
        self.addCleanup(lambda: os.environ.__setitem__("COAL_HOSTED", old) if old is not None
                        else os.environ.pop("COAL_HOSTED", None))

    def test_hosted_tools_refuse_paths_outside_the_sandbox(self):
        from coal_app import tools
        self._hosted(True)
        t = _tmp(self)
        secret = _write(os.path.join(t, "secret.fasta"), ">s1\nMAKL\n")
        dna = _write(os.path.join(t, "dna.fasta"), ">s1\nATGGCTAAACTG\n")
        r = tools.codon_align(secret, dna)
        self.assertFalse(r["ok"])
        self.assertIn("not allowed", r["error"])
        self.assertNotIn("output_path", {k for k, v in r.items() if v})
        for bad in (secret, "/etc/passwd", os.path.join(core.EXAMPLES_DIR, "..", "pyproject.toml")):
            self.assertIn("not allowed", tools.explain_report(bad))
        # a symlink inside the sandbox that points outside is refused too
        link = os.path.join(core.RUNS_DIR, "link.fasta")
        try:
            os.symlink(secret, link)
        except (OSError, NotImplementedError):
            self.skipTest("symlinks not available")
        self.assertIn("not allowed", tools.codon_align(link, dna)["error"])

    def test_hosted_tools_still_work_on_examples_and_their_own_reports(self):
        from coal_app import tools
        self._hosted(True)
        p, d = core.example_paths("clean")
        r = tools.codon_align(p, d)
        self.assertTrue(r["ok"], r)
        self.assertIn("Nothing went wrong", tools.explain_report(r["report_path"]))

    def test_local_tools_accept_any_path(self):
        from coal_app import tools
        self._hosted(False)
        t = _tmp(self)
        p = _write(os.path.join(t, "p.fasta"), ">s1\nMAKL\n")
        d = _write(os.path.join(t, "d.fasta"), ">s1\nATGGCTAAACTG\n")
        self.assertTrue(tools.codon_align(p, d)["ok"])

    def test_intro_depends_on_where_it_runs(self):
        self._hosted(False)
        self.assertIn("Nothing leaves your computer", core.intro_text())
        self._hosted(True)
        text = core.intro_text()
        self.assertNotIn("Nothing leaves your computer", text)
        self.assertIn("Hugging Face", text)

    def test_space_app_marks_itself_hosted_and_turns_mcp_off(self):
        with open(os.path.join(ROOT, "src", "coal_app", "space_app.py"), encoding="utf-8") as fh:
            src = fh.read()
        self.assertIn('os.environ["COAL_HOSTED"] = "1"', src)
        self.assertIn("mcp=False", src)
        with open(os.path.join(ROOT, "space", "README.md"), encoding="utf-8") as fh:
            readme = fh.read()
        self.assertNotIn("Space as an MCP server", readme)


class TestFinding12TempDirs(_RunsDir):
    def test_probe_outputs_are_removed(self):
        p, d = core.example_paths("mito")
        t = _tmp(self)
        r = core.run_coal(p, d, out_dir=t)
        self.assertEqual([h["suggested_codontable"] for h in r.hints if h["code"] == "genetic_code"], [2])
        self.assertEqual([f for f in os.listdir(t) if f.startswith("probe")], [])

    def test_default_out_dir_is_under_runs_dir_not_tmp_root(self):
        import glob
        before = set(glob.glob(os.path.join(tempfile.gettempdir(), "coal_*")))
        p, d = core.example_paths("clean")
        r = core.run_coal(p, d, probe_codes=False)
        self.assertTrue(os.path.realpath(r.output_path).startswith(os.path.realpath(core.RUNS_DIR) + os.sep))
        self.assertEqual(set(glob.glob(os.path.join(tempfile.gettempdir(), "coal_*"))) - before, set())

    def test_prune_runs_by_age(self):
        import time
        old = os.path.join(core.RUNS_DIR, "run_old")
        new = os.path.join(core.RUNS_DIR, "run_new")
        for x in (old, new):
            os.makedirs(x)
            _write(os.path.join(x, "f.txt"), "x")
            _write(os.path.join(x, core.RUN_MARKER), "")  # R2-A8
        past = time.time() - 25 * 3600
        os.utime(old, (past, past))
        os.utime(os.path.join(old, "f.txt"), (past, past))
        removed = core.prune_runs(max_age_hours=24)
        self.assertEqual(removed, 1)
        self.assertFalse(os.path.exists(old))
        self.assertTrue(os.path.exists(new))


class TestFinding13ProbeCost(_RunsDir):
    def test_probes_are_small_and_time_limited(self):
        calls = []
        real = core._run_engine

        def spy(protein, dna, out, fmt, table, nogap, nomismatch, keep_going, report_path, timeout=None):
            calls.append((protein, len(core.read_ids(protein)), timeout))
            return real(protein, dna, out, fmt, table, nogap, nomismatch, keep_going, report_path, timeout=timeout)

        core._run_engine = spy
        self.addCleanup(setattr, core, "_run_engine", real)
        # 250 sequences, every one with a TGA (W in code 2, stop in code 1): more than the 200-sequence sample
        t = _tmp(self)
        prot, dna = [], []
        for i in range(250):
            prot.append(f">s{i}\nMKWL\n")
            dna.append(f">s{i}\nATGAAATGACTG\n")
        p = _write(os.path.join(t, "p.fasta"), "".join(prot))
        d = _write(os.path.join(t, "d.fasta"), "".join(dna))
        r = core.run_coal(p, d, out_dir=t)
        probes = calls[1:]
        self.assertTrue(probes)
        for protein, n, timeout in probes:
            self.assertNotEqual(protein, p)  # a sample file, not the whole input
            self.assertLessEqual(n, 200)
            self.assertGreater(n, 0)
            self.assertIsNotNone(timeout)
            self.assertLessEqual(timeout, core.PROBE_TIMEOUT)
        self.assertLessEqual(core.PROBE_TIMEOUT, 30)
        hint = [h for h in r.hints if h["code"] == "genetic_code"]
        self.assertEqual(hint[0]["suggested_codontable"], 2)
        self.assertIn("200 of 200 sequences checked (a sample of the 250", hint[0]["text"])
        self.assertEqual([f for f in os.listdir(t) if f.startswith("probe")], [])


class TestFinding14MalformedReports(unittest.TestCase):
    def _report(self, obj=None, text=None, bom=False):
        import json
        t = _tmp(self)
        p = os.path.join(t, "r.json")
        data = text if text is not None else json.dumps(obj)
        with open(p, "w", encoding="utf-8-sig" if bom else "utf-8") as fh:
            fh.write(data)
        return p

    def test_wrong_types_give_a_clear_message(self):
        from coal_app import tools
        for obj in ({"schema": "coal.report/1", "summary": []},
                    {"schema": "coal.report/1", "output": "x"},
                    {"schema": "coal.report/1", "sequences": [1, 2]},
                    {"schema": "coal.report/1", "sequences": [{"id": "a", "messages": "oops"}]},
                    {"schema": "coal.report/1", "hints": "x"},
                    {"schema": "coal.report/1", "errors": [3]},
                    {"schema": "coal.report/1", "options": []},
                    {"schema": "coal.report/1", "inputs": {}},
                    [1, 2]):
            text = tools.explain_report(self._report(obj))
            self.assertTrue(text.startswith("Could not read the report:"), (obj, text))
            self.assertIn("not a valid COAL report", text)

    def test_bom_is_accepted(self):
        import json
        from coal_app import tools
        with open(core.SAMPLE_REPORT, encoding="utf-8") as fh:
            obj = json.load(fh)
        text = tools.explain_report(self._report(obj, bom=True))
        self.assertFalse(text.startswith("Could not read"), text)

    def test_not_json(self):
        from coal_app import tools
        self.assertIn("Could not read the report", tools.explain_report(self._report(text="{nope")))


class TestFinding15SampleReport(unittest.TestCase):
    def test_sample_is_a_real_run_labelled_demo(self):
        rep = core.load_sample_report()
        self.assertEqual(rep["coal_version"], "demo")
        self.assertTrue(rep["demo"])
        self.assertIn("DEMO SAMPLE", rep["source"])
        self.assertEqual(rep["options"]["codontable"], 1)
        self.assertEqual([h["suggested_codontable"] for h in rep["hints"]], [2])
        # every mismatch is one that code 2 explains (TGA = W, ATA = M)
        from coal_app import gcode
        for s in rep["sequences"]:
            for m in s["messages"]:
                if m["code"] == "mismatch":
                    self.assertEqual(gcode.translate(m["codon"], 2), m["residue"], m)
        self.assertIn("Demo sample", core.sample_banner(rep))
        self.assertIn("not your data", core.sample_banner(rep))

    def test_sample_matches_a_fresh_run(self):
        mk = _load("make_sample_report_r", os.path.join(ROOT, "demo", "make_sample_report.py"))
        self.assertEqual(mk.build(), core.load_sample_report())


class TestFinding17DemoWording(unittest.TestCase):
    def test_location_is_attributed_to_the_app_estimate(self):
        for rel in ("demo/mcp_demo.md", "demo/make_mcp_demo.py", "pitch/PITCH.md"):
            with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
                text = fh.read()
            self.assertNotIn("says *where*: a missing base", text, rel)
            self.assertNotIn("what is wrong and where: a missing base", text, rel)
        with open(os.path.join(ROOT, "demo", "mcp_demo.md"), encoding="utf-8") as fh:
            self.assertIn("COAL's app estimates", fh.read())


class TestFinding19NoDrift(unittest.TestCase):
    def test_app_genetic_codes_equal_the_engine_tables(self):
        from coal import tables
        from coal_app import gcode
        self.assertEqual(gcode._RAW, tables._RAW)

    def test_viewer_defines_everything_before_main_block(self):
        import ast
        with open(os.path.join(ROOT, "src", "coal_app", "viewer.py"), encoding="utf-8") as fh:
            body = ast.parse(fh.read()).body
        main_at = next(i for i, n in enumerate(body) if isinstance(n, ast.If))
        self.assertEqual([n.name for n in body[main_at:] if isinstance(n, ast.FunctionDef)], [])
        self.assertEqual(main_at, len(body) - 1)

    def test_test_classes_defined_before_unittest_main(self):
        import ast
        import glob
        for path in glob.glob(os.path.join(ROOT, "tests_app", "test_*.py")):
            with open(path, encoding="utf-8") as fh:
                body = ast.parse(fh.read()).body
            ifs = [i for i, n in enumerate(body) if isinstance(n, ast.If)]
            if ifs:
                self.assertEqual(ifs[-1], len(body) - 1, path)  # finding 16: nothing after the main block


class TestFinding16AndC4C5Ci(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(ROOT, ".github", "workflows", "ci.yml"), encoding="utf-8") as fh:
            self.ci = fh.read()

    def _job(self, name):
        m = re.search(rf"^  {name}:\n(.*?)(?=^  [\w-]+:\n|\Z)", self.ci, re.S | re.M)
        self.assertIsNotNone(m, name)
        return m.group(1)

    def test_extras_installed_and_required_in_separate_jobs(self):
        gui, mcp = self._job("gui"), self._job("mcp")
        self.assertIn('pip install -e ".[gui]"', gui)
        self.assertIn("COAL_REQUIRE_EXTRAS: gui", gui)
        self.assertNotIn("[mcp]", gui.replace("gradio[mcp]", ""))
        self.assertIn('pip install -e ".[mcp]"', mcp)
        self.assertIn("COAL_REQUIRE_EXTRAS: mcp", mcp)
        self.assertNotIn("[gui]", mcp)

    def test_cpp_job_cannot_go_green_without_exercising_cpp(self):
        cpp = self._job("cpp")
        self.assertIn("not applicable", cpp.lower())
        self.assertIn("tests/test_cpp_vs_python.py", cpp)
        self.assertIn("skipped", cpp)
        self.assertNotIn("COAL_IMPL:", cpp)  # an env var nothing reads

    def test_require_switch_turns_skips_into_failures(self):
        import subprocess
        import sys
        env = dict(os.environ, COAL_REQUIRE_EXTRAS="gui")
        r = subprocess.run([sys.executable, "-m", "unittest", "test_gui"], capture_output=True, text=True,
                           encoding="utf-8", cwd=os.path.join(ROOT, "tests_app"), env=env, timeout=300)
        try:
            import gradio  # noqa: F401
            has = True
        except ImportError:
            has = False
        if has:
            self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        else:
            self.assertNotEqual(r.returncode, 0)
            self.assertNotIn("skipped", r.stderr)


class TestC2CompatClaims(unittest.TestCase):
    def test_docstring_and_html_error_match_what_works(self):  # DEC-4: clustal is supported now
        from coal_compat import cli
        doc = " ".join(cli.__doc__.split())
        self.assertIn("An unchanged PAL2NAL command line works", doc)
        self.assertIn("default output format is clustal", doc)
        with self.assertRaises(cli.UsageError) as cm:
            cli.translate(["p", "n", "-html"])
        self.assertIn("not supported", str(cm.exception))

    def test_nostderr_is_accepted_as_quiet(self):
        from coal_compat import cli
        self.assertEqual(cli.translate(["p", "n", "-output", "fasta", "-nostderr"]),
                         ["-output", "fasta", "-q", "--", "p", "n"])  # DEC-4

    def test_qa_says_which_pipelines_can_switch(self):
        with open(os.path.join(ROOT, "pitch", "QA.md"), encoding="utf-8") as fh:
            qa = " ".join(fh.read().split())
        self.assertNotIn("a pipeline can switch without being rewritten", qa)
        self.assertIn("an unchanged PAL2NAL command line works", qa)  # DEC-4: was "-output fasta"


class TestC10C11Skill(unittest.TestCase):
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
        for later in ("--keep-going", "--report", "normalised"):  # DEC-4: -output paml is current now
            self.assertNotIn(later, current, later)
        self.assertIn("usage error", current)
        self.assertIn("coal --help", current)

    def test_extras_in_separate_environments(self):
        self.assertIn("separate environments", self.text)


if __name__ == "__main__":
    unittest.main()
