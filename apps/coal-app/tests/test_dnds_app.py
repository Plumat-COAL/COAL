"""Ka/Ks and the dN/dS pipeline in the app: the logic (coal_app.dnds), the GUI's "Ka/Ks" and "dN/dS pipeline" tabs
(against the Gradio stand-in) and the MCP tools dnds_pairs / pipeline_pairs as plain functions.

The CI's app jobs install the pure-Python engine (no coal-cpp, no DIAMOND), so most tests run COAL through
tests/coal_stub.py, which answers `coal --cpp-path`, `coal dnds` and `coal pipeline pairs` with the same files and
progress lines (made-up numbers). TestRealEngine runs the real ones where coal-cpp (and DIAMOND) are found.
"""
import csv
import json
import os
import sys
import threading
import time
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from helpers import HERE, ROOT, RunsDirCase, UiCase, _Comp, _write  # noqa: E402

from coal_app import core, dnds, tools  # noqa: E402

STUB = os.path.join(HERE, "coal_stub.py")
EX = os.path.join(ROOT, "examples", "dnds")
EX_CODON = os.path.join(EX, "codon.fasta")
EX_A = os.path.join(EX, "species_a_cds.fasta")
EX_B = os.path.join(EX, "species_b_cds.fasta")

FASTA4 = ">a\nATGAAACCCGGG\n>b\nATGAAGCCAGGT\n>c\nATGAAACCTGGA\n>d\nATGAGACCCGGC\n"
AXT2 = "# made by hand\na&b\nATGAAACCC\nATGAAGCCA\n\nc&d\nATGAAA\nATGAAG\n"


def _env(test, **kv):
    """Set (value) or remove (None) environment variables for one test."""
    for k, v in kv.items():
        old = os.environ.get(k)
        test.addCleanup(lambda k=k, old=old: os.environ.__setitem__(k, old) if old is not None
                        else os.environ.pop(k, None))
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v


def use_stub(test, diamond="/stub/diamond"):
    """Run `coal` as tests/coal_stub.py, with DIAMOND 'found' at `diamond` (None: not found)."""
    core.detect_features()  # cached from the real engine first (the stub hands --help to it anyway)
    orig_cmd, orig_dia = core._coal_cmd, dnds.diamond_path
    core._coal_cmd = lambda: [sys.executable, STUB]
    dnds.diamond_path = lambda: diamond
    test.addCleanup(setattr, core, "_coal_cmd", orig_cmd)
    test.addCleanup(setattr, dnds, "diamond_path", orig_dia)
    _env(test, COAL_STUB_NO_CPP=None, COAL_STUB_SLEEP=None, COAL_STUB_FAIL=None, COAL_STUB_PIDFILE=None,
         COAL_HOSTED=None)


def _cds(d, name, seqs):
    return _write(os.path.join(d, name), "".join(f">{n}\n{s}\n" for n, s in seqs))


def U(u, key="value"):
    """A gr.update's field, from real Gradio (a dict) or the stand-in (a _Comp)."""
    return u[key] if isinstance(u, dict) else u._kw.get(key)


def two_species(d):
    a = _cds(d, "a.fasta", [("a1", "ATGAAACCCGGGTAA"), ("a2", "ATGTTTCCCTAA"), ("a3", "ATGCCCAAATAA")])
    b = _cds(d, "b.fasta", [("b1", "ATGAAGCCCGGTTAA"), ("b2", "ATGTTCCCGTAA"), ("b3", "ATGCCAAAGTAA")])
    return a, b


# ---------------------------------------------------------------- logic

class TestInputsAndLimits(RunsDirCase):
    def test_methods_are_the_engines(self):
        from coal import dnds_pipeline
        self.assertEqual(dnds.METHODS, tuple(dnds_pipeline.METHODS))
        self.assertEqual(set(dnds.METHOD_NAMES), set(dnds.METHODS))
        self.assertIn("YN", dnds.DEFAULT_METHODS)
        # engine order, any case, the C95 alias, duplicates once
        self.assertEqual(dnds.validate_methods(["yn", "NG", "c95", "YN"], 1), ["NG", "YN", "Comeron"])
        self.assertEqual(dnds.validate_methods("MYN", 11), ["MYN"])
        for bad, code, words in (([], 1, "at least one"), (["GY"], 1, "unknown Ka/Ks method"),
                                 (["Comeron"], 2, "standard genetic code"), (["Li", "Comeron"], 11, "are defined"),
                                 (["YN"], 7, "unknown genetic code")):
            with self.subTest(bad=bad, code=code), self.assertRaises(ValueError) as cm:
                dnds.validate_methods(bad, code)
            self.assertIn(words, str(cm.exception))

    def test_records_pairs_and_format(self):
        d = self.tmp()
        fa, axt = _write(os.path.join(d, "x.fa"), FASTA4), _write(os.path.join(d, "x.axt"), AXT2)
        self.assertEqual(dnds.count_records(fa), ("fasta", 4))
        self.assertEqual(dnds.count_records(axt), ("axt", 2))  # the '#' line before a record is not one
        self.assertEqual(dnds.n_pairs("fasta", 4, "all"), 6)
        self.assertEqual(dnds.n_pairs("fasta", 5, "consecutive"), 2)
        self.assertEqual(dnds.n_pairs("axt", 2, "all"), 2)
        self.assertEqual(dnds.check_kaks_input(fa, "all", 2), ("fasta", 6))
        for text, words in (("", "empty"), ("\n\n", "empty"), (">only\nATG\n", "at least two")):
            with self.subTest(text=text), self.assertRaises(ValueError) as cm:
                dnds.check_kaks_input(_write(os.path.join(d, "e.fa"), text), "all", 1)
            self.assertIn(words, str(cm.exception))

    def test_limits_are_said_and_smaller_when_hosted(self):
        d = self.tmp()
        fa = _write(os.path.join(d, "x.fa"), FASTA4)
        old = dnds.KAKS_MAX_BYTES, dnds.KAKS_MAX_RESULTS
        self.addCleanup(lambda: (setattr(dnds, "KAKS_MAX_BYTES", old[0]), setattr(dnds, "KAKS_MAX_RESULTS", old[1])))
        dnds.KAKS_MAX_BYTES = 10
        with self.assertRaises(ValueError) as cm:
            dnds.check_kaks_input(fa, "all", 1)
        self.assertIn("up to 10 bytes", str(cm.exception))
        self.assertIn("coal dnds", str(cm.exception))
        dnds.KAKS_MAX_BYTES, dnds.KAKS_MAX_RESULTS = old[0], 11
        with self.assertRaises(ValueError) as cm:
            dnds.check_kaks_input(fa, "all", 2)  # 6 pairs x 2 methods = 12 > 11
        self.assertIn("12 results", str(cm.exception))
        self.assertEqual(dnds.check_kaks_input(fa, "consecutive", 2), ("fasta", 2))
        _env(self, COAL_HOSTED="1")
        self.assertEqual(dnds.kaks_limits(), (dnds.KAKS_MAX_BYTES_HOSTED, dnds.KAKS_MAX_RESULTS_HOSTED))
        self.assertLess(dnds.KAKS_MAX_BYTES_HOSTED, old[0])
        self.assertLess(dnds.KAKS_MAX_RESULTS_HOSTED, old[1])
        with self.assertRaises(ValueError):
            dnds.check_cds_input(_write(os.path.join(d, "z.fa"), ""), "first")
        with self.assertRaises(ValueError) as cm:
            dnds.check_cds_input(fa, "first", max_bytes=5)
        self.assertIn("coal pipeline pairs", str(cm.exception))

    def test_pasted_text(self):
        d = self.tmp()
        self.assertTrue(dnds.write_pasted(FASTA4, d).endswith("pasted.fasta"))
        p = dnds.write_pasted(AXT2.rstrip("\n"), d)
        self.assertTrue(p.endswith("pasted.axt"))
        with open(p, encoding="utf-8") as fh:
            self.assertTrue(fh.read().endswith("\n"))
        old = dnds.KAKS_MAX_PASTE
        self.addCleanup(setattr, dnds, "KAKS_MAX_PASTE", old)
        dnds.KAKS_MAX_PASTE = 10
        with self.assertRaises(ValueError) as cm:
            dnds.write_pasted(FASTA4, d)
        self.assertIn("upload the file instead", str(cm.exception))

    def test_every_status_is_explained(self):
        self.assertEqual(dnds.STATUS_ORDER, ("ok", "warning", "not_converged", "saturated", "no_data", "error"))
        for s in dnds.STATUS_ORDER:
            self.assertTrue(dnds.STATUS_TEXT[s])
            self.assertIn(s, dnds.STATUS_ICON)
        counts = dnds.status_counts([{"status": "ok"}, {"status": "saturated"}, {"status": "ok"}])
        self.assertEqual((counts["total"], counts["ok"], counts["saturated"], counts["error"]), (3, 2, 1, 0))
        self.assertEqual([h["status"] for h in dnds.status_help(counts)], ["ok", "saturated"])
        self.assertEqual(len(dnds.status_help(counts, all_statuses=True)), 6)
        md = dnds.status_help_md()
        for s in dnds.STATUS_ORDER:
            self.assertIn(f"**{s}**", md)
        self.assertIn("(2)", dnds.status_help_md(counts))


class TestKaksRun(RunsDirCase):
    def setUp(self):
        super().setUp()
        use_stub(self)

    def test_fasta_pairs_and_methods(self):
        d = self.tmp()
        fa = _write(os.path.join(d, "x.fa"), FASTA4)
        r = dnds.run_kaks(fa, d, ["NG", "YN"], 1, "consecutive")
        self.assertTrue(r.ok)
        self.assertEqual(r.exit_code, 0)
        self.assertEqual((r.counts["total"], r.counts["ok"]), (4, 4))  # 2 pairs x 2 methods
        self.assertEqual([x["method"] for x in r.rows[:2]], ["NG", "YN"])
        self.assertTrue(os.path.isfile(r.tsv_path))
        self.assertEqual(r.command, "coal dnds x.fa -o x.dnds.tsv -m NG -m YN -c 1 --pairs consecutive")
        r = dnds.run_kaks(fa, d, ["YN"], 1, "all", min_codons=50, max_ds=2)
        self.assertEqual(r.counts["total"], 6)
        self.assertIn("--min-codons 50 --max-ds 2", r.command)

    def test_axt_and_unaligned_input(self):
        d = self.tmp()
        r = dnds.run_kaks(_write(os.path.join(d, "x.axt"), AXT2), d, ["YN"])
        self.assertEqual((r.fmt, r.counts["ok"]), ("axt", 2))
        self.assertNotIn("--pairs", r.command)  # AXT is read pair by pair
        bad = _write(os.path.join(d, "bad.fa"), ">a\nATGAAA\n>b\nATGAAACCC\n")
        r = dnds.run_kaks(bad, d, ["YN"])
        self.assertTrue(r.ok)  # the table is written ...
        self.assertEqual(r.exit_code, 2)  # ... but no pair has an estimate
        self.assertIn("no pair has an estimate", r.error)
        html = dnds.kaks_banner(r)
        self.assertIn("No pair has an estimate", html)
        self.assertIn("coal-failed", html)

    def test_refused_before_the_engine(self):
        d = self.tmp()
        fa = _write(os.path.join(d, "x.fa"), FASTA4)
        with self.assertRaises(ValueError):
            dnds.run_kaks(fa, d, ["YN"], 1, "every")
        with self.assertRaises(ValueError):
            dnds.run_kaks(fa, d, ["YN"], 1, "all", min_codons=-1)
        _env(self, COAL_STUB_NO_CPP="1")
        self.assertIsNone(dnds.cpp_path())
        with self.assertRaises(ValueError) as cm:
            dnds.run_kaks(fa, d, ["YN"])
        self.assertIn("coal-cpp", str(cm.exception))
        self.assertFalse(dnds.tool_status().pipeline_ready)


class TestPipelineRun(RunsDirCase):
    def setUp(self):
        super().setUp()
        use_stub(self)

    def test_run_outputs_zip_and_cleanup(self):
        d = self.tmp()
        a, b = two_species(d)
        out = core.new_run_dir()
        cmd = dnds.pipeline_command(a, b, out, ["YN", "NG"], 1, longest_isoform=True, diamond="/stub/diamond",
                                    threads=2)
        self.assertEqual(cmd[2:4], ["pipeline", "pairs"])
        for flag in ("--longest-isoform", "--ortho", "--diamond", "-j"):
            self.assertIn(flag, cmd)
        lines = []
        r = dnds.run_pipeline(cmd, out, on_line=lines.append, hide=[a, b])
        self.assertEqual((r.status, r.exit_code), ("completed", 0))
        self.assertTrue(lines[0].startswith("coal pipeline: translated 3 query"))
        self.assertEqual((r.counts["total"], r.counts["ok"]), (6, 6))
        self.assertTrue(all(os.path.isabs(p) and os.path.isfile(p) for p in r.outputs.values()), r.outputs)
        self.assertEqual(r.report_path, os.path.join(out, "report.json"))
        self.assertFalse(os.path.exists(os.path.join(out, "work")))  # intermediate files removed
        with zipfile.ZipFile(r.zip_path) as z:
            self.assertEqual(sorted(z.namelist()), ["dnds.tsv", "failed.tsv", "genes.tsv", "report.json"])
        self.assertNotIn(d, r.command)  # file names only
        self.assertIn("--query a.fasta --subject b.fasta", r.command)
        self.assertIn("✔", dnds.pipeline_banner(r))

    def test_stages_follow_the_progress_lines(self):
        st = 0
        for line, want in (("coal pipeline: translated 3 query and 3 subject records", 1),
                           ("coal pipeline: 12 pairs (rbh)", 2), ("coal pipeline: warning: x", 2),
                           ("coal pipeline: summary: 6 records", 3)):
            st = dnds.stage_of(line, st)
            self.assertEqual(st, want, line)
        html = dnds.progress_html(1, elapsed=3)
        self.assertEqual(html.count("✔"), 1)
        self.assertIn("▶", html)
        self.assertIn("(running)", html)
        self.assertIn("✖", dnds.progress_html(2, ended=True))
        self.assertEqual(dnds.progress_html(0, done=True).count("✔"), len(dnds.STAGES))

    def test_failure_is_reported(self):
        d = self.tmp()
        a, b = two_species(d)
        _env(self, COAL_STUB_FAIL="1")
        out = core.new_run_dir()
        r = dnds.run_pipeline(dnds.pipeline_command(a, b, out), out)
        self.assertEqual((r.status, r.exit_code), ("failed", 1))
        self.assertIn("DIAMOND failed", r.error)
        self.assertIn("did not finish", dnds.pipeline_banner(r))

    def _slow_run(self, key, timeout=dnds.PIPELINE_TIMEOUT):
        d = self.tmp()
        a, b = two_species(d)
        pidfile = os.path.join(d, "pid")
        _env(self, COAL_STUB_SLEEP="30", COAL_STUB_PIDFILE=pidfile)
        out = core.new_run_dir()
        return dnds.iter_pipeline(dnds.pipeline_command(a, b, out), out, timeout, key=key), out, pidfile

    def _gone(self, pid, limit=15):
        end = time.monotonic() + limit
        while time.monotonic() < end:
            try:
                os.kill(pid, 0)
            except OSError:
                return True
            if os.name == "nt":  # os.kill(pid, 0) on Windows is not a probe
                return True
            time.sleep(0.2)
        return False

    def test_stop_ends_the_run_and_its_processes(self):
        gen, out, pidfile = self._slow_run("session-1")
        kind, line = next(gen)  # the first progress line: the run is under way
        self.assertEqual(kind, "line")
        self.assertFalse(dnds.stop_run("another-session"))
        t0 = time.monotonic()
        threading.Timer(0.3, dnds.stop_run, ("session-1",)).start()
        result = [v for k, v in gen if k == "done"][0]
        self.assertLess(time.monotonic() - t0, 20)
        self.assertEqual(result.status, "interrupted")
        self.assertIn("Stop", result.error)
        self.assertIn("stopped", dnds.pipeline_banner(result))
        self.assertFalse(os.path.exists(os.path.join(out, "work")))
        with open(pidfile) as fh:
            self.assertTrue(self._gone(int(fh.read())))
        self.assertNotIn("session-1", dnds._RUNNING)

    def test_time_limit(self):
        gen, out, pidfile = self._slow_run(None, timeout=1)
        result = [v for k, v in gen if k == "done"][0]
        self.assertEqual(result.status, "interrupted")
        self.assertIn("did not finish within 1 seconds", result.error)

    def test_closing_the_generator_stops_the_run(self):
        gen, out, pidfile = self._slow_run("session-2")
        next(gen)
        gen.close()
        with open(pidfile) as fh:
            self.assertTrue(self._gone(int(fh.read())))
        self.assertFalse(os.path.exists(os.path.join(out, "work")))


class TestHtml(unittest.TestCase):
    def test_table_escapes_and_limits(self):
        rows = [{"seq1": "<b>x</b>", "seq2": "y", "method": "YN", "status": "ok", "note": "<script>"}] * 5
        html = dnds.table_html(rows, limit=3)
        self.assertNotIn("<script>", html)
        self.assertIn("&lt;b&gt;x&lt;/b&gt;", html)
        self.assertEqual(html.count("<tr class="), 3)
        self.assertIn("The first 3 of 5 rows", html)
        self.assertIn("✔", html)  # status as icon + word
        self.assertIn("No result rows", dnds.table_html([]))

    def test_tools_line(self):
        ready = dnds.ToolStatus("/x/coal-cpp", "/x/diamond", "2.1.9")
        self.assertTrue(ready.pipeline_ready)
        self.assertIn("DIAMOND 2.1.9 found", dnds.tools_md(ready))
        none = dnds.ToolStatus(None, None)
        text = dnds.tools_md(none)
        self.assertIn("coal pipeline install diamond", text)
        self.assertIn("coal-cpp", text)
        self.assertFalse(none.kaks_ready or none.pipeline_ready)

    def test_install_command_in_the_windows_app(self):
        self.assertEqual(dnds.install_command(), "coal pipeline install diamond")
        sys.frozen = True  # the frozen COAL.exe has no `coal` on PATH: the command goes through COAL.exe
        try:
            self.assertEqual(dnds.install_command(), f'"{sys.executable}" -m coal pipeline install diamond')
        finally:
            del sys.frozen


# ---------------------------------------------------------------- the GUI (Gradio stand-in)

class TestDndsTabs(UiCase):
    def setUp(self):
        super().setUp()
        use_stub(self)

    def tmp(self):
        import shutil
        import tempfile
        d = tempfile.mkdtemp(prefix="coal_test_")
        self.addCleanup(shutil.rmtree, d, True)
        return d

    def _events(self):
        self.ui.build()
        return list(_Comp.REG)

    def test_tabs_wiring_and_no_path_boxes(self):
        events = self._events()
        import gc
        # the four tabs, in this order: the work (align, Ka/Ks, pipeline), then Settings
        with open(self.ui.__file__, encoding="utf-8") as fh:
            src = fh.read()
        at = [src.index(f'gr.Tab("{t}"') for t in ("Align", "Ka/Ks", "dN/dS pipeline", "Settings")]
        self.assertEqual(at, sorted(at))
        fns = {a[0].__name__: (c, ev, a, kw) for c, ev, a, kw in events if a and callable(a[0])}
        self.assertIn("kaks_ui", fns)
        self.assertIn("pipeline_ui", fns)
        self.assertEqual(fns["pipeline_ui"][3].get("concurrency_limit"), 1)  # one pipeline at a time
        self.assertIn("stop_pipeline_ui", fns)
        self.assertIn("pipeline_tools_ui", fns)
        self.assertEqual(len(fns["kaks_ui"][2][1]), 7)
        self.assertEqual(len(fns["pipeline_ui"][2][2]), 10)
        # the new tabs take uploads and pasted text only: no typed paths from the web page
        boxes = [c for c in gc.get_objects() if isinstance(c, _Comp) and c._name == "Textbox"
                 and str(c._kw.get("elem_id", "")).startswith(("coal-kaks", "coal-pipe"))]
        self.assertEqual(sorted(c._kw["elem_id"] for c in boxes), ["coal-kaks-text", "coal-pipe-log"])
        self.assertFalse(any("path" in str(c._kw.get("label", "")).lower() for c in boxes))
        # the API / MCP endpoints: the four tools
        for name in ("codon_align", "explain_report", "dnds_pairs", "pipeline_pairs"):
            self.assertIn(f'gr.api(tools.{name}, api_name="{name}")', src)

    def test_kaks_ui(self):
        d = self.tmp()
        fa = _write(os.path.join(d, "upload.fa"), FASTA4)
        out = self.ui.kaks_ui(fa, "", "all", ["YN"], 1, 0, 0)
        self.assertIn("All 6 results estimated", U(out[0], "value"))
        self.assertIn("File used: upload.fa", U(out[0], "value"))
        self.assertTrue(U(out[2], "visible") and os.path.isfile(U(out[2], "value")))
        self.assertTrue(U(out[2], "value").startswith(core.RUNS_DIR))
        self.assertIn("<table", U(out[1], "value"))
        self.assertIn("**ok** (6)", U(out[3], "value"))
        self.assertIn("[coal]", U(out[5], "value"))  # what to cite
        self.assertEqual(U(out[6], "value"), "coal dnds upload.fa -o upload.dnds.tsv -m YN -c 1 --pairs all")
        self.assertTrue(U(out[7], "visible"))
        # pasted text, when nothing is uploaded
        out = self.ui.kaks_ui(None, AXT2, "all", ["YN", "NG"], 1)
        self.assertIn("Pasted text", U(out[0], "value"))
        self.assertIn("4 results", U(out[0], "value"))
        # nothing given; a bad method choice
        self.assertIn("Please add an alignment", U(self.ui.kaks_ui(None, "  ", "all", ["YN"], 1)[0], "value"))
        out = self.ui.kaks_ui(fa, "", "all", ["Comeron"], 2)
        self.assertIn("Ka/Ks did not run", U(out[0], "value"))
        self.assertIn("standard genetic code", U(out[0], "value"))
        self.assertFalse(U(out[2], "visible"))
        self.assertIn("Please add", U(self.ui.kaks_ui(None, "", "all", [], 1)[0], "value"))

    def test_kaks_example(self):
        path, text, pairs, methods, code = self.ui.kaks_example_ui()
        self.assertEqual(path, EX_CODON if os.path.isfile(EX_CODON) else path)
        self.assertEqual((text, pairs, methods, code), ("", "consecutive", ["YN", "NG"], 1))
        if path:
            out = self.ui.kaks_ui(path, text, pairs, methods, code)
            self.assertIn("All 24 results estimated", U(out[0], "value"))

    def test_pipeline_ui_streams_and_offers_the_results(self):
        d = self.tmp()
        a, b = two_species(d)
        steps = list(self.ui.pipeline_ui(a, b, ["YN"], 1))
        self.assertGreaterEqual(len(steps), 4)  # start, a progress line or more, the end
        self.assertIn("Running the dN/dS pipeline", U(steps[0][0], "value"))
        self.assertIn("▶", U(steps[0][1], "value"))
        self.assertTrue(any("translated 3 query" in U(s[2]) for s in steps))
        last = steps[-1]
        self.assertIn("✔", U(last[0], "value"))
        self.assertEqual(U(last[1], "value").count("✔"), len(dnds.STAGES))
        for i, name in ((4, "dnds.tsv"), (5, "report.json"), (6, "coal_dnds_results.zip")):
            self.assertTrue(U(last[i], "visible"), name)
            self.assertEqual(os.path.basename(U(last[i], "value")), name)
        with open(U(last[5], "value"), encoding="utf-8") as fh:
            self.assertEqual(json.load(fh)["schema"], "coal.pipeline/1")
        self.assertIn("coal pipeline pairs --query a.fasta --subject b.fasta", U(last[7], "value"))
        self.assertIn("[buchfink2021]", U(last[8], "value"))

    def test_pipeline_ui_refusals(self):
        d = self.tmp()
        a, b = two_species(d)
        out = list(self.ui.pipeline_ui(a, None, ["YN"], 1))
        self.assertEqual(len(out), 1)
        self.assertIn("Please add both files", U(out[0][0], "value"))
        self.assertIn("the second species", U(out[0][0], "value"))
        out = U(list(self.ui.pipeline_ui(a, b, ["Li"], 11))[0][0])
        self.assertIn("standard genetic code", out)
        dnds.diamond_path = lambda: None
        out = U(list(self.ui.pipeline_ui(a, b, ["YN"], 1))[0][0])
        self.assertIn("A tool is missing", out)
        self.assertIn("coal pipeline install diamond", out)
        tools_md, run_btn = self.ui.pipeline_tools_ui()
        self.assertIn("coal pipeline install diamond", U(tools_md))
        self.assertFalse(U(run_btn, "interactive"))
        dnds.diamond_path = lambda: "/stub/diamond"
        self.assertTrue(U(self.ui.pipeline_tools_ui()[1], "interactive"))
        _env(self, COAL_HOSTED="1")
        out = U(list(self.ui.pipeline_ui(a, b, ["YN"], 1))[0][0])
        self.assertIn("Not on this shared server", out)

    def test_stop_button_stops_this_sessions_run(self):
        d = self.tmp()
        a, b = two_species(d)
        _env(self, COAL_STUB_SLEEP="30")

        class Req:
            session_hash = "browser-1"
        gen = self.ui.pipeline_ui(a, b, ["YN"], 1, request=Req())
        next(gen)
        next(gen)  # the first progress line
        threading.Timer(0.3, self.ui.stop_pipeline_ui, (Req(),)).start()
        last = list(gen)[-1]
        self.assertIn("The pipeline was stopped", U(last[0], "value"))
        self.assertIn("you pressed Stop", U(last[0], "value"))
        self.assertIn("✖", U(last[1], "value"))

    def test_about_says_what_the_app_does(self):
        about = self.ui.ABOUT_SCOPE
        for words in ("codon alignments", "Ka/Ks", "pipeline", "DIAMOND", "coal pipeline clusters"):
            self.assertIn(words, about)
        self.assertNotIn("not in this app", about)


# ---------------------------------------------------------------- the MCP tools as plain functions

class TestDndsTools(RunsDirCase):
    def setUp(self):
        super().setUp()
        use_stub(self)

    def test_dnds_pairs(self):
        d = self.tmp()
        recs = "".join(f">s{i}\nATGAAACC{'ACGT'[i % 4]}GGG\n" for i in range(8))
        fa = _write(os.path.join(d, "aln.fasta"), recs)
        r = tools.dnds_pairs(fa, methods=["YN", "NG"])
        self.assertEqual(set(r), set(tools.DndsPairsResult.__annotations__))
        self.assertTrue(r["completed"] and r["all_estimated"])
        self.assertEqual((r["pairs"], r["counts"]["total"], r["rows_total"]), (28, 56, 56))
        self.assertEqual(len(r["rows"]), tools.MAX_ROWS)
        self.assertEqual(r["methods"], ["NG", "YN"])
        self.assertEqual((r["flagged"], r["flagged_total"]), ([], 0))
        self.assertEqual([s["status"] for s in r["statuses"]], ["ok"])
        self.assertIn("coal", r["cite"])
        self.assertIn("Every result is estimated", r["next_step"])
        self.assertEqual(r["limits"]["max_bytes"], dnds.KAKS_MAX_BYTES)
        # save_as: a copy where the person wants it, never over a file
        dest = os.path.join(d, "out.tsv")
        r = tools.dnds_pairs(fa, pairs="consecutive", save_as=dest)
        self.assertEqual(r["output_path"], dest)
        with open(dest, encoding="utf-8") as fh:
            self.assertEqual(len(list(csv.DictReader(fh, delimiter="\t"))), 4)
        with self.assertRaises(tools.CoalToolError):
            tools.dnds_pairs(fa, save_as=dest)
        with self.assertRaises(tools.CoalToolError):
            tools.dnds_pairs(fa, save_as="relative.tsv")

    def test_dnds_pairs_flags_and_errors(self):
        d = self.tmp()
        bad = _write(os.path.join(d, "bad.fa"), ">a\nATGAAA\n>b\nATGAAACCC\n>c\nATGAAG\n")
        r = tools.dnds_pairs(bad)
        self.assertTrue(r["completed"])
        self.assertFalse(r["all_estimated"])
        self.assertEqual((r["counts"]["ok"], r["counts"]["error"], r["flagged_total"]), (1, 2, 2))
        self.assertEqual(r["flagged"][0]["status"], "error")
        self.assertIn("error", [s["status"] for s in r["statuses"]])
        self.assertIn("flagged", r["next_step"])
        for kw, words in (({"methods": ["GY"]}, "unknown Ka/Ks method"),
                          ({"methods": ["Comeron"], "codontable": 2}, "standard genetic code"),
                          ({"codontable": 8}, "unknown genetic code"),
                          ({"pairs": "every"}, "'all' or 'consecutive'"),
                          ({"min_codons": -3}, "negative")):
            with self.subTest(kw=kw), self.assertRaises(tools.CoalToolError) as cm:
                tools.dnds_pairs(bad, **kw)
            self.assertIn(words, str(cm.exception))
        with self.assertRaises(tools.CoalToolError):
            tools.dnds_pairs(os.path.join(d, "missing.fa"))
        old = dnds.KAKS_MAX_RESULTS
        self.addCleanup(setattr, dnds, "KAKS_MAX_RESULTS", old)
        dnds.KAKS_MAX_RESULTS = 2
        with self.assertRaises(tools.CoalToolError) as cm:
            tools.dnds_pairs(bad)
        self.assertIn("up to 2 at a time", str(cm.exception))

    def test_dnds_pairs_when_hosted(self):
        d = self.tmp()
        fa = _write(os.path.join(d, "x.fa"), FASTA4)
        _env(self, COAL_HOSTED="1")
        with self.assertRaises(tools.CoalToolError) as cm:  # outside the sandbox
            tools.dnds_pairs(fa)
        self.assertIn("not allowed on this shared server", str(cm.exception))
        inside = _write(os.path.join(core.RUNS_DIR, "x.fa"), FASTA4)
        with self.assertRaises(tools.CoalToolError):
            tools.dnds_pairs(inside, save_as=os.path.join(d, "y.tsv"))
        self.assertEqual(tools.dnds_pairs(inside)["limits"]["max_results"], dnds.KAKS_MAX_RESULTS_HOSTED)

    def test_pipeline_pairs(self):
        d = self.tmp()
        a, b = two_species(d)
        r = tools.pipeline_pairs(a, b, methods=["YN"])
        self.assertEqual(set(r), set(tools.PipelinePairsResult.__annotations__))
        self.assertTrue(r["completed"])
        self.assertEqual((r["status"], r["exit_code"], r["pair_counts"]["ok"]), ("completed", 0, 3))
        self.assertTrue(r["out_dir"].startswith(core.RUNS_DIR))
        self.assertTrue(os.path.isfile(r["outputs"]["dnds"]) and os.path.isfile(r["report_path"]))
        self.assertTrue(os.path.isfile(r["zip_path"]))
        self.assertEqual(r["diamond"], "/stub/diamond")
        self.assertIn("buchfink2021", r["cite"])
        self.assertTrue(r["progress"][0].startswith("coal pipeline: translated"))
        self.assertIn("Every pair is estimated", r["next_step"])
        # a folder the person names: made if missing; a folder with files in it is refused
        mine = os.path.join(d, "results")
        r = tools.pipeline_pairs(a, b, out_dir=mine)
        self.assertEqual(r["out_dir"], mine)
        self.assertTrue(os.path.isfile(os.path.join(mine, "dnds.tsv")))
        self.assertFalse(os.path.exists(os.path.join(mine, "work")))
        with self.assertRaises(tools.CoalToolError) as cm:
            tools.pipeline_pairs(a, b, out_dir=mine)
        self.assertIn("not empty", str(cm.exception))
        for bad in ("relative/dir", os.path.join(d, "no", "such", "parent"), a):
            with self.subTest(out_dir=bad), self.assertRaises(tools.CoalToolError):
                tools.pipeline_pairs(a, b, out_dir=bad)

    def test_pipeline_pairs_without_diamond(self):
        d = self.tmp()
        a, b = two_species(d)
        dnds.diamond_path = lambda: None
        with self.assertRaises(tools.CoalToolError) as cm:
            tools.pipeline_pairs(a, b)
        self.assertIn("coal pipeline install diamond", str(cm.exception))
        self.assertIn("pairs_path", str(cm.exception))
        pairs = _write(os.path.join(d, "pairs.tsv"), "a1\tb1\na2\tb2\n")
        r = tools.pipeline_pairs(a, b, pairs_path=pairs)  # the given pairs need no DIAMOND
        self.assertTrue(r["completed"])
        self.assertIsNone(r["diamond"])
        self.assertEqual(r["pair_counts"]["total"], 2)
        self.assertNotIn("--diamond", r["coal_command"])
        self.assertIn(f"--pairs {core._q(pairs)}", r["coal_command"])  # the agent gets a command it can run as it is

    def test_pipeline_pairs_refusals(self):
        d = self.tmp()
        a, b = two_species(d)
        _env(self, COAL_STUB_NO_CPP="1")
        with self.assertRaises(tools.CoalToolError) as cm:
            tools.pipeline_pairs(a, b)
        self.assertIn("coal-cpp", str(cm.exception))
        _env(self, COAL_STUB_NO_CPP=None)
        old = dnds.PIPELINE_MAX_BYTES_MCP
        self.addCleanup(setattr, dnds, "PIPELINE_MAX_BYTES_MCP", old)
        dnds.PIPELINE_MAX_BYTES_MCP = 10
        with self.assertRaises(tools.CoalToolError) as cm:
            tools.pipeline_pairs(a, b)
        self.assertIn("up to 10 bytes", str(cm.exception))
        dnds.PIPELINE_MAX_BYTES_MCP = old
        for kw in ({"ortho": "all"}, {"methods": []}, {"threads": -1}):
            with self.subTest(kw=kw), self.assertRaises(tools.CoalToolError):
                tools.pipeline_pairs(a, b, **kw)
        _env(self, COAL_STUB_FAIL="1")
        r = tools.pipeline_pairs(a, b)
        self.assertFalse(r["completed"])
        self.assertEqual(r["status"], "failed")
        self.assertIn("DIAMOND failed", r["error"])
        self.assertIn("did not finish", r["next_step"])
        _env(self, COAL_HOSTED="1")
        with self.assertRaises(tools.CoalToolError) as cm:
            tools.pipeline_pairs(a, b)
        self.assertIn("not available on this shared server", str(cm.exception))


# ---------------------------------------------------------------- the real engine, where it is installed

_CPP = dnds.cpp_path()
# CI's app jobs build coal-cpp and set COAL_REQUIRE_CPP, and COAL_REQUIRE_DIAMOND after installing DIAMOND: there a
# missing engine or DIAMOND fails these tests instead of skipping them (a pure-Python install still skips them).
_REQUIRE_CPP = bool(os.environ.get("COAL_REQUIRE_CPP"))
_REQUIRE_DIAMOND = bool(os.environ.get("COAL_REQUIRE_DIAMOND"))


@unittest.skipUnless((_CPP or _REQUIRE_CPP) and os.path.isfile(EX_CODON),
                     "coal-cpp not found (`coal --cpp-path`): pure-Python install")
class TestRealEngine(RunsDirCase):
    def setUp(self):
        super().setUp()
        if not _CPP:
            self.fail("coal-cpp not found (`coal --cpp-path`) and COAL_REQUIRE_CPP is set")

    def test_kaks_on_the_example(self):
        r = dnds.run_kaks(EX_CODON, self.tmp(), ["YN", "NG"], 1, "consecutive")
        self.assertEqual((r.exit_code, r.counts["ok"], r.counts["total"]), (0, 24, 24))
        first = r.rows[0]  # docs/python-api.md shows this row
        self.assertEqual((first["seq1"], first["seq2"], first["method"], first["dN"], first["dS"]),
                         ("spA_g01", "spB_g01", "NG", "0.0096066871", "0.46871109"))
        r = tools.dnds_pairs(EX_CODON)  # every pair: the two species' genes differ in length
        self.assertTrue(r["completed"])
        self.assertEqual(r["counts"]["ok"], 12)
        self.assertGreater(r["counts"]["error"], 0)
        self.assertTrue(any("differ in length" in m for m in r["engine_messages"]), r["engine_messages"])

    @unittest.skipUnless(dnds.diamond_path() or _REQUIRE_DIAMOND, "DIAMOND not found")
    def test_pipeline_on_the_example(self):
        if not dnds.diamond_path():
            self.fail("DIAMOND not found and COAL_REQUIRE_DIAMOND is set")
        r = tools.pipeline_pairs(EX_A, EX_B, methods=["YN"], threads=1)
        self.assertEqual((r["status"], r["exit_code"], r["pair_counts"]["ok"]), ("completed", 0, 12))
        self.assertEqual(set(r["counts"]["pairs"]) - {"total"}, set(dnds.STATUS_ORDER))
        self.assertFalse(os.path.exists(os.path.join(r["out_dir"], "work")))

    def test_pipeline_with_given_pairs(self):
        r = tools.pipeline_pairs(EX_A, EX_B, pairs_path=os.path.join(EX, "orthologs.tsv"), methods=["NG"])
        self.assertEqual((r["status"], r["pair_counts"]["ok"]), ("completed", 12))
        self.assertEqual(r["rows"][0]["dN"], "0.0096066871")


if __name__ == "__main__":
    unittest.main()
