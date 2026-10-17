"""APP-L: local users give file PATHS so large inputs are read in place, never copied (~50 MB input)."""
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

from coal_app import core, tools, viewer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class TestLargeInputByPath(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = tempfile.mkdtemp(prefix="coal_test_large_")
        r = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks", "gen_data.py"), "22000", "600",
                            cls.data], capture_output=True, text=True, encoding="utf-8", timeout=300)
        assert r.returncode == 0, r.stderr
        cls.pep = os.path.join(cls.data, "pep.fasta")
        cls.nuc = os.path.join(cls.data, "nuc.fasta")
        cls.sizes = {os.path.getsize(cls.pep), os.path.getsize(cls.nuc)}
        cls.shas = {_sha(cls.pep), _sha(cls.nuc)}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.data, True)

    def setUp(self):
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)

    def _copies(self, since, allowed):
        """Files created since `since` in RUNS_DIR (recursive) or the temp folder's top level that look like a
        copy of an input (same size or same content)."""
        found = []
        places = [os.path.join(dp, f) for dp, _, fs in os.walk(core.RUNS_DIR) for f in fs]
        tmp = tempfile.gettempdir()
        places += [os.path.join(tmp, f) for f in os.listdir(tmp)]
        for p in places:
            try:
                if not os.path.isfile(p) or p in allowed or os.path.getmtime(p) < since:
                    continue
                if os.path.getsize(p) in self.sizes or (os.path.getsize(p) > 1_000_000 and _sha(p) in self.shas):
                    found.append(p)
            except OSError:
                continue
        return found

    def test_input_is_about_50_mb(self):
        total = sum(os.path.getsize(p) for p in (self.pep, self.nuc))
        self.assertGreater(total, 45_000_000)

    def test_codon_align_by_path_makes_no_copy_and_output_is_right(self):
        t0 = time.time() - 1
        r = tools.codon_align(self.pep, self.nuc)
        self.assertTrue(r["ok"], r.get("error"))
        self.assertEqual(r["summary"]["ok"], 22000)
        self.assertEqual(self._copies(t0, {r["output_path"], r["report_path"]}), [])
        # correct by construction: each residue's codon in order, '---' for each gap
        prot = viewer.read_alignment(self.pep)
        dna = dict(viewer.read_alignment(self.nuc))
        out = dict(viewer.read_alignment(r["output_path"]))
        for sid, row in prot[:300] + prot[-300:]:
            d, k, want = dna[sid], 0, []
            for a in row:
                if a == "-":
                    want.append("---")
                else:
                    want.append(d[3 * k:3 * k + 3])
                    k += 1
            self.assertEqual(out[sid], "".join(want), sid)
        self.assertEqual(len(out), 22000)

    def test_gui_rendering_is_fast_and_small(self):  # R2-A3
        import tracemalloc
        out_dir = tempfile.mkdtemp(prefix="coal_test_render_")
        self.addCleanup(shutil.rmtree, out_dir, True)
        r = core.run_coal(self.pep, self.nuc, out_dir=out_dir, probe_codes=False)
        tracemalloc.start()
        t = time.perf_counter()
        try:
            frag, page_path = viewer.render_run(self.pep, r.output_path, r.report, self.nuc, 1, out_dir)
            table = core.health_table_html(r.report, links=True)
            elapsed = time.perf_counter() - t
            peak = tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()
        self.assertLess(elapsed, 20.0)  # ~5 s here under tracemalloc; was ~46 s without it
        self.assertLess(peak, 500_000_000)
        self.assertTrue(os.path.isfile(page_path))
        # CHANGED(GUI-10): was "showing the first 200 of 22000 sequences"; now an explicit preview (default size
        # 200 here: render_run without limits keeps the old caps)
        self.assertIn("Preview: 200 of 22,000 sequences", frag)
        self.assertLess(len(table), 300_000)
        # CHANGED(GUI-4): was "21,800 more sequences" (200 OK rows + a summary row); OK genes are now one collapsed
        # row naming the first core.HEALTH_OK_NAMES of them
        self.assertIn("All 22,000 genes match codon for codon", table)
        self.assertIn(f"and {22000 - core.HEALTH_OK_NAMES:,} more", table)

    def test_gui_path_resolution_reads_in_place(self):
        p, info = core.resolve_input(None, f'  "{self.nuc}"  ')  # pasted with quotes
        self.assertEqual(p, self.nuc)
        self.assertIn("MB", info)


class TestResolveInput(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.mkdtemp(prefix="coal_test_")
        self.addCleanup(shutil.rmtree, self.t, True)
        self.f = os.path.join(self.t, "a.fasta")
        with open(self.f, "w", encoding="utf-8") as fh:
            fh.write(">a\nMK\n")
        self._hosted = os.environ.pop("COAL_HOSTED", None)
        self.addCleanup(lambda: self._hosted is not None and os.environ.__setitem__("COAL_HOSTED", self._hosted))

    def test_path_wins_over_upload(self):
        self.assertEqual(core.resolve_input("/uploaded/x.fasta", self.f)[0], self.f)
        self.assertEqual(core.resolve_input("/uploaded/x.fasta", "  ")[0], "/uploaded/x.fasta")
        self.assertEqual(core.resolve_input(None, "")[0], None)

    def test_plain_errors(self):
        for bad, words in ((os.path.join(self.t, "nope.fasta"), "does not exist"), (self.t, "is a folder"),):
            with self.assertRaises(ValueError) as cm:
                core.resolve_input(None, bad)
            self.assertIn(words, str(cm.exception))

    @unittest.skipIf(os.name == "nt" or (hasattr(os, "geteuid") and os.geteuid() == 0),
                     "permission bits are not enforced for root / on Windows")
    def test_unreadable(self):
        os.chmod(self.f, 0)
        self.addCleanup(os.chmod, self.f, 0o644)
        with self.assertRaises(ValueError) as cm:
            core.resolve_input(None, self.f)
        self.assertIn("cannot be read", str(cm.exception))

    def test_hosted_refuses_paths(self):
        os.environ["COAL_HOSTED"] = "1"
        self.addCleanup(os.environ.pop, "COAL_HOSTED", None)
        with self.assertRaises(ValueError) as cm:
            core.resolve_input(None, self.f)
        self.assertIn("upload", str(cm.exception))
        self.assertFalse(core.path_inputs_enabled())

    def test_query_parameters_fill_the_boxes(self):
        self.assertTrue(core.path_inputs_enabled())
        self.assertEqual(core.paths_from_query({"protein": self.f, "dna": self.f, "x": "1"}), (self.f, self.f))
        self.assertEqual(core.paths_from_query({}), ("", ""))
        os.environ["COAL_HOSTED"] = "1"
        self.addCleanup(os.environ.pop, "COAL_HOSTED", None)
        self.assertEqual(core.paths_from_query({"protein": self.f}), ("", ""))


class TestMcpDefaults(unittest.TestCase):
    """APP-M: coal-gui serves MCP unless --no-mcp; the hosted Space does not."""

    def _launch_kwargs(self, argv):
        import types
        from coal_app import __main__ as entry
        seen = {}
        fake = types.ModuleType("coal_app.ui")
        fake.launch = lambda **kw: seen.update(kw)
        old = sys.modules.get("coal_app.ui")
        sys.modules["coal_app.ui"] = fake
        import coal_app
        old_attr = getattr(coal_app, "ui", None)
        coal_app.ui = fake
        try:
            self.assertEqual(entry.main(argv), 0)
        finally:
            if old is not None:
                sys.modules["coal_app.ui"] = old
            else:
                sys.modules.pop("coal_app.ui", None)
            if old_attr is not None:
                coal_app.ui = old_attr
            else:
                delattr(coal_app, "ui")
        return seen

    def test_local_gui_serves_mcp_by_default(self):
        self.assertIs(self._launch_kwargs([])["mcp"], True)
        self.assertIs(self._launch_kwargs(["--no-mcp"])["mcp"], False)
        self.assertEqual(self._launch_kwargs([])["host"], "127.0.0.1")

    def test_remote_host_needs_allow_remote(self):
        import contextlib
        import io
        from coal_app import __main__ as entry
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            self.assertEqual(entry.main(["--host", "0.0.0.0"]), 2)   # refused before Gradio is imported
        self.assertIn("--allow-remote", err.getvalue())
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(self._launch_kwargs(["--host", "0.0.0.0", "--allow-remote"])["host"], "0.0.0.0")
        self.assertEqual(self._launch_kwargs(["--host", "localhost"])["host"], "localhost")

    def test_space_does_not(self):
        import types
        from coal_app import space_app
        seen = {}
        fake = types.ModuleType("coal_app.ui")
        fake.launch = lambda **kw: seen.update(kw)
        old = sys.modules.get("coal_app.ui")
        sys.modules["coal_app.ui"] = fake
        import coal_app
        old_attr = getattr(coal_app, "ui", None)
        coal_app.ui = fake
        try:
            space_app.main()
        finally:
            if old is not None:
                sys.modules["coal_app.ui"] = old
            else:
                sys.modules.pop("coal_app.ui", None)
            if old_attr is not None:
                coal_app.ui = old_attr
            else:
                delattr(coal_app, "ui")
            os.environ.pop("COAL_HOSTED", None)
        self.assertIs(seen["mcp"], False)


if __name__ == "__main__":
    unittest.main()
