"""Engine runs started by the app: timeouts and signals leave no processes or temp files, run folders are pruned safely, genetic-code probes are small and cleaned up."""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from coal_app import core  # noqa: E402
from helpers import ROOT, RunsDirCase, _procs_mentioning, _tmp, _wait_gone, _write  # noqa: E402


def _load_demo(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(f"r3_{name}", os.path.join(ROOT, "demo", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@unittest.skipUnless(os.path.isdir("/proc/self") and os.name != "nt", "R2-A4: needs Linux /proc and POSIX signals")
class TestNoOrphanProcesses(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import subprocess, sys
        cls.data = tempfile.mkdtemp(prefix="coal_test_orphans_")
        r = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks", "gen_data.py"), "40000", "600",
                            cls.data], capture_output=True, text=True, timeout=300)
        assert r.returncode == 0, r.stderr
        cls.pep, cls.nuc = os.path.join(cls.data, "pep.fasta"), os.path.join(cls.data, "nuc.fasta")

    @classmethod
    def tearDownClass(cls):
        for pid in _procs_mentioning(cls.data):
            try:
                os.kill(pid, 9)
            except OSError:
                pass
        shutil.rmtree(cls.data, True)

    def test_app_timeout_leaves_no_processes_or_temp_files(self):
        out_dir = tempfile.mkdtemp(prefix="coal_test_orphan_out_")
        self.addCleanup(shutil.rmtree, out_dir, True)
        code, err, _ = core._run_engine(self.pep, self.nuc, os.path.join(out_dir, "o.fasta"), "fasta", 1,
                                        False, True, False, None, timeout=0.8)
        self.assertNotEqual(code, 0)
        self.assertIn("did not finish within", err)
        self.assertEqual(_wait_gone(self.data), [])
        self.assertEqual([n for n in os.listdir(out_dir) if n.startswith(".coal_")], [])

    def _compat(self, sig):
        import signal, subprocess, sys, time
        p = subprocess.Popen([sys.executable, "-m", "coal_compat", self.pep, self.nuc, "-output", "fasta",
                              "-nomismatch"], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        try:
            time.sleep(0.8)  # let the engine start its work (a full run takes ~2.5 s)
            p.send_signal(sig)
            err = p.communicate(timeout=60)[1].decode()
        finally:
            if p.poll() is None:
                p.kill()
                p.communicate(timeout=60)
        return p.returncode, err

    def test_compat_ctrl_c(self):
        import signal
        rc, err = self._compat(signal.SIGINT)
        self.assertEqual(rc, 130, err)
        self.assertNotIn("Traceback", err)
        self.assertIn("coal: interrupted", err)
        self.assertEqual(_wait_gone(self.data), [])

    def test_compat_sigterm_stops_the_engine(self):
        import signal
        rc, err = self._compat(signal.SIGTERM)
        self.assertEqual(rc, 128 + signal.SIGTERM, err)
        self.assertNotIn("Traceback", err)
        self.assertEqual(_wait_gone(self.data), [])


@unittest.skipUnless(os.path.isdir("/proc/self") and os.name != "nt", "R3-A4: needs Linux /proc and POSIX signals")
class TestDemoTimeoutsStopWorkers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = tempfile.mkdtemp(prefix="coal_test_r3a4_")
        r = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks", "gen_data.py"), "40000", "600",
                            cls.data], capture_output=True, text=True, timeout=300)
        assert r.returncode == 0, r.stderr
        cls.pep, cls.nuc = os.path.join(cls.data, "pep.fasta"), os.path.join(cls.data, "nuc.fasta")

    @classmethod
    def tearDownClass(cls):
        for pid in _procs_mentioning(cls.data):
            try:
                os.kill(pid, 9)
            except OSError:
                pass
        shutil.rmtree(cls.data, True)

    def out_dir(self):
        d = tempfile.mkdtemp(prefix="coal_test_r3a4_out_", dir=self.data)
        return d

    def test_diamond_run_batch_timeout(self):
        import argparse
        dp = _load_demo("diamond_pipeline")
        recs = dp.read_fasta(self.pep)
        dna = {rid: seq for rid, _, seq in dp.read_fasta(self.nuc)}
        items, cds = [], {"a": {}, "b": {}}
        for (x, _, rx), (y, _, ry) in zip(recs[0::2], recs[1::2]):
            items.append((x, y, rx, ry))
            cds["a"][x], cds["b"][y] = dna[x], dna[y]
        out = self.out_dir()
        a = argparse.Namespace(a_name="a", b_name="b", out=out, codontable=1, timeout=0.8)
        res = dp.run_batch(items, a, "", "batch", [], cds)
        self.assertIsInstance(res, tuple)
        self.assertIn("did not finish", res[1])
        self.assertEqual(_wait_gone(self.data), [])
        self.assertEqual([n for n in os.listdir(out) if n.startswith(".coal_")], [])


# A stand-in engine that behaves like a stuck one: it ignores SIGTERM and has a worker (a forked child, holding
# the same stdout/stderr pipes, as the engine's process-pool workers do) that outlives it. argv[1]: file for the
# worker's pid.
_STUCK_ENGINE = """
import os, signal, sys, time
signal.signal(signal.SIGTERM, signal.SIG_IGN)
pid = os.fork()
if pid == 0:
    time.sleep(600)
    os._exit(0)
with open(sys.argv[1], "w") as fh:
    fh.write(str(pid))
time.sleep(600)
"""


def _alive(pid: int) -> bool:
    try:
        with open(f"/proc/{pid}/stat", "rb") as fh:
            return fh.read().decode().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return False


@unittest.skipUnless(os.path.isdir("/proc/self") and os.name != "nt", "needs Linux /proc and POSIX signals")
class TestStopEngineNeverHangs(unittest.TestCase):
    """stop_engine after a timeout: SIGTERM, grace, then SIGKILL of the engine's whole process group. Before the
    group kill, an engine that outlived the grace (a slow runner: SIGTERM waits for the in-flight batches) was
    SIGKILLed alone, its orphaned workers kept its pipes open and the final communicate() waited for them forever
    (the CI app suite hung until the job timeout)."""

    def test_timed_out_engine_and_its_worker_are_gone_quickly(self):
        import subprocess
        import time
        d = tempfile.mkdtemp(prefix="coal_test_stuck_")
        self.addCleanup(shutil.rmtree, d, True)
        pidfile = os.path.join(d, "worker.pid")
        t0 = time.monotonic()
        with self.assertRaises(subprocess.TimeoutExpired):
            core.run_engine_command([sys.executable, "-c", _STUCK_ENGINE, pidfile], timeout=1.5, temp_dir=d)
        elapsed = time.monotonic() - t0
        self.assertLess(elapsed, 30, "stop_engine took too long (hung on the orphaned worker?)")
        with open(pidfile, encoding="utf-8") as fh:
            worker = int(fh.read())
        self.addCleanup(lambda: os.kill(worker, 9) if _alive(worker) else None)
        self.assertFalse(_alive(worker), f"the engine's worker {worker} survived stop_engine")


class TestPruneRunFolders(RunsDirCase):
    def old(self, path):
        t = __import__("time").time() - 48 * 3600
        os.utime(path, (t, t))

    def test_only_coal_run_folders_are_pruned(self):
        mine = core.new_run_dir()
        user = os.path.join(core.RUNS_DIR, "thesis_data")
        os.makedirs(user)
        with open(os.path.join(user, "chapter1.fa"), "w", encoding="utf-8") as fh:
            fh.write(">x\nA\n")
        fake = os.path.join(core.RUNS_DIR, "run_not_ours")
        os.makedirs(fake)
        for p in (mine, user, fake):
            self.old(p)
        core.prune_runs()
        self.assertFalse(os.path.exists(mine))
        self.assertTrue(os.path.isfile(os.path.join(user, "chapter1.fa")))
        self.assertTrue(os.path.isdir(fake))

    def test_fresh_run_folder_kept(self):
        mine = core.new_run_dir()
        core.prune_runs()
        self.assertTrue(os.path.isdir(mine))


class TestTempDirs(RunsDirCase):
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
            _write(os.path.join(x, core.RUN_MARKER), "")
        past = time.time() - 25 * 3600
        os.utime(old, (past, past))
        os.utime(os.path.join(old, "f.txt"), (past, past))
        removed = core.prune_runs(max_age_hours=24)
        self.assertEqual(removed, 1)
        self.assertFalse(os.path.exists(old))
        self.assertTrue(os.path.exists(new))


class TestProbeCost(RunsDirCase):
    def test_probes_are_small_and_time_limited(self):
        calls = []
        real = core._run_engine

        def spy(protein, dna, out, fmt, table, nogap, nomismatch, keep_going, report_path, timeout=None, **kw):
            calls.append((protein, len(core.read_ids(protein)), timeout))  # **kw: strict, blockonly
            return real(protein, dna, out, fmt, table, nogap, nomismatch, keep_going, report_path, timeout=timeout,
                        **kw)

        core._run_engine = spy
        self.addCleanup(setattr, core, "_run_engine", real)
        # the app probes only in its fallback path (an engine without --report); the engine
        # writes the report itself (and has --keep-going, whose "skipped" lines the fallback does not parse), so this
        # test of the fallback's probe cost pins that path
        import dataclasses
        feats = dataclasses.replace(core.detect_features(), report=False, keep_going=False)
        real_df = core.detect_features
        core.detect_features = lambda: feats
        self.addCleanup(setattr, core, "detect_features", real_df)
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
        # only TGA = Trp is needed, so code 4 (which differs from code 1 only there) is the closest
        self.assertEqual(hint[0]["suggested_codontable"], 4)
        self.assertIn("COAL tested 200 of the 250 sequences", hint[0]["text"])  # what was tested
        self.assertIn("200 of them match their protein", hint[0]["text"])
        self.assertEqual([f for f in os.listdir(t) if f.startswith("probe")], [])


if __name__ == "__main__":
    unittest.main()


class TestUnusableRunsDir(unittest.TestCase):
    """COAL_RUNS_DIR pointing somewhere no folder can be made: one plain sentence, never a traceback."""

    def _with_runs_dir(self, path):
        old = core.RUNS_DIR
        core.RUNS_DIR = path
        self.addCleanup(setattr, core, "RUNS_DIR", old)

    def test_new_run_dir_under_an_impossible_path(self):
        self._with_runs_dir("/proc/nonexistent/runs" if os.name != "nt" else "?:\\nonexistent\\runs")
        with self.assertRaisesRegex(ValueError, r"Could not create the results folder .*nonexistent.*runs: .+") as cm:
            core.new_run_dir()
        self.assertNotIn("Traceback", str(cm.exception))
        self.assertIsInstance(cm.exception.__cause__, OSError)

    def test_new_run_dir_when_a_file_is_in_the_way(self):
        d = _tmp(self)
        f = os.path.join(d, "runs")
        _write(f, "not a folder\n")
        self._with_runs_dir(f)
        with self.assertRaisesRegex(ValueError, "Could not create the results folder"):
            core.new_run_dir()
        with self.assertRaisesRegex(ValueError, "Could not create the results folder"):
            core.ensure_runs_dir()

    def test_mcp_tool_reports_a_tool_error(self):
        from coal_app import tools
        self._with_runs_dir("/proc/nonexistent/runs" if os.name != "nt" else "?:\\nonexistent\\runs")
        p, d = core.example_paths("clean")
        with self.assertRaisesRegex(tools.CoalToolError, "Could not create the results folder"):
            tools.codon_align(p, d)

