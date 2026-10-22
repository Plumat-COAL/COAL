"""Regression tests for the robustness campaign (devices as inputs, the report path, dying workers, checkpoint
corruption, closed pipes, unusable TMPDIR, disk-full messages ...).

Every case runs through the Python engine and, when it is built (build/coal-cpp or COAL_CPP_BIN), through coal-cpp,
and the two must agree on stdout, stderr and the exit code. Every input is built inside a temporary directory
created by the test; -o is only ever pointed at paths inside that directory (or at device nodes that discard data).
"""
from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
TIMEOUT = 60
WINDOWS = os.name == "nt"
POSIX = not WINDOWS

PEP = ">s1\nMKV-AL\n>s2\nMKVSAL\n"
NUC = ">s1\nATGAAAGTTGCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n"
EXPECTED = ">s1\nATGAAAGTT---GCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n"


def cpp_binary() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


CPP = cpp_binary()
PY_CMD = [sys.executable, "-m", "coal"]


def engines(cpp: bool = True):
    """(name, command prefix) of the engines to run: Python always, coal-cpp when built."""
    out = [("python", PY_CMD)]
    if cpp and CPP is not None:
        out.append(("cpp", [str(CPP)]))
    return out


class EngineCase(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory(prefix="coal_robust_")
        self.d = Path(self._td.name)
        (self.d / "p.fa").write_text(PEP)
        (self.d / "n.fa").write_text(NUC)

    def tearDown(self):
        self._td.cleanup()

    def env(self, **extra):
        env = dict(os.environ, PYTHONPATH=str(SRC), COAL_ENGINE="python")
        env.update(extra)
        return env

    def run_engine(self, cmd, args, stdin=subprocess.DEVNULL, timeout=TIMEOUT, **kw):
        if isinstance(stdin, (bytes, str)):
            kw["input"] = stdin
            stdin = None
        return subprocess.run([*cmd, *args], cwd=self.d, env=kw.pop("env", None) or self.env(), stdin=stdin,
                              capture_output=True, timeout=timeout, **kw)

    def both(self, args, stdin=subprocess.DEVNULL, cpp=True, **kw):
        """Run both engines; assert they agree on (rc, stdout, stderr); return that triple (bytes)."""
        results = {}
        for name, cmd in engines(cpp):
            r = self.run_engine(cmd, args, stdin=stdin, **kw)
            results[name] = (r.returncode, r.stdout, r.stderr)
        if len(results) == 2:
            self.assertEqual(results["python"], results["cpp"], f"engines differ on {args}")
        return results["python"]

    def leftovers(self, where=None):
        where = where or self.d
        return sorted(p.name for p in where.iterdir() if p.name.startswith(".coal_") or ".tmp" in p.name)


# ---------------------------------------------------------------- 1. devices as inputs
@unittest.skipUnless(POSIX and os.path.exists("/dev/zero"), "needs POSIX device nodes")
class TestDeviceInputs(EngineCase):
    def test_device_path_is_refused(self):
        for args in (["/dev/zero", "n.fa"], ["p.fa", "/dev/zero"], ["p.fa", "n.fa", "/dev/zero"]):
            with self.subTest(args=args):
                rc, out, err = self.both(args, timeout=10)
                self.assertEqual((rc, out), (1, b""))
                self.assertEqual(err, b"coal: error: cannot read /dev/zero: it is a device, not a file\n")
        self.assertEqual(self.leftovers(), [])

    def test_stdin_device_is_refused(self):
        with open("/dev/zero", "rb") as zero:
            rc, out, err = self.both(["-", "n.fa"], stdin=zero, timeout=10)
        self.assertEqual((rc, out), (1, b""))
        self.assertEqual(err, b"coal: error: cannot read standard input: it is a device, not a file\n")

    def test_null_device_is_still_an_empty_input(self):
        rc, out, err = self.both(["-", "n.fa"], stdin=subprocess.DEVNULL)
        self.assertEqual(rc, 1)
        self.assertIn(b"standard input: no sequences found", err)
        rc, out, err = self.both(["p.fa", "/dev/null"])
        self.assertEqual(rc, 1)
        self.assertIn(b"/dev/null: no sequences found", err)

    def test_pipe_input_still_works(self):
        rc, out, err = self.both(["-", "n.fa"], stdin=PEP.encode())
        self.assertEqual((rc, out, err), (0, EXPECTED.encode(), b""))


# ---------------------------------------------------------------- 2. the report path
class TestReportPath(EngineCase):
    def test_report_naming_an_input_is_a_usage_error(self):
        for args, what in ((["p.fa", "n.fa", "--report", "n.fa"], "n.fa is an input file"),
                           (["p.fa", "n.fa", "--report", "./p.fa"], "./p.fa is an input file"),
                           (["p.fa", "n.fa", "-o", "o.fa", "--report", "o.fa"], "o.fa is the output file"),
                           (["p.fa", "n.fa", "-o", "o.fa", "--report", "./o.fa"], "./o.fa is the output file")):
            with self.subTest(args=args):
                rc, out, err = self.both(args)
                self.assertEqual((rc, out), (2, b""))
                self.assertTrue(err.startswith(b"usage: coal "), err)
                self.assertEqual(err.splitlines()[-1], f"coal: error: the report path {what}".encode())
        self.assertEqual((self.d / "n.fa").read_text(), NUC)
        self.assertEqual((self.d / "p.fa").read_text(), PEP)
        self.assertFalse((self.d / "o.fa").exists())

    def test_batch_job_report_naming_its_output_is_refused(self):
        (self.d / "jobs.txt").write_text("p.fa n.fa -o o.fa\np.fa n.fa -o o2.fa --report o2.fa\n")
        r = self.run_engine(PY_CMD, ["--batch", "jobs.txt"])
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertEqual(r.stderr.splitlines()[-1],
                         b"coal: error: jobs.txt line 2: the report path o2.fa is the output file")
        self.assertFalse((self.d / "o.fa").exists())
        (self.d / "jobs.txt").write_text("p.fa n.fa -o o.fa\n")
        r = self.run_engine(PY_CMD, ["--batch", "jobs.txt", "--report", "n.fa"])
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertEqual(r.stderr.splitlines()[-1], b"coal: error: the report path n.fa is an input file")
        self.assertEqual((self.d / "n.fa").read_text(), NUC)

    def test_report_is_written_through_a_symlink(self):
        target = self.d / "target.json"
        target.write_text("old\n")
        link = self.d / "link.json"
        try:
            link.symlink_to(target)
        except (OSError, NotImplementedError):
            self.skipTest("cannot create symlinks here")
        for name, cmd in engines():
            with self.subTest(engine=name):
                target.write_text("old\n")
                r = self.run_engine(cmd, ["p.fa", "n.fa", "-o", "o.fa", "--report", "link.json"])
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertTrue(link.is_symlink())
                self.assertEqual(json.loads(target.read_text())["schema"], "coal.report/1")
                self.assertEqual(self.leftovers(), [])

    @unittest.skipUnless(POSIX and os.path.exists("/dev/null"), "needs /dev/null")
    def test_report_to_a_device_does_not_replace_it(self):
        for name, cmd in engines():
            with self.subTest(engine=name):
                r = self.run_engine(cmd, ["p.fa", "n.fa", "-o", "o.fa", "--report", "/dev/null"])
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertTrue(stat.S_ISCHR(os.stat("/dev/null").st_mode), "/dev/null was replaced by a file!")

    @unittest.skipUnless(hasattr(os, "mkfifo"), "needs FIFOs")
    def test_report_to_a_fifo_is_written_not_replaced(self):
        for name, cmd in engines():
            with self.subTest(engine=name):
                fifo = self.d / f"{name}.fifo"
                os.mkfifo(fifo)
                reader = subprocess.Popen(["cat", str(fifo)], stdout=subprocess.PIPE)
                try:
                    r = self.run_engine(cmd, ["p.fa", "n.fa", "-o", "o.fa", "--report", str(fifo)])
                    data, _ = reader.communicate(timeout=TIMEOUT)
                finally:
                    reader.kill()
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertEqual(json.loads(data)["schema"], "coal.report/1")
                self.assertTrue(stat.S_ISFIFO(fifo.stat().st_mode))


# ---------------------------------------------------------------- 3 + 4. worker processes
try:
    import psutil
except ImportError:  # pragma: no cover
    psutil = None


def _big_input(where: Path, n_seqs: int = 100_000):
    """A protein alignment + DNA large enough for a -j 2 run to use two worker processes for a few seconds."""
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from gen_data import generate
    generate(n_seqs, 150, where)
    return where / "pep.fasta", where / "nuc.fasta"


def _children(pid: int):
    try:
        return psutil.Process(pid).children()
    except psutil.Error:
        return []


def _running(proc) -> bool:
    try:
        return proc.is_running() and proc.status() != psutil.STATUS_ZOMBIE
    except psutil.Error:
        return False


@unittest.skipUnless(POSIX and psutil is not None, "POSIX signals and psutil")
class TestWorkerProcesses(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._td = tempfile.TemporaryDirectory(prefix="coal_workers_")
        cls.d = Path(cls._td.name)
        cls.pep, cls.nuc = _big_input(cls.d)

    @classmethod
    def tearDownClass(cls):
        cls._td.cleanup()

    def start(self, args):
        env = dict(os.environ, PYTHONPATH=str(SRC), COAL_ENGINE="python")
        p = subprocess.Popen([*PY_CMD, str(self.pep), str(self.nuc), "-j", "2", *args], cwd=self.d, env=env,
                             stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        t0 = time.monotonic()
        kids = []
        while time.monotonic() - t0 < 30 and len(kids) < 2:
            kids = _children(p.pid)
            time.sleep(0.05)
        self.assertEqual(len(kids), 2, "the run did not start two workers")
        time.sleep(0.3)
        return p, kids

    def test_a_killed_worker_ends_the_run_with_an_error(self):
        out = self.d / "out.fa"
        p, kids = self.start(["-o", str(out)])
        kids[0].kill()
        try:
            _, err = p.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            p.kill()
            for k in kids:
                try:
                    k.kill()
                except psutil.Error:
                    pass
            self.fail("coal hung after a worker was killed")
        self.assertEqual(p.returncode, 1, err)
        self.assertEqual(err.decode().splitlines()[0], "coal: error: a worker process died (out of memory?)")
        self.assertFalse(out.exists())
        self.assertEqual([x.name for x in self.d.iterdir() if x.name.startswith(("out.fa", ".coal_"))], [])
        time.sleep(0.5)
        self.assertEqual([k.pid for k in kids if _running(k)], [])

    def test_workers_die_with_the_main_process(self):
        p, kids = self.start([])
        os.kill(p.pid, 9)
        p.wait()
        t0 = time.monotonic()
        p.stdout.read()  # EOF only once every process holding the pipe is gone
        self.assertLess(time.monotonic() - t0, 3.0, "a worker kept coal's stdout open")
        deadline = time.monotonic() + 3.0
        while time.monotonic() < deadline and any(_running(k) for k in kids):
            time.sleep(0.1)
        self.assertEqual([k.pid for k in kids if _running(k)], [], "workers survived the main process")
        p.stderr.close()


# ---------------------------------------------------------------- 5 + 6 + 7. checkpoints
def interrupt(cmd, cwd, after=2, extra_env=None):
    """Run `cmd` with frequent checkpoints and stop it (SIGTERM) after the `after`-th one, leaving the checkpoint."""
    import threading
    env = {**os.environ, "PYTHONPATH": str(SRC), "COAL_ENGINE": "python", "COAL_BATCH_NT": "20000",
           "COAL_CHECKPOINT_EVERY": "0", "COAL_TEST_HOOKS": "1", "COAL_TEST_PAUSE_AFTER_CHECKPOINT": str(after),
           **(extra_env or {})}
    p = subprocess.Popen(cmd, cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                         stderr=subprocess.PIPE, text=True)
    dog = threading.Timer(60, p.kill)
    dog.start()
    try:
        for line in p.stderr:
            if line.startswith("coal: test: paused"):
                break
        p.terminate()
        p.wait(timeout=30)
    finally:
        dog.cancel()
        p.stderr.close()


UNREADABLE = b"coal: notice: found an unreadable checkpoint, so nothing can be resumed for out.fa; starting a normal run"
INCOMPLETE = (b"coal: notice: found a checkpoint whose partial output is missing or incomplete, so nothing can be "
              b"resumed for out.fa; starting a normal run")

# edits of a valid checkpoint (a dict) that must make it "unreadable"
CORRUPTIONS = {
    "completed_huge": lambda st: st.update(completed=10 ** 9),
    "completed_99999": lambda st: st.update(completed=99999),
    "completed_negative": lambda st: st.update(completed=-5),
    "completed_bool": lambda st: st.update(completed=True),
    "completed_str": lambda st: st.update(completed="12"),
    "top_level_list": lambda st: [st],
    "sink_bytes_str": lambda st: st.update(sink_bytes="abc"),
    "sink_bytes_float": lambda st: st.update(sink_bytes=1.5e300),
    "sink_bytes_negative": lambda st: st.update(sink_bytes=-1),
    "no_skip_idx": lambda st: st.pop("skip_idx"),
    "skipped_wrong": lambda st: st.update(skipped=7),
    "skip_idx_out_of_range": lambda st: st.update(skipped=1, skip_idx=[10 ** 6]),
    "rec_int": lambda st: st.update(rec=5),
    "rec_kinds_short": lambda st: st["rec"].update(kinds="o"),
    "nseq_wrong": lambda st: st.update(nseq=4),
    "drop_str": lambda st: st.update(drop=["x"]),
}


@unittest.skipIf(WINDOWS, "POSIX signals")
class TestCheckpointCorruption(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._td = tempfile.TemporaryDirectory(prefix="coal_ckpt_")
        cls.d = Path(cls._td.name)
        sys.path.insert(0, str(ROOT / "benchmarks"))
        from gen_data import generate
        generate(3000, 15, cls.d, seed=5)
        env = dict(os.environ, PYTHONPATH=str(SRC), COAL_ENGINE="python")
        cls.ref = subprocess.run([*PY_CMD, "pep.fasta", "nuc.fasta"], cwd=cls.d, env=env, capture_output=True,
                                 timeout=TIMEOUT).stdout
        cls.snapshots = {}
        for name, cmd in engines():
            cls.clear()
            interrupt([*cmd, "pep.fasta", "nuc.fasta", "-o", "out.fa"], cls.d)
            cls.snapshots[name] = {f.name: f.read_bytes() for f in cls.d.iterdir() if f.name.startswith("out.fa")}
            assert "out.fa.coal-checkpoint.json" in cls.snapshots[name], name

    @classmethod
    def tearDownClass(cls):
        cls._td.cleanup()

    @classmethod
    def clear(cls):
        import shutil
        for f in cls.d.iterdir():
            if f.name.startswith("out.fa"):
                shutil.rmtree(f) if f.is_dir() else f.unlink()

    def restore(self, name):
        self.clear()
        for k, v in self.snapshots[name].items():
            (self.d / k).write_bytes(v)
        return self.d / "out.fa.coal-checkpoint.json"

    def resume(self, cmd, env_extra=None):
        env = {**os.environ, "PYTHONPATH": str(SRC), "COAL_ENGINE": "python", "COAL_BATCH_NT": "20000",
               **(env_extra or {})}
        return subprocess.run([*cmd, "pep.fasta", "nuc.fasta", "-o", "out.fa", "--resume"], cwd=self.d, env=env,
                              stdin=subprocess.DEVNULL, capture_output=True, timeout=120)

    def check_normal_run(self, r, notice):
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stderr.splitlines()[0], notice)
        self.assertEqual((self.d / "out.fa").read_bytes(), self.ref)
        self.assertEqual([f.name for f in self.d.iterdir() if f.name.startswith("out.fa.")], [])

    def test_edited_checkpoints_start_a_normal_run(self):
        for name, cmd in engines():
            for case, edit in CORRUPTIONS.items():
                with self.subTest(engine=name, case=case):
                    ck = self.restore(name)
                    st = json.loads(ck.read_text())
                    new = edit(st)
                    ck.write_text(json.dumps(st if new is None else new))
                    self.check_normal_run(self.resume(cmd), UNREADABLE)

    def test_deeply_nested_checkpoint(self):
        for name, cmd in engines():
            with self.subTest(engine=name):
                ck = self.restore(name)
                ck.write_text("[" * 100000 + "]" * 100000)
                self.check_normal_run(self.resume(cmd), UNREADABLE)

    def test_truncated_and_garbage_checkpoint(self):
        for name, cmd in engines():
            for case, data in (("truncated", None), ("garbage", b"\x00\xff{"), ("empty", b"")):
                with self.subTest(engine=name, case=case):
                    ck = self.restore(name)
                    ck.write_bytes(ck.read_bytes()[:50] if data is None else data)
                    self.check_normal_run(self.resume(cmd), UNREADABLE)

    def test_missing_entries_file(self):
        for name, cmd in engines():
            with self.subTest(engine=name):
                self.restore(name)
                (self.d / "out.fa.coal-entries.jsonl").unlink()
                self.check_normal_run(self.resume(cmd), INCOMPLETE)

    def test_missing_partial_file(self):
        for name, cmd in engines():
            with self.subTest(engine=name):
                self.restore(name)
                (self.d / "out.fa.coal-partial").unlink()
                self.check_normal_run(self.resume(cmd), INCOMPLETE)

    def test_checkpoint_path_is_a_folder(self):
        for name, cmd in engines():
            with self.subTest(engine=name):
                ck = self.restore(name)
                ck.unlink()
                ck.mkdir()
                r = self.resume(cmd, {"COAL_CHECKPOINT_EVERY": "0"})
                self.assertEqual(r.returncode, 0, r.stderr)
                lines = r.stderr.splitlines()
                self.assertEqual(lines[0], UNREADABLE)
                self.assertEqual(lines[1], b"coal: notice: cannot write the checkpoint " + bytes(ck) +
                                 b": is a directory; the run cannot be resumed if it is interrupted")
                self.assertEqual(len(lines), 2, r.stderr)
                self.assertEqual((self.d / "out.fa").read_bytes(), self.ref)
                self.assertTrue(ck.is_dir())
                self.assertEqual([f.name for f in self.d.iterdir() if f.name.startswith("out.fa.")],
                                 ["out.fa.coal-checkpoint.json"])

    def test_valid_checkpoint_still_resumes(self):
        for name, cmd in engines():
            with self.subTest(engine=name):
                self.restore(name)
                r = self.resume(cmd)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertTrue(r.stderr.startswith(b"coal: notice: resuming out.fa: "), r.stderr)
                self.assertEqual((self.d / "out.fa").read_bytes(), self.ref)

    @unittest.skipIf(CPP is None, "coal-cpp not built (build/coal-cpp or COAL_CPP_BIN)")
    def test_other_engines_checkpoint(self):
        # the Python engine on a coal-cpp checkpoint, and coal-cpp on a Python one: a clear error, files kept
        self.restore("cpp")
        r = self.resume(PY_CMD)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(r.stderr.splitlines()[0],
                         b"coal: error: out.fa was interrupted by the C++ engine (coal-cpp), which must resume it: run "
                         b"again with COAL_ENGINE=cpp, or start a new run without --resume")
        self.assertTrue((self.d / "out.fa.coal-checkpoint.json").exists())
        self.restore("python")
        r = self.resume([str(CPP)])
        self.assertEqual(r.returncode, 1)
        self.assertEqual(r.stderr.splitlines()[0],
                         b"coal: error: out.fa was interrupted by the Python engine (coal-python), which must resume "
                         b"it: run again with COAL_ENGINE=python, or start a new run without --resume")
        # the command line picks the engine that wrote the checkpoint (COAL_ENGINE unset = auto)
        self.restore("cpp")
        r = self.resume(PY_CMD, {"COAL_ENGINE": "", "COAL_CPP": str(CPP)})
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stderr.splitlines()[0],
                         b"coal: notice: out.fa was interrupted by the C++ engine (coal-cpp); resuming with it")
        self.assertEqual((self.d / "out.fa").read_bytes(), self.ref)
        self.restore("python")
        r = self.resume(PY_CMD, {"COAL_ENGINE": "", "COAL_CPP": str(CPP)})
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(r.stderr.startswith(b"coal: notice: resuming out.fa: "), r.stderr)
        # a coal-cpp checkpoint but no coal-cpp any more: the Python engine says what to do
        self.restore("cpp")
        r = self.resume(PY_CMD, {"COAL_ENGINE": "", "COAL_CPP": str(self.d / "nonexistent")})
        self.assertEqual(r.returncode, 1)
        self.assertIn(b"interrupted by the C++ engine (coal-cpp), which must resume it", r.stderr)


# ---------------------------------------------------------------- 8. -o /dev/fd/1, 9. --json | head
@unittest.skipUnless(POSIX and os.path.exists("/dev/fd/1"), "needs /dev/fd")
class TestOutputToStdoutDevice(EngineCase):
    def test_dev_fd_1_into_a_pipe(self):
        for name, cmd in engines():
            for path in ("/dev/fd/1", "/dev/stdout"):
                with self.subTest(engine=name, path=path):
                    # stdout is a pipe here (capture_output): its real path is "pipe:[N]", which no folder holds
                    r = self.run_engine(cmd, ["p.fa", "n.fa", "-o", path])
                    self.assertEqual((r.returncode, r.stderr), (0, b""), r.stderr)
                    self.assertEqual(r.stdout, EXPECTED.encode())
                    self.assertEqual(self.leftovers(), [])

    def test_dev_stdout_resume_is_a_notice(self):
        rc, out, err = self.both(["p.fa", "n.fa", "-o", "/dev/stdout", "--resume"])
        self.assertEqual((rc, out), (0, EXPECTED.encode()))
        self.assertEqual(err, b"coal: notice: /dev/stdout cannot be resumed (the output is not a regular file); "
                              b"starting a normal run\n")


@unittest.skipIf(WINDOWS, "POSIX pipes")
class TestJsonIntoClosedPipe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._td = tempfile.TemporaryDirectory(prefix="coal_jsonpipe_")
        cls.d = Path(cls._td.name)
        sys.path.insert(0, str(ROOT / "benchmarks"))
        from gen_data import generate
        generate(20000, 15, cls.d, seed=5)  # a report larger than a pipe buffer (20,000 IDs)

    @classmethod
    def tearDownClass(cls):
        cls._td.cleanup()

    def closed_pipe(self, cmd):
        env = dict(os.environ, PYTHONPATH=str(SRC), COAL_ENGINE="python")
        p = subprocess.Popen(cmd, cwd=self.d, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        p.stdout.read(1)
        p.stdout.close()  # like `| head`
        _, err = p.communicate(timeout=TIMEOUT)
        return p.returncode, err

    def test_json_report_into_head(self):
        for name, cmd in engines():
            with self.subTest(engine=name):
                rc, err = self.closed_pipe([*cmd, "pep.fasta", "nuc.fasta", "-o", "o.fa", "--json"])
                self.assertEqual((rc, err), (141, b""))
                self.assertTrue((self.d / "o.fa").exists())

    def test_batch_json_report_into_head(self):
        (self.d / "jobs.txt").write_text("pep.fasta nuc.fasta -o o1.fa\npep.fasta nuc.fasta -o o2.fa\n")
        rc, err = self.closed_pipe([*PY_CMD, "--batch", "jobs.txt", "--json", "-q"])
        self.assertEqual((rc, err), (141, b""))


# ---------------------------------------------------------------- 10. TMPDIR unusable
class TestUnusableTmpdir(EngineCase):
    def test_missing_or_file_tmpdir_falls_back(self):
        for tmpdir in (str(self.d / "nonexistent"), str(self.d / "p.fa")):
            for args, stdin in ((["p.fa", "n.fa"], subprocess.DEVNULL),            # no temp file needed
                                (["-", "n.fa"], PEP.encode()),                      # stdin is spooled
                                (["p.fa", "n.fa", "-output", "clustal"], subprocess.DEVNULL)):  # staged rows
                with self.subTest(tmpdir=tmpdir, args=args):
                    rc, out, err = self.both(args, stdin=stdin, env=self.env(TMPDIR=tmpdir, TEMP=tmpdir, TMP=tmpdir))
                    self.assertEqual((rc, err), (0, b""), err)
                    self.assertTrue(out.startswith(b">s1\n") or out.startswith(b"CLUSTAL"), out)
        self.assertEqual(self.leftovers(), [])


# ---------------------------------------------------------------- 11. disk full
def _mount_tiny_tmpfs(where: Path) -> bool:
    """A 1 MB tmpfs at `where` (needs root on Linux); False when that is not possible here."""
    if not sys.platform.startswith("linux") or os.geteuid() != 0:
        return False
    try:
        r = subprocess.run(["mount", "-t", "tmpfs", "-o", "size=1m", "tmpfs", str(where)], capture_output=True,
                           timeout=30)
    except (OSError, subprocess.SubprocessError):
        return False
    return r.returncode == 0


@unittest.skipIf(WINDOWS, "POSIX device nodes")
class TestDiskFull(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._td = tempfile.TemporaryDirectory(prefix="coal_full_")
        cls.d = Path(cls._td.name)
        sys.path.insert(0, str(ROOT / "benchmarks"))
        from gen_data import generate
        generate(40000, 15, cls.d, seed=5)  # DNA and output of 2-6 MB, more than the tiny file system holds
        cls.full = cls.d / "full"
        cls.full.mkdir()
        cls.mounted = _mount_tiny_tmpfs(cls.full)

    @classmethod
    def tearDownClass(cls):
        if cls.mounted:
            subprocess.run(["umount", str(cls.full)], capture_output=True, timeout=30)
        cls._td.cleanup()

    def run_both(self, args, stdin=subprocess.DEVNULL, **env_extra):
        results = {}
        for name, cmd in engines():
            env = dict(os.environ, PYTHONPATH=str(SRC), COAL_ENGINE="python", **env_extra)
            kw = {"input": stdin} if isinstance(stdin, bytes) else {"stdin": stdin}
            r = subprocess.run([*cmd, *args], cwd=self.d, env=env, capture_output=True, timeout=120, **kw)
            results[name] = (r.returncode, r.stdout, r.stderr)
            for f in self.full.iterdir():  # no temp file left behind on the full file system
                self.assertFalse(f.name.startswith(".coal"), f"{name} left {f.name}")
        if len(results) == 2:
            self.assertEqual(results["python"], results["cpp"], f"engines differ on {args}")
        return results["python"]

    @unittest.skipUnless(os.path.exists("/dev/full"), "needs /dev/full")
    def test_output_to_dev_full(self):
        for fmt in ("fasta", "clustal"):
            with self.subTest(fmt=fmt):
                rc, out, err = self.run_both(["pep.fasta", "nuc.fasta", "-o", "/dev/full", "-output", fmt])
                self.assertEqual((rc, out), (1, b""))
                self.assertEqual(err, b"coal: error: cannot write /dev/full: no space left on device\n")

    def test_output_on_a_full_file_system(self):
        if not self.mounted:
            self.skipTest("cannot mount a tiny tmpfs here (needs root on Linux)")
        out = self.full / "o.fa"
        for fmt, jobs in (("fasta", "1"), ("clustal", "1"), ("fasta", "2")):
            with self.subTest(fmt=fmt, jobs=jobs):
                rc, _, err = self.run_both(["pep.fasta", "nuc.fasta", "-o", str(out), "-output", fmt, "-j", jobs])
                self.assertEqual(rc, 1)
                self.assertEqual(err, f"coal: error: cannot write {out}: no space left on device\n".encode())
                self.assertEqual(sorted(f.name for f in self.full.iterdir()), [])

    def test_full_temp_folder(self):
        if not self.mounted:
            self.skipTest("cannot mount a tiny tmpfs here (needs root on Linux)")
        tmp = {"TMPDIR": str(self.full), "TEMP": str(self.full), "TMP": str(self.full)}
        # standard input is copied to a temporary file first
        rc, _, err = self.run_both(["pep.fasta", "-", "-o", "o.fa"], stdin=(self.d / "nuc.fasta").read_bytes(), **tmp)
        self.assertEqual(rc, 1)
        self.assertEqual(err, f"coal: error: cannot copy standard input to a temporary file in {self.full}: no space "
                              f"left on device\n".encode())
        self.assertFalse((self.d / "o.fa").exists())
        # the alignment for stdout goes through a temporary file
        rc, out, err = self.run_both(["pep.fasta", "nuc.fasta"], **tmp)
        self.assertEqual((rc, out), (1, b""))
        self.assertEqual(err, f"coal: error: cannot write a temporary file in {self.full}: no space left on "
                              f"device\n".encode())


# ---------------------------------------------------------------- 13. smaller items
CONTROL = b"control character: sequence files must be plain text (is this a binary or word-processor file?)"


class TestControlBytes(EngineCase):
    def test_control_bytes_are_an_input_format_error(self):
        cases = {
            "nul_in_protein.fa": (b">s1\nMKV-AL\n>s2\nMKV\x00SAL\n", "n.fa", 4),
            "nul_in_protein_id.fa": (b">s1\x00\nMKV-AL\n>s2\nMKVSAL\n", "n.fa", 1),
            "ctrl_z_at_end.fa": (PEP.encode() + b"\x1a", "n.fa", 5),
            "bell_in_clustal.aln": (b"CLUSTAL\n\ns1 MKV\x07-AL\ns2 MKVSAL\n", "n.fa", 3),
            "utf16.fa": (("\ufeff" + PEP).encode("utf-16-le"), "n.fa", 1),
        }
        for name, (data, nuc, line) in cases.items():
            with self.subTest(file=name):
                (self.d / name).write_bytes(data)
                rc, out, err = self.both([name, nuc])
                self.assertEqual((rc, out), (1, b""))
                self.assertEqual(err, b"coal: error: " + name.encode() + b":" + str(line).encode() + b": " + CONTROL
                                 + b"\n")
        dna = {
            "nul_in_dna.fa": (b">s1\nATGAAAGTT\x00GCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n", 2),
            "nul_after_dna.fa": (NUC.encode() + b"\x00", 5),
            "nul_in_dna_header.fa": (b">s1\x01\nATGAAAGTTGCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n", 1),
            "nul_in_dna_cr.fa": (b">s1\rATGAAAGTTGCTCTG\r>s2\rATGAAGGTC\x00TCAGCACTA\r", 4),
        }
        for name, (data, line) in dna.items():
            with self.subTest(file=name):
                (self.d / name).write_bytes(data)
                rc, out, err = self.both(["p.fa", name])
                self.assertEqual((rc, out), (1, b""))
                self.assertEqual(err, b"coal: error: " + name.encode() + b":" + str(line).encode() + b": " + CONTROL
                                 + b"\n")

    def test_tab_ff_vt_are_still_whitespace(self):
        (self.d / "ws.fa").write_bytes(b">s1\tdesc\x0b\nMKV-AL\x0c\n>s2\nMKVSAL\n")
        rc, out, err = self.both(["ws.fa", "n.fa"])
        self.assertEqual((rc, out, err), (0, EXPECTED.encode(), b""))


class TestLongIds(EngineCase):
    LONG = "L" * 300
    SHORT = "L" * 100 + "..." + "L" * 20

    def test_long_ids_are_shortened_on_stderr_not_in_the_report(self):
        (self.d / "lp.fa").write_text(f">{self.LONG}\nMKV-AL\n>s2\nMKVSAL\n")
        (self.d / "ln.fa").write_text(f">{self.LONG}\nCCCCCCCCC\n>s2\nATGAAGGTCTCAGCACTA\n")
        for name, cmd in engines():  # (--json: the runtime differs between the engines, so one at a time)
            r = self.run_engine(cmd, ["lp.fa", "ln.fa", "-o", "o.fa", "--json"])
            rc, out, err = r.returncode, r.stdout, r.stderr
            self.assertEqual(rc, 1, err)
            self.assertNotIn(self.LONG.encode(), err)
            self.assertEqual(err.splitlines()[0], f"coal: error: {self.SHORT}: protein and DNA are inconsistent (no "
                                                  f"match found)".encode())
            self.assertIn(self.SHORT.encode(), err.splitlines()[-1])  # the summary lines too
            rep = json.loads(out)
            self.assertEqual(rep["sequences"][0]["id"], self.LONG)  # the report keeps the whole ID
            self.assertIn(self.LONG, json.dumps(rep["error"]))
        rc, out, err = self.both(["lp.fa", "ln.fa"])
        self.assertEqual(rc, 1)
        self.assertNotIn(self.LONG.encode(), err)
        self.assertEqual(err.splitlines()[0], f"coal: error: {self.SHORT}: protein and DNA are inconsistent (no match "
                                              f"found)".encode())
        self.assertIn(self.SHORT.encode(), err.splitlines()[-1])  # the summary lines too
        # pairing notices (DNA IDs) and the duplicate-ID error (raised while reading the file)
        rc, out, err = self.both(["lp.fa", "n.fa"])
        self.assertEqual(rc, 0)
        self.assertIn(f"protein '{self.SHORT}' gets DNA 's1'".encode(), err)
        (self.d / "dup.fa").write_text(f">{self.LONG}\nATG\n>{self.LONG}\nATG\n")
        rc, out, err = self.both(["p.fa", "dup.fa"])
        self.assertEqual(err, f"coal: error: dup.fa: duplicate sequence ID '{self.SHORT}': every sequence needs its "
                              f"own ID\n".encode())

    def test_ids_of_200_are_whole(self):
        rid = "I" * 200
        (self.d / "p200.fa").write_text(f">{rid}\nMKV-AL\n>s2\nMKVSAL\n")
        rc, out, err = self.both(["p200.fa", "n.fa"])
        self.assertIn(f"protein '{rid}' gets DNA 's1'".encode(), err)


class TestOutputDash(EngineCase):
    def test_o_dash_is_standard_output(self):
        rc, out, err = self.both(["p.fa", "n.fa", "-o", "-"])
        self.assertEqual((rc, out, err), (0, EXPECTED.encode(), b""))
        self.assertFalse((self.d / "-").exists())
        for extra, msg in ((["--json"], b"--json requires -o (the alignment and the report cannot both go to stdout)"),
                           (["--resume"], b"--resume requires -o (only runs that write to a file can be resumed)")):
            rc, out, err = self.both(["p.fa", "n.fa", "-o", "-", *extra])
            self.assertEqual((rc, out), (2, b""))
            self.assertEqual(err.splitlines()[-1], b"coal: error: " + msg)

    def test_o_dash_on_a_batch_job_line(self):
        (self.d / "jobs.txt").write_text("p.fa n.fa -o -\n")
        r = self.run_engine(PY_CMD, ["--batch", "jobs.txt"])
        self.assertEqual(r.returncode, 2)
        self.assertEqual(r.stderr.splitlines()[-1], b"coal: error: jobs.txt line 1: -o - (standard output) is not "
                                                   b"allowed on a job line: every job of a batch writes to a file")


class TestDnaFileChanged(unittest.TestCase):
    def test_shrunk_dna_file_is_an_error_not_another_gene(self):
        sys.path.insert(0, str(SRC))
        from coal import CoalError
        from coal.seqio import DnaIndex
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "n.fa")
            with open(p, "w") as fh:
                for k in range(200):
                    fh.write(f">s{k}\n" + "ATGAAAGTTGCTCTG" * 10 + "\n")
            ix = DnaIndex([p])
            self.assertTrue(ix.get(0).startswith("ATGAAAGTT"))
            with open(p, "w") as fh:  # the file shrinks after it was indexed
                fh.write(">s1\nATG\n")
            with self.assertRaises(CoalError) as cm:
                ix.get(150)
            self.assertEqual(cm.exception.code, "file_error")
            self.assertEqual(str(cm.exception), f"{p}: could not read the DNA of record 151 (the file changed after "
                                                f"it was read, or cannot be read)")
            ix.close()


class TestEmptyDnaFileAmongSeveral(EngineCase):
    def test_empty_file_is_accepted(self):
        (self.d / "empty.fa").write_bytes(b"")
        for args in (["p.fa", "empty.fa", "n.fa"], ["p.fa", "n.fa", "empty.fa"]):
            rc, out, err = self.both(args)
            self.assertEqual((rc, out, err), (0, EXPECTED.encode(), b""))
        rc, out, err = self.both(["p.fa", "empty.fa"])
        self.assertEqual(rc, 1)
        self.assertIn(b"empty.fa: no sequences found", err)


if __name__ == "__main__":
    unittest.main()
