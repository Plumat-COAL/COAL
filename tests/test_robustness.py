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


if __name__ == "__main__":
    unittest.main()
