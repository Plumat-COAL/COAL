"""The output lock `<out>.coal-lock` (coal.pipeline.take_output_lock; used by `coal` with -o and by `coal pipeline`).

Independent review 2026-11-13, finding 3: (a) several runs that started together and found a stale or empty lock
could all take it; (b) a symbolic link to the folder had a lock of its own; (c) a lock written on another machine
was "stale" because its process id is not alive here. Nothing but Python is needed for the class Lock.

CppLock: the same lock in the C++ engine (cpp/src/sysio.cpp take_output_lock), which writes and reads the Python
format: the content, the refusals, the replacement of a stale lock, a race loop of simultaneous `coal-cpp -o`
starts onto a stale lock (200 trials, never two holders), and Python and C++ runs honouring each other's locks.
Needs build/coal-cpp (COAL_CPP_BIN) built with -DCOAL_TEST_HOOKS=ON (tests/cpp_hooks.py); skipped without it,
failed instead when COAL_REQUIRE_CPP is set.

CI, PR #55 (Windows): after simultaneous starts onto a stale lock the guard `<lock>.takeover` stayed behind.
Windows does not remove a file that another process has open, and the runs that lost the race were reading the
guard at that instant; the removal was tried once and its failure ignored. Every removal is now retried
(pipeline._remove, sysio.cpp remove_file), which the tests here drive on every system through a hook that makes
the first removals fail as they do there; a guard that is left all the same is cleared by the next run (its
process is gone, or its age).
"""
from __future__ import annotations

import multiprocessing as mp
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from cpp_hooks import need_test_hooks  # noqa: E402
from coal import pipeline as P  # noqa: E402
from coal.core import CoalError  # noqa: E402


def _starter(out, start_at, q, tried, n):
    while time.time() < start_at:
        pass
    try:
        lock = P.take_output_lock(out, out)
    except CoalError:
        with tried.get_lock():
            tried.value += 1
        q.put(("refused", os.getpid(), ""))
        return
    with tried.get_lock():
        tried.value += 1
    # hold the output until every other start has made its attempt (each of them must meet this lock), not for
    # a fixed time; a second holder would be waiting here too, so the wait ends after a short while in that case
    deadline = time.time() + 20
    while tried.value < n and time.time() < deadline:
        time.sleep(0.005)
    try:
        with open(lock, encoding="ascii") as fh:
            owner = fh.read().split()[0]
    except OSError:
        owner = "missing"
    q.put(("got", os.getpid(), owner))
    P.release_output_lock(lock)


def dead_pid() -> int:
    p = subprocess.Popen([sys.executable, "-c", "pass"])
    p.wait()
    return p.pid


class Lock(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(os.path.realpath(self._td.name))
        self.out = str(self.tmp / "out")
        self.lock = Path(self.out + ".coal-lock")
        self.sleeper = None

    def tearDown(self):
        if self.sleeper is not None:
            self.sleeper.kill()
            self.sleeper.wait()
        self._td.cleanup()

    def live_pid(self) -> int:
        if self.sleeper is None:
            self.sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
        return self.sleeper.pid

    def take(self, out=None):
        lock = P.take_output_lock(out or self.out, "out")
        self.addCleanup(P.release_output_lock, lock)
        return lock

    def test_format_and_release(self):
        lock = self.take()
        self.assertEqual(lock, str(self.lock))
        f = self.lock.read_text(encoding="ascii").split()
        self.assertEqual((int(f[0]), f[1], len(f)), (os.getpid(), P._hostname(), 4))
        self.assertTrue(self.lock.read_text(encoding="ascii").endswith("\n"))
        with self.assertRaisesRegex(CoalError, "this process is already writing out"):
            P.take_output_lock(self.out, "out")       # a second call or thread of this process
        self.assertTrue(self.lock.exists())
        P.release_output_lock(lock)
        self.assertFalse(self.lock.exists())
        self.assertEqual([p.name for p in self.tmp.iterdir()], [], "no temporary or guard file is left")

    def test_lock_in_the_format_of_earlier_versions(self):
        """`<pid>\\n`, as COAL up to 0.5.0 and coal-cpp write it: held while that process lives, replaced after."""
        self.lock.write_text(f"{self.live_pid()}\n", encoding="ascii")
        with self.assertRaisesRegex(CoalError, rf"another COAL run \(process {self.live_pid()}\) is writing out"):
            P.take_output_lock(self.out, "out")
        self.assertEqual(self.lock.read_text(encoding="ascii"), f"{self.live_pid()}\n")
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        self.take()
        self.assertEqual(int(self.lock.read_text(encoding="ascii").split()[0]), os.getpid())

    def test_lock_from_another_machine_is_held(self):
        for pid in (dead_pid(), os.getpid()):      # (our own id on another machine is not us either)
            text = f"{pid} node17.cluster.example 12345 0123abcd\n"
            self.lock.write_text(text, encoding="ascii")
            with self.assertRaises(CoalError) as cm:
                P.take_output_lock(self.out, "out")
            msg = str(cm.exception)
            self.assertIn(f"process {pid} on the machine node17.cluster.example", msg)
            self.assertIn(f"remove {self.lock}", msg)
            self.assertEqual(self.lock.read_text(encoding="ascii"), text, "the other machine's lock was touched")

    def test_our_own_id_from_an_earlier_process_is_stale(self):
        """A container starts every run with the same process id: a lock with our id that this process does not
        hold is a leftover."""
        self.lock.write_text(f"{os.getpid()} {P._hostname()} - 00\n", encoding="ascii")
        self.take()
        self.assertEqual(len(self.lock.read_text(encoding="ascii").split()), 4)

    @unittest.skipUnless(os.path.exists("/proc/self/stat"), "process start times are read from /proc")
    def test_reused_process_id_is_not_the_holder(self):
        pid = self.live_pid()
        real = P._proc_start(pid)
        self.assertNotEqual(real, "-")
        self.lock.write_text(f"{pid} {P._hostname()} {real} aa\n", encoding="ascii")
        with self.assertRaises(CoalError):
            P.take_output_lock(self.out, "out")
        self.lock.write_text(f"{pid} {P._hostname()} {int(real) + 7} aa\n", encoding="ascii")   # another start time
        self.take()

    @unittest.skipIf(os.name == "nt", "symbolic links need a privilege on Windows")
    def test_a_link_to_the_folder_is_the_same_output(self):
        os.mkdir(self.out)
        link = str(self.tmp / "link")
        os.symlink(self.out, link)
        code = ("import sys, time; sys.path.insert(0, sys.argv[1]); from coal.pipeline import take_output_lock; "
                "take_output_lock(sys.argv[2], 'out'); print('held', flush=True); time.sleep(120)")
        other = subprocess.Popen([sys.executable, "-c", code, str(ROOT / "src"), self.out], stdout=subprocess.PIPE)
        try:
            self.assertEqual(other.stdout.readline().strip(), b"held")
            with self.assertRaisesRegex(CoalError, rf"another COAL run \(process {other.pid}\)"):
                P.take_output_lock(link, "link")
            self.assertFalse(os.path.exists(link + ".coal-lock"))
        finally:
            other.kill()
            other.wait()
            other.stdout.close()
        self.assertEqual(self.take(link), str(self.lock))

    def test_a_stale_judgement_cannot_remove_a_newer_lock(self):
        """The race, step by step: run B read the stale lock; run A replaced it and holds the output; B now acts
        on what it read. B must leave A's lock where it is."""
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        stale = P._read_lock(str(self.lock))[0]
        a_lock = f"{self.live_pid()} {P._hostname()} - a1\n"
        self.lock.write_text(a_lock, encoding="ascii")
        P._replace_stale_lock(str(self.lock), stale)
        self.assertEqual(self.lock.read_text(encoding="ascii"), a_lock)
        self.assertEqual(self.names(), ["out.coal-lock"], "the guard of a takeover that found the lock replaced")
        # and while another run is replacing the lock (its guard file exists), nobody else removes anything
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        stale = P._read_lock(str(self.lock))[0]
        guard = Path(str(self.lock) + ".takeover")
        guard.write_text(f"{self.live_pid()} {P._hostname()} -\n", encoding="ascii")
        P._replace_stale_lock(str(self.lock), stale)
        self.assertTrue(self.lock.exists() and guard.exists())
        guard.write_text(f"{dead_pid()} {P._hostname()} -\n", encoding="ascii")    # a guard left by a killed run
        P._replace_stale_lock(str(self.lock), stale)
        self.assertFalse(guard.exists())
        self.take()

    def names(self):
        return sorted(p.name for p in self.tmp.iterdir())

    def patch(self, name, value):
        old = getattr(P, name)
        setattr(P, name, value)
        self.addCleanup(setattr, P, name, old)

    def test_removals_are_tried_again(self):
        """Windows refuses to remove a file that another process has open (a run reading the lock or the guard at
        that instant). The hook makes the next removals fail in that way: lock, guard and temporary files are gone
        all the same."""
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        self.patch("_unlink_fails", 9)
        lock = P.take_output_lock(self.out, "out")
        self.assertEqual(P._unlink_fails, 0, "the removals did not go through the hook")
        self.assertEqual(self.names(), ["out.coal-lock"])
        self.assertEqual(int(self.lock.read_text(encoding="ascii").split()[0]), os.getpid())
        P._unlink_fails = 5
        P.release_output_lock(lock)
        self.assertEqual(self.names(), [])
        self.assertEqual(P._leftovers, set())

    def test_a_guard_that_cannot_be_removed_is_cleared_later(self):
        """The removal fails for as long as it is tried: the guard stays, is remembered, and is removed when the
        lock is released; and a run that finds such a guard of this process (not in a takeover) clears it."""
        self.patch("_REMOVE_TRIES", 3)
        self.addCleanup(P._leftovers.clear)
        guard = Path(str(self.lock) + ".takeover")
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        stale = self.lock.read_bytes()
        self.patch("_in_takeover", lambda path: setattr(P, "_unlink_fails", 10 ** 6))
        P._replace_stale_lock(str(self.lock), stale)
        self.assertEqual(self.names(), ["out.coal-lock", "out.coal-lock.takeover"])
        self.assertEqual(P._leftovers, {str(self.lock), str(guard)})
        f = guard.read_text(encoding="ascii").split()
        self.assertEqual((int(f[0]), f[1], f[2], len(f)), (os.getpid(), P._hostname(), P._proc_start(os.getpid()), 4))
        P._unlink_fails = 0
        self.patch("_in_takeover", None)
        P._leftovers.clear()                      # (as a process that was killed: nothing remembered)
        lock = P.take_output_lock(self.out, "out")   # the guard names this process, which is in no takeover
        self.assertEqual(self.names(), ["out.coal-lock"])
        # remembered files are removed at the release
        guard.write_text("x\n", encoding="ascii")
        P._leftovers.add(str(guard))
        P.release_output_lock(lock)
        self.assertEqual(self.names(), [])
        self.assertEqual(P._leftovers, set())

    def test_the_guard_is_removed_when_the_takeover_fails(self):
        def boom(path):
            raise OSError(5, "input/output error")
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        self.patch("_in_takeover", boom)
        with self.assertRaisesRegex(CoalError, "cannot write out"):
            P.take_output_lock(self.out, "out")
        self.assertEqual(self.names(), ["out.coal-lock"], "the guard was left by the failed takeover")
        with self.assertRaises(KeyboardInterrupt):     # and when the run is interrupted in that instant
            self.patch("_in_takeover", lambda path: (_ for _ in ()).throw(KeyboardInterrupt()))
            P.take_output_lock(self.out, "out")
        self.assertEqual(self.names(), ["out.coal-lock"])
        self.patch("_in_takeover", None)
        self.take()

    def test_leftover_guards_are_judged_like_locks(self):
        """A guard is cleared when its process is gone, or when it is older than the limit (whatever it holds:
        a live process, another machine, nothing); otherwise it is another run's takeover and nothing is touched."""
        guard = Path(str(self.lock) + ".takeover")
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        stale = self.lock.read_bytes()
        long_ago = time.time() - P._GUARD_MAX_AGE - 5
        live = self.live_pid()
        for what, text, cleared_when_new in (
                ("its process is gone", f"{dead_pid()} {P._hostname()} - aa\n", True),
                ("a reused process id", f"{live} {P._hostname()} 1 aa\n", os.path.exists("/proc/self/stat")),
                ("a live process", f"{live} {P._hostname()} - aa\n", False),
                ("another machine", f"{dead_pid()} node17.cluster.example - aa\n", False),
                ("empty", "", False)):
            guard.write_text(text, encoding="ascii")
            P._replace_stale_lock(str(self.lock), stale)
            self.assertEqual(guard.exists(), not cleared_when_new, what)
            self.assertEqual(self.lock.read_bytes(), stale, what)     # this call only cleared (or waited)
            guard.write_text(text, encoding="ascii")
            os.utime(guard, (long_ago, long_ago))
            P._replace_stale_lock(str(self.lock), stale)
            self.assertFalse(guard.exists(), f"{what}: an old guard was not cleared")
            self.assertEqual(self.names(), ["out.coal-lock"], what)
        guard.write_text(f"{live} {P._hostname()} - aa\n", encoding="ascii")     # a run that waits gets through
        os.utime(guard, (long_ago, long_ago))
        self.take()
        self.assertEqual(self.names(), ["out.coal-lock"])

    def test_a_guard_replaced_while_it_was_judged_is_put_back(self):
        """Run C judges a leftover guard; before C moves it aside, run D clears it and starts its own takeover.
        What C then moves is D's guard: C puts it back."""
        guard = Path(str(self.lock) + ".takeover")
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        guard.write_text(f"{dead_pid()} {P._hostname()} - aa\n", encoding="ascii")
        theirs = f"{self.live_pid()} {P._hostname()} - d4\n"
        judge = P._lock_holder

        def judged(pid, host, start):
            who = judge(pid, host, start)
            guard.write_text(theirs, encoding="ascii")      # D, in the meantime
            return who
        self.patch("_lock_holder", judged)
        P._replace_stale_lock(str(self.lock), self.lock.read_bytes())
        self.assertEqual(guard.read_text(encoding="ascii"), theirs)
        self.assertEqual(self.names(), ["out.coal-lock", "out.coal-lock.takeover"])

    def test_host_name(self):
        """One spelling for both engines: printable ASCII as it is, other ASCII "_", anything else "?" per
        character; compared exactly, on Windows without regard to case."""
        self.assertRegex(P._hostname(), r"^[!-~]{1,255}$")
        if os.name != "nt":
            import socket
            real = socket.gethostname
            self.addCleanup(setattr, socket, "gethostname", real)
            for name, want in (("node 7\t.example", "node_7_.example"), ("r\u00e9seau-\u6771\u4eac", "r?seau-??"),
                               ("", "-"), ("x" * 300, "x" * 255), (os.fsdecode(b"a\xffb\x80"), "a?b")):
                socket.gethostname = lambda name=name: name
                self.assertEqual(P._hostname(), want)
            socket.gethostname = real
        self.assertTrue(P._same_host("node7", "node7"))
        self.assertFalse(P._same_host("node7", "node8"))
        self.assertEqual(P._same_host("Node7", "node7"), os.name == "nt")
        other_case = P._hostname().swapcase()
        if other_case != P._hostname():
            self.lock.write_text(f"{self.live_pid()} {other_case} - aa\n", encoding="ascii")
            want = r"another COAL run \(process \d+\) is" if os.name == "nt" else "on the machine"
            with self.assertRaisesRegex(CoalError, want):
                P.take_output_lock(self.out, "out")

    def test_release_leaves_a_lock_that_is_not_ours(self):
        lock = P.take_output_lock(self.out, "out")
        theirs = f"{self.live_pid()} {P._hostname()} - b2\n"
        self.lock.write_text(theirs, encoding="ascii")     # removed by hand, taken by another run
        P.release_output_lock(lock)
        self.assertEqual(self.lock.read_text(encoding="ascii"), theirs)

    def simultaneous(self, prepare, trials: int, n: int = 6) -> None:
        ctx = mp.get_context("fork" if "fork" in mp.get_all_start_methods() else "spawn")
        for trial in range(trials):
            self.lock.unlink(missing_ok=True)
            prepare()
            q = ctx.Queue()
            at = time.time() + (0.3 if ctx.get_start_method() == "fork" else 3.0)
            tried = ctx.Value("i", 0)   # starts that have made their attempt
            ps = [ctx.Process(target=_starter, args=(self.out, at, q, tried, n)) for _ in range(n)]
            for p in ps:
                p.start()
            ev = [q.get(timeout=60) for _ in ps]
            for p in ps:
                p.join()
            got = [e for e in ev if e[0] == "got"]
            self.assertEqual(len(got), 1, f"trial {trial}: {len(got)} runs held the output at once: {ev}")
            self.assertEqual(got[0][2], str(got[0][1]), f"trial {trial}: the holder's lock was replaced under it")
            self.assertEqual(self.names(), [], f"trial {trial}: lock, guard or temporary files were left")

    def test_simultaneous_starts_on_a_stale_lock_have_one_holder(self):
        self.simultaneous(lambda: self.lock.write_text(f"{dead_pid()}\n", encoding="ascii"), 12)

    def test_simultaneous_starts_on_an_empty_lock_have_one_holder(self):
        self.simultaneous(lambda: self.lock.write_text("", encoding="ascii"), 3)


def cpp_binary() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / "coal-cpp.exe"):
        if cand.is_file():
            return cand
    return None


PY_HOLDER = ("import os, sys, time; sys.path.insert(0, sys.argv[1]); from coal import pipeline as P\n"
             "from coal.core import CoalError\n"
             "while time.time() < float(sys.argv[4]): pass\n"
             "try:\n    lock = P.take_output_lock(sys.argv[2], 'out')\n"
             "except CoalError as e:\n    sys.stderr.write(str(e)); sys.exit(1)\n"
             "t = time.time()\n"
             "while not os.path.exists(sys.argv[3]) and time.time() - t < 60: time.sleep(0.005)\n"
             "mine = open(lock, 'rb').read() == P._held_locks[lock]\n"
             "P.release_output_lock(lock); sys.exit(0 if mine else 3)\n")


class CppLock(unittest.TestCase):
    EXAMPLE = [str(ROOT / "examples" / "clean" / "protein.aln.fasta"), str(ROOT / "examples" / "clean" / "dna.fasta")]

    def setUp(self):
        self.cpp = cpp_binary()
        if self.cpp is None:
            if os.environ.get("COAL_REQUIRE_CPP"):
                self.fail("coal-cpp not built (COAL_REQUIRE_CPP is set)")
            self.skipTest("coal-cpp not built (cmake --build build)")
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(os.path.realpath(self._td.name))
        self.out = self.tmp / "out.fa"
        self.lock = Path(str(self.out) + ".coal-lock")
        self.go = self.tmp / "go"
        self.procs = []
        self.env = {k: v for k, v in os.environ.items() if not k.startswith("COAL_TEST_")}
        self.start_at = None
        self.want = subprocess.run([str(self.cpp), *self.EXAMPLE], capture_output=True, text=True, env=self.env).stdout
        self.assertTrue(self.want.startswith(">"))

    def tearDown(self):
        for p in self.procs:
            if p.poll() is None:
                p.kill()
            p.wait()
            for fh in (p.stdout, p.stderr):
                if fh is not None:
                    fh.close()
        self._td.cleanup()

    def run_cpp(self, out=None, **hooks):
        if hooks:
            need_test_hooks(self, self.cpp)
        return subprocess.run([str(self.cpp), *self.EXAMPLE, "-o", str(out or self.out)], capture_output=True,
                              text=True, env={**self.env, **({"COAL_TEST_HOOKS": "1", **hooks} if hooks else {})},
                              timeout=120)

    def names(self):
        return sorted(x.name for x in self.tmp.iterdir() if ".coal-lock" in x.name)

    def start_cpp(self, out=None, hold=True):
        """A coal-cpp run that, once it holds the lock, waits for the file self.go (the engine's test hook)."""
        need_test_hooks(self, self.cpp)
        env = {**self.env, "COAL_TEST_HOOKS": "1"}
        if hold:
            env["COAL_TEST_HOLD_LOCK_UNTIL"] = str(self.go)
        if self.start_at is not None:     # the race loop: every run reaches the lock at the same instant
            env["COAL_TEST_LOCK_START_AT"] = str(int(self.start_at * 1000))
        p = subprocess.Popen([str(self.cpp), *self.EXAMPLE, "-o", str(out or self.out)], env=env,
                             stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        self.procs.append(p)
        return p

    def start_py(self):
        p = subprocess.Popen([sys.executable, "-c", PY_HOLDER, str(ROOT / "src"), str(self.out), str(self.go),
                              str(self.start_at or 0)],
                             stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        self.procs.append(p)
        return p

    def wait_for_lock(self, pid, timeout=30.0):
        t = time.time()
        while time.time() - t < timeout:
            try:
                f = self.lock.read_text(encoding="ascii").split()
            except OSError:
                f = []
            if f and f[0] == str(pid):
                return f
            time.sleep(0.005)
        self.fail(f"process {pid} did not take {self.lock}")

    def live_pid(self) -> int:
        p = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
        self.procs.append(p)
        return p.pid

    def test_format_is_the_python_engines(self):
        p = self.start_cpp()
        f = self.wait_for_lock(p.pid)
        raw = self.lock.read_bytes()
        self.assertTrue(raw.endswith(b"\n") and b"\r" not in raw, raw)
        self.assertEqual(len(f), 4, f)
        self.assertEqual((int(f[0]), f[1], f[2]), (p.pid, P._hostname(), P._proc_start(p.pid)))
        self.assertRegex(f[3], r"^[0-9a-f]{16}$")
        self.assertEqual(P._lock_holder(int(f[0]), f[1], f[2]), "process")     # as the Python engine judges it
        self.assertEqual(sorted(x.name for x in self.tmp.iterdir() if x.name.startswith("out.fa.coal-lock")),
                         ["out.fa.coal-lock"])                                  # no temporary or guard file left
        self.go.write_text("")
        self.assertEqual(p.wait(timeout=60), 0, p.stderr.read())
        self.assertFalse(self.lock.exists())
        self.assertEqual(self.out.read_text(), self.want)

    def test_python_and_cpp_honour_each_others_locks(self):
        # C++ holds: Python and a second C++ run are refused, and the holder's lock stays as it is
        p = self.start_cpp()
        self.wait_for_lock(p.pid)
        before = self.lock.read_bytes()
        with self.assertRaisesRegex(CoalError, rf"another COAL run \(process {p.pid}\) is writing out"):
            P.take_output_lock(str(self.out), "out")
        r = self.run_cpp()
        self.assertEqual(r.returncode, 1)
        self.assertIn(f"another COAL run (process {p.pid}) is writing {self.out}", r.stderr)
        self.assertIn(f"remove {self.lock}", r.stderr)
        self.assertEqual(self.lock.read_bytes(), before)
        self.go.write_text("")
        self.assertEqual(p.wait(timeout=60), 0, p.stderr.read())
        # Python holds (this process): C++ is refused; released: C++ runs
        lock = P.take_output_lock(str(self.out), "out")
        try:
            mine = self.lock.read_bytes()
            r = self.run_cpp()
            self.assertEqual(r.returncode, 1)
            self.assertIn(f"another COAL run (process {os.getpid()}) is writing", r.stderr)
            self.assertEqual(self.lock.read_bytes(), mine)
        finally:
            P.release_output_lock(lock)
        r = self.run_cpp()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertFalse(self.lock.exists())

    def test_lock_from_another_machine_is_held(self):
        for pid in (self.live_pid(), dead_pid()):
            text = f"{pid} node17.cluster.example - 0123456789abcdef\n"
            self.lock.write_text(text, encoding="ascii")
            r = self.run_cpp()
            self.assertEqual(r.returncode, 1)
            self.assertIn(f"process {pid} on the machine node17.cluster.example", r.stderr)
            self.assertIn(f"remove {self.lock}", r.stderr)
            self.assertEqual(self.lock.read_text(encoding="ascii"), text, "the other machine's lock was touched")
            self.assertFalse(self.out.exists())

    def test_stale_locks_are_replaced(self):
        live = self.live_pid()
        cases = {"an older lock (pid only), process gone": f"{dead_pid()}\n",
                 "this machine, process gone": f"{dead_pid()} {P._hostname()} - aa\n",
                 "not a lock, after 2 s": "what is this\n"}
        if os.path.exists("/proc/self/stat"):
            cases["a reused process id (another start time)"] = \
                f"{live} {P._hostname()} {int(P._proc_start(live)) + 7} aa\n"
        for what, text in cases.items():
            self.lock.write_text(text, encoding="ascii")
            r = self.run_cpp()
            self.assertEqual(r.returncode, 0, (what, r.stderr))
            self.assertFalse(self.lock.exists(), what)
            self.assertEqual(self.out.read_text(), self.want, what)
        # and the locks that are held: an older lock of a live process; a full one with its real start time
        held = [f"{live}\n", f"{live} {P._hostname()} {P._proc_start(live)} aa\n", f"{live} {P._hostname()} - aa\n"]
        for text in held:
            self.lock.write_text(text, encoding="ascii")
            r = self.run_cpp()
            self.assertEqual(r.returncode, 1, text)
            self.assertIn(f"another COAL run (process {live})", r.stderr)
            self.assertEqual(self.lock.read_text(encoding="ascii"), text)
        self.lock.unlink()

    def test_a_stale_judgement_cannot_remove_a_newer_lock(self):
        """The race, step by step: run B judges the lock stale and is held up (the engine's test hook); run A
        replaces the stale lock and holds the output; B goes on with what it read. B must leave A's lock where it
        is and be refused. (With an unguarded takeover B removed A's lock and both ran.)"""
        need_test_hooks(self, self.cpp)
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        b_go = self.go.with_name(self.go.name + ".b")
        b = subprocess.Popen([str(self.cpp), *self.EXAMPLE, "-o", str(self.out)], stdout=subprocess.DEVNULL,
                             stderr=subprocess.PIPE, text=True,
                             env={**self.env, "COAL_TEST_HOOKS": "1",
                                  "COAL_TEST_LOCK_BEFORE_TAKEOVER_UNTIL": str(b_go),
                                  "COAL_TEST_HOLD_LOCK_UNTIL": str(self.go)})
        self.procs.append(b)
        waiting = Path(str(b_go) + ".waiting")   # B says when it has read the stale lock and waits
        t0 = time.time()
        while not waiting.exists():
            self.assertIsNone(b.poll(), "B ended before it reached the takeover")
            self.assertLess(time.time() - t0, 60, "B did not reach the takeover")
            time.sleep(0.005)
        a = self.start_cpp()
        self.wait_for_lock(a.pid)
        mine = self.lock.read_bytes()
        b_go.write_text("")                    # A holds the output: B goes on with what it read
        self.assertEqual(b.wait(timeout=60), 1, "the second run also took the output")
        self.assertIn(f"another COAL run (process {a.pid})", b.stderr.read())
        self.assertEqual(self.lock.read_bytes(), mine)
        self.assertIsNone(a.poll())
        self.go.write_text("")
        self.assertEqual(a.wait(timeout=60), 0, a.stderr.read())
        self.assertFalse(self.lock.exists())

    def test_guard_of_a_takeover(self):
        """While another run is replacing the stale lock (its guard exists) nothing is removed; a guard left by a
        killed run is cleared and the stale lock replaced."""
        guard = Path(str(self.lock) + ".takeover")
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        guard.write_text(f"{dead_pid()} {P._hostname()} -\n", encoding="ascii")
        r = self.run_cpp()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertFalse(guard.exists() or self.lock.exists())
        stale = f"{dead_pid()}\n"
        self.lock.write_text(stale, encoding="ascii")
        guard.write_text(f"{self.live_pid()} {P._hostname()} -\n", encoding="ascii")
        p = self.start_cpp(hold=False)
        time.sleep(1.0)
        self.assertIsNone(p.poll(), "the run went on although another run's takeover guard was there")
        self.assertEqual(self.lock.read_text(encoding="ascii"), stale)
        guard.unlink()                                  # the other run's takeover ends
        self.assertEqual(p.wait(timeout=60), 0, p.stderr.read())
        self.assertEqual(self.names(), [])

    def test_leftover_guards_are_judged_like_locks(self):
        """A guard whose process is gone is cleared at once; one that names a live process, another machine or
        nothing is another run's takeover until it is older than the limit, then it is cleared too."""
        guard = Path(str(self.lock) + ".takeover")
        long_ago = time.time() - P._GUARD_MAX_AGE - 5
        live = self.live_pid()
        for what, text in (("a live process", f"{live} {P._hostname()} - aa\n"),
                           ("another machine", f"{dead_pid()} node17.cluster.example - aa\n"), ("empty", "")):
            self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
            guard.write_text(text, encoding="ascii")
            os.utime(guard, (long_ago, long_ago))
            t = time.time()
            r = self.run_cpp()
            self.assertEqual(r.returncode, 0, (what, r.stderr))
            self.assertLess(time.time() - t, 10, what)
            self.assertEqual(self.names(), [], what)
            self.assertEqual(self.out.read_text(), self.want, what)

    def test_removals_are_tried_again(self):
        """Windows refuses to remove a file that another process has open. The engine's hook makes its first
        removals fail in that way (the temporary name, the stale lock, the guard, its own lock): all are gone."""
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        r = self.run_cpp(COAL_TEST_UNLINK_FAILS="60")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.names(), [])
        self.assertEqual(self.out.read_text(), self.want)
        # and a removal that fails for as long as it is tried (2 s): tried a last time when the run ends
        self.lock.write_text(f"{dead_pid()}\n", encoding="ascii")
        r = self.run_cpp(COAL_TEST_UNLINK_FAILS="225")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.names(), [])

    def test_host_name_is_the_python_engines(self):
        """The same name in another case: this machine on Windows (computer names are not case sensitive there),
        another machine elsewhere; both engines agree."""
        other_case = P._hostname().swapcase()
        if other_case == P._hostname():
            self.skipTest("this machine's name has no letters")
        live = self.live_pid()
        self.lock.write_text(f"{live} {other_case} - aa\n", encoding="ascii")
        r = self.run_cpp()
        self.assertEqual(r.returncode, 1)
        same = os.name == "nt"
        self.assertEqual(f"another COAL run (process {live}) is writing" in r.stderr, same, r.stderr)
        self.assertEqual("on the machine" in r.stderr, not same, r.stderr)
        self.assertEqual(P._lock_holder(live, other_case, "-"), "process" if same else "host")

    @unittest.skipIf(os.name == "nt", "symbolic links need a privilege on Windows")
    def test_a_link_to_the_folder_is_the_same_output(self):
        real = self.tmp / "real"
        real.mkdir()
        os.symlink(real, self.tmp / "link")
        self.out = real / "out.fa"
        self.lock = Path(str(self.out) + ".coal-lock")
        p = self.start_cpp()
        self.wait_for_lock(p.pid)
        r = self.run_cpp(self.tmp / "link" / "out.fa")
        self.assertEqual(r.returncode, 1)
        self.assertIn(f"another COAL run (process {p.pid})", r.stderr)
        with self.assertRaises(CoalError):
            P.take_output_lock(str(self.tmp / "link" / "out.fa"), "out")
        self.go.write_text("")
        self.assertEqual(p.wait(timeout=60), 0, p.stderr.read())

    def race(self, trials: int, starters, prepare, lead: float = 0.03) -> dict:
        """`trials` times: prepare the lock, start the runs so that all reach the lock at the same instant (`lead`
        seconds from now: the engine's test hook), wait until all but the holders have ended, count the holders
        (each waits, holding the lock, for self.go), let them finish. Exactly one may hold."""
        seen = {"trials": 0, "holders": 0, "double": 0, "refused": 0}
        for trial in range(trials):
            for f in (self.lock, self.go, self.out):
                f.unlink(missing_ok=True)
            prepare()
            self.start_at = time.time() + lead
            ps = [start() for start in starters]
            self.start_at = None
            t = time.time()
            while sum(p.poll() is None for p in ps) > 1 and time.time() - t < 20:
                time.sleep(0.002)
            time.sleep(0.02)          # a second holder would still be there; a run about to be refused is not
            holders = [p for p in ps if p.poll() is None]
            self.go.write_text("")
            codes = [p.wait(timeout=60) for p in ps]
            errs = [p.stderr.read() for p in ps]
            for p in ps:
                p.stderr.close()
                self.procs.remove(p)
            seen["trials"] += 1
            seen["holders"] += len(holders)
            seen["double"] += len(holders) > 1
            seen["refused"] += codes.count(1)
            self.assertEqual(len(holders), 1, f"trial {trial}: {len(holders)} runs held the output at once: {errs}")
            self.assertEqual(sorted(codes), [0] + [1] * (len(ps) - 1), f"trial {trial}: {codes} {errs}")
            for c, e in zip(codes, errs):
                if c == 1:
                    self.assertIn("another COAL run (process", e, f"trial {trial}")
            self.assertEqual(self.names(), [], f"trial {trial}: lock, guard or temporary files were left")
        return seen

    def test_simultaneous_cpp_starts_on_a_stale_lock_have_one_holder(self):
        trials = int(os.environ.get("COAL_LOCK_RACE_TRIALS", "200"))
        stale = f"{dead_pid()}\n"
        seen = self.race(trials, [self.start_cpp] * 6, lambda: self.lock.write_text(stale, encoding="ascii"))
        self.assertEqual((seen["trials"], seen["holders"], seen["double"]), (trials, trials, 0))
        self.assertEqual(seen["refused"], 5 * trials)

    def test_simultaneous_python_and_cpp_starts_have_one_holder(self):
        stale = f"{dead_pid()} {P._hostname()} - aa\n"
        seen = self.race(15, [self.start_py, self.start_py, self.start_cpp, self.start_cpp],
                         lambda: self.lock.write_text(stale, encoding="ascii"), lead=0.5)
        self.assertEqual((seen["holders"], seen["double"]), (15, 0))

    def test_simultaneous_starts_while_removals_are_refused_leave_nothing(self):
        """The Windows failure of PR #55 on every system: runs of both engines start together onto a stale lock
        while each run's first removals fail as they do when another run has the file open."""
        stale = f"{dead_pid()} {P._hostname()} - aa\n"

        def cpp():
            self.env["COAL_TEST_UNLINK_FAILS"] = "12"
            try:
                return self.start_cpp()
            finally:
                del self.env["COAL_TEST_UNLINK_FAILS"]
        seen = self.race(10, [cpp] * 5, lambda: self.lock.write_text(stale, encoding="ascii"))
        self.assertEqual((seen["holders"], seen["double"]), (10, 0))

    def test_simultaneous_cpp_starts_without_a_lock_have_one_holder(self):
        seen = self.race(20, [self.start_cpp] * 6, lambda: None)
        self.assertEqual((seen["holders"], seen["double"]), (20, 0))


if __name__ == "__main__":
    unittest.main()
