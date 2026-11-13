"""The output lock `<out>.coal-lock` (coal.pipeline.take_output_lock; used by `coal` with -o and by `coal pipeline`).

Independent review 2026-11-13, finding 3: (a) several runs that started together and found a stale or empty lock
could all take it; (b) a symbolic link to the folder had a lock of its own; (c) a lock written on another machine
was "stale" because its process id is not alive here. Nothing but Python is needed.
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
from coal import pipeline as P  # noqa: E402
from coal.core import CoalError  # noqa: E402


def _starter(out, start_at, q):
    while time.time() < start_at:
        pass
    try:
        lock = P.take_output_lock(out, out)
    except CoalError:
        q.put(("refused", os.getpid(), ""))
        return
    time.sleep(0.4)
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
            ps = [ctx.Process(target=_starter, args=(self.out, at, q)) for _ in range(n)]
            for p in ps:
                p.start()
            ev = [q.get(timeout=60) for _ in ps]
            for p in ps:
                p.join()
            got = [e for e in ev if e[0] == "got"]
            self.assertEqual(len(got), 1, f"trial {trial}: {len(got)} runs held the output at once: {ev}")
            self.assertEqual(got[0][2], str(got[0][1]), f"trial {trial}: the holder's lock was replaced under it")
            self.assertFalse(self.lock.exists())

    def test_simultaneous_starts_on_a_stale_lock_have_one_holder(self):
        self.simultaneous(lambda: self.lock.write_text(f"{dead_pid()}\n", encoding="ascii"), 12)

    def test_simultaneous_starts_on_an_empty_lock_have_one_holder(self):
        self.simultaneous(lambda: self.lock.write_text("", encoding="ascii"), 3)


if __name__ == "__main__":
    unittest.main()
