"""Reproductions of the findings of docs/stress-test-2026-10-23.md, as expected failures: each turns into an
"unexpected success" (remove the decorator) once the behaviour is fixed. Fast (< 2 s) unless COAL_SLOW is set."""
from __future__ import annotations

import os
import random
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CPP = os.environ.get("COAL_CPP_BIN") or str(ROOT / "build" / "coal-cpp")
HAVE_CPP = os.path.isfile(CPP) or os.path.isfile(CPP + ".exe")
ENV = dict(os.environ, COAL_ENGINE="python", PYTHONPATH=str(ROOT / "src"))


def run(engine, args, cwd):
    if engine == "python":
        cmd = [sys.executable, "-m", "coal", *args]
    else:
        cmd = [CPP, *args]
    return subprocess.run(cmd, cwd=cwd, capture_output=True, env=ENV, timeout=600)


def mixed_endings_case(d: Path, first_record_bytes: int):
    """Two genes; the DNA's first record is `first_record_bytes` long with LF line ends, the second record uses lone
    CR (classic Mac) line ends. PAL2NAL v14 reads this file; COAL decides the line-end style from the first 64 kB."""
    rng = random.Random(1)
    g1, g2 = "ATG" + "GCT" * 9, "ATG" + "TGG" * 9
    utr = "".join(rng.choice("ACGT") for _ in range(first_record_bytes))
    (d / "pep.fasta").write_bytes(b">a\nMAAAAAAAAA\n>b\nMWWWWWWWWW\n")
    (d / "nuc.fasta").write_bytes(b">a\n" + (g1 + utr).encode() + b"\n>b\r" + g2.encode() + b"\r")
    return b">a\n" + g1.encode() + b"\n>b\n" + g2.encode() + b"\n"


class MixedLineEndings(unittest.TestCase):
    """Stress-test finding: a lone-CR record after the first 64 kB of a DNA file was not seen as a record (fixed:
    a lone CR ends a line wherever it is, in both engines)."""

    def check(self, engine, first_record_bytes):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            expected = mixed_endings_case(d, first_record_bytes)
            p = run(engine, ["pep.fasta", "nuc.fasta", "-output", "fasta"], tmp)
            self.assertEqual((p.returncode, p.stdout), (0, expected), p.stderr[:300])

    def test_lone_cr_record_within_64kb_python(self):
        self.check("python", 1000)

    @unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
    def test_lone_cr_record_within_64kb_cpp(self):
        self.check("cpp", 1000)

    def test_lone_cr_record_after_64kb_python(self):
        self.check("python", 70_000)

    @unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
    def test_lone_cr_record_after_64kb_cpp(self):
        self.check("cpp", 70_000)

    def test_cr_inside_a_record_and_crlf_mixed(self):
        """LF, CRLF and lone CR in one file, the switch inside a record's data: every record is seen."""
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "pep.fasta").write_bytes(b">a\nMAAA\n>b\nMWWW\n>c\nMKKK\n")
            (d / "nuc.fasta").write_bytes(b">a\r\nATGGCT\rGCTGCT\r>b\nATG\r\nTGGTGGTGG\r>c\rATGAAAAAAAAA")
            want = b">a\nATGGCTGCTGCT\n>b\nATGTGGTGGTGG\n>c\nATGAAAAAAAAA\n"
            for engine in ("python", "cpp") if HAVE_CPP else ("python",):
                p = run(engine, ["pep.fasta", "nuc.fasta", "-output", "fasta"], tmp)
                self.assertEqual((engine, p.returncode, p.stdout), (engine, 0, want), p.stderr[:300])


def repeat_case(d: Path, n: int):
    """A poly-A gene of n residues ending in W, placed after n more codons of poly-A in its DNA (an exact fit, found
    after a long low-complexity stretch)."""
    rng = random.Random(2)
    ala = ["GCA", "GCC", "GCG", "GCT"]
    (d / "pep.fasta").write_text(">r\n" + "A" * (n - 1) + "W\n")
    (d / "nuc.fasta").write_text(">r\n" + "".join(rng.choice(ala) for _ in range(n)) + "".join(rng.choice(ala) for _ in range(n - 1)) + "TGG\n")


@unittest.skipUnless(os.environ.get("COAL_SLOW"), "set COAL_SLOW=1 (about 5 s)")
@unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
class QuadraticExactFitScan(unittest.TestCase):
    """Stress-test finding: coal-cpp's exact-fit search tried every start and verified each one codon by codon, so a
    low-complexity gene after a low-complexity stretch cost time quadratic in its length (fixed: a string search in
    the translated frames, linear). 2x the length must cost well under 4x the time, and the answer must be the
    Python engine's."""

    def test_time_grows_linearly_with_gene_length(self):
        times, outs = {}, {}
        for n in (100_000, 200_000):
            with tempfile.TemporaryDirectory() as tmp:
                repeat_case(Path(tmp), n)
                t = time.monotonic()
                p = run("cpp", ["pep.fasta", "nuc.fasta", "-output", "fasta"], tmp)
                times[n] = time.monotonic() - t
                self.assertEqual(p.returncode, 0, p.stderr[:300])
                q = run("python", ["pep.fasta", "nuc.fasta", "-output", "fasta"], tmp)
                self.assertEqual(p.stdout, q.stdout)
        ratio = times[200_000] / max(times[100_000], 0.05)
        self.assertLess(ratio, 3.0, f"2x the gene length took {ratio:.1f}x the time ({times})")
        self.assertLess(times[200_000], 10.0, times)


class OutOfMemoryMessage(unittest.TestCase):
    """Stress-test finding: the main process running out of memory was a raw MemoryError traceback (seen with the
    1,000,000-sequence input in a 256 MB Windows job). It is one `coal: error: out of memory ...` line, exit 1."""

    def test_memory_error_is_one_error_line(self):
        import io
        from unittest import mock
        from coal import cli
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            mixed_endings_case(d, 100)
            err = io.StringIO()
            with mock.patch.object(cli, "run", side_effect=MemoryError()), mock.patch("sys.stderr", err):
                rc = cli.main([str(d / "pep.fasta"), str(d / "nuc.fasta"), "-o", str(d / "out.fa")])
        self.assertEqual(rc, 1)
        lines = err.getvalue().splitlines()
        self.assertTrue(lines and lines[0].startswith("coal: error: out of memory"), lines)
        self.assertNotIn("Traceback", err.getvalue())
        self.assertFalse((d / "out.fa").exists())


def many_files_case(d: Path, nfiles: int, per_file: int = 2):
    """nfiles DNA files, per_file genes each, round robin over the protein order; returns the expected output."""
    rng = random.Random(3)
    ala = ["GCA", "GCC", "GCG", "GCT"]
    n = nfiles * per_file
    genes = ["ATG" + "".join(rng.choice(ala) for _ in range(5)) for _ in range(n)]
    (d / "pep.fasta").write_text("".join(f">s{k}\nMAAAAA\n" for k in range(n)))
    for f in range(nfiles):
        (d / f"nuc_{f:05d}.fasta").write_text("".join(f">s{k}\n{genes[k]}\n" for k in range(f, n, nfiles)))
    return "".join(f">s{k}\n{genes[k]}\n" for k in range(n)).encode()


class ManyDnaFiles(unittest.TestCase):
    """Stress-test finding: both engines kept every DNA file open, so a thousand files hit the open-file limit (and
    cost 129 KB each in the Python engine). Fixed: an open-file budget with least-recently-used closing. Here the
    hard limit is 160 descriptors and 600 files are read round robin, in both engines."""

    @unittest.skipUnless(hasattr(os, "fork"), "POSIX only (resource limits)")
    def check(self, engine):
        import resource
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            expected = many_files_case(d, 600)
            nucs = sorted(x.name for x in d.glob("nuc_*.fasta"))
            cmd = ([sys.executable, "-m", "coal"] if engine == "python" else [CPP]) + ["pep.fasta", *nucs, "-output", "fasta"]
            p = subprocess.run(cmd, cwd=tmp, capture_output=True, env=dict(ENV, COAL_MIN_NT_PER_JOB="0"), timeout=600,
                               preexec_fn=lambda: resource.setrlimit(resource.RLIMIT_NOFILE, (160, 160)))
            self.assertEqual((p.returncode, p.stdout), (0, expected), p.stderr[-400:])

    def test_python_engine(self):
        self.check("python")

    @unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
    def test_cpp_engine(self):
        self.check("cpp")


class ConcurrentSameOutput(unittest.TestCase):
    """Stress-test finding: two runs writing the same -o file at the same time shared one <out>.coal-partial, so a
    run could exit 0 with an empty or another run's output. Fixed: <out>.coal-lock, one run at a time per output
    file; a lock left by a dead process is replaced. The first run is held on standard input while the second one
    tries."""

    def engines(self):
        return ("python", "cpp") if HAVE_CPP else ("python",)

    def cmd(self, engine, args):
        return ([sys.executable, "-m", "coal"] if engine == "python" else [CPP]) + args

    def test_second_run_is_refused_while_the_first_holds_the_lock(self):
        # the first run is held by the engines' test hook after its first checkpoint (one sequence per batch)
        hold = dict(ENV, COAL_TEST_HOOKS="1", COAL_TEST_PAUSE_AFTER_CHECKPOINT="1", COAL_CHECKPOINT_EVERY="0",
                    COAL_BATCH_NT="1", COAL_MIN_NT_PER_JOB="0")
        for engine in self.engines():
            with tempfile.TemporaryDirectory() as tmp:
                d = Path(tmp)
                expected = mixed_endings_case(d, 100)
                (d / "nuc.fasta").write_bytes((d / "nuc.fasta").read_bytes().replace(b"\r", b"\n"))
                p1 = subprocess.Popen(self.cmd(engine, ["pep.fasta", "nuc.fasta", "-o", "same.fa", "-j", "1"]), cwd=tmp,
                                      env=hold, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                try:
                    line = p1.stderr.readline()
                    self.assertIn(b"paused after checkpoint 1", line, (engine, line))
                    self.assertTrue((d / "same.fa.coal-lock").exists(), engine)
                    p2 = subprocess.run(self.cmd(engine, ["pep.fasta", "nuc.fasta", "-o", "same.fa"]), cwd=tmp, env=ENV,
                                        capture_output=True, timeout=60)
                    self.assertEqual(p2.returncode, 1, (engine, p2.stderr[:300]))
                    self.assertTrue(p2.stderr.startswith(
                        f"coal: error: another COAL run (process {p1.pid}) is writing same.fa".encode()),
                        (engine, p2.stderr[:300]))
                    self.assertFalse((d / "same.fa").exists(), engine)
                finally:
                    p1.terminate()
                    p1.communicate(timeout=60)
                self.assertFalse((d / "same.fa.coal-lock").exists(), (engine, "the lock must go with the run"))
                p3 = subprocess.run(self.cmd(engine, ["pep.fasta", "nuc.fasta", "-o", "same.fa"]), cwd=tmp, env=ENV,
                                    capture_output=True, timeout=60)
                self.assertEqual((p3.returncode, p3.stderr), (0, b""), engine)
                self.assertEqual((d / "same.fa").read_bytes(), expected, engine)
                self.assertFalse((d / "same.fa.coal-lock").exists(), engine)

    def test_stale_lock_of_a_dead_process_is_replaced(self):
        dead = subprocess.Popen([sys.executable, "-c", "pass"])
        dead.wait()
        for engine in self.engines():
            with tempfile.TemporaryDirectory() as tmp:
                d = Path(tmp)
                expected = mixed_endings_case(d, 100)
                (d / "nuc.fasta").write_bytes((d / "nuc.fasta").read_bytes().replace(b"\r", b"\n"))
                (d / "same.fa.coal-lock").write_text(f"{dead.pid}\n")
                p = subprocess.run(self.cmd(engine, ["pep.fasta", "nuc.fasta", "-o", "same.fa"]), cwd=tmp, env=ENV,
                                   capture_output=True, timeout=60)
                self.assertEqual((p.returncode, p.stderr), (0, b""), engine)
                self.assertEqual((d / "same.fa").read_bytes(), expected, engine)
                self.assertFalse((d / "same.fa.coal-lock").exists(), engine)


@unittest.skipUnless(hasattr(signal, "SIGHUP"), "POSIX only")
class SighupUnderNohup(unittest.TestCase):
    """Stress-test finding: `nohup coal ...` died with exit 129 on SIGHUP, because both engines installed their SIGHUP
    handler even when the signal was ignored (what nohup does). The engine must keep an ignored SIGHUP ignored.
    Deterministic: the DNA comes from standard input, which is held open until after the signal."""

    def check(self, engine):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            expected = mixed_endings_case(d, 100)
            nuc = (d / "nuc.fasta").read_bytes().replace(b"\r", b"\n")
            cmd = ([sys.executable, "-m", "coal"] if engine == "python" else [CPP]) + ["pep.fasta", "-", "-output", "fasta"]
            p = subprocess.Popen(cmd, cwd=tmp, env=ENV, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, preexec_fn=lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN))
            time.sleep(0.5)  # the engine is now blocked reading standard input
            os.kill(p.pid, signal.SIGHUP)
            time.sleep(0.2)
            out, err = p.communicate(nuc, timeout=60)
            self.assertEqual((p.returncode, out), (0, expected), err[:300])

    def test_python_engine_keeps_ignored_sighup(self):
        self.check("python")

    @unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
    def test_cpp_engine_keeps_ignored_sighup(self):
        self.check("cpp")


if __name__ == "__main__":
    unittest.main()
