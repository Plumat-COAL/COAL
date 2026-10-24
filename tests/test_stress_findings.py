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
                    # (on Windows a venv's python.exe is a launcher whose child holds the lock: the id differs)
                    who = f"process {p1.pid}" if os.name != "nt" else "process "
                    self.assertTrue(p2.stderr.startswith(f"coal: error: another COAL run ({who}".encode())
                                    and b") is writing same.fa" in p2.stderr, (engine, p2.stderr[:300]))
                    self.assertFalse((d / "same.fa").exists(), engine)
                finally:
                    p1.terminate()
                    p1.communicate(timeout=60)
                if os.name != "nt":  # SIGTERM is handled and releases the lock; TerminateProcess cannot (stale lock)
                    self.assertFalse((d / "same.fa.coal-lock").exists(), (engine, "the lock must go with the run"))
                # a stale lock (Windows) or none: the next run goes ahead quietly
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


def wrong_code_gene(d: Path, n: int):
    """One gene of n residues built with code 2 (TGA = Trp, ATA = Met), run under code 1: a few mismatches."""
    rng = random.Random(5)
    aas = [("W", "TGA"), ("M", "ATA"), ("A", "GCT"), ("K", "AAA"), ("L", "CTG"), ("D", "GAC")]
    picks = [rng.choice(aas) for _ in range(n - 1)]
    (d / "pep.fasta").write_text(">g\nM" + "".join(a for a, _ in picks) + "\n")
    (d / "nuc.fasta").write_text(">g\nATG" + "".join(c for _, c in picks) + "\n")


class HintSkipsEnormousGenes(unittest.TestCase):
    """Stress-test finding: the genetic-code hint re-threaded a 1,000,000-residue gene under 17 codes (28 s of a 30 s
    run). Genes over 100,000 residues are left out of the hint in both engines; shorter ones still get it, and the
    two engines' stderr stays identical."""

    def run_both(self, args, tmp):
        outs = []
        for engine in ("python", "cpp") if HAVE_CPP else ("python",):
            p = run(engine, args, tmp)
            self.assertEqual(p.returncode, 0, (engine, p.stderr[-300:]))
            outs.append(p.stderr)
        self.assertEqual(len(set(outs)), 1, outs)
        return outs[0]

    def test_short_gene_gets_the_hint_and_a_huge_one_does_not(self):
        with tempfile.TemporaryDirectory() as tmp:
            wrong_code_gene(Path(tmp), 1_000)
            err = self.run_both(["pep.fasta", "nuc.fasta", "-o", "out.fa"], tmp)
            self.assertIn(b"wrong genetic code?", err)
        with tempfile.TemporaryDirectory() as tmp:
            wrong_code_gene(Path(tmp), 100_001)
            t = time.monotonic()
            err = self.run_both(["pep.fasta", "nuc.fasta", "-o", "out.fa"], tmp)
            self.assertNotIn(b"genetic code", err)
            self.assertIn(b"does not correspond to", err)
            self.assertLess(time.monotonic() - t, 60)


class FileListsInMessages(unittest.TestCase):
    """Stress-test finding: with thousands of DNA files the duplicate-ID and no-sequences messages listed every path."""

    def test_at_most_three_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "pep.fasta").write_text(">a\nMK\n>b\nMK\n")
            for k in range(5):
                (d / f"n{k}.fa").write_text(f">a\nATGAAA\n" if k in (0, 4) else f">x{k}\nATGAAA\n")
            nucs = [f"n{k}.fa" for k in range(5)]
            for engine in ("python", "cpp") if HAVE_CPP else ("python",):
                p = run(engine, ["pep.fasta", *nucs], tmp)
                self.assertEqual(p.returncode, 1, engine)
                self.assertTrue(p.stderr.startswith(b"coal: error: n0.fa, n1.fa, n2.fa and 2 more files: duplicate sequence ID 'a'"),
                                (engine, p.stderr[:200]))


class QuietKeepsNoEntries(unittest.TestCase):
    """Stress-test finding: a million mismatching genes kept every per-codon message in memory (about 150 KB per
    gene of unrelated DNA). With -q and no report nothing reads them, so the run report keeps none; the output and
    exit status are unchanged, and both engines agree."""

    def test_report_object_keeps_no_entries(self):
        """Since the memory improvements of 2026-10-24 no mode keeps the entries in memory: one byte of kind and one
        of flags per sequence, plus the first flagged entry for the summary line (ReportKeepsCounts has the rest)."""
        from coal.report import RunReport
        warns = [f"pepAlnPos {c}: A does not correspond to TTT" for c in range(1, 401)]
        quiet = RunReport([], "p", ["n"], None, "fasta", 1, False, False, keep_entries=False, hint=False)
        loud = RunReport([], "p", ["n"], None, "fasta", 1, False, False, keep_entries=False, hint=True)
        for k in range(50):
            e1 = quiet.add_sequence(k, f"s{k}", warns, [f"s{k}: protein and DNA are probably unrelated: 400 of 400"])
            e2 = loud.add_sequence(k, f"s{k}", warns, [f"s{k}: protein and DNA are probably unrelated: 400 of 400"])
            self.assertEqual(e1, e2)  # the --resume entries file gets the same line
        for rec in (quiet, loud):
            self.assertEqual(len(rec._store), 0)
            self.assertEqual((len(rec.kinds), len(rec.flags), len(rec.codons)), (50, 50, 150))
            self.assertEqual(bytes(rec.kinds), b"m" * 50)
        self.assertEqual(bytes(quiet.kinds), bytes(loud.kinds))
        self.assertEqual(quiet.counts(), loud.counts())
        self.assertEqual(len(loud.summary_lines()), 1)

    def test_quiet_run_writes_no_resume_entries(self):
        """With -q and no report the --resume entries file stays empty (a million unrelated genes wrote 55 GB of
        it); with the hint wanted it is written, in both engines."""
        rng = random.Random(11)
        hold = dict(ENV, COAL_TEST_HOOKS="1", COAL_TEST_PAUSE_AFTER_CHECKPOINT="1", COAL_CHECKPOINT_EVERY="0",
                    COAL_BATCH_NT="1", COAL_MIN_NT_PER_JOB="0")
        for engine in ("python", "cpp") if HAVE_CPP else ("python",):
            for quiet in (True, False):
                with tempfile.TemporaryDirectory() as tmp:
                    d = Path(tmp)
                    n = 4
                    aa = {"GCT": "A", "AAA": "K", "CTG": "L", "GAC": "D"}
                    genes = [[rng.choice(list(aa)) for _ in range(30)] for _ in range(n)]
                    peps = ["".join(aa[c] for c in g) for g in genes]
                    peps[0] = "W" + peps[0][1:]  # one mismatch: a clean sequence has no entry in either mode
                    (d / "pep.fasta").write_text("".join(f">s{k}\n{pep}\n" for k, pep in enumerate(peps)))
                    (d / "nuc.fasta").write_text("".join(f">s{k}\n" + "".join(g) + "\n" for k, g in enumerate(genes)))
                    cmd = ([sys.executable, "-m", "coal"] if engine == "python" else [CPP]) + \
                          ["pep.fasta", "nuc.fasta", "-o", "out.fa", "-j", "1"] + (["-q"] if quiet else [])
                    p = subprocess.Popen(cmd, cwd=tmp, env=hold, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                    try:
                        for _ in range(8):  # the mismatch warning first (default mode), then the pause
                            line = p.stderr.readline()
                            if b"paused after checkpoint 1" in line:
                                break
                        else:
                            self.fail((engine, quiet, line))
                        size = (d / "out.fa.coal-entries.jsonl").stat().st_size
                        self.assertEqual(size == 0, quiet, (engine, quiet, size))
                    finally:
                        p.terminate()
                        p.communicate(timeout=60)

    def test_quiet_output_unchanged_on_unrelated_dna(self):
        rng = random.Random(9)
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            n = 200
            genes = ["".join(rng.choice(["GCT", "AAA", "CTG", "GAC"]) for _ in range(30)) for _ in range(n)]
            (d / "pep.fasta").write_text("".join(f">s{k}\nM" + "W" * 29 + "\n" for k in range(n)))
            (d / "nuc.fasta").write_text("".join(f">s{k}\n{genes[k]}\n" for k in range(n)))
            outs = {}
            for engine in ("python", "cpp") if HAVE_CPP else ("python",):
                p = run(engine, ["pep.fasta", "nuc.fasta", "-output", "fasta"], tmp)
                q = run(engine, ["pep.fasta", "nuc.fasta", "-output", "fasta", "-q"], tmp)
                self.assertEqual((p.returncode, q.returncode, q.stderr), (0, 0, b""), engine)
                self.assertEqual(p.stdout, q.stdout, engine)
                self.assertIn(b"probably unrelated", p.stderr, engine)
                outs[engine] = (p.stdout, p.stderr)
            self.assertEqual(len(set(outs.values())), 1)


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


# ------------------------------------------------------------------ memory improvements of 2026-10-24 (A1 / A2)
# A1: the run report keeps counts and a bounded sample, not every per-codon message; the entries go to a file.
# A2: the Python engine's workers return mismatch columns and codons, not warning strings.
# Goldens (stdout, stderr, exit status, report JSON) were recorded from `main` at fd46a4e, before either change:
# `COAL_WRITE_GOLDENS=1 python -m unittest tests.test_stress_findings.ReportGoldens` rewrites them (Python engine).

GOLDEN_DIR = ROOT / "tests" / "golden_memory"
WRITE_GOLDENS = os.environ.get("COAL_WRITE_GOLDENS") == "1"
VOLATILE = ("runtime_seconds", "engine", "coal_version")  # differ between runs / engines / releases
CODE1 = {"A": "GCT", "K": "AAA", "L": "CTG", "D": "GAC", "V": "GTT", "E": "GAA", "G": "GGC", "P": "CCA", "S": "TCT",
         "T": "ACC", "N": "AAC", "Q": "CAG", "R": "CGT", "I": "ATC", "F": "TTC", "Y": "TAC", "H": "CAC", "C": "TGC",
         "W": "TGG"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]


def _gene(rng, n, aas=CODE1):
    picks = [rng.choice(sorted(aas)) for _ in range(n - 1)]
    return "M" + "".join(picks), "ATG" + "".join(aas[a] for a in picks)


def _fasta(recs):
    return "".join(f">{i}\n{s}\n" for i, s in recs)


def _write(d: Path, peps, nucs, pep_name="pep.fasta"):
    # LF on every platform: the report records the inputs' sha256
    (d / pep_name).write_text(_fasta(peps), newline="\n")
    (d / "nuc.fasta").write_text(_fasta(nucs), newline="\n")


def unrelated_case(d: Path, n=40, ncod=30):
    """Proteins with DNA that is not theirs: every gene is written with per-codon warnings and flagged."""
    rng = random.Random(101)
    peps, nucs = [], []
    for k in range(n):
        p, _ = _gene(rng, ncod)
        peps.append((f"g{k}", p))
        nucs.append((f"g{k}", "".join(rng.choice(SENSE) for _ in range(ncod))))
    _write(d, peps, nucs)


def wrongcode_case(d: Path, n=30, ncod=40):
    """Genes built with the vertebrate mitochondrial code (TGA = W, ATA = M), run under the standard code."""
    rng = random.Random(102)
    aas = dict(CODE1, W="TGA", M="ATA")
    peps, nucs = [], []
    for k in range(n):
        p, s = _gene(rng, ncod, aas)
        peps.append((f"g{k}", p))
        nucs.append((f"g{k}", s))
    _write(d, peps, nucs)


def startcodon_case(d: Path, n=8, ncod=40):
    """Clean genes whose first codon is ATA: the only mismatch is at the start codon."""
    rng = random.Random(103)
    peps, nucs = [], []
    for k in range(n):
        p, s = _gene(rng, ncod)
        peps.append((f"g{k}", p))
        nucs.append((f"g{k}", "ATA" + s[3:]))
    _write(d, peps, nucs)


def run_mismatch_case(d: Path, n=6, ncod=40):
    """One gene whose last eight codons are wrong: a dense run of mismatches to the end of the gene."""
    rng = random.Random(104)
    peps, nucs = [], []
    for k in range(n):
        p, s = _gene(rng, ncod)
        if k == 3:
            tail = "".join(CODE1["W" if p[i] != "W" else "A"] for i in range(ncod - 8, ncod))
            s = s[:3 * (ncod - 8)] + tail
        peps.append((f"g{k}", p))
        nucs.append((f"g{k}", s))
    _write(d, peps, nucs)


def ambiguous_case(d: Path, n=6, ncod=30):
    """NNN codons (not checkable), one gene with a selenocysteine residue, one with an X."""
    rng = random.Random(105)
    peps, nucs = [], []
    for k in range(n):
        p, s = _gene(rng, ncod)
        if k in (1, 4):
            s = s[:9] + "NNN" + s[12:21] + "NNN" + s[24:]
        if k == 2:
            p = p[:5] + "U" + p[6:]
        if k == 5:
            p = p[:7] + "X" + p[8:]
        peps.append((f"g{k}", p))
        nucs.append((f"g{k}", s))
    _write(d, peps, nucs)


def broken_case(d: Path, n=9, ncod=40, frameshift=(2,), unrelated=(), short=()):
    """Clean genes with some broken ones: an extra base after codon 20 (frameshift), unrelated DNA, DNA too short."""
    rng = random.Random(106)
    peps, nucs = [], []
    for k in range(n):
        p, s = _gene(rng, ncod)
        if k in frameshift:
            s = s[:60] + "G" + s[60:]
        if k in unrelated:
            s = "".join(rng.choice(SENSE) for _ in range(ncod))
        if k in short:
            s = s[:3 * (ncod - 10)]
        peps.append((f"g{k}", p))
        nucs.append((f"g{k}", s))
    _write(d, peps, nucs)


def partial_ids_case(d: Path, n=6, ncod=30):
    """Two DNA records whose IDs are not protein IDs: pairing by order with a partial ID match."""
    rng = random.Random(107)
    peps, nucs = [], []
    for k in range(n):
        p, s = _gene(rng, ncod)
        if k == 2:
            s = s[:30] + CODE1["W" if p[10] != "W" else "A"] + s[33:]  # one mismatch in a mispaired gene
        peps.append((f"g{k}", p))
        nucs.append((f"g{k}" if k not in (2, 4) else f"other{k}", s))
    _write(d, peps, nucs)


def blockonly_case(d: Path, n=3, ncod=30):
    """A CLUSTAL alignment with '#' marks under half of the columns; one gene with mismatches in and out of them."""
    rng = random.Random(108)
    peps, nucs = [], []
    for k in range(n):
        p, s = _gene(rng, ncod)
        if k == 1:
            for i in (4, 5, 20, 21):  # columns 5, 6 (marked) and 21, 22 (not marked)
                s = s[:3 * i] + CODE1["W" if p[i] != "W" else "A"] + s[3 * i + 3:]
        peps.append((f"g{k}", p))
        nucs.append((f"g{k}", s))
    marks = "#" * (ncod // 2) + " " * (ncod - ncod // 2)
    (d / "pep.aln").write_text("CLUSTAL W (1.83) multiple sequence alignment\n\n\n"
                               + "".join(f"{i:<10}{s}\n" for i, s in peps) + " " * 10 + marks + "\n", newline="\n")
    (d / "nuc.fasta").write_text(_fasta(nucs), newline="\n")


GOLDEN_CASES = {
    # name: (generator, arguments, stdin file, extra environment)
    "unrelated": (unrelated_case, ["pep.fasta", "nuc.fasta", "-o", "out.fa", "--report", "r.json"], None, {}),
    "unrelated_json_allwarn": (unrelated_case, ["pep.fasta", "nuc.fasta", "-o", "out.fa", "--json", "--max-warnings", "0"],
                               None, {}),
    "unrelated_quiet": (unrelated_case, ["pep.fasta", "nuc.fasta", "-q", "-o", "out.fa", "--report", "r.json"], None, {}),
    "wrongcode": (wrongcode_case, ["pep.fasta", "nuc.fasta", "-o", "out.fa", "--report", "r.json"], None, {}),
    "wrongcode_stdout": (wrongcode_case, ["pep.fasta", "nuc.fasta", "--report", "r.json"], None, {}),
    "wrongcode_nomismatch": (wrongcode_case, ["pep.fasta", "nuc.fasta", "-nomismatch", "-o", "out.fa", "--report",
                                              "r.json"], None, {}),
    "wrongcode_clustal_workers": (wrongcode_case, ["pep.fasta", "nuc.fasta", "-output", "clustal", "-j", "2", "-o",
                                                   "out.aln", "--report", "r.json"], None,
                                  {"COAL_MIN_NT_PER_JOB": "0", "COAL_BATCH_NT": "700"}),
    "wrongcode_stdin": (wrongcode_case, ["-", "nuc.fasta", "--report", "r.json"], "pep.fasta", {}),
    "startcodon": (startcodon_case, ["pep.fasta", "nuc.fasta", "-o", "out.fa", "--report", "r.json"], None, {}),
    "run_mismatch": (run_mismatch_case, ["pep.fasta", "nuc.fasta", "-o", "out.fa", "--report", "r.json"], None, {}),
    "ambiguous": (ambiguous_case, ["pep.fasta", "nuc.fasta", "-o", "out.fa", "--report", "r.json"], None, {}),
    "keepgoing": (lambda d: broken_case(d, frameshift=(2,), unrelated=(5,), short=(7,)),
                  ["pep.fasta", "nuc.fasta", "--keep-going", "-o", "out.fa", "--report", "r.json"], None, {}),
    "keepgoing_paml_workers": (lambda d: broken_case(d, frameshift=(2,), unrelated=(5,), short=(7,)),
                               ["pep.fasta", "nuc.fasta", "--keep-going", "-output", "paml", "-j", "2", "-o", "out.paml",
                                "--report", "r.json"], None, {"COAL_MIN_NT_PER_JOB": "0", "COAL_BATCH_NT": "500"}),
    "aborted": (lambda d: broken_case(d, frameshift=(4,)), ["pep.fasta", "nuc.fasta", "-o", "out.fa", "--report",
                                                           "r.json"], None, {}),
    "aborted_after_warnings": (lambda d: broken_case(d, frameshift=(6,), unrelated=(1,)),
                               ["pep.fasta", "nuc.fasta", "-o", "out.fa", "--report", "r.json"], None, {}),
    "partial_ids": (partial_ids_case, ["pep.fasta", "nuc.fasta", "-o", "out.fa", "--report", "r.json"], None, {}),
    "blockonly": (blockonly_case, ["pep.aln", "nuc.fasta", "-blockonly", "-o", "out.fa", "--report", "r.json"], None, {}),
    "blockonly_codon": (blockonly_case, ["pep.aln", "nuc.fasta", "-output", "codon", "-o", "out.codon", "--report",
                                         "r.json"], None, {}),
}


def _normalised_report(text: bytes):
    import json
    rep = json.loads(text)
    for k in VOLATILE:
        rep.pop(k, None)
    return rep


def run_case(engine, name, tmp):
    """Run one golden case; returns (rc, stdout, stderr, report bytes or None)."""
    gen, args, stdin_name, extra = GOLDEN_CASES[name]
    d = Path(tmp)
    gen(d)
    cmd = ([sys.executable, "-m", "coal"] if engine == "python" else [CPP]) + args
    stdin = (d / stdin_name).read_bytes() if stdin_name else None
    p = subprocess.run(cmd, cwd=tmp, capture_output=True, env=dict(ENV, **extra), input=stdin, timeout=300)
    rep = (d / "r.json").read_bytes() if (d / "r.json").exists() else None
    # stderr is a text stream: CRLF on Windows in both engines (the alignment on stdout is written in binary mode)
    return p.returncode, p.stdout, p.stderr.replace(b"\r\n", b"\n"), rep


class ReportGoldens(unittest.TestCase):
    """The stdout, stderr, exit status and report JSON of warning-heavy runs (mismatches, flagged genes, failures,
    the hint, -nomismatch, -blockonly, workers, --keep-going, an aborted run, standard input) are byte-identical to
    what `main` produced before the run report stopped keeping every message, in both engines."""

    def check(self, engine):
        for name in GOLDEN_CASES:
            with self.subTest(case=name, engine=engine), tempfile.TemporaryDirectory() as tmp:
                rc, out, err, rep = run_case(engine, name, tmp)
                g = GOLDEN_DIR / name
                if WRITE_GOLDENS and engine == "python":
                    g.mkdir(parents=True, exist_ok=True)
                    (g / "rc").write_text(f"{rc}\n", newline="\n")
                    (g / "stdout").write_bytes(out)
                    (g / "stderr").write_bytes(err)
                    if rep is not None:
                        (g / "report.json").write_bytes(rep)
                    continue
                self.assertEqual(rc, int((g / "rc").read_text()), err[-500:])
                self.assertEqual(err, (g / "stderr").read_bytes())
                if "--json" in GOLDEN_CASES[name][1]:
                    self.assertEqual(_normalised_report(out), _normalised_report((g / "stdout").read_bytes()))
                else:
                    self.assertEqual(out, (g / "stdout").read_bytes())
                if (g / "report.json").exists():
                    self.assertIsNotNone(rep, "no report written")
                    self.assertEqual(_normalised_report(rep), _normalised_report((g / "report.json").read_bytes()))

    def test_python_engine(self):
        self.check("python")

    @unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
    def test_cpp_engine(self):
        self.check("cpp")


class ResumeStaysIdentical(unittest.TestCase):
    """A run killed after its fifth checkpoint and resumed gives the same output, summary lines and report (apart
    from `resumed_from`) as an uninterrupted run, in both engines: the report's per-sequence state is rebuilt from
    the entries file."""

    def check(self, engine):
        hold = dict(ENV, COAL_TEST_HOOKS="1", COAL_TEST_PAUSE_AFTER_CHECKPOINT="5", COAL_CHECKPOINT_EVERY="0",
                    COAL_BATCH_NT="1", COAL_MIN_NT_PER_JOB="0")
        base = [sys.executable, "-m", "coal"] if engine == "python" else [CPP]
        args = ["pep.fasta", "nuc.fasta", "-o", "out.fa", "--report", "r.json", "-j", "1"]
        with tempfile.TemporaryDirectory() as ref, tempfile.TemporaryDirectory() as tmp:
            for d in (ref, tmp):
                broken_case(Path(d), n=24, frameshift=(), unrelated=(3, 9), short=())
            p = subprocess.run(base + args, cwd=ref, capture_output=True, env=ENV, timeout=120)
            self.assertEqual(p.returncode, 0, (engine, p.stderr[-300:]))
            k = subprocess.Popen(base + args, cwd=tmp, env=hold, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            try:
                for _ in range(400):
                    line = k.stderr.readline()
                    if b"paused after checkpoint 5" in line or not line:
                        break
                self.assertIn(b"paused after checkpoint 5", line, engine)
                k.kill()
            finally:
                k.communicate(timeout=60)
            self.assertTrue((Path(tmp) / "out.fa.coal-checkpoint.json").exists(), engine)
            r = subprocess.run(base + args + ["--resume"], cwd=tmp, capture_output=True, env=ENV, timeout=120)
            self.assertEqual(r.returncode, 0, (engine, r.stderr[-300:]))
            self.assertTrue(r.stderr.startswith(b"coal: notice: resuming out.fa: "), (engine, r.stderr[:200]))
            self.assertEqual((Path(tmp) / "out.fa").read_bytes(), (Path(ref) / "out.fa").read_bytes(), engine)

            def summary(e):
                return [ln for ln in e.splitlines() if ln.startswith(b"coal: summary:")]
            self.assertEqual(summary(r.stderr), summary(p.stderr), engine)
            a = _normalised_report((Path(tmp) / "r.json").read_bytes())
            b = _normalised_report((Path(ref) / "r.json").read_bytes())
            self.assertIsNotNone(a.pop("resumed_from", None), engine)
            self.assertEqual(a.pop("command_line"), b.pop("command_line") + " --resume", engine)
            self.assertEqual(a, b, engine)
            self.assertEqual(sorted(x.name for x in Path(tmp).iterdir()),
                             ["nuc.fasta", "out.fa", "pep.fasta", "r.json"], engine)

    def test_python_engine(self):
        self.check("python")

    @unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
    def test_cpp_engine(self):
        self.check("cpp")


def _peak_rss_of(cmd, cwd, env):
    """Run a command and sample the resident set size of its process tree (psutil); returns (rc, peak bytes)."""
    import psutil
    p = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    proc, peak = psutil.Process(p.pid), 0
    while p.poll() is None:
        try:
            peak = max(peak, sum(q.memory_info().rss for q in [proc, *proc.children(recursive=True)]))
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
        time.sleep(0.002)
    return p.returncode, peak


def clean_case(d: Path, n, ncod):
    """`n` clean genes of `ncod` codons: the same sizes as unrelated_case(d, n, ncod), without the warnings."""
    rng = random.Random(101)
    peps, nucs = [], []
    for k in range(n):
        p, s = _gene(rng, ncod)
        peps.append((f"g{k}", p))
        nucs.append((f"g{k}", s))
    _write(d, peps, nucs)


class ReportKeepsCounts(unittest.TestCase):
    """A1: in default mode (summary lines and hint wanted, no report) the run report's memory no longer grows with
    the number of per-codon warnings, in both engines, with small work batches so that only the report's own
    per-sequence state can grow. Before the change the Python engine kept about 150 KB and coal-cpp about 55 KB per
    400-codon gene of unrelated DNA (every codon a warning) until the end of the run: here about 30 KB / 8 KB per
    gene of 50 codons, several times the whole clean run."""

    N, NCOD = 1500, 50

    def test_python_engine_default_mode_memory_is_flat(self):
        """Peak memory grows by well under 1 KB per gene of unrelated DNA (the IDs and file positions, about 350 B
        per sequence): 2,000 more genes used to cost about 60 MB."""
        import tracemalloc
        from coal import pipeline
        from coal.report import RunReport
        peaks = {}
        saved = pipeline.BATCH_NT
        pipeline.BATCH_NT = 20_000  # about 130 genes per batch
        try:
            for n in (1000, 3000):
                with tempfile.TemporaryDirectory() as tmp:
                    d = Path(tmp)
                    unrelated_case(d, n, self.NCOD)
                    pep, nuc, out = str(d / "pep.fasta"), str(d / "nuc.fasta"), str(d / "out.fa")
                    rec = RunReport([], pep, [nuc], out, "fasta", 1, False, False, keep_entries=False, hint=True)
                    tracemalloc.start()
                    pipeline.run(pep, [nuc], out_path=out, jobs=1, warn=lambda w: None, notice=lambda t: None,
                                 note=lambda t: None, rec=rec)
                    lines = rec.summary_lines()
                    peaks[n] = tracemalloc.get_traced_memory()[1]
                    tracemalloc.stop()
                    rec.close()
                    self.assertTrue(lines and lines[0].startswith(f"{n:,} sequences of {n:,} were written"), lines[:1])
                    self.assertEqual(sorted(x.name for x in d.iterdir()), ["nuc.fasta", "out.fa", "pep.fasta"])
        finally:
            pipeline.BATCH_NT = saved
        self.assertLess((peaks[3000] - peaks[1000]) / 2000, 1000, peaks)

    @unittest.skipUnless(HAVE_CPP, "build/coal-cpp not built")
    def test_cpp_engine_default_mode_memory_is_flat(self):
        try:
            import psutil  # noqa: F401
        except ImportError:
            self.skipTest("psutil not installed")
        peaks = {}
        env = dict(ENV, COAL_BATCH_NT="20000")
        for name, gen in (("clean", clean_case), ("unrelated", unrelated_case)):
            with tempfile.TemporaryDirectory() as tmp:
                d = Path(tmp)
                gen(d, 2 * self.N, self.NCOD)
                rc, peaks[name] = _peak_rss_of([CPP, "pep.fasta", "nuc.fasta", "-o", "out.fa", "-j", "1"], tmp, env)
                self.assertEqual(rc, 0)
                self.assertEqual(sorted(x.name for x in d.iterdir()), ["nuc.fasta", "out.fa", "pep.fasta"])
        self.assertLess(peaks["unrelated"] / peaks["clean"], 1.5, peaks)

    def test_report_from_resume_entries_file_equals_temp_file_report(self):
        """A report read back from the --resume entries file next to -o equals one read back from the temporary
        entries file of a run to standard output; both files are gone once the report is closed."""
        import contextlib
        import io
        from coal.pipeline import run
        from coal.report import RunReport
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            broken_case(d, n=12, frameshift=(), unrelated=(2, 7))
            pep, nuc = str(d / "pep.fasta"), str(d / "nuc.fasta")
            reps, paths = [], []
            for out in (str(d / "a.fa"), None):
                rec = RunReport([], pep, [nuc], out, "fasta", 1, False, False)
                sink = io.StringIO()
                with contextlib.redirect_stdout(sink):
                    run(pep, [nuc], out_path=out, jobs=1, warn=lambda w: None, notice=lambda t: None,
                        note=lambda t: None, rec=rec)
                if out is None:
                    self.assertEqual(sink.getvalue(), (d / "a.fa").read_text())
                paths.append(rec.entries_file)
                self.assertTrue(paths[-1] and os.path.exists(paths[-1]), paths[-1])
                reps.append(rec.to_dict())
                rec.close()
                self.assertFalse(os.path.exists(paths[-1]))
            self.assertEqual(os.path.basename(paths[0]), "a.fa.coal-entries.jsonl")
            self.assertEqual(len(reps[0]["sequences"]), 12)
            self.assertEqual([s["status"] for s in reps[0]["sequences"]].count("failed"), 2)
            for r in reps:  # the same report apart from the output
                r["output"] = r["command"] = r["runtime_seconds"] = None
                r["methods_text"] = r["methods_text"].split(" Command: ")[0]
            self.assertEqual(reps[0], reps[1])
            self.assertEqual(sorted(x.name for x in d.iterdir()), ["a.fa", "nuc.fasta", "pep.fasta"])


class WorkersReturnColumnsNotStrings(unittest.TestCase):
    """A2: a worker's result holds the mismatch columns and codons of each sequence, not the warning strings (about
    35 KB per 400-codon gene of unrelated DNA, held per batch in every worker and in the main process); the main
    process rebuilds exactly thread_one's strings."""

    def test_batch_result_is_compact_and_rebuilds_the_same_strings(self):
        import pickle
        from coal import pipeline
        from coal.core import mismatch_warnings, thread_one
        rng = random.Random(7)
        genes = []
        for k in range(10):
            p, _ = _gene(rng, 400)
            genes.append((f"g{k}", p, "".join(rng.choice(SENSE) for _ in range(400))))
        pipeline._init_worker(1, [3] * 400, False, False)
        text, items, mcols, bad = pipeline._work(genes)
        self.assertEqual(len(items), 10)
        self.assertLess(len(pickle.dumps(items)) / 10, 6_000, "about 2 KB per gene expected, 35 KB as strings")
        for (pid, prow, pdna), (iid, row, other, cols, codons, notes, failure) in zip(genes, items):
            self.assertEqual((iid, row, failure), (pid, prow, None))
            self.assertGreater(len(cols), 300)
            self.assertEqual(len(codons), 3 * len(cols))
            it = thread_one(pid, prow, pdna, 1)
            self.assertEqual(list(other) + mismatch_warnings(row, cols, codons), it.warnings)
            self.assertEqual(it.warnings[0], f"pepAlnPos {cols[0]}: {row[cols[0] - 1]} does not correspond to {codons[:3]}")
        self.assertTrue(text.startswith(">g0\n"))


@unittest.skipUnless(os.environ.get("COAL_SLOW"), "set COAL_SLOW=1 (about 2 minutes)")
class UnrelatedTwentyThousand(unittest.TestCase):
    """20,000 genes of 400 codons of unrelated DNA in default mode (the stress test's `unrelated_20k`): peak memory
    within about twice the clean-data figures, in both engines (before: coal-cpp 1,149 MB, Python 3,150 MB)."""

    def test_peak_memory(self):
        try:
            import psutil  # noqa: F401
        except ImportError:
            self.skipTest("psutil not installed")
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            subprocess.run([sys.executable, str(ROOT / "tests" / "tools" / "stress_gen.py"), "--nseq", "20000",
                            "--ncodons", "400", "--mismatch", "unrelated", tmp], check=True, timeout=600)
            limits = {"python": 1_000_000_000, "cpp": 150_000_000}
            for engine in ("python", "cpp") if HAVE_CPP else ("python",):
                cmd = ([sys.executable, "-m", "coal"] if engine == "python" else [CPP]) + \
                      ["pep.fasta", "nuc.fasta", "-o", f"out.{engine}", "-j", "1"]
                rc, peak = _peak_rss_of(cmd, tmp, ENV)
                self.assertEqual(rc, 0, engine)
                self.assertLess(peak, limits[engine], (engine, peak))
                (d / f"out.{engine}").unlink()


if __name__ == "__main__":
    unittest.main()
