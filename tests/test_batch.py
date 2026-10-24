"""Batch mode (`coal --batch JOBS`, GitHub issue #15): every job equals a single run byte for byte, defaults and
per-line overrides, --batch-out-dir, the failure policy and exit codes, up-front JOBS checks with line numbers, the
batch report, stderr prefixes, and what coal-cpp and the dispatcher do with --batch.

Run:  python -m unittest tests.test_batch   (the coal-cpp checks need COAL_CPP_BIN or ./build/coal-cpp; else skipped)
"""
from __future__ import annotations

import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from coal import batch  # noqa: E402
from coal.cli import main as coal_main  # noqa: E402

FIX = ROOT / "tests" / "fixtures"
CLEAN = (FIX / "clean" / "pep.fasta", FIX / "clean" / "nuc.fasta")
TABLE2 = (FIX / "codontable2" / "pep.fasta", FIX / "codontable2" / "nuc.fasta")
BROKEN = (FIX / "codon_mismatch" / "pep.fasta", FIX / "codon_mismatch" / "nuc.fasta")  # exit 1; --keep-going: 2
WARN = (FIX / "long_mismatch" / "pep.fasta", FIX / "long_mismatch" / "nuc.fasta")  # exit 0 with a mismatch warning
CLUSTAL = (FIX / "clustal_input" / "pep.aln", FIX / "clustal_input" / "nuc.fasta")
VOLATILE = ("runtime_seconds", "command_line")  # report fields that differ between a batch job and a single run


def cpp_binary() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe")):
        if cand.is_file():
            return cand
    return None


def run_py(args, cwd, stdin: str | None = None):
    """coal.cli.main in this process (as the tests of single runs do); (exit code, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    old, old_in = os.getcwd(), sys.stdin
    os.chdir(cwd)
    if stdin is not None:
        sys.stdin = io.StringIO(stdin)
    try:
        with redirect_stdout(out), redirect_stderr(err):
            try:
                rc = coal_main(list(args))
            except SystemExit as e:
                rc = e.code
    finally:
        os.chdir(old)
        sys.stdin = old_in
    return rc, out.getvalue(), err.getvalue()


def strip(rep: dict) -> dict:
    return {k: v for k, v in rep.items() if k not in VOLATILE}


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="coal_batch_"))
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def jobs(self, lines, name="jobs.txt") -> str:
        (self.tmp / name).write_text("".join(ln + "\n" for ln in lines), encoding="utf-8")
        return name

    def read(self, name) -> bytes:
        return (self.tmp / name).read_bytes()


def run_cmd(args, cwd, timeout=300):
    """`python -m coal ...` as a subprocess (the parallel path starts worker processes, which a test process run
    under `-m unittest` should not do in-process); (exit code, stdout, stderr)."""
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), PYTHONIOENCODING="utf-8")
    env.pop("COAL_IMPL", None)
    r = subprocess.run([sys.executable, "-m", "coal", *args], cwd=cwd, env=env, capture_output=True, text=True,
                       encoding="utf-8", timeout=timeout, stdin=subprocess.DEVNULL)
    return r.returncode, r.stdout, r.stderr


class TestBatchParallel(_Tmp):
    """--batch-parallel N: the same outputs and report as the sequential batch, each job's messages in one block, the
    failure policy and exit codes unchanged, the option refused on a job line and by coal-cpp."""
    LINES = [f"{CLEAN[0]} {CLEAN[1]} -o a.fa --report a.json",
             f"{TABLE2[0]} {TABLE2[1]} -codontable 2 -output clustal -o b.aln --report b.json",
             f"{CLUSTAL[0]} {CLUSTAL[1]} -blockonly -output paml -o c.phy --report c.json",
             f"{WARN[0]} {WARN[1]} -o d.fa --report d.json"]

    def test_parallel_batch_equals_the_sequential_batch(self):
        seq, par = self.tmp / "seq", self.tmp / "par"
        for d in (seq, par):
            d.mkdir()
        jobs = self.jobs(self.LINES)
        rc1, out1, err1 = run_cmd(["--batch", str(self.tmp / jobs), "--report", "batch.json"], seq)
        rc2, out2, err2 = run_cmd(["--batch", str(self.tmp / jobs), "--report", "batch.json", "--batch-parallel", "3"], par)
        self.assertEqual((rc1, out1), (0, ""), err1)
        self.assertEqual((rc2, out2), (0, ""), err2)
        for n in ("a.fa", "b.aln", "c.phy", "d.fa"):
            self.assertEqual((par / n).read_bytes(), (seq / n).read_bytes(), n)
        for n in ("a.json", "b.json", "c.json", "d.json"):
            self.assertEqual(strip(json.loads((par / n).read_text())), strip(json.loads((seq / n).read_text())), n)
        b1, b2 = json.loads((seq / "batch.json").read_text()), json.loads((par / "batch.json").read_text())
        self.assertEqual(b2["summary"], b1["summary"])
        self.assertEqual([j["index"] for j in b2["jobs"]], [1, 2, 3, 4])  # job order, whatever order they ended in
        for j1, j2 in zip(b1["jobs"], b2["jobs"]):
            self.assertEqual({k: v for k, v in j2.items() if k != "report"}, {k: v for k, v in j1.items() if k != "report"})
            self.assertEqual(strip(j2["report"]), strip(j1["report"]))
            self.assertEqual(j2["report"]["command_line"], j1["report"]["command_line"])  # no --batch-parallel in it
        self.assertEqual(sorted(err2.splitlines()), sorted(err1.splitlines()))  # the same lines, maybe another order
        self.assertEqual(err2.splitlines()[-1], "coal: batch: 4 jobs, 4 ok, 0 failed")
        # each job's lines are one contiguous block
        seen = []
        for line in err2.splitlines()[:-1]:
            k = line.split("[", 1)[1].split("/", 1)[0]
            if not seen or seen[-1] != k:
                seen.append(k)
        self.assertEqual(len(seen), len(set(seen)), err2)

    def test_failing_job_exit_codes_and_progress(self):
        lines = [f"{CLEAN[0]} {CLEAN[1]} -o a.fa", f"{BROKEN[0]} {BROKEN[1]} -o b.fa", f"{CLEAN[0]} {CLEAN[1]} -o c.fa"]
        rc, out, err = run_cmd(["--batch", self.jobs(lines), "--batch-parallel", "0", "--progress"], self.tmp)
        self.assertEqual(rc, 1, err)
        self.assertTrue((self.tmp / "a.fa").is_file() and (self.tmp / "c.fa").is_file())
        self.assertFalse((self.tmp / "b.fa").exists())
        self.assertEqual(err.splitlines()[-1], "coal: batch: 3 jobs, 2 ok, 1 failed")
        self.assertEqual(sum(1 for ln in err.splitlines() if ln.startswith("coal: progress: job ")), 3)
        self.assertIn("coal: [2/3 ", err)
        # --keep-going: exit 2 and the skipped count, as the sequential batch
        rc, out, err = run_cmd(["--batch", self.jobs(lines), "--batch-parallel", "2", "--keep-going", "--json"], self.tmp)
        self.assertEqual(rc, 2, err)
        rep = json.loads(out)
        self.assertEqual(rep["summary"]["jobs"], 3)
        self.assertEqual(rep["summary"]["failed"], 0)
        rc1, out1, _ = run_cmd(["--batch", self.jobs(lines), "--keep-going", "--json"], self.tmp)
        self.assertEqual(rc1, 2)
        self.assertEqual(rep["summary"], json.loads(out1)["summary"])

    def test_option_checks(self):
        rc, _, err = run_py(["--batch", self.jobs([f"{CLEAN[0]} {CLEAN[1]} -o a.fa --batch-parallel 2"])], self.tmp)
        self.assertEqual(rc, 2)
        self.assertIn("--batch-parallel is not allowed on a job line", err)
        rc, _, err = run_py(["--batch", self.jobs([f"{CLEAN[0]} {CLEAN[1]} -o a.fa"]), "--batch-parallel", "-1"], self.tmp)
        self.assertEqual(rc, 2)
        self.assertIn("--batch-parallel must be 0", err)
        rc, _, err = run_py(["--batch", self.jobs([f"{CLEAN[0]} {CLEAN[1]} -o a.fa"]), "--batch-parallel", "1"], self.tmp)
        self.assertEqual(rc, 0, err)  # 1 is the sequential path

    @unittest.skipIf(os.name == "nt", "POSIX signals")
    def test_sigterm_stops_every_worker_and_leaves_no_partial_output(self):
        # a slow batch (8 jobs of 4,000 synthetic genes, 2 at a time) killed with SIGTERM: exit 143 as a single run,
        # every worker gone, and no job's half-written output or temporary file left behind
        import signal
        import time
        sys.path.insert(0, str(ROOT / "tests" / "tools"))
        gen = ROOT / "tests" / "tools" / "stress_gen.py"
        data = self.tmp / "data"
        subprocess.run([sys.executable, str(gen), "--nseq", "4000", "--ncodons", "200", "--mismatch", "none", str(data)],
                       check=True, capture_output=True, timeout=120)
        lines = [f"{data / 'pep.fasta'} {data / 'nuc.fasta'} -o job{i}.fa" for i in range(8)]
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
        env.pop("COAL_IMPL", None)
        p = subprocess.Popen([sys.executable, "-m", "coal", "--batch", self.jobs(lines), "--batch-parallel", "2", "-q"],
                             cwd=self.tmp, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        t0 = time.monotonic()
        while time.monotonic() - t0 < 20 and not list(self.tmp.glob("job*.fa")):
            time.sleep(0.05)  # the first job has finished: the next ones are in flight
        p.send_signal(signal.SIGTERM)
        try:
            _, err = p.communicate(timeout=30)
        except subprocess.TimeoutExpired:
            p.kill()
            self.fail("the batch did not stop within 30 s of SIGTERM")
        self.assertEqual(p.returncode, 143, err)
        self.assertLess(len(list(self.tmp.glob("job*.fa"))), 8)
        self.assertEqual([x.name for x in self.tmp.iterdir() if ".coal" in x.name or x.name.startswith(".coal_")], [])
        self.assertNotIn("Traceback", err)
        time.sleep(0.5)
        try:
            import psutil
        except ImportError:  # the worker check needs psutil (CI installs it); the rest of the test stands
            return
        self.assertEqual([c for c in psutil.Process().children(recursive=True) if "coal" in " ".join(c.cmdline())], [])

    @unittest.skipUnless(cpp_binary(), "coal-cpp not built")
    def test_cpp_refuses_batch_parallel(self):
        r = subprocess.run([str(cpp_binary()), str(CLEAN[0]), str(CLEAN[1]), "--batch-parallel", "2"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 2)
        self.assertIn("--batch is only available in the Python command", r.stderr)


class TestBatchEqualsSingleRuns(_Tmp):
    def test_three_jobs_byte_identical(self):
        lines = [f"{CLEAN[0]} {CLEAN[1]} -o a.fa --report a.json",
                 f"# a comment, and blank lines, are skipped", "",
                 f"{TABLE2[0]} {TABLE2[1]} -codontable 2 -output clustal -o b.aln --report b.json",
                 f"{CLUSTAL[0]} {CLUSTAL[1]} -blockonly -output paml -o c.phy --report c.json"]
        rc, out, err = run_py(["--batch", self.jobs(lines)], self.tmp)
        self.assertEqual((rc, out), (0, ""), err)
        self.assertEqual(err.strip().splitlines()[-1], "coal: batch: 3 jobs, 3 ok, 0 failed")
        got = {n: self.read(n) for n in ("a.fa", "b.aln", "c.phy")}
        reps = {n: json.loads(self.read(n)) for n in ("a.json", "b.json", "c.json")}
        single = self.tmp / "single"
        single.mkdir()
        for args in ([str(CLEAN[0]), str(CLEAN[1]), "-o", "a.fa", "--report", "a.json"],
                     [str(TABLE2[0]), str(TABLE2[1]), "-codontable", "2", "-output", "clustal", "-o", "b.aln",
                      "--report", "b.json"],
                     [str(CLUSTAL[0]), str(CLUSTAL[1]), "-blockonly", "-output", "paml", "-o", "c.phy",
                      "--report", "c.json"]):
            self.assertEqual(run_py(args, single)[0], 0)
        for n, data in got.items():
            self.assertEqual(data, (single / n).read_bytes(), n)
        for n, rep in reps.items():
            one = json.loads((single / n).read_text())
            self.assertEqual(strip(rep), strip(one), n)
            self.assertEqual(rep["command_line"], one["command_line"])  # the job line alone is the whole command

    def test_defaults_and_overrides(self):
        # -codontable 2 and -output clustal from the main command line; the third job overrides the format
        lines = [f"{TABLE2[0]} {TABLE2[1]} -o a.aln", f"{TABLE2[0]} {TABLE2[1]} -o b.aln -codontable 2",
                 f"{TABLE2[0]} {TABLE2[1]} -o c.fa -output fasta"]
        rc, out, err = run_py(["-codontable", "2", "-output", "clustal", "--batch", self.jobs(lines), "--report",
                               "batch.json"], self.tmp)
        self.assertEqual(rc, 0, err)
        single = self.tmp / "single"
        single.mkdir()
        run_py([str(TABLE2[0]), str(TABLE2[1]), "-codontable", "2", "-output", "clustal", "-o", "a.aln"], single)
        run_py([str(TABLE2[0]), str(TABLE2[1]), "-codontable", "2", "-o", "c.fa"], single)
        self.assertEqual(self.read("a.aln"), (single / "a.aln").read_bytes())
        self.assertEqual(self.read("b.aln"), (single / "a.aln").read_bytes())
        self.assertEqual(self.read("c.fa"), (single / "c.fa").read_bytes())
        rep = json.loads(self.read("batch.json"))
        self.assertEqual([j["report"]["options"]["codontable"] for j in rep["jobs"]], [2, 2, 2])
        self.assertEqual([j["report"]["options"]["output_format"] for j in rep["jobs"]],
                         ["clustal", "clustal", "fasta"])
        # the per-job command line = the defaults (without the batch-level options), then the job line
        import shlex  # the report quotes as the single run's report does (shlex, on every platform)
        self.assertEqual(rep["jobs"][2]["report"]["command_line"],
                         " ".join(shlex.quote(a) for a in ["coal", "-codontable", "2", "-output", "clustal",
                                                           str(TABLE2[0]), str(TABLE2[1]), "-o", "c.fa",
                                                           "-output", "fasta"]))

    def test_quiet_and_nomismatch_defaults_reach_the_jobs(self):
        lines = [f"{WARN[0]} {WARN[1]} -o a.fa"]
        rc, _, err = run_py(["--batch", self.jobs(lines)], self.tmp)
        self.assertEqual(rc, 0)
        self.assertIn("] warning: ", err)
        self.assertIn("] summary: ", err)
        rc, _, err = run_py(["-q", "--batch", self.jobs(lines)], self.tmp)
        self.assertEqual((rc, err), (0, ""))  # -q: no warnings, no summary lines, and no batch line either
        plain = self.read("a.fa")
        rc, _, err = run_py(["-nomismatch", "--batch", self.jobs(lines)], self.tmp)
        self.assertEqual(rc, 0)
        self.assertNotEqual(self.read("a.fa"), plain)
        single = self.tmp / "single"
        single.mkdir()
        run_py([str(WARN[0]), str(WARN[1]), "-nomismatch", "-o", "a.fa"], single)
        self.assertEqual(self.read("a.fa"), (single / "a.fa").read_bytes())


class TestOutDir(_Tmp):
    def test_names_by_format(self):
        lines = [f"{CLEAN[0]} {CLEAN[1]}", f"{CLUSTAL[0]} {CLUSTAL[1]} -output clustal",
                 f"{TABLE2[0]} {TABLE2[1]} -codontable 2 -output paml -o explicit.phy"]
        rc, _, err = run_py(["--batch", self.jobs(lines), "--batch-out-dir", "out/sub", "--report", "r.json"],
                            self.tmp)
        self.assertEqual(rc, 0, err)
        self.assertEqual(sorted(os.listdir(self.tmp / "out" / "sub")), ["pep.aln", "pep.fasta"])  # created
        self.assertTrue((self.tmp / "explicit.phy").is_file())
        rep = json.loads(self.read("r.json"))
        self.assertEqual([j["output"] for j in rep["jobs"]],
                         [os.path.join("out/sub", "pep.fasta"), os.path.join("out/sub", "pep.aln"), "explicit.phy"])
        self.assertEqual(rep["jobs"][0]["report"]["output"]["path"], os.path.join("out/sub", "pep.fasta"))
        self.assertEqual(batch.EXTENSION, {"fasta": "fasta", "clustal": "aln", "paml": "paml", "codon": "codon"})

    def test_duplicate_outputs_are_refused_before_running(self):
        lines = [f"{CLEAN[0]} {CLEAN[1]} -o x.fa", f"{TABLE2[0]} {TABLE2[1]} -o ./x.fa"]
        rc, _, err = run_py(["--batch", self.jobs(lines)], self.tmp)
        self.assertEqual(rc, 2)
        self.assertIn("coal: error: jobs.txt line 2: writes ./x.fa, the output of the job on line 1", err)
        self.assertFalse((self.tmp / "x.fa").exists())
        # two families whose protein file is both called pep.fasta collide in --batch-out-dir
        lines = [f"{CLEAN[0]} {CLEAN[1]}", f"{TABLE2[0]} {TABLE2[1]}"]
        rc, _, err = run_py(["--batch", self.jobs(lines), "--batch-out-dir", "d"], self.tmp)
        self.assertEqual(rc, 2)
        self.assertIn("jobs.txt line 2: writes d/pep.fasta, the output of the job on line 1".replace("/", os.sep),
                      err)
        self.assertFalse((self.tmp / "d").exists())


class TestFailurePolicy(_Tmp):
    def test_failing_middle_job_does_not_stop_the_batch(self):
        lines = [f"{CLEAN[0]} {CLEAN[1]} -o a.fa", f"{BROKEN[0]} {BROKEN[1]} -o b.fa", f"{WARN[0]} {WARN[1]} -o c.fa"]
        rc, out, err = run_py(["--batch", self.jobs(lines), "-q", "--json"], self.tmp)
        self.assertEqual(rc, 1)
        self.assertTrue((self.tmp / "a.fa").is_file() and (self.tmp / "c.fa").is_file())
        self.assertFalse((self.tmp / "b.fa").exists())
        self.assertIn(f"coal: [2/3 {BROKEN[0]}] error: ", err)
        rep = json.loads(out)
        self.assertEqual([j["exit_code"] for j in rep["jobs"]], [0, 1, 0])
        self.assertEqual(rep["summary"], {"jobs": 3, "ok": 2, "failed": 1, "skipped_sequences": 0})
        self.assertEqual(rep["jobs"][1]["report"]["exit_code"], 1)
        self.assertIn("error", rep["jobs"][1]["report"])

    def test_exit_2_when_a_job_left_sequences_out_and_none_failed(self):
        lines = [f"{CLEAN[0]} {CLEAN[1]} -o a.fa", f"{BROKEN[0]} {BROKEN[1]} -o b.fa"]
        rc, out, err = run_py(["--keep-going", "--batch", self.jobs(lines), "--json"], self.tmp)
        self.assertEqual(rc, 2, err)
        self.assertTrue((self.tmp / "b.fa").is_file())
        rep = json.loads(out)
        self.assertEqual([j["exit_code"] for j in rep["jobs"]], [0, 2])
        self.assertEqual(rep["summary"]["failed"], 0)
        self.assertEqual(rep["summary"]["ok"], 2)
        self.assertGreater(rep["summary"]["skipped_sequences"], 0)
        self.assertIn("coal: batch: 2 jobs, 2 ok, 1 with sequences left out, 0 failed", err)
        # a failed job wins over a partial one
        lines.append(f"missing.fa {BROKEN[1]} -o c.fa")
        rc, out, err = run_py(["--keep-going", "--batch", self.jobs(lines)], self.tmp)
        self.assertEqual(rc, 1, err)
        self.assertIn("coal: batch: 3 jobs, 2 ok, 1 with sequences left out, 1 failed", err)

    def test_exit_0_when_every_job_is_ok(self):
        lines = [f"{CLEAN[0]} {CLEAN[1]} -o a.fa", f"{CLEAN[0]} {CLEAN[1]} -o b.fa -output codon"]
        rc, out, err = run_py(["--batch", self.jobs(lines)], self.tmp)
        self.assertEqual((rc, out), (0, ""))
        self.assertEqual(err, "coal: batch: 2 jobs, 2 ok, 0 failed\n")


class TestJobsValidation(_Tmp):
    def check(self, lines, main=(), contains=""):
        rc, out, err = run_py([*main, "--batch", self.jobs(lines)], self.tmp)
        self.assertEqual((rc, out), (2, ""), err)
        self.assertIn(contains, err)
        self.assertTrue(err.startswith("usage: coal "), err[:80])
        self.assertFalse([f for f in os.listdir(self.tmp) if f != "jobs.txt"], "nothing may run")
        return err

    def test_line_numbers(self):
        good = f"{CLEAN[0]} {CLEAN[1]} -o ok.fa"
        self.check([good, "", "# c", f"{CLEAN[0]} {CLEAN[1]}"],
                   contains="coal: error: jobs.txt line 4: the job has no -o")
        self.check([good, f"{CLEAN[0]} {CLEAN[1]} -o x.fa --bogus"],
                   contains="coal: error: jobs.txt line 2: unrecognized arguments: --bogus")
        self.check([good, f"{CLEAN[0]} -o x.fa"],
                   contains="jobs.txt line 2: the following arguments are required: nuc_fasta")
        self.check([good, f"{CLEAN[0]} {CLEAN[1]} -o x.fa -codontable 7"],
                   contains="jobs.txt line 2: argument -codontable/--codontable: invalid choice: ")
        self.check([good, f'{CLEAN[0]} {CLEAN[1]} -o "x.fa'],
                   contains="jobs.txt line 2: No closing quotation (quote a path that has spaces; backslashes are ordinary characters)")
        for bad in ("--resume", "--json", "--batch other.txt", "--batch-out-dir d"):
            self.check([good, f"{CLEAN[0]} {CLEAN[1]} -o x.fa {bad}"],
                       contains=f"jobs.txt line 2: {bad.split()[0]}")
        self.check([good, f"{CLEAN[0]} {CLEAN[1]} -o x.fa -h"], contains="jobs.txt line 2: -h and --version")

    def test_main_command_line_mistakes(self):
        good = f"{CLEAN[0]} {CLEAN[1]} -o ok.fa"
        self.check([good], main=[str(CLEAN[0]), str(CLEAN[1])],
                   contains="coal: error: --batch takes the input files from the JOBS file")
        self.check([good], main=["--resume"], contains="coal: error: --resume is not supported with --batch yet")
        self.check([], contains="coal: error: no jobs in jobs.txt")
        rc, _, err = run_py(["--batch", "missing.txt"], self.tmp)
        self.assertEqual(rc, 2)
        self.assertIn("coal: error: cannot read the JOBS file missing.txt: No such file", err)
        rc, _, err = run_py([str(CLEAN[0]), str(CLEAN[1]), "--batch-out-dir", "d"], self.tmp)
        self.assertEqual(rc, 2)
        self.assertIn("coal: error: --batch-out-dir needs --batch", err)
        rc, _, err = run_py(["--bat", "x"], self.tmp)  # ambiguous, as any argparse abbreviation
        self.assertEqual(rc, 2)
        self.assertIn("ambiguous option: --bat could match --batch, --batch-out-dir", err)

    def test_shlex_quoting_works(self):
        d = self.tmp / "a dir"
        d.mkdir()
        shutil.copy(CLEAN[0], d / "p.fa")
        shutil.copy(CLEAN[1], d / "n.fa")
        rc, _, err = run_py(["--batch", self.jobs(["'a dir/p.fa' \"a dir/n.fa\" -o \"a dir/out.fa\""])], self.tmp)
        self.assertEqual(rc, 0, err)
        self.assertTrue((d / "out.fa").is_file())

    def test_backslashes_are_ordinary_characters(self):
        """Windows paths (D:\\a\\pep.fasta) must survive the split on every platform: no POSIX escapes."""
        from coal.batch import split_job_line
        self.assertEqual(split_job_line(r"D:\a\COAL\pep.fasta D:\a\nuc.fasta -o D:\out\x.fa"),
                         [r"D:\a\COAL\pep.fasta", r"D:\a\nuc.fasta", "-o", r"D:\out\x.fa"])
        self.assertEqual(split_job_line("\"C:\\my data\\p.aln\" 'n.fa' -codontable 2"),
                         [r"C:\my data\p.aln", "n.fa", "-codontable", "2"])
        d = self.tmp / "bs"
        d.mkdir()
        shutil.copy(CLEAN[0], d / "p.fa")
        shutil.copy(CLEAN[1], d / "n.fa")
        # a literal backslash in a (POSIX) file name is kept as given
        if os.name != "nt":
            shutil.copy(CLEAN[1], d / "n\\1.fa")
            rc, _, err = run_py(["--batch", self.jobs([f"{d}/p.fa {d}/n\\1.fa -o {d}/o.fa"])], self.tmp)
            self.assertEqual(rc, 0, err)
            self.assertTrue((d / "o.fa").is_file())


class TestReportAndStdin(_Tmp):
    def test_batch_report_shape(self):
        lines = [f"{CLEAN[0]} {CLEAN[1]} -o a.fa", f"{BROKEN[0]} {BROKEN[1]} -o b.fa"]
        rc, out, err = run_py(["--batch", self.jobs(lines), "--report", "batch.json", "--json"], self.tmp)
        self.assertEqual(rc, 1)
        rep = json.loads(self.read("batch.json"))
        self.assertEqual(list(rep), ["schema", "version", "command_line", "jobs", "summary"])
        self.assertEqual(rep["schema"], "coal.batch/1")
        from coal import __version__
        self.assertEqual(rep["version"], __version__)
        self.assertEqual(rep["command_line"], "coal --batch jobs.txt --report batch.json --json")
        for k, job in enumerate(rep["jobs"], 1):
            self.assertEqual(list(job), ["index", "line", "exit_code", "output", "report"])
            self.assertEqual(job["index"], k)
            self.assertEqual(job["line"], lines[k - 1])
            self.assertEqual(job["report"]["schema"], "coal.report/1")
            self.assertEqual(job["report"]["exit_code"], job["exit_code"])
            self.assertTrue(job["report"]["sequences"], "the per-job report has its sequence entries")
        self.assertEqual(rep["summary"], {"jobs": 2, "ok": 1, "failed": 1, "skipped_sequences": 0})
        on_stdout = json.loads(out)
        for j in (*rep["jobs"], *on_stdout["jobs"]):
            j["report"].pop("runtime_seconds")
        self.assertEqual(on_stdout, rep)

    def test_no_report_without_report_or_json(self):
        rc, out, err = run_py(["--batch", self.jobs([f"{CLEAN[0]} {CLEAN[1]} -o a.fa"])], self.tmp)
        self.assertEqual((rc, out), (0, ""))
        self.assertEqual(os.listdir(self.tmp).count("batch.json"), 0)

    def test_jobs_from_stdin(self):
        text = f"# from stdin\n{CLEAN[0]} {CLEAN[1]} -o a.fa\n{TABLE2[0]} {TABLE2[1]} -o b.fa -codontable 2\n"
        rc, out, err = run_py(["--batch", "-"], self.tmp, stdin=text)
        self.assertEqual((rc, out), (0, ""), err)
        self.assertTrue((self.tmp / "a.fa").is_file() and (self.tmp / "b.fa").is_file())
        rc, out, err = run_py(["--batch", "-"], self.tmp, stdin="")
        self.assertEqual(rc, 2)
        self.assertIn("coal: error: no jobs in -", err)


class TestStderrPrefixes(_Tmp):
    def test_every_job_line_is_prefixed(self):
        lines = [f"{WARN[0]} {WARN[1]} -o a.fa", f"{TABLE2[0]} {TABLE2[1]} -o b.fa"]  # a warning; a notice + error
        rc, _, err = run_py(["--batch", self.jobs(lines), "--progress"], self.tmp)
        self.assertEqual(rc, 1)
        for ln in err.splitlines():
            self.assertTrue(ln.startswith("coal: "), ln)
            if ln.startswith("coal: progress: ") or ln.startswith("coal: batch: "):
                continue
            self.assertRegex(ln, r"^coal: \[[12]/2 \S+\] (warning|error|notice|summary): ")
        self.assertIn(f"coal: progress: job 1/2 {WARN[0]}\n", err)
        self.assertIn(f"coal: [1/2 {WARN[0]}] warning: ", err)
        self.assertIn(f"coal: [1/2 {WARN[0]}] summary: ", err)
        self.assertIn(f"coal: [2/2 {TABLE2[0]}] notice: ", err)
        self.assertIn(f"coal: [2/2 {TABLE2[0]}] error: ", err)
        # the same lines as the single runs, but for the prefix
        single = self.tmp / "single"
        single.mkdir()
        one = run_py([str(WARN[0]), str(WARN[1]), "-o", "a.fa"], single)[2]
        two = run_py([str(TABLE2[0]), str(TABLE2[1]), "-o", "b.fa"], single)[2]
        want = (one.replace("coal: ", f"coal: [1/2 {WARN[0]}] ") + two.replace("coal: ", f"coal: [2/2 {TABLE2[0]}] "))
        got = "".join(ln + "\n" for ln in err.splitlines()
                      if not ln.startswith(("coal: progress: ", "coal: batch: ")))
        self.assertEqual(got, want)


class TestEngines(_Tmp):
    """--batch never leaves the Python engine; coal-cpp refuses it."""

    @unittest.skipIf(os.name == "nt", "uses a POSIX shell wrapper to see whether coal-cpp ran")
    @unittest.skipIf(cpp_binary() is None, "coal-cpp not built (set COAL_CPP_BIN)")
    def test_nothing_is_dispatched_to_cpp(self):
        mark = self.tmp / "cpp_ran"
        wrapper = self.tmp / "coal-cpp"
        wrapper.write_text(f'#!/bin/sh\ntouch "{mark}"\nexec "{cpp_binary()}" "$@"\n')
        wrapper.chmod(0o755)
        env = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "COAL_CPP": str(wrapper), "COAL_ENGINE": "cpp"}
        name = self.jobs([f"{CLEAN[0]} {CLEAN[1]} -o a.fa", f"{TABLE2[0]} {TABLE2[1]} -o b.fa -codontable 2"])
        r = subprocess.run([sys.executable, "-m", "coal", "--batch", name], cwd=self.tmp, env=env,
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertFalse(mark.exists(), "the batch (or a job of it) was handed to coal-cpp")
        self.assertEqual(r.stderr, "coal: batch: 2 jobs, 2 ok, 0 failed\n")
        # a single run with the same environment does dispatch (the wrapper is called), so the check is meaningful
        r = subprocess.run([sys.executable, "-m", "coal", str(CLEAN[0]), str(CLEAN[1]), "-o", "c.fa"], cwd=self.tmp,
                           env=env, capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(mark.exists())

    @unittest.skipIf(cpp_binary() is None, "coal-cpp not built (set COAL_CPP_BIN)")
    def test_cpp_refuses_batch(self):
        msg = ("coal: error: --batch is only available in the Python command (coal); coal-cpp converts one alignment "
               "per run")
        for args in (["--batch", "jobs.txt"], ["--batch=jobs.txt"], ["--batch", "jobs.txt", "-codontable", "2"],
                     [str(CLEAN[0]), str(CLEAN[1]), "--batch", "j"], [str(CLEAN[0]), str(CLEAN[1]), "--batch-out-dir", "d"],
                     [str(CLEAN[0]), str(CLEAN[1]), "--batch-o", "d"]):
            with self.subTest(args=args):
                r = subprocess.run([str(cpp_binary()), *args], cwd=self.tmp, capture_output=True, text=True,
                                   timeout=60)
                self.assertEqual(r.returncode, 2)
                self.assertTrue(r.stderr.startswith("usage: coal "), r.stderr)
                self.assertEqual(r.stderr.strip().splitlines()[-1], msg)
        h = subprocess.run([str(cpp_binary()), "--help"], capture_output=True, text=True, timeout=60).stdout
        self.assertIn("--batch JOBS", h)
        self.assertIn("(not in coal-cpp)", h)


if __name__ == "__main__":
    unittest.main()
