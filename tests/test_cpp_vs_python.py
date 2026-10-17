"""Random differential test: the C++ engine (coal-cpp) must behave exactly like
the Python engine on 1,000 random inputs.

Each input: 1-4 sequences, 5-120 codons, gaps, wrong codons, frameshift digits
(1/2/4/5), extra DNA before/after the gene, mixed case, ambiguous bases, B/Z/'.'/
'*' in the protein, the occasional invalid character, shuffled or renamed IDs,
IDs that only match after normalising (version suffix, pipe header, case),
any supported genetic code (with the matching or the default -codontable),
any output format, -nogap / -nomismatch, --keep-going, --max-warnings, and --json
run reports.  Compared: stdout, full stderr (warnings, error messages and the
"coal: summary:" lines), exit code, and for --json runs the report (all fields
except runtime_seconds) and the output file.

Skipped if the binary is not built (COAL_CPP_BIN, default ./build/coal-cpp).
COAL_RANDOM_CHUNK=0..3 runs one quarter of the inputs; COAL_RANDOM_N changes the count.

COAL_REF_SRC=<dir> (CPP-9): compare against the Python engine in that src directory instead of this branch's
(e.g. a worktree of the reviewed engine: COAL_REF_SRC=../coal-ref/src). The same random inputs are generated;
options the reference does not have (--keep-going, --max-warnings, --json) are left out, coal-cpp's own
"coal: summary:" and "IDs don't match ... paired by file order" lines are not compared, and the cases whose DNA IDs
only match after normalising use plain IDs (ID normalisation is coal-cpp's own; see LOG "Needs reconciling").
"""
from __future__ import annotations

import io
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
import zlib
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
REF_SRC = os.environ.get("COAL_REF_SRC")  # CPP-9
if REF_SRC:
    sys.path.insert(0, str(Path(REF_SRC).resolve()))
sys.path.insert(0, str(ROOT / "tests"))
from coal.cli import main as coal_main  # noqa: E402
from coal.tables import SUPPORTED, get_code  # noqa: E402
from test_differential import cpp_binary  # noqa: E402

SEED = int(os.environ.get("COAL_RANDOM_SEED", "20261017"))  # CPP-22: second seeds
N = int(os.environ.get("COAL_RANDOM_N", "1000"))
FORMATS = ["fasta", "clustal", "paml", "codon"]
BASES = "ACGT"


def make_case(rng: random.Random):
    """Return (pep_text, nuc_text, options)."""
    gene_table = rng.choice(SUPPORTED)
    code = get_code(gene_table)
    aas = sorted(a for a in code.codons if a != "*")
    nseq = rng.randint(1, 4)
    ncol = rng.randint(5, 120)
    cons = [rng.choice(aas) for _ in range(ncol)]
    cons[0] = "M" if rng.random() < 0.8 else cons[0]
    rows, dnas = [], []
    for s in range(nseq):
        row, dna = [], []
        for c in range(ncol):
            r = rng.random()
            if c > 0 and r < 0.08:
                row.append("-")
                continue
            aa = cons[c] if rng.random() < 0.85 else rng.choice(aas)
            if c == ncol - 1 and rng.random() < 0.05:
                aa = "*"
            if aa == "M" and c == 0 and rng.random() < 0.3 and code.starts:
                codon = rng.choice(sorted(code.starts))
            else:
                codon = rng.choice(sorted(code.codons[aa]))
            if rng.random() < 0.004:                     # wrong codon
                codon = "".join(rng.choice(BASES) for _ in range(3))
            if rng.random() < 0.002:                     # ambiguous base
                k = rng.randrange(3)
                codon = codon[:k] + rng.choice("NRYKMSW") + codon[k + 1:]
            if c > 2 and rng.random() < 0.02:            # frameshift digit
                d = rng.choice([1, 2, 4, 5])
                row.append(str(d))
                dna.append((codon + "".join(rng.choice(BASES) for _ in range(3)))[:d])
                continue
            if rng.random() < 0.005:
                aa = rng.choice("BZ")
            row.append(aa)
            dna.append(codon)
        if all(ch == "-" for ch in row):
            row[0] = cons[0]
            dna.append(rng.choice(sorted(code.codons.get(cons[0], {"ATG"}))))
        seq = "".join(dna)
        if rng.random() < 0.3:
            seq = "".join(rng.choice(BASES) for _ in range(rng.randint(1, 12))) + seq
        if rng.random() < 0.3:
            seq = seq + rng.choice(["TAA", "TGA", "TAG", ""]) + "".join(rng.choice(BASES) for _ in range(rng.randint(0, 9)))
        if rng.random() < 0.15:
            a = rng.randrange(len(seq) + 1)
            seq = seq[:a] + seq[a:].lower()
        if rng.random() < 0.03:
            seq = seq.replace("T", "U")
        prow = "".join(row)
        if rng.random() < 0.1:
            prow = prow.lower()
        if rng.random() < 0.05:
            prow = prow.replace("-", ".")
        if rng.random() < 0.01:
            k = rng.randrange(len(prow))
            prow = prow[:k] + rng.choice("@#!") + prow[k + 1:]
        rows.append(prow)
        dnas.append(seq)
    ids = [f"seq{i}" for i in range(nseq)]
    dna_ids = list(ids)
    order = list(range(nseq))
    r = rng.random()
    if r < 0.15:
        rng.shuffle(order)
    elif r < 0.2:
        dna_ids = [f"d{i}" for i in range(nseq)]
    elif r < 0.25:                                  # version suffixes, shuffled
        dna_ids = [f"{x}.{rng.randint(1, 3)}" for x in ids]
        rng.shuffle(order)
    elif r < 0.3:                                   # pipe headers
        dna_ids = [f"sp|{x}|{x.upper()}_HUMAN" for x in ids]
    elif r < 0.33:                                  # case, one ambiguous pair sometimes
        dna_ids = [x.upper() for x in ids]
        if nseq > 1 and rng.random() < 0.3:
            dna_ids[1] = dna_ids[0].lower() + ".2"
    if REF_SRC and 0.2 <= r < 0.33:  # CPP-15: ID normalisation is coal-cpp's own (safe ID pairing): plain IDs
        dna_ids = list(ids)
    pep = "".join(f">{ids[i]}\n{rows[i]}\n" for i in range(nseq))
    nuc = "".join(f">{dna_ids[i]} desc\n{dnas[i]}\n" for i in order)
    if rng.random() < 0.01:
        nuc = nuc.split("\n", 2)[-1] if nseq > 1 else nuc + ">extra\nATG\n"
    table = gene_table if rng.random() < 0.85 else 1
    opts = ["-output", rng.choice(FORMATS)]
    if table != 1 or rng.random() < 0.2:
        opts += ["-codontable", str(table)]
    if opts[1] != "codon" or rng.random() < 0.05:
        r = rng.random()
        if r < 0.2:
            opts.append("-nogap")
        elif r < 0.4:
            opts.append("-nomismatch")
    if rng.random() < 0.15:
        opts.append("--keep-going")
    if rng.random() < 0.1:
        opts += ["--max-warnings", str(rng.randint(0, 3))]
    if rng.random() < 0.25:
        opts += ["-o", "out.txt", "--json"]       # compare the run reports too
    if REF_SRC:  # CPP-15: --strict on a fifth and -q on a seventh of the cases (chosen without drawing from rng)
        h = zlib.crc32(pep.encode() + nuc.encode())
        opts += ["--strict"] * (h % 5 == 0) + ["-q"] * (h % 7 == 0)
    if REF_SRC:  # CPP-9: the reference has none of the C++-only options
        # (-o out.txt only comes with --json here: dropped too, so stdout is compared)
        opts = [o for k, o in enumerate(opts) if o not in ("--keep-going", "--json", "--max-warnings", "-o")
                and not (k > 0 and opts[k - 1] in ("--max-warnings", "-o"))]
    return pep, nuc, opts


def run_python(d: Path, opts):
    out, err = io.StringIO(), io.StringIO()
    cwd = os.getcwd()
    os.chdir(d)
    try:
        with redirect_stdout(out), redirect_stderr(err):
            try:
                rc = coal_main(["pep.fasta", "nuc.fasta", *opts])
            except SystemExit as e:
                rc = e.code
    finally:
        os.chdir(cwd)
    return rc, out.getvalue(), err.getvalue()


def run_cpp(binary: Path, d: Path, opts):
    r = subprocess.run([str(binary), "pep.fasta", "nuc.fasta", *opts], cwd=d, capture_output=True, timeout=60)
    # stderr is text (CRLF on Windows, as Python's own stderr): compare lines (CPP-8)
    err = r.stderr.decode("ascii", "replace").replace("\r\n", "\n")
    if REF_SRC:  # CPP-9: summary lines are coal-cpp's own
        err = "".join(ln for ln in err.splitlines(True) if not ln.startswith(
            ("coal: summary:", "coal: warning: IDs don't match, so sequences were paired by file order")))
    return r.returncode, r.stdout.decode("ascii", "replace"), err


def normalise(res, d: Path, opts):
    """For --json runs: (rc, report without runtime, stderr, output file); else unchanged."""
    if "--json" not in opts:
        return res
    rc, out, err = res
    rep = json.loads(out) if out.strip() else None
    if rep:
        rep.pop("runtime_seconds", None)
    f = d / "out.txt"
    data = f.read_bytes() if f.exists() else None
    if f.exists():
        f.unlink()
    return rc, rep, err, data


def selected():
    rng = random.Random(SEED)
    cases = [make_case(rng) for _ in range(N)]
    chunk = os.environ.get("COAL_RANDOM_CHUNK")
    idx = range(N)
    if chunk is not None:
        k = int(chunk)
        idx = range(k * N // 4, (k + 1) * N // 4)
    return [(i, cases[i]) for i in idx]


def compare(binary: Path):
    diffs, stats = [], {"ok": 0, "error": 0, "warning": 0}
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        for i, (pep, nuc, opts) in selected():
            (d / "pep.fasta").write_text(pep, newline="\n")
            (d / "nuc.fasta").write_text(nuc, newline="\n")
            py = normalise(run_python(d, opts), d, opts)
            cc = normalise(run_cpp(binary, d, opts), d, opts)
            stats["error" if py[0] == 1 else ("warning" if "warning" in py[2] else "ok")] += 1
            if py != cc:
                what = [n for n, a, b in zip(("rc", "stdout/report", "stderr", "output"), py, cc) if a != b]
                diffs.append((i, opts, what, py, cc))
    return diffs, stats


class TestCppMatchesPython(unittest.TestCase):
    def test_random_inputs(self):
        binary = cpp_binary()
        if binary is None:
            self.skipTest("coal-cpp not built")
        diffs, stats = compare(binary)
        msg = "\n".join(f"#{i} {opts}: {what}\n  python: {py!r:.400}\n  cpp:    {cc!r:.400}"
                        for i, opts, what, py, cc in diffs[:5])
        self.assertEqual(diffs, [], f"{len(diffs)} differences (stats {stats}):\n{msg}")


if __name__ == "__main__":
    b = cpp_binary()
    if b is None:
        sys.exit("coal-cpp not built")
    diffs, stats = compare(b)
    for i, opts, what, py, cc in diffs[:10]:
        print(f"#{i} {opts}: {what}\n  python: {py!r:.600}\n  cpp:    {cc!r:.600}")
    print(f"{len(selected()) - len(diffs)}/{len(selected())} identical; python outcomes {stats}")


class TestCppThreads(unittest.TestCase):
    """-j must not change anything: several 2M-nt batches, every format."""

    def test_jobs_do_not_change_output(self):
        binary = cpp_binary()
        if binary is None:
            self.skipTest("coal-cpp not built")
        sys.path.insert(0, str(ROOT / "benchmarks"))
        from gen_data import generate
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            generate(3000, 400, d, seed=11)        # 3.4M nt -> 2 batches
            # wrong codons early and late in the file: warnings must come out in input order
            lines = (d / "nuc.fasta").read_text().split("\n")
            for k in (1, 2999, 1500):
                seq = lines[2 * k + 1]
                lines[2 * k + 1] = seq[:30] + ("TGG" if seq[30:33] not in ("TGG",) else "AAA") + seq[33:]
            (d / "nuc.fasta").write_text("\n".join(lines), newline="\n")
            for opts in (["-output", "fasta"], ["-output", "clustal"], ["-output", "paml", "-nogap"],
                         ["-output", "codon"], ["-output", "fasta", "-nomismatch"]):
                ref = run_cpp(binary, d, [*opts, "-j", "1"])
                self.assertEqual(ref[0], 0, ref[2][:300])
                # CPP-11 (reference ENG-i2): no mismatch warnings with -nomismatch, as PAL2NAL
                want = 0 if "-nomismatch" in opts else 3
                self.assertEqual(ref[2].count("does not correspond"), want, ref[2][:300])
                for j in ("2", "3"):
                    with self.subTest(opts=opts, j=j):
                        self.assertEqual(run_cpp(binary, d, [*opts, "-j", j]), ref)
                if opts[1] == "fasta" and len(opts) == 2:
                    self.assertEqual(run_python(d, [*opts, "-j", "1"]), ref)


class TestCppResumeAndProgress(unittest.TestCase):
    """--resume and --progress in coal-cpp: interrupted + resumed C++ runs must give
    exactly what an uninterrupted Python run gives."""

    @classmethod
    def setUpClass(cls):
        cls.binary = cpp_binary()
        if cls.binary is None:
            raise unittest.SkipTest("coal-cpp not built")
        cls._tmp = tempfile.TemporaryDirectory()
        cls.d = Path(cls._tmp.name)
        sys.path.insert(0, str(ROOT / "benchmarks"))
        from gen_data import generate
        generate(20000, 15, cls.d, seed=8)
        lines = (cls.d / "nuc.fasta").read_text().split("\n")
        lines[2 * 12000 + 1] = lines[2 * 12000 + 1][:20] + lines[2 * 12000 + 1][21:]   # one broken gene
        (cls.d / "nuc_bad.fasta").write_text("\n".join(lines), newline="\n")

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_resume_matches_python(self):
        from test_engine_extras import interrupt
        for opts, nuc, rc_ok in ((["-output", "fasta", "-j", "2"], "nuc.fasta", 0),
                                 (["-output", "clustal"], "nuc.fasta", 0),
                                 (["-output", "fasta", "-nogap"], "nuc.fasta", 0),
                                 (["-output", "paml", "--keep-going"], "nuc_bad.fasta", 2),
                                 (["-output", "codon", "--keep-going", "-j", "2"], "nuc_bad.fasta", 2)):
            with self.subTest(opts=opts):
                out_py, out_cc = self.d / "py.out", self.d / "cc.out"
                for f in (out_py, out_cc):
                    if f.exists():
                        f.unlink()
                py = subprocess.run([sys.executable, "-m", "coal", "pep.fasta", nuc, "-o", str(out_py), *opts,
                                     "--report", "py.json"], cwd=self.d, capture_output=True, timeout=120,
                                    env={**os.environ, "PYTHONPATH": str(ROOT / "src")})
                self.assertEqual(py.returncode, rc_ok, py.stderr[-300:])
                cmd = [str(self.binary), "pep.fasta", nuc, "-o", str(out_cc), *opts]
                self.assertTrue(interrupt(cmd, self.d, after_checkpoint=3), "C++ run did not reach its checkpoint")
                ck = Path(str(out_cc) + ".coal-checkpoint.json")
                self.assertTrue(ck.exists())
                self.assertFalse(out_cc.exists())
                r = subprocess.run(cmd + ["--resume", "--report", "cc.json"], cwd=self.d, capture_output=True,
                                   timeout=120)
                self.assertEqual(r.returncode, rc_ok, r.stderr[-300:])
                self.assertIn(b"coal: note: resuming", r.stderr)
                self.assertEqual(out_cc.read_bytes(), out_py.read_bytes())
                self.assertFalse(ck.exists())
                a = json.loads((self.d / "py.json").read_text())
                b = json.loads((self.d / "cc.json").read_text())
                self.assertGreater(b.pop("resumed_from"), 0)
                for x in (a, b):
                    for k in ("runtime_seconds", "command", "output"):
                        x.pop(k)
                self.assertEqual(a, b)

    def test_refusal_and_no_checkpoint_messages(self):
        from test_engine_extras import interrupt
        d = self.d / "chg"
        d.mkdir(exist_ok=True)
        for f in ("pep.fasta", "nuc.fasta"):
            (d / f).write_bytes((self.d / f).read_bytes())
        cmd = [str(self.binary), "pep.fasta", "nuc.fasta", "-o", "o.fa"]
        self.assertTrue(interrupt(cmd, d))
        with open(d / "pep.fasta", "a") as fh:
            fh.write("\n")
        r = subprocess.run(cmd + ["--resume"], cwd=d, capture_output=True, timeout=60)
        self.assertEqual(r.returncode, 1)
        self.assertIn(b"The input files or options changed since the interrupted run, so it can't be resumed "
                      b"safely. Start a new run without --resume.", r.stderr)
        r = subprocess.run([str(self.binary), "pep.fasta", "nuc.fasta", "-o", "fresh.fa", "--resume"], cwd=d,
                           capture_output=True, timeout=60)
        self.assertEqual(r.returncode, 0)
        self.assertIn(b"coal: note: no checkpoint to resume for fresh.fa; starting a normal run", r.stderr)
        r = subprocess.run([str(self.binary), "pep.fasta", "nuc.fasta", "--resume"], cwd=d, capture_output=True,
                           timeout=60)
        self.assertEqual(r.returncode, 2)

    def test_progress(self):
        env = {**os.environ, "COAL_PROGRESS_MIN_SEQS": "100", "COAL_PROGRESS_EVERY": "0", "COAL_BATCH_NT": "20000"}
        r = subprocess.run([str(self.binary), "pep.fasta", "nuc.fasta", "--progress"], cwd=self.d,
                           capture_output=True, timeout=60, env=env)
        lines = [l for l in r.stderr.decode().splitlines() if l.startswith("coal: progress:")]
        self.assertGreater(len(lines), 2)
        self.assertTrue(lines[-1].startswith("coal: progress: 100% (20,000 of 20,000 sequences), done in"))
        plain = subprocess.run([str(self.binary), "pep.fasta", "nuc.fasta"], cwd=self.d, capture_output=True,
                               timeout=60)
        self.assertEqual(r.stdout, plain.stdout)
        self.assertNotIn(b"progress", plain.stderr)
        # small input (100 sequences, under the 10,000 threshold): nothing, even reporting "every 0 s"
        lines = (self.d / "pep.fasta").read_text().split("\n")[:200]
        (self.d / "small_pep.fasta").write_text("\n".join(lines) + "\n")
        lines = (self.d / "nuc.fasta").read_text().split("\n")[:200]
        (self.d / "small_nuc.fasta").write_text("\n".join(lines) + "\n")
        small = subprocess.run([str(self.binary), "small_pep.fasta", "small_nuc.fasta", "--progress"], cwd=self.d,
                               capture_output=True, timeout=60, env={**os.environ, "COAL_PROGRESS_EVERY": "0"})
        self.assertEqual(small.returncode, 0)
        self.assertNotIn(b"progress", small.stderr)
