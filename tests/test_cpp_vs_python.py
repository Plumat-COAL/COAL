"""Random differential test: the C++ engine (coal-cpp) must behave exactly like
the Python engine on 1,000 random inputs.

Each input: 1-4 sequences, 5-120 codons, gaps, wrong codons, frameshift digits
(1/2/4/5), extra DNA before/after the gene, mixed case, ambiguous bases, B/Z/'.'/
'*' in the protein, the occasional invalid character, shuffled or renamed IDs,
any supported genetic code (with the matching or the default -codontable),
any output format, -nogap / -nomismatch.  Compared: stdout, full stderr (warnings
and error messages) and exit code.

Skipped if the binary is not built (COAL_CPP_BIN, default ./build/coal-cpp).
COAL_RANDOM_CHUNK=0..3 runs one quarter of the inputs; COAL_RANDOM_N changes the count.
"""
from __future__ import annotations

import io
import os
import random
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal.cli import main as coal_main  # noqa: E402
from coal.tables import SUPPORTED, get_code  # noqa: E402
from test_differential import cpp_binary  # noqa: E402

SEED = 20261017
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
    pep = "".join(f">{ids[i]}\n{rows[i]}\n" for i in range(nseq))
    nuc = "".join(f">{dna_ids[i]} desc\n{dnas[i]}\n" for i in order)
    if rng.random() < 0.01:
        nuc = nuc.split("\n", 2)[-1] if nseq > 1 else nuc + f">extra\nATG\n"
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
    return rc, out.getvalue(), strip_summary(err.getvalue())


def run_cpp(binary: Path, d: Path, opts):
    r = subprocess.run([str(binary), "pep.fasta", "nuc.fasta", *opts], cwd=d, capture_output=True, timeout=60)
    return r.returncode, r.stdout.decode("ascii"), strip_summary(r.stderr.decode("ascii", "replace"))


def strip_summary(err: str) -> str:
    """Drop the plain-language "coal: summary:" lines (not part of PAL2NAL-compatible stderr)."""
    return "".join(l for l in err.splitlines(True) if not l.startswith("coal: summary:"))


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
            py = run_python(d, opts)
            cc = run_cpp(binary, d, opts)
            stats["error" if py[0] else ("warning" if "warning" in py[2] else "ok")] += 1
            if py != cc:
                what = [n for n, a, b in zip(("rc", "stdout", "stderr"), py, cc) if a != b]
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
                self.assertEqual(ref[2].count("does not correspond"), 3, ref[2][:300])
                for j in ("2", "3"):
                    with self.subTest(opts=opts, j=j):
                        self.assertEqual(run_cpp(binary, d, [*opts, "-j", j]), ref)
                if opts[1] == "fasta" and len(opts) == 2:
                    self.assertEqual(run_python(d, [*opts, "-j", "1"]), ref)
