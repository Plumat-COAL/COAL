"""Ka/Ks engine (`coal-cpp kaks`, `coal-cpp dnds`, cpp/src/kaks.cpp; docs/kaks.md).

* Against the references, stored in tests/kaks/golden (made by tests/kaks/make_fixtures.py): KaKs_Calculator 2.0
  (built with 15 significant digits) for NG, LWL, MLWL, LPB, MLPB, YN, MYN on the edge cases and on 48 simulated
  pairs; PAML yn00 4.10.10 (12 digits) for YN00 with weighting 1 and 0. Values must agree to 1e-9; the only
  allowed differences are the documented ones (docs/kaks.md "Differences from KaKs_Calculator"), each checked to
  be exactly that difference.
* Regression: COAL's own outputs (golden/*.coal.*) field by field (numbers to a relative 1e-6).
* Policy and command line: statuses, the header, method names, FASTA pairs, threads, a 20,000-codon gene, `coal`
  (Python) handing the subcommands to coal-cpp.
Skipped if build/coal-cpp (COAL_CPP_BIN) is not built.
"""
from __future__ import annotations

import math
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
KDIR = ROOT / "tests" / "kaks"
GOLD = KDIR / "golden"
METHODS = ["NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN"]
HEADER = ("Sequence\tMethod\tKa\tKs\tKa/Ks\tP-Value(Fisher)\tLength\tS-Sites\tN-Sites\tFold-Sites(0:2:4)\t"
          "Substitutions\tS-Substitutions\tN-Substitutions\tFold-S-Substitutions(0:2:4)\tFold-N-Substitutions(0:2:4)\t"
          "Divergence-Time\tSubstitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)\tGC(1:2:3)\tML-Score\tAICc\t"
          "Akaike-Weight\tModel")
# KaKs_Calculator's NG/LWL family compares the letters of codons, so U in one sequence and T in the other count
# as a difference there (its YN/MYN read U as T); COAL reads U as T everywhere.
RNA_PAIRS = {"lowercase_rna"}


def cpp_bin() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


def num(x: str):
    try:
        v = float(x)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def lfact(n: float) -> float:
    return math.lgamma(n + 1) if n > 0 else 0.0


def p_observed(sd, nd, S, N) -> float:
    m = [sd, S - sd, nd, N - nd]
    R0, R1, C0, C1 = m[0] + m[2], m[1] + m[3], m[0] + m[1], m[2] + m[3]
    return math.exp(lfact(R0) + lfact(R1) + lfact(C0) + lfact(C1) - lfact(R0 + R1) - sum(lfact(x) for x in m))


def read_rows(text: str) -> dict:
    out = {}
    for line in text.splitlines()[1:]:
        f = line.split("\t")
        out[(f[0], f[1])] = f
    return out


class KaKsBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bin = cpp_bin()
        if cls.bin is None:
            raise unittest.SkipTest("coal-cpp not built (cmake --build build)")

    def run_cpp(self, *args, check=True, **kw):
        r = subprocess.run([str(self.bin), *map(str, args)], capture_output=True, text=True, **kw)
        if check:
            self.assertEqual(r.returncode, 0, r.stderr)
        return r

    def kaks15(self, axt: Path) -> dict:
        r = self.run_cpp("kaks", "-i", axt, "-q", "--digits", "15", *sum((["-m", m] for m in METHODS), []))
        return read_rows(r.stdout)

    def dnds(self, axt: Path, *extra) -> dict:
        r = self.run_cpp("dnds", axt, "-q", "--digits", "15", *extra)
        out = {}
        for line in r.stdout.splitlines()[1:]:
            f = line.split("\t")
            name = f[0] if f[1] == "NA" else f[0] + "&" + f[1]
            out[(name, f[2])] = f
        return out


class TestAgainstKaKsCalculator(KaKsBase):
    """Every field of every row equals KaKs_Calculator's, or differs in one of the documented ways."""

    def check_set(self, name: str):
        axt = KDIR / f"{name}.axt"
        ours = self.kaks15(axt)
        status = {k: v[13] for k, v in self.dnds(axt, *sum((["-m", m] for m in METHODS), [])).items()}
        counts = {}
        for m in METHODS:
            ref = read_rows((GOLD / f"{name}.KaKs_Calculator.{m}.kaks").read_text())
            for key, rf in ref.items():
                pair = key[0]
                if pair in RNA_PAIRS and m not in ("YN", "MYN"):
                    continue
                self.assertIn(key, ours, f"{key} missing in coal kaks output")
                of = ours[key]
                st = status[key]
                for i, (x, y) in enumerate(zip(rf, of)):
                    if x == y or i in (0, 1):
                        continue
                    why = None
                    if i == 17 and pair in RNA_PAIRS:
                        why = "lowercase GC"  # KaKs_Calculator counts upper-case G and C only
                    elif x == "-0" and y == "0":
                        why = "negative zero"
                    elif x in ("nan", "-nan") and y == "NA":
                        why = "nan"
                    elif ":" in x or "(" in x:  # composite fields: compare part by part
                        xs = x.replace("(", ":").replace(")", "").split(":")
                        ys = y.replace("(", ":").replace(")", "").split(":")
                        if len(xs) == len(ys) and all(
                                (num(a) is not None and num(b) is not None and close(num(a), num(b), 1e-9))
                                or a == b for a, b in zip(xs, ys)):
                            why = "rounding"
                        elif y == "NA" and st == "no_data":
                            why = "no data"
                    elif num(x) is not None and num(y) is not None and close(num(x), num(y), 1e-9):
                        why = "rounding"
                    elif i == 5 and num(x) is not None and num(y) is not None:
                        sd, nd, S, N = (float(of[j]) for j in (11, 12, 7, 8))
                        if abs(num(y) - num(x) - p_observed(sd, nd, S, N)) <= 1e-6 * num(y):
                            why = "p-obs"
                    elif i == 5 and y == "NA" and num(x) is not None:
                        if float(of[7]) < float(of[11]) or float(of[8]) < float(of[12]):
                            why = "p-neg"
                    if why is None and y == "NA" and st == "saturated":
                        why = st
                    if why is None and st == "no_data":  # nothing to estimate: KaKs_Calculator prints nan/defaults
                        why = st
                    self.assertIsNotNone(why, f"{name} {key} field {i} ({HEADER.split(chr(9))[i]}): "
                                              f"KaKs_Calculator {x!r}, COAL {y!r} (status {st})")
                    counts[why] = counts.get(why, 0) + 1
        return counts

    def test_edge_cases(self):
        c = self.check_set("edge_cases")
        self.assertGreater(c.get("saturated", 0), 0)

    def test_simulated(self):
        self.check_set("sim_small")


class TestAgainstPAML(KaKsBase):
    def test_yn00(self):
        for w in (1, 0):
            ours = self.dnds(KDIR / "sim_small.axt", "-m", "YN", "--yn-flavour", "paml", "--yn-weighting", str(w))
            lines = (GOLD / f"sim_small.yn00_w{w}.tsv").read_text().splitlines()[1:]
            self.assertEqual(len(lines), 48)
            for line in lines:
                pair, dN, dS, S, N, SEdN, SEdS = line.split("\t")
                o = ours[(pair, "YN")]
                for ref, col, tol in ((dN, 3, 1e-9), (dS, 4, 1e-9), (SEdN, 6, 1e-9), (SEdS, 7, 1e-9),
                                      (N, 8, 1e-8), (S, 9, 1e-8)):
                    self.assertTrue(close(float(ref), float(o[col]), tol), f"w={w} {pair} col {col}: {ref} vs {o[col]}")


class TestRegression(KaKsBase):
    def compare(self, got: str, want: str, what: str):
        g, w = got.splitlines(), want.splitlines()
        self.assertEqual(len(g), len(w), what)
        for a, b in zip(g, w):
            fa, fb = a.split("\t"), b.split("\t")
            self.assertEqual(len(fa), len(fb), what)
            for x, y in zip(fa, fb):
                if x == y:
                    continue
                xs = x.replace("(", ":").replace(")", "").split(":")
                ys = y.replace("(", ":").replace(")", "").split(":")
                ok = len(xs) == len(ys) and all(
                    p == q or (num(p) is not None and num(q) is not None and close(num(p), num(q), 1e-6))
                    for p, q in zip(xs, ys))
                self.assertTrue(ok, f"{what}: {x!r} != {y!r} in\n{a}\n{b}")

    def test_goldens(self):
        flags = sum((["-m", m] for m in METHODS), [])
        for name in ("edge_cases", "sim_small"):
            axt = KDIR / f"{name}.axt"
            self.compare(self.run_cpp("kaks", "-i", axt, "-q", *flags).stdout,
                         (GOLD / f"{name}.coal.kaks").read_text(), f"{name} kaks")
            self.compare(self.run_cpp("dnds", axt, "-q", *flags).stdout,
                         (GOLD / f"{name}.coal.dnds.tsv").read_text(), f"{name} dnds")


class TestPolicy(KaKsBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.rows = {}
        r = subprocess.run([str(cls.bin), "dnds", str(KDIR / "edge_cases.axt"), "-q",
                            *sum((["-m", m] for m in METHODS), [])], capture_output=True, text=True, check=True)
        for line in r.stdout.splitlines()[1:]:
            f = line.split("\t")
            cls.rows[(f[0], f[2])] = f

    def status(self, pair, m="YN"):
        return self.rows[(pair, m)][13]

    def test_statuses(self):
        for m in METHODS:
            self.assertEqual(self.status("identical", m), "ok")
            self.assertEqual(self.rows[("identical", m)][3:6], ["0", "0", "NA"])
            self.assertEqual(self.status("internal_stop", m), "warning")
            self.assertIn("internal stop", self.rows[("internal_stop", m)][14])
            self.assertEqual(self.status("terminal_stop", m), "ok")
            self.assertEqual(self.rows[("terminal_stop", m)][12], "120")
            self.assertEqual(self.rows[("terminal_stop_one_side", m)][12], "120")
            self.assertEqual(self.status("length_not_multiple_of_3", m), "error")
            self.assertEqual(self.status("all_gap", m), "no_data")
            self.assertEqual(self.status("met_trp_only", m), "no_data")
            self.assertEqual(self.status("saturated_unrelated", m), "saturated")
            self.assertEqual(self.rows[("gaps", m)][12], "112")
            self.assertEqual(self.rows[("ambiguous", m)][12], "115")
            self.assertEqual(self.rows[("all_synonymous", m)][13], "ok")
        for m in ("NG", "MLWL", "MLPB", "YN", "MYN"):
            self.assertEqual(self.rows[("all_synonymous", m)][3], "0")

    def test_rna_lowercase_same_as_upper_dna(self):
        text = (KDIR / "edge_cases.axt").read_text().split("\n\n")
        rec = [r for r in text if r.startswith("lowercase_rna")][0].splitlines()
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "u.axt"
            p.write_text(f"upper\n{rec[1].upper().replace('U', 'T')}\n{rec[2].upper()}\n\n")
            r = self.run_cpp("dnds", p, "-q", *sum((["-m", m] for m in METHODS), []))
            for line in r.stdout.splitlines()[1:]:
                f = line.split("\t")
                self.assertEqual(f[3:], self.rows[("lowercase_rna", f[2])][3:])


class TestCli(KaKsBase):
    def test_header_and_default_method(self):
        r = self.run_cpp("kaks", "-i", KDIR / "edge_cases.axt", "-q")
        lines = r.stdout.splitlines()
        self.assertEqual(lines[0], HEADER)
        self.assertTrue(all(l.split("\t")[1] == "YN" for l in lines[1:]))
        # KaKs_Calculator leaves out pairs it cannot read; so does coal kaks (and says so on stderr)
        self.assertFalse(any(l.startswith("length_not_multiple_of_3") for l in lines))
        r = self.run_cpp("kaks", "-i", KDIR / "edge_cases.axt")
        self.assertIn("left out", r.stderr)

    def test_kaks_calculator_flags(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "o.kaks"
            self.run_cpp("kaks", "-I", KDIR / "edge_cases.axt", "-O", out, "-M", "ng", "-m", "YN", "-C", "1", "-q")
            rows = out.read_text().splitlines()[1:]
            methods = [l.split("\t")[1] for l in rows]
            self.assertEqual(methods[:2], ["NG", "YN"])  # KaKs_Calculator's order, one row per method and pair

    def test_unimplemented_methods(self):
        for m in ("GY", "MA", "MS", "GYN", "GMYN", "NONE", "ALL"):
            r = self.run_cpp("kaks", "-i", KDIR / "edge_cases.axt", "-m", m, check=False)
            self.assertEqual(r.returncode, 2)
            self.assertIn("does not implement natively", r.stderr)
        r = self.run_cpp("dnds", KDIR / "edge_cases.axt", "-m", "XYZ", check=False)
        self.assertEqual(r.returncode, 2)
        self.assertIn("unknown method", r.stderr)
        r = self.run_cpp("kaks", "-i", KDIR / "edge_cases.axt", "-c", "7", check=False)
        self.assertEqual(r.returncode, 2)
        self.assertIn("not supported", r.stderr)

    def test_missing_input(self):
        r = self.run_cpp("dnds", "/nonexistent/x.axt", check=False)
        self.assertEqual(r.returncode, 1)

    def test_fasta_pairs(self):
        recs = (KDIR / "sim_small.axt").read_text().split("\n\n")[:3]
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "a.fasta"
            with open(fa, "w") as f:
                for i, r in enumerate(recs):
                    l = r.strip().splitlines()
                    f.write(f">s{i}a\n{l[1]}\n>s{i}b\n{l[2]}\n")
            allp = self.run_cpp("dnds", fa, "-q").stdout.splitlines()[1:]
            self.assertEqual(len(allp), 15)  # 6 sequences: every pair
            cons = self.run_cpp("dnds", fa, "-q", "--pairs", "consecutive").stdout.splitlines()[1:]
            self.assertEqual([l.split("\t")[:2] for l in cons], [["s0a", "s0b"], ["s1a", "s1b"], ["s2a", "s2b"]])
            # consecutive FASTA pairs = the AXT records
            axt = Path(td) / "a.axt"
            axt.write_text("\n\n".join(recs) + "\n\n")
            a = self.run_cpp("dnds", axt, "-q").stdout.splitlines()[1:]
            self.assertEqual([l.split("\t")[3:] for l in a], [l.split("\t")[3:] for l in cons])
            k = self.run_cpp("kaks", "-i", fa, "-q", "--pairs", "consecutive").stdout.splitlines()[1]
            self.assertTrue(k.startswith("s0a&s0b\tYN\t"))

    def test_threads_same_output(self):
        flags = sum((["-m", m] for m in METHODS), [])
        one = self.run_cpp("dnds", KDIR / "sim_small.axt", "-q", "-j", "1", *flags).stdout
        many = self.run_cpp("dnds", KDIR / "sim_small.axt", "-q", "-j", "4", *flags).stdout
        self.assertEqual(one, many)

    def test_long_gene(self):
        rng = random.Random(5)
        bases = "ACGT"
        stops = {"TAA", "TAG", "TGA"}
        cod = []
        while len(cod) < 20000:
            c = "".join(rng.choice(bases) for _ in range(3))
            if c not in stops:
                cod.append(c)
        s1 = "".join(cod)
        s2 = list(s1)
        for i in range(0, len(s2), 37):
            s2[i] = rng.choice([b for b in bases if b != s2[i]])
        s2 = "".join(s2)
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "long.axt"
            p.write_text(f"long\n{s1}\n{s2}\n\n")
            r = self.run_cpp("dnds", p, "-q", *sum((["-m", m] for m in METHODS), []), timeout=120)
            rows = [l.split("\t") for l in r.stdout.splitlines()[1:]]
            self.assertEqual(len(rows), 7)
            for f in rows:
                self.assertIn(f[13], ("ok", "warning"), f)  # (random codons may form stop codons in s2)
                self.assertGreater(float(f[3]), 0)
                self.assertGreater(float(f[4]), 0)

    def test_python_hands_over(self):
        env = dict(os.environ, COAL_CPP=str(self.bin), PYTHONPATH=str(ROOT / "src"))
        r = subprocess.run([sys.executable, "-m", "coal", "dnds", str(KDIR / "edge_cases.axt"), "-q"],
                           capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        direct = self.run_cpp("dnds", KDIR / "edge_cases.axt", "-q").stdout
        self.assertEqual(r.stdout, direct)
        env["COAL_CPP"] = str(Path(tempfile.gettempdir()) / "no-such-coal-cpp")
        r = subprocess.run([sys.executable, "-m", "coal", "kaks", "-i", "x"], capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 1)
        self.assertIn("needs coal-cpp", r.stderr)


if __name__ == "__main__":
    unittest.main()
