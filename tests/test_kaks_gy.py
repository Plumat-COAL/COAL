"""Method GY of the Ka/Ks engine: maximum likelihood under Goldman & Yang (1994) (cpp/src/kaks_gy.cpp; docs/kaks.md
"GY: maximum likelihood").

* Against PAML codeml (runmode = -2, CodonFreq = 2), stored in tests/kaks/golden/gy_cases.codeml.tsv (made by
  tests/kaks/make_gy_fixtures.py with PAML 4.10.10): simulated pairs of dS 0.1 to 3 and real pairs, and three
  pairs under the vertebrate mitochondrial code. COAL's log-likelihood is never below codeml's (beyond 1e-4:
  codeml prints six decimals); where it is the same maximum (within 1e-3) and codeml is not at one of its
  bounds, dN, dS, omega, kappa, t, S and N are codeml's to its printed precision.
* The model itself, recomputed here from the definition for a small case (pure Python: a Taylor series of
  exp(Qt)), so that the likelihood COAL reports is that of the estimates it reports.
* Edge cases and statuses: identical sequences, one to five codons, synonymous changes only, nonsynonymous
  changes only, saturation, gaps and ambiguity letters, every genetic code; invariances (swap, permutation);
  threads; the two output formats, --ml-columns, names and aliases.
Skipped if build/coal-cpp (COAL_CPP_BIN) is not built.
"""
from __future__ import annotations

import math
import os
import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_kaks as tk   # noqa: E402  (the module, not its classes: imported test classes would run again here)

KDIR, GOLD = tk.KDIR, tk.GOLD
BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AA1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODES = (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23)
COLS = ["seq1", "seq2", "method", "dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "Nd", "Sd",
        "codons_compared", "status", "note", "kappa", "t", "lnL"]


def val(x: str) -> float:
    return math.nan if x in ("NA", "") else float(x)


def read_axt(path: Path) -> list[tuple[str, str, str]]:
    out, rec = [], []
    for line in path.read_text(encoding="utf-8").splitlines() + [""]:
        if line:
            rec.append(line)
        elif rec:
            out.append((rec[0], rec[1], rec[2]))
            rec = []
    return out


def lnl_reference(a: str, b: str, t: float, kappa: float, omega: float, aa: str = AA1) -> float:
    """The log-likelihood of the pair under GY94 with F3x4 frequencies, from the definition."""
    cod = [(a[i:i + 3], b[i:i + 3]) for i in range(0, len(a), 3)]
    f = [[0.0] * 4 for _ in range(3)]
    for x, y in cod:
        for p in range(3):
            f[p][BASES.index(x[p])] += 1
            f[p][BASES.index(y[p])] += 1
    f = [[v / (2 * len(cod)) for v in row] for row in f]
    sense = [c for c in range(64) if aa[c] != "*"]
    pi = [f[0][c // 16] * f[1][(c // 4) % 4] * f[2][c % 4] for c in sense]
    tot = sum(pi)
    pi = [p / tot for p in pi]
    n = len(sense)
    Q = [[0.0] * n for _ in range(n)]
    for i, ci in enumerate(sense):
        for j, cj in enumerate(sense):
            d = [(p, CODONS[ci][p], CODONS[cj][p]) for p in range(3) if CODONS[ci][p] != CODONS[cj][p]]
            if len(d) != 1:
                continue
            x, y = d[0][1], d[0][2]
            q = pi[j]
            if {x, y} in ({"T", "C"}, {"A", "G"}):
                q *= kappa
            if aa[ci] != aa[cj]:
                q *= omega
            Q[i][j] = q
    for i in range(n):
        Q[i][i] = -sum(Q[i])
    mr = -sum(pi[i] * Q[i][i] for i in range(n))
    # exp(Q t / mr) by scaling and squaring of a Taylor series
    s = 12
    M = [[Q[i][j] * t / mr / 2 ** s for j in range(n)] for i in range(n)]

    def mul(A, B):
        Bt = list(zip(*B))
        return [[sum(x * y for x, y in zip(row, col)) for col in Bt] for row in A]

    P = [[float(i == j) for j in range(n)] for i in range(n)]
    term = [row[:] for row in P]
    for k in range(1, 8):
        term = [[v / k for v in row] for row in mul(term, M)]
        P = [[P[i][j] + term[i][j] for j in range(n)] for i in range(n)]
    for _ in range(s):
        P = mul(P, P)
    idx = {CODONS[c]: i for i, c in enumerate(sense)}
    return sum(math.log(pi[idx[x]] * P[idx[x]][idx[y]]) for x, y in cod)


class GyBase(tk.KaKsBase):
    def gy(self, text: str, *extra, sub="dnds") -> list[dict]:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "in.axt"
            p.write_text(text, encoding="utf-8", newline="\n")
            r = self.run_cpp("dnds", p, "-m", "GY", "--ml-columns", "--digits", "12", "-q", *extra, check=False)
        self.assertIn(r.returncode, (0, 2), r.stderr)
        lines = r.stdout.splitlines()
        self.assertEqual(lines[0].split("\t"), COLS)
        return [dict(zip(COLS, ln.split("\t"))) for ln in lines[1:]]

    def one(self, a: str, b: str, *extra) -> dict:
        return self.gy(f"x&y\n{a}\n{b}\n\n", *extra)[0]


class TestAgainstCodeml(GyBase):
    def test_golden(self):
        gold = {}
        for line in (GOLD / "gy_cases.codeml.tsv").read_text(encoding="utf-8").splitlines()[1:]:
            f = line.split("\t")
            gold[f[0]] = dict(zip(("icode", "lnL", "t", "kappa", "omega", "dN", "dS", "S", "N", "SE_dN", "SE_dS"),
                                  map(float, f[1:])))
        pairs = read_axt(KDIR / "gy_cases.axt")
        self.assertEqual(len(pairs), len(gold))
        self.assertGreaterEqual(len(pairs), 20)
        same = n_se = 0
        for icode, table in ((0, "1"), (1, "2")):
            sub = [p for p in pairs if gold[p[0].split("&")[0]]["icode"] == icode]
            rows = self.gy("".join(f"{n}\n{a}\n{b}\n\n" for n, a, b in sub), "-c", table)
            self.assertEqual(len(rows), len(sub))
            for (name, _, _), r in zip(sub, rows):
                g = gold[name.split("&")[0]]
                lnl = val(r["lnL"])
                with self.subTest(pair=name):
                    self.assertFalse(math.isnan(lnl), r)
                    # never a worse maximum than codeml's
                    self.assertGreaterEqual(lnl, g["lnL"] - 1e-4)
                    at_bound = (g["t"] >= 49.99 or g["kappa"] <= 0.4001 or g["omega"] <= 0.00101 or
                                g["omega"] >= 98.99)
                    if at_bound or lnl > g["lnL"] + 1e-3:
                        continue   # codeml stopped at one of its bounds or at a lower maximum
                    same += 1
                    self.assertAlmostEqual(val(r["S_sites"]), g["S"], delta=0.051)
                    self.assertAlmostEqual(val(r["N_sites"]), g["N"], delta=0.051)
                    if r["status"] == "saturated":
                        # COAL gives no dS above 10 (codeml prints it); dN is still codeml's. The likelihood is
                        # nearly flat in t there (codeml's own t moves in the fourth digit from run to run)
                        self.assertGreater(g["dS"], 10)
                        self.assertEqual(r["dS"], "NA")
                        self.assertAlmostEqual(val(r["dN"]), g["dN"], delta=6e-5 + 1e-3 * g["dN"])
                        self.assertAlmostEqual(val(r["t"]), g["t"], delta=1e-3 * g["t"])
                        continue
                    self.assertAlmostEqual(val(r["t"]), g["t"], delta=2e-5 + 1e-4 * g["t"])
                    self.assertAlmostEqual(val(r["kappa"]), g["kappa"], delta=2e-5 + 1e-4 * g["kappa"])
                    self.assertEqual(r["status"], "ok")
                    self.assertAlmostEqual(val(r["dN"]), g["dN"], delta=6e-5 + 1e-4 * g["dN"])
                    self.assertAlmostEqual(val(r["dS"]), g["dS"], delta=6e-5 + 1e-4 * g["dS"])
                    self.assertAlmostEqual(val(r["dNdS"]), g["omega"], delta=6e-6 + 1e-4 * g["omega"])
                    # standard errors: codeml's Hessian is made by differences of the likelihood and its SEs
                    # are good to about 15 % (measured: up to 16 % from the independent reference of
                    # test_standard_errors below, which COAL matches to four digits)
                    if 0.001 < g["omega"] < 99 and g["kappa"] > 0.4:
                        for c in ("SE_dN", "SE_dS"):
                            self.assertAlmostEqual(val(r[c]), g[c], delta=2e-5 + 0.25 * g[c], msg=c)
                            n_se += 1
        self.assertGreaterEqual(same, 14)   # (most pairs are compared number by number)
        self.assertGreaterEqual(n_se, 20)

    def test_likelihood_is_the_models(self):
        """lnL is the log-likelihood of the reported t, kappa and omega (recomputed from the definition), dN / dS
        is omega, and no nearby parameter value has a higher likelihood."""
        rng = random.Random(5)
        sense = [c for c in CODONS if AA1[CODONS.index(c)] != "*"]
        a = [rng.choice(sense) for _ in range(40)]
        b = []
        for c in a:
            if rng.random() < 0.4:
                p = rng.randrange(3)
                d = c[:p] + rng.choice(BASES) + c[p + 1:]
                c = d if AA1[CODONS.index(d)] != "*" else c
            b.append(c)
        a, b = "".join(a), "".join(b)
        r = self.one(a, b)
        self.assertEqual(r["status"], "ok")
        t, k, w, lnl = val(r["t"]), val(r["kappa"]), val(r["dNdS"]), val(r["lnL"])
        self.assertAlmostEqual(lnl_reference(a, b, t, k, w), lnl, delta=1e-7)
        for dt, dk, dw in ((1.02, 1, 1), (0.98, 1, 1), (1, 1.02, 1), (1, 0.98, 1), (1, 1, 1.02), (1, 1, 0.98)):
            self.assertLess(lnl_reference(a, b, t * dt, k * dk, w * dw), lnl + 1e-9)
        # dS = t rhoS / (3 rhoS1) and S + N = 3 x codons
        self.assertAlmostEqual(val(r["S_sites"]) + val(r["N_sites"]), 120, delta=1e-8)
        self.assertAlmostEqual((val(r["S_sites"]) * val(r["dS"]) + val(r["N_sites"]) * val(r["dN"])) / 120, t / 3,
                               delta=1e-9)


class TestStandardErrors(GyBase):
    def test_standard_errors(self):
        """SE_dN and SE_dS against an independent computation (tests/kaks/make_gy_se_reference.py: the likelihood
        rebuilt with numpy, its Hessian by second differences of the likelihood, the delta method)."""
        ref = {}
        for line in (GOLD / "gy_cases.se_reference.tsv").read_text(encoding="utf-8").splitlines()[1:]:
            f = line.split("\t")
            ref[f[0]] = tuple(map(float, f[1:]))
        self.assertGreaterEqual(len(ref), 12)
        pairs = [p for p in read_axt(KDIR / "gy_cases.axt") if p[0].split("&")[0] in ref]
        rows = self.gy("".join(f"{n}\n{a}\n{b}\n\n" for n, a, b in pairs))
        for (name, _, _), r in zip(pairs, rows):
            lnl, se_n, se_s = ref[name.split("&")[0]]
            with self.subTest(pair=name):
                self.assertAlmostEqual(val(r["lnL"]), lnl, delta=2e-6)
                self.assertAlmostEqual(val(r["SE_dN"]), se_n, delta=1e-3 * se_n)
                self.assertAlmostEqual(val(r["SE_dS"]), se_s, delta=1e-3 * se_s)

    def test_no_standard_error_at_a_bound_or_when_saturated(self):
        pairs = {n.split("&")[0]: (a, b) for n, a, b in read_axt(KDIR / "gy_cases.axt")}
        # no nonsynonymous difference: omega at its lower bound; dS has an SE (about dS / sqrt(4 differences))
        r = self.one(*pairs["real0_7"])
        self.assertEqual((r["status"], r["SE_dN"]), ("ok", "NA"))
        self.assertLess(val(r["dNdS"]), 1e-4)
        self.assertTrue(0.3 * val(r["dS"]) < val(r["SE_dS"]) < val(r["dS"]))
        for name in ("real1_2", "real1_12", "sim60"):
            r = self.one(*pairs[name])
            self.assertEqual((r["status"], r["SE_dN"], r["SE_dS"]), ("saturated", "NA", "NA"))
        r = self.one("ATGGCTAAAGGTCCC", "ATGGCTAAAGGTCCC")
        self.assertEqual((r["SE_dN"], r["SE_dS"]), ("NA", "NA"))


class TestTinyProbabilities(GyBase):
    def test_one_triple_change_among_thousands_of_identical_codons(self):
        """A probability below the rounding error of the eigen sum (codons three substitutions apart at t about
        1e-4: 1e-13 and less) is computed by the series of non-negative terms: the likelihood is the model's
        (recomputed here from the definition), t is not at its lower bound, and dN is NG's size. Before the fix
        this pair had t = 1e-6 (the bound) and a dN 100 times too low, with status ok (audit of #92, MEDIUM-1)."""
        rng = random.Random(8)
        sense = [c for c in CODONS if AA1[CODONS.index(c)] != "*"]
        a = [rng.choice(sense) for _ in range(20000)]
        b = list(a)
        k = next(i for i, c in enumerate(a) if c == "AAA")
        b[k] = "CCC"                      # three substitutions (Lys -> Pro)
        j = next(i for i, c in enumerate(a) if c == "GCT")
        b[j] = "GCC"                      # and one synonymous change
        a, b = "".join(a), "".join(b)
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "in.axt"
            p.write_text(f"x&y\n{a}\n{b}\n\n", encoding="utf-8", newline="\n")
            out = self.run_cpp("dnds", p, "-m", "GY", "-m", "NG", "--ml-columns", "--digits", "12", "-q").stdout
            swapped = Path(td) / "sw.axt"
            swapped.write_text(f"x&y\n{b}\n{a}\n\n", encoding="utf-8", newline="\n")
            out2 = self.run_cpp("dnds", swapped, "-m", "GY", "-m", "NG", "--ml-columns", "--digits", "12", "-q").stdout
        self.assertEqual(out, out2)
        rows = {ln.split("\t")[2]: dict(zip(COLS, ln.split("\t"))) for ln in out.splitlines()[1:]}
        gy, ng = rows["GY"], rows["NG"]
        self.assertEqual(gy["status"], "ok")
        t, kappa, lnl = val(gy["t"]), val(gy["kappa"]), val(gy["lnL"])
        self.assertGreater(t, 5e-5)       # (4 substitutions in 20,000 codons: 2e-4; the lower bound is 1e-6)
        self.assertLess(t, 1e-3)
        self.assertTrue(0.5 < val(gy["dN"]) / val(ng["dN"]) < 2, (gy["dN"], ng["dN"]))
        omega = val(gy["dNdS"])
        ref = lnl_reference(a, b, t, kappa, omega)
        self.assertAlmostEqual(ref, lnl, delta=1e-6 * abs(ref) * 1e-3)
        # and it is a maximum of the true likelihood in t
        for f in (0.8, 1.25):
            self.assertLess(lnl_reference(a, b, t * f, kappa, omega), ref)


class TestEdgeCases(GyBase):
    GENE = "ATGGCTAAAGGTCCCGAACTGTTTACCGGCGTTGTGCCGATTCTGGTTGAACTGGATGGCGATGTTAACGGCCACAAATTC"

    def codonwise(self, change: dict) -> str:
        return "".join(change.get(self.GENE[i:i + 3], self.GENE[i:i + 3]) for i in range(0, len(self.GENE), 3))

    def test_identical(self):
        r = self.one(self.GENE, self.GENE)
        self.assertEqual((r["status"], r["dN"], r["dS"], r["dNdS"], r["Nd"], r["Sd"], r["t"]),
                         ("ok", "0", "0", "NA", "0", "0", "0"))
        self.assertEqual(r["kappa"], "NA")   # no change: nothing to estimate kappa from
        self.assertLess(val(r["lnL"]), 0)

    def test_one_to_five_codons(self):
        cases = [("AAA", "AAG"), ("ATG", "ATG"), ("ATG", "TGG"), ("GCT", "GCA"), ("AAAGCT", "AAGGCA"),
                 ("ATGGCTAAA", "ATGGCAAAG"), ("ATGGCTAAAGGT", "ATGTCAAAGGGA"), ("ATGGCTAAAGGTCCC", "TTGGCAAGAGGACCA"),
                 ("TTT", "GGG"), ("ATGTGG", "TGGATG")]
        rows = self.gy("".join(f"c{i}&x\n{a}\n{b}\n\n" for i, (a, b) in enumerate(cases)))
        self.assertEqual(len(rows), len(cases))
        for r in rows:
            self.assertIn(r["status"], ("ok", "saturated", "no_data", "not_converged"), r)
            for c in COLS[3:13] + COLS[15:]:
                self.assertNotIn(r[c].lower(), ("nan", "inf", "-inf", "-nan", ""), r)
            for c in ("dN", "dS"):
                if r[c] != "NA":
                    self.assertGreaterEqual(val(r[c]), 0)
                    self.assertLessEqual(val(r[c]), 10)
        # Met / Trp codons only: no synonymous site, as for YN
        self.assertEqual((rows[2]["status"], rows[2]["dS"]), ("no_data", "NA"))
        self.assertIn("no synonymous sites", rows[2]["note"])
        # One codon with one synonymous change: every site of the one kind there is differs, the likelihood is
        # flat in t up to its bound and there is no estimate (not the 4.6 at which the optimiser's tolerance
        # would stop); t is reported at the bound
        r = rows[0]
        self.assertEqual((r["status"], r["dN"], r["dS"], r["t"]), ("no_data", "NA", "NA", "50"))
        self.assertIn("the likelihood is flat in the branch length", r["note"])
        # ... and a pair that is not flat keeps its estimate, with the standard error that says how little it is
        self.assertEqual(rows[7]["status"], "ok")
        self.assertGreater(val(rows[7]["SE_dS"]), 0.5 * val(rows[7]["dS"]))

    def test_synonymous_changes_only(self):
        b = self.codonwise({"GCT": "GCA", "GGT": "GGC", "CCC": "CCA", "CTG": "CTA"})
        r = self.one(self.GENE, b)
        self.assertEqual(r["status"], "ok")
        self.assertGreater(val(r["dS"]), 0.05)
        self.assertLess(val(r["dN"]), 1e-5)       # omega at its lower bound, 1e-6
        self.assertLess(val(r["dNdS"]), 1.1e-6)

    def test_nonsynonymous_changes_only(self):
        b = self.codonwise({"AAA": "ACA", "GAA": "GTA", "GAT": "GCT"})
        r = self.one(self.GENE, b)
        self.assertEqual(r["status"], "ok")
        self.assertGreater(val(r["dN"]), 0.05)
        self.assertLess(val(r["dS"]), val(r["dN"]) / 900)   # omega at its upper bound, 999

    def test_saturated(self):
        """An estimate above 10 substitutions per site is NA (status saturated) and the other one is kept; t at
        its bound of 50 always has one. The three such pairs of the fixture (codeml: dS 14.8 with t 10.6; dS 59.9
        with t = 50; and one where both are above 10)."""
        pairs = {n.split("&")[0]: (a, b) for n, a, b in read_axt(KDIR / "gy_cases.axt")}
        for name, dn_kept, at_bound in (("real1_2", True, False), ("real1_12", True, True), ("sim60", False, True)):
            r = self.one(*pairs[name])
            with self.subTest(pair=name):
                self.assertEqual((r["status"], r["dS"], r["dNdS"]), ("saturated", "NA", "NA"))
                self.assertIn("Ks saturated (the maximum-likelihood estimate, ", r["note"])
                self.assertIn("is above 10 substitutions per site", r["note"])
                self.assertEqual("upper bound, t = 50" in r["note"], at_bound)
                self.assertEqual(r["dN"] != "NA", dn_kept)
                self.assertEqual("Ka saturated" in r["note"], not dn_kept)
                if dn_kept:
                    self.assertTrue(0 < val(r["dN"]) < 10)
                self.assertNotEqual(r["lnL"], "NA")   # the likelihood and kappa are still reported
        # --max-ds flags an estimate that is not NA
        r = self.one(*pairs["real1_7"], "--max-ds", "2")
        self.assertEqual(r["status"], "warning")
        self.assertIn("dS above 2", r["note"])

    def test_gaps_ambiguity_and_stops_as_the_other_methods(self):
        """Columns with a gap, an ambiguity letter or a stop codon are skipped, as for YN: the result is that of
        the pair without them (and lower case / U are read as DNA)."""
        b = self.codonwise({"AAA": "AAG", "GAA": "GTA", "GGT": "GGC"})
        clean = self.one(self.GENE, b)
        a2 = self.GENE[:9] + "---" + "NNN" + "A-G" + "TAA" + self.GENE[9:] + "TAA"
        b2 = b[:9] + "GCT" + "GCT" + "GCT" + "GCT" + b[9:] + "TAG"
        dirty = self.one(a2.lower().replace("t", "u"), b2)
        self.assertEqual(dirty["status"], "warning")   # the internal stop codon
        for c in ("dN", "dS", "dNdS", "kappa", "t", "lnL", "N_sites", "S_sites", "codons_compared"):
            self.assertEqual(dirty[c], clean[c], c)

    def test_every_genetic_code(self):
        rng = random.Random(3)
        a = "".join(rng.choice(CODONS) for _ in range(120))
        b = "".join(c if rng.random() < 0.7 else c[:2] + rng.choice(BASES)
                    for c in (a[i:i + 3] for i in range(0, len(a), 3)))
        for code in CODES:
            r = self.one(a, b, "-c", str(code))
            with self.subTest(code=code):
                self.assertIn(r["status"], ("ok", "warning"), r)
                self.assertGreater(val(r["dS"]), 0)
                self.assertAlmostEqual(val(r["dNdS"]), val(r["dN"]) / val(r["dS"]), delta=1e-9)

    def test_swap_and_permutation_give_the_same_bytes(self):
        pairs = read_axt(KDIR / "gy_cases.axt")[:6]
        text = "".join(f"{n}\n{a}\n{b}\n\n" for n, a, b in pairs)
        swapped = "".join(f"{n}\n{b}\n{a}\n\n" for n, a, b in pairs)
        rng = random.Random(2)
        perm = ""
        for n, a, b in pairs:
            order = list(range(len(a) // 3))
            rng.shuffle(order)
            pa = "".join(a[3 * i:3 * i + 3] for i in order)
            pb = "".join(b[3 * i:3 * i + 3] for i in order)
            perm += f"{n}\n{pa}\n{pb}\n\n"
        base = self.gy(text, "--digits", "17")
        self.assertEqual(self.gy(swapped, "--digits", "17"), base)
        self.assertEqual(self.gy(perm, "--digits", "17"), base)


class TestCliAndFormats(GyBase):
    def test_threads_are_deterministic(self):
        outs = set()
        for j in ("1", "2", "4", "4"):
            outs.add(self.run_cpp("dnds", KDIR / "gy_cases.axt", "-m", "GY", "-m", "YN", "--digits", "17",
                                  "--ml-columns", "-q", "-j", j).stdout)
        self.assertEqual(len(outs), 1)

    def test_names_order_and_default_columns(self):
        r = self.run_cpp("dnds", KDIR / "sim_small.axt", "-m", "yn", "-m", "gy", "-m", "NG", "-q")
        lines = r.stdout.splitlines()
        self.assertEqual(lines[0].split("\t"), COLS[:15])   # no extra column without --ml-columns
        self.assertEqual([ln.split("\t")[2] for ln in lines[1:4]], ["NG", "GY", "YN"])   # KaKs_Calculator's order
        for alias in ("GY-HKY", "gy94"):
            a = self.run_cpp("dnds", KDIR / "gy_cases.axt", "-m", alias, "-q").stdout
            self.assertEqual(a, self.run_cpp("dnds", KDIR / "gy_cases.axt", "-m", "GY", "-q").stdout)
        # standard errors: for an estimate inside the bounds, NA for a saturated pair
        out = self.run_cpp("dnds", KDIR / "gy_cases.axt", "-m", "GY", "-q").stdout
        rows = [ln.split("\t") for ln in out.splitlines()[1:]]
        self.assertTrue(0 < float(rows[0][6]) < float(rows[0][3]) and 0 < float(rows[0][7]) < float(rows[0][4]))
        for row in rows:
            if row[13] == "saturated":
                self.assertEqual(row[6:8], ["NA", "NA"], row[0])

    def test_kaks_format(self):
        r = self.run_cpp("kaks", "-i", KDIR / "gy_cases.axt", "-m", "GY", "-m", "YN", "-q", "--digits", "10")
        lines = r.stdout.splitlines()
        self.assertEqual(lines[0] + "\n", tk.HEADER if tk.HEADER.endswith("\n") else lines[0] + "\n")
        head = lines[0].split("\t")
        gy = dict(zip(head, lines[1].split("\t")))
        yn = dict(zip(head, lines[2].split("\t")))
        self.assertEqual((gy["Method"], gy["Model"], yn["Method"], yn["Model"], yn["ML-Score"], yn["AICc"]),
                         ("GY", "HKY", "YN", "NA", "NA", "NA"))
        d = self.run_cpp("dnds", KDIR / "gy_cases.axt", "-m", "GY", "-q", "--digits", "10", "--ml-columns")
        t = dict(zip(COLS, d.stdout.splitlines()[1].split("\t")))
        self.assertEqual((gy["Ka"], gy["Ks"], gy["Ka/Ks"], gy["ML-Score"]), (t["dN"], t["dS"], t["dNdS"], t["lnL"]))
        n, lnl = int(t["codons_compared"]), float(t["lnL"])
        self.assertAlmostEqual(float(gy["AICc"]), -2 * lnl + 6 * n / (n - 4), delta=1e-6)
        self.assertAlmostEqual(float(gy["Divergence-Time"]), float(t["t"]) / 3, delta=1e-9)
        self.assertEqual(gy["Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)"],
                         f"{t['kappa']}:{t['kappa']}:1:1:1:1")
        self.assertAlmostEqual(float(gy["S-Substitutions"]) + float(gy["N-Substitutions"]),
                               float(gy["Substitutions"]), delta=1e-6)

    def test_ml_columns_for_the_other_methods(self):
        r = self.run_cpp("dnds", KDIR / "gy_cases.axt", "-m", "YN", "-m", "NG", "--ml-columns", "-q")
        rows = [dict(zip(COLS, ln.split("\t"))) for ln in r.stdout.splitlines()[1:3]]
        self.assertEqual([x["method"] for x in rows], ["NG", "YN"])
        self.assertEqual((rows[0]["t"], rows[0]["lnL"], rows[1]["t"], rows[1]["lnL"]), ("NA", "NA", "NA", "NA"))
        self.assertNotEqual(rows[1]["kappa"], "NA")
        # the option belongs to `coal dnds`
        self.assertEqual(self.run_cpp("kaks", "-i", KDIR / "gy_cases.axt", "--ml-columns", check=False).returncode, 2)


if __name__ == "__main__":
    unittest.main()
