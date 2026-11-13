"""The methods Comeron and Li of the Ka/Ks engine (cpp/src/kaks_comeron.cpp, kaks_li.cpp; docs/kaks.md).

* Against the references, stored in tests/kaks/golden (made by tests/kaks/make_comeron_li_fixtures.py):
  orthologr 0.4.3's gestimator code (libsequence) for `-m Comeron --comeron-flavour gestimator`, and
  seqinr::kaks for `-m Li --li-flavour seqinr`, both written with 17 significant digits, on 499 cases made for
  these methods (hand-made edge cases; random short pairs with stop codons, gaps, ambiguity letters, lower case
  and U; genes with terminal stops) and on the edge cases and simulated pairs of the other methods. Every value
  must be the reference's to a relative 1e-12 (they are the same doubles on x86-64), and NA exactly where the
  reference gives its "cannot be calculated" value (999; 9.999999).
* The default flavours: the same numbers wherever COAL's policy reads the pair as the reference does and the
  correction is defined; the documented differences otherwise; regression against COAL's own goldens.
* Names, aliases, statuses, the genetic code, the two output formats, threads.
Skipped if build/coal-cpp (COAL_CPP_BIN) is not built.
"""
from __future__ import annotations

import math
import unittest

import test_kaks as tk   # (the module, not its classes: imported test classes would run again here)

GOLD, HEADER, KDIR = tk.GOLD, tk.HEADER, tk.KDIR

SETS = ("comeron_li_cases", "edge_cases", "sim_small")
STOPS = {"TAA", "TAG", "TGA"}


def val(x: str) -> float:
    x = x.strip()
    return math.nan if x in ("NA", "") else float(x)


def read_axt(path):
    out, rec = {}, []
    for line in path.read_text(encoding="utf-8").splitlines() + [""]:
        if line:
            rec.append(line)
        elif rec:
            out[rec[0]] = rec[1:]
            rec = []
    return out


def policy_applies(a: str, b: str) -> bool:
    """COAL's policy reads the pair differently from the references: a compared column holds a stop codon, or a
    column the references skip holds U (which COAL reads as T)."""
    for i in range(0, len(a), 3):
        x, y = a[i:i + 3].upper(), b[i:i + 3].upper()
        plain = set(x + y) <= set("ACGT")
        if plain and (x in STOPS or y in STOPS):
            return True
        if not plain and "U" in x + y and set((x + y).replace("U", "T")) <= set("ACGT"):
            return True
    return False


class ComeronLiBase(tk.KaKsBase):
    def same(self, got: float, want: float, what: str):
        if math.isnan(want):
            self.assertTrue(math.isnan(got), f"{what}: {got} where the reference has no value")
            return
        self.assertFalse(math.isnan(got), f"{what}: NA where the reference has {want}")
        self.assertLessEqual(abs(got - want), 1e-12 * max(abs(got), abs(want)), f"{what}: {got!r} != {want!r}")

    def tidy(self, axt, *extra):
        r = self.run_cpp("dnds", axt, "--format", "axt", "-q", "--digits", "17", *extra)
        lines = r.stdout.splitlines()
        head = lines[0].split("\t")
        out = {}
        for line in lines[1:]:
            d = dict(zip(head, line.split("\t")))
            name = d["seq1"] if d["seq2"] == "NA" else d["seq1"] + "&" + d["seq2"]
            out[(name, d["method"])] = d
        return out


class TestAgainstGestimator(ComeronLiBase):
    def test_gestimator_flavour_gives_gestimators_numbers(self):
        n_values = 0
        for s in SETS:
            got = self.tidy(KDIR / f"{s}.axt", "-m", "Comeron", "--comeron-flavour", "gestimator")
            for line in (GOLD / f"{s}.gestimator_ref.tsv").read_text(encoding="utf-8").splitlines():
                f = line.split("\t")
                row = got[(f[0], "Comeron")]
                if f[1] == "ERR":       # gestimator cannot take the pair (it reads beyond the rows)
                    self.assertEqual(row["status"], "error", f[0])
                    continue
                for k, col in (("dN", 1), ("dS", 2), ("dNdS", 3)):
                    ref = float(f[col])
                    self.same(val(row[k]), math.nan if ref == 999.0 else ref, f"{s} {f[0]} {k}")
                    n_values += ref != 999.0
                # 999 is never a silent number: the row says why
                if float(f[1]) == 999.0 or float(f[2]) == 999.0:
                    self.assertIn(row["status"], ("saturated", "no_data"), f[0])
        self.assertGreater(n_values, 1200)

    def test_default_flavour_differs_only_as_documented(self):
        for s in SETS:
            seqs = read_axt(KDIR / f"{s}.axt")
            dfl = self.tidy(KDIR / f"{s}.axt", "-m", "Comeron")
            leg = self.tidy(KDIR / f"{s}.axt", "-m", "Comeron", "--comeron-flavour", "gestimator")
            same = 0
            for (name, m), d in dfl.items():
                l = leg[(name, m)]
                if d["status"] == "error":
                    continue
                if policy_applies(*seqs[name]):
                    continue                      # stop codons or U: other columns are compared
                for k in ("dN", "dS"):
                    x, y = val(d[k]), val(l[k])
                    if math.isnan(x) and not math.isnan(y):
                        # the default gives NA where gestimator replaced an undefined logarithm by 0
                        self.assertEqual(d["status"], "saturated", name)
                        self.assertIn("gestimator's value", l["note"], name)
                        self.assertEqual(l["status"], "warning" if "saturated" not in l["status"] else "saturated")
                    else:
                        self.same(x, y, f"{s} {name} {k}")
                        same += 1
            self.assertGreater(same, 30)

    def test_stop_codons(self):
        # gestimator compares a stop codon like any other codon; COAL's policy leaves the column out
        a, b = "GCTAAAGGTCTGTTTGATTAA", "GCAAAGGGTTTGTTCGACTAA"
        d = self.pair(a, b, "-m", "Comeron")
        g = self.pair(a, b, "-m", "Comeron", "--comeron-flavour", "gestimator")
        self.assertEqual((d["codons_compared"], g["codons_compared"]), ("6", "7"))
        # the same pair without its last codon, in the gestimator flavour, is the default's result
        g6 = self.pair(a[:-3], b[:-3], "-m", "Comeron", "--comeron-flavour", "gestimator")
        self.assertEqual(d["dN"], g6["dN"])
        # an internal stop: skipped with a warning by default, compared in the gestimator flavour
        a, b = "ATGGCTTAAGGTCTGTTTGAT", "ATGGCATACGGTTTGTTCGAC"
        d = self.pair(a, b, "-m", "Comeron")
        g = self.pair(a, b, "-m", "Comeron", "--comeron-flavour", "gestimator")
        self.assertIn("internal stop", d["note"])
        self.assertEqual((d["codons_compared"], g["codons_compared"]), ("6", "7"))
        self.assertNotIn("internal stop", g["note"])
        self.assertNotEqual(d["dN"], g["dN"])

    def pair(self, a, b, *extra):
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "p.axt"
            p.write_text(f"x\n{a}\n{b}\n\n", encoding="utf-8")
            return self.tidy(p, *extra)[("x", extra[1] if extra[1] != "C95" else "Comeron")]

    def test_saturation(self):
        # every 4-fold site differs by a transversion: 1 - 2Q <= 0. gestimator prints a number (its logarithm
        # is set to 0) or 999; COAL's default: NA, status saturated
        a, b = "GCTGGTCCTACTGTT" * 6, "GCAGGACCGACAGTG" * 6
        d = self.pair(a, b, "-m", "Comeron")
        self.assertEqual(d["status"], "saturated")
        self.assertEqual(d["dS"], "NA")
        self.assertEqual(d["dN"], "0")
        self.assertIn("Ks saturated", d["note"])


class TestAgainstSeqinr(ComeronLiBase):
    def test_seqinr_flavour_gives_seqinrs_numbers(self):
        n_values = 0
        for s in SETS:
            got = self.tidy(KDIR / f"{s}.axt", "-m", "Li", "--li-flavour", "seqinr")
            for line in (GOLD / f"{s}.seqinr_ref.tsv").read_text(encoding="utf-8").splitlines():
                f = [x.strip() for x in line.split("\t")]
                row = got[(f[0], "Li")]
                if f[1] == "ERR":
                    self.assertEqual(row["status"], "error", f[0])
                    continue
                ka, ks, vka, vks = (float(x) for x in f[1:5])
                for k, ref, var in (("dN", ka, vka), ("dS", ks, vks)):
                    sentinel = ref == 9.999999
                    self.same(val(row[k]), math.nan if sentinel else ref, f"{s} {f[0]} {k}")
                    se = math.nan if sentinel or not var >= 0 else math.sqrt(var)
                    self.same(val(row["SE_" + k]), se, f"{s} {f[0]} SE_{k}")
                    n_values += not sentinel
                if ka == 9.999999 or ks == 9.999999:   # the sentinel is never a number in COAL
                    self.assertIn(row["status"], ("saturated", "no_data"), f[0])
        self.assertGreater(n_values, 700)

    def test_default_flavour_differs_only_by_the_policy(self):
        for s in SETS:
            seqs = read_axt(KDIR / f"{s}.axt")
            dfl = self.tidy(KDIR / f"{s}.axt", "-m", "Li")
            leg = self.tidy(KDIR / f"{s}.axt", "-m", "Li", "--li-flavour", "seqinr")
            for (name, m), d in dfl.items():
                if d["status"] == "error" or policy_applies(*seqs[name]):
                    continue
                for k in ("dN", "dS", "SE_dN", "SE_dS"):
                    self.same(val(d[k]), val(leg[(name, m)][k]), f"{s} {name} {k}")

    def test_li_is_not_lpb(self):
        # the reason Li is a method of its own: seqinr's counting (pathways weighted by amino-acid similarity,
        # its own site classes) is not KaKs_Calculator's LPB / MLPB / LWL
        got = self.tidy(KDIR / "sim_small.axt", "-m", "Li", "-m", "LPB", "-m", "MLPB", "-m", "LWL")
        for other in ("LPB", "MLPB", "LWL"):
            differ = sum(abs(val(got[(n, "Li")]["dS"]) - val(got[(n, other)]["dS"])) > 1e-6
                         for (n, m) in got if m == "Li")
            self.assertGreater(differ, 40, other)


class TestComeronLiRegression(ComeronLiBase):
    def test_goldens(self):
        cmp = tk.TestRegression.compare
        for s in SETS:
            axt = KDIR / f"{s}.axt"
            cmp(self, self.run_cpp("dnds", axt, "--format", "axt", "-q", "--digits", "15", "-m", "Comeron", "-m",
                                   "Li").stdout,
                (GOLD / f"{s}.coal.comeron_li.tsv").read_text(encoding="utf-8"), f"{s} dnds")
            cmp(self, self.run_cpp("kaks", "-i", axt, "-q", "-m", "Comeron", "-m", "Li").stdout,
                (GOLD / f"{s}.coal.comeron_li.kaks").read_text(encoding="utf-8"), f"{s} kaks")


class TestComeronLiCli(ComeronLiBase):
    def test_names_and_aliases(self):
        axt = KDIR / "sim_small.axt"
        want = self.run_cpp("dnds", axt, "-q", "-m", "Comeron").stdout
        for name in ("comeron", "COMERON", "C95", "c95"):
            self.assertEqual(self.run_cpp("dnds", axt, "-q", "-m", name).stdout, want, name)
        self.assertIn("\tComeron\t", want)
        want = self.run_cpp("dnds", axt, "-q", "-m", "Li").stdout
        for name in ("li", "LI"):
            self.assertEqual(self.run_cpp("dnds", axt, "-q", "-m", name).stdout, want, name)
        r = self.run_cpp("dnds", axt, "-m", "XYZ", check=False)
        self.assertIn("Comeron, Li", r.stderr)

    def test_kaks_format(self):
        out = self.run_cpp("kaks", "-i", KDIR / "sim_small.axt", "-q", "-m", "Li", "-m", "Comeron", "-m",
                           "YN").stdout.splitlines()
        self.assertEqual(out[0], HEADER)
        # per pair: KaKs_Calculator's methods first, then Comeron, then Li
        self.assertEqual([l.split("\t")[1] for l in out[1:4]], ["YN", "Comeron", "Li"])
        com, li = out[2].split("\t"), out[3].split("\t")
        self.assertEqual(len(com), 22)
        self.assertEqual(len(li), 22)
        for row in (com, li):
            float(row[2]), float(row[3]), float(row[4])       # Ka, Ks, Ka/Ks
            self.assertEqual(row[5], "NA")                    # no Fisher test: no S / N difference counts
            self.assertEqual(row[7:9], ["NA", "NA"])          # S-Sites, N-Sites
            self.assertEqual(row[18:], ["NA"] * 4)
        self.assertEqual(com[9], "NA")                         # Comeron has four classes of sites, not 0:2:4
        self.assertEqual(len(li[9].split(":")), 3)             # Li: L0:L2:L4
        self.assertEqual(len(li[13].split(":")), 3)

    def test_standard_errors(self):
        got = self.tidy(KDIR / "sim_small.axt", "-m", "Li", "-m", "Comeron")
        for (name, m), d in got.items():
            if m == "Li":
                self.assertGreater(val(d["SE_dN"]), 0)
                self.assertGreater(val(d["SE_dS"]), 0)
            else:
                self.assertEqual((d["SE_dN"], d["SE_dS"]), ("NA", "NA"))

    def test_standard_code_only(self):
        for m in ("Comeron", "Li"):
            for cmd in (("kaks", "-i"), ("dnds",)):
                r = self.run_cpp(*cmd, KDIR / "sim_small.axt", "-m", m, "-c", "2", check=False)
                self.assertEqual(r.returncode, 2, r.stderr)
                self.assertIn("standard genetic code only", r.stderr)
        r = self.run_cpp("dnds", KDIR / "sim_small.axt", "-m", "Comeron", "--comeron-flavour", "x", check=False)
        self.assertEqual(r.returncode, 2)
        r = self.run_cpp("dnds", KDIR / "sim_small.axt", "-m", "Li", "--li-flavour", "x", check=False)
        self.assertEqual(r.returncode, 2)

    def test_statuses(self):
        d = self.tidy(KDIR / "comeron_li_cases.axt", "-m", "Comeron", "-m", "Li")
        st = lambda n, m: d[(n, m)]["status"]
        for m in ("Comeron", "Li"):
            self.assertEqual(st("typical", m), "ok")
            self.assertEqual(st("terminal_stop_both", m), "ok")
            self.assertEqual(st("internal_stop_one_sequence", m), "warning")
            self.assertEqual(st("all_gap", m), "no_data")
            self.assertEqual(st("only_stop_codons", m), "no_data")
            self.assertEqual(st("length_not_multiple_of_3", m), "error")
            # RNA and lower case are read as upper-case DNA
            for k in ("dN", "dS"):
                self.assertEqual(d[("rna", m)][k], d[("typical", m)][k])
                self.assertEqual(d[("lower_case", m)][k], d[("typical", m)][k])
                self.assertEqual(d[("rna_lower_both", m)][k], d[("typical", m)][k])
                self.assertEqual(d[("terminal_stop_both", m)][k], d[("typical", m)][k])
        self.assertEqual(d[("identical", "Comeron")]["dN"], "0")
        self.assertEqual(d[("identical", "Comeron")]["dS"], "0")
        self.assertEqual(d[("identical", "Li")]["dN"], "0")
        self.assertEqual(st("met_trp_only", "Comeron"), "no_data")      # no synonymous site
        self.assertEqual(d[("met_trp_only", "Comeron")]["dS"], "NA")
        self.assertEqual(st("one_codon_identical", "Li"), "no_data")    # Li needs 0-, 2- and 4-fold sites
        self.assertEqual(st("unrelated_saturated", "Li"), "saturated")

    def test_threads_same_output(self):
        for flav in ((), ("--comeron-flavour", "gestimator", "--li-flavour", "seqinr")):
            one = self.run_cpp("dnds", KDIR / "comeron_li_cases.axt", "--format", "axt", "-q", "-j", "1", "-m",
                               "Comeron", "-m", "Li", "-m", "YN", *flav).stdout
            many = self.run_cpp("dnds", KDIR / "comeron_li_cases.axt", "--format", "axt", "-q", "-j", "4", "-m",
                                "Comeron", "-m", "Li", "-m", "YN", *flav).stdout
            self.assertEqual(one, many)

    def test_help_names_the_methods(self):
        for sub in ("kaks", "dnds"):
            h = self.run_cpp(sub, "-h").stdout
            for word in ("Comeron", "C95", "Li", "--comeron-flavour", "gestimator", "--li-flavour", "seqinr"):
                self.assertIn(word, h, sub)


if __name__ == "__main__":
    unittest.main()
