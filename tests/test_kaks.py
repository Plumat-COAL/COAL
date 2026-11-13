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
# Pairs where KaKs_Calculator loses a defined estimate to its own NaN arithmetic and COAL gives the number
# (docs/kaks.md, differences 8-10). Only these (pair, method) may have a number where KaKs_Calculator prints nan
# (or, for no_syn_sites_multi_step, the 0 its lost pathway weights leave). Their values are checked in
# TestOneKindOfSite and TestMynNan.
KK_LOSES = {
    ("met_trp_only", "YN"), ("met_trp_only", "MYN"),                        # S = 0, two differences per codon
    ("no_syn_sites_multi_step", "YN"), ("no_syn_sites_multi_step", "MYN"),
    ("no_syn_sites_single_step", "MYN"), ("no_nonsyn_sites_one_codon", "MYN"), ("no_nonsyn_sites", "MYN"),
    ("short_one_codon", "MYN"), ("myn_tn93_zero_division", "MYN"), ("identical_short", "MYN"),   # TN93 0/0
    ("myn_nan_kappa", "MYN"),                                               # 0 x NaN in the mean kappa
    ("myn_sites_in_one_sequence", "MYN"),                                   # 0/0 site base frequencies
}
KK_NAN = ("nan", "-nan", "inf", "-inf")
# KaKs_Calculator-format columns that are NA when the YN / MYN iteration did not converge: Ka, Ks, Ka/Ks, the
# Fisher p, S- and N-Substitutions, Divergence-Time. Everything else in the row does not depend on the iteration.
NOT_CONVERGED_NA = (2, 3, 4, 5, 11, 12, 15)


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
            if os.environ.get("COAL_REQUIRE_CPP"):   # CI's C++ job: the Ka/Ks engine must be tested, not skipped
                raise AssertionError("coal-cpp not built (COAL_REQUIRE_CPP is set)")
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
            ref = read_rows((GOLD / f"{name}.KaKs_Calculator.{m}.kaks").read_text(encoding="utf-8"))
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
                    # the iteration did not converge: no estimate, where KaKs_Calculator prints its last round
                    # (docs/kaks.md, difference 12). Only the fields that come out of the iteration
                    if y == "NA" and st == "not_converged" and i in NOT_CONVERGED_NA:
                        why = st
                    elif i == 17 and pair in RNA_PAIRS:
                        why = "lowercase GC"  # KaKs_Calculator counts upper-case G and C only
                    elif x == "-0" and y == "0":
                        why = "negative zero"
                    elif x in KK_NAN and y == "NA":
                        why = "nan"
                    elif (pair, m) in KK_LOSES:  # the row is lost to NaN there; COAL's values are tested below
                        self.assertTrue(any(v in KK_NAN or "nan" in v for v in rf) or rf[2] == "NA", rf)
                        why = "KaKs_Calculator loses the value"
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
                    # nothing to estimate: KaKs_Calculator prints nan or a default there. A number that
                    # KaKs_Calculator gives and COAL does not is never excused by the status.
                    if why is None and st == "no_data" and not (y == "NA" and num(x) is not None and num(x) != 0):
                        why = st
                    self.assertIsNotNone(why, f"{name} {key} field {i} ({HEADER.split(chr(9))[i]}): "
                                              f"KaKs_Calculator {x!r}, COAL {y!r} (status {st})")
                    counts[why] = counts.get(why, 0) + 1
        return counts

    def test_edge_cases(self):
        c = self.check_set("edge_cases")
        self.assertGreater(c.get("saturated", 0), 0)

    def test_simulated(self):
        c = self.check_set("sim_small")
        # at most one pair (YN, sim_w0.1_t2_k2_3: dS above 4, ten rounds are not enough) has no estimate; its 7 fields
        self.assertLessEqual(c.get("not_converged", 0), len(NOT_CONVERGED_NA))


class TestOneKindOfSite(KaKsBase):
    """YN and MYN when there are no nonsynonymous (or no synonymous) sites: the estimate that has sites is still
    made, the other is NA and the status is no_data. KaKs_Calculator's YN gives the same numbers when every codon
    differs at one position at most (checked field by field in TestAgainstKaKsCalculator: these pairs are not in
    KK_LOSES for YN); the expected values here are PAML yn00's (4.10.10, 12 digits; weighting = 0 for the pair
    with two differences per codon, where the weighted run of both references is lost to a NaN t)."""

    def setUp(self):
        self.rows = self.dnds(KDIR / "edge_cases.axt", "-m", "YN", "-m", "MYN")

    def check(self, pair, dn, se_dn, ds, se_ds, note):
        for m in ("YN", "MYN"):
            f = self.rows[(pair, m)]
            for got, want, what in ((f[3], dn, "dN"), (f[6], se_dn, "SE_dN"), (f[4], ds, "dS"), (f[7], se_ds, "SE_dS")):
                if want is None:
                    self.assertEqual(got, "NA", f"{pair} {m} {what}")
                else:
                    self.assertTrue(close(float(got), want, 1e-9), f"{pair} {m} {what}: {got} != {want}")
            self.assertEqual(f[5], "NA", f"{pair} {m} dN/dS")
            self.assertEqual((f[13], f[14]), ("no_data", note), f"{pair} {m}")

    def test_no_nonsynonymous_sites(self):
        self.check("no_nonsyn_sites_one_codon", None, None, 0.477385626221, 0.612372435696, "no nonsynonymous sites")
        self.check("no_nonsyn_sites", None, None, 0.440839998677, 0.0346410161514, "no nonsynonymous sites")

    def test_no_synonymous_sites(self):
        self.check("no_syn_sites_multi_step", 0.513345885278, 0.357217254156, None, None, "no synonymous sites")
        # one difference per codon: KaKs_Calculator's YN gives this Ka (golden file), and MYN equals YN
        kk = read_rows((GOLD / "edge_cases.KaKs_Calculator.YN.kaks").read_text(encoding="utf-8"))[("no_syn_sites_single_step", "YN")]
        self.check("no_syn_sites_single_step", float(kk[2]), 0.40654667384577, None, None, "no synonymous sites")

    def test_no_sites_at_all(self):
        for m in ("YN", "MYN"):
            f = self.rows[("all_gap", m)]
            self.assertEqual((f[3], f[4], f[13]), ("NA", "NA", "no_data"))

    def test_kaks_format(self):
        r = self.run_cpp("kaks", "-i", KDIR / "edge_cases.axt", "-q", "-m", "YN")
        f = read_rows(r.stdout)[("no_nonsyn_sites", "YN")]
        self.assertEqual(f[2:5], ["NA", "0.44084", "NA"])   # KaKs_Calculator: -nan  0.44084  -nan


class TestMynNan(KaKsBase):
    """MYN rows that KaKs_Calculator prints as nan although the estimate is defined (docs/kaks.md, difference
    10). A KaKs_Calculator built with the two one-line corrections named there agrees with COAL on these pairs
    and on 4,000 random short pairs with one changed position per codon (benchmarks/kaks/random_short_pairs.py
    --mode single; for codons that differ at several positions see TestMultiPositionCodons)."""

    def setUp(self):
        self.rows = self.dnds(KDIR / "edge_cases.axt", "-m", "YN", "-m", "MYN")

    def test_tn93_zero_division_is_the_f84_distance(self):
        # no estimable kappa (both 2) and no pathway to weight: MYN must equal YN exactly
        for pair in ("myn_tn93_zero_division", "short_one_codon", "identical_short", "myn_sites_in_one_sequence"):
            y, m = self.rows[(pair, "YN")], self.rows[(pair, "MYN")]
            self.assertEqual(m[13], "ok", (pair, m))
            self.assertEqual(m[14], "")
            for col in (3, 4, 6, 7, 8, 9, 10, 11):
                self.assertTrue(close(float(y[col]), float(m[col]), 1e-9), f"{pair} col {col}: {y[col]} vs {m[col]}")
        self.assertEqual(float(self.rows[("myn_tn93_zero_division", "MYN")][3]), 0.0)
        self.assertEqual(float(self.rows[("myn_sites_in_one_sequence", "MYN")][4]), 0.0)

    def test_nan_kappa(self):
        # kappa comes from the class of sites that has one: the corrected KaKs_Calculator prints these values
        m = self.rows[("myn_nan_kappa", "MYN")]
        self.assertEqual(m[13], "ok")
        for col, want in ((3, 0.136556085877548), (4, 0.97117192633034), (8, 16.2479079015356), (9, 7.75209209846442)):
            self.assertTrue(close(float(m[col]), want, 1e-9), f"col {col}: {m[col]} != {want}")
        kk = read_rows((GOLD / "edge_cases.KaKs_Calculator.MYN.kaks").read_text(encoding="utf-8"))[("myn_nan_kappa", "MYN")]
        self.assertIn("nan", kk[2])  # the unmodified reference: the whole row is nan


class TestMultiPositionCodons(KaKsBase):
    """Pairs whose codons differ at two or three positions, as random short pairs do (independent review
    2026-11-13, finding 6; docs/kaks.md "Short and saturated pairs"). The reference values are recorded here:
    KaKs_Calculator 2.0 at 15 digits, one pair per process (unpatched: KK; with the two corrected lines: KKFIX)
    and PAML yn00 at 12 digits with weighting 0."""

    PAIRS = {
        # no synonymous site, one codon differs at two positions
        "one_kind_multi": ("CCTCTTTCT", "TTTCCTTTT"),
        # 12 codons: KaKs_Calculator's MYN meets 0/0 in TN93 in a round and from then on drops the codons that
        # differ at more than one position
        "myn_nan_round": ("TGCATGATATCTTTTAATAACCGAGAAGTTTCCGAC", "TCCATGATATCTTTTAATTCTCGAGAATTTTCCGAC"),
        # every synonymous site differs (Sd = S): saturated, and the iteration takes another path than the reference
        "all_syn_sites_differ": ("ATGAAAAAGAAGTGGTGGGCA", "TGGTGGGCTGCTGCTTGGATG"),
        # the iteration does not settle
        "myn_not_converged": ("GGTGTTAGGCTTATAACGACAGCTGGAGTGGGC", "ATCACTCCCCTACTCATCACACAAGGCATAAAA"),
        "yn_not_converged": ("AAA", "ATG"),
        "sat_not_converged": ("TGGATGGCT", "TGGAAGTGG"),
    }

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        axt = Path(self._td.name) / "multi.axt"
        axt.write_text("".join(f"{n}\n{a}\n{b}\n\n" for n, (a, b) in self.PAIRS.items()), encoding="utf-8")
        self.rows = self.dnds(axt, "-m", "YN", "-m", "MYN")
        self.paml0 = self.dnds(axt, "-m", "YN", "--yn-flavour", "paml", "--yn-weighting", "0")
        self.axt = axt

    def tearDown(self):
        self._td.cleanup()

    def test_one_kind_of_site_counts_every_difference(self):
        """Ka is JC69 of all 4 differences at the 9 sites, as yn00 (weighting 0) gives it. KaKs_Calculator prints
        0.263548415128416 for YN and MYN: the two differences of the codon that differs at one position only
        (its pathway weights are NaN, and the codon with two differences is dropped without a word)."""
        want = -0.75 * math.log(1 - (4 / 3) * (4 / 9))
        self.assertTrue(close(want, 0.673456194904, 1e-11))            # yn00, weighting = 0
        kk_wrong = -0.75 * math.log(1 - (4 / 3) * (2 / 9))
        self.assertTrue(close(kk_wrong, 0.263548415128416, 1e-12))     # KaKs_Calculator, YN and MYN
        for m in ("YN", "MYN"):
            r = self.rows[("one_kind_multi", m)]
            self.assertTrue(close(float(r[3]), want, 1e-11), r)
            self.assertEqual((r[4], r[13], r[14]), ("NA", "no_data", "no synonymous sites"))
        self.assertTrue(close(float(self.paml0[("one_kind_multi", "YN")][3]), want, 1e-11))

    def test_myn_keeps_the_codons_with_several_differences(self):
        """KK MYN: Ka 0.0781105978289474, Ks -nan (and COAL before the review fixes: 0.0781, NA). KKFIX MYN:
        Ka 0.167317060572902, Ks 0.123094343743785. KK YN: 0.167188448579, 0.123094343744."""
        m, y = self.rows[("myn_nan_round", "MYN")], self.rows[("myn_nan_round", "YN")]
        self.assertEqual((m[13], m[14]), ("ok", ""))
        self.assertTrue(close(float(m[3]), 0.167317060572902, 1e-9), m)
        self.assertTrue(close(float(m[4]), 0.123094343743785, 1e-9), m)
        self.assertTrue(close(float(y[3]), 0.167188448579, 1e-9), y)
        self.assertGreater(float(m[3]), 2 * 0.0781105978289474)

    def test_saturated_kind_and_the_other_estimate(self):
        """Sd = S: the proportion of differences is 1 up to rounding, and rounding decides whether a round sees
        "P + Q > 1" (t = 99) or JC69's cap. Either way Ks is NA / saturated in COAL (KK and KKFIX print the cap,
        0.688156925315808). Ka is the same fixed point, reached by another path: KK 2.36633343882711, COAL
        within the iteration's tolerance of it, not to 1e-9."""
        m = self.rows[("all_syn_sites_differ", "MYN")]
        self.assertEqual((m[4], m[13]), ("NA", "saturated"))
        self.assertNotIn("did not converge", m[14])
        self.assertTrue(close(float(m[3]), 2.36633343882711, 1e-6), m)

    # tidy columns: dN dS dNdS SE_dN SE_dS (3..7), Nd Sd (10, 11); N_sites, S_sites (8, 9) do not come out of the
    # iteration and stay
    def assert_no_estimate(self, row, method_name):
        self.assertEqual(row[13], "not_converged", row)
        self.assertEqual([row[i] for i in (3, 4, 5, 6, 7, 10, 11)], ["NA"] * 7, row)
        self.assertTrue(float(row[8]) > 0 and float(row[9]) > 0, row)
        rounds = {"YN": "YN00 iteration did not converge in 10 rounds",
                  "MYN": "MYN iteration did not converge in 100 rounds"}
        self.assertRegex(row[14], "^" + rounds[method_name] + r": no estimate \(last round: Ka \S+, Ks \S+; "
                                  r"such values can differ between platforms\)$")

    def test_iteration_that_does_not_converge_gives_no_estimate(self):
        """An iteration that ends at its round limit has no fixed point to report. The values of its last round
        depend on the platform: for yn_not_converged (one codon, Nd/N on the 0.75 limit of JC69) Ka of the last
        round was 3.16884648849903 on Linux (GCC 15, glibc 2.43), 3.17357163988094 on Windows (MSVC),
        3.1735716240756 on macOS (Apple Clang, arm64) and "saturated" with GCC 13 and glibc 2.39 (CI, 2026-11-13).
        So this test asserts what holds everywhere: the status, the NA fields and the form of the note. No digit of
        a last round is asserted. (KaKs_Calculator prints its last round: for myn_not_converged Ka
        1.74822892588737, Ks 2.17994969816181.)"""
        self.assert_no_estimate(self.rows[("myn_not_converged", "MYN")], "MYN")
        self.assertEqual(self.rows[("myn_not_converged", "YN")][13], "ok")
        self.assert_no_estimate(self.rows[("yn_not_converged", "YN")], "YN")
        # not_converged also where the last round was saturated: whether it was differs by platform as well
        # (sat_not_converged: Ka 0.938 and Ks 1.570 with one C library, Ka 0.574 and Ks saturated with another)
        self.assert_no_estimate(self.rows[("sat_not_converged", "MYN")], "MYN")
        for key, r in self.rows.items():
            if "did not converge" in r[14]:
                self.assertEqual(r[13], "not_converged", key)
            else:
                self.assertNotEqual(r[13], "not_converged", key)
        # the unweighted estimate needs no iteration: it is yn00's (weighting 0: dN 0.636574320355)
        self.assertTrue(close(float(self.paml0[("sat_not_converged", "YN")][3]), 0.636574320355, 1e-11))
        for key, r in self.paml0.items():
            self.assertNotIn("did not converge", r[14], key)
            self.assertNotEqual(r[13], "not_converged", key)

    def test_kaks_format_has_no_last_round_values(self):
        """`coal kaks` has no status column: the row of an iteration that did not converge has NA in the columns
        that come out of the iteration (KaKs_Calculator prints its last round there), the rest of the row is as
        for any pair, and the summary says how many."""
        r = self.run_cpp("kaks", "-i", self.axt, "-m", "MYN")
        self.assertIn("did not converge", r.stderr)
        self.assertIn("3 ", r.stderr.split("did not converge")[0].splitlines()[-1])
        self.assertIn("3 not_converged", r.stderr)
        rows = read_rows(r.stdout)
        for pair in ("myn_not_converged", "yn_not_converged", "sat_not_converged"):
            f = rows[(pair, "MYN")]
            self.assertEqual([f[i] for i in NOT_CONVERGED_NA], ["NA"] * len(NOT_CONVERGED_NA), f)
            self.assertTrue(float(f[7]) > 0 and float(f[8]) > 0 and int(f[6]) > 0 and float(f[10]) > 0, f)
        self.assertNotEqual(rows[("myn_nan_round", "MYN")][2], "NA")

    def test_fisher_test_with_a_cell_that_is_zero_up_to_rounding(self):
        """Weighted counts that add up to the number of sites (Sd = S = 2 here, or 2.0000000000000004): the cell
        S - Sd is zero, and the test is made (p = 0 in KaKs_Calculator's version, which skips tables with an empty
        cell) whatever the last bit says. It used to be made or left out by platform (0 or NA)."""
        with tempfile.TemporaryDirectory() as td:
            axt = Path(td) / "t.axt"
            axt.write_text("p682\nGCATGGAAG\nTGGGCAAAG\n\n", encoding="utf-8")
            f = read_rows(self.run_cpp("kaks", "-i", axt, "-q", "-m", "MYN", "--digits", "15").stdout)[("p682", "MYN")]
        self.assertEqual((f[3], f[5]), ("NA", "0"), f)                  # Ks saturated; p of the degenerate table
        self.assertTrue(close(float(f[11]), float(f[7]), 1e-12), f)     # Sd = S


class TestAgainstPAML(KaKsBase):
    def test_yn00(self):
        for w in (1, 0):
            ours = self.dnds(KDIR / "sim_small.axt", "-m", "YN", "--yn-flavour", "paml", "--yn-weighting", str(w))
            lines = (GOLD / f"sim_small.yn00_w{w}.tsv").read_text(encoding="utf-8").splitlines()[1:]
            self.assertEqual(len(lines), 48)
            for line in lines:
                pair, dN, dS, S, N, SEdN, SEdS = line.split("\t")
                o = ours[(pair, "YN")]
                if o[13] == "not_converged":       # no estimate where yn00 prints its tenth round
                    self.assertEqual((w, pair, o[3], o[4]), (1, "sim_w0.1_t2_k2_3&b", "NA", "NA"), o)
                    continue
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
                         (GOLD / f"{name}.coal.kaks").read_text(encoding="utf-8"), f"{name} kaks")
            self.compare(self.run_cpp("dnds", axt, "-q", *flags).stdout,
                         (GOLD / f"{name}.coal.dnds.tsv").read_text(encoding="utf-8"), f"{name} dnds")


class TestPolicy(KaKsBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.rows = {}
        r = subprocess.run([str(cls.bin), "dnds", str(KDIR / "edge_cases.axt"), "-q",
                            *sum((["-m", m] for m in METHODS), [])], capture_output=True, text=True, check=True, encoding="utf-8")
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
        text = (KDIR / "edge_cases.axt").read_text(encoding="utf-8").split("\n\n")
        rec = [r for r in text if r.startswith("lowercase_rna")][0].splitlines()
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "u.axt"
            p.write_text(f"upper\n{rec[1].upper().replace('U', 'T')}\n{rec[2].upper()}\n\n", encoding="utf-8")
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
            rows = out.read_text(encoding="utf-8").splitlines()[1:]
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
        recs = (KDIR / "sim_small.axt").read_text(encoding="utf-8").split("\n\n")[:3]
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "a.fasta"
            with open(fa, "w", encoding="utf-8") as f:
                for i, r in enumerate(recs):
                    l = r.strip().splitlines()
                    f.write(f">s{i}a\n{l[1]}\n>s{i}b\n{l[2]}\n")
            allp = self.run_cpp("dnds", fa, "-q").stdout.splitlines()[1:]
            self.assertEqual(len(allp), 15)  # 6 sequences: every pair
            cons = self.run_cpp("dnds", fa, "-q", "--pairs", "consecutive").stdout.splitlines()[1:]
            self.assertEqual([l.split("\t")[:2] for l in cons], [["s0a", "s0b"], ["s1a", "s1b"], ["s2a", "s2b"]])
            # consecutive FASTA pairs = the AXT records
            axt = Path(td) / "a.axt"
            axt.write_text("\n\n".join(recs) + "\n\n", encoding="utf-8")
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
            p.write_text(f"long\n{s1}\n{s2}\n\n", encoding="utf-8")
            r = self.run_cpp("dnds", p, "-q", *sum((["-m", m] for m in METHODS), []), timeout=120)
            rows = [l.split("\t") for l in r.stdout.splitlines()[1:]]
            self.assertEqual(len(rows), 7)
            for f in rows:
                self.assertIn(f[13], ("ok", "warning"), f)  # (random codons may form stop codons in s2)
                self.assertGreater(float(f[3]), 0)
                self.assertGreater(float(f[4]), 0)

    def test_output_onto_the_input_is_refused(self):
        """The output is opened before the input is read to its end: -o naming the input would empty it."""
        text = "a&b\nATGAAACCC\nATGAAGCCC\n\n"
        with tempfile.TemporaryDirectory() as td:
            axt = Path(td) / "in.axt"
            axt.write_text(text, encoding="utf-8")
            spellings = [axt, Path(td) / "." / "in.axt"]
            if os.name != "nt":
                os.symlink(axt, Path(td) / "link.axt")
                os.link(axt, Path(td) / "hard.axt")
                spellings += [Path(td) / "link.axt", Path(td) / "hard.axt"]
            for out in spellings:
                for cmd in (["kaks", "-i", axt, "-o", out], ["dnds", axt, "-o", out], ["dnds", out, "-o", axt]):
                    r = self.run_cpp(*cmd, check=False)
                    self.assertEqual(r.returncode, 2, (cmd, r.stderr))
                    self.assertIn("is the input file", r.stderr)
                    self.assertEqual(axt.read_text(encoding="utf-8"), text, f"{cmd}: the input was changed")
            fa = Path(td) / "cds.fa"
            fa.write_text(">g\nATGAAA\n", encoding="utf-8")
            for cmd in (["translate", fa, "-o", fa], ["translate", fa, "--genes", fa],
                        ["codon-dnds", "-i", fa, "-o", fa]):
                r = self.run_cpp(*cmd, check=False)
                self.assertEqual(r.returncode, 2, (cmd, r.stderr))
                self.assertIn("is the input file", r.stderr)
                self.assertEqual(fa.read_text(encoding="utf-8"), ">g\nATGAAA\n")
            r = self.run_cpp("kaks", "-i", axt, "-o", Path(td) / "out.kaks", "-q")  # another file: fine
            self.assertEqual((Path(td) / "out.kaks").read_text(encoding="utf-8").count("\n"), 2)

    def test_axt_record_with_one_sequence_line(self):
        """A name and one sequence line is an error for that record; it is never cut in half and compared."""
        with tempfile.TemporaryDirectory() as td:
            axt = Path(td) / "in.axt"
            axt.write_text("a&b\nATGAAACCC\nATGAAGCCC\n\none\nATGAAACCCATGAAGCCC\n\n"
                           "wrapped\nATGAAA\nCCC\nATGAAG\nCCC\n\nc&d\nATGAAA\nATGAAG\n", encoding="utf-8")
            r = self.run_cpp("dnds", axt)
            rows = {f[0]: f for f in (l.split("\t") for l in r.stdout.splitlines()[1:])}
            self.assertEqual(list(rows), ["a", "one", "wrapped", "c"])
            self.assertEqual(rows["one"][13], "error")
            self.assertIn("one sequence line", rows["one"][14])
            self.assertEqual(rows["one"][3:6], ["NA", "NA", "NA"])
            self.assertEqual(rows["a"][3:13], rows["wrapped"][3:13])   # 4 lines: joined and split in half
            self.assertEqual((rows["a"][13], rows["c"][13]), ("ok", "ok"))
            self.assertIn("1 error", r.stderr)
            k = self.run_cpp("kaks", "-i", axt)   # KaKs_Calculator's format: the record is left out, with a warning
            self.assertEqual([l.split("\t")[0] for l in k.stdout.splitlines()[1:]], ["a&b", "wrapped", "c&d"])
            self.assertIn("one: record has one sequence line", k.stderr)

    def test_python_hands_over(self):
        env = dict(os.environ, COAL_CPP=str(self.bin), PYTHONPATH=str(ROOT / "src"))
        r = subprocess.run([sys.executable, "-m", "coal", "dnds", str(KDIR / "edge_cases.axt"), "-q"],
                           capture_output=True, text=True, env=env, encoding="utf-8")
        self.assertEqual(r.returncode, 0, r.stderr)
        direct = self.run_cpp("dnds", KDIR / "edge_cases.axt", "-q").stdout
        self.assertEqual(r.stdout, direct)
        env["COAL_CPP"] = str(Path(tempfile.gettempdir()) / "no-such-coal-cpp")
        r = subprocess.run([sys.executable, "-m", "coal", "kaks", "-i", "x"], capture_output=True, text=True, env=env, encoding="utf-8")
        self.assertEqual(r.returncode, 1)
        self.assertIn("needs coal-cpp", r.stderr)


if __name__ == "__main__":
    unittest.main()
