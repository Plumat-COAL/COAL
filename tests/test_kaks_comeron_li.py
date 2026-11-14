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
import os
import sys
import unittest

# (importable both as `discover -s tests` / `test_kaks_comeron_li` and as `tests.test_kaks_comeron_li`)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_kaks as tk   # noqa: E402  (the module, not its classes: imported test classes would run again here)

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
        # the policy of the default flavour, with the pathways and the NA rule of the gestimator flavour
        # (--comeron-twosubs libsequence --comeron-min-sites 0: COAL's defaults before D8 / D9, 2026-11-14)
        for s in SETS:
            seqs = read_axt(KDIR / f"{s}.axt")
            dfl = self.tidy(KDIR / f"{s}.axt", "-m", "Comeron", "--comeron-twosubs", "libsequence",
                            "--comeron-min-sites", "0")
            leg = self.tidy(KDIR / f"{s}.axt", "-m", "Comeron", "--comeron-flavour", "gestimator",
                            "--comeron-min-sites", "0")
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
            # the two options restore the output of COAL before D8 / D9 byte for byte
            got = self.run_cpp("dnds", axt, "--format", "axt", "-q", "--digits", "15", "-m", "Comeron", "-m", "Li",
                               "--comeron-twosubs", "libsequence", "--comeron-min-sites", "0").stdout
            self.assertEqual(got, (GOLD / f"{s}.coal.comeron_li.libsequence.tsv").read_text(encoding="utf-8"), s)
            cmp(self, self.run_cpp("kaks", "-i", axt, "-q", "-m", "Comeron", "-m", "Li").stdout,
                (GOLD / f"{s}.coal.comeron_li.kaks").read_text(encoding="utf-8"), f"{s} kaks")


def random_pairs(n: int, codons: int, seed: int):
    """Random sense-codon pairs in which many codons differ at two positions (the TwoSubs case)."""
    import random
    rng = random.Random(seed)
    bases = "TCAG"
    sense = [a + b + c for a in bases for b in bases for c in bases if a + b + c not in STOPS]
    out = []
    for k in range(n):
        x, y = [], []
        for _ in range(codons):
            c = rng.choice(sense)
            x.append(c)
            r = rng.random()
            if r < 0.6:
                y.append(c)
                continue
            npos = 2 if r < 0.85 else 1
            while True:
                d = list(c)
                for pos in rng.sample(range(3), npos):
                    d[pos] = rng.choice([b for b in bases if b != c[pos]])
                if "".join(d) not in STOPS:
                    break
            y.append("".join(d))
        out.append((f"p{k}", "".join(x), "".join(y)))
    return out


class TestTwoSubsAndMinSites(ComeronLiBase):
    """D8 (2026-11-14): Comeron's two-position pathways (--comeron-twosubs); D9: the minimum number of sites for
    which a class's undefined correction makes dN or dS NA (--comeron-min-sites). docs/kaks.md."""

    def write_pairs(self, tmp, pairs, swap=False):
        from pathlib import Path
        p = Path(tmp) / ("ba.axt" if swap else "ab.axt")
        with open(p, "w", encoding="utf-8", newline="\n") as f:
            for name, a, b in pairs:
                f.write(f"{name}\n{b if swap else a}\n{a if swap else b}\n\n")
        return p

    def test_symmetric_by_default(self):
        # a distance does not depend on which sequence is first: with Comeron's pathways swapping the rows changes
        # no digit (0 of 300 pairs); libsequence's branches change Ka in most pairs (279 of these 300)
        import tempfile
        pairs = random_pairs(300, 150, 20261114)
        with tempfile.TemporaryDirectory() as tmp:
            ab, ba = self.write_pairs(tmp, pairs), self.write_pairs(tmp, pairs, swap=True)
            for flav in ((), ("--comeron-flavour", "gestimator", "--comeron-twosubs", "comeron")):
                x, y = self.tidy(ab, "-m", "Comeron", *flav), self.tidy(ba, "-m", "Comeron", *flav)
                for key in x:
                    for k in ("dN", "dS", "dNdS", "status"):
                        self.assertEqual(x[key][k], y[key][k], f"{flav} {key} {k}")
            x = self.tidy(ab, "-m", "Comeron", "--comeron-twosubs", "libsequence")
            y = self.tidy(ba, "-m", "Comeron", "--comeron-twosubs", "libsequence")
            self.assertGreater(sum(x[k]["dN"] != y[k]["dN"] for k in x), 200)
            # the default is comeron in the coal flavour, libsequence in the gestimator flavour
            dfl = self.tidy(ab, "-m", "Comeron")
            self.assertEqual(dfl, self.tidy(ab, "-m", "Comeron", "--comeron-twosubs", "comeron"))
            g = self.tidy(ab, "-m", "Comeron", "--comeron-flavour", "gestimator")
            self.assertEqual(g, self.tidy(ab, "-m", "Comeron", "--comeron-flavour", "gestimator",
                                          "--comeron-twosubs", "libsequence"))
            self.assertNotEqual([dfl[k]["dN"] for k in dfl], [g[k]["dN"] for k in g])

    def test_one_position_pairs_unchanged(self):
        # only codons that differ at exactly two positions take the TwoSubs pathways
        import tempfile
        a = "GCTAAAGGTCTGTTTGATCATACC" * 5
        b = "GCAAAGGGTTTGTTCGACCACACT" * 5      # every codon differs at one position at most
        with tempfile.TemporaryDirectory() as tmp:
            p = self.write_pairs(tmp, [("x", a, b)])
            self.assertEqual(self.tidy(p, "-m", "Comeron", "--comeron-twosubs", "comeron"),
                             self.tidy(p, "-m", "Comeron", "--comeron-twosubs", "libsequence"))

    def test_min_sites_never_changes_a_value(self):
        # a value that is not NA under the previous rule (--comeron-min-sites 0) is the same under the default (3)
        # and under a larger threshold; only NA rows can get a value, with status warning and the note
        for s in SETS:
            old = self.tidy(KDIR / f"{s}.axt", "-m", "Comeron", "--comeron-min-sites", "0")
            for ms in ("3", "10"):
                new = self.tidy(KDIR / f"{s}.axt", "-m", "Comeron", "--comeron-min-sites", ms)
                if ms == "3":
                    self.assertEqual(new, self.tidy(KDIR / f"{s}.axt", "-m", "Comeron"))
                for key, o in old.items():
                    n = new[key]
                    for k, which in (("dN", "Ka"), ("dS", "Ks")):
                        if o[k] != "NA":
                            self.assertEqual(n[k], o[k], f"{s} {key} {k}")
                        elif n[k] != "NA":
                            self.assertIn(f"{which}: a site class with fewer than {ms} sites", n["note"], key)
                            self.assertIn(n["status"], ("warning", "saturated", "no_data"), key)

    def test_min_sites_one_codon(self):
        # one codon, nonsynonymous change: the 2-fold transversional class has < 1 site, its correction is
        # undefined. Previous rule: dN NA (saturated); default: gestimator's value, flagged
        d = self.pair_rows("GCT", "GAT")
        self.assertEqual(d["0"]["dN"], "NA")
        self.assertEqual(d["0"]["status"], "saturated")
        self.assertNotEqual(d["3"]["dN"], "NA")
        self.assertEqual(d["3"]["status"], "warning")
        self.assertIn("fewer than 3 sites has an undefined correction", d["3"]["note"])
        g = self.pair_rows("GCT", "GAT", "--comeron-flavour", "gestimator")
        self.assertEqual(g["3"]["dN"], d["3"]["dN"])      # gestimator's number
        self.assertEqual(g["0"]["dN"], d["3"]["dN"])      # the gestimator flavour always gives it
        # only the note says why: orthologr reads "Ka: correction undefined" as "not an estimate"
        self.assertNotIn("Ka: correction undefined", g["3"]["note"])
        self.assertIn("Ka: correction undefined", g["0"]["note"])

    def pair_rows(self, a, b, *extra):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            p = self.write_pairs(tmp, [("x", a, b)])
            return {ms: self.tidy(p, "-m", "Comeron", "--comeron-min-sites", ms, *extra)[("x", "Comeron")]
                    for ms in ("0", "3")}

    def test_options_refused(self):
        r = self.run_cpp("dnds", KDIR / "sim_small.axt", "-m", "Comeron", "--comeron-twosubs", "x", check=False)
        self.assertEqual(r.returncode, 2)
        self.assertIn("expected comeron or libsequence", r.stderr)
        for bad in ("-1", "x", "nan"):
            r = self.run_cpp("dnds", KDIR / "sim_small.axt", "-m", "Comeron", "--comeron-min-sites", bad, check=False)
            self.assertEqual(r.returncode, 2, bad)
            self.assertIn("--comeron-min-sites", r.stderr)
        r = self.run_cpp("kaks", "-i", KDIR / "sim_small.axt", "-q", "-m", "Comeron", "--comeron-twosubs",
                         "libsequence", "--comeron-min-sites", "0")
        self.assertEqual(r.stdout, self.run_cpp("kaks", "-i", KDIR / "sim_small.axt", "-q", "-m", "Comeron",
                                                "--comeron-twosubs=libsequence", "--comeron-min-sites=0").stdout)


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

    def test_not_converged_is_a_status_of_the_iterating_methods_only(self):
        """Status `not_converged` (YN, MYN: no estimate) next to Comeron and Li, which do not iterate: the summary
        line names all six statuses in the engine's order whatever the methods, and no Comeron or Li row has it."""
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "p.axt"
            p.write_text("nc\nAAA\nATG\n\n", encoding="utf-8")
            r = self.run_cpp("dnds", p, "--format", "axt", "-m", "YN", "-m", "Comeron", "-m", "Li")
            rows = {f.split("\t")[2]: f.split("\t") for f in r.stdout.splitlines()[1:]}
            self.assertEqual(rows["YN"][13], "not_converged")
            for m in ("Comeron", "Li"):
                self.assertNotEqual(rows[m][13], "not_converged", m)
                self.assertIn(rows[m][13], ("ok", "warning", "saturated", "no_data"), m)
            counts = r.stderr.splitlines()[0]
            self.assertRegex(counts, r"3 method\(s\): \d+ ok, \d+ warning, 1 not_converged, \d+ saturated, "
                                     r"\d+ no_data, 0 error$")
        for flav in ((), ("--comeron-flavour", "gestimator", "--li-flavour", "seqinr")):
            for name in SETS:
                for (_, m), d in self.tidy(KDIR / (name + ".axt"), "-m", "Comeron", "-m", "Li", *flav).items():
                    self.assertNotEqual(d["status"], "not_converged", (name, m))
                    self.assertIn(d["status"], ("ok", "warning", "saturated", "no_data", "error"), (name, m))

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


class TestComeronLiPipeline(ComeronLiBase):
    """`coal pipeline` and coal.dnds_pipeline pass the two methods and their flavours to the engine."""

    def test_method_names(self):
        from coal import dnds_pipeline as dp
        self.assertEqual(dp.METHODS[-2:], ("Comeron", "Li"))
        for given, want in (("comeron", "Comeron"), ("C95", "Comeron"), ("c95", "Comeron"), ("LI", "Li"),
                            ("yn", "YN"), ("nope", "nope")):
            self.assertEqual(dp._method_name(given), want)
        a = dp._namespace("pairs", ["--query", "q", "--subject", "s", "-o", "o"], {"method": ["c95", "li", "yn"]})
        self.assertEqual(a.method, ["Comeron", "Li", "YN"])
        self.assertEqual((a.comeron_flavour, a.li_flavour), ("coal", "coal"))
        with self.assertRaisesRegex(ValueError, "standard genetic code only"):
            dp._namespace("pairs", ["--query", "q", "--subject", "s", "-o", "o"],
                          {"method": ["Comeron"], "codontable": 2})
        with self.assertRaises(ValueError):
            dp._namespace("pairs", ["--query", "q", "--subject", "s", "-o", "o"], {"comeron_flavour": "x"})

    def test_codon_dnds(self):
        import tempfile
        from pathlib import Path
        from coal import dnds_pipeline as dp
        a = "ATGGCTAAAGGTCTGTTTGATCCGGAATGGTAA"
        b = "ATGGCAAAGGGTTTGTTCGACCCGGAGTGGTGA"
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "pairs.tsv"
            dp.codon_dnds([("g1", [("x", a), ("y", b)])], str(out), method=["Comeron", "Li"], cpp=str(self.bin))
            rows = [l.split("\t") for l in out.read_text(encoding="utf-8").splitlines()]
            self.assertEqual([r[2] for r in rows[1:]], ["Comeron", "Li"])
            default = rows[1][3:5]
            axt = Path(tmp) / "p.axt"

            def direct(x, y, *flags):
                axt.write_text(f"x&y\n{x}\n{y}\n\n", encoding="utf-8")
                lines = self.run_cpp("dnds", axt, "-q", *flags).stdout.splitlines()
                return [l.split("\t")[3:5] for l in lines[1:]]

            def piped(x, y, method, *flags):
                dp.codon_dnds([("g1", [("x", x), ("y", y)])], str(out), method=[method], cpp=str(self.bin),
                              extra_args=list(flags))
                rows = [l.split("\t") for l in out.read_text(encoding="utf-8").splitlines()]
                self.assertEqual(rows[1][2], method)
                return rows[1][3:5]

            gest = ("--comeron-flavour", "gestimator")
            seqi = ("--li-flavour", "seqinr")
            # the reference flavours end the pair before the terminal stop codons, like everything in COAL ...
            self.assertEqual(piped(a, b, "Comeron", *gest), direct(a[:-3], b[:-3], "-m", "Comeron", *gest)[0])
            self.assertEqual(piped(a, b, "Li", *seqi), direct(a[:-3], b[:-3], "-m", "Li", *seqi)[0])
            # ... unless --terminal-stops compare: TAA against TGA is then the last compared column, which is
            # what `coal dnds` gives for the rows with their stop codons (and what orthologr 0.4.3 computes)
            with_stops = piped(a, b, "Comeron", *gest, "--terminal-stops", "compare")
            self.assertEqual(with_stops, direct(a, b, "-m", "Comeron", *gest)[0])
            self.assertNotEqual(with_stops, piped(a, b, "Comeron", *gest))
            self.assertEqual(piped(a, b, "Li", *seqi, "--terminal-stops", "compare"),
                             direct(a, b, "-m", "Li", *seqi)[0])
            # one gene without a terminal stop: no column to restore
            self.assertEqual(piped(a[:-3], b, "Comeron", *gest, "--terminal-stops", "compare"),
                             direct(a[:-3], b[:-3], "-m", "Comeron", *gest)[0])
            # the default flavours never compare stop codons, whatever the option says
            self.assertEqual(piped(a, b, "Comeron", "--terminal-stops", "compare"), default)
            self.assertEqual(piped(a, b, "YN", "--terminal-stops", "compare"), piped(a, b, "YN"))

    def test_terminal_stops_option(self):
        from coal import dnds_pipeline as dp
        a = dp._namespace("pairs", ["--query", "q", "--subject", "s", "-o", "o"], {"terminal_stops": "compare"})
        self.assertEqual(a.terminal_stops, "compare")
        with self.assertRaises(ValueError):
            dp._namespace("pairs", ["--query", "q", "--subject", "s", "-o", "o"], {"terminal_stops": "x"})
        r = self.run_cpp("codon-dnds", "-o", "x.tsv", "--terminal-stops", "maybe", check=False)
        self.assertEqual(r.returncode, 2)
        self.assertIn("--terminal-stops", self.run_cpp("codon-dnds", "-h").stdout)

if __name__ == "__main__":
    unittest.main()
