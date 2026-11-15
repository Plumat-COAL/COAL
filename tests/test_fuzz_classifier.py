"""The classifier of tests/fuzz_kaks_vs_reference.py (classify_kk) on one fixture per documented difference.

The fuzz against KaKs_Calculator runs by hand only (it needs 15-digit builds of the references), so a classifier
that no longer knows a documented difference would show up as "unexplained" weeks later, or, worse, a class that
is too wide would hide a regression. This test needs no binary: each fixture is a KaKs_Calculator row, COAL's row
and COAL's status, and the field must land in the class of its difference. The numbered list of docs/kaks.md
("Differences from KaKs_Calculator") is read, so a new difference fails here until it has a fixture or is listed
in NOT_IN_FUZZ with the reason the fuzz cannot meet it. The negative fixtures check that the classes are not
wider than their difference.
"""
from __future__ import annotations

import re
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fuzz_kaks_vs_reference as F  # noqa: E402
import test_kaks_properties as P  # noqa: E402

DOC = HERE.parent / "docs" / "kaks.md"
SAT_KS = "Ks saturated (multiple-hit correction undefined)"
SAT_KA = "Ka saturated (multiple-hit correction undefined)"
COL = {"Ka": 2, "Ks": 3, "Ka/Ks": 4, "P-Value": 5, "S-Sites": 7, "N-Sites": 8, "Fold-Sites": 9,
       "S-Substitutions": 11, "N-Substitutions": 12, "Rate-Ratio": 16}

# differences the fuzz cannot meet, and why
NOT_IN_FUZZ = {
    2: "Fisher's p (P-Value) is not compared by the fuzz",
    3: "Fisher's p (P-Value) is not compared by the fuzz",
    6: "reference_cases() leaves out RNA and lower-case pairs",
    7: "reference_cases() leaves out pairs whose length is not a multiple of 3 (and unequal ones)",
    17: "YN / MYN's S- and N-Substitutions (Sd, Nd) are not compared by the fuzz (only the LWL family's)",
}

# the MLWL pair of tests/kaks/mlwl_kappa_boundary.axt (code 1, ts = tv = 3.5, kappa exactly 2): COAL takes the
# formula of kappa = 2 ("not above 2"), KaKs_Calculator that of kappa > 2
K2_FOLD = "11:4:3"
K2_OURS = {"Ka": "0.457209261327193", "Ks": "1.04591163912279", "Ka/Ks": "0.437139471658098",
           "S-Sites": "4.71428571428571", "N-Sites": "13.2857142857143", "Fold-Sites": K2_FOLD,
           "Rate-Ratio": "2:2:1:1:1:1"}
K2_REF = {"Ka": "0.446195", "Ks": "1.18", "Ka/Ks": "0.378", "S-Sites": "4.33333333333333",
          "N-Sites": "13.6666666666667", "Fold-Sites": K2_FOLD, "Rate-Ratio": "2:2:1:1:1:1"}
# docs/kaks.md difference 16: code 2 (ATA is Met), ATA -> ATG is synonymous; KaKs_Calculator's MLWL counts it
# as Tzeng et al.'s nonsynonymous Ile <-> Met
CODE2 = ("ATAGCTAAACTGGGT", "ATGGCCAAACTGGGT")


def row(name, method, **fields):
    r = ["NA"] * 22
    r[0], r[1] = name, method
    for k, v in fields.items():
        r[COL[k]] = v
    return r


def fx(method, field, ref, ours, status="ok", note="", case=None):
    """(method, field, KaKs_Calculator's fields, COAL's fields, COAL's status and note, Case or None)"""
    return method, field, ref, ours, status, note, case


# difference -> fixtures (each: the field and the substring its class must contain; None: not counted at all)
FIXTURES = {
    1: [(fx("YN", "Ks", {"Ks": "3.69"}, {}, "saturated", SAT_KS), "difference 1")],
    4: [(fx("NG", "Ka/Ks", {"Ka/Ks": "-0", "Ka": "0.1", "Ks": "0.2"}, {"Ka/Ks": "0", "Ka": "0.1", "Ks": "0.2"}),
         None)],
    5: [(fx("LPB", "S-Sites", {"S-Sites": "-nan"}, {}), "difference 5")],
    8: [(fx("YN", "Ka", {"Ka": "NA"}, {"Ka": "1.64791843300216"}, "no_data", "no synonymous sites"),
         "difference 8")],
    9: [(fx("MYN", "Ka", {"Ka": "nan", "Ks": "nan"}, {"Ka": "0.39", "Ks": "0.25"}), "differences 8-10")],
    10: [(fx("MYN", "Ks", {"Ka": "-nan", "Ks": "-nan"}, {"Ka": "0.39", "Ks": "0.25"}), "differences 8-10")],
    11: [(fx("YN", "Ka", {"Ka": "0.5", "Ks": "NA"}, {"Ka": "0.6"}, "no_data", "no synonymous sites"),
          "difference 11"),
         (fx("MYN", "Ka", {"Ka": "0.166967", "Ks": "nan"}, {"Ka": "1.017615", "Ks": "0.4"}), "difference 11")],
    12: [(fx("YN", "Ka", {"Ka": "4.49007072568802"}, {}, "not_converged",
             "YN00 iteration did not converge in 10 rounds: no estimate"), "difference 12")],
    13: [(fx("MYN", "Ks", {"Ks": "1.89884662745176"}, {"Ks": "1.89884663167766"}, "saturated", SAT_KA),
          "difference 13")],
    # a rounding residue printed as an estimate by KaKs_Calculator, saturated in COAL: such rows are difference 1
    14: [(fx("NG", "Ks", {"Ks": "27.5526004272578"}, {}, "saturated", SAT_KS), "difference 1")],
    # MYN with an infinite kappa: KaKs_Calculator prints nan for the estimates and the site counts
    15: [(fx("MYN", "S-Sites", {"Ka": "nan", "Ks": "nan", "S-Sites": "nan"}, {}, "saturated",
             "MYN: kappa_AG estimated as infinite"), "difference 5")],
    16: [(fx("MLWL", "Ka", {"Ka": "0.154", "N-Substitutions": "1.5"}, {"Ks": "0.1"}, "ok", "",
             P.Case("d16", *CODE2, code=2)), "difference 16"),
         (fx("MLPB", "N-Substitutions", {"Ka": "0.137", "N-Substitutions": "1.5"}, {"Ks": "0.1"}, "ok", "",
             P.Case("d16", *CODE2, code=2)), "difference 16"),
         (fx("MLWL", "Ks", {"Ks": "0.30"}, {"Ks": "0.25"}, "ok", "", P.Case("d16", *CODE2, code=2)),
          "difference 16")],
    18: [(fx("MLWL", "S-Sites", K2_REF, K2_OURS), "difference 18"),
         (fx("MLWL", "Ks", K2_REF, K2_OURS), "difference 18")],
    19: [(fx("LPB", "S-Sites", {"S-Sites": "6", "Ks": "0.366"}, {}, "saturated", SAT_KS), "difference 19"),
         (fx("MLPB", "N-Sites", {"N-Sites": "333", "Ka": "0.131"}, {"Ks": "0.2"}, "saturated", SAT_KA),
          "difference 19"),
         (fx("LWL", "S-Substitutions", {"S-Substitutions": "2.19936276082409"}, {}, "saturated", SAT_KS),
          "difference 19"),
         (fx("MLWL", "N-Substitutions", {"N-Substitutions": "43.75"}, {}, "saturated", f"{SAT_KS}; {SAT_KA}"),
          "difference 19")],
}

# fixtures that must stay unexplained: the classes are no wider than their difference
UNEXPLAINED = [
    # LWL / MLWL's site counts are fold-class counts, never NA beside a saturated correction
    fx("LWL", "S-Sites", {"S-Sites": "5.1667"}, {}, "saturated", SAT_KS),
    # LPB's S-Sites (Sd / Ks) NA although only Ka is saturated and Ks is a number
    fx("LPB", "S-Sites", {"S-Sites": "11.5", "Ks": "0.3"}, {"Ks": "0.3"}, "saturated", SAT_KA),
    # Sd NA on the side that is not saturated
    fx("LWL", "S-Substitutions", {"S-Substitutions": "3.3"}, {"Ks": "0.3"}, "saturated", SAT_KA),
    # the pair of difference 16 with the standard code: the rules are the same, so a difference is a bug
    fx("MLWL", "Ka", {"Ka": "0.154"}, {"Ka": "0.140"}, "ok", "", P.Case("d16", *CODE2, code=1)),
    # an LWL-family difference in a pair that meets no special case of its code
    fx("LWL", "Ka", {"Ka": "0.154"}, {"Ka": "0.140"}, "ok", "", P.Case("d16", "GCTGCT", "GCCGCA", code=2)),
    # MLWL at kappa 2 whose sites are not the other formula
    fx("MLWL", "Ks", dict(K2_REF, **{"S-Sites": K2_OURS["S-Sites"], "N-Sites": K2_OURS["N-Sites"]}), K2_OURS),
]


def classify(f):
    method, field, ref, ours, status, note, case = f
    name = case.name if case else "p"
    rf, of = row(name, method, **ref), row(name, method, **ours)
    counts, _, bad = F.classify_kk({(name, method): rf}, {(name, method): of},
                                   {(name, method): {"status": status, "note": note}}, 1e-9,
                                   {name: case} if case else {})
    kinds = [k[2] for k in counts if k[1] == field]
    return kinds, bad


def documented_differences() -> set[int]:
    text = DOC.read_text(encoding="utf-8")
    sec = text.split("### Differences from KaKs_Calculator", 1)[1].split("### Platform reproducibility", 1)[0]
    return {int(m) for m in re.findall(r"^(\d+)\. \*\*", sec, re.M)}


class TestFuzzClassifier(unittest.TestCase):
    def test_every_documented_difference_has_a_fixture(self):
        docs = documented_differences()
        self.assertGreaterEqual(len(docs), 19)
        self.assertEqual(docs, set(FIXTURES) | set(NOT_IN_FUZZ),
                         "docs/kaks.md's differences and this test's FIXTURES / NOT_IN_FUZZ must name the same ones")
        self.assertFalse(set(FIXTURES) & set(NOT_IN_FUZZ))

    def test_fixtures_land_in_their_class(self):
        for d, items in FIXTURES.items():
            for f, want in items:
                with self.subTest(difference=d, method=f[0], field=f[1]):
                    kinds, bad = classify(f)
                    if want is None:
                        self.assertEqual(kinds, [])
                    else:
                        self.assertEqual(len(kinds), 1, kinds)
                        self.assertIn(want, kinds[0])
                        self.assertNotIn("unexplained", kinds[0])
                    self.assertFalse([b for b in bad if b[1] == f[1]], bad)

    def test_classes_are_not_wider_than_their_difference(self):
        for f in UNEXPLAINED:
            with self.subTest(method=f[0], field=f[1], status=f[4]):
                kinds, bad = classify(f)
                self.assertEqual(len(kinds), 1, kinds)
                self.assertIn("unexplained", kinds[0])
                self.assertTrue(any(b[1] == f[1] for b in bad))

    def test_substitution_counts_compared_for_the_lwl_family_only(self):
        kinds, bad = classify(fx("YN", "S-Substitutions", {"S-Substitutions": "69.3"}, {"S-Substitutions": "70"}))
        self.assertEqual((kinds, bad), ([], []))

    def test_special_case_detector(self):
        a, b = CODE2
        self.assertTrue(F.lwl_code_special_case(P.Case("x", a, b, code=2), "MLWL"))
        self.assertFalse(F.lwl_code_special_case(P.Case("x", a, b, code=2), "LWL"))   # LWL has no Ile rule
        self.assertFalse(F.lwl_code_special_case(P.Case("x", a, b, code=1), "MLWL"))
        # code 13: AGA is Gly, KaKs_Calculator's LWL still splits CGA <-> AGA into half a transition
        self.assertTrue(F.lwl_code_special_case(P.Case("x", "CGAGCT", "AGAGCC", code=13), "LWL"))
        self.assertFalse(F.lwl_code_special_case(P.Case("x", "CGAGCT", "AGAGCC", code=1), "LWL"))
        self.assertFalse(F.lwl_code_special_case(P.Case("x", "CGAGCT", "AGAGCC", code=13), "YN"))
        # the two step rules agree on every step of the standard code
        for x in P.SENSE[1]:
            for pos in range(3):
                for base in "TCAG":
                    y = x[:pos] + base + x[pos + 1:]
                    if y == x or P.AA[1][y] == "*":
                        continue
                    for mod in (False, True):
                        self.assertEqual(F._kk_step_lwl(1, mod, x, y, pos), F._coal_step_lwl(1, mod, x, y, pos),
                                         (x, y, pos, mod))

    def test_mlwl_boundary_is_the_benchmarks_rule(self):
        sys.path.insert(0, str(HERE.parent / "benchmarks" / "kaks"))
        import random_short_pairs as R
        rf, of = row("p", "MLWL", **K2_REF), row("p", "MLWL", **K2_OURS)
        self.assertTrue(F.mlwl_boundary(rf, of))
        self.assertEqual(F.mlwl_boundary(rf, of), R.mlwl_boundary(rf, of))
        of2 = row("p", "MLWL", **dict(K2_OURS, **{"Rate-Ratio": "2.1:2.1:1:1:1:1"}))
        self.assertFalse(F.mlwl_boundary(rf, of2))
        self.assertEqual(F.mlwl_boundary(rf, of2), R.mlwl_boundary(rf, of2))


if __name__ == "__main__":
    unittest.main()
