"""Differential fuzz of the Ka/Ks engine against its references: KaKs_Calculator 2.0 (NG, LWL, MLWL, LPB, MLPB,
YN, MYN, with every genetic code) and PAML yn00 (YN00, weighting 1 and 0). Random pairs from the generator of
tests/test_kaks_properties.py (1 to 2,000 codons, divergence 0 to 100 %, gaps, N / IUPAC, stop codons).

Every field where COAL and the reference print different values is put in a class. The classes of
benchmarks/kaks/random_short_pairs.py are the documented ones for KaKs_Calculator (docs/kaks.md, differences 1
and 5 to 13: nan -> NA, saturated, not converged, sites of one kind only, a NaN in the reference's row, the
saturated kind's other path); for yn00 the same documented policy is applied (saturated: yn00's capped value,
COAL NA; not converged: COAL NA; one kind of site: yn00 nan; dN/dS when dN = 0: yn00's floor 1e-5). Anything
else is "unexplained" and fails the run; the pair is written to the --keep folder.

    python tests/fuzz_kaks_vs_reference.py --kaks KaKs_Calculator_15digits --yn00 yn00_12digits \\
        [--kaks-fixed KaKs_Calculator_corrected] [--n 600] [--seed 1] [--max-codons 2000] [--keep DIR]

As a test (skipped unless KAKS_CALCULATOR and / or YN00 name the binaries; KAKS_CALCULATOR_FIXED optional):
FUZZ_N (default 120), FUZZ_SEED, COAL_SLOW=1 for 1,000 pairs. Both references must print at least 12
significant digits (the builds of docs/kaks.md "Validation"); with the stock 6-digit binaries set
KAKS_REF_DIGITS=6 and the tolerance becomes 1e-5.
"""
from __future__ import annotations

import argparse
import math
import os
import random
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "benchmarks" / "kaks"))
import test_kaks_properties as P  # noqa: E402

KK_METHODS = P.KK_METHODS
# NCBI table -> PAML's icode (PAML's codeml.ctl: 0 universal, 1 vertebrate mt, 2 yeast mt, 3 mold mt,
# 4 invertebrate mt, 5 ciliate, 6 echinoderm mt, 7 euplotid, 8 alternative yeast, 9 ascidian mt, 10 blepharisma)
PAML_ICODE = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 9: 6, 10: 7, 12: 8, 13: 9, 15: 10}
KFIELDS = {2: "Ka", 3: "Ks", 4: "Ka/Ks", 7: "S-Sites", 8: "N-Sites"}
# K03 (a rounding residue printed as an estimate; tests/test_kaks_properties.py) is a known bug: counted in its
# own class, not failed, until it is fixed (then set this to False)
KNOWN_K03_OK = True


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def reference_cases(n: int, seed: int, max_codons: int, codes=None) -> list[P.Case]:
    """Pairs both programs read the same way: no U, no lower case (difference 6), no malformed rows (7)."""
    rng = random.Random(seed)
    out = []
    k = 0
    while len(out) < n:
        code = rng.choice(codes) if codes else None
        c = P.make_case(rng, k, max_codons=max_codons, code=code)
        k += 1
        if c.tags & {"rna", "lower", "unequal", "not_triplet"}:
            continue
        out.append(c)
    return out


# ------------------------------------------------------------------------------------------- KaKs_Calculator

def kk_rows(exe: str, cases, method: str, workdir: Path, per_pair: bool) -> dict:
    """KaKs_Calculator's rows: one process per pair (its MYN carries state from pair to pair, docs/kaks.md)."""
    def one(batch, tag):
        axt, out = workdir / f"{tag}.axt", workdir / f"{tag}.{method}.kaks"
        P.write_axt(batch, axt)
        code = batch[0].code
        try:
            subprocess.run([exe, "-i", str(axt), "-o", str(out), "-m", method, "-c", str(code)],
                           capture_output=True, timeout=600, stdin=subprocess.DEVNULL)
        except subprocess.TimeoutExpired:
            return {(c.name, method): ["TIMEOUT"] * 22 for c in batch}
        text = out.read_text() if out.exists() else ""
        axt.unlink()
        out.unlink(missing_ok=True)
        return P.parse_kaks(text)
    rows = {}
    if per_pair:
        with ThreadPoolExecutor(2) as ex:
            for r in ex.map(lambda c: one([c], c.name), cases):
                rows.update(r)
    else:
        groups = {}
        for c in cases:
            groups.setdefault(c.code, []).append(c)
        for code, g in groups.items():
            rows.update(one(g, f"all{code}"))
    return rows


def coal_rows(cpp, cases, workdir: Path, extra=()) -> tuple[dict, dict]:
    """COAL's kaks rows (15 digits) and its dnds rows (status, note), by genetic code."""
    ours, tidy = {}, {}
    groups = {}
    for c in cases:
        groups.setdefault(c.code, []).append(c)
    for code, g in groups.items():
        axt = workdir / f"coal{code}.axt"
        P.write_axt(g, axt)
        r = P.run([cpp, "kaks", "-i", axt, "-q", "--digits", "15", "-c", code, *P.mflags(KK_METHODS), *extra])
        ours.update(P.parse_kaks(r.out))
        r = P.run([cpp, "dnds", axt, "-q", "--digits", "15", "-c", code, *P.mflags(KK_METHODS), *extra])
        tidy.update(P.parse_dnds(r.out))
    return ours, tidy


def one_sequence_without_synonymous_sites(case) -> bool:
    """Every compared codon of one row has no synonymous neighbour (Met / Trp only): its synonymous-site base
    frequencies are 0/0 in KaKs_Calculator and yn00 (docs/kaks.md difference 10; COAL averages over the other)."""
    if case is None:
        return False
    cols = P.compared_columns(case)
    code = case.code

    def nosyn(c):
        return not any(P.AA[code][x] == P.AA[code][c] for x in P._neighbours(c) if P.AA[code][x] != "*")
    return bool(cols) and (all(nosyn(x) for x, _ in cols) or all(nosyn(y) for _, y in cols))


def classify_kk(ref: dict, ours: dict, tidy: dict, tol: float, cases=None):
    """Field by field, KaKs_Calculator against COAL (the classes of random_short_pairs.classify, generalised to
    any length and to the tolerance of the reference's digits)."""
    counts, examples, bad = Counter(), {}, []
    for key, rf in ref.items():
        of = ours.get(key)
        d = tidy.get(key)
        if of is None or d is None:
            if key[0] in {k[0] for k in tidy} and tidy.get(key, {}).get("status") == "error":
                continue
            bad.append((key, "row", "present", "missing", ""))
            continue
        st, note = d["status"], d.get("note", "")
        if rf[0] == "TIMEOUT":
            counts[(key[1], "*", "reference timeout", st)] += 1
            continue
        ref_nan = any(num(rf[i]) is None and rf[i] != "NA" for i in (2, 3))   # nan / -nan / inf printed
        if key[1] in ("YN", "MYN") and "NA" in (rf[2], rf[3]) and num(of[2]) is not None and num(of[3]) is not None:
            # difference 11: once a NaN is in its iteration, KaKs_Calculator prints a Ka from the remaining
            # codons "with Ks nan or NA beside it" (values above 10, K03, are classed before this)
            ref_nan = True
        if one_sequence_without_synonymous_sites((cases or {}).get(key[0])) and st in P.ESTIMATE:
            ref_nan = True   # the reference's NaN is printed as NA there, and its iteration lost codons (11)
        for i, name in KFIELDS.items():
            x, y = rf[i], of[i]
            if x == y:
                continue
            nx, ny = num(x), num(y)
            if ny is not None and ny > 10 and name in ("Ka", "Ks"):
                kind = "COAL value above 10: rounding residue (K03)"
            elif nx is not None and ny is not None:
                rel = abs(nx - ny) / max(1.0, abs(nx), abs(ny))
                if rel <= tol:
                    continue
                if st == "no_data":
                    kind = "number != number: one kind of site (difference 11)"
                elif st == "saturated" and rel <= 1e-6:
                    kind = "number != number: saturated, other path (difference 13)"
                elif ref_nan:
                    kind = "number != number: reference nan row (difference 11)"
                else:
                    kind = "number != number: unexplained"
            elif nx is not None:
                if st == "not_converged" and name in ("Ka", "Ks", "Ka/Ks"):
                    kind = "number -> NA: not converged (difference 12)"
                elif st == "saturated":
                    kind = "number -> NA: saturated (difference 1)"
                elif abs(nx) < 1e-6 and name in ("Ka", "Ks", "Ka/Ks"):
                    kind = "number -> NA: below 1e-6"
                elif name in ("S-Sites", "N-Sites") and key[1] in ("LPB", "MLPB"):
                    kind = "number -> NA: LPB sites from an NA estimate (difference 5)"
                elif name == "Ka/Ks" and (num(of[3]) is None or num(of[2]) is None):
                    kind = "number -> NA: Ka/Ks of an NA estimate"
                else:
                    kind = "number -> NA: unexplained"
            elif ny is not None:
                if "nan" in x or "inf" in x:
                    kind = "nan -> number (differences 8-10)"
                elif st == "no_data":
                    kind = f"{x} -> number (difference 8)"
                elif ref_nan or (st == "saturated" and key[1] in ("YN", "MYN")):
                    kind = f"{x} -> number: reference NaN (differences 10, 11)"
                elif name in ("Ka", "Ks") and any("nan" in rf[j] for j in range(2, 9)):
                    kind = f"{x} -> number: reference nan row (difference 11)"
                else:
                    kind = f"{x} -> number: unexplained"
            else:
                kind = "nan -> NA (difference 5)"
            k = (key[1], name, kind, st)
            counts[k] += 1
            examples.setdefault(k, f"{key[0]}: {x} / {y}")
            if "unexplained" in kind or "K03" in kind and not KNOWN_K03_OK:
                bad.append((key, name, x, y, st))
    return counts, examples, bad


# ------------------------------------------------------------------------------------------------------ yn00

def yn00_one(exe: str, case: P.Case, weighting: int):
    a, b = [], []
    for x, y in P.compared_columns(case):
        a.append(x)
        b.append(y)
    if not a:
        return None
    with tempfile.TemporaryDirectory() as td:
        Path(td, "in.phy").write_text(f"  2  {3 * len(a)}\ns1  {''.join(a)}\ns2  {''.join(b)}\n")
        Path(td, "yn00.ctl").write_text(
            "seqfile = in.phy\noutfile = out.txt\nverbose = 0\nnoisy = 0\n"
            f"icode = {PAML_ICODE[case.code]}\nweighting = {weighting}\ncommonkappa = 0\ncommonf3x4 = 0\nndata = 1\n")
        try:
            subprocess.run([exe, "yn00.ctl"], cwd=td, capture_output=True, text=True, stdin=subprocess.DEVNULL,
                           timeout=300)
        except subprocess.TimeoutExpired:
            return {"timeout": True}
        rst = Path(td, "rst")
        if not rst.exists():
            return None
        for line in rst.read_text().splitlines():
            if line.strip().startswith("YN:"):
                f = line.split()[1:]
                return dict(zip(["S", "N", "t", "kappa", "omega", "dN", "SEdN", "dS", "SEdS"],
                                (float(x) if x.lower().lstrip("-") not in ("nan", "inf", "-nan(ind)") else math.nan
                                 for x in f)))
    return None


def classify_yn00(exe: str, cpp, cases, workdir: Path, tol: float, weighting: int):
    counts, examples, bad = Counter(), {}, []
    groups = {}
    for c in cases:
        if c.code in PAML_ICODE:
            groups.setdefault(c.code, []).append(c)
    tidy = {}
    for code, g in groups.items():
        axt = workdir / f"yn{code}.axt"
        P.write_axt(g, axt)
        r = P.run([cpp, "dnds", axt, "-q", "--digits", "15", "-c", code, "-m", "YN", "--yn-flavour", "paml",
                   "--yn-weighting", str(weighting)])
        tidy.update(P.parse_dnds(r.out))
    todo = [c for g in groups.values() for c in g]
    with ThreadPoolExecutor(2) as ex:
        refs = list(ex.map(lambda c: yn00_one(exe, c, weighting), todo))
    for c, y in zip(todo, refs):
        d = tidy.get((c.name, "YN"))
        if d is None:
            continue
        st = d["status"]
        if y is None:
            k = ("YN", "*", "yn00 no result" + ("" if st in ("no_data", "error") else ": unexplained"), st)
            counts[k] += 1
            examples.setdefault(k, c.name)
            if "unexplained" in k[2]:
                bad.append((c, "*", "none", d["dN"], st))
            continue
        if y.get("timeout"):
            counts[("YN", "*", "yn00 timeout", st)] += 1
            continue
        for name, col, yk in (("dN", "dN", "dN"), ("dS", "dS", "dS"), ("S", "S_sites", "S"), ("N", "N_sites", "N"),
                              ("SE_dN", "SE_dN", "SEdN"), ("SE_dS", "SE_dS", "SEdS"), ("dN/dS", "dNdS", "omega")):
            x, v = y[yk], num(d[col])
            xn = x if math.isfinite(x) else None
            if v is not None and v > 10 and name in ("dN", "dS"):
                kind = "COAL value above 10: rounding residue (K03)"
            elif xn is not None and v is not None:
                rel = abs(xn - v) / max(1.0, abs(xn), abs(v))
                if rel <= tol:
                    continue
                if name == "dN/dS" and num(d["dN"]) == 0 and abs(xn - 1e-5) < 1e-9:
                    kind = "number != number: omega floor 1e-5 when dN = 0"
                elif st == "no_data":
                    kind = "number != number: one kind of site (difference 8)"
                elif st == "saturated" and rel <= 1e-6:
                    kind = "number != number: saturated, other path (difference 13)"
                elif name.startswith("SE") and st == "saturated":
                    kind = "number != number: SE of a saturated row"
                else:
                    kind = "number != number: unexplained"
            elif xn is not None:
                if st == "not_converged" and name not in ("S", "N"):
                    kind = "number -> NA: not converged (difference 12)"
                elif st == "saturated":
                    kind = "number -> NA: saturated (difference 1)"
                elif st == "no_data":
                    kind = "number -> NA: one kind of site (difference 8)"
                elif name == "dN/dS" and (num(d["dS"]) is not None and num(d["dS"]) < 1e-6):
                    kind = "number -> NA: dN/dS with dS < 1e-6"
                elif name.startswith("SE") and st in ("ok", "warning"):
                    kind = "number -> NA: SE: unexplained"
                else:
                    kind = "number -> NA: unexplained"
            elif v is not None:
                if st == "no_data":
                    kind = "nan -> number: one kind of site (difference 8)"
                elif not math.isfinite(y["t"]):
                    # yn00's t is nan (its F84 has no value and it has no fallback); COAL's YN has the
                    # K80 / JC69 fallbacks of docs/kaks.md: recorded (BUGS.md, reference disagreements)
                    kind = "nan -> number: yn00 t is nan (COAL's fallback distance)"
                else:
                    kind = "nan -> number: unexplained"
            else:
                continue
            k = ("YN", name, kind, st)
            counts[k] += 1
            examples.setdefault(k, f"{c.name}: {x!r} / {d[col]}")
            if "unexplained" in kind or "K03" in kind and not KNOWN_K03_OK:
                bad.append((c, name, x, d[col], st))
    return counts, examples, bad


# ---------------------------------------------------------------------------------------------------- runner

def campaign(cpp, cases, kk=None, kk_fixed=None, yn00=None, tol=1e-9, keep=None, log=print):
    total_bad = []
    report = []
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        ours, tidy = coal_rows(cpp, cases, wd)
        byname = {c.name: c for c in cases}
        for label, exe in (("KaKs_Calculator", kk), ("KaKs_Calculator_corrected", kk_fixed)):
            if not exe:
                continue
            ref = {}
            for m in KK_METHODS:
                ref.update(kk_rows(exe, cases, m, wd, per_pair=True))
            counts, ex, bad = classify_kk(ref, ours, tidy, tol, byname)
            log(f"{label}: {len(cases)} pairs, {len(ref)} rows: {sum(counts.values())} differing fields; "
                f"{len(bad)} unexplained")
            for k in sorted(counts, key=lambda k: -counts[k]):
                report.append((label, *k, counts[k], ex[k]))
            for b in bad:
                total_bad.append((label, byname.get(b[0][0]), *b))
        if yn00:
            for w in (1, 0):
                counts, ex, bad = classify_yn00(yn00, cpp, cases, wd, max(tol, 1e-9), w)
                log(f"yn00 weighting {w}: {sum(counts.values())} differing fields; {len(bad)} unexplained")
                for k in sorted(counts, key=lambda k: -counts[k]):
                    report.append((f"yn00_w{w}", *k, counts[k], ex[k]))
                for b in bad:
                    total_bad.append((f"yn00_w{w}", b[0], (b[0].name, "YN"), *b[1:]))
    for r in report:
        log("  " + " | ".join(map(str, r)))
    if keep:
        kd = Path(keep)
        kd.mkdir(parents=True, exist_ok=True)
        with open(kd / "unexplained.tsv", "w", newline="\n") as fh:
            fh.write("reference\tpair\tmethod\tfield\treference_value\tcoal_value\tcoal_status\tcode\ts1\ts2\n")
            for lab, case, key, field, x, y, st in total_bad:
                fh.write("\t".join(map(str, [lab, key[0], key[1], field, x, y, st, case.code if case else "",
                                             case.s1 if case else "", case.s2 if case else ""])) + "\n")
        with open(kd / "classes.tsv", "w", newline="\n") as fh:
            fh.write("reference\tmethod\tfield\tclass\tcoal_status\tcount\texample\n")
            for r in report:
                fh.write("\t".join(map(str, r)) + "\n")
    return total_bad, report


KK = os.environ.get("KAKS_CALCULATOR")
YN00 = os.environ.get("YN00")


@unittest.skipUnless(KK or YN00, "set KAKS_CALCULATOR=/path/to/KaKs_Calculator and / or YN00=/path/to/yn00 "
                     "(builds printing 15 / 12 digits, docs/kaks.md Validation)")
class FuzzVsReference(unittest.TestCase):
    def test_random_pairs(self):
        cpp = P.cpp_bin()
        if cpp is None:
            self.skipTest("coal-cpp not built")
        n = int(os.environ.get("FUZZ_N", "1000" if P.SLOW else "120"))
        seed = int(os.environ.get("FUZZ_SEED", "20261114"))
        digits = int(os.environ.get("KAKS_REF_DIGITS", "15"))
        tol = 1e-9 if digits >= 12 else 1e-5
        cases = reference_cases(n, seed, 2000 if P.SLOW else 400)
        out = []
        bad, _ = campaign(cpp, cases, KK, os.environ.get("KAKS_CALCULATOR_FIXED"), YN00, tol,
                          os.environ.get("FUZZ_KEEP"), log=out.append)
        print("\n".join(out[:4]))
        self.assertFalse(bad, "\n".join(out) + "\n" + "\n".join(map(str, bad[:20])))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cpp", default=str(P.cpp_bin() or ""))
    ap.add_argument("--kaks")
    ap.add_argument("--kaks-fixed")
    ap.add_argument("--yn00")
    ap.add_argument("--n", type=int, default=600)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-codons", type=int, default=2000)
    ap.add_argument("--codes", help="comma-separated NCBI tables (default: the generator's mix)")
    ap.add_argument("--digits", type=int, default=15, help="significant digits the references print")
    ap.add_argument("--keep")
    a = ap.parse_args(argv)
    codes = [int(x) for x in a.codes.split(",")] if a.codes else None
    cases = reference_cases(a.n, a.seed, a.max_codons, codes)
    bad, _ = campaign(a.cpp, cases, a.kaks, a.kaks_fixed, a.yn00, 1e-9 if a.digits >= 12 else 1e-5, a.keep)
    print(f"fuzz vs reference: {len(cases)} pairs (seed {a.seed}), {len(bad)} unexplained differences")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
