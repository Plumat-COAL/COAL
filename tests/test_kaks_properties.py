"""Property-based fuzz tests of the Ka/Ks engine (`coal-cpp dnds` / `coal-cpp kaks`, `coal dnds` / `coal kaks`;
docs/kaks.md). A seeded generator (standard library only, so that the tests need nothing beyond the engine) makes
random pairs of aligned coding sequences of 1 to 2,000 codons, divergence 0 to 100 %, with gaps (whole and split
codons, `-` and `.`), internal and terminal stop codons, N / IUPAC / X, lower case, U, every genetic code COAL
supports, and a few malformed pairs (rows of different length, a length that is not a multiple of 3). For every
pair and method (NG, LWL, MLWL, LPB, MLPB, YN, MYN, Comeron, Li; the reference flavours too):

P1  no crash, no hang, no sanitizer report: the exit status is 0 or 2 (2 only when no pair has an estimate)
P2  status is one of ok, warning, saturated, not_converged, no_data, error
P3  no nan / inf in either output format; dN or dS NA only with a status that is not ok/warning; dN/dS NA only
    with dN or dS NA or dS < 1e-6 (gestimator: dS = 0); negative values only where documented (Li: warning)
P4  identical sequences: every printed dN and dS is 0, status ok / warning / no_data / error
P5  swapping the two sequences gives the same row (digits 17), for every method and flavour except the ones
    documented as order dependent (`--comeron-twosubs libsequence`, the default of `--comeron-flavour
    gestimator`); that documentation is checked to be there
P6  permuting the codon columns (pairs without stop codons) gives the same values to a relative 1e-9: every
    method here is a function of the multiset of codon pairs (NG, LWL, MLWL, LPB, MLPB, Comeron and Li sum per
    column; YN and MYN use whole-sequence codon frequencies and summed, pathway-weighted counts), so only the
    order of floating-point sums may change; for YN / MYN a status may flip between converged and not_converged
    (docs/kaks.md difference 12: the iteration amplifies the last bit), which is counted, not failed
P7  `coal dnds` / `coal kaks` (the Python command, plain and gzip-compressed input) print the bytes of
    `coal-cpp dnds` / `coal-cpp kaks`; AXT and FASTA (`--pairs consecutive`) input give the same values; the two
    output formats agree (Ka = dN, Ks = dS, Length = 3 x codons compared, the same codons for every method)
P8  the output is the same bytes across runs, for -j 1, 2 and 4, and for a pair alone or among others
P9  synonymous changes only (every compared codon pair identical or one synonymous change apart): dN = Nd = 0
    for NG, YN, MYN, MLWL and MLPB
P10 no dN or dS above 10: with fewer than a million sites only a rounding residue of an exact 0 inside the
    logarithm of a correction gives such a value

The fast subset (a few hundred pairs, fixed seed) takes well under 60 s. COAL_SLOW=1 runs more and longer pairs
(KAKS_FUZZ_N, KAKS_FUZZ_SEED override the count and the seed). Run it as a script for a campaign:

    python tests/test_kaks_properties.py --n 3000 --seed 7 [--max-codons 2000] [--keep DIR] [--cpp build/coal-cpp]

The bugs this fuzz found (K01-K06; tests/fixtures/kaks_fuzz/K0*.axt hold the minimal pairs) have their
reproducers in TestKnownBugs; a bug that is not yet fixed is an expected failure there, and the random tests
leave out exactly the rows it explains (KNOWN below), so that the suite stays green and still checks everything
else. K01, K02, K03 and K06 are fixed (2026-11-14). Skipped if build/coal-cpp
(COAL_CPP_BIN) is not built.
"""
from __future__ import annotations

import argparse
import gzip
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "tests" / "fixtures" / "kaks_fuzz"

METHODS = ["NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN", "Comeron", "Li"]
KK_METHODS = METHODS[:7]                    # the methods that run with every genetic code
STATUSES = {"ok", "warning", "saturated", "not_converged", "no_data", "error"}
ESTIMATE = {"ok", "warning"}
COLS = ["seq1", "seq2", "method", "dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "Nd", "Sd",
        "codons_compared", "status", "note"]
VALUE_COLS = COLS[3:13]
SE_METHODS = {"NG", "YN", "MYN", "Li"}
COUNT_METHODS = {"NG", "LWL", "MLWL", "LPB", "MLPB", "Comeron", "Li"}
CODES = (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23)
BASES = "TCAG"
CODON_LIST = [a + b + c for a in BASES for b in BASES for c in BASES]
# (the NCBI tables of COAL, src/coal/tables.py; copied so that this module needs no installed coal)
_AAS = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    13: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    14: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    15: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    16: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    21: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    22: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    23: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
AA = {n: dict(zip(CODON_LIST, s)) for n, s in _AAS.items()}
SENSE = {n: [c for c in CODON_LIST if AA[n][c] != "*"] for n in CODES}
STOP = {n: [c for c in CODON_LIST if AA[n][c] == "*"] for n in CODES}
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}
IUPAC = "NNNNRYKMSWBDHVX"
SLOW = bool(os.environ.get("COAL_SLOW"))


def cpp_bin() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


# ---------------------------------------------------------------------------------------------------- generator

@dataclass
class Case:
    name: str
    s1: str
    s2: str
    code: int = 1
    tags: set = field(default_factory=set)

    def swapped(self) -> "Case":
        return Case(self.name, self.s2, self.s1, self.code, set(self.tags) | {"swapped"})


def _length(rng: random.Random, max_codons: int) -> int:
    r = rng.random()
    lo, hi = ((1, 3) if r < 0.15 else (4, 12) if r < 0.40 else (13, 60) if r < 0.65 else (61, 300) if r < 0.85
              else (301, max(301, max_codons)))
    hi = min(hi, max_codons)
    lo = min(lo, hi)
    return int(round(math.exp(rng.uniform(math.log(lo), math.log(hi + 0.999)))))


def _pool(rng: random.Random, code: int) -> list[str]:
    """The codons of the first sequence: all sense codons with a random usage, or a reduced alphabet (codons of two or
    three bases, one amino acid's codons, Met/Trp), which reach sites of one kind only and frequencies of zero."""
    sense = SENSE[code]
    r = rng.random()
    if r < 0.6:
        return sense
    if r < 0.8:
        alpha = rng.choice(["TCA", "TCG", "CAG", "TAG", "TC", "AG", "AT", "GC"])
        pool = [c for c in sense if set(c) <= set(alpha)]
        return pool or sense
    if r < 0.9:
        aa = rng.choice(sorted({AA[code][c] for c in sense}))
        return [c for c in sense if AA[code][c] == aa]
    return rng.sample(sense, rng.randint(2, 6))


def _mutate(rng: random.Random, codon: str, d: float, kappa: float, code: int, stop_ok: float) -> str:
    out = list(codon)
    for i in range(3):
        if rng.random() < d:
            if rng.random() < kappa / (kappa + 2):
                out[i] = TRANSITION[out[i]]
            else:
                out[i] = rng.choice([b for b in "ACGT" if b != out[i] and b != TRANSITION[out[i]]])
    c = "".join(out)
    if AA[code][c] == "*" and rng.random() >= stop_ok:
        return codon
    return c


def _neighbours(c: str) -> list[str]:
    return [c[:i] + b + c[i + 1:] for i in range(3) for b in "TCAG" if b != c[i]]


def compared_columns(case: Case) -> list[tuple[str, str]]:
    """The codon columns COAL's policy compares (docs/kaks.md, Policy): both codons A/C/G/T (U read as T, any
    case) and neither a stop codon; a terminal stop is dropped like an internal one here."""
    out = []
    a, b = case.s1.upper().replace("U", "T"), case.s2.upper().replace("U", "T")
    for i in range(0, min(len(a), len(b)) - 2, 3):
        x, y = a[i:i + 3], b[i:i + 3]
        if set(x + y) <= set("ACGT") and AA[case.code][x] != "*" and AA[case.code][y] != "*":
            out.append((x, y))
    return out


def synonymous_single(case: Case) -> bool:
    """Every compared codon pair is identical or differs at one position without changing the amino acid."""
    cols = compared_columns(case)
    return bool(cols) and all(x == y or (sum(p != q for p, q in zip(x, y)) == 1 and AA[case.code][x] ==
                                         AA[case.code][y]) for x, y in cols)


def make_case(rng: random.Random, k: int, *, max_codons: int = 2000, code: int | None = None,
              plain: bool = False) -> Case:
    """One random pair. `plain`: no gaps, stops, ambiguity, case or U, no malformed rows (for the reference fuzz)."""
    code = code if code is not None else (1 if rng.random() < 0.55 else rng.choice(CODES))
    n = _length(rng, max_codons)
    pool = _pool(rng, code)
    weights = [rng.random() ** 2 + 0.01 for _ in pool]
    a = rng.choices(pool, weights, k=n)
    tags = set()
    r = rng.random()
    if r < 0.08:
        b, tags = list(a), {"identical"}
    elif r < 0.14:   # synonymous changes only, one position per codon at most
        b = []
        for c in a:
            alt = [x for x in _neighbours(c) if AA[code][x] == AA[code][c]]
            b.append(rng.choice(alt) if alt and rng.random() < 0.5 else c)
        tags.add("synonymous")
    elif r < 0.22:   # drawn independently: 100 % divergence in expectation
        b = rng.choices(pool if rng.random() < 0.5 else SENSE[code], k=n)
        tags.add("independent")
    else:
        d = rng.choice([rng.uniform(0, 0.02), rng.uniform(0, 0.1), rng.uniform(0, 0.4), rng.uniform(0, 1.0)])
        kappa = rng.choice([0.5, 1, 2, 4, 10])
        b = [_mutate(rng, c, d, kappa, code, 0.0 if plain else 0.02) for c in a]
    if not plain:
        if rng.random() < 0.12:   # internal stop codons
            for _ in range(rng.randint(1, max(1, n // 20))):
                (a if rng.random() < 0.5 else b)[rng.randrange(n)] = rng.choice(STOP[code])
            tags.add("stop")
        if rng.random() < 0.2:    # a terminal stop codon in one or both rows
            s = rng.choice(STOP[code])
            a.append(s)
            b.append(s if rng.random() < 0.7 else (rng.choice(STOP[code]) if rng.random() < 0.5
                                                  else rng.choice(SENSE[code])))
            tags.add("stop")
        n = len(a)
        if rng.random() < 0.2:    # whole-codon gaps, in one row or in both (an alignment column of gaps)
            for _ in range(rng.randint(1, max(1, n // 8))):
                i = rng.randrange(n)
                g = "..." if rng.random() < 0.1 else "---"
                if rng.random() < 0.5:
                    a[i] = g
                else:
                    b[i] = g
                if rng.random() < 0.1:
                    a[i] = b[i] = g
            tags.add("gap")
        if rng.random() < 0.06:   # a codon split by a gap
            i = rng.randrange(n)
            x = list(a[i])
            x[rng.randrange(3)] = "-"
            a[i] = "".join(x)
            tags.add("gap")
        if rng.random() < 0.1:    # N, IUPAC, X
            for _ in range(rng.randint(1, max(1, n // 10))):
                row = a if rng.random() < 0.5 else b
                i = rng.randrange(n)
                x = list(row[i])
                x[rng.randrange(3)] = rng.choice(IUPAC)
                row[i] = "".join(x)
            tags.add("ambiguous")
    s1, s2 = "".join(a), "".join(b)
    if not plain:
        if rng.random() < 0.08:
            s1, s2 = (s1.lower(), s2) if rng.random() < 0.5 else (s1.lower(), s2.lower())
            tags.add("lower")
        if rng.random() < 0.05:
            s1 = s1.replace("T", "U").replace("t", "u")
            tags.add("rna")
        r = rng.random()
        if r < 0.015:
            s2 = s2[:-3] if len(s2) > 3 else s2 + "ACG"
            tags.add("unequal")
        elif r < 0.03:
            cut = rng.choice([1, 2])
            s1, s2 = s1[:-cut] or "A", s2[:-cut] or "A"
            tags.add("not_triplet")
    return Case(f"c{k}", s1, s2, code, tags)


def make_cases(n: int, seed: int, **kw) -> list[Case]:
    rng = random.Random(seed)
    return [make_case(rng, k, **kw) for k in range(n)]


# ------------------------------------------------------------------------------------------------------ running

def write_axt(cases, path: Path) -> None:
    with open(path, "w", newline="\n") as fh:
        for c in cases:
            fh.write(f"{c.name}\n{c.s1}\n{c.s2}\n\n")


def write_fasta(cases, path: Path) -> None:
    with open(path, "w", newline="\n") as fh:
        for c in cases:
            fh.write(f">{c.name}_a\n{c.s1}\n>{c.name}_b\n{c.s2}\n")


@dataclass
class Run:
    cmd: list
    rc: int
    out: str
    err: str
    secs: float


def run(cmd, timeout=300, env=None) -> Run:
    t = time.monotonic()
    try:
        r = subprocess.run([str(x) for x in cmd], capture_output=True, text=True, timeout=timeout, env=env,
                           stdin=subprocess.DEVNULL)
    except subprocess.TimeoutExpired as e:
        return Run(cmd, -999, (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""),
                   "TIMEOUT", time.monotonic() - t)
    return Run(cmd, r.returncode, r.stdout, r.stderr, time.monotonic() - t)


def mflags(methods) -> list[str]:
    return sum((["-m", m] for m in methods), [])


def methods_for(code: int) -> list[str]:
    return METHODS if code == 1 else KK_METHODS


def parse_dnds(text: str) -> dict:
    """(pair name, method) -> {column: text}. AXT names here have no '&', so seq2 is NA and seq1 is the name;
    FASTA rows are keyed by seq1 without the _a suffix."""
    lines = text.splitlines()
    out = {}
    if not lines:
        return out
    head = lines[0].split("\t")
    for line in lines[1:]:
        f = line.split("\t")
        d = dict(zip(head, f))
        d["_nf"] = len(f)
        name = d["seq1"][:-2] if d["seq1"].endswith("_a") else d["seq1"]
        out[(name, d["method"])] = d
    return out


def parse_kaks(text: str) -> dict:
    out = {}
    for line in text.splitlines()[1:]:
        f = line.split("\t")
        out[(f[0], f[1])] = f
    return out


def num(x: str):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v


SAN_MARKERS = ("AddressSanitizer", "LeakSanitizer", "runtime error", "UndefinedBehaviorSanitizer",
               "terminate called", "Assertion", "Segmentation", "core dumped")


def crashed(r: Run) -> str:
    if r.rc == -999:
        return "timeout"
    if r.rc not in (0, 2) or r.rc < 0:
        return f"exit status {r.rc}"
    for m in SAN_MARKERS:
        if m in r.err:
            return f"stderr has {m!r}"
    return ""


# ------------------------------------------------------------------------------------------------------ checks

@dataclass
class Finding:
    prop: str
    case: Case
    method: str
    detail: str
    extra: tuple = ()

    def key(self):
        return (self.prop, self.method, self.detail.split(":")[0])


def check_row(case: Case, method: str, d: dict, flavour: dict) -> list[str]:
    """P2-P4 on one dnds row; returns the violated properties as 'Pn: text'."""
    bad = []
    st = d.get("status", "")
    if d.get("_nf") != len(COLS):
        bad.append(f"P3: {d.get('_nf')} columns")
    if st not in STATUSES:
        bad.append(f"P2: status {st!r}")
    for c in VALUE_COLS:
        v = d.get(c, "")
        low = v.lower()
        if "nan" in low or "inf" in low or v == "":
            bad.append(f"P3: {c} = {v!r}")
        elif v != "NA" and num(v) is None:
            bad.append(f"P3: {c} = {v!r} is not a number")
        elif v != "NA" and not math.isfinite(num(v)):
            bad.append(f"P3: {c} = {v!r}")
    if bad:
        return bad
    dn, ds, r = num(d["dN"]), num(d["dS"]), num(d["dNdS"])
    if st in ESTIMATE and (dn is None or ds is None):
        bad.append(f"P3: status {st} with dN {d['dN']} dS {d['dS']}")
    if st not in ESTIMATE and dn is not None and ds is not None and st != "no_data":
        # saturated / not_converged / error with both numbers: a status without its NA
        bad.append(f"P3: status {st} but dN {d['dN']} and dS {d['dS']} both printed")
    gest = flavour.get("comeron", "").startswith("gestimator") and method == "Comeron"
    if dn is not None and ds is not None:
        want_na = (ds == 0) if gest else (ds < 1e-6)
        if want_na and r is not None:
            bad.append(f"P3: dNdS {d['dNdS']} printed with dS {d['dS']}")
        if not want_na and r is None:
            bad.append(f"P3: dNdS NA with dN {d['dN']} dS {d['dS']}")
        if r is not None and ds != 0:
            if abs(r - dn / ds) > 1e-9 * max(1.0, abs(r)):
                bad.append(f"P3: dNdS {r} != dN/dS {dn / ds}")
    elif r is not None:
        bad.append(f"P3: dNdS {d['dNdS']} printed with dN {d['dN']} dS {d['dS']}")
    for c, v in (("dN", dn), ("dS", ds)):
        li_note = ("Ka" if c == "dN" else "Ks") + " is negative" in d["note"]
        if v is not None and v < 0 and not (method == "Li" and st != "ok" and li_note):
            bad.append(f"P3: negative {c} {v} (status {st})")
    for c in ("SE_dN", "SE_dS"):
        v = num(d[c])
        if v is not None and v < 0:
            bad.append(f"P3: negative {c} {v}")
        if method not in SE_METHODS and d[c] != "NA":
            bad.append(f"P3: {c} = {d[c]} for a method without a variance")
    if st == "not_converged" and (dn is not None or ds is not None):
        bad.append("P3: not_converged with a value")
    if st == "not_converged" and method not in ("YN", "MYN"):
        bad.append("P2: not_converged for a method without iteration")
    if "identical" in case.tags and not (case.tags & {"stop", "gap", "ambiguous", "rna", "lower", "unequal",
                                                     "not_triplet"}):
        if st not in ("ok", "warning", "no_data", "error"):
            bad.append(f"P4: identical sequences have status {st}")
        for c in ("dN", "dS", "Nd", "Sd"):
            v = num(d[c])
            if v is not None and v != 0:
                bad.append(f"P4: identical sequences have {c} = {d[c]}")
    if "unequal" in case.tags or "not_triplet" in case.tags:
        if st != "error":
            bad.append(f"P2: malformed pair has status {st}")
    elif st == "error":
        bad.append("P2: status error for a well-formed pair")
    # P9: no nonsynonymous difference when every compared codon pair is identical or one synonymous change apart
    # (NG, YN and MYN count such a change as synonymous; MLWL / MLPB by Tzeng et al.'s rules, standard code)
    if method in ("NG", "YN", "MYN", "MLWL", "MLPB") and \
            "unequal" not in case.tags and "not_triplet" not in case.tags and synonymous_single(case):
        for c in ("dN", "Nd"):
            v = num(d[c])
            if v is not None and v != 0:
                bad.append(f"P9: synonymous changes only, but {c} = {d[c]}")
    # P10: a multiple-hit correction (Jukes-Cantor, Kimura, F84 / TN93) of a pair of at most a million sites
    # cannot exceed about 8 unless the argument of its logarithm is a rounding residue of an exact 0 (p = 3/4,
    # 2P + Q = 1): such a value is not an estimate, whatever the status of the row. Two exceptions: the gestimator
    # flavour reproduces gestimator's number, residue included (docs/kaks.md "Differences", 14), and LWL / MLWL
    # divide Sd (Nd) by S (N), which with MLWL's kappa-based sites can be a small fraction of a site
    for c in ("dN", "dS"):
        v = num(d[c])
        sites = num(d["S_sites" if c == "dS" else "N_sites"])
        if gest or (method in ("LWL", "MLWL") and sites is not None and sites < 1):
            continue
        if v is not None and v > 10:
            bad.append(f"P10: {c} = {d[c]} with status {st}: a rounding residue at the boundary of the correction")
    cc = num(d["codons_compared"])
    if cc is None or cc < 0 or cc > len(case.s1) // 3 or cc != int(cc):
        if st != "error":
            bad.append(f"P3: codons_compared {d['codons_compared']}")
    return bad


def check_batch(cases, rows: dict, methods_of, flavour: dict) -> list[Finding]:
    out = []
    for c in cases:
        ms = methods_of(c)
        cc = set()
        for m in ms:
            d = rows.get((c.name, m))
            if d is None:
                out.append(Finding("P1", c, m, "row missing"))
                continue
            for b in check_row(c, m, d, flavour):
                out.append(Finding(b.split(":")[0], c, m, b))
            cc.add(d["codons_compared"])
        # the policy is applied once per pair, so the default flavours compare the same columns in every method
        if not flavour and len(cc) > 1:
            out.append(Finding("P7", c, "*", f"P7: codons_compared differs between methods: {sorted(cc)}"))
    return out


def rel_close(a: str, b: str, tol: float) -> bool:
    if a == b:
        return True
    x, y = num(a), num(b)
    if x is None or y is None:
        return False
    return abs(x - y) <= tol * max(1.0, abs(x), abs(y))


def permute(case: Case, rng: random.Random) -> Case:
    n = len(case.s1) // 3
    idx = list(range(n))
    rng.shuffle(idx)
    a = "".join(case.s1[3 * i:3 * i + 3] for i in idx)
    b = "".join(case.s2[3 * i:3 * i + 3] for i in idx)
    return Case(case.name, a, b, case.code, set(case.tags) | {"permuted"})


def has_stop(case: Case) -> bool:
    for i in range(0, len(case.s1) - 2, 3):
        for s in (case.s1, case.s2):
            x = s[i:i + 3].upper().replace("U", "T")
            if x in STOP[case.code]:
                return True
    return False


# Swapping the rows changes the order of some floating-point sums, so the last digits of a 17-digit value may
# change (counted, not failed); beyond these relative tolerances it is a finding. The YN / MYN iteration
# amplifies such a difference (it stops at a tolerance of 5e-4 / 5e-8 on t; docs/kaks.md difference 12).
SWAP_TOL = {m: 1e-9 for m in METHODS} | {"YN": 1e-6, "MYN": 1e-6}

# The flavours whose output depends on which sequence is first (docs/kaks.md, "Two-position codon differences").
ORDER_DEPENDENT = ("--comeron-twosubs libsequence", "--comeron-flavour gestimator")
FLAVOURS = [
    ([], {}),
    (["--yn-flavour", "paml"], {"yn": "paml"}),
    (["--yn-weighting", "0"], {"yn": "w0"}),
    (["--comeron-flavour", "gestimator", "--li-flavour", "seqinr"], {"comeron": "gestimator", "li": "seqinr"}),
    (["--comeron-flavour", "gestimator", "--comeron-twosubs", "comeron"], {"comeron": "gestimator-comeron"}),
    (["--comeron-twosubs", "libsequence", "--comeron-min-sites", "0"], {"comeron": "libsequence"}),
]


class Harness:
    """Runs batches of cases through coal-cpp and collects findings."""

    def __init__(self, cpp: Path, workdir: Path, jobs: int = 2, timeout: int = 600):
        self.cpp, self.dir, self.jobs, self.timeout = cpp, workdir, jobs, timeout
        self.findings: list[Finding] = []
        self.last_digits = []
        self.counts = Counter()
        self.n = 0

    def _file(self, stem: str) -> Path:
        self.n += 1
        return self.dir / f"{self.n:04d}_{stem}"

    def dnds(self, cases, code, extra=(), methods=None, jobs=None, digits=17, fmt="axt") -> tuple[Run, dict]:
        methods = methods or methods_for(code)
        path = self._file(f"c{code}.{fmt}")
        (write_axt if fmt == "axt" else write_fasta)(cases, path)
        cmd = [self.cpp, "dnds", path, "-q", "--digits", str(digits), "-c", str(code), "-j", str(jobs or self.jobs),
               *mflags(methods), *extra]
        if fmt == "fasta":
            cmd += ["--pairs", "consecutive"]
        r = run(cmd, self.timeout)
        why = crashed(r)
        if why:
            for c in cases[:1]:
                self.findings.append(Finding("P1", c, ",".join(methods), f"P1: {why}: {r.err[-400:]}",
                                             (tuple(map(str, cmd)), str(path))))
        if r.rc == 2 and r.out.count("\n") > 1:
            rows = parse_dnds(r.out)
            if any(d["status"] in ESTIMATE | {"saturated", "not_converged"} for d in rows.values()):
                self.findings.append(Finding("P1", cases[0], "*", "P1: exit status 2 although a pair has an estimate"))
        return r, parse_dnds(r.out)

    def by_code(self, cases):
        groups = {}
        for c in cases:
            groups.setdefault(c.code, []).append(c)
        return groups

    def properties(self, cases, flavours=FLAVOURS, swap=True) -> None:
        """P1-P5 for every flavour, by genetic code."""
        for code, group in sorted(self.by_code(cases).items()):
            for extra, fl in flavours:
                if fl and code != 1 and any(k in fl for k in ("comeron", "li")):
                    continue
                r, rows = self.dnds(group, code, extra)
                self.counts["rows"] += len(rows)
                self.findings += check_batch(group, rows, lambda c: methods_for(c.code), fl)
                for d in rows.values():
                    self.counts[f"status {d['status']}"] += 1
                if not swap:
                    continue
                r2, rows2 = self.dnds([c.swapped() for c in group], code, extra)
                dependent = any(" ".join(extra).find(x) >= 0 for x in ORDER_DEPENDENT)
                for c in group:
                    for m in methods_for(code):
                        a, b = rows.get((c.name, m)), rows2.get((c.name, m))
                        if a is None or b is None:
                            continue
                        self.counts["swap compared"] += 1
                        diff = [k for k in VALUE_COLS + ["status"] if a[k] != b[k]]
                        if not diff:
                            continue
                        if dependent and m == "Comeron":
                            self.counts["swap differs (documented order dependence)"] += 1
                            continue
                        if a["status"] != b["status"]:
                            if m in ("YN", "MYN") and "not_converged" in (a["status"], b["status"]):
                                self.counts["swap: YN/MYN convergence flips (difference 12)"] += 1
                                continue
                            big = diff
                        else:
                            big = [k for k in diff if not rel_close(a[k], b[k], SWAP_TOL[m])]
                        if not big:
                            self.counts["swap: last digits only"] += 1
                            self.last_digits.append((c, m, {k: (a[k], b[k]) for k in diff}))
                            continue
                        self.findings.append(Finding(
                            "P5", c, m, f"P5: swap changes {','.join(big)} with {' '.join(extra) or 'defaults'}",
                            (tuple(extra), {k: (a[k], b[k]) for k in diff}, a["status"])))

    def permutation(self, cases, seed: int) -> None:
        rng = random.Random(seed)
        pool = [c for c in cases if not has_stop(c) and not (c.tags & {"unequal", "not_triplet"})
                and len(c.s1) >= 6]
        for code, group in sorted(self.by_code(pool).items()):
            _, rows = self.dnds(group, code)
            perm = [permute(c, rng) for c in group]
            _, rows2 = self.dnds(perm, code)
            for c in group:
                for m in methods_for(code):
                    a, b = rows.get((c.name, m)), rows2.get((c.name, m))
                    if a is None or b is None:
                        continue
                    self.counts["permutation compared"] += 1
                    if a["status"] != b["status"]:
                        if m in ("YN", "MYN") and {a["status"], b["status"]} & {"not_converged"}:
                            self.counts["permutation: YN/MYN convergence flips (difference 12)"] += 1
                            continue
                        self.findings.append(Finding("P6", c, m, f"P6: permutation changes status {a['status']} -> "
                                                     f"{b['status']}", ({"a": a, "b": b}, perm[group.index(c)])))
                        continue
                    big = [k for k in VALUE_COLS if not rel_close(a[k], b[k], SWAP_TOL[m])]
                    if big:
                        self.findings.append(Finding("P6", c, m, f"P6: permutation changes {','.join(big)}",
                                                     ({k: (a[k], b[k]) for k in big}, perm[group.index(c)],
                                                      a["status"])))
                    elif any(a[k] != b[k] for k in VALUE_COLS):
                        self.counts["permutation: last digits only"] += 1

    def determinism(self, cases) -> None:
        """P8: -j 1 / 2 / 4, twice, and each pair alone against the batch, byte for byte (stdout)."""
        for code, group in sorted(self.by_code(cases).items()):
            path = self._file(f"det{code}.axt")
            write_axt(group, path)
            outs = {}
            for j in (1, 2, 4, 4):
                for sub in ("dnds", "kaks"):
                    cmd = ([self.cpp, sub, path] if sub == "dnds" else [self.cpp, sub, "-i", path]) + [
                        "-q", "-c", str(code), "-j", str(j), "--digits", "17", *mflags(methods_for(code))]
                    r = run(cmd, self.timeout)
                    outs.setdefault(sub, []).append((j, r.out))
            for sub, lst in outs.items():
                if len({o for _, o in lst}) != 1:
                    self.findings.append(Finding("P8", group[0], "*", f"P8: {sub} output depends on -j or the run "
                                                 f"(code {code})", (str(path),)))
            # reversed order of the pairs in the file: every row the same
            rev = self._file(f"rev{code}.axt")
            write_axt(group[::-1], rev)
            r = run([self.cpp, "dnds", rev, "-q", "-c", str(code), "-j", "2", "--digits", "17",
                     *mflags(methods_for(code))], self.timeout)
            fwd = outs["dnds"][0][1].splitlines()
            if sorted(fwd[1:]) != sorted(r.out.splitlines()[1:]):
                self.findings.append(Finding("P8", group[0], "*",
                                             f"P8: a row depends on the order of the pairs (code {code})"))
            self.counts["determinism batches"] += 1

    def formats(self, cases) -> None:
        """P7: AXT vs FASTA input, dnds vs kaks output."""
        for code, group in sorted(self.by_code(cases).items()):
            group = [c for c in group if not (c.tags & {"unequal", "not_triplet"})]
            if not group:
                continue
            _, ax = self.dnds(group, code, digits=15)
            _, fa = self.dnds(group, code, digits=15, fmt="fasta")
            path = self._file(f"k{code}.axt")
            write_axt(group, path)
            rk = run([self.cpp, "kaks", "-i", path, "-q", "-c", str(code), "--digits", "15",
                      *mflags(methods_for(code))], self.timeout)
            kk = parse_kaks(rk.out)
            if any(x in rk.out.lower() for x in ("nan", "inf")):
                self.findings.append(Finding("P3", group[0], "*", "P3: nan/inf in the kaks output", (str(path),)))
            for c in group:
                for m in methods_for(code):
                    a, b = ax.get((c.name, m)), fa.get((c.name, m))
                    if a is None or b is None:
                        self.findings.append(Finding("P7", c, m, "P7: row missing in one input format"))
                        continue
                    diff = [k for k in VALUE_COLS + ["status", "note"] if a[k] != b[k]]
                    if diff:
                        self.findings.append(Finding("P7", c, m, f"P7: AXT and FASTA input differ in {','.join(diff)}"))
                    f = kk.get((c.name, m))
                    if f is None:
                        self.findings.append(Finding("P7", c, m, "P7: row missing in the kaks output"))
                        continue
                    self.counts["format compared"] += 1
                    for col, kcol in (("dN", 2), ("dS", 3)):
                        v, w = a[col], f[kcol]
                        x = num(v)
                        if x is None or x < 1e-6:
                            if w != "NA":
                                self.findings.append(Finding("P7", c, m, f"P7: kaks {col}={w} where dnds has {v}"))
                        elif w != v:
                            self.findings.append(Finding("P7", c, m, f"P7: kaks {col}={w} where dnds has {v}"))
                    if a["codons_compared"] != "NA" and f[6] != str(3 * int(a["codons_compared"])):
                        self.findings.append(Finding("P7", c, m, f"P7: kaks Length {f[6]} for "
                                                     f"{a['codons_compared']} codons compared"))

    def python_cli(self, cases, python=sys.executable) -> None:
        """P7: `python -m coal dnds|kaks` prints coal-cpp's bytes, for plain and gzip-compressed input."""
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), COAL_CPP=str(self.cpp))
        for code, group in sorted(self.by_code(cases).items()):
            path = self._file(f"py{code}.axt")
            write_axt(group, path)
            gz = path.with_suffix(".axt.gz")
            gz.write_bytes(gzip.compress(path.read_bytes(), mtime=0))
            for sub in ("dnds", "kaks"):
                base = ([path] if sub == "dnds" else ["-i", path])
                args = ["-q", "-c", str(code), "--digits", "15", *mflags(methods_for(code))]
                want = run([self.cpp, sub, *base, *args], self.timeout)
                for inp in (path, gz):
                    b2 = [inp] if sub == "dnds" else ["-i", inp]
                    got = run([python, "-m", "coal", sub, *b2, *args], self.timeout, env=env)
                    self.counts["python cli compared"] += 1
                    if (got.rc, got.out) != (want.rc, want.out):
                        self.findings.append(Finding("P7", group[0], "*", f"P7: python -m coal {sub} differs from "
                                                     f"coal-cpp ({inp.name}; rc {got.rc} vs {want.rc})",
                                                     (got.err[-300:],)))

    def summary(self) -> str:
        byp = Counter(f.prop for f in self.findings)
        return (f"kaks properties: {self.counts['rows']} rows, {self.counts['swap compared']} swap, "
                f"{self.counts['permutation compared']} permutation, {self.counts['format compared']} format, "
                f"{self.counts['python cli compared']} python-cli comparisons; findings "
                + (", ".join(f"{p} {n}" for p, n in sorted(byp.items())) or "none"))


# ---------------------------------------------------------------------------------------------------- shrinking

def shrink(case: Case, still_fails, max_steps: int = 400) -> Case:
    """Smallest pair (by codon columns, greedy delta debugging) for which still_fails(case) holds."""
    cur = case
    steps = 0
    chunk = max(1, len(cur.s1) // 3 // 2)
    while chunk >= 1 and steps < max_steps:
        n = len(cur.s1) // 3
        progressed = False
        i = 0
        while i < n and steps < max_steps:
            a = cur.s1[:3 * i] + cur.s1[3 * (i + chunk):]
            b = cur.s2[:3 * i] + cur.s2[3 * (i + chunk):]
            steps += 1
            if a and len(a) % 3 == len(cur.s1) % 3:
                cand = Case(cur.name, a, b, cur.code, cur.tags)
                if still_fails(cand):
                    cur, n, progressed = cand, len(a) // 3, True
                    continue
            i += chunk
        if not progressed:
            chunk //= 2
    return cur


# -------------------------------------------------------------------------------------------------------- tests

# Known bugs found by this fuzz (tests/fixtures/kaks_fuzz/K*.axt). A property that a known bug breaks is checked
# on every other row, and on the bug's reproducer as an expected failure (TestKnownBugs).
KNOWN = {
    # K01, K02, K03 and K06 were fixed on 2026-11-14 (branch fix/kaks-accuracy; docs/kaks.md "Differences", 14 to
    # 17): their reproducers in TestKnownBugs pass, and the random tests check every row again. Add an entry here
    # (id -> predicate on a Finding) for a bug that is found and not yet fixed.
}


def known(f: Finding) -> str:
    for kid, pred in KNOWN.items():
        if pred(f):
            return kid
    return ""


class PropertyBase(unittest.TestCase):
    N, SEED, MAXC = (int(os.environ.get("KAKS_FUZZ_N", 1500 if SLOW else 160)),
                     int(os.environ.get("KAKS_FUZZ_SEED", 20261114)), 2000 if SLOW else 400)

    @classmethod
    def setUpClass(cls):
        cls.bin = cpp_bin()
        if cls.bin is None:
            if os.environ.get("COAL_REQUIRE_CPP"):
                raise AssertionError("coal-cpp not built (COAL_REQUIRE_CPP is set)")
            raise unittest.SkipTest("coal-cpp not built (cmake --build build)")
        cls.tmp = tempfile.mkdtemp(prefix="kaks_props_")
        cls.cases = make_cases(cls.N, cls.SEED, max_codons=cls.MAXC)

    @classmethod
    def tearDownClass(cls):
        if getattr(cls, "tmp", None):
            shutil.rmtree(cls.tmp, ignore_errors=True)

    def harness(self) -> Harness:
        d = Path(tempfile.mkdtemp(dir=self.tmp))
        return Harness(self.bin, d, jobs=2)

    def assertNoFindings(self, h: Harness, props):
        new = [f for f in h.findings if f.prop in props and not known(f)]
        lines = [f"{f.detail} [{f.case.name} code {f.case.code} method {f.method} tags {sorted(f.case.tags)}; "
                 f"{len(f.case.s1) // 3} codons: {f.case.s1[:60]} / {f.case.s2[:60]}]" for f in new[:15]]
        if new:
            self.fail(f"{len(new)} findings ({h.summary()}):\n" + "\n".join(lines))


class TestGenerator(unittest.TestCase):
    def test_deterministic_and_covers_the_input_space(self):
        a, b = make_cases(400, 5), make_cases(400, 5)
        self.assertEqual([(c.s1, c.s2, c.code) for c in a], [(c.s1, c.s2, c.code) for c in b])
        tags = Counter(t for c in a for t in c.tags)
        for t in ("identical", "synonymous", "independent", "stop", "gap", "ambiguous", "lower", "rna"):
            self.assertGreater(tags[t], 0, t)
        self.assertGreater(len({c.code for c in a}), 10)
        lengths = [len(c.s1) // 3 for c in make_cases(400, 5, max_codons=2000)]
        self.assertLessEqual(min(lengths), 1)
        self.assertGreater(max(lengths), 1000)

    def test_shrink(self):
        c = Case("x", "AAA" * 10 + "TGG" + "CCC" * 5, "AAG" * 10 + "TGG" + "CCA" * 5)
        small = shrink(c, lambda k: "TGG" in k.s1)
        self.assertEqual(small.s1, "TGG")


class TestProperties(PropertyBase):
    def test_p1_to_p5_every_flavour(self):
        h = self.harness()
        h.properties(self.cases)
        self.assertGreater(h.counts["rows"], 1000)
        self.assertNoFindings(h, {"P1", "P2", "P3", "P4", "P5", "P7", "P9", "P10"})

    def test_p6_codon_column_permutation(self):
        h = self.harness()
        h.permutation(self.cases, self.SEED + 1)
        self.assertGreater(h.counts["permutation compared"], 500)
        self.assertNoFindings(h, {"P1", "P6"})

    def test_p7_formats_and_python_cli(self):
        h = self.harness()
        sub = self.cases[: max(40, len(self.cases) // 4)]
        h.formats(sub)
        h.python_cli(sub[:30])
        self.assertGreater(h.counts["python cli compared"], 0)
        self.assertNoFindings(h, {"P1", "P3", "P7"})

    def test_p8_threads_runs_and_order(self):
        h = self.harness()
        h.determinism(self.cases[: max(60, len(self.cases) // 3)])
        self.assertNoFindings(h, {"P1", "P8"})

    def test_order_dependent_flavours_are_documented(self):
        """P5's exception: `--comeron-twosubs libsequence` (the default of the gestimator flavour) depends on
        which sequence is first. The documentation and the help say so, and it is still true (else the exception
        would hide nothing and should go)."""
        doc = " ".join((ROOT / "docs" / "kaks.md").read_text(encoding="utf-8").split())
        self.assertIn("the estimate depends on which sequence is first", doc)
        self.assertIn("| `libsequence` | libsequence's four branches, as above | `--comeron-flavour gestimator` |", doc)
        helptext = " ".join(run([self.bin, "dnds", "--help"]).out.split())
        self.assertIn("depend on which sequence is first", helptext)
        h = self.harness()
        rng = random.Random(3)   # a fifth of the codons redrawn: many differ at two positions, few saturate
        cases = []
        for k in range(10):
            a = rng.choices(SENSE[1], k=150)
            b = [rng.choice(SENSE[1]) if rng.random() < 0.2 else c for c in a]
            cases.append(Case(f"d{k}", "".join(a), "".join(b)))
        _, a = h.dnds(cases, 1, ["--comeron-twosubs", "libsequence"], methods=["Comeron"])
        _, b = h.dnds([c.swapped() for c in cases], 1, ["--comeron-twosubs", "libsequence"], methods=["Comeron"])
        self.assertTrue(any(a[k]["dN"] != b[k]["dN"] or a[k]["dS"] != b[k]["dS"] for k in a))

    def test_comeron_li_refused_with_other_codes(self):
        d = Path(tempfile.mkdtemp(dir=self.tmp))
        write_axt([Case("x", "ATGAAACCC", "ATGAAGCCA")], d / "x.axt")
        for m in ("Comeron", "Li"):
            for code in (2, 11):
                r = run([self.bin, "dnds", d / "x.axt", "-m", m, "-c", str(code)])
                self.assertEqual(r.rc, 2, r.err)
                self.assertIn("standard genetic code only", r.err)
                self.assertFalse(crashed(Run([], 0, "", r.err, 0)))


def read_fixture(path: Path):
    """(# key=value header, [Case]) of a corpus file: AXT records after '#' lines."""
    meta, recs = {}, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            for tok in line[1:].split():
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    meta.setdefault(k, v)
            continue
        recs.append(line.strip())
    code = int(meta.get("code", 1))
    cases, i = [], 0
    while i < len(recs):
        if not recs[i]:
            i += 1
            continue
        cases.append(Case(recs[i], recs[i + 1], recs[i + 2], code))
        i += 3
    return meta, cases


class TestCorpus(PropertyBase):
    """tests/fixtures/kaks_fuzz/corpus_*.axt: pairs the fuzz found interesting that behave as documented."""
    N = 0

    def test_corpus(self):
        files = sorted(FIXTURES.glob("corpus_*.axt"))
        self.assertTrue(files)
        for f in files:
            meta, cases = read_fixture(f)
            for extra, fl in FLAVOURS:
                with self.subTest(f.name, flavour=" ".join(extra)):
                    h = self.harness()
                    h.properties(cases, flavours=[(extra, fl)])
                    self.assertNoFindings(h, {"P1", "P2", "P3", "P4", "P5", "P7", "P9", "P10"})


class TestKnownBugs(PropertyBase):
    """The minimal reproductions of the bugs this fuzz found (tests/fixtures/kaks_fuzz/K*.axt; the campaign's
    BUGS.md has the analysis). Each test asserts the correct behaviour and is an expected failure: it becomes
    an "unexpected success" (remove the decorator) once the bug is fixed."""
    N = 0

    def rows(self, name, code, methods, extra=()):
        meta, cases = read_fixture(FIXTURES / name)
        h = self.harness()
        r, rows = h.dnds(cases, code, list(extra), methods=methods)
        self.assertFalse(crashed(r), r.err)
        return cases, rows

    def test_K05_an_input_without_pairs_exits_2(self):
        """docs/kaks.md, Exit status: 2 when no pair has an estimate because of the input."""
        d = Path(tempfile.mkdtemp(dir=self.tmp))
        for name, text in (("empty.axt", ""), ("blank.axt", "\n\n\n"), ("comment.axt", "# only a comment\n")):
            (d / name).write_text(text)
            for sub in ("dnds", "kaks"):
                r = run([self.bin, sub, *([] if sub == "dnds" else ["-i"]), d / name, "-q"])
                self.assertFalse(crashed(r), r.err)
                self.assertEqual(r.rc, 2, f"{sub} {name}")

    def test_K06_saturated_myn_counts_do_not_depend_on_the_order(self):
        _, rows = self.rows("K06_myn_saturated_counts_swap_code5.axt", 5, ["MYN"])
        a, b = rows[("f", "MYN")], rows[("r", "MYN")]
        self.assertEqual((a["status"], b["status"]), ("saturated", "saturated"))
        self.assertEqual((a["Nd"], a["Sd"]), (b["Nd"], b["Sd"]))

    def test_K01_lwl_lpb_do_not_depend_on_the_order_of_the_sequences(self):
        _, rows = self.rows("K01_lwl_lpb_swap_code13.axt", 13, ["LWL", "LPB"])
        for m in ("LWL", "LPB"):
            a, b = rows[("f", m)], rows[("r", m)]
            self.assertEqual([a[k] for k in VALUE_COLS], [b[k] for k in VALUE_COLS], m)

    def test_K02_mlwl_mlpb_no_nonsynonymous_difference_for_synonymous_changes(self):
        _, rows = self.rows("K02_mlwl_mlpb_synonymous_ata_atg_code2.axt", 2, ["NG", "MLWL", "MLPB"])
        self.assertEqual(rows[("k2", "NG")]["dN"], "0")
        for m in ("MLWL", "MLPB"):
            self.assertEqual(rows[("k2", m)]["dN"], "0", m)

    def test_K03_an_exact_zero_of_the_correction_is_saturated(self):
        ms = ["NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN", "Comeron"]
        _, rows = self.rows("K03_kimura_boundary_residue.axt", 1, ms)
        for m in ms[1:5] + ["Comeron"]:
            self.assertEqual(rows[("b4", m)]["status"], "saturated", m)   # exact in binary: right already
        found = [(k, d["dN"], d["dS"], d["status"]) for k, d in rows.items()
                 if any((num(d[c]) or 0) > 10 for c in ("dN", "dS"))]
        self.assertEqual(found, [])
        self.assertEqual((rows[("ng1", "NG")]["status"], rows[("ng1", "NG")]["dS"]), ("saturated", "NA"))


# ------------------------------------------------------------------------------------------------------ script

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cpp", default=str(cpp_bin() or ""))
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-codons", type=int, default=2000)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--keep", help="folder for the findings (findings.tsv, one AXT per finding, shrunk)")
    ap.add_argument("--only", default="props,perm,formats,python,det")
    ap.add_argument("--timeout", type=int, default=900)
    a = ap.parse_args(argv)
    if not a.cpp:
        print("coal-cpp not found: --cpp", file=sys.stderr)
        return 2
    cases = make_cases(a.n, a.seed, max_codons=a.max_codons)
    with tempfile.TemporaryDirectory() as td:
        h = Harness(Path(a.cpp), Path(td), a.jobs, a.timeout)
        only = a.only.split(",")
        t = time.monotonic()
        if "props" in only:
            h.properties(cases)
        if "perm" in only:
            h.permutation(cases, a.seed + 1)
        if "formats" in only:
            h.formats(cases)
        if "python" in only:
            h.python_cli(cases[:200])
        if "det" in only:
            h.determinism(cases[:400])
        print(h.summary(), f"({time.monotonic() - t:.0f} s, seed {a.seed}, {a.n} pairs)")
        groups = {}
        for f in h.findings:
            groups.setdefault(f.key(), []).append(f)
        for k, fs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            f = min(fs, key=lambda x: len(x.case.s1))
            print(f"  {len(fs):5d}  {f.detail[:160]}  e.g. {f.case.name} code {f.case.code} {f.method} "
                  f"{len(f.case.s1) // 3} codons tags {sorted(f.case.tags)}{' known ' + known(f) if known(f) else ''}")
        if a.keep:
            keep = Path(a.keep)
            keep.mkdir(parents=True, exist_ok=True)
            with open(keep / "findings.tsv", "w", newline="\n") as fh:
                fh.write("prop\tmethod\tdetail\tcase\tcode\ttags\ts1\ts2\textra\n")
                for f in h.findings:
                    fh.write("\t".join([f.prop, f.method, f.detail.replace("\t", " "), f.case.name,
                                        str(f.case.code), ",".join(sorted(f.case.tags)), f.case.s1, f.case.s2,
                                        repr(f.extra)[:2000].replace("\t", " ")]) + "\n")
        return 1 if any(not known(f) for f in h.findings) else 0


if __name__ == "__main__":
    sys.exit(main())
