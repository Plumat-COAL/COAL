"""Validation of the methods Comeron and Li against their references (docs/kaks.md, "Comeron and Li").

  python benchmarks/kaks/validate_comeron_li.py --coal-cpp build/coal-cpp \
      --gestimator-ref work/gestimator_ref --rscript Rscript \
      --set ortho=val/ortho/pairs.axt --set sim=val/sim/pairs.axt ... --out work/

References (benchmarks/kaks/ref/):
  gestimator_ref   orthologr 0.4.3's unmodified copies of the libsequence headers behind its gestimator(),
                   compiled with a driver that prints 17 significant digits (gestimator_ref.cpp)
  seqinr_ref.R     seqinr::kaks() per pair, 17 significant digits

Generated sets (always added unless --no-generated): `short_single` and `short_multi` are the generators of
random_short_pairs.py; `short_stops` draws 1 to 12 codons from all 64 codons (stop codons included) with gaps,
N, IUPAC letters, lower case and U sprinkled in.

Compared, for every pair of every set:
  1. `coal dnds -m Comeron --comeron-flavour gestimator` with gestimator_ref: dN, dS, dN/dS. A value is
     "identical" when the 17-digit texts give the same double; gestimator's 999 must be NA in COAL.
  2. `coal dnds -m Li --li-flavour seqinr` with seqinr::kaks: ka, ks, and SE against sqrt(vka), sqrt(vks);
     seqinr's 9.999999 must be NA in COAL.
  3. the default flavours with the reference flavours: how many pairs change, why (columns removed by COAL's
     policy; a correction that is undefined) and by how much.
  4. sanity: Comeron and Li against COAL's MLPB, LPB, LWL and YN (which are KaKs_Calculator's numbers).
Exit status 1 when a value differs from its reference by more than 1e-12 (relative) or an NA does not match.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import random
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load_generators():
    spec = importlib.util.spec_from_file_location("rsp", os.path.join(HERE, "random_short_pairs.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_pairs_stops(n: int, seed: int):
    """Codons from all 64 (stops included), with gaps, ambiguity letters, lower case and U."""
    rng = random.Random(seed)
    codons = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
    stops = ["TAA", "TAG", "TGA"]
    out = []
    for k in range(n):
        ncod = rng.randint(1, 12)
        s1, s2 = [], []
        for i in range(ncod):
            c1 = rng.choice(codons)
            r = rng.random()
            if r < 0.35:
                c2 = c1
            elif r < 0.65:
                pos = rng.randrange(3)
                c2 = c1[:pos] + rng.choice("TCAG") + c1[pos + 1:]
            else:
                c2 = rng.choice(codons)
            if rng.random() < 0.06:
                c2 = rng.choice(stops)
            r = rng.random()
            if r < 0.05:
                c1 = "---"
            elif r < 0.08:
                c2 = c2[0] + "-" + c2[2]
            elif r < 0.11:
                pos = rng.randrange(3)
                c1 = c1[:pos] + rng.choice("NRYKMSWBDHVX") + c1[pos + 1:]
            elif r < 0.13:
                c2 = c2.replace("T", "U")
            elif r < 0.15:
                c1 = c1.lower()
            elif r < 0.16:
                c2 = "..."
            s1.append(c1)
            s2.append(c2)
        if rng.random() < 0.3:
            s1.append(rng.choice(stops))
            s2.append(rng.choice(stops) if rng.random() < 0.7 else "---")
        out.append((f"st{k}", "".join(s1), "".join(s2)))
    return out


def write_axt(path, pairs):
    with open(path, "w") as f:
        for name, a, b in pairs:
            f.write(f"{name}\n{a}\n{b}\n\n")


STOPS = {"TAA", "TAG", "TGA"}


def policy_applies(a, b):
    """Does COAL's policy read this pair differently from the reference programs? True when a column that the
    references compare (both codons A/C/G/T only) holds a stop codon (COAL drops or skips it), or when a column
    they skip holds U (COAL reads U as T)."""
    if len(a) != len(b) or len(a) % 3:
        return False
    for i in range(0, len(a), 3):
        x, y = a[i:i + 3].upper(), b[i:i + 3].upper()
        plain = set(x + y) <= set("ACGT")
        if plain and (x in STOPS or y in STOPS):
            return True
        if not plain and "U" in x + y and set((x + y).replace("U", "T")) <= set("ACGT"):
            return True
    return False


def read_axt(path):
    """name -> does the policy apply (dict keeps the input order)"""
    out, rec = {}, []

    def flush():
        if rec:
            out[rec[0]] = policy_applies(rec[1], rec[2]) if len(rec) == 3 else False

    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line:
                flush()
                rec = []
            else:
                rec.append(line)
    flush()
    return out


def run_coal(coal, axt, methods, extra, out):
    cmd = [coal, "dnds", axt, "--format", "axt", "-o", out, "--digits", "17", "-q", "-j", "2"]
    for m in methods:
        cmd += ["-m", m]
    subprocess.run(cmd + extra, check=True)
    rows = {}
    with open(out) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            v = line.rstrip("\n").split("\t")
            d = dict(zip(head, v))
            # AXT names: seq1 is the name up to the first '&', seq2 the rest (or NA)
            name = d["seq1"] if d["seq2"] == "NA" else d["seq1"] + "&" + d["seq2"]
            rows.setdefault(d["method"], {})[name] = d
    return rows


def fl(s):
    s = s.strip()
    return math.nan if s in ("NA", "NaN", "nan", "") else float(s)


def reldiff(a, b):
    if a == b:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


class Cmp:
    """One field compared over a set."""

    def __init__(self):
        self.ref_numbers = self.identical = self.within = self.na_both = 0
        self.max_rel = 0.0
        self.bad = []

    def add(self, name, ref, got, ref_is_na):
        if ref_is_na:
            if math.isnan(got):
                self.na_both += 1
            else:
                self.bad.append((name, "reference NA", got))
            return
        self.ref_numbers += 1
        if math.isnan(got):
            self.bad.append((name, ref, "NA"))
            return
        d = reldiff(ref, got)
        self.max_rel = max(self.max_rel, d)
        if d == 0.0:
            self.identical += 1
        elif d <= 1e-12:
            self.within += 1
        else:
            self.bad.append((name, ref, got))

    def cells(self):
        return [str(self.ref_numbers), str(self.identical), str(self.within), f"{self.max_rel:.1e}",
                str(self.na_both), str(len(self.bad))]


def pearson(x, y):
    if len(x) < 3:
        return math.nan
    mx, my = statistics.fmean(x), statistics.fmean(y)
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    if sx == 0 or sy == 0:
        return math.nan
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coal-cpp", required=True)
    ap.add_argument("--gestimator-ref", help="the binary built from ref/gestimator_ref.cpp")
    ap.add_argument("--rscript", help="Rscript with seqinr installed (omit to skip the seqinr comparison)")
    ap.add_argument("--set", action="append", default=[], metavar="NAME=FILE.axt")
    ap.add_argument("--no-generated", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sets = []
    for s in a.set:
        name, path = s.split("=", 1)
        sets.append((name, path))
    if not a.no_generated:
        gen = load_generators()
        for name, pairs in (("short_single", gen.make_pairs(4000, 5)),
                            ("short_multi", gen.make_pairs_multi(3000, 20261112)),
                            ("short_stops", make_pairs_stops(4000, 20261113))):
            path = os.path.join(a.out, name + ".axt")
            write_axt(path, pairs)
            sets.append((name, path))

    report = {"sets": {}}
    failed = False
    lines = []
    t1 = ["| Set | Pairs | Field | Reference numbers | Identical doubles | Within 1e-12 | Max rel. diff | NA in both | Mismatches |",
          "|---|---:|---|---:|---:|---:|---:|---:|---:|"]
    t2 = list(t1)
    t3 = ["| Set | Pairs | Pairs the policy reads differently | dN changes | dN > 1% | dN > 20% | dN number -> NA | dS changes | dS > 1% | dS > 20% | dS number -> NA | Unexplained |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    t3li = list(t3)
    t4 = ["| Set | Pair of methods | n | r (dN) | r (dS) | median ratio dN | median ratio dS | max abs diff dN | max abs diff dS | equal within 1e-6 (dN) | equal within 1e-6 (dS) |",
          "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for name, axt in sets:
        policy = read_axt(axt)
        names = list(policy)
        rep = report["sets"].setdefault(name, {"pairs": len(names)})
        leg = run_coal(a.coal_cpp, axt, ["Comeron", "Li"], ["--comeron-flavour", "gestimator", "--li-flavour", "seqinr"],
                       os.path.join(a.out, name + ".coal.ref_flavours.tsv"))
        dfl = run_coal(a.coal_cpp, axt, ["Comeron", "Li", "LWL", "LPB", "MLPB", "YN"], [],
                       os.path.join(a.out, name + ".coal.default.tsv"))
        # 1. gestimator
        if a.gestimator_ref:
            ref_path = os.path.join(a.out, name + ".gestimator_ref.tsv")
            with open(ref_path, "w") as f:
                subprocess.run([a.gestimator_ref, axt], check=True, stdout=f)
            cmp = {"dN": Cmp(), "dS": Cmp(), "dNdS": Cmp()}
            nerr = 0
            with open(ref_path) as f:
                for line in f:
                    v = line.rstrip("\n").split("\t")
                    row = leg["Comeron"].get(v[0])
                    if v[1] == "ERR":
                        nerr += 1
                        if row is None or row["status"] != "error":
                            cmp["dN"].bad.append((v[0], "reference refuses the pair", row and row["status"]))
                        continue
                    for k, col in (("dN", 1), ("dS", 2), ("dNdS", 3)):
                        r = float(v[col])
                        cmp[k].add(v[0], r, fl(row[k]), r == 999.0)
            for k in ("dN", "dS", "dNdS"):
                t1.append(f"| {name} | {len(names)} | {k} | " + " | ".join(cmp[k].cells()) + " |")
                if cmp[k].bad:
                    failed = True
                    lines.append(f"MISMATCH Comeron/gestimator {name} {k}: {cmp[k].bad[:10]}")
            rep["gestimator"] = {k: dict(zip(("ref_numbers", "identical", "within_1e-12", "max_rel", "na_both",
                                             "mismatches"), cmp[k].cells())) for k in cmp}
            rep["gestimator"]["refused_pairs"] = nerr
        # 2. seqinr
        if a.rscript:
            ref_path = os.path.join(a.out, name + ".seqinr_ref.tsv")
            subprocess.run([a.rscript, os.path.join(HERE, "ref", "seqinr_ref.R"), axt, ref_path], check=True,
                           stdout=subprocess.DEVNULL)
            cmp = {"dN": Cmp(), "dS": Cmp(), "SE_dN": Cmp(), "SE_dS": Cmp()}
            with open(ref_path) as f:
                for line in f:
                    v = [x.strip() for x in line.rstrip("\n").split("\t")]
                    row = leg["Li"].get(v[0])
                    if v[1] == "ERR":
                        if row is None or row["status"] != "error":
                            cmp["dN"].bad.append((v[0], "reference refuses the pair", row and row["status"]))
                        continue
                    ka, ks, vka, vks = (float(x) for x in v[1:5])
                    cmp["dN"].add(v[0], ka, fl(row["dN"]), ka == 9.999999)
                    cmp["dS"].add(v[0], ks, fl(row["dS"]), ks == 9.999999)
                    # the variance can be negative or NaN in seqinr (it is a formula, not a sum of squares):
                    # COAL's SE is then NA
                    for k, var, sat in (("SE_dN", vka, ka == 9.999999), ("SE_dS", vks, ks == 9.999999)):
                        na = sat or not (var >= 0)
                        cmp[k].add(v[0], math.sqrt(var) if not na else math.nan, fl(row[k]), na)
            for k in cmp:
                t2.append(f"| {name} | {len(names)} | {k} | " + " | ".join(cmp[k].cells()) + " |")
                if cmp[k].bad:
                    failed = True
                    lines.append(f"MISMATCH Li/seqinr {name} {k}: {cmp[k].bad[:10]}")
            rep["seqinr"] = {k: dict(zip(("ref_numbers", "identical", "within_1e-12", "max_rel", "na_both",
                                         "mismatches"), cmp[k].cells())) for k in cmp}
        # 3. default flavour against the reference flavour
        for meth, table in (("Comeron", t3), ("Li", t3li)):
            c = dict(cols=0, unexplained=0)
            for k in ("dN", "dS"):
                c[k] = dict(changed=0, p1=0, p20=0, to_na=0)
            examples = []
            for n in names:
                d, l = dfl[meth].get(n), leg[meth].get(n)
                if d is None or l is None or d["status"] == "error":
                    continue
                cols = policy[n]
                c["cols"] += cols
                for k in ("dN", "dS"):
                    x, y = fl(d[k]), fl(l[k])
                    if (math.isnan(x) and math.isnan(y)) or x == y:
                        continue
                    c[k]["changed"] += 1
                    if math.isnan(x) and not math.isnan(y):
                        c[k]["to_na"] += 1
                        # explained by the policy (columns removed) or by the saturation rule (the reference
                        # flavour's row says its value comes from an undefined correction)
                        if not cols and "undefined" not in l["note"] and "saturated" not in d["note"]:
                            c["unexplained"] += 1
                            examples.append((n, k, d[k], l[k]))
                    elif math.isnan(y):
                        if not cols:
                            c["unexplained"] += 1
                            examples.append((n, k, d[k], l[k]))
                    else:
                        rd = abs(x - y) / max(abs(y), 1e-300)
                        c[k]["p1"] += rd > 0.01
                        c[k]["p20"] += rd > 0.20
                        if not cols:
                            c["unexplained"] += 1
                            examples.append((n, k, d[k], l[k]))
            table.append(f"| {name} | {len(names)} | {c['cols']} | " + " | ".join(
                f"{c[k]['changed']} | {c[k]['p1']} | {c[k]['p20']} | {c[k]['to_na']}" for k in ("dN", "dS")) +
                f" | {c['unexplained']} |")
            rep[meth + "_default_vs_reference_flavour"] = c
            if c["unexplained"]:
                failed = True
                lines.append(f"UNEXPLAINED default/reference flavour {meth} {name}: {examples[:10]}")
        # 4. sanity against the other methods (pairs with status ok or warning in both)
        for m1, m2 in (("Comeron", "MLPB"), ("Comeron", "YN"), ("Comeron", "Li"), ("Li", "LPB"), ("Li", "MLPB"),
                       ("Li", "LWL"), ("Li", "YN")):
            xs = {"dN": ([], []), "dS": ([], [])}
            for n in names:
                r1, r2 = dfl[m1].get(n), dfl[m2].get(n)
                if not r1 or not r2:
                    continue
                for k in ("dN", "dS"):
                    x, y = fl(r1[k]), fl(r2[k])
                    if not math.isnan(x) and not math.isnan(y):
                        xs[k][0].append(x)
                        xs[k][1].append(y)
            cells = [str(min(len(xs["dN"][0]), len(xs["dS"][0])))]
            for k in ("dN", "dS"):
                cells.append(f"{pearson(*xs[k]):.4f}")
            for k in ("dN", "dS"):
                ratios = [x / y for x, y in zip(*xs[k]) if y > 0 and x > 0]
                cells.append(f"{statistics.median(ratios):.4f}" if ratios else "NA")
            for k in ("dN", "dS"):
                cells.append(f"{max((abs(x - y) for x, y in zip(*xs[k])), default=math.nan):.3g}")
            for k in ("dN", "dS"):
                cells.append(f"{sum(abs(x - y) <= 1e-6 for x, y in zip(*xs[k]))} / {len(xs[k][0])}")
            t4.append(f"| {name} | {m1} vs {m2} | " + " | ".join(cells) + " |")
            rep[f"{m1}_vs_{m2}"] = cells

    text = "\n".join(
        ["## Comeron, flavour gestimator, against orthologr 0.4.3's gestimator code (17 digits)", ""] + t1 +
        ["", "## Li, flavour seqinr, against seqinr::kaks (17 digits)", ""] + t2 +
        ["", "## Comeron: default flavour against flavour gestimator", ""] + t3 +
        ["", "## Li: default flavour against flavour seqinr", ""] + t3li +
        ["", "## Sanity: default flavours against COAL's other methods", ""] + t4 + [""] + lines)
    with open(os.path.join(a.out, "report.md"), "w") as f:
        f.write(text + "\n")
    with open(os.path.join(a.out, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(text)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
