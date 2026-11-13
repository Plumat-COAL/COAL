"""Compare `coal kaks` with KaKs_Calculator 2.0 and PAML yn00 on an AXT file (docs/kaks.md "Validation").

  python benchmarks/kaks/validate_kaks.py --axt pairs.axt --coal-cpp build/coal-cpp \
      --kaks KaKs_Calculator [--kaks-hp KaKs_Calculator_15digits] [--yn00 yn00_hp] --out results/

* --kaks: the stock KaKs_Calculator 2.0 (6 significant digits). Every row of `coal kaks` (default --digits 6)
  is compared as text with KaKs_Calculator's row for the same pair and method: "identical rows".
* --kaks-hp: a KaKs_Calculator built with 15 significant digits (base.h CONVERT: stream.precision(15)), for the
  numeric comparison (max / 95th percentile absolute difference, Pearson r) against `coal kaks --digits 15`.
* --yn00: a PAML yn00 whose rst line prints S N t kappa omega dN SEdN dS SEdS with %.12g (the stock binary
  prints 4 decimals). Run per pair (yn00 carries t over between data sets), with weighting = 1 and 0, against
  `coal dnds --yn-flavour paml --yn-weighting 1|0`. Columns with a gap, an ambiguous base or a stop codon in
  either sequence are removed first (yn00 refuses stop codons; COAL drops those columns too).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import math
import os
import subprocess
import tempfile
from pathlib import Path

METHODS = ["NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN"]
STOPS1 = {"TAA", "TAG", "TGA"}


def read_axt(path: Path):
    lines = path.read_text().splitlines()
    i, out = 0, []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        name = lines[i].strip()
        rows = []
        i += 1
        while i < len(lines) and lines[i].strip():
            rows.append(lines[i].strip())
            i += 1
        out.append((name, rows[0], rows[1]))
    return out


def read_kaks(path: Path) -> dict[tuple[str, str], list[str]]:
    rows = {}
    for line in path.read_text().splitlines()[1:]:
        f = line.split("\t")
        if len(f) > 2:
            rows[(f[0], f[1])] = f
    return rows


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise SystemExit(f"{cmd} failed: {r.stderr}")
    return r


def num(x: str):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except ValueError:
        return None


def stats(pairs):
    """pairs: list of (a, b) floats or None. -> dict"""
    both = [(a, b) for a, b in pairs if a is not None and b is not None]
    na_agree = sum(1 for a, b in pairs if a is None and b is None)
    na_diff = [(i, a, b) for i, (a, b) in enumerate(pairs) if (a is None) != (b is None)]
    d = sorted(abs(a - b) for a, b in both)
    r = float("nan")
    if len(both) > 2:
        xs, ys = zip(*both)
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        sxy = sum((x - mx) * (y - my) for x, y in both)
        sxx = sum((x - mx) ** 2 for x in xs)
        syy = sum((y - my) ** 2 for y in ys)
        r = sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")
    return {
        "n": len(both), "na_agree": na_agree, "na_diff": len(na_diff), "na_diff_list": na_diff,
        "max": d[-1] if d else 0.0, "p95": d[int(0.95 * (len(d) - 1))] if d else 0.0, "r": r,
        "over_1e4": sum(1 for x in d if x > 1e-4),
    }


def filter_cols(r1: str, r2: str) -> tuple[str, str]:
    a, b = [], []
    for i in range(0, len(r1) - len(r1) % 3, 3):
        c1, c2 = r1[i:i + 3].upper().replace("U", "T"), r2[i:i + 3].upper().replace("U", "T")
        if any(ch not in "ACGT" for ch in c1 + c2) or c1 in STOPS1 or c2 in STOPS1:
            continue
        a.append(c1)
        b.append(c2)
    return "".join(a), "".join(b)


def lfact(n: float) -> float:
    return math.lgamma(n + 1) if n > 0 else 0.0


def p_observed(sd: float, nd: float, S: float, N: float) -> float:
    m = [sd, S - sd, nd, N - nd]
    R0, R1, C0, C1 = m[0] + m[2], m[1] + m[3], m[0] + m[1], m[2] + m[3]
    return math.exp(lfact(R0) + lfact(R1) + lfact(C0) + lfact(C1) - lfact(R0 + R1) - sum(lfact(x) for x in m))


def p_observed_explains(row, kc: str, co: str) -> bool:
    """COAL's Fisher p = KaKs_Calculator's + P(observed table) (to 6 digits): KaKs_Calculator lost the observed
    table to rounding in its 'no more probable than observed' test."""
    try:
        sd, nd, S, N = float(row[11]), float(row[12]), float(row[7]), float(row[8])
        kcv, cov = float(kc), float(co)
    except ValueError:
        return False
    return abs(cov - kcv - p_observed(sd, nd, S, N)) <= 1e-5 * cov


def negative_cell(row) -> bool:
    try:
        return float(row[7]) < float(row[11]) or float(row[8]) < float(row[12])
    except ValueError:
        return False


def yn00_one(yn00: str, r1: str, r2: str, weighting: int):
    s1, s2 = filter_cols(r1, r2)
    if not s1:
        return None
    with tempfile.TemporaryDirectory() as td:
        Path(td, "in.phy").write_text(f"  2  {len(s1)}\ns1  {s1}\ns2  {s2}\n")
        Path(td, "yn00.ctl").write_text(
            "seqfile = in.phy\noutfile = out.txt\nverbose = 0\nnoisy = 0\nicode = 0\n"
            f"weighting = {weighting}\ncommonkappa = 0\ncommonf3x4 = 0\nndata = 1\n")
        r = subprocess.run([yn00, "yn00.ctl"], cwd=td, capture_output=True, text=True, stdin=subprocess.DEVNULL)
        rst = Path(td, "rst")
        if not rst.exists():
            return None
        for line in rst.read_text().splitlines():
            if line.strip().startswith("YN:"):
                f = line.split()[1:]
                return dict(zip(["S", "N", "t", "kappa", "omega", "dN", "SEdN", "dS", "SEdS"], map(float, f)))
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--axt", required=True)
    ap.add_argument("--coal-cpp", required=True)
    ap.add_argument("--kaks")
    ap.add_argument("--kaks-hp")
    ap.add_argument("--yn00")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--methods", default=",".join(METHODS))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    axt = Path(a.axt).resolve()
    methods = a.methods.split(",")
    mflags = sum((["-m", m] for m in methods), [])
    report = []
    pairs = read_axt(axt)
    report.append(f"# {axt.name}: {len(pairs)} pairs\n")

    status = {}
    if a.kaks or a.kaks_hp:
        tsv = out / "coal_status.tsv"
        run([a.coal_cpp, "dnds", str(axt), "-o", str(tsv), "-q", *mflags])
        for line in tsv.read_text().splitlines()[1:]:
            f = line.split("\t")
            status[(f[0] + ("&" + f[1] if f[1] != "NA" else ""), f[2])] = f[13]
        run([a.coal_cpp, "kaks", "-i", str(axt), "-o", str(out / "coal15.kaks"), "-q", "--digits", "15", *mflags])
        ours15 = read_kaks(out / "coal15.kaks")

    if a.kaks:
        run([a.coal_cpp, "kaks", "-i", str(axt), "-o", str(out / "coal6.kaks"), "-q", *mflags])
        ours = read_kaks(out / "coal6.kaks")
        report.append("## coal kaks vs stock KaKs_Calculator 2.0 (text of each row, 6 significant digits)\n")
        report.append("Rows that differ are classified field by field: `-0` (KaKs_Calculator prints a negative zero), "
                      "`nan` (it prints nan where COAL prints NA), `p-obs` (Fisher p: KaKs_Calculator lost the observed "
                      "table to rounding; COAL - KaKs_Calculator = P(observed table)), `p-neg` (Fisher p: COAL gives NA because the "
                      "table has a negative cell, S < Sd or N < Nd, as LPB93's S = Sd/Ks once Ks > 1), `saturated` (COAL: NA with status "
                      "saturated; KaKs_Calculator prints a capped or truncated value), `not converged` (COAL: NA with status "
                      "not_converged, the YN / MYN iteration did not converge; KaKs_Calculator prints its last round), "
                      "`other` (unexplained).\n")
        report.append("| Method | rows | identical | different | -0 | nan | p-obs | p-neg | saturated | not converged | other |")
        report.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for m in methods:
            run([a.kaks, "-i", str(axt), "-o", str(out / f"kc_{m}.kaks"), "-m", m])
            ref = {k: v for k, v in read_kaks(out / f"kc_{m}.kaks").items()}
            mine = {k: v for k, v in ours.items() if k[1] == m}
            same = sum(1 for k in ref if k in mine and ref[k] == mine[k])
            diff = [k for k in ref if k in mine and ref[k] != mine[k]]
            cls = {"-0": 0, "nan": 0, "p-obs": 0, "p-neg": 0, "saturated": 0, "not converged": 0, "other": 0}
            with open(out / f"diff6_{m}.txt", "w") as f:
                for k in diff:
                    kinds = set()
                    for i, (x, y) in enumerate(zip(ref[k], mine[k])):
                        if x == y:
                            continue
                        if y == "NA" and status.get(k) == "not_converged" and i in (2, 3, 4, 5, 11, 12, 15):
                            kinds.add("not converged")
                        elif x == "-0" and y == "0":
                            kinds.add("-0")
                        elif x in ("nan", "-nan") and y == "NA":
                            kinds.add("nan")
                        elif i == 5 and p_observed_explains(ours15[k], x, ours15[k][5]):
                            kinds.add("p-obs")
                        elif i == 5 and y == "NA" and negative_cell(ours15[k]):
                            kinds.add("p-neg")
                        elif y == "NA" and status.get(k) == "saturated":
                            kinds.add("saturated")
                        else:
                            kinds.add("other")
                    for kd in kinds:
                        cls[kd] += 1
                    f.write(f"[{','.join(sorted(kinds))}]\nKC  " + "\t".join(ref[k]) + "\nCOAL" + "\t".join(mine[k]) + "\n")
            only = len(set(ref) ^ set(mine))
            report.append(f"| {m} | {len(ref)} | {same} | {len(diff)} | {cls['-0']} | {cls['nan']} | {cls['p-obs']} | {cls['p-neg']} | "
                          f"{cls['saturated']} | {cls['not converged']} | {cls['other']} |" + (f" ({only} rows in one tool only)" if only else ""))
        report.append("")

    if a.kaks_hp:
        ours = ours15
        report.append("## coal kaks vs KaKs_Calculator 2.0 at 15 significant digits\n")
        report.append("| Method | Value | both numeric | both NA | NA in one | max abs diff | 95th pct | r | > 1e-4 |")
        report.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
        for m in methods:
            run([a.kaks_hp, "-i", str(axt), "-o", str(out / f"kchp_{m}.kaks"), "-m", m])
            ref = read_kaks(out / f"kchp_{m}.kaks")
            keys = [k for k in ref if k in ours]
            for col, name in ((2, "Ka"), (3, "Ks"), (4, "Ka/Ks"), (5, "P-Value")):
                vals = [(num(ref[k][col]), num(ours[k][col])) for k in keys]
                if col == 5:  # leave out the p-values explained by KaKs_Calculator losing the observed table
                    nexp = sum(1 for k, (x, y) in zip(keys, vals) if x is not None and y is not None and
                               abs(x - y) > 1e-9 * max(abs(y), 1e-300) and p_observed_explains(ours[k], ref[k][col], ours[k][col]))
                    vals = [(x, y) for k, (x, y) in zip(keys, vals) if not (x is not None and y is not None and
                            abs(x - y) > 1e-9 * max(abs(y), 1e-300) and p_observed_explains(ours[k], ref[k][col], ours[k][col]))]
                    nneg = sum(1 for k in keys if negative_cell(ours[k]))
                    name = f"P-Value (excl. {nexp} p-obs)" + (f" ({nneg} p-neg are the NA-in-one)" if nneg else "")
                st = stats(vals)
                report.append(f"| {m} | {name} | {st['n']} | {st['na_agree']} | {st['na_diff']} | {st['max']:.2e} | "
                              f"{st['p95']:.2e} | {st['r']:.9f} | {st['over_1e4']} |")
                with open(out / f"hp_nadiff_{m}_{name.replace('/', '')}.txt", "w") as f:
                    for i, x, y in st["na_diff_list"]:
                        f.write(f"{keys[i][0]}\tKC={ref[keys[i]][col]}\tCOAL={ours[keys[i]][col]}\n")
        report.append("")

    if a.yn00:
        report.append("## coal dnds (YN, --yn-flavour paml) vs PAML yn00 (12 significant digits)\n")
        report.append("| weighting | Value | both numeric | both NA | NA in one | max abs diff | 95th pct | r | > 1e-4 |")
        report.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
        for w in (1, 0):
            tsv = out / f"coal_paml_w{w}.tsv"
            run([a.coal_cpp, "dnds", str(axt), "-o", str(tsv), "-q", "--digits", "15", "-m", "YN",
                 "--yn-flavour", "paml", "--yn-weighting", str(w)])
            mine = [l.split("\t") for l in tsv.read_text().splitlines()[1:]]
            with cf.ThreadPoolExecutor(a.threads) as ex:
                refs = list(ex.map(lambda p: yn00_one(a.yn00, p[1], p[2], w), pairs))
            with open(out / f"yn00_w{w}.tsv", "w") as f:
                f.write("pair\tdN\tdS\tomega\tS\tN\tSEdN\tSEdS\n")
                for p, r in zip(pairs, refs):
                    if r:
                        f.write(f"{p[0]}\t{r['dN']}\t{r['dS']}\t{r['omega']}\t{r['S']}\t{r['N']}\t{r['SEdN']}\t{r['SEdS']}\n")
            for col, key, name in ((3, "dN", "dN"), (4, "dS", "dS"), (9, "S", "S sites"), (8, "N", "N sites"),
                                   (6, "SEdN", "SE dN"), (7, "SEdS", "SE dS")):
                st = stats([(r[key] if r else None, num(mi[col])) for r, mi in zip(refs, mine)])
                report.append(f"| {w} | {name} | {st['n']} | {st['na_agree']} | {st['na_diff']} | {st['max']:.2e} | "
                              f"{st['p95']:.2e} | {st['r']:.9f} | {st['over_1e4']} |")
            # omega: PAML's own omega
            st = stats([(r["omega"] if r and r["dS"] > 1e-6 else None, num(mi[5])) for r, mi in zip(refs, mine)])
            report.append(f"| {w} | dN/dS | {st['n']} | {st['na_agree']} | {st['na_diff']} | {st['max']:.2e} | "
                          f"{st['p95']:.2e} | {st['r']:.9f} | {st['over_1e4']} |")
        report.append("")
    text = "\n".join(report) + "\n"
    (out / "report.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
