#!/usr/bin/env python3
"""Compare the Ka/Ks methods of one COAL dnds.tsv that has several methods (-m YN -m NG ...).

Used by tutorial 5 (docs/tutorials/05-choosing-a-method.md). Standard library only.

  python3 compare_methods.py dnds.tsv [--reference YN]

It prints: the statuses per method; for the pairs where every method has a dN/dS, the median dS and dN/dS of
each method in bins of the reference method's dS; and, for the reference method, how the relative standard
error of dS and dN shrinks with the number of compared codons (methods with standard errors: YN, MYN, NG, Li).
"""
import argparse
import csv
import statistics


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dnds")
    ap.add_argument("--reference", default="YN", help="the method whose dS defines the bins (default YN)")
    a = ap.parse_args()

    with open(a.dnds, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    methods, by = [], {}
    for r in rows:
        if r["method"] not in methods:
            methods.append(r["method"])
        by.setdefault((r.get("cluster_id"), r["seq1"], r["seq2"]), {})[r["method"]] = r
    ref = a.reference
    if ref not in methods:
        raise SystemExit(f"method {ref} is not in {a.dnds} (it has {', '.join(methods)})")

    statuses = ["ok", "warning", "saturated", "not_converged", "no_data", "error"]
    print("statuses per method:")
    print("  " + f"{'method':<8}" + "".join(f"{s:>14}" for s in statuses))
    for m in methods:
        c = {}
        for d in by.values():
            if m in d:
                c[d[m]["status"]] = c.get(d[m]["status"], 0) + 1
        print("  " + f"{m:<8}" + "".join(f"{c.get(s, 0):>14}" for s in statuses))

    def has(r):
        return r["status"] in ("ok", "warning") and r["dNdS"] != "NA"

    common = [d for d in by.values() if all(m in d and has(d[m]) for m in methods)]
    print(f"\n{len(common)} pairs have a dN/dS from every method; median dS / median dN/dS, by {ref}'s dS:")
    print("  " + f"{ref + ' dS':<10}{'pairs':>6}" + "".join(f"{m:>18}" for m in methods))
    for lo, hi in [(0, 0.5), (0.5, 1), (1, 1.5), (1.5, 2), (2, float("inf"))]:
        sel = [d for d in common if lo <= float(d[ref]["dS"]) < hi]
        if not sel:
            continue
        cells = []
        for m in methods:
            ds = statistics.median(float(d[m]["dS"]) for d in sel)
            w = statistics.median(float(d[m]["dNdS"]) for d in sel)
            cells.append(f"{ds:.3f} / {w:.4f}")
        label = f"{lo}-{hi}" if hi != float("inf") else f">= {lo}"
        print("  " + f"{label:<10}{len(sel):>6}" + "".join(f"{c:>18}" for c in cells))

    se = [d[ref] for d in by.values() if ref in d and d[ref]["status"] == "ok" and d[ref]["SE_dS"] != "NA"
          and 0.5 <= float(d[ref]["dS"]) <= 1.5 and float(d[ref]["dN"]) > 0]
    if se:
        print(f"\n{ref}, pairs with status ok and 0.5 <= dS <= 1.5: median relative standard error, by compared codons:")
        for lo, hi in [(0, 100), (100, 200), (200, 400), (400, None)]:
            sel = [r for r in se if int(r["codons_compared"]) >= lo and (hi is None or int(r["codons_compared"]) < hi)]
            if not sel:
                continue
            rs = statistics.median(float(r["SE_dS"]) / float(r["dS"]) for r in sel)
            rn = statistics.median(float(r["SE_dN"]) / float(r["dN"]) for r in sel)
            label = f"{lo}-{hi}" if hi else f">= {lo}"
            print(f"  {label:<9} codons: {len(sel):>5} pairs   SE_dS/dS {rs:.2f}   SE_dN/dN {rn:.2f}")


if __name__ == "__main__":
    main()
