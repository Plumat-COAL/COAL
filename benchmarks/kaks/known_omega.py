"""Recovery of a known omega: `coal dnds` on simulate_pairs.py output, summarised per true omega and t.

  python benchmarks/kaks/simulate_pairs.py --out sim.axt ...
  coal dnds sim.axt -m NG -m LWL -m LPB -m YN -m MYN -o sim.tsv
  python benchmarks/kaks/known_omega.py sim.tsv
Prints, per method, true omega and t: pairs, saturated/NA, median and mean of dN/dS, and the median dS.
"""
from __future__ import annotations

import re
import statistics
import sys
from collections import defaultdict


def main() -> int:
    rows = defaultdict(list)
    for line in open(sys.argv[1]).read().splitlines()[1:]:
        f = line.split("\t")
        m = re.match(r"sim_w([0-9.]+)_t([0-9.]+)_k([0-9.]+)_\d+", f[0])
        if not m:
            continue
        rows[(f[2], float(m.group(1)), float(m.group(2)))].append(f)
    print("| Method | true omega | t (subst/codon) | pairs | NA | median dN/dS | mean dN/dS | median dS |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    for (meth, w, t), fs in sorted(rows.items()):
        vals = [float(f[5]) for f in fs if f[5] != "NA"]
        ds = [float(f[4]) for f in fs if f[4] != "NA"]
        med = statistics.median(vals) if vals else float("nan")
        mean = statistics.fmean(vals) if vals else float("nan")
        print(f"| {meth} | {w:g} | {t:g} | {len(fs)} | {len(fs) - len(vals)} | {med:.3f} | {mean:.3f} | "
              f"{statistics.median(ds) if ds else float('nan'):.3f} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
