"""Plots for results/windows_<commit>/ from runs.csv (warm medians of timed runs; warm-up runs excluded).
usage: plots_windows.py RESULTS_DIR"""
import csv
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

D = Path(sys.argv[1])
rows = list(csv.DictReader(open(D / "runs.csv", encoding="utf-8")))
TOOLS = ["p2n14", "p2n16", "coal-py@j1", "coal-py@jdefault", "coal-cpp@j1", "coal-cpp@jdefault"]
STYLE = {"p2n14": ("#7f7f7f", "o"), "p2n16": ("#bcbd22", "s"), "coal-py@j1": ("#1f77b4", "^"),
         "coal-py@jdefault": ("#1f77b4", "v"), "coal-cpp@j1": ("#d62728", "^"), "coal-cpp@jdefault": ("#d62728", "v")}


def cfg(r):
    return r["tool"] + ("" if r["tool"].startswith("p2n") else "@j" + r["jobs"])


def med(tier, case, c, field="wall_s", allow_single=False):
    v = [float(r[field]) for r in rows if r["tier"] == tier and r["case"].endswith(case) and cfg(r) == c
         and r["cold"] == "0" and r["status"] == "OK" and r["mode"] in ("single", "single_cap_test")
         and (r["run_idx"] != "0")]
    if not v and allow_single:
        v = [float(r[field]) for r in rows if r["tier"] == tier and r["case"].endswith(case) and cfg(r) == c
             and r["status"] == "OK"]
    return st.median(v) if v else None


def sweep(cases, xs, xlabel, fname, field="wall_s", ylabel="wall time (s, median)"):
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for c in TOOLS:
        pts = [(x, med("B", case, c, field)) for x, case in zip(xs, cases)]
        pts = [(x, y) for x, y in pts if y]
        if pts:
            col, mk = STYLE[c]
            ax.plot(*zip(*pts), marker=mk, color=col, ls="-" if c.endswith(("j1", "4", "6")) else "--", label=c)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(D / fname, dpi=130); plt.close(fig)


ns = [1000, 10000, 100000, 1000000]
sweep([f"n{n}_L400" for n in ns], ns, "sequences (L=400 codons)", "tierB_nsweep_time.png")
sweep([f"n{n}_L400" for n in ns], ns, "sequences (L=400 codons)", "tierB_nsweep_mem.png", field="peak_mem_mb",
      ylabel="peak committed memory (MB, median)")
Ls = [100, 400, 1000, 10000]
sweep([f"n1000_L{L}" for L in Ls], Ls, "codons per sequence (n=1000)", "tierB_lsweep_time.png")

# Tier E: coal-cpp speed-up and efficiency
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
for case, col in (("n100000_L400", "#d62728"), ("n1000000_L400", "#9467bd")):
    js = [1, 2, 4, 8, 12, 20]
    t = {j: st.median([float(r["wall_s"]) for r in rows if r["tier"] == "E" and r["case"].endswith(case)
                       and r["tool"] == "coal-cpp" and r["jobs"] == str(j) and r["run_idx"] != "0" and r["status"] == "OK"] or [0])
         for j in js}
    if not t[1]:
        continue
    sp = [t[1] / t[j] for j in js]
    a1.plot(js, sp, marker="o", color=col, label=f"coal-cpp {case}")
    a2.plot(js, [s / j for s, j in zip(sp, js)], marker="o", color=col, label=f"coal-cpp {case}")
a1.plot([1, 20], [1, 20], color="k", lw=.6, ls=":", label="ideal")
for ax in (a1, a2):
    for x, lab in ((6, "6 P-cores"), (12, "12 P-threads"), (14, "14 cores")):
        ax.axvline(x, color="gray", lw=.7, ls="--"); ax.text(x, ax.get_ylim()[1] * .95 if ax is a2 else 17, lab,
                                                             rotation=90, fontsize=7, va="top", ha="right")
    ax.set_xlabel("-j"); ax.grid(alpha=.3); ax.legend(fontsize=8)
a1.set_ylabel("speed-up vs -j1"); a2.set_ylabel("parallel efficiency (speed-up / j)")
fig.suptitle("Tier E (reduced): coal-cpp thread scaling")
fig.tight_layout(); fig.savefig(D / "tierE_coalcpp_scaling.png", dpi=130); plt.close(fig)

# Tier C (reduced): per-file time
c = {cfg(r): [] for r in rows if r["case"] == "ncbi:sample100"}
for r in rows:
    if r["case"] == "ncbi:sample100" and r["run_idx"] != "0" and r["status"] == "OK":
        c[cfg(r)].append(float(r["wall_s"]) / 100 * 1000)
labels = [k for k in TOOLS if c.get(k)]
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.barh(labels, [st.median(c[k]) for k in labels], color=[STYLE[k][0] for k in labels])
for i, k in enumerate(labels):
    ax.text(st.median(c[k]), i, f" {st.median(c[k]):.0f} ms", va="center", fontsize=8)
ax.set_xlabel("ms per file (median of 3; one process per file; 100 NCBI families)")
fig.tight_layout(); fig.savefig(D / "tierC_ncbi100_per_file.png", dpi=130); plt.close(fig)

print("plots written:", sorted(p.name for p in D.glob("*.png")))
