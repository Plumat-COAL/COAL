"""Final merge for results/windows_<commit>/: correctness.csv (Tiers A, B, C, D; annotations applied) and the Tier D
timing table (tierD_timing.md). usage: wrapup_windows.py RESULTS_DIR"""
import csv
import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tiers import CORR_FIELDS  # noqa: E402

D = Path(sys.argv[1])
ann = {(r["tier"], r["case"], r["tool"]): r["reason_note"] for r in csv.DictReader(open(D / "annotations.csv", encoding="utf-8"))}
out = []
for t in ("A", "C", "D"):
    out += list(csv.DictReader(open(D / f"correctness_{t}.csv", encoding="utf-8")))
runs = list(csv.DictReader(open(D / "runs.csv", encoding="utf-8")))
# Tier B and the 1M cases: the graded warm-up run (or the first run of a single-run test) of each (tool, jobs, case)
seen = set()
for r in runs:
    if r["tier"] in ("B",) and r["cold"] == "0" and r["mode"] in ("single", "single_cap_test") and r.get("harness_commit"):
        k = (r["tool"], r["jobs"], r["case"])
        if k in seen or not (r["grade"] or r["status"] != "OK"):
            continue
        seen.add(k)
        cls = r["grade"] if r["status"] == "OK" else r["status"]
        out.append({"tool": r["tool"], "tool_commit": r["tool_commit"], "jobs": "" if r["jobs"] == "default" else r["jobs"],
                    "tier": "B", "case": r["case"], "scenario": r["case"].split("/")[-1], "class": cls,
                    "reason": r["reason"], "n_seqs": r["n_seqs"], "aln_cols": r["aln_cols"], "cds_nt": r["cds_nt"],
                    "exit_code": r["exit_code"], "wall_s": r["wall_s"], "peak_mem_mb": r["peak_mem_mb"]})
for r in out:
    note = ann.get((r["tier"], r["case"], r["tool"]))
    if note and note not in (r.get("reason") or ""):
        r["reason"] = f"{note} | {r.get('reason') or ''}"
with open(D / "correctness.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=CORR_FIELDS, extrasaction="ignore")
    w.writeheader()
    w.writerows(out)
print("correctness.csv rows:", len(out), {t: sum(r["tier"] == t for r in out) for t in "ABCD"})

# Tier D timing: per (group, config) over trees every config finished (like-for-like)
cfgname = lambda r: r["tool"] + ("" if r["tool"].startswith("p2n") else "@j" + r["jobs"])
med = defaultdict(list)
meta = {}
for r in runs:
    if r["tier"] == "D" and r["status"] == "OK" and r["run_idx"] != "0":
        g = r["mode"].split(":", 1)[1]
        med[(g, r["case"], cfgname(r))].append(float(r["wall_s"]))
        meta[r["case"]] = int(r["cds_nt"])
plan = json.load(open(D / "tierD_timing_plan.json"))
cfgs = ["p2n14", "p2n16", "coal-py@j1", "coal-py@jdefault", "coal-cpp@j1", "coal-cpp@jdefault"]
lines = ["| group | trees timed (all 6 configs) / planned | " + " | ".join(cfgs) + " |", "|---|---|" + "---|" * len(cfgs)]
for g, trees in plan.items():
    full = [c for c in trees if all(len(med.get((g, c, k), [])) == 3 for k in cfgs)]
    if not full:
        lines.append(f"| {g} | 0/{len(trees)} **not run (time)** |" + " — |" * len(cfgs))
        continue
    nt = sum(meta[c] for c in full)
    cells = []
    for k in cfgs:
        tot = sum(st.median(med[(g, c, k)]) for c in full)
        cells.append(f"{tot:.2f} s ({nt / tot / 1e6:.2f} Mnt/s)")
    tag = "" if len(full) == len(trees) else " (partial: rest not run (time))"
    lines.append(f"| {g} | {len(full)}/{len(trees)}{tag} | " + " | ".join(cells) + " |")
txt = ("Total of per-tree medians over the trees all 6 configs finished, with throughput in Mnt/s "
       "(CDS nucleotides per second).\n\n" + "\n".join(lines) + "\n")
(D / "tierD_timing.md").write_text(txt, encoding="utf-8")
print(txt)
