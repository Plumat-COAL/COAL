#!/usr/bin/env python3
"""Run a codon-alignment tool over test cases, timing it and scoring every output.

    python validation/run.py --tool coal accuracy                 # simulated edge cases
    python validation/run.py --tool coal speed --repeats 3
    python validation/run.py --tool pal2nal accuracy --keep-outputs validation/data/ref/pal2nal
    python validation/run.py --tool coal compara --limit 500
    python validation/run.py --cmd "mytool {pep} {nuc} -o {out} {opts}" accuracy

Tiers: accuracy, speed, families, compara, ncbi (any folder under validation/data/*/ with a
manifest.jsonl), or a path to a manifest.jsonl. Placeholders: {pep} {nuc} (space-separated
files) {out} {opts}. If {out} is absent, stdout is captured as the output.

--compare-to DIR also checks each output is byte-identical to DIR/<case>.fasta (e.g. PAL2NAL's).
Writes validation/results/<tool>_<tier>_<timestamp>.csv and prints a per-scenario summary.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import shlex
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
from check import check_case  # noqa: E402

TOOLS = {
    "coal": f"{shlex.quote(sys.executable)} -m coal {{pep}} {{nuc}} {{opts}} -o {{out}}",
    "pal2nal": f"perl {shlex.quote(str(REPO / 'PAL2NAL' / 'pal2nal.v14' / 'pal2nal.pl'))} {{pep}} {{nuc}} -output fasta {{opts}}",
}


def find_manifests(tier: str) -> list[Path]:
    p = Path(tier)
    if p.suffix == ".jsonl" and p.exists():
        return [p]
    hits = sorted((HERE / "data").glob(f"*/{tier}/manifest.jsonl")) + sorted((HERE / "data").glob(f"{tier}/manifest.jsonl"))
    if not hits:
        sys.exit(f"no manifest for tier {tier!r} under {HERE / 'data'}; build it first")
    return hits


def run_one(cmd_t: str, case_dir: Path, rec: dict, out: Path, err: Path, timeout: float):
    nuc = " ".join(shlex.quote(str(case_dir / f)) for f in rec["nuc"])
    cmd = cmd_t.format(pep=shlex.quote(str(case_dir / rec["pep"])), nuc=nuc, out=shlex.quote(str(out)),
                       opts=" ".join(rec.get("opts", [])))
    capture = "{out}" not in cmd_t
    env = {**os.environ, "PYTHONPATH": str(REPO / "src") + os.pathsep + os.environ.get("PYTHONPATH", "")}
    with open(err, "wb") as fe, open(out if capture else os.devnull, "wb") as fo:
        t0 = time.perf_counter()
        p = subprocess.Popen(cmd, shell=True, stdout=fo, stderr=fe, env=env, cwd=case_dir)
        peak = 0.0
        try:
            import psutil
            proc = psutil.Process(p.pid)
        except Exception:  # psutil missing or process already gone
            proc = None
        deadline = t0 + timeout
        rc = None
        while rc is None:
            rc = p.poll()
            if rc is not None:
                break
            if proc is not None:
                try:
                    tree = [proc] + proc.children(recursive=True)
                    peak = max(peak, sum(q.memory_info().rss for q in tree) / 2**20)
                except Exception:
                    pass
            if time.perf_counter() > deadline:
                p.kill()
                p.wait()
                rc = -9
                break
            time.sleep(0.005)
        if proc is None:
            peak = None
        dt = time.perf_counter() - t0
    return rc, dt, peak


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("tier")
    ap.add_argument("--tool", default="coal", help="coal | pal2nal (ignored if --cmd given)")
    ap.add_argument("--cmd", help="command template")
    ap.add_argument("--name", help="label for results (default: tool)")
    ap.add_argument("--scenario", nargs="*", help="only these scenarios")
    ap.add_argument("--limit", type=int, help="first N cases only")
    ap.add_argument("--repeats", type=int, default=1, help="timing repeats (median reported)")
    ap.add_argument("--timeout", type=float, default=3600)
    ap.add_argument("--keep-outputs", type=Path, help="save outputs here as <case>.fasta (+ .stderr)")
    ap.add_argument("--compare-to", type=Path, help="directory of reference outputs <case>.fasta")
    ap.add_argument("--no-check", action="store_true", help="time only")
    a = ap.parse_args()

    cmd_t = a.cmd or TOOLS[a.tool]
    if a.keep_outputs:
        a.keep_outputs = a.keep_outputs.resolve()
    if a.compare_to:
        a.compare_to = a.compare_to.resolve()
    name = a.name or (a.tool if not a.cmd else "custom")
    recs = []
    for m in find_manifests(a.tier):
        for line in open(m):
            r = json.loads(line)
            if a.scenario and r["scenario"] not in a.scenario:
                continue
            recs.append(r)
    if a.limit:
        recs = recs[: a.limit]
    stamp = time.strftime("%Y%m%d_%H%M%S")
    resdir = HERE / "results"
    resdir.mkdir(exist_ok=True)
    csv_path = resdir / f"{name}_{Path(a.tier).stem if a.tier.endswith('.jsonl') else a.tier}_{stamp}.csv"
    tmp = Path(tempfile.mkdtemp(prefix="coalval_"))
    rows = []
    fields = ["case", "tier", "scenario", "n_seqs", "aln_cols", "cds_nt", "rc", "seconds", "peak_mb",
              "pass", "rows_total", "rows_refree_ok", "rows_truth_exact", "col_agreement",
              "warnings_expected", "warnings_found", "identical_to_ref", "problems", "example"]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for i, rec in enumerate(recs):
            case_dir = HERE / "data" / rec["dir"]
            key = rec["dir"].replace("/", "__")
            out = (a.keep_outputs / f"{key}.fasta") if a.keep_outputs else tmp / "out.fasta"
            err = (a.keep_outputs / f"{key}.stderr") if a.keep_outputs else tmp / "err.txt"
            out.parent.mkdir(parents=True, exist_ok=True)
            times, peaks = [], []
            for _ in range(a.repeats):
                if out.exists():
                    out.unlink()
                rc, dt, peak = run_one(cmd_t, case_dir, rec, out, err, a.timeout)
                times.append(dt)
                if peak is not None:
                    peaks.append(peak)
            times.sort()
            row = {"case": rec["dir"], "tier": rec["tier"], "scenario": rec["scenario"], "n_seqs": rec["n_seqs"],
                   "aln_cols": rec["aln_cols"], "cds_nt": rec["cds_nt"], "rc": rc,
                   "seconds": round(times[len(times) // 2], 4), "peak_mb": round(max(peaks), 1) if peaks else ""}
            if not a.no_check:
                c = check_case(case_dir, out, err.read_text(errors="replace"), rec)
                row.update({k: c.get(k, "") for k in ("pass", "rows_total", "rows_refree_ok", "rows_truth_exact",
                                                     "warnings_expected", "warnings_found")})
                row["col_agreement"] = round(c["col_agreement"], 6) if "col_agreement" in c else ""
                row["problems"] = ";".join(f"{k}={v}" for k, v in c["problem_counts"].items())
                row["example"] = c["examples"][0] if c["examples"] else ""
            if a.compare_to:
                ref = a.compare_to / f"{key}.fasta"
                row["identical_to_ref"] = ref.exists() and out.exists() and ref.read_bytes().replace(b"\r\n", b"\n") == out.read_bytes().replace(b"\r\n", b"\n")
            w.writerow(row)
            fh.flush()
            rows.append(row)
            if (i + 1) % 100 == 0:
                print(f"  {i + 1}/{len(recs)}", file=sys.stderr, flush=True)
    (csv_path.with_suffix(".txt")).write_text(
        f"command: {cmd_t}\nplatform: {platform.platform()}\npython: {sys.version.split()[0]}\n"
        f"cpus: {os.cpu_count()}\nrepeats: {a.repeats}\n")
    summarise(rows, a)
    print(f"\nresults: {csv_path}")


def summarise(rows, a):
    by = defaultdict(list)
    for r in rows:
        by[r["scenario"]].append(r)
    print(f"\n{'scenario':<20}{'cases':>6}{'pass':>6}{'truth rows':>14}{'errors':>8}{'seconds':>10}{'peak MB':>9}" +
          ("  identical" if a.compare_to else ""))
    for s, rs in by.items():
        npass = sum(1 for r in rs if r.get("pass") is True)
        te = sum(int(r["rows_truth_exact"] or 0) for r in rs if r.get("rows_truth_exact") != "")
        tt = sum(int(r["rows_total"] or 0) for r in rs if r.get("rows_truth_exact") != "")
        errs = sum(1 for r in rs if r["rc"] != 0)
        secs = sum(r["seconds"] for r in rs)
        pk = max([r["peak_mb"] for r in rs if r["peak_mb"] != ""] or [0])
        line = f"{s:<20}{len(rs):>6}{npass:>6}{(f'{te}/{tt}' if tt else '-'):>14}{errs:>8}{secs:>10.2f}{pk:>9.0f}"
        if a.compare_to:
            line += f"  {sum(1 for r in rs if r['identical_to_ref'] is True)}/{len(rs)}"
        print(line)


if __name__ == "__main__":
    main()
