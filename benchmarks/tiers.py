"""Manifest-driven benchmark tiers for the Windows release benchmark (run via run_bench.py --tier X).

Tools: p2n14 (PAL2NAL v14 Perl), p2n16 (grp-bork Python port), coal-py, coal-cpp. Every run goes through
jobrun.run (Job Object: wall time, peak committed memory, CPU time, 24 GB cap). Correctness passes may run
cases concurrently; timed passes never do (see run_timed).
"""
from __future__ import annotations

import csv
import json
import os
import random
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import grade
import jobrun

ROOT = Path(__file__).resolve().parent.parent
BENCH = Path(os.environ.get("COALBENCH", r"C:\coalbench"))
DATA = BENCH / "data"
WORK = BENCH / "tmp"
COMMIT = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
OUTDIR = ROOT / "benchmarks" / "results" / f"windows_{COMMIT[:7]}"

PERL = r"C:\Strawberry\perl\bin\perl.exe"
P2N14 = ROOT / "PAL2NAL" / "pal2nal.v14" / "pal2nal.pl"
P2N16 = BENCH / "venv-p2n16" / "Scripts" / "pal2nal.exe"
COAL_PY = BENCH / "venv-coal" / "Scripts" / "python.exe"
COAL_CPP = ROOT / "build" / "coal-cpp.exe"
TOOL_COMMIT = {"p2n14": "v14 sha256:2b419417e304", "p2n16": "da572f796bb5",
               "coal-py": COMMIT[:12], "coal-cpp": COMMIT[:12]}


def tool_cmd(tool: str, jobs, meta: dict, case_dir: Path) -> list[str]:
    """jobs: None = tool default; int = -j N (COAL only)."""
    files = [str(case_dir / meta["pep"])] + [str(case_dir / n) for n in meta["nuc"]]
    opts = ["-output", "fasta"] + list(meta.get("opts", []))
    if tool == "p2n14":
        return [PERL, str(P2N14), *files, *opts]
    if tool == "p2n16":
        return [str(P2N16), *files, *opts]
    j = [] if jobs is None else ["-j", str(jobs)]
    if tool == "coal-py":
        return [str(COAL_PY), "-m", "coal", *files, *opts, *j]
    if tool == "coal-cpp":
        return [str(COAL_CPP), *files, *opts, *j]
    raise ValueError(tool)


def config_name(tool, jobs):
    return tool if jobs is None and tool.startswith("p2n") else f"{tool}@j{'default' if jobs is None else jobs}"


def load_manifest(rel: str) -> list[dict]:
    with open(DATA / rel / "manifest.jsonl") as fh:
        return [json.loads(l) for l in fh if l.strip()]


def log(msg: str):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with open(OUTDIR / "status.log", "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def append_csv(path: Path, rows: list[dict], fields: list[str]):
    new = not path.exists()
    with open(path, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerows(rows)


CORR_FIELDS = ["tool", "tool_commit", "jobs", "tier", "case", "scenario", "class", "reason", "diff_path",
               "features", "opts", "n_seqs", "aln_cols", "cds_nt", "exit_code", "wall_s", "peak_mem_mb",
               "warn_emitted", "warn_expected", "warn_set_match"]


# ------------------------------------------------------------------ correctness pass
def correctness_pass(tier: str, rows: list[dict], configs: list[tuple], workers: int = 10,
                     serial_cases: set = frozenset(), timeout: float = 7200, ref_tool="p2n14") -> list[dict]:
    """Run every (config, case) once, grade it, delete the output (diffs kept). Returns correctness rows."""
    work = WORK / tier
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)

    def one(cfg, m):
        tool, jobs = cfg
        cname = config_name(tool, jobs)
        safe = m["dir"].replace("/", "__")
        o, e = work / f"{cname}__{safe}.out", work / f"{cname}__{safe}.err"
        r = jobrun.run(tool_cmd(tool, jobs, m, DATA / m["dir"]), o, e, timeout)
        return (cname, m["dir"]), (r, o, e)

    results = {}
    par = [m for m in rows if m["dir"] not in serial_cases]
    ser = [m for m in rows if m["dir"] in serial_cases]
    # the reference first (REF grading needs it), then the others
    order = sorted(configs, key=lambda c: c[0] != ref_tool)
    with ThreadPoolExecutor(workers) as ex:
        for k, v in ex.map(lambda a: one(*a), [(c, m) for c in order for m in par]):
            results[k] = v
    for m in ser:
        for c in order:
            k, v = one(c, m)
            results[k] = v

    out_rows = []
    for m in rows:
        case_dir = DATA / m["dir"]
        ref_cfg = config_name(ref_tool, None)
        ref_recs = None
        if (ref_cfg, m["dir"]) in results and results[(ref_cfg, m["dir"])][0].status == "OK":
            ref_recs = grade.read_fasta(results[(ref_cfg, m["dir"])][1])
        exp = grade.expected_warning_set(case_dir)
        for tool, jobs in configs:
            cname = config_name(tool, jobs)
            r, o, e = results[(cname, m["dir"])]
            dpath = OUTDIR / "diffs" / tier / cname / (m["dir"].replace("/", "__") + ".diff")
            cls, reason, diff = grade.classify(tool, m, case_dir, r, o, ref_recs, dpath)
            w = grade.warnings_from_stderr(e)
            out_rows.append({
                "tool": tool, "tool_commit": TOOL_COMMIT[tool], "jobs": "" if jobs is None else jobs,
                "tier": tier, "case": m["dir"], "scenario": m.get("scenario", ""), "class": cls,
                "reason": reason, "diff_path": os.path.relpath(diff, OUTDIR) if diff else "",
                "features": ";".join(m.get("features", [])), "opts": " ".join(m.get("opts", [])),
                "n_seqs": m["n_seqs"], "aln_cols": m["aln_cols"], "cds_nt": m["cds_nt"],
                "exit_code": r.exit_code, "wall_s": round(r.wall_s, 4), "peak_mem_mb": round(r.peak_mem_mb or 0, 1),
                "warn_emitted": len(w), "warn_expected": m.get("expected_warnings", ""),
                "warn_set_match": "" if exp is None else (set(w) == exp),
            })
            if cls == "MISMATCH":  # keep stderr next to the diff for triage
                shutil.copy(e, dpath.with_suffix(".stderr"))
    shutil.rmtree(work, ignore_errors=True)
    return out_rows


def class_table(rows, key="scenario") -> str:
    classes = ["TRUTH", "REF", "KNOWN_DIV", "MISMATCH", "ERROR", "TIMEOUT", "OOM"]
    cfgs = list(dict.fromkeys(config_name(r["tool"], None if r["jobs"] == "" else r["jobs"]) for r in rows))
    groups = list(dict.fromkeys(r[key] for r in rows))
    abbr = {"TRUTH": "T", "REF": "R", "KNOWN_DIV": "K", "MISMATCH": "M", "ERROR": "E", "TIMEOUT": "TO", "OOM": "OOM"}
    def cell(rs):
        c = {k: sum(1 for r in rs if r["class"] == k) for k in classes}
        return " ".join(f"{abbr[k]}{v}" for k, v in c.items() if v) or "-"
    lines = [f"| {key} | " + " | ".join(cfgs) + " |", "|---|" + "---|" * len(cfgs)]
    for g in groups + ["**all**"]:
        rs_g = rows if g == "**all**" else [r for r in rows if r[key] == g]
        lines.append(f"| {g} | " + " | ".join(
            cell([r for r in rs_g if config_name(r["tool"], None if r["jobs"] == "" else r["jobs"]) == c])
            for c in cfgs) + " |")
    return "\n".join(lines) + "\n\nT=TRUTH R=REF K=KNOWN_DIV M=MISMATCH E=ERROR TO=TIMEOUT OOM=OOM\n"


ALL_CONFIGS = [("p2n14", None), ("p2n16", None), ("coal-py", 1), ("coal-py", None), ("coal-cpp", 1), ("coal-cpp", None)]


def tier_a(args):
    rows = load_manifest("sim/accuracy")
    log(f"Tier A start: {len(rows)} cases x {len(ALL_CONFIGS)} configs, workers={args.workers}")
    t0 = time.time()
    res = correctness_pass("A", rows, ALL_CONFIGS, workers=args.workers)
    path = OUTDIR / "correctness_A.csv"
    if path.exists():
        path.unlink()
    append_csv(path, res, CORR_FIELDS)
    table = class_table(res)
    (OUTDIR / "tierA_table.md").write_text(table)
    log(f"Tier A done in {time.time() - t0:.0f}s; wrote {path.name}, tierA_table.md")
    print(table)


TIERS = {"A": tier_a}


def main(args):
    assert "COAL_REF_SRC" not in os.environ or not os.environ["COAL_REF_SRC"], "unset COAL_REF_SRC first"
    os.environ.pop("COAL_REF_SRC", None)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    TIERS[args.tier](args)
    return 0
