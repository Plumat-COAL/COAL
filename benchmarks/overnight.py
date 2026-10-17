"""Overnight 1M PAL2NAL runs (user instruction, 2026-10-17 ~21:35), unattended, strictly sequential:
coal-cpp drift check; p2n14 n1000000_L400 (1 run, no warm-up, 3 h timeout); p2n16 (1 run, no warm-up, 6 h timeout);
then, only if each is expected to finish before 08:00, up to 2 more p2n14 runs; drift check. 24 GB cap throughout;
OOM/TIMEOUT are results. Every completed output is graded against truth. A background thread logs the average clock
speed every 30 min (PDH counter only, no busy-loop probe, so it does not load a core). Rows go to runs.csv (tier B).
At the end: a morning report is put at the top of status.log and the results are committed locally (never pushed).
usage: overnight.py"""
from __future__ import annotations

import ctypes
import os
import subprocess
import sys
import threading
import time
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tiers  # noqa: E402
from tiers import (HARNESS_COMMIT, OUTDIR, ROOT, RUN_FIELDS, Clock, append_csv, config_name, drift_check,  # noqa: E402
                   load_manifest, log, timed_case)

SEED = 20261017
END = datetime.combine(datetime.now().date() + timedelta(days=1), datetime.strptime("08:00", "%H:%M").time())
RUNS_CSV = OUTDIR / "runs.csv"
done: list[dict] = []
drift: dict = {}


def keep_awake():
    # ES_CONTINUOUS | ES_SYSTEM_REQUIRED: no idle sleep while this process runs (does not stop Windows Update restarts)
    ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)


def clock_logger():
    c = Clock()
    while True:
        time.sleep(1800)
        log(f"overnight clock (avg over last 30 min, all cores): {c.sample()}% of nominal")


def one_run(m, tool, timeout, idx, clock, label):
    tiers.TIMEOUT = timeout
    rows = timed_case("B", (tool, None), m, 1, clock, SEED, warmup=False,
                      extra={"mode": "single_overnight"})
    for r in rows:
        r["run_idx"] = idx
    append_csv(RUNS_CSV, rows, RUN_FIELDS)
    r = rows[-1]
    log(f"overnight {label} {tool} n1000000_L400 run {idx}: {r['status']} {r.get('wall_s', '')}s "
        f"peak {r.get('peak_mem_mb', '')} MB grade={r.get('grade', '')} {r.get('reason') or ''}".rstrip())
    done.append(r)
    return r


def drift_med(when, clock):
    drift_check("overnight", clock, when, SEED)
    import csv
    rows = [r for r in csv.DictReader(open(RUNS_CSV, encoding="utf-8"))
            if r["tier"] == "overnight" and r["mode"] == f"drift_{when}" and r["run_idx"] != "0" and r["status"] == "OK"]
    w = sorted(float(r["wall_s"]) for r in rows[-5:])
    drift[when] = w[len(w) // 2] if w else None


def morning_report(started):
    lines = [f"===== MORNING REPORT (overnight 1M PAL2NAL runs, {started:%Y-%m-%d %H:%M} -> {datetime.now():%Y-%m-%d %H:%M}) ====="]
    for r in done:
        lines.append(f"  {r['tool']} run {r['run_idx']}: {r['status']}  wall {r.get('wall_s', '')} s "
                     f"({float(r['wall_s']) / 60:.1f} min)  peak committed {r.get('peak_mem_mb', '')} MB  "
                     f"grade {r.get('grade') or '-'}  {r.get('reason') or ''}".rstrip()
                     if r.get("wall_s") not in (None, "") else f"  {r['tool']} run {r['run_idx']}: {r['status']} {r.get('reason') or ''}")
    if not done:
        lines.append("  no run completed")
    d0, d1 = drift.get("start"), drift.get("end")
    lines.append(f"  drift (coal-cpp -j1 n10000_L400 median): start {d0} s, end {d1} s"
                 + (f" ({(d1 / d0 - 1) * 100:+.1f}%)" if d0 and d1 else ""))
    lines.append("  not pushed; commit is local only. SUMMARY.md / PR body still need the measured figures if Claude "
                 "was not around to update them.")
    lines.append("=" * 60)
    p = OUTDIR / "status.log"
    old = p.read_text(encoding="utf-8") if p.exists() else ""
    p.write_text("\n".join(lines) + "\n\n" + old, encoding="utf-8")


def commit():
    subprocess.run(["git", "-C", str(ROOT), "add", "benchmarks/results/" + OUTDIR.name], check=False)
    subprocess.run(["git", "-C", str(ROOT), "add", "-f", f"benchmarks/results/{OUTDIR.name}/status.log"], check=False)
    subprocess.run(["git", "-C", str(ROOT), "commit", "-q", "-m",
                    "Overnight 1M PAL2NAL runs (p2n14, p2n16): runs.csv, status.log morning report\n\n"
                    "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"], check=False)


def main():
    assert not os.environ.get("COAL_REF_SRC"), "unset COAL_REF_SRC first"
    os.environ.pop("COAL_REF_SRC", None)
    started = datetime.now()
    keep_awake()
    log(f"OVERNIGHT start: n1000000_L400 p2n14 (3 h timeout) then p2n16 (6 h timeout), 1 run each, no warm-up "
        f"(a ~40-min warm-up is not worth it), 24 GB cap; extra p2n14 runs only if they can finish before {END:%H:%M}")
    threading.Thread(target=clock_logger, daemon=True).start()
    clock = Clock()
    clock.check(force=True)
    m = next(r for r in load_manifest("sim/speed") if r["scenario"] == "n1000000_L400")
    try:
        drift_med("start", clock)
        p14 = one_run(m, "p2n14", 3 * 3600, 1, clock, "main")
        one_run(m, "p2n16", 6 * 3600, 1, clock, "main")
        if p14["status"] == "OK":
            est = float(p14["wall_s"]) * 1.15 + 120  # margin for grading
            for idx in (2, 3):
                if datetime.now() + timedelta(seconds=est) >= END:
                    log(f"overnight: extra p2n14 run {idx} not started (would not finish before {END:%H:%M})")
                    break
                one_run(m, "p2n14", 3 * 3600, idx, clock, "extra")
        else:
            log("overnight: extra p2n14 runs not started (first p2n14 run did not finish OK)")
        drift_med("end", clock)
        log(f"STEP overnight FINISHED harness_commit={HARNESS_COMMIT}")
    finally:
        morning_report(started)
        commit()


if __name__ == "__main__":
    main()
