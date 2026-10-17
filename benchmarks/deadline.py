"""Reduced queue run under the 2026-10-17 deadline (user instruction at 19:5x): one step per invocation, strictly
sequential, 1 warm-up + 3 timed runs throughout. A watchdog stops any step at HARD_STOP (rows already written stay;
the run in flight is lost and logged). usage: deadline.py STEP   (STEP = 1m | E | C | Dcorr | Dtime)"""
from __future__ import annotations

import csv
import json
import os
import random
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tiers  # noqa: E402
from tiers import (DATA, HARNESS_COMMIT, HEADLINE, NCPU, OUTDIR, RUN_FIELDS, STATE, TOOL_COMMIT, Clock,  # noqa: E402
                   append_csv, c_loop_run, c_samples, class_table, config_name, correctness_pass, drift_check,
                   load_manifest, log, timed_case, _median_txt)

TODAY = datetime.now().date()
HARD_STOP = datetime.combine(TODAY, datetime.strptime("21:30", "%H:%M").time())  # user moved it from 21:15 (20:4x)
D_CORR_STOP = datetime.combine(TODAY, datetime.strptime("21:00", "%H:%M").time())
RUNS = 3
SEED = 20261017
RUNS_CSV = OUTDIR / "runs.csv"


def watchdog(step):
    def stop():
        while datetime.now() < HARD_STOP:
            time.sleep(5)
        log(f"STEP {step} HARD STOP at {HARD_STOP:%H:%M} (deadline): run in flight discarded, rows already written kept")
        os._exit(0)  # jobrun's Job Objects are KILL_ON_JOB_CLOSE: the tool in flight dies with this process
    threading.Thread(target=stop, daemon=True).start()


def speed_case(name):
    return next(r for r in load_manifest("sim/speed") if r["scenario"] == name)


def run_and_log(tier, cfg, m, clock, label, runs=RUNS, warmup=True, extra=None):
    rows = timed_case(tier, cfg, m, runs, clock, SEED, warmup=warmup, extra=extra)
    append_csv(RUNS_CSV, rows, RUN_FIELDS)
    log(f"{tier} {label} {config_name(*cfg)}: {_median_txt(rows) if runs > 1 else (rows[-1]['status'] + ' ' + str(rows[-1].get('wall_s', '')) + 's peak ' + str(rows[-1].get('peak_mem_mb', '')) + ' MB ' + str(rows[-1].get('reason') or ''))}")


def step_1m(clock):
    m = speed_case("n1000000_L400")
    log("1M (COAL only, no PAL2NAL): 1 warm-up + 3 timed; then coal-py default as a single run under the 24 GB cap")
    for cfg in [("coal-cpp", 1), ("coal-cpp", None), ("coal-py", 1)]:
        run_and_log("B", cfg, m, clock, "n1000000_L400")
    run_and_log("B", ("coal-py", None), m, clock, "n1000000_L400 (single run, 24 GB cap test)", runs=1,
                warmup=False, extra={"mode": "single_cap_test"})


def step_1m_p2n16(clock):
    """User addition (20:1x): one p2n16 run on 1M, no warm-up, 24 GB cap, 15-min timeout; OOM/TIMEOUT are results."""
    m = speed_case("n1000000_L400")
    tiers.TIMEOUT = 900
    log("1M p2n16: single run, no warm-up, 24 GB cap, 15-min timeout")
    run_and_log("B", ("p2n16", None), m, clock, "n1000000_L400 (single run, 24 GB cap, 900 s timeout)", runs=1,
                warmup=False, extra={"mode": "single_cap_test"})


def step_e(clock):
    log("Tier E (reduced): coal-cpp -j 1,2,4,8,12,20 on n100000_L400 and n1000000_L400; coal-py -j1/default on n100000_L400")
    for case in ("n100000_L400", "n1000000_L400"):
        m = speed_case(case)
        for j in (1, 2, 4, 8, 12, 20):
            run_and_log("E", ("coal-cpp", j), m, clock, case)
    m = speed_case("n100000_L400")
    for cfg in (("coal-py", 1), ("coal-py", None)):
        run_and_log("E", cfg, m, clock, "n100000_L400")


def step_c(clock):
    rows = load_manifest("real/ncbi")
    by_dir = {r["dir"]: r for r in rows}
    state = json.load(open(STATE / "C_ncbi.json"))
    ok = ("TRUTH", "REF", "KNOWN_DIV")
    valid = {d for d in by_dir if all(state[config_name(*c)][d]["class"] in ok for c in HEADLINE)}
    sample = c_samples(rows, valid, SEED)[100]  # same nested fixed-seed sample design as the full Tier C
    json.dump({"100": sample}, open(OUTDIR / "tierC_ncbi_samples_reduced.json", "w"))
    metas = [by_dir[d] for d in sample]
    log(f"Tier C (reduced): ncbi, one fixed-seed batch of 100 from {len(valid)} cases valid for all 6 configs "
        f"(INPUT_INVALID and p2n14-failing SELD/FDHF excluded); per-file loop; 1 warm-up + 3 timed")
    cfgs = HEADLINE[:]
    random.Random(SEED).shuffle(cfgs)
    for cfg in cfgs:
        tool, jobs = cfg
        hashes = {d: state[config_name(*cfg)][d]["hash"] for d in sample}
        base = {"tool": tool, "tool_commit": TOOL_COMMIT[tool], "tier": "C", "case": "ncbi:sample100", "seed": SEED,
                "n_seqs": sum(x["n_seqs"] for x in metas), "aln_cols": "", "cds_nt": sum(x["cds_nt"] for x in metas),
                "batch_size": len(metas), "mode": "per_file_loop", "jobs": "default" if jobs is None else jobs,
                "jobs_effective": tiers.jobs_effective(tool, jobs, max(x["cds_nt"] for x in metas)) if tool != "coal-cpp"
                else (jobs or NCPU), "cold": 0, "cold_method": "", "cpu_temp_c": "", "harness_commit": HARNESS_COMMIT}
        out = []
        for idx in range(RUNS + 1):
            clock.check()
            started = datetime.now().isoformat(timespec="milliseconds")
            res = c_loop_run("C", cfg, metas, hashes, clock)
            out.append({**base, "run_idx": idx, "started_at": started, **res, "grade": "hash_checked"})
            if res["status"] != "OK":
                out += [{**base, "run_idx": j, "status": "SKIPPED", "reason": f"after {res['status']} at run {idx}"}
                        for j in range(idx + 1, RUNS + 1)]
                break
        append_csv(RUNS_CSV, out, RUN_FIELDS)
        w = sorted(r["wall_s"] for r in out if r["run_idx"] and r["status"] == "OK")
        log(f"C ncbi batch=100 per_file_loop {config_name(*cfg)}: "
            + (f"median {w[len(w) // 2]:.2f}s ({w[len(w) // 2] / 100 * 1000:.0f} ms/file)" if w else out[0]["status"]))


def d_bins(rows):
    bins = {}
    for r in rows:
        bins.setdefault(r["scenario"], []).append(r)
    return dict(sorted(bins.items()))


def step_dcorr(clock):
    rows = load_manifest("real/compara")
    largest = {r["dir"] for r in sorted(rows, key=lambda r: -r["cds_nt"])[:10]}
    path = OUTDIR / "correctness_D.csv"
    path.unlink(missing_ok=True)
    done_rows, status = [], {}
    log(f"Tier D correctness: {len(rows)} trees x {len(HEADLINE)} configs, workers=10, 10 largest serial; "
        f"stops between chunks at 21:00")
    for b, rs in d_bins(rows).items():
        n_done = 0
        for k in range(0, len(rs), 50):
            if datetime.now() >= D_CORR_STOP:
                break
            chunk = rs[k:k + 50]
            res = correctness_pass("D", chunk, HEADLINE, workers=10, serial_cases=largest & {r["dir"] for r in chunk})
            append_csv(path, res, tiers.CORR_FIELDS)
            done_rows += res
            n_done += len(chunk)
        status[b] = (n_done, len(rs))
        log(f"D correctness {b}: {n_done}/{len(rs)} trees" + ("" if n_done == len(rs) else " (stopped at 21:00: rest not run (time))"))
    json.dump({b: {"done": d, "total": t} for b, (d, t) in status.items()}, open(OUTDIR / "tierD_correctness_coverage.json", "w"))
    if done_rows:
        (OUTDIR / "tierD_table.md").write_text(class_table(done_rows, key="scenario"))
        print(class_table(done_rows, key="scenario"))


def step_dtime(clock):
    """User update (20:4x): 20 trees per bin + the 10 largest, all 6 configs, 1 warm-up + 3; runs until HARD_STOP
    (21:30). Bins in order, the 10 largest last; what is not reached is "not run (time)"."""
    rows = load_manifest("real/compara")
    rng = random.Random(SEED)
    largest = sorted(rows, key=lambda r: -r["cds_nt"])[:10]
    big = {r["dir"] for r in largest}
    groups = [(b, rng.sample([r for r in rs if r["dir"] not in big], min(20, len(rs)))) for b, rs in d_bins(rows).items()]
    groups.append(("largest10", largest))
    plan = {g: [m["dir"] for m in ms] for g, ms in groups}
    json.dump(plan, open(OUTDIR / "tierD_timing_plan.json", "w"), indent=1)
    log(f"Tier D timing: {sum(len(v) for v in plan.values())} trees (20 per bin + the 10 largest), all 6 configs, "
        f"1 warm-up + 3; stops at {HARD_STOP:%H:%M}")
    for g, ms in groups:
        for i, m in enumerate(ms, 1):
            cfgs = HEADLINE[:]
            rng.shuffle(cfgs)
            for cfg in cfgs:
                run_and_log("D", cfg, m, clock, f"{g} {m['dir'].split('/')[-1]}", extra={"mode": f"single:{g}"})
            log(f"D timing {g}: {i}/{len(ms)} trees done")


STEPS = {"1m": step_1m, "1m_p2n16": step_1m_p2n16,"E": step_e, "C": step_c, "Dcorr": step_dcorr, "Dtime": step_dtime}

if __name__ == "__main__":
    step = sys.argv[1]
    assert not os.environ.get("COAL_REF_SRC"), "unset COAL_REF_SRC first"
    os.environ.pop("COAL_REF_SRC", None)
    if datetime.now() >= HARD_STOP:
        log(f"STEP {step}: not run (time; after 21:15)")
        sys.exit(0)
    watchdog(step)
    clock = Clock()
    clock.check(force=True)
    timed = step != "Dcorr"
    if timed:
        drift_check(step, clock, "start", SEED)
    STEPS[step](clock)
    if timed:
        drift_check(step, clock, "end", SEED)
    log(f"STEP deadline-{step} FINISHED harness_commit={HARNESS_COMMIT}")
