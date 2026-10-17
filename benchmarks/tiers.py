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
# the frozen COAL commit under test (not HEAD: the benchmark branch adds harness/results commits on top);
# the tool sources must be unchanged from it
COMMIT = os.environ.get("COAL_FROZEN_COMMIT", "47f57f176e78d8d66eed0edc98c2b17bc3cbe6a2")
assert subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", COMMIT, "--", "src", "cpp", "CMakeLists.txt",
                       "pyproject.toml"]).returncode == 0, f"COAL tool sources differ from frozen commit {COMMIT}"
OUTDIR = ROOT / "benchmarks" / "results" / f"windows_{COMMIT[:7]}"

PERL = r"C:\Strawberry\perl\bin\perl.exe"
P2N14 = ROOT / "PAL2NAL" / "pal2nal.v14" / "pal2nal.pl"
P2N16 = BENCH / "venv-p2n16" / "Scripts" / "pal2nal.exe"
COAL_PY = BENCH / "venv-coal" / "Scripts" / "python.exe"
COAL_CPP = ROOT / "build" / "coal-cpp.exe"
HARNESS_COMMIT = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short=12", "HEAD"], capture_output=True,
                                text=True).stdout.strip() + ("" if subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--",
                                                                                  "benchmarks", ":!benchmarks/results"]).returncode == 0 else "-dirty")
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
            recs = grade.read_fasta(o) if r.status == "OK" else None
            out_rows.append({
                "tool": tool, "tool_commit": TOOL_COMMIT[tool], "jobs": "" if jobs is None else jobs,
                "tier": tier, "case": m["dir"], "scenario": m.get("scenario", ""), "class": cls,
                "reason": reason, "diff_path": os.path.relpath(diff, OUTDIR) if diff else "",
                "features": ";".join(m.get("features", [])), "opts": " ".join(m.get("opts", [])),
                "n_seqs": m["n_seqs"], "aln_cols": m["aln_cols"], "cds_nt": m["cds_nt"],
                "exit_code": r.exit_code, "wall_s": round(r.wall_s, 4), "peak_mem_mb": round(r.peak_mem_mb or 0, 1),
                "warn_emitted": len(w), "warn_expected": m.get("expected_warnings", ""),
                "warn_set_match": "" if exp is None else (set(w) == exp),
                "_hash": grade.records_hash(recs) if recs else "",
            })
            if cls == "MISMATCH":  # keep stderr next to the diff for triage
                shutil.copy(e, dpath.with_suffix(".stderr"))
    shutil.rmtree(work, ignore_errors=True)
    return out_rows


def class_table(rows, key="scenario") -> str:
    classes = ["TRUTH", "REF", "KNOWN_DIV", "MISMATCH", "ERROR", "TIMEOUT", "OOM", "INPUT_INVALID"]
    cfgs = list(dict.fromkeys(config_name(r["tool"], None if r["jobs"] == "" else r["jobs"]) for r in rows))
    groups = list(dict.fromkeys(r[key] for r in rows))
    abbr = {"TRUTH": "T", "REF": "R", "KNOWN_DIV": "K", "MISMATCH": "M", "ERROR": "E", "TIMEOUT": "TO", "OOM": "OOM",
            "INPUT_INVALID": "I"}
    def cell(rs):
        c = {k: sum(1 for r in rs if r["class"] == k) for k in classes}
        return " ".join(f"{abbr[k]}{v}" for k, v in c.items() if v) or "-"
    lines = [f"| {key} | " + " | ".join(cfgs) + " |", "|---|" + "---|" * len(cfgs)]
    for g in groups + ["**all**"]:
        rs_g = rows if g == "**all**" else [r for r in rows if r[key] == g]
        lines.append(f"| {g} | " + " | ".join(
            cell([r for r in rs_g if config_name(r["tool"], None if r["jobs"] == "" else r["jobs"]) == c])
            for c in cfgs) + " |")
    return "\n".join(lines) + "\n\nT=TRUTH R=REF K=KNOWN_DIV M=MISMATCH E=ERROR TO=TIMEOUT OOM=OOM I=INPUT_INVALID\n"


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


# ------------------------------------------------------------------ timed runs
RAMMAP = Path(os.environ.get("RAMMAP", r"C:\Users\chloe\AppData\Local\Microsoft\WinGet\Packages"
                             r"\Microsoft.Sysinternals.RAMMap_Microsoft.Winget.Source_8wekyb3d8bbwe\RAMMap64.exe"))
TIMEOUT = 7200
RUN_FIELDS = ["tool", "tool_commit", "harness_commit", "started_at", "tier", "case", "n_seqs", "aln_cols", "cds_nt", "batch_size", "mode", "jobs",
              "jobs_effective", "run_idx", "cold", "cold_method", "wall_s", "user_s", "kernel_s", "peak_mem_mb",
              "mem_source", "cpu_temp_c", "cpu_perf_pct", "exit_code", "status", "reason", "grade", "seed"]
HEADLINE = [("p2n14", None), ("p2n16", None), ("coal-py", 1), ("coal-py", None), ("coal-cpp", 1), ("coal-cpp", None)]
EXTRA_J12 = [("coal-py", 12), ("coal-cpp", 12)]
NCPU = os.cpu_count() or 1


def jobs_effective(tool, jobs, cds_nt):
    if tool.startswith("p2n"):
        return 1
    if tool == "coal-cpp":
        return jobs or NCPU  # coal-cpp: default = hardware_concurrency, independent of input size
    return max(1, min(jobs or NCPU, cds_nt // 4_000_000))  # src/coal/pipeline.py effective_jobs


class Clock:
    """% Processor Performance (actual/nominal frequency x100), averaged over the interval since last sample.
    Package temperature is not logged (user decision: no LibreHardwareMonitor driver on this desktop)."""
    def __init__(self):
        import win32pdh
        self.pdh = win32pdh
        self.q = win32pdh.OpenQuery()
        self.c = win32pdh.AddEnglishCounter(self.q, r"\Processor Information(_Total)\% Processor Performance")
        win32pdh.CollectQueryData(self.q)
        self.ref = None
        self.last_check = 0.0

    def sample(self):
        """Blank when the interval since the last sample is too short for PDH to form a rate
        ("negative denominator"); a missing clock reading must never crash a run."""
        try:
            self.pdh.CollectQueryData(self.q)
            return round(self.pdh.GetFormattedCounterValue(self.c, self.pdh.PDH_FMT_DOUBLE)[1], 1)
        except Exception:
            return ""

    def probe(self) -> float:
        """Single-thread 2 s busy loop; the clock it ran at."""
        self.sample()
        t = time.perf_counter()
        while time.perf_counter() - t < 2:
            pass
        return self.sample()

    def check(self, force=False):
        """At tier start (force) and every ~30 min: log the probe clock; pause while >10% below the start value."""
        if not force and time.time() - self.last_check < 1800:
            return
        self.last_check = time.time()
        v = self.probe()
        if self.ref is None:
            self.ref = v
            log(f"clock reference: {v}% of nominal (single-thread probe)")
            return
        waited = 0
        while v < 0.9 * self.ref and waited < 900:
            log(f"clock {v}% < 90% of start {self.ref}%: pausing 60 s")
            time.sleep(60)
            waited += 60
            v = self.probe()
        log(f"clock check: {v}% (start {self.ref}%)" + (f", paused {waited}s" if waited else ""))


def standby_mb() -> float:
    import win32pdh
    q = win32pdh.OpenQuery()
    cs = [win32pdh.AddEnglishCounter(q, rf"\Memory\Standby Cache {n} Bytes") for n in ("Core", "Normal Priority", "Reserve")]
    win32pdh.CollectQueryData(q)
    v = sum(win32pdh.GetFormattedCounterValue(c, win32pdh.PDH_FMT_LARGE)[1] for c in cs) / 2**20
    win32pdh.CloseQuery(q)
    return v


def empty_standby() -> str:
    """Empty the standby list with RAMMap and verify it emptied (RAMMap exits 0 even when it cannot)."""
    before = standby_mb()
    r = subprocess.run([str(RAMMAP), "-accepteula", "-Et"], capture_output=True)
    time.sleep(2)
    after = standby_mb()
    if r.returncode != 0 or after > 256:
        raise RuntimeError(f"RAMMap -Et did not empty the standby list (exit {r.returncode}, "
                           f"{before:.0f} -> {after:.0f} MB); needs an elevated process")
    return f"RAMMap64 -Et (standby {before:.0f}->{after:.0f} MB)"


def _sha(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def timed_case(tier, cfg, m, runs, clock, seed, warmup=True, cold=False, extra=None) -> list[dict]:
    """1 warm-up (graded; gates the rest) + `runs` timed runs, strictly sequential. cold=True: one cold run only.
    Every run's output is checked: byte-identical to the graded run, else re-graded."""
    tool, jobs = cfg
    case_dir = DATA / m["dir"]
    work = WORK / tier
    work.mkdir(parents=True, exist_ok=True)
    o, e = work / "out.fa", work / "err.txt"
    base = {"tool": tool, "tool_commit": TOOL_COMMIT[tool], "tier": tier, "case": m["dir"], "n_seqs": m["n_seqs"],
            "aln_cols": m["aln_cols"], "cds_nt": m["cds_nt"], "batch_size": 1, "mode": "single",
            "jobs": "default" if jobs is None else jobs, "jobs_effective": jobs_effective(tool, jobs, m["cds_nt"]),
            "cpu_temp_c": "", "seed": seed, "harness_commit": HARNESS_COMMIT, **(extra or {})}
    rows, good_hash = [], None
    plan = ([("cold", 0)] if cold else ([("warmup", 0)] if warmup else []) + [("timed", i) for i in range(1, runs + 1)])
    for k, (kind, idx) in enumerate(plan):
        clock.check()
        method = empty_standby() if kind == "cold" else ""
        started = datetime.now().isoformat(timespec="milliseconds")
        clock.sample()
        r = jobrun.run(tool_cmd(tool, jobs, m, case_dir), o, e, TIMEOUT)
        perf = clock.sample()
        g = ""
        if r.status == "OK":
            h = _sha(o)
            if h == good_hash:
                g = "same_as_graded"
            else:
                g, why, _ = grade.classify(tool, m, case_dir, r, o, None,
                                           OUTDIR / "diffs" / tier / config_name(tool, jobs) / (m["dir"].replace("/", "__") + ".diff"))
                if g in ("TRUTH", "REF", "KNOWN_DIV"):
                    good_hash = h
                else:
                    r.reason = why
        rows.append({**base, "run_idx": idx, "started_at": started, "cold": int(kind == "cold"), "cold_method": method,
                     "wall_s": round(r.wall_s, 4), "user_s": round(r.user_s or 0, 4), "kernel_s": round(r.kernel_s or 0, 4),
                     "peak_mem_mb": round(r.peak_mem_mb or 0, 1), "mem_source": r.mem_source, "cpu_perf_pct": perf,
                     "exit_code": r.exit_code, "status": r.status if r.status != "OK" or g in ("TRUTH", "REF", "KNOWN_DIV", "same_as_graded") else "INVALID",
                     "reason": r.reason, "grade": g})
        if rows[-1]["status"] != "OK":
            # correctness gates timing; TIMEOUT/OOM early stop: remaining runs recorded as skipped
            why = f"after {rows[-1]['status']} at run {idx}"
            for kind2, idx2 in plan[k + 1:]:
                rows.append({**base, "run_idx": idx2, "cold": int(kind2 == "cold"), "cold_method": "",
                             "status": "SKIPPED", "reason": why})
            break
    shutil.rmtree(work, ignore_errors=True)
    return rows


def drift_check(tier, clock, when, seed):
    m = next(r for r in load_manifest("sim/speed") if r["scenario"] == "n10000_L400")
    rows = timed_case(tier, ("coal-cpp", 1), m, 5, clock, seed, extra={"mode": f"drift_{when}"})
    append_csv(OUTDIR / "runs.csv", rows, RUN_FIELDS)
    log(f"drift check {tier}/{when}: coal-cpp -j1 n10000_L400 {_median_txt(rows)}")


def _median_txt(rows) -> str:
    """Summary for the log; never raises when every run failed (failures are recorded rows, not crashes)."""
    w = sorted(r["wall_s"] for r in rows if r["run_idx"] and r["status"] == "OK")
    if not w:
        return f"no valid timed runs ({rows[0]['status']} {rows[0].get('reason', '')})"
    mm = max(r["peak_mem_mb"] for r in rows if r["status"] == "OK")
    return f"median {w[len(w) // 2]:.4f}s peak {mm} MB"


def baseline(args):
    """Fixed start-up cost of each tool on a trivial 2-sequence input (sim/accuracy/pairs/fam001)."""
    clock = Clock()
    clock.check(force=True)
    m = next(r for r in load_manifest("sim/accuracy") if r["dir"] == "sim/accuracy/pairs/fam001")
    seed = args.seed
    cfgs = HEADLINE[:]
    random.Random(seed).shuffle(cfgs)
    for cfg in cfgs:
        rows = timed_case("baseline", cfg, m, 10, clock, seed, extra={"mode": "startup"})
        append_csv(OUTDIR / "runs.csv", rows, RUN_FIELDS)
        log(f"baseline {config_name(*cfg)}: {_median_txt(rows)}")


def tier_b(args):
    clock = Clock()
    clock.check(force=True)
    rows = [r for r in load_manifest("sim/speed") if args.include_1m or r["scenario"] != "n1000000_L400"]
    seed = args.seed
    rng = random.Random(seed)
    rng.shuffle(rows)
    if args.cold_only:
        log(f"Tier B cold block start (seed {seed}): {len(rows)} cases x {len(HEADLINE)} configs")
        for m in rows:
            cfgs = HEADLINE[:]
            rng.shuffle(cfgs)
            for cfg in cfgs:
                try:
                    out = timed_case("B", cfg, m, 0, clock, seed, cold=True)
                except RuntimeError as ex:  # standby list could not be emptied: no valid cold run possible
                    log(f"COLD BLOCK ABORTED at {m['scenario']} {config_name(*cfg)}: {ex}")
                    return "ABORTED"
                append_csv(OUTDIR / "runs.csv", out, RUN_FIELDS)
                log(f"cold {m['scenario']} {config_name(*cfg)}: {out[0]['status']} {out[0]['wall_s']}s {out[0].get('reason', '')}")
        log("Tier B cold block done")
        return
    log(f"Tier B start (seed {seed}): {len(rows)} cases")
    drift_check("B", clock, "start", seed)
    for m in rows:
        cfgs = HEADLINE + EXTRA_J12
        rng.shuffle(cfgs)
        for cfg in cfgs:
            n = 3 if m["cds_nt"] >= 1e9 else 5
            out = timed_case("B", cfg, m, n, clock, seed)
            append_csv(OUTDIR / "runs.csv", out, RUN_FIELDS)
            ok = sorted(r["wall_s"] for r in out if r["run_idx"] and r["status"] == "OK")
            log(f"B {m['scenario']} {config_name(*cfg)}: " + (f"median {ok[len(ok) // 2]:.3f}s" if ok else
                f"{out[0]['status']} {out[0].get('reason', '')}"))
    drift_check("B", clock, "end", seed)
    log("Tier B done")


# ------------------------------------------------------------------ Tier C: batch size
C_SETS = {"families": "sim/families", "ncbi": "real/ncbi"}
C_SIZES = [1, 10, 100, 1000, "all"]
STATE = WORK.parent / "state"   # grading hashes (outside the repo)
# coal-py: every Tier C file is < 8 Mnt, so its automatic -j resolves to 1 and "-j1" and "default" are the same
# code path; it is timed once (jobs=default, jobs_effective=1) and that row serves both headline columns.
C_CONFIGS = [("p2n14", None), ("p2n16", None), ("coal-py", None), ("coal-cpp", 1), ("coal-cpp", None)]
C_GRADE_CONFIGS = HEADLINE


# Input features whose v14 failure is a documented divergence (case meta "features" -> where it is documented).
DOCUMENTED_V14_FAILURES = {
    "selenocysteine": "documented v14 defect: p2n16 tests/divergences.tsv fix_selenocysteine; COAL README 'U (selenocysteine)'",
}


def reclassify(res, metas):
    """Cross-tool pass after grading (per-tool classify() cannot see the other tools):
    - INPUT_INVALID: case listed in input_invalid.csv (triage_inputs.py) and every config failed on it.
    - no truth and p2n14 failed: if every other config produced the same output and the case has a feature whose
      v14 failure is documented -> those configs KNOWN_DIV "v14 fails; others agree"; if they disagree -> MISMATCH.
    Rewrites class/reason in place; hashes are unchanged."""
    inval = {}
    p = OUTDIR / "input_invalid.csv"
    if p.exists():
        inval = {r["case"]: r for r in csv.DictReader(open(p, encoding="utf-8"))}
    by_case = {}
    for r in res:
        by_case.setdefault(r["case"], []).append(r)
    for case, rs in by_case.items():
        if case in inval and all(r["class"] == "ERROR" for r in rs):
            for r in rs:
                r["class"] = "INPUT_INVALID"
                r["reason"] = f"{inval[case]['category']}: {inval[case]['reason']} | tool said: {r['reason']}"
            continue
        p14 = [r for r in rs if r["tool"] == "p2n14"]
        others = [r for r in rs if r["tool"] != "p2n14"]
        if metas[case].get("truth") or not p14 or p14[0]["class"] != "ERROR" or not others:
            continue
        if any(r["class"] in ("ERROR", "TIMEOUT", "OOM") for r in others):
            continue
        agree = len({r["_hash"] for r in others}) == 1
        doc = next((DOCUMENTED_V14_FAILURES[f] for f in metas[case].get("features", []) if f in DOCUMENTED_V14_FAILURES), None)
        if agree and doc:
            if not p14[0]["reason"].startswith(doc):
                p14[0]["reason"] = f"{doc} | {p14[0]['reason']}"
            for r in others:
                r["class"], r["reason"] = "KNOWN_DIV", f"v14 fails; others agree ({doc})"
        else:
            for r in others:
                r["class"] = "MISMATCH"
                r["reason"] = "no reference (p2n14 failed); " + ("others agree, but no documented divergence" if agree
                                                                 else "other tools disagree")
    return res


def tier_c_grade(args):
    STATE.mkdir(parents=True, exist_ok=True)
    path = OUTDIR / "correctness_C.csv"
    if getattr(args, "reclassify", False):
        # re-apply reclassify() to the recorded grading (tool outputs are not rerun); hashes come from the state files
        old = list(csv.DictReader(open(path, encoding="utf-8")))
        st = {n: json.load(open(STATE / f"C_{n}.json")) for n in C_SETS}
        sets = {r["dir"]: n for n in C_SETS for r in load_manifest(C_SETS[n])}
        for r in old:
            r["_hash"] = st[sets[r["case"]]][config_name(r["tool"], None if r["jobs"] == "" else r["jobs"])][r["case"]]["hash"]
        graded = {n: [r for r in old if sets[r["case"]] == n] for n in C_SETS}
        log("Tier C reclassify: re-applying input_invalid.csv and documented-v14-failure rules to recorded grading")
    path.unlink(missing_ok=True)
    for name, rel in C_SETS.items():
        rows = load_manifest(rel)
        if getattr(args, "reclassify", False):
            res = graded[name]
        else:
            log(f"Tier C grading {name}: {len(rows)} cases x {len(C_GRADE_CONFIGS)} configs, workers={args.workers}")
            res = correctness_pass("C", rows, C_GRADE_CONFIGS, workers=args.workers)
        res = reclassify(res, {r["dir"]: r for r in rows})
        append_csv(path, res, CORR_FIELDS)
        state = {}
        for r in res:
            state.setdefault(config_name(r["tool"], None if r["jobs"] == "" else r["jobs"]), {})[r["case"]] = \
                {"class": r["class"], "hash": r["_hash"]}
        json.dump(state, open(STATE / f"C_{name}.json", "w"))
        (OUTDIR / f"tierC_{name}_table.md").write_text(class_table(res, key="scenario"))
        log(f"Tier C grading {name} done")
        print(class_table(res, key="scenario"))


def c_samples(rows, valid: set, seed: int) -> dict:
    """Nested fixed-seed samples drawn from the cases every timed config passed (like-for-like file sets)."""
    pool = sorted(r["dir"] for r in rows if r["dir"] in valid)
    random.Random(seed).shuffle(pool)
    return {k: pool[: (len(pool) if k == "all" else k)] for k in C_SIZES}


def c_loop_run(tier, cfg, metas, hashes, clock):
    """One timed per-file loop: one process per alignment, sequential. Returns aggregate row fields."""
    tool, jobs = cfg
    work = WORK / tier
    work.mkdir(parents=True, exist_ok=True)
    o, e = work / "out.fa", work / "err.txt"
    clock.sample()
    t0 = time.perf_counter()
    user = kern = peak = 0.0
    bad = []
    for m in metas:
        r = jobrun.run(tool_cmd(tool, jobs, m, DATA / m["dir"]), o, e, TIMEOUT)
        user += r.user_s or 0
        kern += r.kernel_s or 0
        peak = max(peak, r.peak_mem_mb or 0)
        recs = grade.read_fasta(o) if r.status == "OK" else None
        if r.status != "OK" or not recs or grade.records_hash(recs) != hashes[m["dir"]]:
            bad.append(f"{m['dir']}:{r.status}")
    wall = time.perf_counter() - t0
    return {"wall_s": round(wall, 4), "user_s": round(user, 4), "kernel_s": round(kern, 4), "peak_mem_mb": round(peak, 1),
            "mem_source": "job_peak_committed (max over files)", "cpu_perf_pct": clock.sample(),
            "exit_code": 0 if not bad else 1, "status": "OK" if not bad else "INVALID",
            "reason": "" if not bad else f"{len(bad)} files differ from graded output: {';'.join(bad[:5])}"}


def c_inproc_run(tier, tool, metas, hashes, clock):
    work = WORK / tier
    work.mkdir(parents=True, exist_ok=True)
    jobs_file, hash_out = work / "inproc_jobs.json", work / "inproc_hashes.json"
    calls = [[m["dir"], tool_cmd(tool, None, m, DATA / m["dir"])[3 if tool == "coal-py" else 1:]] for m in metas]
    json.dump(calls, open(jobs_file, "w"))
    py = COAL_PY if tool == "coal-py" else P2N16.parent / "python.exe"
    cmd = [str(py), str(Path(__file__).with_name("inproc_loop.py")), "coal" if tool == "coal-py" else "p2n16",
           str(jobs_file), str(hash_out)]
    clock.sample()
    r = jobrun.run(cmd, work / "out.txt", work / "err.txt", TIMEOUT)
    got = json.load(open(hash_out)) if r.status == "OK" and hash_out.exists() else {}
    bad = [m["dir"] for m in metas if got.get(m["dir"]) != hashes[m["dir"]]]
    return {"wall_s": round(r.wall_s, 4), "user_s": round(r.user_s or 0, 4), "kernel_s": round(r.kernel_s or 0, 4),
            "peak_mem_mb": round(r.peak_mem_mb or 0, 1), "mem_source": r.mem_source, "cpu_perf_pct": clock.sample(),
            "exit_code": r.exit_code, "status": r.status if r.status != "OK" else ("OK" if not bad else "INVALID"),
            "reason": r.reason or (f"{len(bad)} outputs differ: {bad[:5]}" if bad else "")}


def tier_c(args):
    if args.grade_only:
        return tier_c_grade(args)
    clock = Clock()
    clock.check(force=True)
    seed = args.seed
    rng = random.Random(seed)
    log(f"Tier C timing start (seed {seed}); runs: 1 warm-up + 5 for batch 1/10/100, 1 warm-up + 3 for 1000/all")
    drift_check("C", clock, "start", seed)
    for name, rel in C_SETS.items():
        rows = load_manifest(rel)
        by_dir = {r["dir"]: r for r in rows}
        state = json.load(open(STATE / f"C_{name}.json"))
        ok = ("TRUTH", "REF", "KNOWN_DIV")
        valid = {d for d in by_dir if all(state[config_name(*c)][d]["class"] in ok for c in C_CONFIGS)}
        log(f"Tier C {name}: {len(valid)}/{len(rows)} cases passed grading for every timed config "
            f"(excluded {len(rows) - len(valid)}, not timed for any tool)")
        samples = c_samples(rows, valid, seed)
        json.dump(samples, open(OUTDIR / f"tierC_{name}_samples.json", "w"))
        sizes = C_SIZES[:]
        rng.shuffle(sizes)
        for k in sizes:
            metas = [by_dir[d] for d in samples[k]]
            nt = sum(m["cds_nt"] for m in metas)
            n_timed = 5 if k in (1, 10, 100) else 3
            modes = [("per_file_loop", c) for c in C_CONFIGS] + ([("inproc_api_loop", ("coal-py", None)),
                                                                   ("inproc_api_loop", ("p2n16", None))] if args.inproc else [])
            rng.shuffle(modes)
            for mode, cfg in modes:
                tool, jobs = cfg
                hashes = {d: state[config_name(*cfg)][d]["hash"] for d in samples[k]}
                base = {"tool": tool, "tool_commit": TOOL_COMMIT[tool], "tier": "C", "case": f"{name}:sample{k}",
                        "n_seqs": sum(m["n_seqs"] for m in metas), "aln_cols": "", "cds_nt": nt, "batch_size": len(metas),
                        "mode": mode, "jobs": "default" if jobs is None else jobs,
                        "jobs_effective": 1 if tool != "coal-cpp" else (jobs or NCPU), "cold": 0, "cold_method": "",
                        "cpu_temp_c": "", "seed": seed, "harness_commit": HARNESS_COMMIT}
                out = []
                for idx in range(0, n_timed + 1):
                    clock.check()
                    started = datetime.now().isoformat(timespec="milliseconds")
                    res = (c_loop_run("C", cfg, metas, hashes, clock) if mode == "per_file_loop"
                           else c_inproc_run("C", tool, metas, hashes, clock))
                    out.append({**base, "run_idx": idx, "started_at": started, **res, "grade": "hash_checked"})
                    if res["status"] != "OK":
                        for j in range(idx + 1, n_timed + 1):
                            out.append({**base, "run_idx": j, "status": "SKIPPED", "reason": f"after {res['status']} at run {idx}"})
                        break
                append_csv(OUTDIR / "runs.csv", out, RUN_FIELDS)
                w = sorted(r["wall_s"] for r in out if r["run_idx"] and r["status"] == "OK")
                log(f"C {name} batch={len(metas)} {mode} {config_name(*cfg)}: " +
                    (f"median {w[len(w) // 2]:.2f}s" if w else f"{out[0]['status']} {out[0]['reason']}"))
    drift_check("C", clock, "end", seed)
    log("Tier C timing done")


TIERS = {"A": tier_a, "B": tier_b, "C": tier_c, "baseline": baseline}


def main(args):
    assert "COAL_REF_SRC" not in os.environ or not os.environ["COAL_REF_SRC"], "unset COAL_REF_SRC first"
    os.environ.pop("COAL_REF_SRC", None)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    step = f"{args.tier}{' cold-only' if getattr(args, 'cold_only', False) else ''}{' grade-only' if getattr(args, 'grade_only', False) else ''}"
    status = TIERS[args.tier](args) or "OK"
    # recorded TIMEOUT/OOM/ERROR runs are results, not failures: exit 0. Only an uncaught exception (a harness
    # crash) exits non-zero and breaks a && queue.
    log(f"STEP {step} FINISHED status={status} harness_commit={HARNESS_COMMIT}")
    return 0
