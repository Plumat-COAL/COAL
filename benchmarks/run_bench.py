#!/usr/bin/env python3
"""Benchmark PAL2NAL v14 vs COAL on synthetic datasets of increasing size.

For each dataset: generate it (cached in benchmarks/data/, git-ignored), run each
tool REPEATS times as a command-line process (wall-clock includes interpreter
start-up, i.e. what a user experiences), record median time, peak memory
(if psutil is installed), and whether COAL's FASTA output is byte-identical to
PAL2NAL's. Results go to benchmarks/results/<name>.csv plus a platform note.

Usage:
    python benchmarks/run_bench.py                     # default size grid
    python benchmarks/run_bench.py --quick             # tiny grid, smoke test
    python benchmarks/run_bench.py --sizes 100x500 1000x500 --repeats 5
    python benchmarks/run_bench.py --impl both --sizes 10000x500   # COAL Python vs C++ only

--impl python|cpp|both times COAL's engines only (no PAL2NAL): "coal" = Python -j 1,
"coal_cpp_j1" = C++ on one thread, "coal_cpp" = C++ on all hardware threads (its default).
The C++ binary comes from --cpp-bin (default build/coal-cpp).
"""
from __future__ import annotations

import argparse
import csv
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gen_data import generate  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "benchmarks" / "data"
RESULTS = ROOT / "benchmarks" / "results"
PAL2NAL = ROOT / "PAL2NAL" / "pal2nal.v14" / "pal2nal.pl"

# (n_seqs, n_codons). Two axes: more sequences, and longer sequences.
DEFAULT_SIZES = [(10, 500), (100, 500), (1000, 500), (10000, 500),
                 (100, 100), (100, 2000), (100, 10000)]
QUICK_SIZES = [(10, 100), (100, 300)]

try:
    import psutil
except ImportError:
    psutil = None


def run_timed(cmd, stdout_path, timeout):
    """Run cmd, stdout -> file. Return (seconds, peak_rss_MB or None, returncode)."""
    with open(stdout_path, "wb") as out:
        t0 = time.perf_counter()
        env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}  # always benchmark this checkout
        p = subprocess.Popen(cmd, stdout=out, stderr=subprocess.DEVNULL, env=env)
        peak = 0
        if psutil:
            try:
                ps = psutil.Process(p.pid)
                while p.poll() is None:
                    peak = max(peak, ps.memory_info().rss)
                    if time.perf_counter() - t0 > timeout:
                        p.kill()
                        break
                    time.sleep(0.005)
            except psutil.NoSuchProcess:
                pass
        try:
            p.wait(timeout=max(1, timeout - (time.perf_counter() - t0)))
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
        dt = time.perf_counter() - t0
    if dt >= timeout:
        return None, None, "timeout"
    return dt, (peak / 2**20 if psutil and peak else None), p.returncode


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sizes", nargs="*", help="NSEQSxNCODONS, e.g. 1000x500")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--timeout", type=float, default=900, help="seconds per run")
    ap.add_argument("--perl", default="perl")
    ap.add_argument("--name", default=None, help="results file stem")
    ap.add_argument("--coal-jobs", type=int, nargs="*", default=[],
                    help="also time COAL with these -j values (0 = all cores), e.g. --coal-jobs 2 4 0")
    ap.add_argument("--no-pal2nal", action="store_true",
                    help="time COAL only (for sizes where PAL2NAL takes too long)")
    ap.add_argument("--impl", choices=["python", "cpp", "both"],
                    help="time COAL's Python and/or C++ engine only (implies --no-pal2nal)")
    ap.add_argument("--cpp-bin", default=str(ROOT / "build" / "coal-cpp"), help="coal-cpp executable")
    ap.add_argument("--data-dir", type=Path, default=DATA,
                    help="where generated datasets are cached (default benchmarks/data)")
    ap.add_argument("--tier", choices=["A", "B", "C", "D", "E", "baseline"],
                    help="run a manifest-driven tier of the Windows release benchmark (see tiers.py); "
                         "other options are ignored")
    ap.add_argument("--workers", type=int, default=10, help="--tier: concurrent cases in correctness passes")
    ap.add_argument("--seed", type=int, default=20261017, help="--tier: seed for case/tool order and sampling")
    ap.add_argument("--cold-only", action="store_true", help="--tier B: one cold-cache run per case/tool (needs admin)")
    ap.add_argument("--include-1m", action="store_true", help="--tier B: include the n1000000_L400 case")
    ap.add_argument("--grade-only", action="store_true", help="--tier C: correctness pass only (run before timing)")
    ap.add_argument("--inproc", action="store_true",
                    help="--tier C: also time the supplementary in-process API loop (coal-py, p2n16)")
    a = ap.parse_args()
    if a.tier:
        import tiers
        return tiers.main(a)
    if a.impl:
        a.no_pal2nal = True
    cpp_bin = Path(a.cpp_bin)
    if a.impl in ("cpp", "both") and not cpp_bin.exists() and cpp_bin.with_suffix(".exe").exists():
        cpp_bin = cpp_bin.with_suffix(".exe")
    if a.impl in ("cpp", "both") and not cpp_bin.exists():
        sys.exit(f"{cpp_bin} not found; build it: cmake -S . -B build && cmake --build build")

    perl = None
    if not a.no_pal2nal:
        if not PAL2NAL.exists():
            sys.exit("PAL2NAL not found; run: python scripts/fetch_pal2nal.py")
        perl = shutil.which(a.perl) or sys.exit("perl not found")
    sizes = ([tuple(map(int, s.lower().split("x"))) for s in a.sizes] if a.sizes
             else QUICK_SIZES if a.quick else DEFAULT_SIZES)

    tools = {}
    if not a.no_pal2nal:
        tools["pal2nal"] = lambda d: [perl, str(PAL2NAL), str(d / "pep.fasta"), str(d / "nuc.fasta"), "-output", "fasta"]
    # "coal" = single process (COAL's default is automatic -j); use --coal-jobs for more
    if a.impl != "cpp":
        tools["coal"] = lambda d: [sys.executable, "-m", "coal", str(d / "pep.fasta"), str(d / "nuc.fasta"), "-j", "1"]
    if a.impl in ("cpp", "both"):
        tools["coal_cpp_j1"] = lambda d: [str(cpp_bin), str(d / "pep.fasta"), str(d / "nuc.fasta"), "-j", "1"]
        tools["coal_cpp"] = lambda d: [str(cpp_bin), str(d / "pep.fasta"), str(d / "nuc.fasta")]
    for j in a.coal_jobs:
        tools[f"coal_j{j}"] = (lambda j: lambda d: [sys.executable, "-m", "coal", str(d / "pep.fasta"),
                                                    str(d / "nuc.fasta"), "-j", str(j)])(j)

    RESULTS.mkdir(parents=True, exist_ok=True)
    stem = a.name or f"bench_{platform.system().lower()}_{datetime.now():%Y%m%d_%H%M}"
    csv_path = RESULTS / f"{stem}.csv"
    rows = []
    tmp = Path(tempfile.mkdtemp(prefix="coal_bench_"))
    try:
        for n, L in sizes:
            d = a.data_dir / f"n{n}_L{L}"
            if not (d / "nuc.fasta").exists():
                generate(n, L, d)
            nuc_mb = (d / "nuc.fasta").stat().st_size / 2**20
            res = {}
            for tool, mk in tools.items():
                times, mems, status = [], [], "ok"
                for r in range(a.repeats):
                    t, m, rc = run_timed(mk(d), tmp / f"{tool}.fa", a.timeout)
                    if t is None or rc != 0:
                        status = "timeout" if rc == "timeout" else f"rc={rc}"
                        break
                    times.append(t)
                    if m:
                        mems.append(m)
                res[tool] = (times, mems, status)
            ref_tool = "pal2nal" if "pal2nal" in res else "coal"
            outs = {t: (tmp / f"{t}.fa").read_bytes().replace(b"\r\n", b"\n") for t in res
                    if res[t][2] == "ok"}
            ref = outs.get(ref_tool)
            identical = (len(outs) == len(res) and bool(ref)
                         and all(o == ref for o in outs.values()))
            row = {"n_seqs": n, "n_codons": L, "total_codons": n * L, "nuc_MB": round(nuc_mb, 2)}
            for tool, (times, mems, status) in res.items():
                med = statistics.median(times) if times else None
                row[f"{tool}_status"] = status
                row[f"{tool}_median_s"] = round(med, 4) if med else ""
                row[f"{tool}_min_s"] = round(min(times), 4) if times else ""
                row[f"{tool}_peak_MB"] = round(max(mems), 1) if mems else ""
                row[f"{tool}_Mcodons_per_s"] = round(n * L / med / 1e6, 3) if med else ""
            py = row.get("coal_median_s")
            for tool in ("coal_cpp_j1", "coal_cpp"):
                c = row.get(f"{tool}_median_s")
                if py and c:
                    row[f"speedup_{tool}_vs_python"] = round(py / c, 2)
            p = row.get("pal2nal_median_s")
            for tool in res:
                if tool != "pal2nal":
                    c = row[f"{tool}_median_s"]
                    row[f"speedup_{tool}"] = round(p / c, 2) if p and c else ""
            row["outputs_identical"] = identical
            rows.append(row)
            times_txt = "  ".join(f"{t} {row[f'{t}_median_s'] or res[t][2]!s:>8}s" for t in res)
            sp = row.get("speedup_coal", "")
            print(f"{n:>7} x {L:<6} {times_txt}  speedup {sp!s:>6}  identical={identical}", flush=True)
            with open(csv_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    pv = subprocess.run([perl, "-e", "print $^V"], capture_output=True, text=True).stdout if perl else "n/a"
    note = (f"date: {datetime.now():%Y-%m-%d %H:%M}\nplatform: {platform.platform()}\n"
            f"cpu: {platform.processor() or platform.machine()}\npython: {platform.python_version()}\n"
            f"perl: {pv}\nrepeats: {a.repeats}\npsutil: {bool(psutil)}\ncpu_count: {os.cpu_count()}\n"
            f"coal_jobs: {a.coal_jobs}\nimpl: {a.impl}\ndata_dir: {a.data_dir}\n")
    csv_path.with_suffix(".txt").write_text(note)
    print(f"\nwrote {csv_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
