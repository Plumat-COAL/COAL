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
    a = ap.parse_args()

    if not PAL2NAL.exists():
        sys.exit("PAL2NAL not found; run: python scripts/fetch_pal2nal.py")
    perl = shutil.which(a.perl) or sys.exit("perl not found")
    sizes = ([tuple(map(int, s.lower().split("x"))) for s in a.sizes] if a.sizes
             else QUICK_SIZES if a.quick else DEFAULT_SIZES)

    tools = {
        "pal2nal": lambda d: [perl, str(PAL2NAL), str(d / "pep.fasta"), str(d / "nuc.fasta"), "-output", "fasta"],
        "coal": lambda d: [sys.executable, "-m", "coal", str(d / "pep.fasta"), str(d / "nuc.fasta")],
    }

    RESULTS.mkdir(parents=True, exist_ok=True)
    stem = a.name or f"bench_{platform.system().lower()}_{datetime.now():%Y%m%d_%H%M}"
    csv_path = RESULTS / f"{stem}.csv"
    rows = []
    tmp = Path(tempfile.mkdtemp(prefix="coal_bench_"))
    try:
        for n, L in sizes:
            d = DATA / f"n{n}_L{L}"
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
            pal_out = (tmp / "pal2nal.fa").read_bytes().replace(b"\r\n", b"\n")
            coal_out = (tmp / "coal.fa").read_bytes().replace(b"\r\n", b"\n")
            identical = (res["pal2nal"][2] == res["coal"][2] == "ok"
                         and len(pal_out) > 0 and pal_out == coal_out)
            row = {"n_seqs": n, "n_codons": L, "total_codons": n * L, "nuc_MB": round(nuc_mb, 2)}
            for tool, (times, mems, status) in res.items():
                med = statistics.median(times) if times else None
                row[f"{tool}_status"] = status
                row[f"{tool}_median_s"] = round(med, 4) if med else ""
                row[f"{tool}_min_s"] = round(min(times), 4) if times else ""
                row[f"{tool}_peak_MB"] = round(max(mems), 1) if mems else ""
                row[f"{tool}_Mcodons_per_s"] = round(n * L / med / 1e6, 3) if med else ""
            p, c = row["pal2nal_median_s"], row["coal_median_s"]
            row["speedup"] = round(p / c, 2) if p and c else ""
            row["outputs_identical"] = identical
            rows.append(row)
            print(f"{n:>6} x {L:<6} pal2nal {p or res['pal2nal'][2]!s:>9}s  coal {c or res['coal'][2]!s:>9}s  "
                  f"speedup {row['speedup']!s:>6}  identical={identical}", flush=True)
            with open(csv_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    pv = subprocess.run([perl, "-e", "print $^V"], capture_output=True, text=True).stdout
    note = (f"date: {datetime.now():%Y-%m-%d %H:%M}\nplatform: {platform.platform()}\n"
            f"cpu: {platform.processor() or platform.machine()}\npython: {platform.python_version()}\n"
            f"perl: {pv}\nrepeats: {a.repeats}\npsutil: {bool(psutil)}\n")
    csv_path.with_suffix(".txt").write_text(note)
    print(f"\nwrote {csv_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
