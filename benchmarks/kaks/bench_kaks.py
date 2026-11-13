"""Speed of `coal kaks` against KaKs_Calculator 2.0 on one AXT file (docs/kaks.md "Benchmark").

  python benchmarks/kaks/bench_kaks.py --axt pairs.axt --coal-cpp build/coal-cpp --kaks KaKs_Calculator \
      [--methods YN,NG,LWL,LPB,MYN] [--repeats 3] [--jobs 0]
Per method: wall time and CPU time (user + sys, of the child) for KaKs_Calculator (single-threaded),
coal kaks -j 1 and coal kaks -j N (default: all hardware threads); best of --repeats; pairs per second.
"""
from __future__ import annotations

import argparse
import os
import resource
import subprocess
import tempfile
import time
from pathlib import Path


def timed(cmd) -> tuple[float, float]:
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.perf_counter()
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    wall = time.perf_counter() - t0
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    return wall, cpu


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--axt", required=True)
    ap.add_argument("--coal-cpp", required=True)
    ap.add_argument("--kaks", required=True)
    ap.add_argument("--methods", default="YN,NG,LWL,LPB,MYN")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=0)
    a = ap.parse_args()
    npairs = sum(1 for l in Path(a.axt).read_text().split("\n\n") if l.strip())
    jobs = a.jobs or os.cpu_count() or 1
    print(f"{npairs} pairs, {os.cpu_count()} hardware threads; best of {a.repeats}\n")
    print("Wall times depend on what else the machine runs; CPU time (user + sys) is the fairer single-thread "
          "comparison.\n")
    print(f"| Method | KaKs_Calculator wall / CPU s | coal -j 1 wall / CPU s | coal -j {jobs} wall s | "
          f"pairs per CPU-s KaKs_Calculator | pairs per CPU-s coal | CPU speed-up | wall speed-up -j {jobs} |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    with tempfile.TemporaryDirectory() as td:
        for m in a.methods.split(","):
            kc = min((timed([a.kaks, "-i", a.axt, "-o", f"{td}/kc", "-m", m]) for _ in range(a.repeats)),
                     key=lambda x: x[1])
            c1 = min((timed([a.coal_cpp, "kaks", "-i", a.axt, "-o", f"{td}/c1", "-m", m, "-j", "1", "-q"])
                      for _ in range(a.repeats)), key=lambda x: x[1])
            cn = min(timed([a.coal_cpp, "kaks", "-i", a.axt, "-o", f"{td}/cn", "-m", m, "-j", str(jobs), "-q"])
                     for _ in range(a.repeats))
            print(f"| {m} | {kc[0]:.2f} / {kc[1]:.2f} | {c1[0]:.2f} / {c1[1]:.2f} | {cn[0]:.2f} | "
                  f"{npairs / kc[1]:.0f} | {npairs / c1[1]:.0f} | {kc[1] / c1[1]:.1f}x | {kc[0] / cn[0]:.1f}x |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
