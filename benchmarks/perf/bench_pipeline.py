"""Where the time of the dN/dS pipeline goes, and before/after comparisons of two builds (benchmarks/perf).

Two commands:

  python benchmarks/perf/bench_pipeline.py pipeline --repeats 3 -- pairs --query A.fna --subject B.fna -j 2 -m YN
      runs `coal pipeline ...` (the arguments after --; -o is added) --repeats times and prints, per stage of
      report.json, the median wall seconds, with the run's wall, CPU (user + sys of all child processes) and peak
      RSS (the largest child)

  python benchmarks/perf/bench_pipeline.py methods --aln codon_alignments.fasta --coal-cpp OLD NEW [...] \
      --methods YN,MYN --repeats 7 [--jobs 1] [--extra --kaks-baseline]
      `coal-cpp dnds ALN --pairs consecutive -m M` per method and binary, the binaries interleaved run by run;
      prints the minimum and the median CPU seconds of each and checks that every binary writes the same bytes

The machine should be otherwise idle; on a shared one the CPU time (and its minimum over interleaved runs) is the
number to compare, not the wall time. The pipeline uses COAL_CPP / COAL_DIAMOND as usual (docs/pipeline.md).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def run(cmd: list[str]) -> tuple[float, float, float]:
    """(wall, cpu, peak RSS in MB) of cmd and its children."""
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.perf_counter()
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    wall = time.perf_counter() - t0
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    rss = after.ru_maxrss / 1024 if sys.platform != "darwin" else after.ru_maxrss / 2**20
    return wall, cpu, rss


def pipeline(a: argparse.Namespace) -> int:
    rows: dict[str, list[float]] = {}
    walls, cpus, rsss = [], [], []
    with tempfile.TemporaryDirectory() as td:
        for k in range(a.repeats):
            out = Path(td) / f"run{k}"
            w, c, r = run([sys.executable, "-m", "coal", "pipeline", *a.args, "-o", str(out), "-q"])
            walls.append(w)
            cpus.append(c)
            rsss.append(r)
            rep = json.loads((out / "report.json").read_text(encoding="utf-8"))
            for s in rep["stages"]:
                rows.setdefault(s["name"], []).append(s["seconds"])
                for key in ("msa_seconds", "align_codon_seconds", "kaks_seconds"):
                    if key in s:
                        rows.setdefault(f"  {key}", []).append(s[key])
    med = statistics.median
    print(f"coal pipeline {' '.join(a.args)}: median of {a.repeats}")
    for name, v in rows.items():
        print(f"  {name:28s} {med(v):9.3f} s")
    print(f"  {'run wall':28s} {med(walls):9.3f} s\n  {'run CPU (user+sys)':28s} {med(cpus):9.3f} s\n"
          f"  {'peak RSS':28s} {max(rsss):9.0f} MB")
    return 0


def methods(a: argparse.Namespace) -> int:
    ms = a.methods.split(",")
    fail = 0
    with tempfile.TemporaryDirectory() as td:
        for m in ms:
            cpu: dict[str, list[float]] = {b: [] for b in a.coal_cpp}
            digest: dict[str, set] = {b: set() for b in a.coal_cpp}
            for _ in range(a.repeats):
                for k, b in enumerate(a.coal_cpp):
                    out = Path(td) / f"o{k}.tsv"
                    _, c, _ = run([b, "dnds", a.aln, "--pairs", "consecutive", "-m", m, "-j", str(a.jobs),
                                   "--digits", "17", "-q", "-o", str(out), *a.extra])
                    cpu[b].append(c)
                    digest[b].add(hashlib.sha256(out.read_bytes()).hexdigest())
            same = len(set().union(*digest.values())) == 1
            fail |= not same
            for b in a.coal_cpp:
                print(f"{m:8s} {b:50s} CPU min {min(cpu[b]):7.3f} s  median {statistics.median(cpu[b]):7.3f} s")
            print(f"{m:8s} outputs {'identical' if same else 'DIFFER'}")
    return 1 if fail else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pipeline")
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("args", nargs=argparse.REMAINDER)
    q = sub.add_parser("methods")
    q.add_argument("--aln", required=True)
    q.add_argument("--coal-cpp", nargs="+", required=True)
    q.add_argument("--methods", default="YN,NG,LWL,MLWL,LPB,MLPB,MYN,Comeron,Li")
    q.add_argument("--repeats", type=int, default=5)
    q.add_argument("--jobs", type=int, default=1)
    q.add_argument("--extra", nargs=argparse.REMAINDER, default=[])
    a = ap.parse_args()
    if a.cmd == "pipeline":
        if a.args and a.args[0] == "--":
            a.args = a.args[1:]
        return pipeline(a)
    return methods(a)


if __name__ == "__main__":
    raise SystemExit(main())
