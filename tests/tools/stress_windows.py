#!/usr/bin/env python3
"""Windows-side stress cases: paths with spaces, UNC paths, both engines, a closed stdout pipe, a 256 MB job limit.

    venv-win\\Scripts\\python.exe tests\\tools\\stress_windows.py --data data --cpp COAL-win\\coal-cpp.exe

Prints one line per case; the memory-limit case uses stress_winjob.py (a Job Object). Standard library + psutil.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def run(label, cmd, env, timeout=1800, stdin=None):
    t = time.monotonic()
    try:
        p = subprocess.run(cmd, env=env, capture_output=True, timeout=timeout, stdin=stdin)
        rc, err = p.returncode, p.stderr.decode("utf-8", "replace")
    except subprocess.TimeoutExpired:
        rc, err = "timeout", ""
    dt = time.monotonic() - t
    first = err.strip().splitlines()[0][:200] if err.strip() else ""
    print(f"[{label}] rc={rc} {dt:.1f}s {'TRACEBACK ' if 'Traceback' in err else ''}:: {first}", flush=True)
    return rc, err


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--cpp", required=True)
    ap.add_argument("--coal", default="coal")
    a = ap.parse_args(argv)
    data = Path(a.data).resolve()
    base_env = dict(os.environ, COAL_CPP=str(Path(a.cpp).resolve()))
    envs = {"python": dict(base_env, COAL_ENGINE="python"), "cpp": dict(base_env, COAL_ENGINE="cpp")}
    sp = data / "dir with spaces" / "sub dir"
    out = data / "out"
    out.mkdir(exist_ok=True)
    truth = data / "m100k" / "truth.fasta"
    # 1. paths with spaces, both engines, -o and stdout
    hashes = {}
    for eng, env in envs.items():
        o = sp / f"out put {eng}.fa"
        run(f"win_spaces_o_{eng}", [a.coal, str(sp / "pep.fasta"), str(sp / "nuc.fasta"), "-o", str(o)], env)
        hashes[eng] = sha(o) if o.exists() else None
        print(f"    output equals truth: {hashes[eng] == sha(truth)}")
    print(f"    engines identical: {hashes['python'] == hashes['cpp']}")
    # 2. forward slashes, relative paths, a UNC path (\\wsl$)
    os.chdir(data)
    for eng, env in envs.items():
        run(f"win_forward_slashes_{eng}", [a.coal, "m100k/pep.fasta", "m100k/nuc.fasta", "-o", f"out/fs_{eng}.fa"], env)
        print(f"    equals truth: {sha(out / f'fs_{eng}.fa') == sha(truth) if (out / f'fs_{eng}.fa').exists() else None}")
        unc = r"\\wsl$\Ubuntu\home\keanu\work\stress\m100k"
        if os.path.isdir(unc):
            run(f"win_unc_input_{eng}", [a.coal, unc + r"\pep.fasta", unc + r"\nuc.fasta", "-o", f"out/unc_{eng}.fa"], env)
            print(f"    equals truth: {sha(out / f'unc_{eng}.fa') == sha(truth) if (out / f'unc_{eng}.fa').exists() else None}")
        # a report and clustal output in the folder with spaces
        run(f"win_spaces_clustal_report_{eng}", [a.coal, str(sp / "pep.fasta"), str(sp / "nuc.fasta"), "-o", str(sp / f"out {eng}.aln"),
                                                  "-output", "clustal", "--report", str(sp / f"rep {eng}.json")], env)
    # 3. closed stdout pipe: read one line and close
    for eng, env in envs.items():
        t = time.monotonic()
        p = subprocess.Popen([a.coal, "m100k/pep.fasta", "m100k/nuc.fasta", "-output", "fasta"], env=env,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.stdout.readline()
        p.stdout.close()
        err = p.stderr.read().decode("utf-8", "replace")
        rc = p.wait()
        print(f"[win_closed_pipe_{eng}] rc={rc} {time.monotonic() - t:.1f}s {'TRACEBACK ' if 'Traceback' in err else ''}:: {err.strip()[:200]}", flush=True)
    # 4. stdin
    for eng, env in envs.items():
        with open("m100k/nuc.fasta", "rb") as fh:
            run(f"win_stdin_{eng}", [a.coal, "m100k/pep.fasta", "-", "-o", f"out/stdin_{eng}.fa"], env, stdin=fh)
        print(f"    equals truth: {sha(out / f'stdin_{eng}.fa') == sha(truth) if (out / f'stdin_{eng}.fa').exists() else None}")
    # 5. memory limit: 1 GB input in a 256 MB job (then 512 MB, 1 GB)
    if (data / "m1M" / "nuc.fasta").exists():
        for eng, env in envs.items():
            for mb in (256, 512, 1024):
                for jobs in ("", "1"):
                    cmd = [sys.executable, str(HERE / "stress_winjob.py"), "--limit-mb", str(mb), "--", a.coal,
                           "m1M/pep.fasta", "m1M/nuc.fasta", "-o", f"out/mem_{eng}_{mb}.fa"] + (["-j", jobs] if jobs else [])
                    t = time.monotonic()
                    p = subprocess.run(cmd, env=env, capture_output=True, timeout=1800)
                    err = p.stderr.decode("utf-8", "replace")
                    so = p.stdout.decode("utf-8", "replace").strip().splitlines()
                    lines = [l for l in err.strip().splitlines() if l.strip()]
                    print(f"[win_memlimit_{mb}MB{'_j' + jobs if jobs else ''}_{eng}] rc={p.returncode} {time.monotonic() - t:.1f}s "
                          f"{'TRACEBACK ' if 'Traceback' in err else ''}:: {so[-1][:160] if so else ''} :: {lines[-1][:200] if lines else ''}", flush=True)
                    o = out / f"mem_{eng}_{mb}.fa"
                    print(f"    output present: {o.exists()} ({o.stat().st_size if o.exists() else 0} bytes); leftovers: "
                          f"{[x.name for x in out.glob('.coal_*')] + [x.name for x in out.glob('mem_*.coal-*')]}")
                    if o.exists():
                        o.unlink()
    return 0


if __name__ == "__main__":
    sys.exit(main())
