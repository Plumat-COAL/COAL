#!/usr/bin/env python3
"""Run one stress case through the Python and C++ engines (and optionally PAL2NAL v14), measure, compare, log.

    python tests/tools/stress_run.py --name long_10k --dir /tmp/c1 -- -output fasta
    python tests/tools/stress_run.py --name many_100k --dir /tmp/c2 --pal2nal --use-o -- -output clustal -nogap

For each engine (COAL_ENGINE=python / cpp; `--engines` picks) the case's `pep.*` and `nuc*.fasta` are run with the
options after `--`, under a timeout (default 1800 s; the process group is killed and the case recorded as a
timeout). Recorded per engine: wall time, peak memory of the whole process tree (psutil, resident set sizes summed
over the tree, sampled every 50 ms), exit status, stdout / -o file size and sha256, stderr line counts. Then:
  * the two engines are compared: stdout (or the -o file) bytes, stderr bytes, exit status, and the --report JSON
    (when `--report rep.{E}.json` is among the options; `{E}` is replaced by the engine name) without its
    engine-naming fields;
  * `--truth` (default: DIR/truth.fasta if it exists) compares the output bytes with the simulated truth
    (fasta output) and names the first differing record;
  * `--pal2nal` runs PAL2NAL v14 (perl, PAL2NAL env var or PAL2NAL/pal2nal.v14/pal2nal.pl) with the PAL2NAL subset
    of the options and compares its stdout with COAL's.
A JSON line per case goes to --log (default tests/tools/stress_results.jsonl), with everything above.
`--stdin nuc|pep` pipes that file to standard input (`-`), `--gzip-stdin` pipes `gzip -dc nuc.fasta.gz`.
`--env K=V` sets environment variables for the COAL runs; `--cwd` runs in another folder (temp-file tests).
Standard library + psutil.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path

try:
    import psutil
except ImportError:  # measured without memory then
    psutil = None

ROOT = Path(__file__).resolve().parent.parent.parent
VOLATILE = ("runtime_seconds", "engine", "command", "command_line", "source", "methods_text", "coal_version")
PAL2NAL_OPTS = {"-output": 1, "-codontable": 1, "-nogap": 0, "-nomismatch": 0, "-blockonly": 0}


def sha256_file(p: Path) -> tuple[str, int]:
    h = hashlib.sha256()
    n = 0
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
            n += len(chunk)
    return h.hexdigest(), n


def files_equal(a: Path, b: Path) -> tuple[bool, int]:
    """(equal, first differing byte offset or -1)."""
    off = 0
    with open(a, "rb") as fa, open(b, "rb") as fb:
        while True:
            x, y = fa.read(1 << 20), fb.read(1 << 20)
            if x != y:
                n = min(len(x), len(y))
                for i in range(n):
                    if x[i] != y[i]:
                        return False, off + i
                return False, off + n
            if not x:
                return True, -1
            off += len(x)


def fasta_records(p: Path):
    rid, chunks = None, []
    with open(p, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if rid is not None:
                    yield rid, b"".join(chunks)
                rid, chunks = line[1:].split()[0] if line[1:].split() else b"", []
            else:
                chunks.append(line.strip())
    if rid is not None:
        yield rid, b"".join(chunks)


def first_record_diff(a: Path, b: Path) -> str:
    """Describe the first FASTA record that differs between two files."""
    ia, ib = fasta_records(a), fasta_records(b)
    k = 0
    while True:
        ra, rb = next(ia, None), next(ib, None)
        if ra is None and rb is None:
            return "no record differs (files differ only in layout?)"
        if ra is None or rb is None:
            return f"record {k}: {'output' if ra is None else 'truth'} has no more records"
        if ra != rb:
            if ra[0] != rb[0]:
                return f"record {k}: ID {ra[0][:60]!r} vs truth {rb[0][:60]!r}"
            sa, sb = ra[1], rb[1]
            n = min(len(sa), len(sb))
            i = next((i for i in range(n) if sa[i] != sb[i]), n)
            return (f"record {k} ({ra[0][:40]!r}): lengths {len(sa)} vs {len(sb)}, first difference at base {i}: "
                    f"out {sa[max(0, i - 12):i + 12]!r} truth {sb[max(0, i - 12):i + 12]!r}")
        k += 1


class Measured:
    def __init__(self):
        self.rc = None
        self.wall = 0.0
        self.peak_mb = 0.0      # peak of the summed PSS over the process tree (true memory; RSS if PSS is unavailable)
        self.peak_rss_mb = 0.0  # peak of the summed RSS (counts pages shared between forked workers once per process)
        self.timeout = False
        self.cpu = 0.0
        self.max_procs = 0


def run_measured(cmd, cwd, env, stdin, stdout_path, stderr_path, timeout, sample=0.05) -> Measured:
    m = Measured()
    t0 = time.monotonic()
    with open(stdout_path, "wb") as so, open(stderr_path, "wb") as se:
        p = subprocess.Popen(cmd, cwd=cwd, env=env, stdin=stdin, stdout=so, stderr=se, start_new_session=True)
        proc = psutil.Process(p.pid) if psutil else None
        peak_pss = peak_rss = 0
        have_pss = True
        cpu = 0.0
        while True:
            try:
                rc = p.wait(timeout=sample)
                m.rc = rc
                break
            except subprocess.TimeoutExpired:
                pass
            if proc is not None:
                try:
                    procs = [proc, *proc.children(recursive=True)]
                    m.max_procs = max(m.max_procs, len(procs))
                    rss = pss = 0
                    c = 0.0
                    for q in procs:
                        try:
                            if have_pss:
                                try:
                                    mf = q.memory_full_info()
                                    rss += mf.rss
                                    pss += mf.pss
                                except (AttributeError, psutil.AccessDenied):
                                    have_pss = False
                                    rss += q.memory_info().rss
                            else:
                                rss += q.memory_info().rss
                            t = q.cpu_times()
                            c += t.user + t.system
                        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                            pass
                    peak_rss = max(peak_rss, rss)
                    peak_pss = max(peak_pss, pss)
                    cpu = max(cpu, c)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            if time.monotonic() - t0 > timeout:
                m.timeout = True
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except OSError:
                    pass
                m.rc = p.wait()
                break
        m.wall = time.monotonic() - t0
        m.peak_rss_mb = peak_rss / 1e6
        m.peak_mb = (peak_pss if have_pss and peak_pss else peak_rss) / 1e6
        m.cpu = cpu
    return m


def stderr_stats(p: Path) -> dict:
    n = w = e = no = 0
    first = []
    last = b""
    with open(p, "rb") as fh:
        for line in fh:
            n += 1
            if n <= 3:
                first.append(line.decode("utf-8", "replace").rstrip()[:200])
            if b"warning:" in line[:40]:
                w += 1
            elif b"error:" in line[:40]:
                e += 1
            elif b"notice:" in line[:40]:
                no += 1
            last = line
    return {"lines": n, "warnings": w, "errors": e, "notices": no, "first": first,
            "last": last.decode("utf-8", "replace").rstrip()[:300], "bytes": p.stat().st_size,
            "traceback": b"Traceback" in open(p, "rb").read(1 << 20)}


def strip_report(path: Path):
    try:
        obj = json.loads(path.read_bytes())
    except (OSError, ValueError) as e:
        return f"BAD {e}"
    for k in VOLATILE:
        obj.pop(k, None)
    if isinstance(obj.get("output"), dict):
        obj["output"].pop("path", None)  # the two engines write to differently named -o files here
    return obj


def plain_fasta(args) -> bool:
    """True when the options give the plain fasta output that truth.fasta describes."""
    if any(x in args for x in ("-nogap", "-nomismatch", "-blockonly", "--json")):
        return False
    if "-output" in args:
        return args[args.index("-output") + 1] == "fasta"
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--dir", required=True, help="case folder from stress_gen.py (pep.*, nuc*.fasta, truth.fasta)")
    ap.add_argument("--engines", default="python,cpp")
    ap.add_argument("--truth", help="expected output (default DIR/truth.fasta if present; 'none' to skip)")
    ap.add_argument("--pal2nal", action="store_true", help="also run PAL2NAL v14 and compare stdout")
    ap.add_argument("--use-o", action="store_true", help="write the output with -o instead of stdout")
    ap.add_argument("--stdin", choices=("nuc", "pep"), help="pipe that input to standard input ('-')")
    ap.add_argument("--gzip-stdin", action="store_true", help="pipe gzip -dc DIR/nuc.fasta.gz to standard input")
    ap.add_argument("--timeout", type=float, default=1800)
    ap.add_argument("--pal2nal-timeout", type=float, default=1800)
    ap.add_argument("--log", default=str(ROOT / "tests" / "tools" / "stress_results.jsonl"))
    ap.add_argument("--env", action="append", default=[], metavar="K=V")
    ap.add_argument("--cwd", help="run the engines in this folder (default: the case folder)")
    ap.add_argument("--note", default="", help="free text stored with the result")
    ap.add_argument("--keep-out", action="store_true", help="keep the output files (default: deleted when > 50 MB)")
    ap.add_argument("--shape", default="", help="input shape for the table, e.g. '1 x 1,000,000'")
    ap.add_argument("--relative", action="store_true", help="pass the input files as names relative to the case folder (long command lines)")
    ap.add_argument("args", nargs="*", help="coal options (after --)")
    a = ap.parse_args(argv)

    d = Path(a.dir).resolve()
    peps = sorted(glob.glob(str(d / "pep.*")))
    if not peps:
        raise SystemExit(f"{d}: no pep.* file")
    pep = peps[0]
    nucs = sorted(p for p in glob.glob(str(d / "nuc*.fasta")))
    if not nucs:
        raise SystemExit(f"{d}: no nuc*.fasta")
    truth = None
    if a.truth != "none":
        t = Path(a.truth) if a.truth else d / "truth.fasta"
        truth = t if t.exists() and (a.truth or plain_fasta(a.args)) else None
    in_bytes = sum(os.path.getsize(p) for p in [pep, *nucs])
    cwd = a.cwd or str(d)
    if a.relative:
        pep = os.path.basename(pep)
        nucs = [os.path.basename(p) for p in nucs]
    base_env = dict(os.environ)
    for kv in a.env:
        k, _, v = kv.partition("=")
        if v == "":
            base_env.pop(k, None)
        else:
            base_env[k] = v
    cpp = base_env.get("COAL_CPP") or str(ROOT / "build" / "coal-cpp")
    base_env["COAL_CPP"] = cpp
    result = {"name": a.name, "shape": a.shape, "dir": str(d), "args": a.args, "input_mb": round(in_bytes / 1e6, 2),
              "note": a.note, "engines": {}, "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    outs = {}
    for eng in a.engines.split(","):
        env = dict(base_env, COAL_ENGINE=eng)
        args = [x.replace("{E}", eng) for x in a.args]
        pep_arg, nuc_args = pep, list(nucs)
        stdin = None
        gz = None
        if a.stdin == "pep":
            stdin = open(pep, "rb")
            pep_arg = "-"
        elif a.stdin == "nuc":
            stdin = open(nucs[0], "rb")
            nuc_args = ["-", *nucs[1:]]
        elif a.gzip_stdin:
            gz = subprocess.Popen(["gzip", "-dc", str(d / "nuc.fasta.gz")], stdout=subprocess.PIPE)
            stdin = gz.stdout
            nuc_args = ["-"]
        out_path = d / f"out.{a.name}.{eng}"
        so_path = d / f"stdout.{a.name}.{eng}"
        se_path = d / f"stderr.{a.name}.{eng}"
        for p in (out_path, so_path, se_path):
            if p.exists():
                p.unlink()
        cmd = ["coal", pep_arg, *nuc_args, *args]
        if a.use_o:
            cmd += ["-o", str(out_path)]
        m = run_measured(cmd, cwd, env, stdin, so_path, se_path, a.timeout)
        if stdin is not None and gz is None:
            stdin.close()
        if gz is not None:
            gz.wait()
        produced = out_path if a.use_o else so_path
        r = {"cmd": " ".join(shlex.quote(c) for c in cmd), "rc": m.rc, "wall_s": round(m.wall, 3),
             "peak_mb": round(m.peak_mb, 1), "peak_rss_mb": round(m.peak_rss_mb, 1), "cpu_s": round(m.cpu, 2),
             "processes": m.max_procs, "timeout": m.timeout, "stderr": stderr_stats(se_path)}
        if produced.exists():
            h, n = sha256_file(produced)
            r["out_sha256"], r["out_bytes"] = h, n
        else:
            r["out_sha256"], r["out_bytes"] = None, 0
        rep = None
        for x in args:
            if x.endswith(".json") and "rep" in x:
                rp = Path(cwd) / x if not os.path.isabs(x) else Path(x)
                if rp.exists():
                    rep = strip_report(rp)
                    r["report_bytes"] = rp.stat().st_size
        r["report"] = rep
        if truth is not None and produced.exists():
            eq, off = files_equal(produced, truth)
            r["truth_equal"] = eq
            if not eq:
                r["truth_first_diff_offset"] = off
                r["truth_first_record_diff"] = first_record_diff(produced, truth)
        outs[eng] = (produced, se_path, r)
        result["engines"][eng] = r
        status = "TIMEOUT" if m.timeout else f"rc={m.rc}"
        tr = "" if truth is None else (" truth=OK" if r.get("truth_equal") else f" truth=DIFF({r.get('truth_first_record_diff', '')[:120]})")
        print(f"[{a.name}] {eng:6s} {status} {m.wall:8.2f}s {m.peak_mb:9.1f} MB (rss sum {m.peak_rss_mb:.0f}, "
              f"{m.max_procs} procs) out={r['out_bytes']:,}B stderr={r['stderr']['lines']} lines{tr}", flush=True)
        if r["stderr"]["traceback"]:
            print(f"    TRACEBACK in stderr: {r['stderr']['last'][:200]}", flush=True)
    # engine vs engine
    engs = list(outs)
    if len(engs) == 2:
        (o1, e1, r1), (o2, e2, r2) = outs[engs[0]], outs[engs[1]]
        cmp = {"rc": r1["rc"] == r2["rc"], "stdout": r1["out_sha256"] == r2["out_sha256"],
               "stderr": files_equal(e1, e2)[0] if e1.exists() and e2.exists() else None,
               "report": (r1["report"] == r2["report"]) if (r1["report"] is not None or r2["report"] is not None) else None}
        if cmp["stderr"] is False:
            # the first differing stderr line
            with open(e1, "rb") as f1, open(e2, "rb") as f2:
                for n, (x, y) in enumerate(zip(f1, f2), 1):
                    if x != y:
                        cmp["stderr_first_diff"] = {"line": n, engs[0]: x.decode("utf-8", "replace")[:300],
                                                    engs[1]: y.decode("utf-8", "replace")[:300]}
                        break
                else:
                    cmp["stderr_first_diff"] = {"line": "different lengths"}
        if cmp["report"] is False:
            cmp["report_diff"] = diff_json(r1["report"], r2["report"])
        result["engines_identical"] = cmp
        same = all(v in (True, None) for v in cmp.values() if not isinstance(v, dict))
        print(f"[{a.name}] engines identical: {'YES' if same else 'NO ' + json.dumps({k: v for k, v in cmp.items() if v is not True and v is not None})[:600]}", flush=True)
    # PAL2NAL
    if a.pal2nal:
        perl = os.environ.get("PAL2NAL") or str(ROOT / "PAL2NAL" / "pal2nal.v14" / "pal2nal.pl")
        popts = []
        it = iter(a.args)
        for x in it:
            if x in PAL2NAL_OPTS:
                popts.append(x)
                for _ in range(PAL2NAL_OPTS[x]):
                    popts.append(next(it))
            elif x.startswith("-") and x not in ("-q",) and not x.startswith("--"):
                pass
        if "-output" not in popts:
            popts += ["-output", "fasta"]
        cmd = ["perl", perl, pep, *nucs, *popts]
        so = d / f"stdout.{a.name}.pal2nal"
        se = d / f"stderr.{a.name}.pal2nal"
        m = run_measured(cmd, cwd, dict(os.environ), None, so, se, a.pal2nal_timeout)
        h, n = sha256_file(so) if so.exists() else (None, 0)
        pr = {"cmd": " ".join(shlex.quote(c) for c in cmd), "rc": m.rc, "wall_s": round(m.wall, 3),
              "peak_mb": round(m.peak_mb, 1), "timeout": m.timeout, "out_bytes": n, "out_sha256": h,
              "stderr_lines": stderr_stats(se)["lines"] if se.exists() else 0}
        for eng, (o, e, r) in outs.items():
            pr[f"same_as_{eng}"] = (h is not None and r["out_sha256"] == h and n > 0)
            if not pr[f"same_as_{eng}"] and n > 0 and o.exists():
                pr[f"diff_vs_{eng}"] = first_record_diff(o, so)
        result["pal2nal"] = pr
        print(f"[{a.name}] pal2nal rc={m.rc} {m.wall:8.2f}s {m.peak_mb:9.1f} MB out={n:,}B "
              f"same={[pr[f'same_as_{e}'] for e in outs]}", flush=True)
    with open(a.log, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(result) + "\n")
    if not a.keep_out:
        for eng, (o, e, r) in outs.items():
            if o.exists() and o.stat().st_size > 50_000_000:
                o.unlink()
    return 0


def diff_json(x, y, path="$"):
    """First difference between two JSON values, as text."""
    if type(x) is not type(y):
        return f"{path}: {str(x)[:80]!r} vs {str(y)[:80]!r}"
    if isinstance(x, dict):
        for k in sorted(set(x) | set(y)):
            if k not in x or k not in y:
                return f"{path}.{k}: only in {'second' if k not in x else 'first'}"
            r = diff_json(x[k], y[k], f"{path}.{k}")
            if r:
                return r
        return ""
    if isinstance(x, list):
        if len(x) != len(y):
            return f"{path}: list lengths {len(x)} vs {len(y)}"
        for i, (p, q) in enumerate(zip(x, y)):
            r = diff_json(p, q, f"{path}[{i}]")
            if r:
                return r
        return ""
    return "" if x == y else f"{path}: {str(x)[:120]!r} vs {str(y)[:120]!r}"


if __name__ == "__main__":
    sys.exit(main())
