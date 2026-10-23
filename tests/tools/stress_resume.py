#!/usr/bin/env python3
"""Kill a `coal -o` run at random points and resume it; the result must be byte-identical to an uninterrupted run.

    python tests/tools/stress_resume.py --dir /tmp/m1M --engine cpp --rounds 5 -- -output fasta
    python tests/tools/stress_resume.py --dir /tmp/m1M --engine python --corrupt -- -output clustal

Each round: start `coal PEP NUC... OPTS -o out`, wait a random fraction of the reference run's time, send SIGTERM or
SIGKILL (alternating), then run the same command with --resume (as many times as it takes, each time killed again at
a random point, the last one running to the end). The final output must equal the reference output. --corrupt
instead interrupts once (SIGKILL), damages the checkpoint (truncated JSON, `completed` out of range, `sink_bytes`
past the partial file, the partial file deleted, the other engine's name) and checks that `--resume` refuses it
with a notice and still produces the right output. Standard library only. Prints one line per attempt and a JSON
summary line; exit 1 if any final output differs or a traceback appears on stderr.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import random
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def run_to_end(cmd, env, cwd, log):
    t = time.monotonic()
    with open(log, "wb") as se:
        p = subprocess.run(cmd, env=env, cwd=cwd, stdout=subprocess.DEVNULL, stderr=se)
    return p.returncode, time.monotonic() - t


def kill_after(cmd, env, cwd, delay, sig, log):
    """Start the command, kill it with `sig` after `delay` s. Returns (killed?, exit code, stderr text)."""
    with open(log, "wb") as se:
        p = subprocess.Popen(cmd, env=env, cwd=cwd, stdout=subprocess.DEVNULL, stderr=se, start_new_session=True)
        try:
            rc = p.wait(timeout=delay)
            return False, rc, Path(log).read_text(errors="replace")
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, sig)
            rc = p.wait()
            time.sleep(0.2)
            return True, rc, Path(log).read_text(errors="replace")


def leftovers(out):
    return sorted(os.path.basename(p) for p in glob.glob(str(out) + ".coal-*") + glob.glob(str(out.parent / ".coal_*")))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--engine", choices=("python", "cpp"), required=True)
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--corrupt", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--log", default=str(ROOT / "tests" / "tools" / "stress_results.jsonl"))
    ap.add_argument("--name", default="")
    ap.add_argument("--min-delay", type=float, default=0.5, help="earliest kill, seconds (a checkpoint needs about 2 s)")
    ap.add_argument("args", nargs="*")
    a = ap.parse_args(argv)
    rng = random.Random(a.seed)
    d = Path(a.dir).resolve()
    pep = sorted(glob.glob(str(d / "pep.*")))[0]
    nucs = sorted(glob.glob(str(d / "nuc*.fasta")))
    env = dict(os.environ, COAL_ENGINE=a.engine)
    env.setdefault("COAL_CPP", str(ROOT / "build" / "coal-cpp"))
    tag = a.name or f"resume_{a.engine}_{'_'.join(x.strip('-') for x in a.args)}"
    out = d / f"res.{tag}.out"
    ref = d / f"res.{tag}.ref"
    base = ["coal", pep, *nucs, *a.args]
    for p in (out, ref):
        if p.exists():
            p.unlink()
    rc, t_ref = run_to_end(base + ["-o", str(ref)], env, str(d), d / f"res.{tag}.ref.err")
    if rc != 0:
        print(f"reference run failed rc={rc}", file=sys.stderr)
        return 1
    ref_sha = sha(ref)
    print(f"[{tag}] reference: {t_ref:.1f}s sha {ref_sha[:12]}", flush=True)
    results = {"name": tag, "engine": a.engine, "args": a.args, "ref_seconds": round(t_ref, 2), "rounds": [], "ok": True}
    problems = []
    for r in range(a.rounds):
        for p in [out, *[d / x for x in leftovers(out)]]:
            if p.exists():
                p.unlink()
        attempts = []
        sig = signal.SIGTERM if r % 2 == 0 else signal.SIGKILL
        cmd = base + ["-o", str(out)]
        delay = rng.uniform(a.min_delay, max(a.min_delay + 0.5, 0.7 * t_ref))
        killed, rc, err = kill_after(cmd, env, str(d), delay, sig, d / f"res.{tag}.r{r}.err")
        attempts.append({"delay": round(delay, 2), "signal": sig.name, "killed": killed, "rc": rc,
                         "leftovers": leftovers(out), "traceback": "Traceback" in err,
                         "checkpoint_after_kill": os.path.exists(str(out) + ".coal-checkpoint.json")})
        if a.corrupt:
            ck = Path(str(out) + ".coal-checkpoint.json")
            kinds = ["truncate", "completed_huge", "sink_bytes_huge", "partial_deleted", "other_engine", "not_json"]
            kind = kinds[r % len(kinds)]
            had = ck.exists()
            if had:
                st = json.loads(ck.read_text())
                if kind == "truncate":
                    ck.write_text(ck.read_text()[: len(ck.read_text()) // 2])
                elif kind == "completed_huge":
                    st["completed"] = st["nseq"] + 5
                    ck.write_text(json.dumps(st))
                elif kind == "sink_bytes_huge":
                    st["sink_bytes"] = st["sink_bytes"] * 100 + 12345
                    ck.write_text(json.dumps(st))
                elif kind == "partial_deleted":
                    for p in glob.glob(str(out) + ".coal-partial") + glob.glob(str(out) + ".coal-stage"):
                        os.unlink(p)
                elif kind == "other_engine":
                    st["engine"] = "coal-cpp" if a.engine == "python" else "coal-python"
                    ck.write_text(json.dumps(st))
                elif kind == "not_json":
                    ck.write_bytes(b"\x00\xff{{{{ not json")
            rc2, t2 = run_to_end(cmd + ["--resume"], env, str(d), d / f"res.{tag}.r{r}.resume.err")
            err2 = (d / f"res.{tag}.r{r}.resume.err").read_text(errors="replace")
            first = err2.strip().splitlines()[0][:220] if err2.strip() else ""
            same = out.exists() and sha(out) == ref_sha
            # the other-engine checkpoint must be refused with an error naming the engine (exit 1), everything else
            # with a notice and a normal run
            expect_err = kind == "other_engine"
            good = (rc2 == 1 and not out.exists()) if expect_err else (rc2 == 0 and same)
            attempts.append({"corrupt": kind, "had_checkpoint": had, "resume_rc": rc2, "seconds": round(t2, 1),
                             "same_as_reference": same, "first_stderr": first, "ok": good,
                             "traceback": "Traceback" in err2, "leftovers": leftovers(out)})
            print(f"[{tag}] round {r} {sig.name} after {delay:.1f}s (killed={killed}, checkpoint={had}) -> corrupt {kind}: "
                  f"resume rc={rc2} same={same} ok={good} :: {first[:150]}", flush=True)
            if not good or "Traceback" in err2:
                problems.append(attempts[-1])
                results["ok"] = False
            results["rounds"].append(attempts)
            continue
        # resume (possibly several times)
        n = 0
        while killed and n < 6:
            n += 1
            last = rng.random() < 0.5 or n >= 5
            delay = rng.uniform(0.5, max(1.0, 0.7 * t_ref))
            sig = rng.choice([signal.SIGTERM, signal.SIGKILL])
            if last:
                rc2, t2 = run_to_end(cmd + ["--resume"], env, str(d), d / f"res.{tag}.r{r}.{n}.err")
                err2 = (d / f"res.{tag}.r{r}.{n}.err").read_text(errors="replace")
                killed = False
                attempts.append({"resume": n, "rc": rc2, "seconds": round(t2, 1), "traceback": "Traceback" in err2,
                                 "notice": next((l[:200] for l in err2.splitlines() if "notice" in l), "")})
            else:
                killed, rc2, err2 = kill_after(cmd + ["--resume"], env, str(d), delay, sig, d / f"res.{tag}.r{r}.{n}.err")
                attempts.append({"resume": n, "delay": round(delay, 2), "signal": sig.name, "killed": killed, "rc": rc2,
                                 "traceback": "Traceback" in err2,
                                 "notice": next((l[:200] for l in err2.splitlines() if "notice" in l), "")})
        same = out.exists() and sha(out) == ref_sha
        left = leftovers(out)
        tb = any(x.get("traceback") for x in attempts)
        ok = same and not left and not tb
        print(f"[{tag}] round {r}: {len(attempts)} attempts, final same={same} leftovers={left} traceback={tb}", flush=True)
        for x in attempts:
            print(f"    {x}", flush=True)
        results["rounds"].append(attempts)
        if not ok:
            results["ok"] = False
            problems.append({"round": r, "same": same, "leftovers": left, "attempts": attempts})
    results["problems"] = problems
    with open(a.log, "a", encoding="utf-8") as fh:
        fh.write(json.dumps({"name": tag, "resume": results}) + "\n")
    print(f"[{tag}] {'OK' if results['ok'] else 'PROBLEMS'}: {len(problems)} problem(s)", flush=True)
    return 0 if results["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
