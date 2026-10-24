#!/usr/bin/env python3
"""Phase 4 of the verification: every flag and contract, both engines.

    COAL_CPP_BIN=build/coal-cpp python verify/phase4.py OUTDIR

Writes OUTDIR/results.json and prints one line per check: OK / FINDING-CANDIDATE with the reason.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adversarial import CODE1, ENV, PY, CPP, fasta, gene, genes, write  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append({"name": name, "ok": bool(ok), "detail": str(detail)[:600]})
    print(("OK   " if ok else "FINDING-CANDIDATE ") + name + ("" if ok else f": {str(detail)[:300]}"), flush=True)


def run(cmd, cwd, env=None, stdin=None, timeout=600, **kw):
    return subprocess.run(cmd, cwd=cwd, input=stdin, capture_output=True, timeout=timeout, env=env or ENV, **kw)


def engines():
    return (("python", PY), ("cpp", [CPP]))


def sha(p: Path):
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ---------------------------------------------------------------- 4.1 report
def report_checks(out: Path):
    import random
    schema = json.loads((ROOT / "docs" / "report.schema.json").read_text(encoding="utf-8"))
    try:
        import jsonschema
    except ImportError:
        jsonschema = None
    rng = random.Random(41)
    d = out / "report"
    p, n = genes(rng, 6, 30)
    j, s = n[1]; n[1] = (j, s[:15] + CODE1["W" if p[1][1][5] != "W" else "A"] + s[18:])  # one mismatch
    j, s = n[3]; n[3] = (j, s[:40] + "G" + s[40:])  # frameshift -> fails
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    good = [(i, r) for k, (i, r) in enumerate(p) if k != 3]
    write(d, "pep_ok.fasta", fasta(good)); write(d, "nuc_ok.fasta", fasta([x for k, x in enumerate(n) if k != 3]))
    cases = {
        "success": (["pep_ok.fasta", "nuc_ok.fasta", "-o", "s.fa", "--report", "s.json"], 0),
        "warnings_json": (["pep_ok.fasta", "nuc_ok.fasta", "-o", "w.fa", "--json"], 0),
        "failure": (["pep.fasta", "nuc.fasta", "-o", "f.fa", "--report", "f.json"], 1),
        "keep_going": (["pep.fasta", "nuc.fasta", "-o", "k.fa", "--report", "k.json", "--keep-going", "-j", "2"], 2),
        "keep_going_paml": (["pep.fasta", "nuc.fasta", "-o", "k.paml", "-output", "paml", "--report", "kp.json", "--keep-going"], 2),
    }
    reps = {}
    for name, (args, want_rc) in cases.items():
        for eng, cmd in engines():
            env = dict(ENV, COAL_MIN_NT_PER_JOB="0")
            r = run(cmd + args, d, env=env)
            check(f"report/{name}/{eng}: exit {want_rc}", r.returncode == want_rc, f"rc={r.returncode} {r.stderr[-300:]!r}")
            text = r.stdout if "--json" in args else (d / args[args.index("--report") + 1]).read_bytes()
            rep = json.loads(text)
            reps[(name, eng)] = rep
            if jsonschema:
                try:
                    jsonschema.validate(rep, schema)
                    check(f"report/{name}/{eng}: schema valid", True)
                except jsonschema.ValidationError as e:
                    check(f"report/{name}/{eng}: schema valid", False, e.message)
            check(f"report/{name}/{eng}: schema id", rep.get("schema") == "coal.report/1", rep.get("schema"))
            # sha256 of inputs and output
            for inp in rep["inputs"]:
                check(f"report/{name}/{eng}: sha256 of {inp['path']}", inp["sha256"] == sha(d / inp["path"]), inp["sha256"])
            o = rep.get("output")
            if o:
                check(f"report/{name}/{eng}: output sha256", o["sha256"] == sha(d / o["path"]), o)
                nrec = (d / o["path"]).read_text().count(">") if o["format"] == "fasta" else None
                if nrec is not None:
                    check(f"report/{name}/{eng}: output.sequences == records written", o["sequences"] == nrec, (o["sequences"], nrec))
            sm = rep["summary"]
            total = sm["total"]
            check(f"report/{name}/{eng}: counts add up", sm["ok"] + sm["warning"] + sm["failed"] + sm.get("not_checked", 0) == total, sm)
            check(f"report/{name}/{eng}: sequences length == total", len(rep["sequences"]) == total, len(rep["sequences"]))
            st = [x["status"] for x in rep["sequences"]]
            if name == "keep_going":
                check(f"report/{name}/{eng}: failed gene is g3 and others in input order",
                      [x["id"] for x in rep["sequences"]] == [i for i, _ in p] and st.count("failed") == 1 and rep["sequences"][3]["status"] == "failed", st)
            if name == "failure":
                check(f"report/{name}/{eng}: error recorded", rep.get("error", {}).get("code") == "inconsistent", rep.get("error"))
                check(f"report/{name}/{eng}: not_checked present and output null", "not_checked" in sm and rep["output"] is None, sm)
        a, b = reps[(name, "python")], reps[(name, "cpp")]
        for k in ("runtime_seconds", "engine"):
            a.pop(k, None); b.pop(k, None)
        check(f"report/{name}: engines identical (minus engine, runtime)", a == b, "differs")
    # every field of the schema's top level present
    top = set(schema.get("properties", {}))
    have = set(reps[("success", "python")])
    check("report: every top-level schema field appears in a success report (or is optional)",
          top - have <= {"error", "resumed_from", "not_checked"} and not (have - top), f"missing={sorted(top - have)} extra={sorted(have - top)}")


# ---------------------------------------------------------------- 4.2 --strict, -q, --max-warnings, --progress
def flag_checks(out: Path):
    import random
    rng = random.Random(42)
    d = out / "flags"
    # strict rule 1: unrelated DNA
    p, n = genes(rng, 4, 40)
    n[1] = (n[1][0], "".join(rng.choice(["GCT", "AAA", "CTG", "GAC"]) for _ in range(40)))
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    # strict rule 3: partial IDs
    write(d, "nuc_partial.fasta", fasta([(i if k != 2 else "other", s) for k, (i, s) in enumerate(n[:1] + genes(rng, 3, 40)[1])]))
    for eng, cmd in engines():
        r = run(cmd + ["pep.fasta", "nuc.fasta", "-o", "o1.fa"], d)
        check(f"strict-off/{eng}: unrelated gene written with a notice, exit 0", r.returncode == 0 and b"probably unrelated" in r.stderr and (d / "o1.fa").exists(), r.stderr[-200:])
        r = run(cmd + ["pep.fasta", "nuc.fasta", "-o", "o2.fa", "--strict"], d)
        check(f"strict/{eng}: unrelated gene -> exit 1, nothing written", r.returncode == 1 and not (d / "o2.fa").exists() and b"error" in r.stderr, (r.returncode, r.stderr[-200:]))
        r = run(cmd + ["pep.fasta", "nuc_partial.fasta", "-o", "o3.fa", "--strict"], d)
        check(f"strict/{eng}: partial ID match -> exit 1", r.returncode == 1 and not (d / "o3.fa").exists(), (r.returncode, r.stderr[-300:]))
        r2 = run(cmd + ["pep.fasta", "nuc.fasta", "-o", "o4.fa"], d)
        rq = run(cmd + ["pep.fasta", "nuc.fasta", "-o", "o5.fa", "-q"], d)
        check(f"-q/{eng}: silent stderr, same output and exit", rq.stderr == b"" and rq.returncode == r2.returncode and sha(d / "o4.fa") == sha(d / "o5.fa"), rq.stderr[:200])
        # --max-warnings
        r = run(cmd + ["pep.fasta", "nuc.fasta", "-o", "o6.fa"], d)
        nwarn = sum(1 for l in r.stderr.splitlines() if b"pepAlnPos" in l)
        more = [l for l in r.stderr.splitlines() if b"more warnings" in l]
        check(f"--max-warnings default/{eng}: 20 per-codon lines then a count line", nwarn == 20 and len(more) == 1, (nwarn, more[:1]))
        r0 = run(cmd + ["pep.fasta", "nuc.fasta", "-o", "o7.fa", "--max-warnings", "0"], d)
        n0 = sum(1 for l in r0.stderr.splitlines() if b"pepAlnPos" in l)
        check(f"--max-warnings 0/{eng}: all warnings, no count line", n0 > 20 and not any(b"more warnings" in l for l in r0.stderr.splitlines()), n0)
        r5 = run(cmd + ["pep.fasta", "nuc.fasta", "-o", "o8.fa", "--max-warnings", "5"], d)
        check(f"--max-warnings 5/{eng}", sum(1 for l in r5.stderr.splitlines() if b"pepAlnPos" in l) == 5, r5.stderr[-200:])
        rneg = run(cmd + ["pep.fasta", "nuc.fasta", "--max-warnings", "-1"], d)
        check(f"--max-warnings -1/{eng}: usage error exit 2", rneg.returncode == 2, rneg.returncode)
    # pal2nal drop-in prints all warnings (PAL2NAL compatibility)
    r = run([sys.executable, "-m", "coal_compat", "pep.fasta", "nuc.fasta"], d)
    check("pal2nal drop-in: prints every per-codon warning (no cap)", sum(1 for l in r.stderr.splitlines() if b"pepAlnPos" in l) > 20 and not any(b"more warnings" in l for l in r.stderr.splitlines()), r.stderr[-200:])
    # --progress: small input -> no progress lines; large input -> lines on stderr only
    p, n = genes(rng, 50, 20)
    write(d, "small_p.fasta", fasta(p)); write(d, "small_n.fasta", fasta(n))
    for eng, cmd in engines():
        r = run(cmd + ["small_p.fasta", "small_n.fasta", "--progress"], d)
        check(f"--progress small/{eng}: no progress lines, stdout clean", b"progress" not in r.stderr and b"progress" not in r.stdout and r.returncode == 0, r.stderr[:200])
    big_p, big_n = genes(rng, 12_000, 20)
    write(d, "big_p.fasta", fasta(big_p)); write(d, "big_n.fasta", fasta(big_n))
    for eng, cmd in engines():
        env = dict(ENV, COAL_PROGRESS_EVERY="0")
        r = run(cmd + ["big_p.fasta", "big_n.fasta", "--progress", "-o", f"big.{eng}.fa"], d, env=env)
        check(f"--progress 12,000 seqs/{eng}: progress lines on stderr, none on stdout, exit 0",
              b"coal: progress:" in r.stderr and r.stdout == b"" and r.returncode == 0, r.stderr[-200:])
        rq = run(cmd + ["big_p.fasta", "big_n.fasta", "--progress", "-q", "-o", f"bigq.{eng}.fa"], d, env=env)
        check(f"-q --progress/{eng}: progress still printed", b"coal: progress:" in rq.stderr, rq.stderr[-200:])


# ---------------------------------------------------------------- 4.3 --resume
def resume_checks(out: Path):
    import random
    rng = random.Random(43)
    d = out / "resume"
    p, n = genes(rng, 20_000, 15)
    # a few warnings in the mix (mismatches) so the report has entries
    for k in range(0, 20_000, 997):
        j, s = n[k]; n[k] = (j, s[:9] + CODE1["W" if p[k][1][3] != "W" else "A"] + s[12:])
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    hold = dict(ENV, COAL_TEST_HOOKS="1", COAL_CHECKPOINT_EVERY="0", COAL_BATCH_NT="30000", COAL_MIN_NT_PER_JOB="0",
                COAL_TEST_PAUSE_AFTER_CHECKPOINT="4")
    for fmt in ("fasta", "clustal", "paml", "codon"):
        for eng, cmd in engines():
            ext = {"fasta": "fa", "clustal": "aln", "paml": "paml", "codon": "codon"}[fmt]
            ref = d / f"ref.{eng}.{ext}"
            args = ["pep.fasta", "nuc.fasta", "-output", fmt, "-o", f"int.{eng}.{ext}", "--report", f"int.{eng}.{ext}.json", "-j", "2"]
            r = run(cmd + ["pep.fasta", "nuc.fasta", "-output", fmt, "-o", ref.name, "--report", ref.name + ".json", "-j", "2"], d, env=dict(ENV, COAL_MIN_NT_PER_JOB="0"))
            check(f"resume/{fmt}/{eng}: reference run ok", r.returncode == 0, r.stderr[-200:])
            k = subprocess.Popen(cmd + args, cwd=d, env=hold, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            paused = False
            for _ in range(2000):
                line = k.stderr.readline()
                if not line:
                    break
                if b"paused after checkpoint 4" in line:
                    paused = True
                    break
            k.kill(); k.communicate(timeout=60)
            check(f"resume/{fmt}/{eng}: run paused at checkpoint 4 then killed", paused, "no pause marker")
            r = run(cmd + args + ["--resume"], d, env=dict(ENV, COAL_MIN_NT_PER_JOB="0"))
            check(f"resume/{fmt}/{eng}: resumed run exit 0 with the resuming notice", r.returncode == 0 and b"resuming" in r.stderr, r.stderr[:300])
            check(f"resume/{fmt}/{eng}: output byte-identical to the uninterrupted run", (d / args[args.index("-o") + 1]).read_bytes() == ref.read_bytes())
            a = json.loads((d / (args[args.index("-o") + 1] + ".json")).read_bytes()); b = json.loads((ref.with_name(ref.name + ".json")).read_bytes())
            for x in (a, b):
                for key in ("runtime_seconds", "command_line", "command", "output", "resumed_from"):
                    x.pop(key, None)
                x["methods_text"] = x.get("methods_text", "").split(" Command: ")[0]
            check(f"resume/{fmt}/{eng}: report identical (minus runtime, command, output path, resumed_from)", a == b)
            left = sorted(x.name for x in d.iterdir() if x.name.startswith((".coal_", f"int.{eng}.{ext}.coal")))
            check(f"resume/{fmt}/{eng}: no checkpoint or temp files left", left == [], left)
        # inputs changed -> refusal (Python engine only needs one format)
    for eng, cmd in engines():
        k = subprocess.Popen(cmd + ["pep.fasta", "nuc.fasta", "-o", f"chg.{eng}.fa"], cwd=d, env=hold, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        for _ in range(2000):
            line = k.stderr.readline()
            if not line or b"paused after checkpoint 4" in line:
                break
        k.kill(); k.communicate(timeout=60)
        with open(d / "nuc.fasta", "a") as fh:
            fh.write("\n")
        r = run(cmd + ["pep.fasta", "nuc.fasta", "-o", f"chg.{eng}.fa", "--resume"], d)
        check(f"resume changed inputs/{eng}: refused with the documented message, exit 1",
              r.returncode == 1 and b"changed since the interrupted run" in r.stderr, (r.returncode, r.stderr[:300]))
        for x in d.glob(f"chg.{eng}.fa*"):
            x.unlink()
    for x in list(d.glob("ref.*")) + list(d.glob("int.*")):
        x.unlink()


# ---------------------------------------------------------------- 4.4 signals and exit codes
def signal_checks(out: Path):
    import random
    rng = random.Random(44)
    d = out / "signals"
    p, n = genes(rng, 20_000, 15)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    hold = dict(ENV, COAL_TEST_HOOKS="1", COAL_CHECKPOINT_EVERY="0", COAL_BATCH_NT="30000", COAL_MIN_NT_PER_JOB="0",
                COAL_TEST_PAUSE_AFTER_CHECKPOINT="2")
    for eng, cmd in engines():
        for sig, want in ((signal.SIGINT, 130), (signal.SIGTERM, 143)):
            k = subprocess.Popen(cmd + ["pep.fasta", "nuc.fasta", "-o", f"sig.{eng}.fa"], cwd=d, env=hold, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, start_new_session=True)
            for _ in range(2000):
                line = k.stderr.readline()
                if not line or b"paused after checkpoint 2" in line:
                    break
            os.kill(k.pid, sig)
            try:
                _, err = k.communicate(timeout=60)
            except subprocess.TimeoutExpired:
                k.kill(); _, err = k.communicate()
            check(f"signal {sig.name}/{eng}: exit {want}", k.returncode == want, (k.returncode, err[-200:]))
            # the checkpoint files are kept for --resume (a checkpoint was written); the lock must go
            left = sorted(x.name for x in d.iterdir() if x.name.startswith((".coal_", f"sig.{eng}.fa")))
            check(f"signal {sig.name}/{eng}: only resumable files left (no lock, no .coal_ temp)", all(("coal-partial" in x or "coal-checkpoint" in x or "coal-entries" in x or "coal-stage" in x) for x in left), left)
            for x in d.glob(f"sig.{eng}.fa*"):
                x.unlink()
        # closed pipe -> 141, silent
        r = subprocess.run(f"{' '.join(repr(c) for c in cmd)} pep.fasta nuc.fasta | head -c 100 > /dev/null; echo ${{PIPESTATUS[0]}}", shell=True, cwd=d, env=ENV, capture_output=True, timeout=600, executable="/bin/bash")
        check(f"closed pipe/{eng}: exit 141 and silent stderr", r.stdout.strip() == b"141" and r.stderr == b"", (r.stdout, r.stderr[:200]))
        left = sorted(x.name for x in d.iterdir() if x.name.startswith(".coal_"))
        check(f"closed pipe/{eng}: no temp files", left == [], left)
    # temp folder leftovers
    tmp = Path(os.environ.get("TMPDIR", "/tmp"))
    left = sorted(x.name for x in tmp.iterdir() if x.name.startswith(".coal_"))
    check("no .coal_ files left in the temp folder after the phase", left == [], left[:10])


# ---------------------------------------------------------------- 4.5 pairing, hint, methods text
def pairing_hint_methods(out: Path):
    import random
    rng = random.Random(45)
    d = out / "pairing"
    p, n = genes(rng, 4, 30)
    # normalised: AT1G01.1 vs AT1G01, DNA reversed
    pn = [(f"AT1G0{k}.1", r) for k, (i, r) in enumerate(p)]
    nn = [(f"AT1G0{k}", s) for k, (i, s) in enumerate(n)][::-1]
    write(d, "pep_norm.fasta", fasta(pn)); write(d, "nuc_norm.fasta", fasta(nn))
    # ambiguous: two DNA records normalise to the same protein ID -> must not pair by normalised IDs
    na = [("AT1G00", n[0][1]), ("AT1G00.2", n[1][1]), ("AT1G02", n[2][1]), ("AT1G03", n[3][1])]
    write(d, "nuc_amb.fasta", fasta(na))
    # order fallback: no ID matches
    no = [(f"x{k}", s) for k, (i, s) in enumerate(n)]
    write(d, "nuc_order.fasta", fasta(no))
    for eng, cmd in engines():
        r = run(cmd + ["pep_norm.fasta", "nuc_norm.fasta", "-o", f"norm.{eng}.fa", "--report", f"norm.{eng}.json"], d)
        rep = json.loads((d / f"norm.{eng}.json").read_bytes())
        check(f"pairing normalised/{eng}: exit 0, pairing 'normalised', notice printed, every gene ok",
              r.returncode == 0 and rep["summary"]["pairing"] == "normalised" and rep["summary"]["ok"] == 4 and b"normalis" in r.stderr, (r.returncode, rep["summary"], r.stderr[:200]))
        r = run(cmd + ["pep_norm.fasta", "nuc_amb.fasta", "-o", f"amb.{eng}.fa", "--report", f"amb.{eng}.json"], d)
        rep = json.loads((d / f"amb.{eng}.json").read_bytes())
        check(f"pairing ambiguous/{eng}: no normalised pairing when two DNA IDs tidy to one protein ID (pairing != normalised, or every pair verified)",
              rep["summary"]["pairing"] != "normalised" or rep["summary"]["ok"] == 4, (rep["summary"], r.stderr[:300]))
        r = run(cmd + ["pep_norm.fasta", "nuc_order.fasta", "-o", f"ord.{eng}.fa", "--report", f"ord.{eng}.json"], d)
        rep = json.loads((d / f"ord.{eng}.json").read_bytes())
        check(f"pairing by order/{eng}: pairing 'order', id_mismatch hint in the report", rep["summary"]["pairing"] == "order" and any(h["code"] == "id_mismatch" for h in rep["hints"]), (rep["summary"], [h["code"] for h in rep["hints"]]))
    # the mito example: code 2 suggested
    ex = ROOT / "examples" / "mito"
    pep = next(ex.glob("*.aln*")) if list(ex.glob("*.aln*")) else next(ex.glob("protein*"))
    dna = next(ex.glob("dna*"))
    for eng, cmd in engines():
        r = run(cmd + [str(pep), str(dna), "-o", str(d / f"mito.{eng}.fa"), "--report", str(d / f"mito.{eng}.json")], d)
        rep = json.loads((d / f"mito.{eng}.json").read_bytes())
        gc = [h for h in rep["hints"] if h["code"] == "genetic_code"]
        check(f"mito hint/{eng}: suggests code 2", gc and gc[0]["suggested_codontable"] == 2 and gc[0]["fit"] > gc[0]["fit_current"], gc)
        check(f"mito hint/{eng}: early notice on stderr names code 2", b"Re-run with -codontable 2" in r.stderr, r.stderr[:300])
        r2 = run(cmd + [str(pep), str(dna), "-codontable", "2", "-o", str(d / f"mito2.{eng}.fa"), "--report", str(d / f"mito2.{eng}.json")], d)
        rep2 = json.loads((d / f"mito2.{eng}.json").read_bytes())
        check(f"mito with code 2/{eng}: all ok, no hint", rep2["summary"]["ok"] == rep2["summary"]["total"] and not [h for h in rep2["hints"] if h["code"] == "genetic_code"], rep2["summary"])
        m = rep2["methods_text"]
        from coal import __version__
        check(f"methods_text/{eng}: version, URL, code name, counts", f"COAL v{__version__}" in m and "github.com/plumat-coal/COAL" in m and "vertebrate mitochondrial genetic code (NCBI translation table 2)" in m and f"Of {rep2['summary']['total']} sequences, {rep2['summary']['ok']} matched" in m, m[:300])
    g = subprocess.run(["grep", "-rIl", "chloe-newt", str(ROOT), "--exclude-dir=.git", "--exclude-dir=.venv", "--exclude-dir=.venv-mcp", "--exclude-dir=build", "--exclude-dir=PAL2NAL", "--exclude-dir=verify"], capture_output=True, text=True)
    check("no 'chloe-newt' anywhere in the repository", g.stdout.strip() == "", g.stdout[:500])
    # a case where another code must NOT be suggested: standard-code genes with one wrong codon
    p, n = genes(rng, 10, 40)
    j, s = n[2]; n[2] = (j, s[:30] + CODE1["W" if p[2][1][10] != "W" else "A"] + s[33:])
    write(d, "pep_one.fasta", fasta(p)); write(d, "nuc_one.fasta", fasta(n))
    for eng, cmd in engines():
        r = run(cmd + ["pep_one.fasta", "nuc_one.fasta", "-o", f"one.{eng}.fa", "--report", f"one.{eng}.json"], d)
        rep = json.loads((d / f"one.{eng}.json").read_bytes())
        gc = [h for h in rep["hints"] if h["code"] == "genetic_code"]
        check(f"no false hint/{eng}: one wrong codon under the right code gives no genetic-code hint", not gc, gc)


def main(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    for fn in (report_checks, flag_checks, resume_checks, signal_checks, pairing_hint_methods):
        try:
            fn(out)
        except Exception as e:  # noqa: BLE001 - a harness error is recorded, the other parts still run
            check(f"{fn.__name__}: harness error", False, repr(e))
    (out / "results.json").write_text(json.dumps(CHECKS, indent=1), encoding="utf-8")
    bad = [c for c in CHECKS if not c["ok"]]
    print(f"\n{len(CHECKS)} checks, {len(bad)} finding candidates")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "phase4_out"))
