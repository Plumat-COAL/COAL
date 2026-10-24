#!/usr/bin/env python3
"""Phase 5: the pal2nal drop-in on typical old-pipeline command lines, against the real PAL2NAL v14.

    PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl python verify/phase5.py OUTDIR
"""
from __future__ import annotations

import json
import os
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adversarial import CODE1, ENV, PAL2NAL, clustal, fasta, gene, genes, write  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DROPIN = [sys.executable, "-m", "coal_compat"]
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append({"name": name, "ok": bool(ok), "detail": str(detail)[:600]})
    print(("OK   " if ok else "FINDING-CANDIDATE ") + name + ("" if ok else f": {str(detail)[:300]}"), flush=True)


def run(cmd, cwd, stdin=None, env=None, shell=False):
    return subprocess.run(cmd, cwd=cwd, input=stdin, capture_output=True, timeout=600, env=env or ENV, shell=shell)


def main(outdir):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(55)
    p, n = genes(rng, 5, 70)
    j, s = n[1]; n[1] = (j, s[:30] + CODE1["W" if p[1][1][10] != "W" else "A"] + s[33:])  # one mismatch
    write(out, "pep.fasta", fasta(p)); write(out, "nuc.fasta", fasta(n))
    marks = "".join("#" if 10 <= c < 60 else " " for c in range(70))
    write(out, "pep.aln", clustal(p, marks))
    aas2 = dict(CODE1, W="TGA", M="ATA")
    p2, n2 = [], []
    for k in range(4):
        pp, ss = gene(rng, 50, aas2); p2.append((f"m{k}", pp)); n2.append((f"m{k}", ss))
    write(out, "pep2.fasta", fasta(p2)); write(out, "nuc2.fasta", fasta(n2))
    lines = [
        ("default (no -output)", ["pep.fasta", "nuc.fasta"]),
        ("-output clustal", ["pep.fasta", "nuc.fasta", "-output", "clustal"]),
        ("-output paml", ["pep.fasta", "nuc.fasta", "-output", "paml"]),
        ("-output fasta", ["pep.fasta", "nuc.fasta", "-output", "fasta"]),
        ("-output codon", ["pep.fasta", "nuc.fasta", "-output", "codon"]),
        ("-codontable 2", ["pep2.fasta", "nuc2.fasta", "-codontable", "2"]),
        ("-codontable 2 -output paml -nogap", ["pep2.fasta", "nuc2.fasta", "-codontable", "2", "-output", "paml", "-nogap"]),
        ("-nogap", ["pep.fasta", "nuc.fasta", "-nogap"]),
        ("-nomismatch", ["pep.fasta", "nuc.fasta", "-nomismatch"]),
        ("-nomismatch -output fasta", ["pep.fasta", "nuc.fasta", "-nomismatch", "-output", "fasta"]),
        ("-blockonly (CLUSTAL with marks)", ["pep.aln", "nuc.fasta", "-blockonly"]),
        ("-blockonly -output paml", ["pep.aln", "nuc.fasta", "-blockonly", "-output", "paml"]),
        ("-nostderr", ["pep.fasta", "nuc.fasta", "-nostderr"]),
        ("options before files", ["-output", "fasta", "-codontable", "1", "pep.fasta", "nuc.fasta"]),
        ("dangling -codontable", ["pep.fasta", "nuc.fasta", "-codontable"]),
        ("dangling -output", ["pep.fasta", "nuc.fasta", "-output"]),
        ("-codontable 2x (means 2)", ["pep2.fasta", "nuc2.fasta", "-codontable", "2x"]),
        ("two DNA files", ["pep.fasta", "nuc.fasta", "nuc.fasta"]),
    ]
    results = {}
    for name, args in lines:
        a = run(DROPIN + args, out)
        b = run(["perl", PAL2NAL] + args, out, env=dict(os.environ))
        pal_ok = b.returncode == 0 and len(b.stdout) > 0
        same = a.stdout == b.stdout
        both_fail = (not pal_ok) and a.returncode != 0
        results[name] = {"coal_rc": a.returncode, "pal_rc": b.returncode, "pal_output": pal_ok, "same_stdout": same,
                         "coal_stderr": a.stderr.decode("utf-8", "replace").splitlines()[:5],
                         "pal_stderr": b.stderr.decode("utf-8", "replace").splitlines()[:3]}
        # PAL2NAL's -nostderr prints nothing; the drop-in maps it to -q
        check(f"drop-in {name}: stdout identical to PAL2NAL or both fail", (pal_ok and same and a.returncode == 0) or both_fail,
              f"coal rc={a.returncode} pal rc={b.returncode} pal_output={pal_ok} same={same} coal_err={results[name]['coal_stderr'][:2]} pal_err={results[name]['pal_stderr'][:2]}")
        if name == "-nostderr":
            check("drop-in -nostderr: no warnings on stderr", a.stderr == b"", a.stderr[:200])
        if name == "-nomismatch":
            check("drop-in -nomismatch: no mismatch warnings (as PAL2NAL)", b"does not correspond" not in a.stderr and b"does not correspond" not in b.stderr, a.stderr[:200])
        if name == "default (no -output)":
            check("drop-in default output is CLUSTAL", a.stdout.startswith(b"CLUSTAL"), a.stdout[:40])
            # stderr wording vs PAL2NAL: both print the per-codon warning for g1
            check("drop-in prints the per-codon warning PAL2NAL prints (gene g1)", b"pepAlnPos" in a.stderr and b"g1" in a.stderr, a.stderr[:300])
    # unsupported options: clear error, non-zero exit
    for bad in (["-html"], ["--strict"], ["-j", "2"], ["-o", "x.fa"], ["--report", "r.json"], ["-foo"]):
        a = run(DROPIN + ["pep.fasta", "nuc.fasta"] + bad, out)
        check(f"drop-in refuses {' '.join(bad)} with a clear error and exit 1", a.returncode == 1 and a.stdout == b"" and b"error" in a.stderr.lower(),
              f"rc={a.returncode} out={a.stdout[:40]!r} err={a.stderr[:200]!r}")
    a = run(DROPIN, out)
    check("drop-in without arguments: usage text, exit 1", a.returncode == 1 and (b"pal2nal" in a.stdout or b"pal2nal" in a.stderr), (a.returncode, a.stdout[:80], a.stderr[:80]))
    a = run(DROPIN + ["-h"], out)
    check("drop-in -h: usage text, exit 0", a.returncode == 0 and b"pal2nal" in a.stdout, (a.returncode, a.stdout[:80]))
    # stdin
    a = run(DROPIN + ["-", "nuc.fasta", "-output", "fasta"], out, stdin=(out / "pep.fasta").read_bytes())
    b = run(DROPIN + ["pep.fasta", "nuc.fasta", "-output", "fasta"], out)
    check("drop-in `-` reads the protein from standard input, same output", a.returncode == 0 and a.stdout == b.stdout, (a.returncode, a.stderr[:200]))
    # closed pipe: 141 silent, also for the usage text
    for what in ("pep.fasta nuc.fasta", ""):
        r = subprocess.run(f"{sys.executable} -m coal_compat {what} | head -c 10 > /dev/null; echo ${{PIPESTATUS[0]}}", shell=True, cwd=out, env=ENV, capture_output=True, timeout=600, executable="/bin/bash")
        check(f"drop-in closed pipe ({what or 'usage text'}): exit 141, silent", r.stdout.strip() == b"141" and r.stderr == b"", (r.stdout, r.stderr[:200]))
    # bad input: exit 1 (not 2)
    write(out, "bad.fasta", ">a\nMKV\n")
    a = run(DROPIN + ["bad.fasta", "nuc.fasta"], out)
    check("drop-in bad input (count mismatch): exit 1 with a message", a.returncode == 1 and b"error" in a.stderr.lower(), (a.returncode, a.stderr[:200]))
    # entry point from a fresh shell: the venv's bin/pal2nal
    exe = Path(sys.executable).parent / "pal2nal"
    r = subprocess.run(["/bin/bash", "-lc", f"cd {out} && {exe} pep.fasta nuc.fasta -output fasta | head -c 20"], capture_output=True, timeout=120)
    check("drop-in entry point works from a fresh login shell", r.returncode == 0 and r.stdout.startswith(b">"), (r.returncode, r.stdout[:40], r.stderr[-200:]))
    (out / "results.json").write_text(json.dumps({"checks": CHECKS, "lines": results}, indent=1), encoding="utf-8")
    bad = [c for c in CHECKS if not c["ok"]]
    print(f"\n{len(CHECKS)} checks, {len(bad)} finding candidates")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "phase5_out")
