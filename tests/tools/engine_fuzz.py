#!/usr/bin/env python3
"""Random engine-vs-engine differential: coal-cpp against the Python engine on a wide input space.

Every case is a random protein alignment (FASTA or CLUSTAL; gaps, '.', B/Z/J/O/U/X, '*', frameshift digits 0-9,
tfastx '\\' and '/', lower case, invalid characters, unequal row lengths, CRLF, BOM, long / Unicode / spaced /
tabbed IDs, duplicate and near-duplicate IDs) plus coding DNA (UTRs, stop tails, poly-A, truncated or empty genes,
mixed case, RNA, stray characters; DNA IDs equal, shuffled, versioned, database-prefixed, upper-cased, trimmed or
unrelated; a second DNA file, an extra or missing record) and a random option set (-output, -codontable, -nogap,
-nomismatch, -blockonly, --keep-going, --strict, -q, --max-warnings, -j, --report, -o, --json). Both engines run
with the same command line in the same directory; stdout bytes, stderr bytes, the exit code, the --report JSON and
the -o file (the volatile report fields removed) must be identical. Any difference, and any sanitizer report in
coal-cpp's stderr (AddressSanitizer, LeakSanitizer, UBSan "runtime error"), makes the exit status 1.

    python tests/tools/engine_fuzz.py --n 1500 --seed 42 --keep fuzz-keep     # COAL_CPP_BIN or build/coal-cpp
    python tests/tools/engine_fuzz.py --n 300 --sanitizer                      # coal-cpp only (an ASan/UBSan build)

--sanitizer runs coal-cpp alone and only looks for sanitizer reports (the Python engine is not run, so the
comparison is skipped: what matters is that no input makes the instrumented binary report anything). --keep DIR
copies each differing case (and each sanitizer hit) to DIR/c<N>/ (DIR/san<N>/) with the command line in args.txt;
re-run one with `coal-cpp $(cat args.txt)` in that folder. Needs only the standard library and this repo's src/.
Run by .github/workflows/nightly.yml; tests/test_engine_fuzz.py runs it with --n 20.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "src"))
from coal.tables import SUPPORTED, get_code  # noqa: E402

FORMATS = ["fasta", "clustal", "paml", "codon"]
BASES = "ACGT"
VOLATILE = ("runtime_seconds", "engine", "command", "command_line", "source", "methods_text", "coal_version")
SANITIZER = re.compile(r"AddressSanitizer|LeakSanitizer|runtime error")
# Both engines see the same environment: the Python one as the Python engine (COAL_ENGINE=python keeps `coal`
# from delegating to a coal-cpp on PATH), UTF-8 I/O, the C locale, this repo's src/ first.
ENV = dict(os.environ, COAL_ENGINE="python", PYTHONIOENCODING="utf-8", PYTHONUTF8="1", LC_ALL="C.UTF-8",
           PYTHONPATH=os.pathsep.join(p for p in (str(ROOT / "src"), os.environ.get("PYTHONPATH")) if p))


def cpp_binary(given: str | None = None) -> Path | None:
    """coal-cpp: --cpp, else COAL_CPP_BIN, else build/coal-cpp (.exe and build/Release/ tried too); None if absent."""
    p = Path(given or os.environ.get("COAL_CPP_BIN") or ROOT / "build" / "coal-cpp")
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file():
            return cand
    return None


def rid(rng, k):
    r = rng.random()
    if r < 0.5:
        return f"seq{k}"
    if r < 0.6:
        return f"sp|P{k:05d}|PROT{k}_HUMAN"
    if r < 0.7:
        return f"AT1G0{k}010.{rng.randint(1, 3)}"
    if r < 0.73:
        return "x" * rng.randint(100, 10000) + str(k)
    if r < 0.75:
        return f"gène{k}_αβ"
    if r < 0.78:
        return f"id{k} with spaces"
    if r < 0.8:
        return f"id{k}\ttab"
    if r < 0.95:
        return "".join(rng.choice("abcXYZ019_.|/:-") for _ in range(rng.randint(1, 20)))
    return f"seq{k}"


def make_case(rng):
    """(files: name -> text, options, info tags) of one random case."""
    table = rng.choice(SUPPORTED) if rng.random() < 0.5 else 1
    code = get_code(table)
    aas = sorted(a for a in code.codons if a != "*")
    nseq = rng.choice([1, 1, 2, 2, 3, 4, 6, 10])
    r = rng.random()
    ncol = 0 if r < 0.03 else rng.randint(1, 3) if r < 0.15 else rng.randint(4, 40) if r < 0.8 else rng.randint(41, 200)
    cons = [rng.choice(aas) for _ in range(ncol)]
    if ncol and rng.random() < 0.7:
        cons[0] = "M"
    special = rng.random()
    rows, dnas, info = [], [], set()
    for _ in range(nseq):
        row, dna = [], []
        for c in range(ncol):
            x = rng.random()
            if c > 0 and x < 0.08:
                row.append("-" if rng.random() < 0.9 else ".")
                continue
            aa = cons[c] if rng.random() < 0.85 else rng.choice(aas)
            if rng.random() < 0.02:
                aa = rng.choice("BZJOUX")
            if rng.random() < 0.02 and c > 0:
                aa = "*"
            if aa == "M" and c == 0 and code.starts and rng.random() < 0.3:
                codon = rng.choice(sorted(code.starts))
            elif aa in code.codons:
                codon = rng.choice(sorted(code.codons[aa]))
            else:
                codon = "".join(rng.choice(BASES) for _ in range(3))
            if rng.random() < 0.01:
                codon = "".join(rng.choice(BASES) for _ in range(3))
            if rng.random() < 0.01:
                k = rng.randrange(3)
                codon = codon[:k] + rng.choice("NRYKMSWBDHV") + codon[k + 1:]
            if c > 0 and rng.random() < 0.03:  # frameshift digit, 3/6/7/8/9/0 included
                d = rng.choice([1, 2, 4, 5, 1, 2, 3, 6, 7, 8, 9, 0])
                row.append(str(d))
                dna.append("".join(rng.choice(BASES) for _ in range(d)))
                info.add("digit")
                continue
            if c > 0 and rng.random() < 0.02:  # tfastx marks
                if rng.random() < 0.5:
                    row.append("\\")
                    dna.append(rng.choice(BASES))
                else:
                    row.append("/")
                    dna.append("".join(rng.choice(BASES) for _ in range(2)))
                info.add("tfastx")
                continue
            row.append(aa)
            dna.append(codon)
        if special < 0.03 and ncol:
            row, dna = ["-"] * len(row), []
            info.add("allgap")
        elif special < 0.05 and ncol:
            row = ["X"] * len(row)
            info.add("allX")
        elif special < 0.07 and ncol:
            row = ["*"] * len(row)
            dna = [rng.choice(sorted(code.codons.get("*", ["TAA"]))) for _ in row]
            info.add("allstop")
        elif special < 0.09 and ncol:
            row = [rng.choice("124578") for _ in row]
            dna = ["".join(rng.choice(BASES) for _ in range(int(d))) for d in row]
            info.add("alldigit")
        seq = "".join(dna)
        if rng.random() < 0.25:
            seq = "".join(rng.choice(BASES) for _ in range(rng.randint(1, 40))) + seq
        if rng.random() < 0.25:
            seq = seq + rng.choice(["TAA", "TGA", "TAG", ""]) + "".join(rng.choice(BASES) for _ in range(rng.randint(0, 30)))
        if rng.random() < 0.05:
            seq = seq + "A" * rng.randint(10, 40)
        if rng.random() < 0.05 and seq:
            seq = seq[:rng.randrange(len(seq))]  # DNA shorter than the protein
        if rng.random() < 0.03:
            seq = ""
        if rng.random() < 0.1:
            a = rng.randrange(len(seq) + 1)
            seq = seq[:a] + seq[a:].lower()
        if rng.random() < 0.05:
            seq = seq.replace("T", "U").replace("t", "u")
        if rng.random() < 0.03:
            seq = seq[:5] + rng.choice("@#!- 1") + seq[5:]
        prow = "".join(row)
        if rng.random() < 0.08:
            prow = prow.lower()
        if rng.random() < 0.02 and prow:
            k = rng.randrange(len(prow))
            prow = prow[:k] + rng.choice("@#!_ ") + prow[k + 1:]
        if rng.random() < 0.02 and prow:
            prow = prow[:-1]  # unequal row lengths
        rows.append(prow)
        dnas.append(seq)
    ids = [rid(rng, i) for i in range(nseq)]
    if rng.random() < 0.05 and nseq > 1:
        ids[1] = ids[0]
        info.add("dupid")
    if rng.random() < 0.03 and nseq > 1:
        ids[1] = ids[0] + ".2"
        info.add("neardup")
    dna_ids, order = list(ids), list(range(nseq))
    r = rng.random()
    if r < 0.15:
        rng.shuffle(order)
    elif r < 0.2:
        dna_ids = [f"d{i}" for i in range(nseq)]
    elif r < 0.25:
        dna_ids = [f"{x}.{rng.randint(1, 3)}" for x in ids]
        rng.shuffle(order)
    elif r < 0.3:
        dna_ids = [f"sp|{x}|{x.upper()}_HUMAN" for x in ids]
    elif r < 0.33:
        dna_ids = [x.upper() for x in ids]
    elif r < 0.36:
        dna_ids = [x.split(".")[0] for x in ids]
        rng.shuffle(order)
    elif r < 0.38:
        dna_ids = ids[:1] + [f"o{i}" for i in range(1, nseq)]
    # the protein file: CLUSTAL (IDs without blanks) or FASTA
    if rng.random() < 0.3 and all("\t" not in i and " " not in i for i in ids):
        idw = max(len(i) for i in ids) + rng.choice((1, 2, 6))
        blk = rng.choice((60, 50, 7, 10))
        length = max((len(x) for x in rows), default=0)
        mask = "".join(rng.choice("# ") for _ in range(length)) if rng.random() < 0.5 else None
        pep = rng.choice(["CLUSTAL W (1.83) multiple sequence alignment\n\n\n", "CLUSTAL multiple sequence alignment\n\n",
                          "MUSCLE (3.8) multiple sequence alignment\n\n"])
        b = 0
        while b < length or (b == 0 and length == 0):
            w = blk if rng.random() < 0.9 else rng.randint(1, blk)  # inconsistent block widths
            for i, s in zip(ids, rows):
                pep += i.ljust(idw) + s[b:b + w] + ("" if rng.random() < 0.9 else " 12") + "\n"
            if mask is not None:
                m = mask[b:b + w]
                if rng.random() < 0.1:
                    m = m + "##"
                if rng.random() < 0.1:
                    m = m[:-2]
                pep += " " * idw + m.rstrip() + "\n"
            elif rng.random() < 0.5:
                pep += " " * idw + "".join(rng.choice(" *:.") for _ in range(w)).rstrip() + "\n"
            pep += "\n" * rng.choice((0, 1, 1, 2))
            b += w
            if length == 0:
                break
        pf = "pep.aln"
    else:
        w = rng.choice((0, 0, 60, 7))

        def wrap(s):
            return (s + "\n") if not w or not s else "\n".join(s[i:i + w] for i in range(0, len(s), w)) + "\n"
        pep = "".join(f">{i}{rng.choice(('', ' desc', '  d1 d2'))}\n{wrap(s)}" for i, s in zip(ids, rows))
        if rng.random() < 0.02 and rows and rows[0]:
            pep = pep.replace(rows[0][:1], ">" + rows[0][:1], 1)  # '>' inside a sequence
        if rng.random() < 0.02:
            pep = "\n" + pep
        if rng.random() < 0.01:
            pep = "﻿" + pep
        if rng.random() < 0.02:
            pep = pep.replace("\n", "\r\n")
        pf = "pep.fasta"
    w = rng.choice((0, 60, 13))

    def wrapn(s):
        return (s + "\n") if not w or not s else "\n".join(s[i:i + w] for i in range(0, len(s), w)) + "\n"
    nuc = "".join(f">{dna_ids[i]}{rng.choice(('', ' desc'))}\n{wrapn(dnas[i])}" for i in order)
    if rng.random() < 0.02:
        nuc = nuc + ">extra\nATGAAA\n"
    if rng.random() < 0.02 and nseq > 1:
        nuc = nuc.split("\n", 2)[-1]  # headerless first record
    if rng.random() < 0.02:
        nuc = nuc.replace("\n", "\r\n")
    files = {pf: pep, "nuc.fasta": nuc}
    if rng.random() < 0.1 and nseq > 1:  # the DNA split over two files
        recs = nuc.split(">")[1:]
        k = rng.randint(1, len(recs) - 1) if len(recs) > 1 else 1
        files["nuc.fasta"] = "".join(">" + x for x in recs[:k])
        files["nuc2.fasta"] = "".join(">" + x for x in recs[k:])
    opts = ["-output", rng.choice(FORMATS)]
    if table != 1 or rng.random() < 0.2:
        opts += ["-codontable", str(table if rng.random() < 0.85 else rng.choice(SUPPORTED))]
    x = rng.random()
    if x < 0.2:
        opts.append("-nogap")
    elif x < 0.35:
        opts.append("-nomismatch")
    elif x < 0.4:
        opts += ["-nogap", "-nomismatch"]
    if rng.random() < 0.3:
        opts.append("-blockonly")
    if rng.random() < 0.15:
        opts.append("--keep-going")
    if rng.random() < 0.15:
        opts.append("--strict")
    if rng.random() < 0.1:
        opts.append("-q")
    if rng.random() < 0.1:
        opts += ["--max-warnings", str(rng.randint(0, 3))]
    if rng.random() < 0.2:
        opts += ["-j", rng.choice(("1", "2"))]
    rep = rng.random() < 0.35
    if rep:
        opts += ["--report", "rep.json"]
    if rep and rng.random() < 0.4:
        opts += ["-o", "out.txt", "--json"]
    return files, opts, info


def _strip(obj):
    for k in VOLATILE:
        obj.pop(k, None)
    return obj


def run(cmd, d, files, opts, timeout):
    """(exit code, stdout, stderr, report JSON or None, -o file bytes or None) of one engine on a case."""
    pep = [f for f in files if f.startswith("pep")][0]
    args = [*cmd, pep, *sorted(f for f in files if f.startswith("nuc")), *opts]
    p = subprocess.run(args, cwd=d, capture_output=True, timeout=timeout, env=ENV)
    rep = out = None
    rp = Path(d) / "rep.json"
    if rp.exists():
        try:
            rep = _strip(json.loads(rp.read_bytes()))
        except ValueError as e:
            rep = f"BADJSON {e}"
        rp.unlink()
    op = Path(d) / "out.txt"
    if op.exists():
        out = op.read_bytes()
        op.unlink()
    so = p.stdout
    if "--json" in opts and so.strip():
        try:
            so = json.dumps(_strip(json.loads(so)), sort_keys=True).encode()
        except ValueError:
            so = b"BADJSON"
    return p.returncode, so, p.stderr, rep, out


def case_args(files, opts):
    pep = [f for f in files if f.startswith("pep")][0]
    return [pep, *sorted(f for f in files if f.startswith("nuc")), *opts]


def keep_case(keep, d, name, files, opts):
    dst = Path(keep) / name
    shutil.copytree(d, dst, dirs_exist_ok=True)
    (dst / "args.txt").write_text(" ".join(case_args(files, opts)) + "\n", encoding="utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=500, help="number of random cases (default 500)")
    ap.add_argument("--seed", type=int, default=1, help="random seed (default 1); the summary line repeats it")
    ap.add_argument("--keep", metavar="DIR", help="copy each differing case to DIR/c<N>/ (sanitizer hits: DIR/san<N>/)")
    ap.add_argument("--sanitizer", action="store_true",
                    help="run coal-cpp only and fail on any AddressSanitizer / LeakSanitizer / UBSan report")
    ap.add_argument("--cpp", metavar="BIN", help="coal-cpp binary (default: COAL_CPP_BIN, else build/coal-cpp)")
    ap.add_argument("--python", metavar="EXE", default=sys.executable,
                    help="interpreter that runs the Python engine as `-m coal` (default: this one)")
    ap.add_argument("--timeout", type=float, default=300, help="seconds per engine run (default 300)")
    a = ap.parse_args(argv)

    cpp = cpp_binary(a.cpp)
    if cpp is None:
        print("engine_fuzz: coal-cpp not found (build it, or set COAL_CPP_BIN / --cpp)", file=sys.stderr)
        return 2
    py = [a.python, "-m", "coal"]
    rng = random.Random(a.seed)
    base = tempfile.mkdtemp(prefix="coal_engine_fuzz_")
    ndiff = nsan = 0
    rc_stats, errs = {}, {}
    try:
        for i in range(a.n):
            files, opts, info = make_case(rng)
            d = os.path.join(base, f"c{i}")
            os.mkdir(d)
            for name, text in files.items():
                with open(os.path.join(d, name), "w", newline="", encoding="utf-8") as fh:
                    fh.write(text)
            cc = run([str(cpp)], d, files, opts, a.timeout)
            err = cc[2].decode("latin-1")
            if SANITIZER.search(err):
                nsan += 1
                print(f"SANITIZER #{i} opts={opts} info={sorted(info)}\n{err[-1500:]}", flush=True)
                if a.keep:
                    keep_case(a.keep, d, f"san{i}", files, opts)
            if not a.sanitizer:
                pr = run(py, d, files, opts, a.timeout)
                rc_stats[pr[0]] = rc_stats.get(pr[0], 0) + 1
                if pr[0] != 0:
                    m = re.search(rb"coal: error: ([^\n]{0,60})", pr[2])
                    k = (m.group(1).decode("latin-1") if m else "?")[:50]
                    errs[k] = errs.get(k, 0) + 1
                if pr != cc:
                    ndiff += 1
                    what = [n for n, x, y in zip(("rc", "stdout", "stderr", "report", "outfile"), pr, cc) if x != y]
                    print(f"DIFF #{i} opts={opts} info={sorted(info)} what={what}\n"
                          f"  py: rc={pr[0]} out={pr[1][:300]!r} err={pr[2][-400:]!r}\n"
                          f"  cc: rc={cc[0]} out={cc[1][:300]!r} err={cc[2][-400:]!r}", flush=True)
                    if "report" in what:
                        print(f"  report py={json.dumps(pr[3])[:600]}\n  report cc={json.dumps(cc[3])[:600]}", flush=True)
                    if a.keep:
                        keep_case(a.keep, d, f"c{i}", files, opts)
            shutil.rmtree(d, ignore_errors=True)
    finally:
        shutil.rmtree(base, ignore_errors=True)
    mode = "sanitizer run (coal-cpp only)" if a.sanitizer else "coal-cpp vs python"
    print(f"engine fuzz, {mode}, seed {a.seed}: {a.n} cases, {ndiff} differences, {nsan} sanitizer reports"
          + ("" if a.sanitizer else f", python exit codes {rc_stats}"))
    for k, v in sorted(errs.items(), key=lambda x: -x[1]):
        print(f"  {v:4d} {k}")
    return 1 if ndiff or nsan else 0


if __name__ == "__main__":
    sys.exit(main())
