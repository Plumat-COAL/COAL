"""Race: PAL2NAL (once per gene, as typical pipelines run it) vs one batched COAL run.

    python demo/race.py                         # 300 synthetic genes x 400 codons
    python demo/race.py --genes 1000 --codons 500
    python demo/race.py --protein pep.aln --dna nuc.fasta   # your own data

PAL2NAL is NOT included (GPL-2.0). This script runs *your* local copy, by default
PAL2NAL/pal2nal.v14/pal2nal.pl as fetched by `python scripts/fetch_pal2nal.py`.
It never imports or ships PAL2NAL. The C++ COAL binary is raced too if it is built
(cpp/build/...) or given with --cpp.
"""
from __future__ import annotations

import argparse
import glob
import importlib.util
import os
import random
import shutil
import subprocess
import sys
import tempfile
import threading
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_PAL2NAL = os.path.join(ROOT, "PAL2NAL", "pal2nal.v14", "pal2nal.pl")

HOW_TO_FETCH = f"""\
PAL2NAL was not found at: {{path}}

This demo races COAL against your own local copy of PAL2NAL (it is GPL-2.0, so it is not
included in COAL). To fetch it (needs Perl 5 to run):

    python scripts/fetch_pal2nal.py

If the download is blocked, download pal2nal.v14.tar.gz from https://www.bork.embl.de/pal2nal/
in a browser, put it in {ROOT} and run:

    python scripts/fetch_pal2nal.py --archive pal2nal.v14.tar.gz

Or point this script at an existing copy:  python demo/race.py --pal2nal /path/to/pal2nal.pl
"""


# ---------------------------------------------------------------- data

def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0] if line[1:].split() else "", []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_fasta(path, recs):
    with open(path, "w", newline="\n") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def make_synthetic(workdir, genes, codons, seed=7):
    """Synthetic aligned proteins + their DNA, using examples/make_examples.py helpers."""
    spec = importlib.util.spec_from_file_location("make_examples", os.path.join(ROOT, "examples", "make_examples.py"))
    mx = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mx)
    rng = random.Random(seed)
    cmap = mx.codon_map(mx.CODE1)
    rows = mx.make_family(rng, genes, codons)
    names = [f"gene{i + 1:05d}" for i in range(genes)]
    prot = os.path.join(workdir, "race.pep.fasta")
    dna = os.path.join(workdir, "race.nuc.fasta")
    write_fasta(prot, list(zip(names, rows)))
    write_fasta(dna, [(n, mx.back_translate(rng, r.replace("-", ""), cmap)) for n, r in zip(names, rows)])
    return prot, dna


# ---------------------------------------------------------------- racers

class Racer:
    def __init__(self, label, total):
        self.label, self.total = label, total
        self.done, self.t0, self.t1, self.error = 0, None, None, None
        self.records = None

    @property
    def elapsed(self):
        if self.t0 is None:
            return 0.0
        return (self.t1 or time.perf_counter()) - self.t0


def run_pal2nal_per_gene(r, perl, pal2nal, prot_recs, dna_recs, workdir, deadline):
    dna = dict(dna_recs)
    out = []
    r.t0 = time.perf_counter()
    try:
        for name, row in prot_recs:
            p = os.path.join(workdir, "one.pep")
            d = os.path.join(workdir, "one.nuc")
            write_fasta(p, [(name, row)])
            write_fasta(d, [(name, dna[name])])
            left = deadline - time.perf_counter()
            if left <= 0:
                raise TimeoutError("time limit reached")
            res = subprocess.run([perl, pal2nal, p, d, "-output", "fasta"], capture_output=True, text=True,
                                 timeout=left)
            if res.returncode != 0 or not res.stdout.startswith(">"):
                raise RuntimeError(f"PAL2NAL failed on {name}: {res.stderr.strip()[:200]}")
            fa = os.path.join(workdir, "one.out")
            with open(fa, "w") as fh:
                fh.write(res.stdout)
            out.extend(read_fasta(fa))
            r.done += 1
        r.records = out
    except Exception as e:  # shown in the summary
        r.error = str(e)
    r.t1 = time.perf_counter()


def run_batched(r, cmd, out_path, deadline):
    r.t0 = time.perf_counter()
    try:
        res = subprocess.run(cmd + ["-o", out_path], capture_output=True, text=True,
                             timeout=max(1.0, deadline - time.perf_counter()))
        if res.returncode not in (0, 2):
            raise RuntimeError(res.stderr.strip()[:300])
        r.records = read_fasta(out_path)
        r.done = r.total
    except Exception as e:
        r.error = str(e)
    r.t1 = time.perf_counter()


# ---------------------------------------------------------------- display

def bar(r, width=28):
    frac = r.done / r.total if r.total else 1
    fill = int(frac * width)
    state = "error" if r.error else ("done" if r.t1 else "")
    return f"  {r.label:<26} [{'#' * fill}{'.' * (width - fill)}] {r.done:>6}/{r.total} genes {r.elapsed:8.2f} s  {state}"


def find_cpp():
    pats = ["cpp/build/coal", "cpp/build/coal.exe", "cpp/build/Release/coal.exe", "cpp/build/*/coal*"]
    for pat in pats:
        for p in glob.glob(os.path.join(ROOT, pat)):
            if os.path.isfile(p) and os.access(p, os.X_OK):
                return p
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pal2nal", default=DEFAULT_PAL2NAL, help="path to your local pal2nal.pl")
    ap.add_argument("--perl", default=shutil.which("perl") or "perl", help="Perl interpreter")
    ap.add_argument("--cpp", help="path to the C++ COAL binary (default: look in cpp/build/)")
    ap.add_argument("--protein", help="your protein alignment (FASTA); default: synthetic data")
    ap.add_argument("--dna", help="your coding DNA (FASTA)")
    ap.add_argument("--genes", type=int, default=300)
    ap.add_argument("--codons", type=int, default=400)
    ap.add_argument("--timeout", type=float, default=900, help="give up after this many seconds")
    args = ap.parse_args(argv)

    if not os.path.isfile(args.pal2nal):
        print(HOW_TO_FETCH.format(path=args.pal2nal))
        return 0
    if not shutil.which(args.perl) and not os.path.isfile(args.perl):
        print("Perl was not found. PAL2NAL needs Perl 5 (Windows: Strawberry Perl). Use --perl /path/to/perl.")
        return 0
    if bool(args.protein) != bool(args.dna):
        ap.error("give both --protein and --dna, or neither")
    if os.name == "nt":
        os.system("")  # enable ANSI escape codes in the Windows console

    work = tempfile.mkdtemp(prefix="coal_race_")
    try:
        if args.protein:
            prot, dna = args.protein, args.dna
        else:
            print(f"Making {args.genes} synthetic genes x {args.codons} codons ...")
            prot, dna = make_synthetic(work, args.genes, args.codons)
        prot_recs, dna_recs = read_fasta(prot), read_fasta(dna)
        n = len(prot_recs)
        missing = [p for p, _ in prot_recs if p not in dict(dna_recs)]
        if missing:
            print(f"Protein and DNA IDs must match for the per-gene run (e.g. {missing[0]} has no DNA).")
            return 1
        deadline = time.perf_counter() + args.timeout
        racers = []
        threads = []
        r_p = Racer("PAL2NAL (once per gene)", n)
        racers.append(r_p)
        threads.append(threading.Thread(target=run_pal2nal_per_gene, daemon=True,
                                        args=(r_p, args.perl, args.pal2nal, prot_recs, dna_recs,
                                              os.path.join(work, "p"), deadline)))
        os.makedirs(os.path.join(work, "p"))
        r_c = Racer("COAL (one batched run)", n)
        racers.append(r_c)
        threads.append(threading.Thread(target=run_batched, daemon=True, args=(
            r_c, [sys.executable, "-m", "coal", prot, dna, "-j", "1"], os.path.join(work, "coal.fa"), deadline)))
        cpp = args.cpp or find_cpp()
        if cpp:
            r_x = Racer("COAL C++ (one batched run)", n)
            racers.append(r_x)
            threads.append(threading.Thread(target=run_batched, daemon=True,
                                            args=(r_x, [cpp, prot, dna], os.path.join(work, "cpp.fa"), deadline)))

        tty = sys.stdout.isatty()
        print(f"\nRace: {n} genes\n")
        for t in threads:
            t.start()
        stop = threading.Event()
        first = True
        while any(t.is_alive() for t in threads) and time.perf_counter() < deadline + 5:
            if tty:
                if not first:
                    sys.stdout.write(f"\x1b[{len(racers)}A")
                sys.stdout.write("\n".join("\x1b[2K" + bar(r) for r in racers) + "\n")
                sys.stdout.flush()
                first = False
            stop.wait(0.1)
        for t in threads:
            t.join(timeout=5)
        if tty and not first:
            sys.stdout.write(f"\x1b[{len(racers)}A")
        print("\n".join(("\x1b[2K" if tty else "") + bar(r) for r in racers))

        print("\nResults")
        base = racers[0]
        for r in racers:
            if r.error:
                print(f"  {r.label}: ERROR {r.error}")
            else:
                speed = f"  ({base.elapsed / r.elapsed:.1f}x faster than PAL2NAL)" if r is not base and r.elapsed else ""
                print(f"  {r.label}: {r.elapsed:.2f} s{speed}")
        if base.records is not None:
            ref = dict(base.records)
            for r in racers[1:]:
                if r.records is not None:
                    same = dict(r.records) == ref and len(r.records) == len(base.records)
                    print(f"  {r.label} output {'MATCHES' if same else 'DIFFERS FROM'} PAL2NAL's.")
        return 0 if not any(r.error for r in racers) else 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
