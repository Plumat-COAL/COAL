"""Race: PAL2NAL once per gene family (as pipelines call it) vs one batched COAL run on all families.

    python demo/race.py                              # 300 synthetic gene families x 4 sequences x 400 codons
    python demo/race.py --families 1000 --seqs 2     # e.g. 1000 ortholog pairs
    python demo/race.py --family pep1.aln nuc1.fa --family pep2.aln nuc2.fa   # your own families

Pipelines (OrthoFinder-style ortholog groups, pairwise dN/dS) call PAL2NAL once per gene family:
one protein alignment + its DNA per call. COAL can do the same, but it can also take all families in
one run (rows padded with trailing gaps to one width), which is what is raced here. Outputs are checked
family by family against PAL2NAL's (COAL's padding columns removed).

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
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def make_synthetic(workdir, families, seqs, codons, seed=7):
    """Synthetic gene families (aligned proteins + their DNA), using examples/make_examples.py helpers.
    Returns [(protein_path, dna_path)], one pair of files per family."""
    spec = importlib.util.spec_from_file_location("make_examples", os.path.join(ROOT, "examples", "make_examples.py"))
    mx = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mx)
    rng = random.Random(seed)
    cmap = mx.codon_map(mx.CODE1)
    out = []
    for f in range(families):
        rows = mx.make_family(rng, seqs, codons)
        names = [f"fam{f + 1:05d}_sp{k + 1}" for k in range(seqs)]
        prot = os.path.join(workdir, f"fam{f + 1:05d}.pep.fasta")
        dna = os.path.join(workdir, f"fam{f + 1:05d}.nuc.fasta")
        write_fasta(prot, list(zip(names, rows)))
        write_fasta(dna, [(n, mx.back_translate(rng, r.replace("-", ""), cmap)) for n, r in zip(names, rows)])
        out.append((prot, dna))
    return out


def batch_families(fams, workdir):
    """All families in one protein alignment + one DNA file (rows padded with trailing gaps).
    Returns (protein_path, dna_path, {id: protein row width})."""
    prot_all, dna_all, width = [], [], {}
    for p, d in fams:
        prot_all += read_fasta(p)
        dna_all += read_fasta(d)
    ids = [n for n, _ in prot_all]
    if len(set(ids)) != len(ids):
        raise ValueError("sequence IDs must be unique across families to batch them")
    w = max((len(r) for _, r in prot_all), default=0)
    for n, r in prot_all:
        width[n] = len(r)
    prot = os.path.join(workdir, "all.pep.fasta")
    dna = os.path.join(workdir, "all.nuc.fasta")
    write_fasta(prot, [(n, r + "-" * (w - len(r))) for n, r in prot_all])
    write_fasta(dna, dna_all)
    return prot, dna, width


def unpad(records, ref_len):
    """Cut COAL's rows from the batch to the length of PAL2NAL's row for the same ID (the rest is the
    padding added to batch the families). Returns None if anything but gaps was cut off, or an ID is unknown."""
    out = []
    for n, s in records:
        if n not in ref_len or s[ref_len[n]:].strip("-"):
            return None
        out.append((n, s[:ref_len[n]]))
    return out


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


# CHANGED(APP-8): was PAL2NAL once per sequence of one family (measured process start-up); now once per gene
# family vs one COAL run on all families. See overnight/app/CHANGES.md
def run_pal2nal_per_family(r, perl, pal2nal, fams, workdir, deadline):
    out = []
    r.t0 = time.perf_counter()
    try:
        for p, d in fams:
            left = deadline - time.perf_counter()
            if left <= 0:
                raise TimeoutError("time limit reached")
            res = subprocess.run([perl, pal2nal, p, d, "-output", "fasta"], capture_output=True, text=True,
                                 encoding="utf-8", errors="replace", timeout=left)
            if res.returncode != 0 or not res.stdout.startswith(">"):
                raise RuntimeError(f"PAL2NAL failed on {os.path.basename(p)}: {res.stderr.strip()[:200]}")
            fa = os.path.join(workdir, "one.out")
            with open(fa, "w", encoding="utf-8") as fh:
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
        res = subprocess.run(cmd + ["-o", out_path], capture_output=True, text=True, encoding="utf-8",
                             errors="replace", timeout=max(1.0, deadline - time.perf_counter()))
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
    return f"  {r.label:<30} [{'#' * fill}{'.' * (width - fill)}] {r.done:>6}/{r.total} families {r.elapsed:8.2f} s  {state}"


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
    ap.add_argument("--family", nargs=2, action="append", metavar=("PEP_ALN", "NUC_FASTA"),
                    help="one of your gene families (protein alignment FASTA + its DNA); repeat it. "
                         "Default: synthetic families")
    ap.add_argument("--families", type=int, default=300, help="synthetic gene families (default 300)")
    ap.add_argument("--seqs", type=int, default=4, help="sequences per synthetic family (default 4)")
    ap.add_argument("--codons", type=int, default=400, help="codons per synthetic gene (default 400)")
    ap.add_argument("--timeout", type=float, default=900, help="give up after this many seconds")
    args = ap.parse_args(argv)

    if not os.path.isfile(args.pal2nal):
        print(HOW_TO_FETCH.format(path=args.pal2nal))
        return 0
    if not shutil.which(args.perl) and not os.path.isfile(args.perl):
        print("Perl was not found. PAL2NAL needs Perl 5 (Windows: Strawberry Perl). Use --perl /path/to/perl.")
        return 0
    if os.name == "nt":
        os.system("")  # enable ANSI escape codes in the Windows console

    work = tempfile.mkdtemp(prefix="coal_race_")
    try:
        if args.family:
            fams = [tuple(x) for x in args.family]
        else:
            print(f"Making {args.families} synthetic gene families x {args.seqs} sequences x {args.codons} codons ...")
            fams = make_synthetic(work, args.families, args.seqs, args.codons)
        try:
            prot, dna, width = batch_families(fams, work)
        except ValueError as e:
            print(f"Cannot batch these families: {e}")
            return 1
        n = len(fams)
        deadline = time.perf_counter() + args.timeout
        racers = []
        threads = []
        r_p = Racer("PAL2NAL (once per family)", n)
        racers.append(r_p)
        os.makedirs(os.path.join(work, "p"))
        threads.append(threading.Thread(target=run_pal2nal_per_family, daemon=True,
                                        args=(r_p, args.perl, args.pal2nal, fams, os.path.join(work, "p"),
                                              deadline)))
        r_c = Racer("COAL (one run, all families)", n)
        racers.append(r_c)
        threads.append(threading.Thread(target=run_batched, daemon=True, args=(
            r_c, [sys.executable, "-m", "coal", prot, dna, "-j", "1"], os.path.join(work, "coal.fa"), deadline)))
        cpp = args.cpp or find_cpp()
        if cpp:
            r_x = Racer("COAL C++ (one run)", n)
            racers.append(r_x)
            threads.append(threading.Thread(target=run_batched, daemon=True,
                                            args=(r_x, [cpp, prot, dna], os.path.join(work, "cpp.fa"), deadline)))

        tty = sys.stdout.isatty()
        nseq = len(width)
        print(f"\nRace: {n} gene families ({nseq} sequences). PAL2NAL is called once per gene family, "
              f"as pipelines call it; COAL runs once on all of them.\n")
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
                    mine = unpad(r.records, {k: len(v) for k, v in ref.items()})
                    same = mine is not None and dict(mine) == ref and len(mine) == len(base.records)
                    print(f"  {r.label} output {'MATCHES' if same else 'DIFFERS FROM'} PAL2NAL's.")
        return 0 if not any(r.error for r in racers) else 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
