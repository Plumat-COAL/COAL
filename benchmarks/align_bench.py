"""Benchmark of the protein-alignment layer (docs/alignment.md "Benchmark"): global NW (coal-cpp-align) vs the
DIAMOND gapped-alignment route, on orthologr's 1,000-gene A. thaliana / A. lyrata example.

    python benchmarks/align_bench.py \
        --cds-a orthologr/inst/seqs/ortho_thal_cds_1000.fasta --cds-b orthologr/inst/seqs/ortho_lyra_cds_1000.fasta \
        --diamond diamond --out bench_align/ [--variant scalar=/path/to/coal-cpp-align ...]

Steps: translate the CDS (standard code, trailing stop removed) -> DIAMOND blastp both ways with
--outfmt 6 qseqid sseqid qstart qend sstart send qseq_gapped sseq_gapped bitscore -> reciprocal best hits ->
protein alignment of every pair by each method -> COAL codon alignment (coal-cpp) -> NG86 dN, dS
(Biopython's Bio.codonalign; only to compare the methods with each other) -> results.md.

Needs Biopython (NG86; also timed: PairwiseAligner with the same scoring); parasail (pip install parasail) is timed
too when importable. Every timing is the median of --repeats runs, wall clock, on this machine only.
"""
from __future__ import annotations

import argparse
import math
import os
import platform
import shutil
import statistics
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from coal import align as A  # noqa: E402
from coal.tables import get_code  # noqa: E402


def read_fasta(path):
    recs, cur = [], None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                cur = [line[1:].split()[0], []]
                recs.append(cur)
            elif cur is not None:
                cur[1].append(line)
    return [(i, "".join(s)) for i, s in recs]


def write_fasta(path, recs):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for i, s in recs:
            fh.write(f">{i}\n{s}\n")


def translate(cds):
    code = get_code(1).code
    s = cds.upper().replace("U", "T")
    p = "".join(code.get(s[k:k + 3], "X") for k in range(0, len(s) - len(s) % 3, 3))
    return p[:-1] if p.endswith("*") else p


def _cpu():
    import resource
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time() + c.ru_utime + c.ru_stime


def timed(fn, repeats):
    """(median wall seconds, median CPU seconds of this process and its children, result of the last run).
    CPU time is the steadier number on a shared machine."""
    ts, cs, out = [], [], None
    for _ in range(repeats):
        t, c = time.perf_counter(), _cpu()
        out = fn()
        ts.append(time.perf_counter() - t)
        cs.append(_cpu() - c)
    TIMES.append((statistics.median(ts), statistics.median(cs)))
    return statistics.median(ts), out


TIMES: list = []


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode not in (0, 2):
        raise RuntimeError(f"{cmd[0]} failed ({r.returncode}): {r.stderr[-2000:]}")
    return r


MASKING: list = []  # extra blastp options (--masking)


def blastp(diamond, q, db, out, threads, extra_cols=True):
    cols = list(A.DIAMOND_COLUMNS) + ["bitscore"] if extra_cols else ["qseqid", "sseqid", "bitscore"]
    run([diamond, "blastp", "-q", q, "-d", db, "-o", out, "--outfmt", "6", *cols, "--max-target-seqs", "5",
         "--evalue", "1e-5", "--threads", str(threads), "--quiet", *MASKING])


def best(path, bitcol):
    b = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] not in b or float(f[bitcol]) > b[f[0]][1]:
                b[f[0]] = (f[1], float(f[bitcol]))
    return {q: s for q, (s, _) in b.items()}


def read_pairs_fasta(path):
    rows, name = [], None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
            else:
                rows.append((name, line))
    return [(rows[k][1], rows[k + 1][1]) for k in range(0, len(rows), 2)]


def codon_align(coal_cpp, alns, cds_a, cds_b, work, tag):
    """Codon rows for each pair: one coal-cpp run on the batch (pairs padded with trailing gaps)."""
    L = max(len(a) for a, _ in alns)
    pep, nuc = [], []
    for k, ((ra, rb), (ia, ib)) in enumerate(alns_ids(alns)):
        pep += [(f"p{k}a", ra + "-" * (L - len(ra))), (f"p{k}b", rb + "-" * (L - len(rb)))]
        nuc += [(f"p{k}a", cds_a[ia]), (f"p{k}b", cds_b[ib])]
    pp, nn, oo = (os.path.join(work, f"{tag}.{x}") for x in ("pep.fasta", "nuc.fasta", "codon.fasta"))
    write_fasta(pp, pep)
    write_fasta(nn, nuc)
    run([coal_cpp, pp, nn, "-o", oo, "--keep-going", "-q"])
    got = dict(read_fasta(oo))
    out = []
    for k, ((ra, _), _) in enumerate(alns_ids(alns)):
        x, y = got.get(f"p{k}a"), got.get(f"p{k}b")
        out.append(None if x is None or y is None else (x[:3 * len(ra)], y[:3 * len(ra)]))
    return out


def alns_ids(alns):
    return [((a, b), ids) for (a, b), ids in zip(alns, ALN_IDS)]


ALN_IDS: list = []


def ng86(rows):
    import warnings
    from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
    a, b = rows
    # codons with a gap in either sequence, or a stop, are dropped (pairwise deletion)
    ka, kb, n = [], [], 0
    stops = {"TAA", "TAG", "TGA"}
    for k in range(0, len(a) - 2, 3):
        x, y = a[k:k + 3].upper(), b[k:k + 3].upper()
        if "-" in x or "-" in y or x in stops or y in stops or not set(x + y) <= set("ACGT"):
            continue
        ka.append(x)
        kb.append(y)
        n += 1
    if n == 0:
        return math.nan, math.nan, 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            dn, ds = cal_dn_ds(CodonSeq("".join(ka)), CodonSeq("".join(kb)), method="NG86")
        except Exception:  # noqa: BLE001
            return math.nan, math.nan, n
    return (dn if dn >= 0 else math.nan), (ds if ds >= 0 else math.nan), n


def pearson(x, y):
    pts = [(a, b) for a, b in zip(x, y) if a == a and b == b and math.isfinite(a) and math.isfinite(b)]
    if len(pts) < 3:
        return math.nan, len(pts)
    mx = sum(a for a, _ in pts) / len(pts)
    my = sum(b for _, b in pts) / len(pts)
    sxy = sum((a - mx) * (b - my) for a, b in pts)
    sxx = sum((a - mx) ** 2 for a, _ in pts)
    syy = sum((b - my) ** 2 for _, b in pts)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else math.nan, len(pts)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cds-a", required=True)
    ap.add_argument("--cds-b", required=True)
    ap.add_argument("--diamond", default=shutil.which("diamond"))
    ap.add_argument("--coal-cpp", default=os.path.join(ROOT, "build", "coal-cpp"))
    ap.add_argument("--coal-cpp-align", default=os.path.join(ROOT, "build", "coal-cpp-align"))
    ap.add_argument("--variant", action="append", default=[], metavar="NAME=PATH",
                    help="another coal-cpp-align build to time (e.g. scalar=build-scalar/coal-cpp-align)")
    ap.add_argument("--jobs", type=int, nargs="+", default=[1, 0])
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 1, help="DIAMOND threads")
    ap.add_argument("--clustalo", help="Clustal Omega binary: also align every pair as a 2-sequence cluster")
    ap.add_argument("--mafft", help="MAFFT binary: the same with MAFFT")
    ap.add_argument("--masking", help="DIAMOND --masking (default: DIAMOND's own, tantan); 0 = off")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.masking is not None:
        MASKING[:] = ["--masking", a.masking]
    os.makedirs(a.out, exist_ok=True)
    work = a.out
    lines = []

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    cds_a, cds_b = dict(read_fasta(a.cds_a)), dict(read_fasta(a.cds_b))
    prot_a = {i: translate(s) for i, s in cds_a.items()}
    prot_b = {i: translate(s) for i, s in cds_b.items()}
    fa, fb = os.path.join(work, "a.faa"), os.path.join(work, "b.faa")
    write_fasta(fa, prot_a.items())
    write_fasta(fb, prot_b.items())
    da, db = os.path.join(work, "a_db"), os.path.join(work, "b_db")
    run([a.diamond, "makedb", "--in", fa, "-d", da, "--quiet"])
    run([a.diamond, "makedb", "--in", fb, "-d", db, "--quiet"])
    ab, ba = os.path.join(work, "ab.tsv"), os.path.join(work, "ba.tsv")
    t_dmd_gapped, _ = timed(lambda: (blastp(a.diamond, fa, db, ab, a.threads), blastp(a.diamond, fb, da, ba,
                                                                                    a.threads)), a.repeats)
    t_dmd_plain, _ = timed(lambda: (blastp(a.diamond, fa, db, ab + ".plain", a.threads, False),
                                    blastp(a.diamond, fb, da, ba + ".plain", a.threads, False)), a.repeats)
    bab, bba = best(ab, 8), best(ba, 8)
    rbh = sorted((q, s) for q, s in bab.items() if bba.get(s) == q)
    ALN_IDS[:] = rbh
    pairs = [(prot_a[q], prot_b[s]) for q, s in rbh]
    cells = sum(len(x) * len(y) for x, y in pairs)
    # the RBH hit lines (first line per pair in ab.tsv) for the DIAMOND route
    rbhset = set(rbh)
    hits_tsv = os.path.join(work, "rbh_hits.tsv")
    with open(ab) as fh, open(hits_tsv, "w") as out:
        seen = set()
        for line in fh:
            f = line.split("\t")
            if (f[0], f[1]) in rbhset and (f[0], f[1]) not in seen:
                seen.add((f[0], f[1]))
                out.write(line)
    hits = {(h.qseqid, h.sseqid): h for h in A.read_diamond_hits(hits_tsv)}
    cover = sorted(min(h.qend - h.qstart + 1, h.send - h.sstart + 1) / min(len(prot_a[q]), len(prot_b[s]))
                   for (q, s), h in hits.items())

    say(f"# Protein-alignment layer benchmark\n")
    say(f"Machine: {platform.platform()}, {os.cpu_count()} hardware threads, Python {platform.python_version()}. "
        f"DIAMOND: {run([a.diamond, 'version']).stdout.strip()}.")
    say(f"Data: {len(cds_a)} + {len(cds_b)} CDS -> {len(rbh)} DIAMOND reciprocal best hits; "
        f"{cells / 1e6:.1f} M DP cells; median protein length "
        f"{statistics.median([len(x) for x, _ in pairs] + [len(y) for _, y in pairs]):.0f}.")
    say(f"DIAMOND hit coverage of the shorter protein: median {statistics.median(cover):.3f}, "
        f"{sum(c < 0.9 for c in cover)} pairs < 0.9, {sum(c < 0.5 for c in cover)} pairs < 0.5.\n")

    pf = os.path.join(work, "pairs.faa")
    write_fasta(pf, [(f"{k}{x}", s) for k, (p, q) in enumerate(pairs) for x, s in (("a", p), ("b", q))])
    say("## Time (protein alignment of all pairs; median of %d runs)\n" % a.repeats)
    say("| Method | Threads | Wall (s) | CPU (s) | Pairs/s (wall) | GCUPS (CPU) |")
    say("|---|---:|---:|---:|---:|---:|")
    methods = {}

    def nw(binary, params_args, jobs, tag):
        o = os.path.join(work, f"{tag}.aln.fasta")
        run([binary, pf, "-o", o, "-j", str(jobs), *params_args])
        return o
    def row(name, th, t):
        cpu = TIMES[-1][1]
        say(f"| {name} | {th} | {t:.3f} | {cpu:.3f} | {len(pairs) / t:,.0f} | {cells / cpu / 1e9:.2f} |")
    for jobs in a.jobs:
        t, o = timed(lambda: nw(a.coal_cpp_align, [], jobs, f"nw_needle_j{jobs}"), a.repeats)
        row("coal-cpp-align NW, needle defaults", jobs or os.cpu_count(), t)
        methods["NW needle (reference)"] = read_pairs_fasta(o)
    t, o = timed(lambda: nw(a.coal_cpp_align, ["--preset", "pwalign"], 1, "nw_pwalign"), a.repeats)
    row("coal-cpp-align NW, pwalign preset", 1, t)
    methods["NW pwalign preset"] = read_pairs_fasta(o)
    for v in a.variant:
        name, path = v.split("=", 1)
        t, _ = timed(lambda: nw(path, [], 1, f"nw_{name}"), a.repeats)
        row(f"coal-cpp-align NW, build '{name}'", 1, t)
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
        al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                             open_gap_score=-10, extend_gap_score=-0.5)
        al.end_gap_score = 0
        t, _ = timed(lambda: [al.align(x, y)[0] for x, y in pairs], 1)
        row("Biopython PairwiseAligner (same scoring)", 1, t)
    except ImportError:
        pass
    try:
        import parasail
        t, _ = timed(lambda: [parasail.sg_trace_striped_16(x, y, 10, 1, parasail.blosum62).traceback
                              for x, y in pairs], a.repeats)
        row("parasail sg_trace_striped_16 (open 10, extend 1)", 1, t)
        t, _ = timed(lambda: [parasail.sg_trace_scan_16(x, y, 10, 1, parasail.blosum62).traceback
                              for x, y in pairs], a.repeats)
        row("parasail sg_trace_scan_16 (open 10, extend 1)", 1, t)
    except ImportError:
        pass
    t, o = timed(lambda: nw(a.coal_cpp_align, ["--gap-open", "11", "--gap-extend", "1"], 1, "nw_diamondgaps"),
                 a.repeats)
    row("coal-cpp-align NW, DIAMOND's gaps (open 11, extend 1, free ends)", 1, t)
    methods["NW open 11 / extend 1"] = read_pairs_fasta(o)
    from coal import msa
    for tool, exe in (("clustalo", a.clustalo), ("mafft", a.mafft)):
        if not exe:
            continue
        clusters = [[("a", x), ("b", y)] for x, y in pairs]
        jobs = a.jobs[-1]  # one process per pair: slow, timed once with the most jobs
        t, res = timed(lambda: msa.align_clusters(clusters, tool=tool, jobs=jobs, threads=1, pairs_with="msa",
                                                  exe=exe), 1)
        row(f"{tool} {msa.tool_version(exe)}, one process per pair", jobs or os.cpu_count(), t)
        methods[f"{tool} (pair as a cluster)"] = [(r.rows[0][1], r.rows[1][1]) for r in res.values()]
    for fl in ("gap", "align"):
        o = os.path.join(work, f"diamond_{fl}.aln.fasta")
        t, _ = timed(lambda: run([a.coal_cpp_align, fa, fb, "--diamond", hits_tsv, "--flanks", fl, "-o", o,
                                  "-j", "1"]), a.repeats)
        row(f"DIAMOND hit expanded, flanks={fl}", 1, t)
        got = read_pairs_fasta(o)
        order = [(h.qseqid, h.sseqid) for h in A.read_diamond_hits(hits_tsv)]
        m = dict(zip(order, got))
        methods[f"DIAMOND gapped, flanks={fl}"] = [m[p] for p in rbh]
    say(f"\nDIAMOND blastp both ways ({a.threads} threads), RBH input: {t_dmd_plain:.2f} s with "
        f"`qseqid sseqid bitscore`, {t_dmd_gapped:.2f} s with the gapped columns added.\n")

    say("## dN/dS (NG86, Biopython) against the needle-default NW\n")
    say("Codon alignments by coal-cpp; codons with a gap or stop in either sequence dropped. r = Pearson over pairs "
        "with finite values in both; |dw| = |dN/dS - reference dN/dS|.\n")
    res = {}
    for name, alns in methods.items():
        import re
        rows = codon_align(a.coal_cpp, alns, cds_a, cds_b, work, re.sub(r"\W+", "_", name).strip("_"))
        vals = []
        for r in rows:
            vals.append(ng86(r) if r else (math.nan, math.nan, 0))
        res[name] = vals
    ref = res["NW needle (reference)"]
    say("| Method | codons compared (sum) | r dN | r dS | r dN/dS | median abs dw | pairs abs dw > 0.05 | NA dS |")
    say("|---|---:|---:|---:|---:|---:|---:|---:|")
    for name, vals in res.items():
        dn = [v[0] for v in vals]
        ds = [v[1] for v in vals]
        w = [v[0] / v[1] if v[1] and v[1] == v[1] and v[0] == v[0] else math.nan for v in vals]
        rw = [v[0] / v[1] if v[1] and v[1] == v[1] and v[0] == v[0] else math.nan for v in ref]
        dw = [abs(x - y) for x, y in zip(w, rw) if x == x and y == y]
        say(f"| {name} | {sum(v[2] for v in vals):,} | {pearson(dn, [v[0] for v in ref])[0]:.4f} | "
            f"{pearson(ds, [v[1] for v in ref])[0]:.4f} | {pearson(w, rw)[0]:.4f} | "
            f"{statistics.median(dw) if dw else math.nan:.4f} | {sum(d > 0.05 for d in dw)} | "
            f"{sum(1 for x in ds if x != x)} |")
    # where the DIAMOND route differs: by hit coverage of the shorter protein
    cov = [min(h.qend - h.qstart + 1, h.send - h.sstart + 1) / min(len(prot_a[q]), len(prot_b[s]))
           for (q, s), h in ((p, hits[p]) for p in rbh)]
    say("\nBy DIAMOND hit coverage of the shorter protein (r dN/dS against the reference; pairs):\n")
    say("| Method | coverage >= 0.9 | coverage < 0.9 |")
    say("|---|---:|---:|")
    rw = [v[0] / v[1] if v[1] and v[1] == v[1] and v[0] == v[0] else math.nan for v in ref]
    for name, vals in res.items():
        w = [v[0] / v[1] if v[1] and v[1] == v[1] and v[0] == v[0] else math.nan for v in vals]
        cells_ = []
        for keep in (lambda c: c >= 0.9, lambda c: c < 0.9):
            idx = [k for k, c in enumerate(cov) if keep(c)]
            r, n = pearson([w[k] for k in idx], [rw[k] for k in idx])
            cells_.append(f"{r:.4f} ({n})")
        say(f"| {name} | {cells_[0]} | {cells_[1]} |")
    with open(os.path.join(a.out, "results.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
