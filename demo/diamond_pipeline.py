"""Ortholog codon alignments for two species in one go: DIAMOND reciprocal best hits -> COAL.

    python demo/diamond_pipeline.py \
        --a-protein human.pep.fa --a-cds human.cds.fa --a-name human \
        --b-protein mouse.pep.fa --b-cds mouse.cds.fa --b-name mouse \
        --out orthologs/

Steps
  1. Keep the longest protein isoform per gene (gene from a 'gene:' / 'gene=' tag in the FASTA header,
     else the ID without an isoform suffix such as .1, -RA, -PA, _i1, .p1, -T1).
  2. Find reciprocal best hits with YOUR local `diamond` binary (blastp both ways).
  3. Align each pair of proteins (built-in Needleman-Wunsch, BLOSUM62, affine gaps; pure Python).
  4. Write ONE batched protein alignment + DNA file (pairs padded with gaps to a common length) and run
     COAL once on the batch, then split the result into one codon alignment per pair.

Outputs in --out: rbh.tsv, batch.pep.fasta, batch.cds.fasta, batch.codon.fasta, pairs/<A>__<B>.codon.fasta,
and coal.stderr.txt. DIAMOND is not included: install it from https://github.com/bbuchfink/diamond
(or `conda install -c bioconda diamond`).
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

# ---------------------------------------------------------------- FASTA


def read_fasta(path):
    recs, name, header, buf = [], None, "", []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, header, "".join(buf)))
                header = line[1:]
                name, buf = (header.split()[0] if header.split() else ""), []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, header, "".join(buf)))
    return recs


def write_fasta(path, recs):
    with open(path, "w", newline="\n") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


# ---------------------------------------------------------------- 1. longest isoform

_SUFFIX = re.compile(r"(\.\d+|-R[A-Z]+|-P[A-Z]+|_i\d+|\.p\d+|-T\d+|\.t\d+)$")


def gene_of(seq_id: str, header: str) -> str:
    m = re.search(r"\bgene[:=]([^\s;]+)", header)
    return m.group(1) if m else _SUFFIX.sub("", seq_id)


def longest_isoforms(proteins, cds_ids):
    """{gene: (transcript_id, protein)} keeping the longest protein that has a CDS."""
    best = {}
    for tid, header, seq in proteins:
        if tid not in cds_ids:
            continue
        seq = seq.rstrip("*")
        g = gene_of(tid, header)
        if g not in best or len(seq) > len(best[g][1]):
            best[g] = (tid, seq)
    return best


# ---------------------------------------------------------------- 2. DIAMOND RBH

def diamond_cmd(diamond: str) -> list[str]:
    # a .py stand-in (used in tests) runs with this Python
    return [sys.executable, diamond] if diamond.endswith(".py") else [diamond]


def find_diamond(diamond: str | None) -> str | None:
    if diamond:
        return diamond if (os.path.isfile(diamond) or shutil.which(diamond)) else None
    return shutil.which("diamond")


def best_hits(diamond, query, db_fasta, work, tag, threads, timeout):
    db = os.path.join(work, f"{tag}_db")
    out = os.path.join(work, f"{tag}.tsv")
    run = lambda args: subprocess.run(diamond_cmd(diamond) + args, capture_output=True, text=True,  # noqa: E731
                                      timeout=timeout)
    r = run(["makedb", "--in", db_fasta, "-d", db, "--quiet"])
    if r.returncode != 0:
        raise RuntimeError(f"diamond makedb failed: {r.stderr.strip()[:300]}")
    r = run(["blastp", "-q", query, "-d", db, "-o", out, "--outfmt", "6", "qseqid", "sseqid", "bitscore",
             "evalue", "--max-target-seqs", "5", "--evalue", "1e-5", "--threads", str(threads), "--quiet"])
    if r.returncode != 0:
        raise RuntimeError(f"diamond blastp failed: {r.stderr.strip()[:300]}")
    best: dict[str, tuple[str, float]] = {}
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            q, s, score = f[0], f[1], float(f[2])
            if q not in best or score > best[q][1]:
                best[q] = (s, score)
    return {q: s for q, (s, _) in best.items()}


def reciprocal_best_hits(ab: dict, ba: dict) -> list[tuple[str, str]]:
    return sorted((a, b) for a, b in ab.items() if ba.get(b) == a)


# ---------------------------------------------------------------- 3. Needleman-Wunsch (Gotoh)

_B62_ORDER = "ARNDCQEGHILKMFPSTWYVBZX*"
_B62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""
_IDX = {a: i for i, a in enumerate(_B62_ORDER)}
_B62 = [[int(x) for x in row.split()] for row in _B62_ROWS.strip().splitlines()]


def _score(a, b):
    return _B62[_IDX.get(a, 22)][_IDX.get(b, 22)]


def needleman_wunsch(s: str, t: str, gap_open: int = 11, gap_extend: int = 1) -> tuple[str, str]:
    """Global alignment, BLOSUM62, affine gaps (Gotoh: a gap of length k costs open + (k-1)*extend).
    Returns the two gapped rows."""
    s, t = s.upper(), t.upper()
    n, m = len(s), len(t)
    NEG = -10 ** 9
    # three states: 0 = s[i] aligned to t[j]; 1 = s[i] against a gap; 2 = t[j] against a gap
    V = [[[NEG] * (m + 1) for _ in range(n + 1)] for _ in range(3)]
    P = [[[0] * (m + 1) for _ in range(n + 1)] for _ in range(3)]  # previous state for traceback
    V[0][0][0] = 0
    for i in range(1, n + 1):
        V[1][i][0] = -gap_open - (i - 1) * gap_extend
        P[1][i][0] = 0 if i == 1 else 1
    for j in range(1, m + 1):
        V[2][0][j] = -gap_open - (j - 1) * gap_extend
        P[2][0][j] = 0 if j == 1 else 2
    M, X, Y = V
    PM, PX, PY = P
    for i in range(1, n + 1):
        si = s[i - 1]
        for j in range(1, m + 1):
            c = (M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1])
            k = 0 if c[0] >= c[1] and c[0] >= c[2] else (1 if c[1] >= c[2] else 2)
            M[i][j] = c[k] + _score(si, t[j - 1])
            PM[i][j] = k
            c = (M[i - 1][j] - gap_open, X[i - 1][j] - gap_extend, Y[i - 1][j] - gap_open)
            k = 0 if c[0] >= c[1] and c[0] >= c[2] else (1 if c[1] >= c[2] else 2)
            X[i][j] = c[k]
            PX[i][j] = k
            c = (M[i][j - 1] - gap_open, X[i][j - 1] - gap_open, Y[i][j - 1] - gap_extend)
            k = 0 if c[0] >= c[1] and c[0] >= c[2] else (1 if c[1] >= c[2] else 2)
            Y[i][j] = c[k]
            PY[i][j] = k
    a, b = [], []
    i, j = n, m
    state = max(range(3), key=lambda k: V[k][n][m])
    while i > 0 or j > 0:
        prev = P[state][i][j]
        if state == 0:
            a.append(s[i - 1]); b.append(t[j - 1]); i -= 1; j -= 1  # noqa: E702
        elif state == 1:
            a.append(s[i - 1]); b.append("-"); i -= 1  # noqa: E702
        else:
            a.append("-"); b.append(t[j - 1]); j -= 1  # noqa: E702
        state = prev
    return "".join(reversed(a)), "".join(reversed(b))


# ---------------------------------------------------------------- 4. batch + COAL

def coal_help() -> str:
    try:
        return subprocess.run([sys.executable, "-m", "coal", "--help"], capture_output=True, text=True,
                              timeout=60).stdout
    except (OSError, subprocess.SubprocessError):
        return ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for s in ("a", "b"):
        ap.add_argument(f"--{s}-protein", required=True, help=f"species {s.upper()} proteins (FASTA)")
        ap.add_argument(f"--{s}-cds", required=True, help=f"species {s.upper()} coding DNA (FASTA, same IDs)")
        ap.add_argument(f"--{s}-name", default=s.upper(), help=f"short name for species {s.upper()}")
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--diamond", help="path to the diamond binary (default: `diamond` on PATH)")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--codontable", type=int, default=1, help="NCBI genetic code for COAL (default 1)")
    ap.add_argument("--max-length", type=int, default=5000,
                    help="skip pairs with a protein longer than this (pure-Python alignment is O(n*m))")
    ap.add_argument("--timeout", type=float, default=3600, help="time limit for each DIAMOND / COAL step (s)")
    a = ap.parse_args(argv)

    diamond = find_diamond(a.diamond)
    if diamond is None:
        print("DIAMOND was not found" + (f" at {a.diamond}" if a.diamond else " on your PATH") + ".\n"
              "Install it (https://github.com/bbuchfink/diamond/releases, or `conda install -c bioconda "
              "diamond`) and run again, or give its path with --diamond /path/to/diamond.", file=sys.stderr)
        return 1
    if "__" in a.a_name + a.b_name:
        ap.error("species names must not contain '__'")
    os.makedirs(os.path.join(a.out, "pairs"), exist_ok=True)

    species = {}
    for s, name in (("a", a.a_name), ("b", a.b_name)):
        prot = read_fasta(getattr(a, f"{s}_protein"))
        cds = {i: seq for i, _, seq in read_fasta(getattr(a, f"{s}_cds"))}
        iso = longest_isoforms(prot, cds)
        species[s] = {"name": name, "iso": iso, "cds": cds}
        print(f"{name}: {len(prot)} proteins -> {len(iso)} genes (longest isoform with a CDS)")

    with tempfile.TemporaryDirectory(prefix="coal_rbh_") as work:
        fa = {}
        for s in ("a", "b"):
            fa[s] = os.path.join(work, f"{s}.faa")
            write_fasta(fa[s], [(tid, seq) for tid, seq in species[s]["iso"].values()])
        try:
            ab = best_hits(diamond, fa["a"], fa["b"], work, "ab", a.threads, a.timeout)
            ba = best_hits(diamond, fa["b"], fa["a"], work, "ba", a.threads, a.timeout)
        except (RuntimeError, OSError, subprocess.SubprocessError) as e:
            print(f"error: {e}", file=sys.stderr)
            return 1
    pairs = reciprocal_best_hits(ab, ba)
    print(f"reciprocal best hits: {len(pairs)}")
    with open(os.path.join(a.out, "rbh.tsv"), "w", newline="\n") as fh:
        fh.write(f"{a.a_name}\t{a.b_name}\n")
        for x, y in pairs:
            fh.write(f"{x}\t{y}\n")
    if not pairs:
        print("No reciprocal best hits, nothing to align.")
        return 0

    pa = {tid: seq for tid, seq in species["a"]["iso"].values()}
    pb = {tid: seq for tid, seq in species["b"]["iso"].values()}
    aligned, skipped = [], []
    for x, y in pairs:
        if max(len(pa[x]), len(pb[y])) > a.max_length:
            skipped.append((x, y))
            continue
        ra, rb = needleman_wunsch(pa[x], pb[y])
        aligned.append((x, y, ra, rb))
    if skipped:
        print(f"skipped {len(skipped)} pairs longer than --max-length {a.max_length}")

    # one batch: every pair padded with trailing gaps to the same length, unique IDs per row
    width = max(len(ra) for _, _, ra, _ in aligned)
    pep, nuc, ids = [], [], []
    for x, y, ra, rb in aligned:
        ia, ib = f"{a.a_name}__{x}", f"{a.b_name}__{y}"
        pep += [(ia, ra.ljust(width, "-")), (ib, rb.ljust(width, "-"))]
        nuc += [(ia, species["a"]["cds"][x]), (ib, species["b"]["cds"][y])]
        ids.append((ia, ib, len(ra)))
    batch_pep = os.path.join(a.out, "batch.pep.fasta")
    batch_cds = os.path.join(a.out, "batch.cds.fasta")
    batch_out = os.path.join(a.out, "batch.codon.fasta")
    write_fasta(batch_pep, pep)
    write_fasta(batch_cds, nuc)

    cmd = [sys.executable, "-m", "coal", batch_pep, batch_cds, "-codontable", str(a.codontable), "-o", batch_out]
    h = coal_help()
    if "--keep-going" in h:
        cmd.append("--keep-going")
    if "--report" in h:
        cmd += ["--report", os.path.join(a.out, "batch.report.json")]
    print(f"running COAL once on {len(pep)} sequences ({len(aligned)} pairs)")
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=a.timeout)
    except subprocess.TimeoutExpired:
        print("error: COAL did not finish in time", file=sys.stderr)
        return 1
    with open(os.path.join(a.out, "coal.stderr.txt"), "w") as fh:
        fh.write(r.stderr)
    if r.returncode not in (0, 2) or not os.path.exists(batch_out):
        print(f"COAL failed (exit {r.returncode}):\n{r.stderr.strip()[-1500:]}", file=sys.stderr)
        print("Tip: open the batch in the COAL GUI or run `python -m coal_app.viewer "
              f"{batch_pep} --dna {batch_cds} -o view.html` to see which pair breaks.", file=sys.stderr)
        return 1
    codon = {i: seq for i, _, seq in read_fasta(batch_out)}
    written = 0
    for ia, ib, n in ids:
        if ia not in codon or ib not in codon:
            continue  # failed pair (with --keep-going)
        name = f"{ia}__{ib.split('__', 1)[1]}".replace("/", "_")
        write_fasta(os.path.join(a.out, "pairs", f"{name}.codon.fasta"),
                    [(ia, codon[ia][:3 * n]), (ib, codon[ib][:3 * n])])
        written += 1
    warn = sum(1 for line in r.stderr.splitlines() if "warning" in line)
    print(f"done: {written} pairwise codon alignments in {os.path.join(a.out, 'pairs')}"
          + (f" ({warn} COAL warnings, see coal.stderr.txt)" if warn else ""))
    return 0 if r.returncode == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
