"""Check a `coal pipeline pairs` run against the same steps done by hand (docs/pipeline.md "Validation").

For the pairs of RUN/orthologs.tsv: translate in Python (an independent translation, standard code) ->
`coal-cpp-align` (the pairwise aligner as a separate program) -> `coal-cpp` (the codon engine, one run per pair) ->
`coal-cpp kaks` and `coal-cpp dnds` on the resulting codon alignments, and, when given, KaKs_Calculator 2.0 on the
same AXT file. Compared with RUN/dnds.tsv, RUN/codon_alignments.fasta and RUN/dnds.kaks.

  python benchmarks/pipeline/validate_pairs.py --query q_cds.fa --subject s_cds.fa --run RUN \
      --coal-cpp build/coal-cpp --coal-cpp-align build/coal-cpp-align [--kaks-calculator KaKs_Calculator] [-m YN]

The run must have been made with --keep-alignments --kaks-format and the default aligner (nw).
"""
from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
import tempfile

BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AAS[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
        for k, c in enumerate(BASES)}


def read_fasta(path):
    """{id: sequence}; of a repeated ID the first record, as the pipeline."""
    recs, cur = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = (line[1:].split() or [""])[0]
                cur = None if name in recs else recs.setdefault(name, [])
            elif cur is not None:
                cur.append(line.strip())
    return {k: "".join(v) for k, v in recs.items()}


def read_fasta_list(path):
    out, name = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out.append([name, ""])
            elif name is not None:
                out[-1][1] += line.strip()
    return out


def translate(cds):
    cds = cds.upper().replace("U", "T")
    cds = cds[:len(cds) - len(cds) % 3]
    p = "".join(CODE.get(cds[i:i + 3], "X") for i in range(0, len(cds), 3))
    return p[:-1] if p.endswith("*") else p


def num(x):
    return math.nan if x == "NA" else float(x)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True)
    ap.add_argument("--subject", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--coal-cpp", required=True)
    ap.add_argument("--coal-cpp-align", required=True)
    ap.add_argument("--kaks-calculator")
    ap.add_argument("-m", "--method", action="append")
    a = ap.parse_args()
    methods = a.method or ["YN"]
    q, s = read_fasta(a.query), read_fasta(a.subject)
    pairs = [line.split("\t")[:2] for line in open(os.path.join(a.run, "orthologs.tsv")).read().splitlines()[1:]]
    ours = {}
    for line in open(os.path.join(a.run, "dnds.tsv")).read().splitlines()[1:]:
        f = line.split("\t")
        ours[(f[0], f[1], f[2])] = f
    bad = 0
    with tempfile.TemporaryDirectory() as td:
        A, B = os.path.join(td, "a.faa"), os.path.join(td, "b.faa")
        with open(A, "w") as fa, open(B, "w") as fb:
            for x, y in pairs:
                fa.write(f">{x}\n{translate(q[x])}\n")
                fb.write(f">{y}\n{translate(s[y])}\n")
        aln = os.path.join(td, "aln.fasta")
        subprocess.run([a.coal_cpp_align, A, B, "-o", aln], check=True)
        rows = read_fasta_list(aln)
        assert len(rows) == 2 * len(pairs)
        kept = read_fasta_list(os.path.join(a.run, "codon_alignments.fasta"))
        kept_by = {(kept[i][0], kept[i + 1][0]): (kept[i][1], kept[i + 1][1]) for i in range(0, len(kept), 2)}
        axt = os.path.join(td, "pairs.axt")
        n_codon_same = 0
        with open(axt, "w") as fx:
            for k, (x, y) in enumerate(pairs):
                pep, nuc = os.path.join(td, "pep.fa"), os.path.join(td, "nuc.fa")
                with open(pep, "w") as f:
                    f.write(f">a\n{rows[2 * k][1]}\n>b\n{rows[2 * k + 1][1]}\n")
                with open(nuc, "w") as f:
                    f.write(f">a\n{q[x]}\n>b\n{s[y]}\n")
                r = subprocess.run([a.coal_cpp, pep, nuc, "-q"], capture_output=True, text=True)
                if r.returncode != 0:
                    print(f"by hand: codon step failed for {x} {y}: {r.stderr.strip()}")
                    bad += 1
                    continue
                parts = r.stdout.split(">")[1:]
                r1, r2 = ("".join(p.splitlines()[1:]) for p in parts)
                fx.write(f"{x}&{y}\n{r1}\n{r2}\n\n")
                if kept_by.get((x, y)) == (r1, r2):
                    n_codon_same += 1
                else:
                    bad += 1
                    print(f"codon alignment differs: {x} {y}")
        print(f"{len(pairs)} pairs; codon alignments identical to the by-hand ones: {n_codon_same}")
        margs = [x for m in methods for x in ("-m", m)]
        tsv, kk = os.path.join(td, "hand.tsv"), os.path.join(td, "hand.kaks")
        subprocess.run([a.coal_cpp, "dnds", axt, "-o", tsv, "-q", *margs], check=True)
        subprocess.run([a.coal_cpp, "kaks", "-i", axt, "-o", kk, "-q", *margs], check=True)
        n_same = n_rows = 0
        for line in open(tsv).read().splitlines()[1:]:
            f = line.split("\t")
            o = ours.get((f[0], f[1], f[2]))
            n_rows += 1
            # by hand: seq1 seq2 method dN dS dNdS SE_dN SE_dS N_sites S_sites Nd Sd codons_compared status note
            # the pipeline adds what it knows about the CDS (length not a multiple of 3) to status and note
            cds_note = o is not None and o[15] == "warning" and f[13] == "ok" and o[16].startswith("seq")
            if o is not None and o[3:13] == f[3:13] and (o[15] == f[13] or cds_note):
                n_same += 1
            else:
                bad += 1
                print("differs:", f, o)
        print(f"dnds.tsv rows equal to `coal dnds` on the by-hand codon alignments (text of dN dS dNdS SEs sites "
              f"differences codons status): {n_same} of {n_rows}")
        same_kaks = open(kk).read() == open(os.path.join(a.run, "dnds.kaks")).read()
        print(f"dnds.kaks byte-identical to `coal kaks` on the by-hand codon alignments: {same_kaks}")
        bad += 0 if same_kaks else 1
        if a.kaks_calculator:
            for m in methods:
                out = os.path.join(td, f"kc_{m}.kaks")
                subprocess.run([a.kaks_calculator, "-i", axt, "-o", out, "-m", m], check=True,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                n = both = 0
                mx = [0.0, 0.0, 0.0]
                na_ours, na_theirs, over = [], [], 0
                for line in open(out).read().splitlines()[1:]:
                    f = line.split("\t")
                    x, y = f[0].split("&")
                    o = ours[(x, y, m)]
                    n += 1
                    for k, (mine, theirs) in enumerate(((o[3], f[2]), (o[4], f[3]), (o[5], f[4]))):
                        u, v = num(mine), num(theirs)
                        if math.isnan(v) and (math.isnan(u) or u < 1e-6):
                            continue   # KaKs_Calculator prints values below 1e-6 as NA
                        if math.isnan(u) and not math.isnan(v):
                            na_ours.append((x, y, "Ka Ks Ka/Ks".split()[k], theirs, o[15]))
                            continue
                        if math.isnan(v):
                            na_theirs.append((x, y, k, mine))
                            continue
                        d = abs(u - v)
                        tol = 5.1e-6 * max(abs(v), 1e-300)   # KaKs_Calculator prints 6 significant digits
                        if d > tol:
                            over += 1
                            print("KaKs_Calculator differs:", x, y, m, mine, theirs)
                        mx[k] = max(mx[k], d / max(abs(v), 1e-300))
                    both += 1
                print(f"KaKs_Calculator -m {m}: {n} pairs; max relative difference Ka {mx[0]:.1e}, Ks {mx[1]:.1e}, "
                      f"Ka/Ks {mx[2]:.1e} (its 6 printed digits allow 5e-6); beyond that: {over}; "
                      f"NA in COAL only: {len(na_ours)} {na_ours[:6]}; NA in KaKs_Calculator only: {len(na_theirs)}")
                bad += over + len(na_theirs)
    print("OK" if not bad else f"{bad} difference(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
