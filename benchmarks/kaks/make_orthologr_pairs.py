"""Build the Ka/Ks validation set from orthologr's example genes (docs/kaks.md "Validation").

A. thaliana and A. lyrata CDS (orthologr inst/seqs/ortho_{thal,lyra}_cds_1000.fasta) -> translate -> DIAMOND
reciprocal best hits -> global protein alignment (Biopython PairwiseAligner, BLOSUM62, gap open 11, extend 1,
full length, no end truncation) -> codon alignment with coal-cpp -> pairs.axt (KaKs_Calculator's input),
pairs.fasta (consecutive pairs) and pairs.tsv (ids). The sequences are not part of this repository; the script
rebuilds the set from an orthologr checkout.

  python benchmarks/kaks/make_orthologr_pairs.py --orthologr ../orthologr --diamond diamond \
      --coal-cpp build/coal-cpp --out /tmp/kaks-val
"""
from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq


def read_cds(path: Path) -> dict[str, str]:
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(str(path), "fasta")}


def translate(cds: str) -> str:
    cds = cds[: len(cds) - len(cds) % 3]
    p = str(Seq(cds).translate())
    return p[:-1] if p.endswith("*") else p


def best_hits(diamond: str, query: Path, db: Path, work: Path, threads: int) -> dict[str, str]:
    out = work / f"{query.stem}_vs_{db.stem}.tsv"
    subprocess.run([diamond, "makedb", "--in", str(db), "-d", str(work / db.stem), "--quiet"], check=True)
    subprocess.run([diamond, "blastp", "-q", str(query), "-d", str(work / db.stem), "-o", str(out),
                    "--max-target-seqs", "1", "--evalue", "1e-5", "--quiet", "-p", str(threads),
                    "--outfmt", "6", "qseqid", "sseqid", "bitscore"], check=True)
    best: dict[str, tuple[float, str]] = {}
    for line in out.read_text().splitlines():
        q, s, b = line.split("\t")
        if q not in best or float(b) > best[q][0]:
            best[q] = (float(b), s)
    return {q: s for q, (b, s) in best.items()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthologr", required=True)
    ap.add_argument("--diamond", default="diamond")
    ap.add_argument("--coal-cpp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args()
    seqs = Path(a.orthologr) / "inst" / "seqs"
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    thal = read_cds(seqs / "ortho_thal_cds_1000.fasta")
    lyra = read_cds(seqs / "ortho_lyra_cds_1000.fasta")
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        for name, d in (("thal", thal), ("lyra", lyra)):
            with open(work / f"{name}.faa", "w") as f:
                for k, v in d.items():
                    f.write(f">{k}\n{translate(v)}\n")
        ab = best_hits(a.diamond, work / "thal.faa", work / "lyra.faa", work, a.threads)
        ba = best_hits(a.diamond, work / "lyra.faa", work / "thal.faa", work, a.threads)
        rbh = [(q, s) for q, s in ab.items() if ba.get(s) == q]
        rbh.sort()
        aligner = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                  open_gap_score=-11, extend_gap_score=-1)
        axt = open(out / "pairs.axt", "w")
        fa = open(out / "pairs.fasta", "w")
        tsv = open(out / "pairs.tsv", "w")
        n = 0
        for q, s in rbh:
            p1, p2 = translate(thal[q]), translate(lyra[s])
            aln = aligner.align(p1, p2)[0]
            pep = work / "pep.fa"
            nuc = work / "nuc.fa"
            pep.write_text(f">a\n{aln[0]}\n>b\n{aln[1]}\n")
            nuc.write_text(f">a\n{thal[q]}\n>b\n{lyra[s]}\n")
            r = subprocess.run([a.coal_cpp, str(pep), str(nuc), "-q"], capture_output=True, text=True)
            if r.returncode != 0:
                print(f"skip {q} {s}: {r.stderr.strip()}")
                continue
            rows = [x for x in r.stdout.split(">")[1:]]
            r1 = "".join(rows[0].splitlines()[1:])
            r2 = "".join(rows[1].splitlines()[1:])
            axt.write(f"{q}&{s}\n{r1}\n{r2}\n\n")
            fa.write(f">{q}\n{r1}\n>{s}\n{r2}\n")
            tsv.write(f"{q}\t{s}\n")
            n += 1
        for f in (axt, fa, tsv):
            f.close()
    print(f"{len(rbh)} reciprocal best hits, {n} codon-aligned pairs -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
