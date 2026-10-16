#!/usr/bin/env python3
"""Generate a synthetic, correct-by-construction benchmark dataset.

Model: one random ancestral coding sequence of L codons. Each of N descendants
copies it, mutates ~10% of codons (any sense codon), and deletes ~5% of codons.
The protein alignment has L columns with '-' at deleted positions; the DNA file
holds the ungapped codons (no stop codons, no UTRs -- within COAL v0.1 scope).

Usage: python benchmarks/gen_data.py N_SEQS N_CODONS OUT_DIR [--seed S]
Writes OUT_DIR/pep.fasta and OUT_DIR/nuc.fasta (sequences on one line each).
"""
from __future__ import annotations

import argparse
import random
from itertools import product
from pathlib import Path

_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {"".join(c): _AA[i] for i, c in enumerate(product(_B, repeat=3))}
SENSE = [c for c, a in CODE.items() if a != "*"]


def generate(n_seqs: int, n_codons: int, out: Path, seed: int = 1, p_mut=0.10, p_del=0.05):
    rng = random.Random(seed)
    anc = [rng.choice(SENSE) for _ in range(n_codons)]
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "pep.fasta", "w", newline="\n") as pf, open(out / "nuc.fasta", "w", newline="\n") as nf:
        for i in range(n_seqs):
            aln, dna = [], []
            for c in anc:
                r = rng.random()
                if r < p_del:
                    aln.append("-")
                    continue
                if r < p_del + p_mut:
                    c = rng.choice(SENSE)
                aln.append(CODE[c])
                dna.append(c)
            pf.write(f">seq{i}\n{''.join(aln)}\n")
            nf.write(f">seq{i}\n{''.join(dna)}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("n_seqs", type=int)
    ap.add_argument("n_codons", type=int)
    ap.add_argument("out", type=Path)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    generate(a.n_seqs, a.n_codons, a.out, a.seed)


if __name__ == "__main__":
    main()
