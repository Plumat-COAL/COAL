"""Simulate pairs of coding sequences with a known omega (dN/dS) under a Goldman-Yang codon model, as AXT.

Codon model (GY94 / M0): rate i->j (one nucleotide apart, no stop codons) = pi_j * (kappa if a transition) *
(omega if nonsynonymous), scaled to one expected nucleotide substitution per codon per unit time; codon
frequencies F3x4 from position-specific base frequencies. The first sequence is drawn from pi, the second is the
first evolved for time t (P = exp(Qt), numpy eigendecomposition). Under this model the YN00 / ML definitions of
dN/dS equal omega, so the estimators should recover it on average (docs/kaks.md "Known omega").

  python benchmarks/kaks/simulate_pairs.py --out sim.axt --codons 500 --reps 50 \
      --omega 0.1 0.5 1 2 --t 0.1 0.3 0.6 --kappa 2 --seed 1
Names: sim_w<omega>_t<t>_k<kappa>_<rep>&b
"""
from __future__ import annotations

import argparse

import numpy as np

BASES = "TCAG"
CODE1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
# position-specific base frequencies (T C A G), unequal on purpose
F3X4 = np.array([[0.20, 0.22, 0.28, 0.30], [0.30, 0.22, 0.30, 0.18], [0.30, 0.25, 0.20, 0.25]])


def codon_model(kappa: float, omega: float):
    sense = [c for c in range(64) if CODE1[c] != "*"]
    pi = np.array([F3X4[0][c // 16] * F3X4[1][(c // 4) % 4] * F3X4[2][c % 4] for c in sense])
    pi /= pi.sum()
    n = len(sense)
    Q = np.zeros((n, n))
    for a, i in enumerate(sense):
        for b, j in enumerate(sense):
            if a == b:
                continue
            bi = [i // 16, (i // 4) % 4, i % 4]
            bj = [j // 16, (j // 4) % 4, j % 4]
            d = [k for k in range(3) if bi[k] != bj[k]]
            if len(d) != 1:
                continue
            x, y = bi[d[0]], bj[d[0]]
            q = pi[b]
            if x + y in (1, 5):
                q *= kappa
            if CODE1[i] != CODE1[j]:
                q *= omega
            Q[a, b] = q
    np.fill_diagonal(Q, -Q.sum(axis=1))
    # one expected nucleotide substitution per codon per unit time
    Q /= -(pi * np.diag(Q)).sum()
    return sense, pi, Q


def pmatrix(pi, Q, t):
    s = np.sqrt(pi)
    A = (Q * s[:, None]) / s[None, :]
    A = (A + A.T) / 2
    w, V = np.linalg.eigh(A)
    P = (V * np.exp(w * t)) @ V.T
    P = P / s[:, None] * s[None, :]
    P = np.clip(P, 0, None)
    return P / P.sum(axis=1, keepdims=True)


def codon_str(c: int) -> str:
    return BASES[c // 16] + BASES[(c // 4) % 4] + BASES[c % 4]


def simulate(rng, ncod: int, kappa: float, omega: float, t: float):
    sense, pi, Q = codon_model(kappa, omega)
    P = pmatrix(pi, Q, t)
    a = rng.choice(len(sense), size=ncod, p=pi)
    cum = np.cumsum(P, axis=1)
    u = rng.random(ncod)
    b = np.array([min(int(np.searchsorted(cum[x], u[k])), len(sense) - 1) for k, x in enumerate(a)])
    s1 = "".join(codon_str(sense[x]) for x in a)
    s2 = "".join(codon_str(sense[x]) for x in b)
    return s1, s2


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--codons", type=int, default=500)
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--omega", type=float, nargs="+", default=[0.1, 0.5, 1.0, 2.0])
    ap.add_argument("--t", type=float, nargs="+", default=[0.1, 0.3, 0.6])
    ap.add_argument("--kappa", type=float, default=2.0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    with open(a.out, "w") as f:
        for w in a.omega:
            for t in a.t:
                for r in range(a.reps):
                    s1, s2 = simulate(rng, a.codons, a.kappa, w, t)
                    f.write(f"sim_w{w:g}_t{t:g}_k{a.kappa:g}_{r}&b\n{s1}\n{s2}\n\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
