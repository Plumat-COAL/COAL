"""Codon-pair simulator for the SIM-VALIDATE study (GY94 / M0 with F3x4 codon frequencies).

Model: instantaneous rate i->j (codons one nucleotide apart, no stop codons)
    q_ij = pi_j * (kappa if transition) * (omega if nonsynonymous)
scaled so that the expected number of nucleotide substitutions per codon per unit time is 1
(t = substitutions per codon, the PAML convention).

True dN and dS follow Yang & Nielsen (2000) / Yang (2006, sec. 2.5.3), the definition that YN00 and codeml
estimate: rho_S = proportion of substitutions that are synonymous under the model; S = 3 * rho_S1 where rho_S1 is
that proportion under the same model with omega = 1 (sites are defined by the mutation process: kappa and codon
frequencies); N = 3 - S. Then dS = t*rho_S/S and dN = t*(1-rho_S)/N, so dN/dS = omega exactly.

Two samplers:
  sample_pair_P   : root codon ~ pi, descendant ~ P(t)[root]  (fast; used for the study)
  sample_pair_path: Gillespie path simulation; also returns realised synonymous / nonsynonymous counts
                    (used only to validate the fast sampler and the dN/dS definitions)
"""
from __future__ import annotations

import numpy as np

BASES = "TCAG"
CODE1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
SENSE = [c for c in range(64) if CODE1[c] != "*"]
CODON_STR = [BASES[c // 16] + BASES[(c // 4) % 4] + BASES[c % 4] for c in SENSE]

FREQS = {
    # position-specific base frequencies (T, C, A, G) for codon positions 1, 2, 3
    "equal": np.full((3, 4), 0.25),
    # GC3-rich skew (GC3 = 0.80), unequal at positions 1 and 2 as well
    "skewed": np.array([[0.20, 0.20, 0.30, 0.30], [0.30, 0.25, 0.30, 0.15], [0.10, 0.40, 0.10, 0.40]]),
}


def is_transition(x: int, y: int) -> bool:
    # T=0 C=1 A=2 G=3: T<->C and A<->G
    return {x, y} in ({0, 1}, {2, 3})


class CodonModel:
    def __init__(self, kappa: float, omega: float, f3x4: np.ndarray):
        self.kappa, self.omega = kappa, omega
        pi = np.array([f3x4[0][c // 16] * f3x4[1][(c // 4) % 4] * f3x4[2][c % 4] for c in SENSE])
        self.pi = pi / pi.sum()
        self.n = len(SENSE)
        self.Q, self.syn_mask = self._build(omega)
        Q1, _ = self._build(1.0)
        self.rhoS = float((self.pi[:, None] * self.Q * self.syn_mask).sum())  # Q scaled: total rate 1
        self.rhoS1 = float((self.pi[:, None] * Q1 * self.syn_mask).sum())
        self.S = 3 * self.rhoS1  # synonymous sites per codon
        self.N = 3 - self.S
        self._eig()

    def _build(self, omega):
        n = self.n
        Q = np.zeros((n, n))
        syn = np.zeros((n, n), dtype=bool)
        for a, i in enumerate(SENSE):
            bi = (i // 16, (i // 4) % 4, i % 4)
            for b, j in enumerate(SENSE):
                if a == b:
                    continue
                bj = (j // 16, (j // 4) % 4, j % 4)
                d = [k for k in range(3) if bi[k] != bj[k]]
                if len(d) != 1:
                    continue
                q = self.pi[b]
                if is_transition(bi[d[0]], bj[d[0]]):
                    q *= self.kappa
                if CODE1[i] == CODE1[j]:
                    syn[a, b] = True
                else:
                    q *= omega
                Q[a, b] = q
        np.fill_diagonal(Q, -Q.sum(axis=1))
        Q /= -(self.pi * np.diag(Q)).sum()
        return Q, syn

    def _eig(self):
        s = np.sqrt(self.pi)
        A = (self.Q * s[:, None]) / s[None, :]
        A = (A + A.T) / 2
        self._w, self._V = np.linalg.eigh(A)
        self._s = s

    def P(self, t: float) -> np.ndarray:
        P = (self._V * np.exp(self._w * t)) @ self._V.T
        P = P / self._s[:, None] * self._s[None, :]
        P = np.clip(P, 0, None)
        return P / P.sum(axis=1, keepdims=True)

    # true distances for a time t
    def dS_of_t(self, t):
        return t * self.rhoS / self.S

    def dN_of_t(self, t):
        return t * (1 - self.rhoS) / self.N

    def t_for_dS(self, dS):
        return dS * self.S / self.rhoS


def sample_pair_P(model: CodonModel, t: float, ncod: int, rng: np.random.Generator):
    root = rng.choice(model.n, size=ncod, p=model.pi)
    cum = np.cumsum(model.P(t), axis=1)
    cum[:, -1] = 1.0
    u = rng.random(ncod)
    desc = (cum[root] < u[:, None]).sum(axis=1)
    return root, desc


def sample_pair_path(model: CodonModel, t: float, ncod: int, rng: np.random.Generator):
    """Gillespie simulation along a branch of length t; returns root, descendant, #syn, #nonsyn events."""
    Q = model.Q
    rates = -np.diag(Q)
    jump = np.where(np.eye(model.n, dtype=bool), 0.0, Q) / rates[:, None]
    cumj = np.cumsum(jump, axis=1)
    cumj[:, -1] = 1.0
    root = rng.choice(model.n, size=ncod, p=model.pi)
    cur = root.copy()
    nsyn = np.zeros(ncod, dtype=np.int64)
    nnon = np.zeros(ncod, dtype=np.int64)
    time = rng.exponential(1.0 / rates[cur])
    active = time < t
    while active.any():
        idx = np.nonzero(active)[0]
        c = cur[idx]
        u = rng.random(len(idx))
        nxt = (cumj[c] < u[:, None]).sum(axis=1)
        s = model.syn_mask[c, nxt]
        nsyn[idx] += s
        nnon[idx] += ~s
        cur[idx] = nxt
        time[idx] += rng.exponential(1.0 / rates[nxt])
        active = time < t
    return root, cur, nsyn, nnon


def to_str(codons) -> str:
    return "".join(CODON_STR[c] for c in codons)
