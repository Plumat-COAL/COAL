#!/usr/bin/env python3
"""Simulate codon evolution down a tree, with codon indels, keeping the true codon alignment.

Independent of COAL's own code: genetic codes come from Biopython's NCBI tables.

Model
-----
* Tree: Kingman coalescent with N leaves, scaled to a root-to-tip height H
  (expected substitutions per codon).
* Substitution: Goldman-Yang (GY94) codon model with transition/transversion ratio
  kappa, dN/dS omega and F3x4 codon frequencies (random per family). Gamma rate
  variation across sites (4 discrete categories, shape alpha).
* Indels: whole codons only. Events per branch ~ Poisson(indel_rate * t * length),
  half insertions, half deletions; lengths geometric (mean `indel_mean`, max 30).
  Insertions create new alignment columns, so the true alignment has realistic,
  tree-structured gap patterns. Codon 0 is the start Met (ATG) and is never
  substituted or deleted, so every sequence still starts with a Met.
* "clades" mode, for very large N: a full tree of K clade ancestors as above,
  then each leaf is a short independent branch from its clade ancestor
  (substitutions + deletions only). Much faster; used for the speed tier.

A family is returned as a `Family`: leaf IDs, the ordered alignment columns, and per
leaf a list of codon strings (None = gap). Scenario edits (UTRs, case, frameshifts,
...) are applied afterwards by build_sim.py.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from Bio.Data import CodonTable

BASES = "TCAG"
_TRANSITIONS = {frozenset("AG"), frozenset("CT")}


# ---------------------------------------------------------------- genetic codes
@dataclass
class Code:
    table: int
    sense: list[str]            # sense codons, model state order
    aa: dict[str, str]          # codon -> amino acid (stops -> '*')
    stops: list[str]
    starts: list[str]


_CODES: dict[int, Code] = {}


def code(table: int) -> Code:
    if table not in _CODES:
        t = CodonTable.unambiguous_dna_by_id[table]
        aa = dict(t.forward_table)
        for s in t.stop_codons:
            aa[s] = "*"
        sense = sorted(t.forward_table)
        _CODES[table] = Code(table, sense, aa, list(t.stop_codons), list(t.start_codons))
    return _CODES[table]


def translate(codon: str, table: int) -> str:
    """Translate one (possibly ambiguous, any-case, RNA) codon; 'X' if ambiguous."""
    c = code(table)
    s = codon.upper().replace("U", "T")
    if s in c.aa:
        return c.aa[s]
    expansions = _expand(s)
    aas = {c.aa.get(e) for e in expansions}
    return aas.pop() if len(aas) == 1 and None not in aas else "X"


_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
          "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def _expand(c: str) -> list[str]:
    out = [""]
    for ch in c:
        out = [o + b for o in out for b in _IUPAC.get(ch, "")]
    return out


# ---------------------------------------------------------------- codon model
@dataclass
class Model:
    table: int
    states: list[str]
    pi: np.ndarray
    evals: np.ndarray
    left: np.ndarray            # P(t) = left @ diag(exp(evals t)) @ right
    right: np.ndarray
    rates: np.ndarray           # gamma category rates
    _cache: dict = field(default_factory=dict)

    def P(self, t: float) -> np.ndarray:
        key = round(float(t), 6)
        p = self._cache.get(key)
        if p is None:
            p = (self.left * np.exp(self.evals * t)) @ self.right
            np.clip(p, 0, None, out=p)
            p /= p.sum(1, keepdims=True)
            p = np.cumsum(p, 1)
            if len(self._cache) < 4096:
                self._cache[key] = p
        return p


def gy94(table: int, rng: np.random.Generator, kappa=2.5, omega=0.2, alpha=0.8) -> Model:
    c = code(table)
    states = c.sense
    n = len(states)
    f3x4 = rng.dirichlet([20.0] * 4, size=3)              # position-specific base freqs
    idx = {b: i for i, b in enumerate(BASES)}
    pi = np.array([np.prod([f3x4[k][idx[s[k]]] for k in range(3)]) for s in states])
    pi /= pi.sum()
    Q = np.zeros((n, n))
    for i, a in enumerate(states):
        for j, b in enumerate(states):
            diff = [k for k in range(3) if a[k] != b[k]]
            if len(diff) != 1:
                continue
            k = diff[0]
            r = pi[j]
            if frozenset((a[k], b[k])) in _TRANSITIONS:
                r *= kappa
            if c.aa[a] != c.aa[b]:
                r *= omega
            Q[i, j] = r
    np.fill_diagonal(Q, -Q.sum(1))
    Q /= -(pi * np.diag(Q)).sum()                         # 1 substitution per codon per unit
    d = np.sqrt(pi)
    S = (Q * d[:, None]) / d[None, :]
    S = (S + S.T) / 2
    evals, U = np.linalg.eigh(S)
    left = U / d[:, None]
    right = U.T * d[None, :]
    # 4 discrete gamma categories (median of each quartile), mean-normalised
    from scipy.stats import gamma as G
    q = G.ppf([0.125, 0.375, 0.625, 0.875], alpha, scale=1 / alpha)
    rates = q / q.mean()
    return Model(table, states, pi, evals, left, right, rates)


def evolve_subst(model: Model, seq: np.ndarray, cats: np.ndarray, t: float, rng, protect_first=True):
    """Return a new codon-state array after a branch of length t."""
    out = seq.copy()
    for k, r in enumerate(model.rates):
        m = cats == k
        if protect_first:
            m[0] = False
        if not m.any():
            continue
        cum = model.P(t * r)
        u = rng.random(int(m.sum()))
        new = (cum[seq[m]] < u[:, None]).sum(1)
        out[m] = np.minimum(new, len(model.states) - 1)
    return out


# ---------------------------------------------------------------- tree
def coalescent(n: int, height: float, rng) -> tuple[list[int], list[float], list[list[int]]]:
    """Return (parent, branch_length, children) for 2n-1 nodes; leaves are 0..n-1, root 2n-2."""
    parent = [-1] * (2 * n - 1)
    times = [0.0] * (2 * n - 1)
    children: list[list[int]] = [[] for _ in range(2 * n - 1)]
    active = list(range(n))
    t = 0.0
    nxt = n
    while len(active) > 1:
        k = len(active)
        t += rng.exponential(2.0 / (k * (k - 1)))
        i, j = rng.choice(k, 2, replace=False)
        a, b = active[i], active[j]
        for x in (a, b):
            parent[x] = nxt
            children[nxt].append(x)
        times[nxt] = t
        active = [x for x in active if x not in (a, b)] + [nxt]
        nxt += 1
    scale = height / t if t > 0 else 0.0
    bl = [0.0 if parent[x] < 0 else (times[parent[x]] - times[x]) * scale for x in range(2 * n - 1)]
    return parent, bl, children


# ---------------------------------------------------------------- family simulation
@dataclass
class Family:
    table: int
    ids: list[str]
    states: list[str]           # codon for each state index
    mat: np.ndarray             # n x ncol int16 state indices, -1 = gap

    @property
    def ncol(self) -> int:
        return self.mat.shape[1]

    def rows(self) -> list[list]:
        """Per leaf: list of codon strings (None = gap)."""
        st = self.states
        return [[st[x] if x >= 0 else None for x in r] for r in self.mat.tolist()]


def _geom(rng, mean):
    return min(30, int(rng.geometric(1.0 / mean)))


def _indels(cols, seq, cats, t, rate, mean, rng, order, model):
    """Apply codon indels on a branch. `order` is the global column order (list), mutated."""
    n_ev = rng.poisson(rate * t * max(len(seq), 1))
    for _ in range(n_ev):
        L = len(seq)
        if L <= 2:
            break
        ln = _geom(rng, mean)
        if rng.random() < 0.5:
            p = int(rng.integers(1, L + 1))           # insert after position p-1 (never before the Met)
            pred = int(cols[p - 1])
            new_ids = list(range(len(order), len(order) + ln))
            gi = order.index(pred) + 1
            order[gi:gi] = new_ids
            cols = np.insert(cols, p, new_ids)
            seq = np.insert(seq, p, rng.choice(len(model.states), ln, p=model.pi))
            cats = np.insert(cats, p, rng.integers(0, 4, ln))
        else:
            ln = min(ln, L - 2)
            p = int(rng.integers(1, L - ln + 1))      # never delete codon 0
            cols = np.delete(cols, slice(p, p + ln))
            seq = np.delete(seq, slice(p, p + ln))
            cats = np.delete(cats, slice(p, p + ln))
    return cols, seq, cats


def _tree_leaves(n, L, rng, table, height, kappa, omega, alpha, indel_rate, indel_mean):
    model = gy94(table, rng, kappa, omega, alpha)
    root_seq = rng.choice(len(model.states), L, p=model.pi)
    root_seq[0] = model.states.index("ATG")
    order = list(range(L))
    _, bl, children = coalescent(n, height, rng)
    root = 2 * n - 2
    state = {root: (np.arange(L), root_seq, rng.integers(0, 4, L))}
    stack, leaves = [root], {}
    while stack:
        v = stack.pop()
        cols, seq, cats = state.pop(v)
        if not children[v]:
            leaves[v] = (cols, seq)
            continue
        for ch in children[v]:
            s2 = evolve_subst(model, seq, cats, bl[ch], rng)
            state[ch] = _indels(cols, s2, cats, bl[ch], indel_rate, indel_mean, rng, order, model)
            stack.append(ch)
    return model, order, [leaves[i] for i in range(n)]


def _matrix(order, leaves):
    """(global column order, [(cols, seq)]) -> n x ncol matrix, all-gap columns dropped."""
    pos = np.full(len(order), -1, np.int64)
    pos[np.array(order)] = np.arange(len(order))
    mat = np.full((len(leaves), len(order)), -1, np.int16)
    for i, (cols, seq) in enumerate(leaves):
        mat[i, pos[cols]] = seq
    return mat[:, (mat >= 0).any(0)]


def simulate_family(n: int, L: int, rng: np.random.Generator, table: int = 1, height: float = 0.6,
                    kappa: float = 2.5, omega: float = 0.2, alpha: float = 0.8,
                    indel_rate: float = 0.03, indel_mean: float = 2.0, id_prefix: str = "seq") -> Family:
    """Full-tree simulation (for up to a few thousand sequences)."""
    n = max(n, 2)
    model, order, leaves = _tree_leaves(n, L, rng, table, height, kappa, omega, alpha, indel_rate, indel_mean)
    return Family(table, [f"{id_prefix}{i + 1}" for i in range(n)], model.states, _matrix(order, leaves))


def simulate_clades(n: int, L: int, rng: np.random.Generator, table: int = 1, n_clades: int = 256,
                    height: float = 0.5, leaf_t: float = 0.05, leaf_del: float = 0.01,
                    kappa: float = 2.5, omega: float = 0.2, alpha: float = 0.8,
                    indel_rate: float = 0.03, indel_mean: float = 2.0, id_prefix: str = "seq") -> Family:
    """Fast large-N variant: full tree for K clade ancestors (with indels), then each leaf is an
    independent short branch from its clade ancestor (substitutions + single-codon deletions)."""
    k = max(2, min(n_clades, n))
    model, order, anc = _tree_leaves(k, L, rng, table, height, kappa, omega, alpha, indel_rate, indel_mean)
    base = _matrix(order, anc)
    ncol = base.shape[1]
    sizes = np.full(k, n // k)
    sizes[: n % k] += 1
    nstate = len(model.states)
    # per category: cumulative P(leaf_t * rate) and probability of no change
    cums = [model.P(leaf_t * r) for r in model.rates]
    stays = [np.diag(np.diff(c, axis=1, prepend=0.0)) for c in cums]
    mat = np.empty((n, ncol), np.int16)
    row = 0
    for ci in range(k):
        m = int(sizes[ci])
        if m == 0:
            continue
        blk = np.repeat(base[ci][None, :], m, 0)
        cats = rng.integers(0, 4, ncol)
        for cat in range(4):
            sel = (blk >= 0) & (cats[None, :] == cat)
            sel[:, 0] = False
            vals = blk[sel].astype(np.int64)
            change = rng.random(vals.size) > stays[cat][vals]
            if change.any():
                c = vals[change]
                new = np.minimum((cums[cat][c] < rng.random(c.size)[:, None]).sum(1), nstate - 1)
                vals[change] = new
                blk[sel] = vals
        dele = (rng.random(blk.shape) < leaf_del) & (blk >= 0)
        dele[:, 0] = False
        blk[dele] = -1
        mat[row:row + m] = blk
        row += m
    mat = mat[:, (mat >= 0).any(0)]
    return Family(table, [f"{id_prefix}{i + 1}" for i in range(n)], model.states, mat)
