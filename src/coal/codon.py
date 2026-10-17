"""Codon acceptance rules, built on the NCBI tables in tables.py.

A DNA triplet "encodes" a residue if every concrete codon it can stand for
(IUPAC ambiguity codes expanded; U treated as T) is in the residue's codon set.
DNA is normalised to upper case with U->T before any of these are used.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from .tables import get_code

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
ANCHOR_STARTS = frozenset({"ATG", "CTG", "GTG"})  # PAL2NAL's first-Met rule in fallback mode (FIN-E1: + the code's)

# Residue "modes": how the first Met of a sequence is treated.
NORMAL, START, ANCHOR_START = 0, 1, 2


def _allowed(table: int, aa: str, mode: int) -> frozenset[str]:
    code = get_code(table)
    if aa == "M" and mode == START:
        return code.starts
    if aa == "M" and mode == ANCHOR_START:
        # CHANGED(FIN-E1): was ATG/CTG/GTG only (PAL2NAL's fallback rule) whatever the table, so code 2's ATT/ATA start
        # was a mismatch (and -nomismatch dropped column 1) whenever a gene needed the fallback; now also the selected
        # code's start and Met codons, as in the exact-fit path. Deliberate difference from PAL2NAL v14.
        # teammate's engine/C++ must mirror. See overnight/integration/CHANGES.md
        return ANCHOR_STARTS | code.starts | code.codons.get("M", frozenset())
    return code.codons.get(aa, frozenset())


_BIT = {"A": 1, "C": 2, "G": 4, "T": 8}
_MASK = {k: sum(_BIT[b] for b in v) for k, v in IUPAC.items()}  # IUPAC code -> bit set of its bases


@lru_cache(maxsize=None)
def accepted(table: int, aa: str, mode: int = NORMAL) -> frozenset[str]:
    """All upper-case IUPAC triplets that encode `aa` (no 'X' handling here).
    CHANGED(IMP-E-H): was a brute-force expansion of all 15^3 triplets into their concrete codons (3.2 ms per
    residue/table/mode, ~1 s for the 17 codes the genetic-code hint tries); now from per-position base sets: for
    each pair of first two IUPAC bases, the third bases every concrete prefix allows. Same sets.
    See overnight/integration/CHANGES.md"""
    allowed = _allowed(table, aa, mode)
    third: dict[str, int] = {}  # concrete two-base prefix -> bit set of third bases that complete it
    for c in allowed:
        if len(c) == 3 and c[2] in _BIT and c[0] in _BIT and c[1] in _BIT:
            third[c[:2]] = third.get(c[:2], 0) | _BIT[c[2]]
    out = []
    for b1, s1 in IUPAC.items():
        if not all(any(p[0] == x for p in third) for x in s1):
            continue
        for b2, s2 in IUPAC.items():
            inter = 15
            for x in s1:
                for y in s2:
                    inter &= third.get(x + y, 0)
                    if not inter:
                        break
                if not inter:
                    break
            if inter:
                for b3, m in _MASK.items():
                    if not m & ~inter:
                        out.append(b1 + b2 + b3)
    return frozenset(out)


@lru_cache(maxsize=None)
def fast_code(table: int) -> dict[str, str]:
    """Plain upper-case codon -> amino acid ('*' for stop)."""
    return dict(get_code(table).code)


@lru_cache(maxsize=None)
def residue_regex(table: int, aa: str, mode: int = NORMAL) -> str:
    """Regex (upper-case IUPAC alphabet) matching exactly `accepted(table, aa, mode)`."""
    if aa == "X":
        return "[A-Z]{3}"
    trie: dict = {}
    for t in sorted(accepted(table, aa, mode)):
        node = trie
        for ch in t:
            node = node.setdefault(ch, {})

    def emit(node) -> str:
        # Group children that share an identical subtree into one character class.
        groups: dict[str, list[str]] = {}
        for ch, sub in node.items():
            groups.setdefault(emit(sub) if sub else "", []).append(ch)
        alts = []
        for tail, chars in groups.items():
            cls = chars[0] if len(chars) == 1 else "[" + "".join(sorted(chars)) + "]"
            alts.append(cls + tail)
        return alts[0] if len(alts) == 1 else "(?:" + "|".join(alts) + ")"

    return emit(trie) if trie else "(?!)"


def codon_ok(table: int, codon: str, aa: str, mode: int = NORMAL) -> bool:
    if aa == "X":
        return len(codon) == 3
    return codon in accepted(table, aa, mode)


