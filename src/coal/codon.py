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
ANCHOR_STARTS = frozenset({"ATG", "CTG", "GTG"})  # PAL2NAL's first-Met rule in fallback mode

# Residue "modes": how the first Met of a sequence is treated.
NORMAL, START, ANCHOR_START = 0, 1, 2


def _allowed(table: int, aa: str, mode: int) -> frozenset[str]:
    code = get_code(table)
    if aa == "M" and mode == START:
        return code.starts
    if aa == "M" and mode == ANCHOR_START:
        return ANCHOR_STARTS
    return code.codons.get(aa, frozenset())


@lru_cache(maxsize=None)
def accepted(table: int, aa: str, mode: int = NORMAL) -> frozenset[str]:
    """All upper-case IUPAC triplets that encode `aa` (no 'X' handling here)."""
    allowed = _allowed(table, aa, mode)
    out = set()
    for trip in product(IUPAC, repeat=3):
        if all("".join(c) in allowed for c in product(*(IUPAC[b] for b in trip))):
            out.add("".join(trip))
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


