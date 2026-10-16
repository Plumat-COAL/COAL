"""Genetic code (NCBI table 1) with IUPAC-ambiguity-aware codon checking."""
from __future__ import annotations

from functools import lru_cache
from itertools import product

_BASES = "TCAG"
_AA1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
STANDARD = {a + b + c: _AA1[i] for i, (a, b, c) in enumerate(product(_BASES, repeat=3))}

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


@lru_cache(maxsize=None)
def translations(codon: str) -> frozenset[str]:
    """All amino acids a (possibly ambiguous) codon can encode. Empty set if invalid."""
    codon = codon.upper().replace("U", "T")
    try:
        options = [IUPAC[b] for b in codon]
    except KeyError:
        return frozenset()
    return frozenset(STANDARD["".join(p)] for p in product(*options))


def codon_matches(codon: str, aa: str) -> bool:
    """True if `codon` is consistent with residue `aa` ('X' accepts any valid codon)."""
    t = translations(codon)
    if not t:
        return False
    if aa == "X":
        return True
    return t == {aa}
