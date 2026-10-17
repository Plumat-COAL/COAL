"""NCBI genetic code tables (data from NCBI "The Genetic Codes",
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for the viewer's own
codon checks. Kept separate from the engine on purpose (the app never imports coal).

Each table: 64 amino-acid letters and 64 start flags, codons ordered TCAG x TCAG x TCAG.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

CODONS = ["".join(c) for c in product("TCAG", repeat=3)]

_RAW = {
    1:  ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "---M------**--*----M---------------M----------------------------"),
    2:  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
         "----------**--------------------MMMM----------**---M------------"),
    3:  ("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "----------**----------------------MM----------------------------"),
    4:  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "--MM------**-------M------------MMMM---------------M------------"),
    5:  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
         "---M------**--------------------MMMM---------------M------------"),
    6:  ("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "--------------*--------------------M----------------------------"),
    9:  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "----------**-----------------------M---------------M------------"),
    10: ("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "----------**-----------------------M----------------------------"),
    11: ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "---M------**--*----M------------MMMM---------------M------------"),
    12: ("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "----------**--*----M---------------M----------------------------"),
    13: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
         "---M------**----------------------MM---------------M------------"),
    14: ("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "-----------*-----------------------M----------------------------"),
    15: ("FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "----------*---*--------------------M----------------------------"),
    16: ("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "----------*---*--------------------M----------------------------"),
    21: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "----------**-----------------------M---------------M------------"),
    22: ("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "------*---*---*--------------------M----------------------------"),
    23: ("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "--*-------**--*-----------------M--M---------------M------------"),
}



@lru_cache(maxsize=None)
def table(number: int) -> tuple[dict[str, str], frozenset[str]]:
    """(codon -> amino acid letter, start codons) for an NCBI code number."""
    aas, starts = _RAW[int(number)]
    return dict(zip(CODONS, aas)), frozenset(c for c, s in zip(CODONS, starts) if s == "M")


IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


@lru_cache(maxsize=4096)
def expand(codon: str) -> tuple[str, ...]:
    """All plain codons an IUPAC triplet stands for (U read as T); () if it is not a triplet of IUPAC letters."""
    c = codon.upper().replace("U", "T")
    if len(c) != 3 or any(b not in IUPAC for b in c):
        return ()
    return tuple("".join(x) for x in product(*(IUPAC[b] for b in c)))


def translate(codon: str, number: int = 1) -> str:
    """Amino acid for a codon. An ambiguous codon (IUPAC, e.g. GCN) gives its amino acid only if every
    codon it stands for gives the same one (the engine's rule), else 'X'. Non-IUPAC characters give 'X'."""
    code = table(number)[0]
    aas = {code[c] for c in expand(codon)}
    return aas.pop() if len(aas) == 1 else "X"
