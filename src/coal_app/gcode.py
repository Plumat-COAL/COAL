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


def translate(codon: str, number: int = 1) -> str:
    """Amino acid for a codon ('X' if it contains anything but A/C/G/T/U)."""
    return table(number)[0].get(codon.upper().replace("U", "T"), "X")
