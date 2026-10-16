"""NCBI genetic codes supported by PAL2NAL v14 (transl_table numbers).

Each table: 64 amino-acid letters and 64 start flags in the standard NCBI
codon order (first base T,C,A,G; then second; then third). Source: NCBI
"The Genetic Codes" (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

BASES = "TCAG"
CODONS = ["".join(c) for c in product(BASES, repeat=3)]

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

SUPPORTED = tuple(sorted(_RAW))


class GeneticCode:
    def __init__(self, number: int):
        if number not in _RAW:
            raise ValueError(f"unsupported codon table {number}; supported: {SUPPORTED}")
        aas, starts = _RAW[number]
        self.number = number
        self.code = dict(zip(CODONS, aas))
        self.starts = frozenset(c for c, s in zip(CODONS, starts) if s == "M")
        # aa -> set of plain codons; '*' = stops
        self.codons: dict[str, frozenset[str]] = {}
        for c, a in self.code.items():
            self.codons[a] = self.codons.get(a, frozenset()) | {c}

    def __repr__(self):
        return f"GeneticCode({self.number})"


@lru_cache(maxsize=None)
def get_code(number: int) -> GeneticCode:
    return GeneticCode(number)
