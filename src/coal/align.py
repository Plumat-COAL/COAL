"""Pairwise protein alignment for the codon step (docs/alignment.md).

Global Needleman-Wunsch/Gotoh alignment with affine gaps, full length (overhangs are gapped, never cut off), and the
expansion of a local DIAMOND hit to a full-length pairwise alignment without realignment. The same scoring as the C++
library (cpp/src/align.hpp) and its command-line tool ``coal-cpp-align``:

* a gap of length k costs ``gap_open + (k - 1) * gap_extend`` (EMBOSS needle / Biopython convention);
* ``free_end_gaps``: gaps before the first or after the last residue of either sequence cost nothing (EMBOSS
  needle's default); off, they cost the same as internal gaps (pwalign ``type = "global"``);
* residues: letters (either case) and ``*``; letters missing from the matrix (J, O, U) score as X.

:func:`align_pairs` runs ``coal-cpp-align`` (multithreaded, ~1 GCUPS per thread) when it is found
(:func:`find_align_binary`), else this module's pure-Python kernel (~1,000x slower; fine for a few short pairs).
:func:`expand_local_hit`, :func:`score_alignment` and :func:`read_diamond_hits` are pure Python and cheap.

Not part of the stable ``import coal`` API yet (python-api.md): it may change before 1.0.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from typing import Iterable, Iterator, Sequence

__all__ = ["AlignParams", "PairAlignment", "LocalHit", "align_pair", "align_pairs", "score_alignment",
           "expand_local_hit", "read_diamond_hits", "find_align_binary", "MATRICES", "DIAMOND_COLUMNS"]

ALPHABET = "ARNDCQEGHILKMFPSTWYVBZX*"
# NCBI matrices (the values of Biopython's / NCBI's files; the same as cpp/src/align.cpp)
_MATRIX_TEXT = {
    "BLOSUM62": """
  4  -1  -2  -2   0  -1  -1   0  -2  -1  -1  -1  -1  -2  -1   1   0  -3  -2   0  -2  -1   0  -4
 -1   5   0  -2  -3   1   0  -2   0  -3  -2   2  -1  -3  -2  -1  -1  -3  -2  -3  -1   0  -1  -4
 -2   0   6   1  -3   0   0   0   1  -3  -3   0  -2  -3  -2   1   0  -4  -2  -3   3   0  -1  -4
 -2  -2   1   6  -3   0   2  -1  -1  -3  -4  -1  -3  -3  -1   0  -1  -4  -3  -3   4   1  -1  -4
  0  -3  -3  -3   9  -3  -4  -3  -3  -1  -1  -3  -1  -2  -3  -1  -1  -2  -2  -1  -3  -3  -2  -4
 -1   1   0   0  -3   5   2  -2   0  -3  -2   1   0  -3  -1   0  -1  -2  -1  -2   0   3  -1  -4
 -1   0   0   2  -4   2   5  -2   0  -3  -3   1  -2  -3  -1   0  -1  -3  -2  -2   1   4  -1  -4
  0  -2   0  -1  -3  -2  -2   6  -2  -4  -4  -2  -3  -3  -2   0  -2  -2  -3  -3  -1  -2  -1  -4
 -2   0   1  -1  -3   0   0  -2   8  -3  -3  -1  -2  -1  -2  -1  -2  -2   2  -3   0   0  -1  -4
 -1  -3  -3  -3  -1  -3  -3  -4  -3   4   2  -3   1   0  -3  -2  -1  -3  -1   3  -3  -3  -1  -4
 -1  -2  -3  -4  -1  -2  -3  -4  -3   2   4  -2   2   0  -3  -2  -1  -2  -1   1  -4  -3  -1  -4
 -1   2   0  -1  -3   1   1  -2  -1  -3  -2   5  -1  -3  -1   0  -1  -3  -2  -2   0   1  -1  -4
 -1  -1  -2  -3  -1   0  -2  -3  -2   1   2  -1   5   0  -2  -1  -1  -1  -1   1  -3  -1  -1  -4
 -2  -3  -3  -3  -2  -3  -3  -3  -1   0   0  -3   0   6  -4  -2  -2   1   3  -1  -3  -3  -1  -4
 -1  -2  -2  -1  -3  -1  -1  -2  -2  -3  -3  -1  -2  -4   7  -1  -1  -4  -3  -2  -2  -1  -2  -4
  1  -1   1   0  -1   0   0   0  -1  -2  -2   0  -1  -2  -1   4   1  -3  -2  -2   0   0   0  -4
  0  -1   0  -1  -1  -1  -1  -2  -2  -1  -1  -1  -1  -2  -1   1   5  -2  -2   0  -1  -1   0  -4
 -3  -3  -4  -4  -2  -2  -3  -2  -2  -3  -2  -3  -1   1  -4  -3  -2  11   2  -3  -4  -3  -2  -4
 -2  -2  -2  -3  -2  -1  -2  -3   2  -1  -1  -2  -1   3  -3  -2  -2   2   7  -1  -3  -2  -1  -4
  0  -3  -3  -3  -1  -2  -2  -3  -3   3   1  -2   1  -1  -2  -2   0  -3  -1   4  -3  -2  -1  -4
 -2  -1   3   4  -3   0   1  -1   0  -3  -4   0  -3  -3  -2   0  -1  -4  -3  -3   4   1  -1  -4
 -1   0   0   1  -3   3   4  -2   0  -3  -3   1  -1  -3  -1   0  -1  -3  -2  -2   1   4  -1  -4
  0  -1  -1  -1  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2   0   0  -2  -1  -1  -1  -1  -1  -4
 -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4   1
""",
    "BLOSUM45": """
  5  -2  -1  -2  -1  -1  -1   0  -2  -1  -1  -1  -1  -2  -1   1   0  -2  -2   0  -1  -1   0  -5
 -2   7   0  -1  -3   1   0  -2   0  -3  -2   3  -1  -2  -2  -1  -1  -2  -1  -2  -1   0  -1  -5
 -1   0   6   2  -2   0   0   0   1  -2  -3   0  -2  -2  -2   1   0  -4  -2  -3   4   0  -1  -5
 -2  -1   2   7  -3   0   2  -1   0  -4  -3   0  -3  -4  -1   0  -1  -4  -2  -3   5   1  -1  -5
 -1  -3  -2  -3  12  -3  -3  -3  -3  -3  -2  -3  -2  -2  -4  -1  -1  -5  -3  -1  -2  -3  -2  -5
 -1   1   0   0  -3   6   2  -2   1  -2  -2   1   0  -4  -1   0  -1  -2  -1  -3   0   4  -1  -5
 -1   0   0   2  -3   2   6  -2   0  -3  -2   1  -2  -3   0   0  -1  -3  -2  -3   1   4  -1  -5
  0  -2   0  -1  -3  -2  -2   7  -2  -4  -3  -2  -2  -3  -2   0  -2  -2  -3  -3  -1  -2  -1  -5
 -2   0   1   0  -3   1   0  -2  10  -3  -2  -1   0  -2  -2  -1  -2  -3   2  -3   0   0  -1  -5
 -1  -3  -2  -4  -3  -2  -3  -4  -3   5   2  -3   2   0  -2  -2  -1  -2   0   3  -3  -3  -1  -5
 -1  -2  -3  -3  -2  -2  -2  -3  -2   2   5  -3   2   1  -3  -3  -1  -2   0   1  -3  -2  -1  -5
 -1   3   0   0  -3   1   1  -2  -1  -3  -3   5  -1  -3  -1  -1  -1  -2  -1  -2   0   1  -1  -5
 -1  -1  -2  -3  -2   0  -2  -2   0   2   2  -1   6   0  -2  -2  -1  -2   0   1  -2  -1  -1  -5
 -2  -2  -2  -4  -2  -4  -3  -3  -2   0   1  -3   0   8  -3  -2  -1   1   3   0  -3  -3  -1  -5
 -1  -2  -2  -1  -4  -1   0  -2  -2  -2  -3  -1  -2  -3   9  -1  -1  -3  -3  -3  -2  -1  -1  -5
  1  -1   1   0  -1   0   0   0  -1  -2  -3  -1  -2  -2  -1   4   2  -4  -2  -1   0   0   0  -5
  0  -1   0  -1  -1  -1  -1  -2  -2  -1  -1  -1  -1  -1  -1   2   5  -3  -1   0   0  -1   0  -5
 -2  -2  -4  -4  -5  -2  -3  -2  -3  -2  -2  -2  -2   1  -3  -4  -3  15   3  -3  -4  -2  -2  -5
 -2  -1  -2  -2  -3  -1  -2  -3   2   0   0  -1   0   3  -3  -2  -1   3   8  -1  -2  -2  -1  -5
  0  -2  -3  -3  -1  -3  -3  -3  -3   3   1  -2   1   0  -3  -1   0  -3  -1   5  -3  -3  -1  -5
 -1  -1   4   5  -2   0   1  -1   0  -3  -3   0  -2  -3  -2   0   0  -4  -2  -3   4   2  -1  -5
 -1   0   0   1  -3   4   4  -2   0  -3  -2   1  -1  -3  -1   0  -1  -2  -2  -3   2   4  -1  -5
  0  -1  -1  -1  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -2  -1  -1  -1  -1  -1  -5
 -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5  -5   1
""",
    "BLOSUM80": """
  7  -3  -3  -3  -1  -2  -2   0  -3  -3  -3  -1  -2  -4  -1   2   0  -5  -4  -1  -3  -2  -1  -8
 -3   9  -1  -3  -6   1  -1  -4   0  -5  -4   3  -3  -5  -3  -2  -2  -5  -4  -4  -2   0  -2  -8
 -3  -1   9   2  -5   0  -1  -1   1  -6  -6   0  -4  -6  -4   1   0  -7  -4  -5   5  -1  -2  -8
 -3  -3   2  10  -7  -1   2  -3  -2  -7  -7  -2  -6  -6  -3  -1  -2  -8  -6  -6   6   1  -3  -8
 -1  -6  -5  -7  13  -5  -7  -6  -7  -2  -3  -6  -3  -4  -6  -2  -2  -5  -5  -2  -6  -7  -4  -8
 -2   1   0  -1  -5   9   3  -4   1  -5  -4   2  -1  -5  -3  -1  -1  -4  -3  -4  -1   5  -2  -8
 -2  -1  -1   2  -7   3   8  -4   0  -6  -6   1  -4  -6  -2  -1  -2  -6  -5  -4   1   6  -2  -8
  0  -4  -1  -3  -6  -4  -4   9  -4  -7  -7  -3  -5  -6  -5  -1  -3  -6  -6  -6  -2  -4  -3  -8
 -3   0   1  -2  -7   1   0  -4  12  -6  -5  -1  -4  -2  -4  -2  -3  -4   3  -5  -1   0  -2  -8
 -3  -5  -6  -7  -2  -5  -6  -7  -6   7   2  -5   2  -1  -5  -4  -2  -5  -3   4  -6  -6  -2  -8
 -3  -4  -6  -7  -3  -4  -6  -7  -5   2   6  -4   3   0  -5  -4  -3  -4  -2   1  -7  -5  -2  -8
 -1   3   0  -2  -6   2   1  -3  -1  -5  -4   8  -3  -5  -2  -1  -1  -6  -4  -4  -1   1  -2  -8
 -2  -3  -4  -6  -3  -1  -4  -5  -4   2   3  -3   9   0  -4  -3  -1  -3  -3   1  -5  -3  -2  -8
 -4  -5  -6  -6  -4  -5  -6  -6  -2  -1   0  -5   0  10  -6  -4  -4   0   4  -2  -6  -6  -3  -8
 -1  -3  -4  -3  -6  -3  -2  -5  -4  -5  -5  -2  -4  -6  12  -2  -3  -7  -6  -4  -4  -2  -3  -8
  2  -2   1  -1  -2  -1  -1  -1  -2  -4  -4  -1  -3  -4  -2   7   2  -6  -3  -3   0  -1  -1  -8
  0  -2   0  -2  -2  -1  -2  -3  -3  -2  -3  -1  -1  -4  -3   2   8  -5  -3   0  -1  -2  -1  -8
 -5  -5  -7  -8  -5  -4  -6  -6  -4  -5  -4  -6  -3   0  -7  -6  -5  16   3  -5  -8  -5  -5  -8
 -4  -4  -4  -6  -5  -3  -5  -6   3  -3  -2  -4  -3   4  -6  -3  -3   3  11  -3  -5  -4  -3  -8
 -1  -4  -5  -6  -2  -4  -4  -6  -5   4   1  -4   1  -2  -4  -3   0  -5  -3   7  -6  -4  -2  -8
 -3  -2   5   6  -6  -1   1  -2  -1  -6  -7  -1  -5  -6  -4   0  -1  -8  -5  -6   6   0  -3  -8
 -2   0  -1   1  -7   5   6  -4   0  -6  -5   1  -3  -6  -2  -1  -2  -5  -4  -4   0   6  -1  -8
 -1  -2  -2  -3  -4  -2  -2  -3  -2  -2  -2  -2  -2  -3  -3  -1  -1  -5  -3  -2  -3  -1  -2  -8
 -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8   1
""",
    "PAM250": """
  2  -2   0   0  -2   0   0   1  -1  -1  -2  -1  -1  -3   1   1   1  -6  -3   0   0   0   0  -8
 -2   6   0  -1  -4   1  -1  -3   2  -2  -3   3   0  -4   0   0  -1   2  -4  -2  -1   0  -1  -8
  0   0   2   2  -4   1   1   0   2  -2  -3   1  -2  -3   0   1   0  -4  -2  -2   2   1   0  -8
  0  -1   2   4  -5   2   3   1   1  -2  -4   0  -3  -6  -1   0   0  -7  -4  -2   3   3  -1  -8
 -2  -4  -4  -5  12  -5  -5  -3  -3  -2  -6  -5  -5  -4  -3   0  -2  -8   0  -2  -4  -5  -3  -8
  0   1   1   2  -5   4   2  -1   3  -2  -2   1  -1  -5   0  -1  -1  -5  -4  -2   1   3  -1  -8
  0  -1   1   3  -5   2   4   0   1  -2  -3   0  -2  -5  -1   0   0  -7  -4  -2   3   3  -1  -8
  1  -3   0   1  -3  -1   0   5  -2  -3  -4  -2  -3  -5   0   1   0  -7  -5  -1   0   0  -1  -8
 -1   2   2   1  -3   3   1  -2   6  -2  -2   0  -2  -2   0  -1  -1  -3   0  -2   1   2  -1  -8
 -1  -2  -2  -2  -2  -2  -2  -3  -2   5   2  -2   2   1  -2  -1   0  -5  -1   4  -2  -2  -1  -8
 -2  -3  -3  -4  -6  -2  -3  -4  -2   2   6  -3   4   2  -3  -3  -2  -2  -1   2  -3  -3  -1  -8
 -1   3   1   0  -5   1   0  -2   0  -2  -3   5   0  -5  -1   0   0  -3  -4  -2   1   0  -1  -8
 -1   0  -2  -3  -5  -1  -2  -3  -2   2   4   0   6   0  -2  -2  -1  -4  -2   2  -2  -2  -1  -8
 -3  -4  -3  -6  -4  -5  -5  -5  -2   1   2  -5   0   9  -5  -3  -3   0   7  -1  -4  -5  -2  -8
  1   0   0  -1  -3   0  -1   0   0  -2  -3  -1  -2  -5   6   1   0  -6  -5  -1  -1   0  -1  -8
  1   0   1   0   0  -1   0   1  -1  -1  -3   0  -2  -3   1   2   1  -2  -3  -1   0   0   0  -8
  1  -1   0   0  -2  -1   0   0  -1   0  -2   0  -1  -3   0   1   3  -5  -3   0   0  -1   0  -8
 -6   2  -4  -7  -8  -5  -7  -7  -3  -5  -2  -3  -4   0  -6  -2  -5  17   0  -6  -5  -6  -4  -8
 -3  -4  -2  -4   0  -4  -4  -5   0  -1  -1  -4  -2   7  -5  -3  -3   0  10  -2  -3  -4  -2  -8
  0  -2  -2  -2  -2  -2  -2  -1  -2   4   2  -2   2  -1  -1  -1   0  -6  -2   4  -2  -2  -1  -8
  0  -1   2   3  -4   1   3   0   1  -2  -3   1  -2  -4  -1   0   0  -5  -3  -2   3   2  -1  -8
  0   0   1   3  -5   3   3   0   2  -2  -3   0  -2  -5   0   0  -1  -6  -4  -2   2   3  -1  -8
  0  -1   0  -1  -3  -1  -1  -1  -1  -1  -1  -1  -1  -2  -1   0   0  -4  -2  -1  -1  -1  -1  -8
 -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8   1
""",
}
MATRICES: dict[str, list[list[int]]] = {
    name: [[int(v) for v in line.split()] for line in text.strip().splitlines()] for name, text in _MATRIX_TEXT.items()
}
_CODE = {}
for _c in range(ord("A"), ord("Z") + 1):
    _ch = chr(_c)
    _CODE[_ch] = _CODE[_ch.lower()] = ALPHABET.index(_ch) if _ch in ALPHABET else ALPHABET.index("X")
_CODE["*"] = 23

# DIAMOND --outfmt 6 columns read by read_diamond_hits / coal-cpp-align --diamond (in this order, first)
DIAMOND_COLUMNS = ("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qseq_gapped", "sseq_gapped")


class AlignError(ValueError):
    """Bad residues, parameters or local hit (code: align_input, align_params, local_hit, align_too_large)."""

    def __init__(self, msg: str, code: str = "align_input"):
        super().__init__(msg)
        self.code = code


@dataclass(frozen=True)
class AlignParams:
    matrix: str = "BLOSUM62"
    gap_open: float = 10.0
    gap_extend: float = 0.5
    free_end_gaps: bool = True

    @classmethod
    def needle(cls) -> "AlignParams":
        """EMBOSS needle's defaults (COAL's default): BLOSUM62, open 10, extend 0.5, end gaps free."""
        return cls()

    @classmethod
    def pwalign(cls) -> "AlignParams":
        """pwalign::pairwiseAlignment(type="global", substitutionMatrix="BLOSUM62") with its default gapOpening=10,
        gapExtension=4: pwalign charges gapOpening + k*gapExtension, i.e. open 14, extend 4; end gaps penalised."""
        return cls(gap_open=14.0, gap_extend=4.0, free_end_gaps=False)

    def check(self) -> None:
        if self.matrix.upper() not in MATRICES:
            raise AlignError(f"unknown substitution matrix {self.matrix!r} (known: {', '.join(MATRICES)})",
                             "align_params")
        o, e = float(self.gap_open), float(self.gap_extend)
        if not (0 <= e <= o <= 1000):
            raise AlignError(f"gap costs must satisfy 1000 >= gap_open >= gap_extend >= 0 (given open {o}, "
                             f"extend {e})", "align_params")

    def args(self) -> list[str]:
        """The coal-cpp-align options for these parameters."""
        return ["--matrix", self.matrix, "--gap-open", repr(float(self.gap_open)), "--gap-extend",
                repr(float(self.gap_extend)), "--end-gaps", "free" if self.free_end_gaps else "penalized"]


@dataclass
class PairAlignment:
    a: str = ""
    b: str = ""
    score: float = 0.0
    length: int = 0
    identities: int = 0
    similarities: int = 0
    gaps: int = 0
    id_a: str = ""
    id_b: str = ""
    error: str | None = None


@dataclass
class LocalHit:
    """One DIAMOND hit: 1-based inclusive coordinates and the aligned parts with '-' for gaps."""
    qstart: int
    qend: int
    sstart: int
    send: int
    qseq_gapped: str
    sseq_gapped: str
    qseqid: str = ""
    sseqid: str = ""
    extra: list[str] = field(default_factory=list)


def _scale(p: AlignParams) -> tuple[int, int, int]:
    for k in (1, 2, 4, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000):
        o, e = p.gap_open * k, p.gap_extend * k
        if abs(o - round(o)) < 1e-6 and abs(e - round(e)) < 1e-6:
            return k, int(round(o)), int(round(e))
    raise AlignError(f"gap costs must be multiples of 0.001 (given {p.gap_open}, {p.gap_extend})", "align_params")


def _encode(seq: str, which: str) -> list[int]:
    out = []
    for i, ch in enumerate(seq):
        c = _CODE.get(ch)
        if c is None:
            raise AlignError(f"{which} sequence: {ch!r} at position {i + 1} is not an amino acid letter or '*' "
                             "(give ungapped protein sequences)")
        out.append(c)
    return out


def _stats(res: PairAlignment, mat) -> PairAlignment:
    res.length = len(res.a)
    res.identities = res.similarities = res.gaps = 0
    for x, y in zip(res.a, res.b):
        if x == "-" and y == "-":
            continue
        if x == "-" or y == "-":
            res.gaps += 1
            continue
        res.identities += x == y
        res.similarities += mat[_CODE[x]][_CODE[y]] > 0
    return res


def score_alignment(row_a: str, row_b: str, params: AlignParams = AlignParams()) -> float:
    """Score of two gapped rows under params (columns that are gaps in both rows are ignored)."""
    params.check()
    if len(row_a) != len(row_b):
        raise AlignError(f"alignment rows differ in length ({len(row_a)} and {len(row_b)})")
    k, o, e = _scale(params)
    mat = MATRICES[params.matrix.upper()]
    for which, row in (("first", row_a), ("second", row_b)):
        for i, ch in enumerate(row):
            if ch != "-" and ch not in _CODE:
                raise AlignError(f"{which} row: {ch!r} at column {i + 1} is not an amino acid letter, '*' or '-'")
    ia = [i for i, c in enumerate(row_a) if c != "-"]
    ib = [i for i, c in enumerate(row_b) if c != "-"]
    total, prev = 0, 0
    for c, (x, y) in enumerate(zip(row_a, row_b)):
        if x == "-" and y == "-":
            continue
        if x != "-" and y != "-":
            total += mat[_CODE[x]][_CODE[y]] * k
            prev = 0
            continue
        st = 1 if y == "-" else 2
        idx = ib if st == 1 else ia
        end = not idx or c < idx[0] or c > idx[-1]
        if not (end and params.free_end_gaps):
            total -= e if prev == st else o
        prev = st
    return total / k


def align_pair(a: str, b: str, params: AlignParams = AlignParams()) -> PairAlignment:
    """Global alignment of a against b in pure Python (O(len(a) * len(b)) time; use align_pairs for many or long
    pairs). Same kernel, tie-breaking and result as coal-cpp-align."""
    params.check()
    k, o, e = _scale(params)
    mat = MATRICES[params.matrix.upper()]
    ea, eb = _encode(a, "first"), _encode(b, "second")
    n, m = len(ea), len(eb)
    free = params.free_end_gaps
    gap = (lambda length: o + (length - 1) * e if length else 0)
    if n == 0 or m == 0:
        res = PairAlignment(a.upper() + "-" * m, "-" * n + b.upper(), score=(0 if free else -gap(n + m)) / k)
        return _stats(res, mat)
    NEG = -(1 << 40)
    hp = [0] + [0 if free else -gap(j) for j in range(1, m + 1)]
    xp = [NEG] * (m + 1)
    tb = []
    lastcol = [hp[m]]
    for i in range(1, n + 1):
        row = [s * k for s in mat[ea[i - 1]]]
        h0 = 0 if free else -gap(i)
        hc = [h0] + [0] * m
        xc = [NEG] * (m + 1)
        t = bytearray(m + 1)
        y, dprev = NEG, h0
        for j in range(1, m + 1):
            xo, xe = hp[j] - o, xp[j] - e
            x = xe if xe > xo else xo
            mm = hp[j - 1] + row[eb[j - 1]]
            d = x if x > mm else mm
            yo, ye = dprev - o, y - e
            y = ye if ye > yo else yo
            h = y if y > d else d
            xc[j], hc[j] = x, h
            t[j] = (xe > xo) | ((x > mm) << 1) | ((ye > yo) << 2) | ((y > d) << 3)
            dprev = d
        tb.append(t)
        lastcol.append(hc[m])
        hp, xp = hc, xc
    bi, bj, best = n, m, hp[m]
    if free:
        for j in range(m):
            if hp[j] > best:
                best, bi, bj = hp[j], n, j
        for i in range(n):
            if lastcol[i] > best:
                best, bi, bj = lastcol[i], i, m
    ra, rb = [], []
    for j in range(m, bj, -1):
        ra.append("-"); rb.append(b[j - 1])  # noqa: E702
    for i in range(n, bi, -1):
        ra.append(a[i - 1]); rb.append("-")  # noqa: E702
    i, j, state = bi, bj, 0
    while i > 0 and j > 0:
        c = tb[i - 1][j]
        if state == 0:
            if c & 8:
                state = 2
            elif c & 2:
                state = 1
            else:
                ra.append(a[i - 1]); rb.append(b[j - 1]); i -= 1; j -= 1  # noqa: E702
                continue
        if state == 1:
            ra.append(a[i - 1]); rb.append("-"); i -= 1  # noqa: E702
            state = 1 if c & 1 else 0
        else:
            ra.append("-"); rb.append(b[j - 1]); j -= 1  # noqa: E702
            state = 2 if c & 4 else 0
    while i > 0:
        ra.append(a[i - 1]); rb.append("-"); i -= 1  # noqa: E702
    while j > 0:
        ra.append("-"); rb.append(b[j - 1]); j -= 1  # noqa: E702
    res = PairAlignment("".join(reversed(ra)).upper(), "".join(reversed(rb)).upper(), score=best / k)
    return _stats(res, mat)


def expand_local_hit(query: str, subject: str, hit: LocalHit, flanks: str = "gap",
                     params: AlignParams = AlignParams()) -> PairAlignment:
    """The full-length pairwise alignment of query and subject that contains the hit's aligned block unchanged
    (C++: coal::expand_local_hit). Only the gap pattern of the hit is used; residues come from the sequences, and
    the hit must agree with them ('X' in the hit, e.g. DIAMOND's masking, matches anything).
    flanks="gap": the unaligned ends go against gaps, staggered (no realignment at all); "align": they are aligned
    with align_pair (pure Python; short)."""
    params.check()
    if flanks not in ("gap", "align"):
        raise AlignError(f"flanks: 'gap' or 'align', not {flanks!r}", "align_params")

    def fail(msg):
        raise AlignError("local hit: " + msg, "local_hit")
    nq, ns = len(query), len(subject)
    if not (1 <= hit.qstart <= hit.qend <= nq):
        fail(f"query coordinates {hit.qstart}..{hit.qend} do not fit a query of {nq} residues")
    if not (1 <= hit.sstart <= hit.send <= ns):
        fail(f"subject coordinates {hit.sstart}..{hit.send} do not fit a subject of {ns} residues")
    qg, sg = hit.qseq_gapped, hit.sseq_gapped
    if len(qg) != len(sg):
        fail(f"qseq_gapped and sseq_gapped differ in length ({len(qg)} and {len(sg)})")
    rq, rs = len(qg) - qg.count("-"), len(sg) - sg.count("-")
    if rq != hit.qend - hit.qstart + 1:
        fail(f"qseq_gapped has {rq} residues, qstart..qend spans {hit.qend - hit.qstart + 1}")
    if rs != hit.send - hit.sstart + 1:
        fail(f"sseq_gapped has {rs} residues, sstart..send spans {hit.send - hit.sstart + 1}")
    ca, cb = [], []
    qi, si = hit.qstart - 1, hit.sstart - 1
    for x, y in zip(qg, sg):
        if x == "-" and y == "-":
            continue
        if x != "-":
            if x.upper() not in ("X", query[qi].upper()):
                fail(f"qseq_gapped has {x!r} where the query has {query[qi]!r} (position {qi + 1})")
            ca.append(query[qi]); qi += 1  # noqa: E702
        else:
            ca.append("-")
        if y != "-":
            if y.upper() not in ("X", subject[si].upper()):
                fail(f"sseq_gapped has {y!r} where the subject has {subject[si]!r} (position {si + 1})")
            cb.append(subject[si]); si += 1  # noqa: E702
        else:
            cb.append("-")
    ql, sl = query[:hit.qstart - 1], subject[:hit.sstart - 1]
    qr, sr = query[hit.qend:], subject[hit.send:]
    if flanks == "gap":
        ra = ql + "-" * len(sl) + "".join(ca) + qr + "-" * len(sr)
        rb = "-" * len(ql) + sl + "".join(cb) + "-" * len(qr) + sr
    else:
        left, right = align_pair(ql, sl, params), align_pair(qr, sr, params)
        ra = left.a + "".join(ca) + right.a
        rb = left.b + "".join(cb) + right.b
    res = PairAlignment(ra.upper(), rb.upper(), id_a=hit.qseqid, id_b=hit.sseqid)
    res.score = score_alignment(res.a, res.b, params)
    return _stats(res, MATRICES[params.matrix.upper()])


def read_diamond_hits(path: str, first_only: bool = True) -> Iterator[LocalHit]:
    """Hits from DIAMOND --outfmt 6 whose first columns are DIAMOND_COLUMNS (later columns go to .extra).
    first_only: only the first line of each query-subject pair (DIAMOND writes the best HSP first)."""
    seen = set()
    with open(path, encoding="utf-8", errors="backslashreplace") as fh:
        for ln, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 8:
                raise AlignError(f"{path} line {ln}: expected the columns {' '.join(DIAMOND_COLUMNS)}", "local_hit")
            if first_only:
                if (f[0], f[1]) in seen:
                    continue
                seen.add((f[0], f[1]))
            try:
                q0, q1, s0, s1 = (int(v) for v in f[2:6])
            except ValueError:
                raise AlignError(f"{path} line {ln}: qstart qend sstart send must be whole numbers",
                                 "local_hit") from None
            yield LocalHit(q0, q1, s0, s1, f[6], f[7], qseqid=f[0], sseqid=f[1], extra=f[8:])


# ---------------------------------------------------------------- coal-cpp-align

_NAMES = ("coal-cpp-align.exe", "coal-cpp-align") if os.name == "nt" else ("coal-cpp-align",)


def find_align_binary() -> str | None:
    """coal-cpp-align: $COAL_CPP_ALIGN; inside the package (coal/bin); next to the Python; next to the coal-cpp
    that coal would use (coal --cpp-path); PATH."""
    env = os.environ.get("COAL_CPP_ALIGN")
    if env:
        return env if os.path.isfile(env) else None
    from . import dispatch
    dirs = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin"), *dispatch._dirs()]
    cpp = dispatch.find_cpp()
    if cpp:
        dirs.append(os.path.dirname(os.path.abspath(cpp)))
    for d in dirs:
        for n in _NAMES:
            p = os.path.join(d, n)
            if os.path.isfile(p) and os.access(p, os.X_OK):
                return p
    for n in _NAMES:
        p = shutil.which(n)
        if p:
            return p
    return None


def _write_fasta(path: str, recs: Iterable[tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n{seq}\n")


def align_pairs(pairs: Sequence[tuple[str, str]] | Sequence[tuple[str, str, str, str]],
                params: AlignParams = AlignParams(), jobs: int = 0, binary: str | None = None,
                python_fallback: bool = True, timeout: float | None = None) -> list[PairAlignment]:
    """Align many pairs: (seq_a, seq_b) or (id_a, seq_a, id_b, seq_b) tuples. Results in input order; a pair that
    cannot be aligned has .error set (the others are unaffected). Runs coal-cpp-align (binary, or
    find_align_binary()) on `jobs` threads (0 = all); without it, the pure-Python align_pair if python_fallback,
    else RuntimeError."""
    params.check()
    norm = []
    for i, p in enumerate(pairs):
        if len(p) == 2:
            norm.append((f"a{i}", p[0], f"b{i}", p[1]))
        elif len(p) == 4:
            norm.append(tuple(p))
        else:
            raise ValueError("pairs: (seq_a, seq_b) or (id_a, seq_a, id_b, seq_b) tuples")
    exe = binary or find_align_binary()
    if exe is None:
        if not python_fallback:
            raise RuntimeError("coal-cpp-align was not found (build it: cmake -S . -B build && cmake --build build; "
                               "or set COAL_CPP_ALIGN)")
        out = []
        for ia, sa, ib, sb in norm:
            try:
                r = align_pair(sa, sb, params)
            except AlignError as e:
                r = PairAlignment(error=str(e))
            r.id_a, r.id_b = ia, ib
            out.append(r)
        return out
    with tempfile.TemporaryDirectory(prefix="coal_align_") as tmp:
        inp, outp, sc = (os.path.join(tmp, n) for n in ("in.fasta", "out.fasta", "scores.tsv"))
        # record names are positions, so IDs with spaces or duplicates are no problem
        _write_fasta(inp, ((f"{2 * i}" if k == 0 else f"{2 * i + 1}", s)
                           for i, (_, sa, _, sb) in enumerate(norm) for k, s in ((0, sa), (1, sb))))
        cmd = [exe, inp, "-o", outp, "--scores", sc, "-j", str(jobs), *params.args()]
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="backslashreplace", timeout=timeout)
        if r.returncode not in (0, 2):
            raise RuntimeError(f"coal-cpp-align failed (exit {r.returncode}): {r.stderr.strip()[-1000:]}")
        rows = {}
        with open(outp, encoding="utf-8", errors="backslashreplace") as fh:
            name = None
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    name = int(line[1:])
                else:
                    rows[name] = line
        res = []
        with open(sc, encoding="utf-8", errors="backslashreplace") as fh:
            next(fh)
            for i, line in enumerate(fh):
                f = line.rstrip("\n").split("\t")
                ia, _, ib, _ = norm[i]
                if f[10] != "ok":
                    res.append(PairAlignment(id_a=ia, id_b=ib, error=f[10]))
                    continue
                res.append(PairAlignment(rows[2 * i], rows[2 * i + 1], float(f[6]), int(f[5]), int(f[7]), int(f[8]),
                                         int(f[9]), ia, ib))
        return res
