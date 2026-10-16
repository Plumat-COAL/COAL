"""Minimal readers for FASTA and CLUSTAL protein alignments and FASTA DNA."""
from __future__ import annotations

import re
from pathlib import Path

_NON_LETTER = re.compile(r"[^A-Za-z]")


class InputFormatError(ValueError):
    pass


def read_fasta(path: str | Path, letters_only: bool = False) -> list[tuple[str, str]]:
    """Return [(id, sequence)]. ID = first whitespace token of the header.
    Whitespace inside sequence lines is removed. With letters_only (used for DNA,
    as PAL2NAL does), every non-letter is removed (gaps, digits, '*', ...) and
    lines starting with '#' are skipped."""
    records: list[tuple[str, list[str]]] = []
    with open(path, encoding="ascii", errors="strict") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                parts = line[1:].split()
                if not parts:
                    raise InputFormatError(f"{path}:{lineno}: empty FASTA header")
                records.append((parts[0], []))
            elif line.strip():
                if letters_only and line.startswith("#"):
                    continue
                if not records:
                    raise InputFormatError(f"{path}:{lineno}: sequence data before first '>' header")
                if letters_only:
                    chunk = line if line.isalpha() else _NON_LETTER.sub("", line)
                else:
                    chunk = "".join(line.split())
                records[-1][1].append(chunk)
    _check_unique([r[0] for r in records], path)
    return [(rid, "".join(chunks)) for rid, chunks in records]


def read_clustal(path: str | Path) -> list[tuple[str, str]]:
    """Return [(id, aligned sequence)] from a CLUSTAL file (interleaved blocks)."""
    order: list[str] = []
    seqs: dict[str, list[str]] = {}
    with open(path, encoding="ascii", errors="strict") as fh:
        first = fh.readline()
        if not first.startswith("CLUSTAL"):
            raise InputFormatError(f"{path}: not a CLUSTAL file")
        for line in fh:
            line = line.rstrip("\r\n")
            # Skip blank lines and conservation lines (start with whitespace).
            if not line.strip() or line[0] in " \t":
                continue
            parts = line.split()
            if len(parts) < 2:
                raise InputFormatError(f"{path}: malformed CLUSTAL line: {line!r}")
            name, chunk = parts[0], parts[1]
            if name not in seqs:
                order.append(name)
                seqs[name] = []
            seqs[name].append(chunk)
    return [(n, "".join(seqs[n])) for n in order]


def read_alignment(path: str | Path) -> list[tuple[str, str]]:
    """Read a protein alignment, auto-detecting CLUSTAL vs FASTA."""
    with open(path, encoding="ascii", errors="strict") as fh:
        head = fh.read(7)
    return read_clustal(path) if head == "CLUSTAL" else read_fasta(path)


def _check_unique(ids: list[str], path) -> None:
    seen = set()
    for i in ids:
        if i in seen:
            raise InputFormatError(f"{path}: duplicate sequence ID {i!r}")
        seen.add(i)


def read_dna(paths) -> list[tuple[str, str]]:
    """Read one or more DNA FASTA files (records concatenated in order)."""
    out: list[tuple[str, str]] = []
    for p in paths:
        out += read_fasta(p, letters_only=True)
    _check_unique([r[0] for r in out], ", ".join(map(str, paths)))
    return out
