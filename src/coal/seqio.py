"""Minimal readers for FASTA and CLUSTAL protein alignments and FASTA DNA."""
from __future__ import annotations

import re
from array import array
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


# ---------------------------------------------------------------- streaming

def iter_fasta(path: str | Path):
    """Yield (id, sequence) one record at a time (whitespace removed)."""
    rid, chunks = None, []
    seen = set()
    with open(path, encoding="ascii", errors="strict") as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(chunks)
                parts = line[1:].split()
                if not parts:
                    raise InputFormatError(f"{path}:{lineno}: empty FASTA header")
                rid, chunks = parts[0], []
                if rid in seen:
                    raise InputFormatError(f"{path}: duplicate sequence ID {rid!r}")
                seen.add(rid)
            elif line.strip():
                if rid is None:
                    raise InputFormatError(f"{path}:{lineno}: sequence data before first '>' header")
                chunks.append("".join(line.split()))
    if rid is not None:
        yield rid, "".join(chunks)


def is_clustal(path: str | Path) -> bool:
    with open(path, encoding="ascii", errors="strict") as fh:
        return fh.read(7) == "CLUSTAL"


_LETTERS = bytes(range(65, 91)) + bytes(range(97, 123))
_NON_LETTER_BYTES = bytes(b for b in range(256) if b not in _LETTERS)


class DnaIndex:
    """Byte-offset index of one or more DNA FASTA files.

    Only IDs and offsets are kept in memory; each sequence is read from disk
    when needed. Sequences in the order they will be requested (the usual case)
    are read sequentially.
    """

    def __init__(self, paths):
        self.paths = [str(p) for p in paths]
        self.ids: list[str] = []
        # compact per-record location: file number, byte start, byte end
        self.fno, self.start, self.end = array("H"), array("q"), array("q")
        self.pos: dict[str, int] = {}
        for fno, path in enumerate(self.paths):
            self._scan(fno, path)
        self._fh: dict[int, object] = {}

    def _scan(self, fno: int, path: str) -> None:
        with open(path, "rb") as fh:
            offset, start = 0, None
            for line in fh:
                if line.startswith(b">"):
                    if start is not None:
                        self._add(fno, start, offset)
                    parts = line[1:].split()
                    if not parts:
                        raise InputFormatError(f"{path}: empty FASTA header")
                    rid = parts[0].decode("ascii")
                    if rid in self.pos:
                        raise InputFormatError(f"{', '.join(self.paths)}: duplicate sequence ID {rid!r}")
                    self.pos[rid] = len(self.ids)
                    self.ids.append(rid)
                    start = offset + len(line)
                elif start is None and line.strip():
                    raise InputFormatError(f"{path}: sequence data before first '>' header")
                offset += len(line)
            if start is not None:
                self._add(fno, start, offset)

    def _add(self, fno: int, a: int, b: int) -> None:
        self.fno.append(fno)
        self.start.append(a)
        self.end.append(b)

    def forget_ids(self) -> None:
        """Free the ID tables once sequences have been paired (records stay readable)."""
        self.ids, self.pos = [], {}

    def __len__(self) -> int:
        return len(self.start)

    def total_bytes(self) -> int:
        """Size of all sequence data (including line breaks): a cheap work estimate."""
        return sum(self.end) - sum(self.start)

    def get(self, i: int) -> str:
        """Sequence number i: non-letters removed, '#' lines skipped (as PAL2NAL)."""
        fno, a, b = self.fno[i], self.start[i], self.end[i]
        fh = self._fh.get(fno)
        if fh is None:
            fh = self._fh[fno] = open(self.paths[fno], "rb")
        fh.seek(a)
        raw = fh.read(b - a)
        if b"#" in raw:
            raw = b"".join(l for l in raw.splitlines(True) if not l.startswith(b"#"))
        return raw.translate(None, _NON_LETTER_BYTES).decode("ascii")

    def close(self) -> None:
        for fh in self._fh.values():
            fh.close()
        self._fh.clear()
