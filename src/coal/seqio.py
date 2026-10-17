"""Minimal readers for FASTA and CLUSTAL protein alignments and FASTA DNA."""
from __future__ import annotations

import codecs
import errno
import os
import re
from array import array
from contextlib import contextmanager
from pathlib import Path

_NON_LETTER = re.compile(r"[^A-Za-z]")


class InputFormatError(ValueError):
    pass


# CHANGED(ENG-UX-11): the reader messages say what is wrong in plain words and what a valid file looks like (were
# e.g. "p.fa:1: sequence data before first '>' header", "p.fa:1: non-ASCII character"). teammate's C++ mirrors them
# (cpp/src/seqio.cpp). See overnight/integration/CHANGES.md
NON_ASCII = "non-ASCII character: sequence files must be plain text (is this a binary or word-processor file?)"
NON_ASCII_ID = "non-ASCII character in the sequence ID: use only ASCII letters, digits and punctuation in IDs"
EMPTY_HEADER = "empty FASTA header: a '>' line without a sequence ID (every record needs a name after '>')"
_BEFORE = "sequence data before the first '>' header: not a FASTA file? (each FASTA record starts with a '>name' line"
BEFORE_HEADER_PROTEIN = _BEFORE + "; a CLUSTAL file must start with 'CLUSTAL')"
BEFORE_HEADER_DNA = _BEFORE + ")"


def duplicate_message(where, rid: str) -> str:
    return f"{where}: duplicate sequence ID {rid!r}: every sequence needs its own ID"


def open_error_message(path: str, e: OSError) -> str:
    """CHANGED(ENG-UX-4): new. "cannot open 'x.fa': no such file" (or: permission denied, it is a folder, or the
    system's reason) for an input that cannot be opened (was "[Errno 2] No such file or directory: 'x.fa'")."""
    if os.path.isdir(path):
        why = "it is a folder, not a file"
    elif e.errno == errno.ENOENT:
        why = "no such file"
    elif e.errno in (errno.EACCES, errno.EPERM):
        why = "permission denied"
    else:
        why = e.strerror or str(e)
    return f"cannot open '{path}': {why}"


# CHANGED(ENG-h): text inputs were opened as strict ASCII (a UTF-8 BOM failed with "'ascii' codec
# can't decode byte 0xef", no file name); now a BOM is skipped and other non-ASCII input is an
# InputFormatError naming the file. See overnight/engine-fixes/CHANGES.md
@contextmanager
def _text(path):
    """Lines of a text input: BOM skipped; \\n, \\r\\n and \\r all end a line."""
    # CHANGED(R2-E6): was errors="strict": a bad byte failed the decoder on its 8 kB chunk before the line was
    # counted ("near line 1"); now bad bytes become surrogates and the line check names the right line.
    # See overnight/integration/CHANGES.md
    fh = open(path, encoding="utf-8-sig", errors="surrogateescape")
    try:
        yield _ascii_lines(fh, path)
    finally:
        fh.close()


def _ascii_lines(fh, path):
    lineno = 0
    try:
        for line in fh:
            lineno += 1
            if not line.isascii():
                # CHANGED(FZ-2): was an error for any non-ASCII byte; now a header's description (after the ID) may
                # hold any text, as in PAL2NAL (e.g. ">a Krüppel-like factor"). Non-ASCII IDs and sequence data stay
                # errors. See overnight/integration/CHANGES.md
                if line.startswith(">"):
                    head = line[1:].split(None, 1)
                    if head and head[0].isascii():
                        yield ">" + head[0] + "\n"
                        continue
                    raise InputFormatError(f"{path}:{lineno}: {NON_ASCII_ID}")
                raise InputFormatError(f"{path}:{lineno}: {NON_ASCII}")
            yield line
    except UnicodeDecodeError:
        raise InputFormatError(f"{path}: near line {lineno + 1}: not an ASCII text file") from None


def read_fasta(path: str | Path, letters_only: bool = False) -> list[tuple[str, str]]:
    """Return [(id, sequence)]. ID = first whitespace token of the header.
    Whitespace inside sequence lines is removed. With letters_only (used for DNA,
    as PAL2NAL does), every non-letter is removed (gaps, digits, '*', ...) and
    lines starting with '#' are skipped."""
    records: list[tuple[str, list[str]]] = []
    with _text(path) as fh:  # CHANGED(ENG-h)
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                parts = line[1:].split()
                if not parts:
                    raise InputFormatError(f"{path}:{lineno}: {EMPTY_HEADER}")
                records.append((parts[0], []))
            elif line.strip():
                if letters_only and line.startswith("#"):
                    continue
                if not records:
                    raise InputFormatError(f"{path}:{lineno}: "
                                           + (BEFORE_HEADER_DNA if letters_only else BEFORE_HEADER_PROTEIN))
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
    with _text(path) as fh:  # CHANGED(ENG-h)
        first = next(fh, "")
        lineno = 1
        while first and not first.strip():  # CHANGED(FZ-1): leading blank lines (PAL2NAL reads such files as CLUSTAL)
            first = next(fh, "")
            lineno += 1
        if not first.startswith("CLUSTAL"):
            raise InputFormatError(f"{path}: not a CLUSTAL file")
        for lineno, line in enumerate(fh, lineno + 1):
            line = line.rstrip("\r\n")
            # Skip blank lines and conservation lines (start with whitespace).
            # CHANGED(FZ-1): also lines starting with '#' (PAL2NAL's block-mark lines; was read as a sequence named
            # '#') and further CLUSTAL header lines, as PAL2NAL. See overnight/integration/CHANGES.md
            if not line.strip() or line[0] in " \t#" or line.startswith("CLUSTAL"):
                continue
            parts = line.split()
            if len(parts) < 2:
                # CHANGED(ENG-UX-11): was "<file>: malformed CLUSTAL line: '...'" (no line number, no explanation)
                raise InputFormatError(f"{path}:{lineno}: malformed CLUSTAL line {line!r}: each line of a block must "
                                       f"hold a sequence name, then (after spaces) its aligned sequence")
            name, chunk = parts[0], parts[1]
            if name not in seqs:
                order.append(name)
                seqs[name] = []
            seqs[name].append(chunk)
    return [(n, "".join(seqs[n])) for n in order]


# CHANGED(DEC-4): new. PAL2NAL's '#' block marks (used by -blockonly and shown in clustal/codon output); was not read.
# See overnight/integration/CHANGES.md
def clustal_mask(path: str | Path) -> str:
    """The mark line of a CLUSTAL alignment, one character per alignment column, read the way PAL2NAL v14 reads it:
    after each sequence line, the next line that is not a sequence line (the conservation line, or a blank line) is
    padded with spaces to that sequence line's length and the part under its sequence is appended. Only '#' has a
    meaning (PAL2NAL's -blockonly keeps the '#' columns). A block whose last sequence line ends the file adds
    nothing, as in PAL2NAL."""
    with _text(path) as fh:
        text = "".join(fh)
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    while lines and lines[-1] == "":
        lines.pop()  # Perl's split drops trailing empty fields
    mask: list[str] = []
    get, tmplen, idspc, sublen = False, 0, 0, 0
    for line in lines:
        if not line.startswith("CLUSTAL") and line[:1] not in ("", " ", "\t", "#") and not line[:1].isspace():
            line = line.rstrip()
            parts = line.split()
            tmplen = len(line)
            m = re.match(r"\S+\s+", line)
            if m:
                idspc = m.end()
            sublen = len(parts[1]) if len(parts) > 1 else 0
            get = True
        elif get:
            line += " " * (tmplen - len(line))
            mask.append(line[idspc:idspc + sublen])
            get = False
    return "".join(mask)


def read_alignment(path: str | Path) -> list[tuple[str, str]]:
    """Read a protein alignment, auto-detecting CLUSTAL vs FASTA."""
    return read_clustal(path) if is_clustal(path) else read_fasta(path)  # CHANGED(FZ-1): same rule as is_clustal


def _check_unique(ids: list[str], path) -> None:
    seen = set()
    for i in ids:
        if i in seen:
            raise InputFormatError(duplicate_message(path, i))
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
    with _text(path) as fh:  # CHANGED(ENG-h)
        for lineno, line in enumerate(fh, 1):
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(chunks)
                parts = line[1:].split()
                if not parts:
                    raise InputFormatError(f"{path}:{lineno}: {EMPTY_HEADER}")
                rid, chunks = parts[0], []
                if rid in seen:
                    raise InputFormatError(duplicate_message(path, rid))
                seen.add(rid)
            elif line.strip():
                if rid is None:
                    raise InputFormatError(f"{path}:{lineno}: {BEFORE_HEADER_PROTEIN}")
                chunks.append("".join(line.split()))
    if rid is not None:
        yield rid, "".join(chunks)


def is_clustal(path: str | Path) -> bool:
    """True if the first non-blank line starts with CLUSTAL (CHANGED(FZ-1): was the first line)."""
    with _text(path) as fh:  # CHANGED(ENG-h)
        for line in fh:
            if line.strip():
                return line.startswith("CLUSTAL")
    return False


_LETTERS = bytes(range(65, 91)) + bytes(range(97, 123))
# CHANGED(ENG-h): new helpers so a lone \r also ends a line in DNA files (classic Mac).
_LINES = re.compile(rb"[^\r\n]*(?:\r\n|\r|\n)|[^\r\n]+")


def _has_lone_cr(fh) -> bool:
    """True if the first 64 kB of the (seekable) binary file have a \\r not followed by
    \\n (old Mac line endings); position kept. Checking only the start keeps \\n and
    \\r\\n files on the fast path."""
    pos = fh.tell()
    try:
        c = fh.read(1 << 16)
        if c.endswith(b"\r"):
            c += fh.read(1)
        return c.count(b"\r") != c.count(b"\r\n") + c.endswith(b"\r")
    finally:
        fh.seek(pos)


def _cr_lines(fh):
    for raw in fh:
        yield from _LINES.findall(raw)
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
        # CHANGED(ENG-h): was \n-only lines, BOM not skipped, '#' before the first header an
        # error; now a BOM is skipped, a lone \r also ends a line, and '#' lines before the
        # first header are skipped (as read_fasta does). See overnight/engine-fixes/CHANGES.md
        with open(path, "rb") as fh:
            offset, start = 0, None
            if fh.read(len(codecs.BOM_UTF8)) == codecs.BOM_UTF8:
                offset = len(codecs.BOM_UTF8)
            else:
                fh.seek(0)
            # CHANGED(R3-E3): was "<file>: empty FASTA header" / "non-ASCII sequence ID" / "sequence data before
            # first '>' header" with no line number; now "<file>:<line>: ..." as for protein files.
            # See overnight/integration/CHANGES.md
            for lineno, line in enumerate(_cr_lines(fh) if _has_lone_cr(fh) else fh, 1):
                if line.startswith(b">"):
                    if start is not None:
                        self._add(fno, start, offset)
                    parts = line[1:].split()
                    if not parts:
                        raise InputFormatError(f"{path}:{lineno}: {EMPTY_HEADER}")
                    try:
                        rid = parts[0].decode("ascii")
                    except UnicodeDecodeError:
                        raise InputFormatError(f"{path}:{lineno}: {NON_ASCII_ID}") from None
                    if rid in self.pos:
                        raise InputFormatError(duplicate_message(", ".join(self.paths), rid))
                    self.pos[rid] = len(self.ids)
                    self.ids.append(rid)
                    start = offset + len(line)
                elif start is None and line.strip() and not line.startswith(b"#"):
                    raise InputFormatError(f"{path}:{lineno}: {BEFORE_HEADER_DNA}")
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
