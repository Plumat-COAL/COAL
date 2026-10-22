"""Minimal readers for FASTA and CLUSTAL protein alignments and FASTA DNA."""
from __future__ import annotations

import codecs
import errno
import os
import re
from array import array
from contextlib import contextmanager
from pathlib import Path

from .core import CoalError

_NON_LETTER = re.compile(r"[^A-Za-z]")


class InputFormatError(CoalError, ValueError):
    """An input that is not FASTA or CLUSTAL (or not ASCII text): a CoalError with code `input_format`, as the
    `coal` command reports it; still a ValueError for callers that catch that."""

    def __init__(self, message: str):
        CoalError.__init__(self, message, code="input_format")

    def __reduce__(self):
        return (InputFormatError, (str(self),), dict(self.__dict__))


# The reader messages say what is wrong in plain words and what a valid file looks like. The C++ engine mirrors them
# (cpp/src/seqio.cpp).
NON_ASCII = "non-ASCII character: sequence files must be plain text (is this a binary or word-processor file?)"
NON_ASCII_ID = "non-ASCII character in the sequence ID: use only ASCII letters, digits and punctuation in IDs"
# NUL and the other C0 control bytes (except tab, LF, VT, FF, CR) never occur in a sequence file: a binary file, or
# UTF-16 text
CONTROL_CHAR = "control character: sequence files must be plain text (is this a binary or word-processor file?)"
_CONTROL = re.compile(r"[\x00-\x08\x0e-\x1f]")
_CONTROL_BYTES = re.compile(rb"[\x00-\x08\x0e-\x1f]")
EMPTY_HEADER = "empty FASTA header: a '>' line without a sequence ID (every record needs a name after '>')"
_BEFORE = "sequence data before the first '>' header: not a FASTA file? (each FASTA record starts with a '>name' line"
BEFORE_HEADER_PROTEIN = _BEFORE + "; a CLUSTAL file must start with 'CLUSTAL')"
BEFORE_HEADER_DNA = _BEFORE + ")"


def duplicate_message(where, rid: str) -> str:
    if len(rid) > 200:  # (raised before the IDs are known to the message layer: shortened here)
        rid = rid[:100] + "..." + rid[-20:]
    return f"{where}: duplicate sequence ID {rid!r}: every sequence needs its own ID"


def open_error_message(path: str, e: OSError) -> str:
    """The error for an input that cannot be opened: "cannot open 'x.fa': no such file" (or: permission denied, it
    is a folder, or the system's reason)."""
    if os.path.isdir(path):
        why = "it is a folder, not a file"
    elif e.errno == errno.ENOENT:
        why = "no such file"
    elif e.errno in (errno.EACCES, errno.EPERM):
        why = "permission denied"
    else:
        why = e.strerror or str(e)
    return f"cannot open '{path}': {why}"


# Text inputs: a UTF-8 BOM is skipped and other non-ASCII input is an InputFormatError naming the file.
@contextmanager
def _text(path):
    """Lines of a text input: BOM skipped; \\n, \\r\\n and \\r all end a line."""
    # Not errors="strict": that fails the decoder on its 8 kB chunk before the line is counted ("near line 1").
    # Bad bytes become surrogates instead, so the line check names the right line.
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
            if _CONTROL.search(line):
                raise InputFormatError(f"{path}:{lineno}: {CONTROL_CHAR}")
            if not line.isascii():
                # A header's description (after the ID) may hold any text, as in PAL2NAL (e.g. ">a Krüppel-like
                # factor"). Non-ASCII IDs and sequence data are errors.
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
    with _text(path) as fh:
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
    with _text(path) as fh:
        first = next(fh, "")
        lineno = 1
        while first and not first.strip():  # leading blank lines (PAL2NAL reads such files as CLUSTAL)
            first = next(fh, "")
            lineno += 1
        if not first.startswith("CLUSTAL"):
            raise InputFormatError(f"{path}: not a CLUSTAL file")
        for lineno, line in enumerate(fh, lineno + 1):
            line = line.rstrip("\r\n")
            # Skip blank lines and conservation lines (start with whitespace).
            # Also lines starting with '#' (PAL2NAL's block-mark lines) and further CLUSTAL header lines, as PAL2NAL.
            if not line.strip() or line[0] in " \t#" or line.startswith("CLUSTAL"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise InputFormatError(f"{path}:{lineno}: malformed CLUSTAL line {line!r}: each line of a block must "
                                       f"hold a sequence name, then (after spaces) its aligned sequence")
            name, chunk = parts[0], parts[1]
            if name not in seqs:
                order.append(name)
                seqs[name] = []
            seqs[name].append(chunk)
    return [(n, "".join(seqs[n])) for n in order]


# PAL2NAL's '#' block marks (used by -blockonly and shown in clustal/codon output).
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
    return read_clustal(path) if is_clustal(path) else read_fasta(path)


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
    with _text(path) as fh:
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
    """True if the first non-blank line starts with CLUSTAL."""
    with _text(path) as fh:
        for line in fh:
            if line.strip():
                return line.startswith("CLUSTAL")
    return False


_LETTERS = bytes(range(65, 91)) + bytes(range(97, 123))
# Helpers so a lone \r also ends a line in DNA files (classic Mac).
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
        # A BOM is skipped, a lone \r also ends a line, and '#' lines before the first header are skipped (as
        # read_fasta does).
        with open(path, "rb") as fh:
            offset, start = 0, None
            if fh.read(len(codecs.BOM_UTF8)) == codecs.BOM_UTF8:
                offset = len(codecs.BOM_UTF8)
            else:
                fh.seek(0)
            # Errors are "<file>:<line>: ..." as for protein files.
            for lineno, line in enumerate(_cr_lines(fh) if _has_lone_cr(fh) else fh, 1):
                if line.startswith(b">"):
                    if start is not None:
                        self._add(fno, start, offset)
                    if _CONTROL_BYTES.search(line):  # (the sequence data are checked when they are read: get())
                        raise InputFormatError(f"{path}:{lineno}: {CONTROL_CHAR}")
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
        if len(raw) != b - a:
            # Never another record's bytes: a file that shrank or changed after it was indexed is an error (the
            # same words as coal-cpp)
            raise CoalError(f"{self.paths[fno]}: could not read the DNA of record {i + 1} (the file changed after "
                            f"it was read, or cannot be read)", code="file_error")
        m = _CONTROL_BYTES.search(raw)
        if m:
            raise InputFormatError(f"{self.paths[fno]}:{self._line_of(fno, a + m.start())}: {CONTROL_CHAR}")
        if b"#" in raw:
            raw = b"".join(l for l in raw.splitlines(True) if not l.startswith(b"#"))
        return raw.translate(None, _NON_LETTER_BYTES).decode("ascii")

    def _line_of(self, fno: int, offset: int) -> int:
        """The 1-based line number of a byte offset (for an error message only: reads the file up to there)."""
        with open(self.paths[fno], "rb") as fh:
            head = fh.read(offset)
        n = head.count(b"\n")
        if n == 0:
            n = head.count(b"\r")  # old Mac line ends
        return n + 1

    def close(self) -> None:
        for fh in self._fh.values():
            fh.close()
        self._fh.clear()
