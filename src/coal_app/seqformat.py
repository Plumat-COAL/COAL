"""FASTA / CLUSTAL reading rules shared by the app's own readers.

CHANGED(R3-A1): new. The app read files as plain UTF-8 and wanted the first non-blank line to start with '>' or
'CLUSTAL', so files the engine accepts (UTF-8 BOM, '#' lines before the first header, lone CR line ends) were
"not FASTA or CLUSTAL" or gave no IDs. This module mirrors src/coal/seqio.py (the app must not import engine
internals); tests_app/test_round3.py runs the same edge-case files through both and compares.
See overnight/integration/CHANGES.md

Rules (as src/coal/seqio.py):
- a UTF-8 BOM is skipped; \\n, \\r\\n and a lone \\r all end a line;
- the file is CLUSTAL if its first line starts with "CLUSTAL", otherwise FASTA;
- FASTA: blank lines are ignored, the ID is the first word after '>', lines starting with '#' are skipped
  (the engine skips them in DNA files; in a protein file it rejects them, so skipping changes nothing for
  any input the engine accepts); DNA keeps only letters, protein drops whitespace;
- CLUSTAL: lines that are blank or start with whitespace are skipped; name = first word, data = second word.
"""
from __future__ import annotations

import re

_NON_LETTER = re.compile(r"[^A-Za-z]")


def _lines(path: str):
    """Lines without their line end; BOM skipped, any line end (\\n, \\r\\n, \\r) accepted."""
    with open(path, encoding="utf-8-sig", errors="replace", newline=None) as fh:
        for line in fh:
            yield line.rstrip("\r\n")


def is_clustal_line(first: str) -> bool:
    return first.startswith("CLUSTAL")


def iter_records(path: str, dna: bool = False, ids_only: bool = False):
    """Yield (id, sequence) in file order (sequence is "" with ids_only). Raises ValueError for a file that is
    neither FASTA nor CLUSTAL (no header at all)."""
    it = _lines(path)
    first = next(it, None)
    if first is None:
        raise ValueError(f"{path}: not FASTA or CLUSTAL")
    if is_clustal_line(first):
        seqs: dict[str, list[str]] = {}
        for line in it:
            if not line.strip() or line[0] in " \t":
                continue
            parts = line.split()
            if len(parts) >= 2 or ids_only:
                seqs.setdefault(parts[0], []).append(parts[1] if len(parts) >= 2 else "")
        for k, v in seqs.items():
            yield k, ("" if ids_only else "".join(v))
        return
    name, buf = None, []

    def lines():
        yield first
        yield from it

    for line in lines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf)
            parts = line[1:].split()
            name, buf = (parts[0] if parts else ""), []
        elif not line.strip() or line.startswith("#"):
            continue
        elif name is None:
            raise ValueError(f"{path}: not FASTA or CLUSTAL")
        elif not ids_only:
            buf.append(_NON_LETTER.sub("", line) if dna else "".join(line.split()))
    if name is None:
        raise ValueError(f"{path}: not FASTA or CLUSTAL")
    yield name, "".join(buf)


def read_records(path: str, dna: bool = False) -> list[tuple[str, str]]:
    return list(iter_records(path, dna=dna))


def read_ids(path: str) -> list[str]:
    """IDs in file order; [] for a file that is neither FASTA nor CLUSTAL (as before)."""
    ids: list[str] = []
    try:
        for i, _ in iter_records(path, ids_only=True):
            ids.append(i)
    except ValueError:
        pass
    return ids
