"""FASTA / CLUSTAL reading rules shared by the app's own readers.

This module mirrors src/coal/seqio.py (the app must not import engine internals), so the app accepts every file
the engine accepts (UTF-8 BOM, '#' lines before the first header, lone CR line ends); tests_app/test_round3.py runs
the same edge-case files through both and compares.

Rules (as src/coal/seqio.py):
- a UTF-8 BOM is skipped; \\n, \\r\\n and a lone \\r all end a line;
- the file is CLUSTAL if its first non-blank line starts with "CLUSTAL", otherwise FASTA;
- FASTA: blank lines are ignored, the ID is the first word after '>', lines starting with '#' are skipped
  (the engine skips them in DNA files; in a protein file it rejects them, so skipping changes nothing for
  any input the engine accepts); DNA keeps only letters, protein drops whitespace;
- CLUSTAL: lines that are blank, start with whitespace or '#', or start with CLUSTAL are skipped;
  name = first word, data = second word.
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


FAST_READ_MAX = 512 * 1024 * 1024  # FASTA files up to this size are read in one go (fast path)


def _fasta_whole(text: str, path: str, dna: bool, ids_only: bool, want):
    """The FASTA rules below, on a whole file at once: records are cut at every line that starts with '>' by one
    str.split, and only the records wanted are assembled (no Python work per line). The
    line-by-line reader takes ~2 s for a 125 MB codon alignment of 100,000 genes, of which the preview needs 100."""
    if text.startswith(">"):
        start = 0
    else:
        j = text.find("\n>")
        if j < 0:
            raise ValueError(f"{path}: not FASTA or CLUSTAL")
        start = j + 1
        for line in text[:start].split("\n"):
            if line.strip() and not line.startswith("#"):
                raise ValueError(f"{path}: not FASTA or CLUSTAL")
    for piece in text[start + 1:].split("\n>"):
        nl = piece.find("\n")
        parts = (piece if nl < 0 else piece[:nl]).split(None, 1)
        name = parts[0] if parts else ""
        if want is not None and name not in want:
            continue
        if ids_only:
            yield name, ""
            continue
        body = "" if nl < 0 else piece[nl:]
        if "\n#" in body:
            body = "\n".join(ln for ln in body.split("\n") if not ln.startswith("#"))
        if dna:
            yield name, _NON_LETTER.sub("", body)
        else:
            flat = body.replace("\n", "")
            yield name, ("".join(flat.split()) if (" " in flat or "\t" in flat or "\f" in flat or "\v" in flat)
                         else flat)


def iter_records(path: str, dna: bool = False, ids_only: bool = False, want=None):
    """Yield (id, sequence) in file order (sequence is "" with ids_only). Raises ValueError for a file that is
    neither FASTA nor CLUSTAL (no header at all).
    want (a set of IDs, or None for all): only those records are yielded, and the others are
    skipped without building their sequences (the alignment preview reads 100 of 100,000 records)."""
    import os
    try:
        small = os.path.getsize(path) <= FAST_READ_MAX
    except OSError:
        small = False
    if small:
        with open(path, encoding="utf-8-sig", errors="replace", newline=None) as fh:
            text = fh.read()
        head = text[:65536].lstrip()
        first_line = head.split("\n", 1)[0]
        if first_line and not is_clustal_line(first_line):
            yield from _fasta_whole(text, path, dna, ids_only, want)
            return
        del text, head
    it = _lines(path)
    first = next(it, None)
    while first is not None and not first.strip():  # skip leading blank lines
        first = next(it, None)
    if first is None:
        raise ValueError(f"{path}: not FASTA or CLUSTAL")
    if is_clustal_line(first):
        seqs: dict[str, list[str]] = {}
        for line in it:
            if not line.strip() or line[0] in " \t#" or line.startswith("CLUSTAL"):
                continue
            parts = line.split()
            if len(parts) >= 2 or ids_only:
                seqs.setdefault(parts[0], []).append(parts[1] if len(parts) >= 2 else "")
        for k, v in seqs.items():
            if want is None or k in want:
                yield k, ("" if ids_only else "".join(v))
        return
    name, buf = None, []
    keep = True

    def lines():
        yield first
        yield from it

    for line in lines():
        if line.startswith(">"):
            if name is not None and keep:
                yield name, "".join(buf)
            parts = line[1:].split()
            name, buf = (parts[0] if parts else ""), []
            keep = want is None or name in want
        elif not keep:
            continue
        elif not line.strip() or line.startswith("#"):
            continue
        elif name is None:
            raise ValueError(f"{path}: not FASTA or CLUSTAL")
        elif not ids_only:
            buf.append(_NON_LETTER.sub("", line) if dna else "".join(line.split()))
    if name is None:
        raise ValueError(f"{path}: not FASTA or CLUSTAL")
    if keep:
        yield name, "".join(buf)


def read_records(path: str, dna: bool = False, want=None) -> list[tuple[str, str]]:
    return list(iter_records(path, dna=dna, want=want))


def read_ids(path: str) -> list[str]:
    """IDs in file order; [] for a file that is neither FASTA nor CLUSTAL (as before)."""
    ids: list[str] = []
    try:
        for i, _ in iter_records(path, ids_only=True):
            ids.append(i)
    except ValueError:
        pass
    return ids
