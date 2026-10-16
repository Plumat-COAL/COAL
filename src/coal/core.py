"""Core conversion: protein alignment + coding DNA -> codon alignment.

Behaviour follows PAL2NAL v14 (used as a black-box/behavioural reference; this
is an independent implementation):

1. Locate: find the leftmost position in the DNA (any frame) where the whole
   protein fits codon by codon. Extra sequence before/after (UTRs, stop codon,
   poly-A) is therefore ignored. The first Met may use any start codon of the
   chosen genetic code. A digit N in the alignment is a frameshift: N bases.
2. Fallback: if there is no exact fit, split the protein into ~10-residue
   anchors; anchors found anywhere in the DNA are kept, the others are replaced
   by wildcards; locate that pattern instead and report every residue whose
   codon does not match as a warning.
3. Thread: residues -> their codons, gaps -> '-', each column padded to the
   widest entry (3, or 3*ceil(N/3) for frameshift digits).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import lru_cache

from .codon import ANCHOR_START, NORMAL, START, accepted, codon_ok, fast_code, residue_regex
from .tables import get_code

GAP = "-"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY*X")
AS_X = set("BZJOU")  # accepted but treated as 'X' (with a warning), like PAL2NAL
_VALID = re.compile(r"[^-ACDEFGHIKLMNPQRSTVWY*XBZJOU0-9]")
_DIGIT = re.compile(r"[0-9]")
_AS_X = re.compile("[" + "".join(sorted(AS_X)) + "]")
_CODON = re.compile(r"...", re.S)
_STOP_FOR_NOGAP = re.compile(r"[UT]A[AGR]|[TU]GA")  # PAL2NAL -nogap rule (case-sensitive)
ANCHOR = 10
PREFIX_TOKENS = 8
RUN_WINDOW, RUN_MIN, RUN_MAX = 24, 4, 8
_NON_ACGT = re.compile("[^ACGT]")


class CoalError(Exception):
    pass


@dataclass
class Threaded:
    """Per-sequence result: the matched DNA and the column layout."""
    seq_id: str
    row: str                 # aligned protein row (upper case, '.' -> '-')
    dna: str                 # matched DNA, original case
    warnings: list[str] = field(default_factory=list)
    mismatch_cols: set[int] = field(default_factory=set)   # 1-based alignment columns


# ---------------------------------------------------------------- pairing

def pair_sequences(prot_ids: list[str], dna_ids: list[str]) -> list[int]:
    """Counts must match; match by ID if every protein ID is a DNA ID, else by order."""
    if len(prot_ids) != len(dna_ids):
        raise CoalError(
            f"number of sequences differs (protein: {len(prot_ids)}, DNA: {len(dna_ids)})"
        )
    if not prot_ids:
        raise CoalError("no sequences in input")
    pos = {d: i for i, d in enumerate(dna_ids)}
    if all(p in pos for p in prot_ids):
        return [pos[p] for p in prot_ids]
    return list(range(len(prot_ids)))


# ---------------------------------------------------------------- tokens
# A token is (kind, value, column):  ("r", aa, col, mode) residue, ("n", N, col) N bases,
# ("w", 3|N, col) wildcard of that many bases (fallback mode).

def _tokens(row: str, start_mode: int, col0: int = 1):
    toks = []
    seen_specific = False   # PAL2NAL: "first Met" = an M before any non-X residue
    for i, ch in enumerate(row):
        col = col0 + i
        if ch == GAP:
            continue
        if ch.isdigit():
            toks.append(("n", int(ch), col))
            continue
        aa = "X" if ch in AS_X else ch
        mode = start_mode if (aa == "M" and not seen_specific) else NORMAL
        if aa != "X":
            seen_specific = True
        toks.append(("r", aa, col, mode))
    return toks


def _length(toks) -> int:
    return sum(3 if t[0] == "r" else t[1] for t in toks)


def _token_regex(table: int, t) -> str:
    if t[0] == "r":
        return residue_regex(table, t[1], t[2])
    return "[A-Z]{%d}" % t[1] if t[1] else ""


def _matches_at(table: int, dna: str, s: int, toks) -> bool:
    pos = s
    for t in toks:
        if t[0] == "r":
            if t[1] != "X" and dna[pos:pos + 3] not in accepted(table, t[1], t[3]):
                return False
            pos += 3
        else:
            pos += t[1]
    return pos <= len(dna)


class _Target:
    """A normalised DNA sequence plus lazily built 3-frame translations."""
    __slots__ = ("table", "dna", "plain", "_frames")

    def __init__(self, table: int, dna: str):
        self.table, self.dna = table, dna
        self.plain = _NON_ACGT.search(dna) is None
        self._frames = None

    def frames(self) -> list[str]:
        if self._frames is None:
            get = fast_code(self.table).__getitem__
            self._frames = ["".join(map(get, _CODON.findall(self.dna, f))) for f in range(3)]
        return self._frames


def _plain_run(toks):
    """Longest run (<= RUN_MAX) of ordinary residues among the first tokens.
    Returns (dna_offset_of_run, run_string) or None."""
    best, cur, cur_off, off = None, [], 0, 0
    for t in toks[:RUN_WINDOW]:
        ok = t[0] == "r" and t[1] != "X" and t[3] == NORMAL
        if ok:
            if not cur:
                cur_off = off
            cur.append(t[1])
            if best is None or len(cur) > len(best[1]):
                best = (cur_off, "".join(cur))
            if len(cur) >= RUN_MAX:
                break
        else:
            cur = []
        off += 3 if t[0] == "r" else t[1]
    return best if best and len(best[1]) >= RUN_MIN else None


@lru_cache(maxsize=4096)
def _prefix_regex(table: int, key: tuple):
    return re.compile("(?=" + "".join(_token_regex(table, t) for t in key) + ")")


def _find(target: _Target, toks) -> int:
    """Leftmost start where all tokens match consecutively, or -1."""
    dna, table = target.dna, target.table
    total = _length(toks)
    last = len(dna) - total
    if last < 0:
        return -1
    run = _plain_run(toks) if target.plain else None
    if run:
        # Every token has a fixed length, so the run sits at a fixed offset from
        # the start: find the run in the translated frames, then verify.
        off, word = run
        cands = []
        for f, tr in enumerate(target.frames()):
            k = tr.find(word)
            while k >= 0:
                st = f + 3 * k - off
                if 0 <= st <= last:
                    cands.append(st)
                k = tr.find(word, k + 1)
        for st in sorted(cands):
            if _matches_at(table, dna, st, toks):
                return st
        return -1
    key = tuple((t[0], t[1], t[3] if t[0] == "r" else None) for t in toks[:PREFIX_TOKENS])
    rx = _prefix_regex(table, key)
    for m in rx.finditer(dna):
        st = m.start()
        if st > last:
            break
        if _matches_at(table, dna, st, toks):
            return st
    return -1


# ---------------------------------------------------------------- locate

def _anchors(row: str) -> list[str]:
    """Split the aligned row into chunks of 10 non-gap characters (PAL2NAL rule);
    a final chunk shorter than 10 *characters* is merged into the previous one."""
    pre, cur, n = [], [], 0
    for i, ch in enumerate(row):
        cur.append(ch)
        if ch != GAP:
            n += 1
        if n == ANCHOR or i == len(row) - 1:
            pre.append("".join(cur))
            cur, n = [], 0
    if pre and len(pre[-1]) < ANCHOR:
        pre = pre[:-2] + ["".join(pre[-2:])] if len(pre) >= 2 else []
    return pre


def _quick_offset0(table: int, dna: str, residues: str) -> bool:
    """True if every codon from offset 0 is a plain codon for its residue.
    (False just means "use the general search", never "no match".)"""
    need = 3 * len(residues)
    if len(dna) < need:
        return False
    if residues[0] == "M" and dna[:3] not in get_code(table).starts:
        return False
    try:
        return "".join(map(fast_code(table).__getitem__, _CODON.findall(dna, 0, need))) == residues
    except KeyError:
        return False


def thread_one(seq_id: str, row: str, raw_dna: str, table: int) -> Threaded:
    bad = _VALID.search(row)
    if bad:
        chars = sorted(set(_VALID.findall(row)))
        raise CoalError(f"{seq_id}: invalid protein alignment characters {chars}")

    res = Threaded(seq_id, row, "")
    for m in _AS_X.finditer(row):
        res.warnings.append(f"pepAlnPos {m.start() + 1}: {m.group()} unknown AA type. Taken as 'X'")

    dna = raw_dna.upper().replace("U", "T")
    residues = row.replace(GAP, "")
    if not residues:
        raise CoalError(f"{seq_id}: protein sequence is empty")

    # 1. exact fit (quick check at offset 0 first: the common, clean case)
    if not _DIGIT.search(residues) and _quick_offset0(table, dna, residues):
        res.dna = raw_dna[:3 * len(residues)]
        return res
    target = _Target(table, dna)
    toks = _tokens(row, START)
    s = _find(target, toks)
    if s >= 0:
        res.dna = raw_dna[s:s + _length(toks)]
        return res

    # 2. anchor fallback
    anchors = _anchors(row)
    if not anchors:
        raise CoalError(f"{seq_id}: protein and DNA are inconsistent (no match found)")
    pattern, col = [], 1
    for k, a in enumerate(anchors):
        at = _tokens(a, ANCHOR_START if k == 0 else NORMAL, col)
        if _find(target, at) < 0:
            at = [("w", 3 if t[0] == "r" else t[1], t[2]) for t in at]
        pattern += at
        col += len(a)
    s = _find(target, [("n", t[1], t[2]) if t[0] == "w" else t for t in pattern])
    if s < 0:
        raise CoalError(f"{seq_id}: protein and DNA are inconsistent (no match found)")
    res.dna = raw_dna[s:s + _length(toks)]

    # mismatch report: first residue Met -> ATG/CTG/GTG, others -> genetic code
    pos, first = 0, True
    for t in toks:
        if t[0] == "n":
            pos += t[1]
            continue
        aa, c = t[1], t[2]
        codon = res.dna[pos:pos + 3]
        mode = ANCHOR_START if (first and aa == "M") else NORMAL
        if not codon_ok(table, codon.upper().replace("U", "T"), aa, mode):
            res.warnings.append(f"pepAlnPos {c}: {row[c - 1]} does not correspond to {codon}")
            res.mismatch_cols.add(c)
        first = False
        pos += 3
    return res


# ---------------------------------------------------------------- assemble

def build_alignment(items: list[Threaded], nogap: bool = False, nomismatch: bool = False) -> list[str]:
    ncol = len(items[0].row)
    widths = [3] * ncol
    for it in items:
        for m in _DIGIT.finditer(it.row):
            i = m.start()
            widths[i] = max(widths[i], (max(int(m.group()) - 1, 0) // 3 + 1) * 3)
    drop = set()
    if nomismatch:
        for it in items:
            drop |= {c - 1 for c in it.mismatch_cols}

    wide = any(w != 3 for w in widths)
    rows = []
    for it in items:
        out, pos, row, dna = [], 0, it.row, it.dna
        simple = not drop and not wide and not _DIGIT.search(row)
        if simple:
            # Fast path: copy runs of residues as whole DNA slices.
            for m in re.finditer(r"-+|[^-]+", row):
                n = m.end() - m.start()
                if row[m.start()] == GAP:
                    out.append("---" * n)
                else:
                    out.append(dna[pos:pos + 3 * n])
                    pos += 3 * n
        else:
            for i, ch in enumerate(row):
                w = widths[i]
                if ch == GAP:
                    piece = "-" * w
                elif ch.isdigit():
                    d = int(ch)
                    piece = dna[pos:pos + d] + "-" * (w - d)
                    pos += d
                else:
                    piece = dna[pos:pos + 3] + "-" * (w - 3)
                    pos += 3
                if i not in drop:
                    out.append(piece)
        rows.append("".join(out))

    if nogap:
        keep = []
        for k in range(0, len(rows[0]), 3):
            trip = [r[k:k + 3] for r in rows]
            if not any("-" in t or _STOP_FOR_NOGAP.search(t) for t in trip):
                keep.append(k)
        rows = ["".join(r[k:k + 3] for k in keep) for r in rows]
    return rows


def codon_align(
    proteins: list[tuple[str, str]],
    dnas: list[tuple[str, str]],
    table: int = 1,
    nogap: bool = False,
    nomismatch: bool = False,
) -> tuple[list[tuple[str, str]], list[str]]:
    """Return ([(protein_id, codon_aligned_dna)], warnings)."""
    get_code(table)  # validates the table number
    proteins = [(pid, aln.upper().replace(".", GAP)) for pid, aln in proteins]
    lengths = {len(s) for _, s in proteins}
    if len(lengths) > 1:
        raise CoalError(f"protein alignment rows have different lengths: {sorted(lengths)}")
    idx = pair_sequences([p[0] for p in proteins], [d[0] for d in dnas])
    items = [thread_one(pid, row, dnas[j][1], table) for (pid, row), j in zip(proteins, idx)]
    rows = build_alignment(items, nogap=nogap, nomismatch=nomismatch)
    warnings = [f"{it.seq_id} {w}" for it in items for w in it.warnings]
    return [(it.seq_id, r) for it, r in zip(items, rows)], warnings
