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

import os
import re
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import accumulate

from .codon import ANCHOR_START, NORMAL, START, accepted, codon_ok, fast_code, residue_regex
from .tables import get_code

GAP = "-"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY*X")
AS_X = set("BZJOU")  # accepted but treated as 'X' (with a warning), like PAL2NAL
_VALID = re.compile(r"[^-ACDEFGHIKLMNPQRSTVWY*XBZJOU0-9]")
_DIGIT = re.compile(r"[0-9]")
_AS_X = re.compile("[" + "".join(sorted(AS_X)) + "]")
_CODON = re.compile(r"...", re.S)
_GAPCHAR = re.compile("-")
_TIMES3 = (3).__mul__
ANCHOR = 10
PREFIX_TOKENS = 8
RUN_WINDOW, RUN_MIN, RUN_MAX = 24, 4, 8
_NON_ACGT = re.compile("[^ACGT]")


MISMATCH_WARNING = " does not correspond to "  # CHANGED(ENG-i2): new constant, used to filter


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
    # CHANGED(ENG-c): new field; notices that -q does not suppress. See overnight/engine-fixes/CHANGES.md
    notices: list[str] = field(default_factory=list)


# ---------------------------------------------------------------- pairing

def pair_sequences(prot_ids: list[str], dna_ids: list[str], dna_pos: dict | None = None,
                   warn=None, strict: bool = False) -> list[int]:
    """Counts must match; match by ID if every protein ID is a DNA ID, else by order.
    `dna_pos` (ID -> index) may be passed in to avoid building it again.
    `warn` is called when some (but not all) protein IDs are DNA IDs; with `strict` that is an error."""
    if len(prot_ids) != len(dna_ids):
        raise CoalError(
            f"number of sequences differs (protein: {len(prot_ids)}, DNA: {len(dna_ids)})"
        )
    if not prot_ids:
        raise CoalError("no sequences in input")
    pos = dna_pos if dna_pos is not None else {d: i for i, d in enumerate(dna_ids)}
    if all(p in pos for p in prot_ids):
        return [pos[p] for p in prot_ids]
    # CHANGED(ENG-d): was silent order pairing; now warns when some but not all IDs match.
    # See overnight/engine-fixes/CHANGES.md
    hits = sum(p in pos for p in prot_ids)
    if hits and (warn or strict):
        p, d = next((p, d) for p, d in zip(prot_ids, dna_ids) if p != d)
        text = (f"only {hits} of {len(prot_ids)} protein IDs are DNA IDs; pairing sequences by order "
                f"instead (e.g. protein '{p}' gets DNA '{d}')")
        # CHANGED(DEC-1): was always a warning (ENG-d); now an error with --strict, and the warning says so.
        # See overnight/integration/CHANGES.md
        if strict:
            raise CoalError(f"only {hits} of {len(prot_ids)} protein IDs are DNA IDs, so sequences would be paired "
                            f"by order (e.g. protein '{p}' with DNA '{d}'); make the IDs match, or run without "
                            f"--strict to pair by order as PAL2NAL does")
        warn(text + " (--strict makes this an error)")
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
            if self.plain:
                get = fast_code(self.table).__getitem__
                self._frames = ["".join(map(get, _CODON.findall(self.dna, f))) for f in range(3)]
            else:
                # CHANGED(R2-E9): new. DNA with N/IUPAC codes: a triplet translates to the one residue whose
                # accepted() set contains it, else '?', so the frame search stays exact (candidates are verified
                # by _matches_at as before). Was: no frame search at all for such DNA (regex scans, ~quadratic
                # when many anchors fail). See overnight/integration/CHANGES.md
                get = _iupac_code(self.table).get
                self._frames = ["".join([get(c, "?") for c in _CODON.findall(self.dna, f)]) for f in range(3)]
        return self._frames


@lru_cache(maxsize=None)
def _iupac_code(table: int) -> dict[str, str]:
    """Every upper-case IUPAC triplet that encodes exactly one residue (NORMAL mode) -> that residue."""
    out: dict[str, str] = {}
    for aa in get_code(table).codons:
        for trip in accepted(table, aa, NORMAL):
            out[trip] = aa
    return out


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
    run = _plain_run(toks)  # CHANGED(R2-E9): was only for DNA without N/IUPAC codes
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
    if len(residues) >= NP_MIN_CODONS:
        # numpy costs ~0.1 s to import; switch to it only once this run has
        # processed enough codons for the import to pay for itself.
        _NP_STATE["seen"] = _NP_STATE.get("seen", 0) + len(residues)
        if _NP_STATE["seen"] >= NP_AFTER_CODONS:
            np = _numpy()
            if np is not None:
                return _np_quick(np, table, dna, residues)
    try:
        return "".join(map(fast_code(table).__getitem__, _CODON.findall(dna, 0, need))) == residues
    except KeyError:
        return False


# Optional numpy acceleration of the quick check (same answer, ~10x faster on long
# sequences). Imported lazily so short runs don't pay numpy's start-up time.
_NP_STATE: dict = {}
NP_MIN_CODONS = 300        # per sequence
NP_AFTER_CODONS = 2_000_000  # per process, cumulative


def _numpy():
    if "np" not in _NP_STATE:
        if os.environ.get("COAL_NO_NUMPY"):
            _NP_STATE["np"] = None
        else:
            try:
                import numpy as np
                base = np.full(256, 255, dtype=np.uint8)
                for i, b in enumerate(b"TCAG"):
                    base[b] = i
                _NP_STATE.update(np=np, base=base)
            except ImportError:
                _NP_STATE["np"] = None
    return _NP_STATE["np"]


@lru_cache(maxsize=None)
def _np_table(table: int):
    np = _NP_STATE["np"]
    return np.frombuffer("".join(get_code(table).code.values()).encode("ascii"), dtype=np.uint8).copy()


def _np_quick(np, table: int, dna: str, residues: str) -> bool:
    need = 3 * len(residues)
    b = _NP_STATE["base"][np.frombuffer(dna[:need].encode("ascii"), dtype=np.uint8)]
    if (b == 255).any():
        return False  # non-ACGT character: let the general search decide
    b = b.reshape(-1, 3).astype(np.intp)
    aa = _np_table(table)[b[:, 0] * 16 + b[:, 1] * 4 + b[:, 2]]
    return bool((aa == np.frombuffer(residues.encode("ascii"), dtype=np.uint8)).all())


# CHANGED(R3-E1): new. Fewest checkable (A/C/G/T-only) codons that make a placement without any matching anchor
# trustworthy when the DNA is longer than the gene. See overnight/integration/CHANGES.md
MIN_CHECKABLE = 10


def thread_one(seq_id: str, row: str, raw_dna: str, table: int, strict: bool = False) -> Threaded:
    bad = _VALID.search(row)
    if bad:
        chars = sorted(set(_VALID.findall(row)))
        raise CoalError(f"{seq_id}: invalid protein alignment characters {chars}")

    res = Threaded(seq_id, row, "")
    for m in _AS_X.finditer(row):
        res.warnings.append(f"pepAlnPos {m.start() + 1}: {m.group()} unknown AA type. Taken as 'X'")
        # CHANGED(DEC-4f): was kept by -nomismatch; now removed, as PAL2NAL v14 removes every column it printed a
        # message for (found by the differential fuzz). See overnight/integration/CHANGES.md
        res.mismatch_cols.add(m.start() + 1)

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
    found_any = False
    expect = -1  # where this anchor starts if the previous anchors' matches are the real ones
    for k, a in enumerate(anchors):
        at = _tokens(a, ANCHOR_START if k == 0 else NORMAL, col)
        # CHANGED(ENG-e): was a full _find (regex compile + scan from the DNA start) for every
        # anchor (quadratic); now tries the expected position first, which proves the anchor
        # occurs, and only searches the whole DNA if it is not there. Same result.
        # See overnight/engine-fixes/CHANGES.md
        if expect >= 0 and _matches_at(table, dna, expect, at):
            hit = expect
        else:
            hit = _find(target, at)
        if hit < 0:
            at = [("w", 3 if t[0] == "r" else t[1], t[2]) for t in at]
            expect = expect + _length(at) if expect >= 0 else -1
        else:
            found_any = True
            expect = hit + _length(at)
        pattern += at
        col += len(a)
    s = _find(target, [("n", t[1], t[2]) if t[0] == "w" else t for t in pattern])
    if s < 0:
        raise CoalError(f"{seq_id}: protein and DNA are inconsistent (no match found)")
    res.dna = raw_dna[s:s + _length(toks)]

    # mismatch report: first residue Met -> ATG/CTG/GTG, others -> genetic code
    pos, first = 0, True
    nchk = nbad_chk = 0  # CHANGED(R2-E3): codons with only A/C/G/T(U) bases, and those of them that mismatch
    for t in toks:
        if t[0] == "n":
            pos += t[1]
            continue
        aa, c = t[1], t[2]
        codon = res.dna[pos:pos + 3]
        mode = ANCHOR_START if (first and aa == "M") else NORMAL
        cu = codon.upper().replace("U", "T")
        checkable = not cu.strip("ACGT")
        nchk += checkable
        if not codon_ok(table, cu, aa, mode):
            res.warnings.append(f"pepAlnPos {c}: {row[c - 1]}{MISMATCH_WARNING}{codon}")
            res.mismatch_cols.add(c)
            nbad_chk += checkable
        first = False
        pos += 3
    # CHANGED(ENG-c): was success whenever the wildcard pattern fit (unrelated protein and DNA
    # gave exit 0); now > half the codons mismatching is an error, and no anchor matching
    # gives a notice. See overnight/engine-fixes/CHANGES.md
    # CHANGED(R2-E3): was all codons counted, so a correctly placed gene with > half NNN (masked / low
    # coverage) codons was rejected and the run aborted; now codons containing any base other than
    # A/C/G/T(U) count in neither the numerator nor the denominator (no rule if none is checkable).
    # See overnight/integration/CHANGES.md
    nres = sum(1 for t in toks if t[0] == "r")
    namb = nres - nchk
    # CHANGED(DEC-1): was always an error (ENG-c, R2-E3, R3-E1); now only with --strict (strict=True). By default
    # COAL writes the same row as PAL2NAL (which never refuses on mismatches) and adds a notice saying the data is
    # probably wrong. See overnight/integration/CHANGES.md
    if nchk and 2 * nbad_chk > nchk:
        if strict:
            raise CoalError(f"{seq_id}: protein and DNA are inconsistent "
                            f"({nbad_chk} of {nchk} checkable codons do not match the protein)")
        res.notices.append(f"{seq_id}: protein and DNA are probably unrelated: {nbad_chk} of {nchk} checkable codons "
                           f"do not match the protein. The row was written as PAL2NAL writes it, but it is "
                           f"probably the wrong DNA or the wrong genetic code (--strict makes this an error)")
        return res
    # CHANGED(R3-E1): was accepted (exit 0, notice only) when no anchor matched and almost every codon was
    # ambiguous: the all-wildcard pattern then sits at offset 0, so a 5' UTR ended up inside the codons (output
    # shifted). Now, when no anchor matched, the DNA is longer than the gene (so the offset is not forced) and
    # fewer than MIN_CHECKABLE codons can be checked, it is an error. See overnight/integration/CHANGES.md
    if not found_any and len(dna) > _length(toks) and nchk < MIN_CHECKABLE:
        if strict:  # CHANGED(DEC-1): was always an error; see above
            raise CoalError(f"{seq_id}: protein and DNA are inconsistent (no 10-residue segment of the protein was "
                            f"found in the DNA and only {nchk} of {nres} codons have no ambiguous bases: too few to "
                            f"tell where the gene starts in the longer DNA)")
        res.notices.append(f"{seq_id}: position of the gene in the DNA is a guess: no 10-residue segment of the "
                           f"protein was found in the DNA and only {nchk} of {nres} codons have no ambiguous bases, "
                           f"too few to tell where the gene starts in the longer DNA. The codons were taken from the "
                           f"start of the DNA as PAL2NAL does, so the row may be shifted (--strict makes this an error)")
        return res
    if not found_any:
        # CHANGED(R3-E1): was "(K of N do not match)" counting ambiguous codons as mismatches; now the same
        # checkable-codon counts as the error rule, plus the number of ambiguous codons. See overnight/integration/CHANGES.md
        amb = f"; {namb} contain N or other ambiguous bases" if namb else ""
        res.notices.append(f"{seq_id}: no 10-residue segment of the protein was found in the DNA; "
                           f"codons were assigned by position ({nbad_chk} of {nchk} checkable codons do not "
                           f"match{amb})")
    return res


# ---------------------------------------------------------------- assemble

def column_widths(rows, ncol: int) -> list[int]:
    """Width (in nt) of each alignment column: 3, or wider where a frameshift digit
    N needs it (3 * ceil(N / 3)), as in PAL2NAL."""
    widths = [3] * ncol
    for row in rows:
        for m in _DIGIT.finditer(row):
            i = m.start()
            widths[i] = max(widths[i], (max(int(m.group()) - 1, 0) // 3 + 1) * 3)
    return widths


def expand(it: Threaded, widths: list[int], wide: bool) -> str:
    """Codon-aligned row for one sequence (before any column removal)."""
    out, pos, row, dna = [], 0, it.row, it.dna
    if not wide and not _DIGIT.search(row):
        # Fast path, all at C level: split the row at gaps; residue segment k
        # takes DNA [3*a_k, 3*b_k); joining with '---' puts one gap codon back
        # per '-' (consecutive gaps give empty segments).
        bounds = list(map(_TIMES3, accumulate(map(len, row.split(GAP)), initial=0)))
        return "---".join(map(dna.__getitem__, map(slice, bounds, bounds[1:])))
    for i, ch in enumerate(row):
        w = widths[i]
        if ch == GAP:
            out.append("-" * w)
        elif ch.isdigit():
            d = int(ch)
            out.append(dna[pos:pos + d] + "-" * (w - d))
            pos += d
        else:
            out.append(dna[pos:pos + 3] + "-" * (w - 3))
            pos += 3
    return "".join(out)


def nogap_bad_triplets(expanded: str, table: int) -> set[int]:
    """Indices of 3-nt output columns that -nogap removes for this row: a gap, or a stop codon of genetic code
    `table` (any case, U read as T; IUPAC triplets whose every reading is a stop count too)."""
    # CHANGED(FZ-5): was PAL2NAL's rule, the standard code's stops TAA/TAG/TGA (and TAR), upper case only, whatever
    # -codontable said (table 2 removed TGA = Trp and kept AGA = stop; lower-case 'taa' was never removed). Now the
    # selected code's stops, case-insensitive: a deliberate difference. See overnight/integration/CHANGES.md
    bad = {m.start() // 3 for m in _GAPCHAR.finditer(expanded)}
    up = expanded.upper().replace("U", "T")
    bad.update(m.start() // 3 for m in _stop_lookahead(table).finditer(up) if m.start() % 3 == 0)
    return bad


@lru_cache(maxsize=None)
def _stop_lookahead(table: int):
    return re.compile("(?=" + residue_regex(table, "*") + ")")


def keep_slices(widths: list[int], drop_cols, bad_triplets) -> list[tuple[int, int]]:
    """Character ranges of the expanded row that survive -nomismatch / -nogap.
    Columns are multiples of 3 nt wide, so both rules work on 3-nt units."""
    drop_cols, bad_triplets = set(drop_cols), set(bad_triplets)
    slices: list[tuple[int, int]] = []
    pos = 0
    for i, w in enumerate(widths):
        if i not in drop_cols:
            for k in range(pos, pos + w, 3):
                if k // 3 not in bad_triplets:
                    if slices and slices[-1][1] == k:
                        slices[-1] = (slices[-1][0], k + 3)
                    else:
                        slices.append((k, k + 3))
        pos += w
    return slices


def apply_slices(s: str, slices: list[tuple[int, int]]) -> str:
    return "".join(s[a:b] for a, b in slices)


def build_alignment(items: list[Threaded], nogap: bool = False, nomismatch: bool = False,
                    table: int = 1) -> list[str]:  # CHANGED(FZ-5): table is new (the -nogap stops)
    """In-memory assembly (library API). The streaming pipeline uses the same pieces."""
    widths = column_widths((it.row for it in items), len(items[0].row))
    wide = any(w != 3 for w in widths)
    rows = [expand(it, widths, wide) for it in items]
    if not (nogap or nomismatch):
        return rows
    drop = {c - 1 for it in items for c in it.mismatch_cols} if nomismatch else set()
    bad = set().union(*(nogap_bad_triplets(r, table) for r in rows)) if nogap else set()
    sl = keep_slices(widths, drop, bad)
    return [apply_slices(r, sl) for r in rows]


_TFASTX_FWD = re.compile(r"/(-*)[A-Z*]")


def prepare_row(pid: str, aln: str) -> str:
    # CHANGED(FZ-3): new: tfastx/tfasty frameshift notation, as PAL2NAL v14 reads it: '\' = 1 base (digit 1), and
    # '/' + gaps + residue = a gap, the gaps, then 2 bases (digit 2). Was rejected as invalid characters.
    # See overnight/integration/CHANGES.md
    row = aln.upper()
    if "\\" in row or "/" in row:
        row = _TFASTX_FWD.sub(r"-\g<1>2", row.replace("\\", "1"))
    return row.replace(".", GAP)


def codon_align(
    proteins: list[tuple[str, str]],
    dnas: list[tuple[str, str]],
    table: int = 1,
    nogap: bool = False,
    nomismatch: bool = False,
    strict: bool = False,
) -> tuple[list[tuple[str, str]], list[str]]:
    """Return ([(protein_id, codon_aligned_dna)], warnings). `strict`: see thread_one / pair_sequences."""
    get_code(table)  # validates the table number
    proteins = [(pid, prepare_row(pid, aln)) for pid, aln in proteins]
    lengths = {len(s) for _, s in proteins}
    if len(lengths) > 1:
        raise CoalError(f"protein alignment rows have different lengths: {sorted(lengths)}")
    pair_warn: list[str] = []  # CHANGED(ENG-d): partial-ID-match warning added to warnings
    idx = pair_sequences([p[0] for p in proteins], [d[0] for d in dnas], warn=pair_warn.append, strict=strict)
    items = [thread_one(pid, row, dnas[j][1], table, strict) for (pid, row), j in zip(proteins, idx)]
    rows = build_alignment(items, nogap=nogap, nomismatch=nomismatch, table=table)
    warnings = [f"{it.seq_id} {w}" for it in items for w in it.warnings
                if not (nomismatch and MISMATCH_WARNING in w)]  # CHANGED(ENG-i2): as PAL2NAL
    warnings += [n for it in items for n in it.notices]  # CHANGED(ENG-c): notices included
    warnings = pair_warn + warnings
    return [(it.seq_id, r) for it, r in zip(items, rows)], warnings
