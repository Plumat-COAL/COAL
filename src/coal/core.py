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
from functools import lru_cache
from itertools import accumulate

from .codon import ANCHOR_START, NORMAL, START, accepted, codon_ok, fast_code, residue_regex
from .tables import SUPPORTED, get_code

GAP = "-"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY*X")
AS_X = set("BZJOU")  # accepted but treated as 'X' (with a warning), like PAL2NAL
AS_X_NAMES = {"U": "selenocysteine", "O": "pyrrolysine", "B": "D or N", "Z": "E or Q", "J": "I or L"}
_VALID = re.compile(r"[^-ACDEFGHIKLMNPQRSTVWY*XBZJOU0-9]")
_DIGIT = re.compile(r"[0-9]")
_MORE_THAN_3_RES = re.compile(r"[^-0-9](?:[-0-9]*[^-0-9]){3}")  # _note_ambiguous's quick test
_AS_X = re.compile("[" + "".join(sorted(AS_X)) + "]")
_AS_X_TABLE = str.maketrans({c: "X" for c in AS_X})
_CODON = re.compile(r"...", re.S)
_GAPCHAR = re.compile("-")
_TIMES3 = (3).__mul__
ANCHOR = 10
PREFIX_TOKENS = 8
RUN_WINDOW, RUN_MIN, RUN_MAX = 24, 4, 8
_NON_ACGT = re.compile("[^ACGT]")


MISMATCH_WARNING = " does not correspond to "  # used to filter the mismatch warnings (-nomismatch)


class CoalError(Exception):
    """An error COAL reports to the user. `code` is a stable machine-readable reason (used in the run report);
    `seq_id` is set when one sequence failed (so --keep-going can skip it); `index` is its position in the input."""

    def __init__(self, message: str, code: str = "error", seq_id: str | None = None,
                 index: int | None = None):
        super().__init__(message)
        self.code, self.seq_id, self.index = code, seq_id, index

    def __reduce__(self):  # keep code/seq_id/index (and any extra attributes) across worker processes
        return (CoalError, (str(self), self.code, self.seq_id, self.index), dict(self.__dict__))


def check_table(table) -> None:
    """Raise CoalError `unknown_table` unless `table` is a supported NCBI genetic-code number (tables.SUPPORTED).
    The public API's check: tables.get_code itself raises ValueError, as an internal helper."""
    if table not in SUPPORTED:
        raise CoalError(f"unknown genetic code table {table!r}: the supported NCBI tables are "
                        f"{', '.join(map(str, SUPPORTED))}", code="unknown_table")


class Threaded:
    """Per-sequence result: the matched DNA and the column layout. (A plain class, not a dataclass: importing
    dataclasses costs ~20 ms of start-up, which shows on pipelines that run COAL once per small file.)"""
    __slots__ = ("seq_id", "row", "dna", "warnings", "mismatch_cols", "notices")

    def __init__(self, seq_id: str, row: str, dna: str, warnings: list[str] | None = None,
                 mismatch_cols: set[int] | None = None, notices: list[str] | None = None):
        self.seq_id = seq_id
        self.row = row            # aligned protein row (upper case, '.' -> '-')
        self.dna = dna            # matched DNA, original case
        self.warnings = [] if warnings is None else warnings
        self.mismatch_cols = set() if mismatch_cols is None else mismatch_cols  # 1-based alignment columns
        self.notices = [] if notices is None else notices  # probably-wrong data, apart from per-codon warnings

    def __repr__(self):
        return f"Threaded({', '.join(f'{k}={getattr(self, k)!r}' for k in self.__slots__)})"

    def __eq__(self, other):
        if not isinstance(other, Threaded):
            return NotImplemented
        return all(getattr(self, k) == getattr(other, k) for k in self.__slots__)

    __hash__ = None


# ---------------------------------------------------------------- pairing

DB_TAGS = frozenset({"sp", "tr", "gi", "ref", "gb", "emb", "dbj", "pdb", "lcl", "gnl", "pir", "prf",
                     "tpg", "tpe", "tpd", "bbs", "bbm", "pat"})


def _strip_version(x: str) -> str:
    """AT1G01010.1 -> AT1G01010 (a trailing '.<digits>')."""
    head, dot, tail = x.rpartition(".")
    return head if dot and head and tail.isdigit() else x


def id_keys(rid: str) -> set[str]:
    """Normalised forms of an ID used for safe pairing: version suffix removed, lower case; for pipe headers
    (sp|P12345|NAME_HUMAN) every field except database tags is tried as well."""
    keys = {_strip_version(rid).lower()}
    fields = [f for f in rid.split("|") if f]
    if len(fields) > 1:
        keys.update(_strip_version(f).lower() for f in fields if f.lower() not in DB_TAGS)
    return keys


def normalised_pairing(prot_ids: list[str], dna_ids: list[str], pos: dict) -> tuple[list[int], list[str]] | None:
    """DNA record per protein by normalised ID (id_keys), or None. Accepted only if every protein has exactly one
    candidate (an exact ID match wins) and the pairing is one-to-one. Returns (indices, per-protein 'id' or
    'normalised')."""
    by_key: dict[str, set[int]] = {}
    for i, d in enumerate(dna_ids):
        for k in id_keys(d):
            by_key.setdefault(k, set()).add(i)
    idx, modes, used = [], [], set()
    for p in prot_ids:
        if p in pos:
            j, m = pos[p], "id"
        else:
            cands = set().union(*(by_key.get(k, ()) for k in id_keys(p)))
            if len(cands) != 1:
                return None
            j, m = next(iter(cands)), "normalised"
        if j in used:
            return None
        used.add(j)
        idx.append(j)
        modes.append(m)
    return idx, modes


def _records(n: int, what: str) -> str:
    return f"{n} {what}" + ("" if n == 1 else "s")


def count_mismatch_message(prot_ids: list[str], dna_ids: list[str], pos: dict) -> str:
    """The count error with what to do; names an example record when the IDs show which
    records are extra (every protein ID is a DNA ID) or which proteins have no DNA (exactly that many are missing)."""
    p, d = len(prot_ids), len(dna_ids)
    head = f"number of sequences differs (protein: {p}, DNA: {d}): each protein needs exactly one DNA record"
    if d > p:
        eg = ""
        if all(x in pos for x in prot_ids):
            pset = set(prot_ids)
            extra = next(x for x in dna_ids if x not in pset)
            eg = f" (e.g. '{extra}', which is not in the protein alignment)"
        return (f"{head}; remove the {_records(d - p, 'extra DNA record')}{eg} or check you used the right files")
    missing = [x for x in prot_ids if x not in pos]
    eg = f" (e.g. '{missing[0]}')" if missing and len(missing) == p - d else ""
    return (f"{head}; {p - d} protein" + (" has" if p - d == 1 else "s have") + f" no DNA record{eg}: add "
            f"the missing DNA, remove those proteins from the alignment, or check you used the right files")


_NUC = frozenset("ACGTUN")
_PROTEIN_ONLY = frozenset("EFIJLOPQZ")  # amino-acid letters that are not IUPAC nucleotide codes (X excluded: masking)


def _letters(s: str) -> list[str]:
    return [c for c in s.upper() if "A" <= c <= "Z"]


SNIFF_MIN = 20  # letters needed before COAL guesses what kind of sequence a text is


def looks_like_dna(text: str) -> bool:
    """True if the text has at least SNIFF_MIN letters and at least 90% of them are
    A/C/G/T/U/N (a protein alignment never is)."""
    letters = _letters(text)
    return len(letters) >= SNIFF_MIN and 10 * sum(c in _NUC for c in letters) >= 9 * len(letters)


def looks_like_protein(text: str) -> bool:
    """True if the text has at least SNIFF_MIN letters and at least 10% of them are amino acids that are not
    nucleotide codes."""
    letters = _letters(text)
    return len(letters) >= SNIFF_MIN and 10 * sum(c in _PROTEIN_ONLY for c in letters) >= len(letters)


def row_length_message(label: str, ncol: int, n: int, pid: str, rows: str) -> str:
    """The error for protein rows of different lengths: the usual causes are the two files
    given the other way round (the first file is DNA) and a protein file that is not aligned. `rows`: the text of
    the rows seen (the first and the offending one)."""
    if looks_like_dna(rows):
        return (f"{label} looks like DNA, not a protein alignment: give the protein alignment first and the DNA "
                f"second (coal PROTEIN_ALIGNMENT DNA_FASTA)")
    return (f"protein alignment rows have different lengths ({ncol} and {n}, at {pid}): is {label} an aligned "
            f"protein file, e.g. from MAFFT or MUSCLE? Every row of an alignment has the same length, gaps included")


def inconsistent_message(seq_id: str, detail: str, row: str, dna: str) -> str:
    """The per-sequence inconsistency error, plus a hint when the inputs look swapped."""
    msg = f"{seq_id}: protein and DNA are inconsistent ({detail})"
    if looks_like_dna(row):
        msg += ("; its protein row looks like DNA: are the two files the other way round? "
                "(coal PROTEIN_ALIGNMENT DNA_FASTA)")
    elif looks_like_protein(dna):
        msg += "; its DNA looks like protein: the second file must hold the coding DNA, not protein sequences"
    return msg


def order_pairing_warning(prot_ids: list[str], dna_ids: list[str]) -> str:
    """The file-order notice. Not printed (when no protein ID is a DNA ID, COAL stays silent as PAL2NAL does); the
    run report's id_mismatch hint says it instead. Kept for the C++ comparison."""
    pairs = ", ".join(f"{p} <-> DNA {d}" if i == 0 else f"{p} <-> {d}"
                      for i, (p, d) in enumerate(zip(prot_ids[:2], dna_ids[:2])))
    more = ", ..." if len(prot_ids) > 2 else ""
    return (f"IDs don't match, so sequences were paired by file order (protein {pairs}{more}). "
            f"Check this is what you want.")


def pair_ids(prot_ids: list[str], dna_ids: list[str], dna_pos: dict | None = None, warn=None,
             strict: bool = False):
    """Pair proteins with DNA records. Returns (DNA index per protein, mode, per-protein modes or None).

    mode "id": every protein ID is a DNA ID. "order": file order, PAL2NAL's rule. "normalised": by normalised ID
    (version suffix, case, pipe fields; see normalised_pairing), used ONLY where PAL2NAL's file order would pair
    a protein with a different DNA record than the one its normalised ID names (otherwise the order pairing is
    the same pairing, and COAL keeps PAL2NAL's output and messages exactly).
    `warn` receives the notices; `strict` makes a partial exact-ID match an error."""
    pos = dna_pos if dna_pos is not None else {d: i for i, d in enumerate(dna_ids)}
    if len(prot_ids) != len(dna_ids):
        raise CoalError(count_mismatch_message(prot_ids, dna_ids, pos), code="count_mismatch")
    if not prot_ids:
        raise CoalError("no sequences in input", code="no_sequences")
    if all(p in pos for p in prot_ids):
        return [pos[p] for p in prot_ids], "id", None
    hits = sum(p in pos for p in prot_ids)
    if hits and strict:
        # An error with --strict, also when normalised IDs would pair them (otherwise a warning, below)
        p, d = next((p, d) for p, d in zip(prot_ids, dna_ids) if p != d)
        raise CoalError(f"only {hits} of {len(prot_ids)} protein IDs are DNA IDs, so sequences would be paired "
                        f"by order (e.g. protein '{p}' with DNA '{d}'); make the IDs match, or run without "
                        f"--strict to pair by order as PAL2NAL does", code="partial_id_match")
    norm = normalised_pairing(prot_ids, dna_ids, pos)
    if norm is not None and norm[0] != list(range(len(prot_ids))):
        # Safe ID pairing, limited to where it matters: PAL2NAL's file order would give this protein the DNA of
        # another gene. Deliberate difference from PAL2NAL v14 (definition of a codon alignment: a protein's codons
        # must come from its own DNA). The C++ engine mirrors it.
        idx, modes = norm
        moved = [k for k, j in enumerate(idx) if j != k]
        k = moved[0]
        if warn:
            nn = modes.count("normalised")
            warn(f"{nn} of {len(prot_ids)} protein IDs match a DNA ID only after normalising (version suffix, case "
                 f"or database fields), so sequences were paired by those IDs, not by file order as PAL2NAL "
                 f"pairs them: by order, {len(moved)} of {len(prot_ids)} proteins would get the DNA of another gene "
                 f"(e.g. by file order, protein '{prot_ids[k]}' would get DNA '{dna_ids[k]}'; COAL paired it with "
                 f"'{dna_ids[idx[k]]}')")
        return idx, "normalised", modes
    if hits and warn:
        # Some but not all IDs match: pair by order, with a warning.
        p, d = next((p, d) for p, d in zip(prot_ids, dna_ids) if p != d)
        warn(f"only {hits} of {len(prot_ids)} protein IDs are DNA IDs; pairing sequences by order "
             f"instead (e.g. protein '{p}' gets DNA '{d}') (--strict makes this an error)")
    return list(range(len(prot_ids))), "order", None


def pair_sequences(prot_ids: list[str], dna_ids: list[str], dna_pos: dict | None = None,
                   warn=None, strict: bool = False) -> list[int]:
    """Counts must match; DNA record number for each protein (see pair_ids).
    `warn` is called with the pairing notices; with `strict` a partial exact-ID match is an error."""
    return pair_ids(prot_ids, dna_ids, dna_pos, warn=warn, strict=strict)[0]


# ---------------------------------------------------------------- tokens
# A token is (kind, value, column):  ("r", aa, col, mode) residue, ("n", N, col) N bases,
# ("w", 3|N, col) wildcard of that many bases (fallback mode).

def _tokens(row: str, start_mode: int, col0: int = 1):
    if not (_DIGIT.search(row) or _AS_X.search(row)):
        # Only ordinary residues: one comprehension, then the first Met's mode
        toks = [("r", ch, col0 + i, NORMAL) for i, ch in enumerate(row) if ch != GAP]
        for q, t in enumerate(toks):
            if t[1] != "X":
                if t[1] == "M" and start_mode != NORMAL:
                    toks[q] = ("r", "M", t[2], start_mode)
                break
        return toks
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
    n = 3 * len(toks)
    for t in toks:
        if t[0] != "r":
            n += t[1] - 3
    return n


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
    """A normalised DNA sequence plus lazily built 3-frame translations and an index of their 10-residue windows."""
    __slots__ = ("table", "dna", "plain", "_frames", "_windows")

    def __init__(self, table: int, dna: str):
        self.table, self.dna = table, dna
        self.plain = _NON_ACGT.search(dna) is None
        self._frames = None
        self._windows = None

    def windows(self) -> list[dict[str, int]]:
        """Per frame: every ANCHOR-residue window of the translation -> the DNA position of its leftmost occurrence
        (built once per sequence, in O(len(dna))). The anchor fallback looks its chunks up here instead of scanning
        the DNA for each of them (which is quadratic in the gene length)."""
        if self._windows is None:
            out = []
            for f, tr in enumerate(self.frames()):
                nw = len(tr) - ANCHOR + 1
                # Windows in decreasing order of position, so the dict keeps the leftmost occurrence of each word
                out.append({tr[k:k + ANCHOR]: f + 3 * k for k in range(nw - 1, -1, -1)})
            self._windows = out
        return self._windows

    def find_plain(self, word: str) -> int:
        """Leftmost DNA position (any frame) whose next len(word) codons translate to `word`, or -1; len(word) ==
        ANCHOR. Same answer as _find for tokens that are all ordinary residues (no X, digit or first-Met mode): a
        triplet matches such a residue exactly when it translates to it (frames())."""
        best = -1
        for d in self.windows():
            s = d.get(word, -1)
            if s >= 0 and (best < 0 or s < best):
                best = s
        return best

    def frames(self) -> list[str]:
        if self._frames is None:
            if self.plain:
                get = fast_code(self.table).__getitem__
                self._frames = ["".join(map(get, _CODON.findall(self.dna, f))) for f in range(3)]
            else:
                # DNA with N/IUPAC codes: a triplet translates to the one residue whose accepted() set contains it,
                # else '?', so the frame search stays exact (candidates are verified by _matches_at). Without it,
                # such DNA would need regex scans (~quadratic when many anchors fail).
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


def _plain_word(toks) -> str | None:
    """The residues of an anchor's tokens as a word for _Target.find_plain, or None when a token is not an ordinary
    residue (a frameshift digit, X, or a first Met with its start-codon mode) or there are not exactly ANCHOR."""
    if len(toks) != ANCHOR:
        return None
    for t in toks:
        if t[0] != "r" or t[1] == "X" or t[3] != NORMAL:
            return None
    return "".join([t[1] for t in toks])


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
    run = _plain_run(toks)  # also for DNA with N/IUPAC codes
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

_CHUNK = re.compile("(?:-*[^-]){%d}" % ANCHOR)  # a chunk: up to and including its 10th non-gap character


def _anchors(row: str) -> list[str]:
    """Split the aligned row into chunks of 10 non-gap characters (PAL2NAL rule); what is left at the end of the row
    (fewer than 10 non-gap characters, possibly none) is a chunk too, and a final chunk shorter than 10 *characters*
    is merged into the previous one."""
    pre = _CHUNK.findall(row)
    tail = row[sum(map(len, pre)):]
    if tail:
        pre.append(tail)
    if pre and len(pre[-1]) < ANCHOR:
        pre = pre[:-2] + ["".join(pre[-2:])] if len(pre) >= 2 else []
    return pre


def _anchor_tokens(toks, anchors: list[str]):
    """The tokens of each anchor (_tokens(anchor, ANCHOR_START if first else NORMAL, its first column)), taken from
    the row's tokens: every non-gap character is one token, and the only mode that differs from the row's is the
    first Met's (START in the row; ANCHOR_START in the first anchor, NORMAL in a later one)."""
    out, i = [], 0
    for k, a in enumerate(anchors):
        j = i + len(a) - a.count(GAP)
        at = toks[i:j]  # a new list
        for q, t in enumerate(at):
            if t[0] == "r" and t[3] != NORMAL:
                at[q] = (t[0], t[1], t[2], ANCHOR_START if k == 0 else NORMAL)
        out.append(at)
        i = j
    return out


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


def limit_blas_threads() -> None:
    """COAL uses no BLAS. Importing numpy otherwise starts a BLAS thread pool sized to the machine, each thread with
    its own buffers, in every worker process: memory that grows with workers x cores. Called only in processes COAL
    owns (the command line and its workers), never on import, so a program using COAL as a library is unaffected.
    A value the user set wins."""
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(var, "1")


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


# Fewest checkable (A/C/G/T-only) codons that make a placement without any matching anchor trustworthy when the
# DNA is longer than the gene.
MIN_CHECKABLE = 10


AMBIGUOUS_MAX_RESIDUES = 3  # proteins this short get a notice when they fit in several places


def _note_ambiguous(res: Threaded, table: int, dna: str, row: str, toks=None) -> None:
    """A protein of at most 3 residues often fits in several places of a longer DNA (e.g. inside a 5' UTR); like
    PAL2NAL, COAL takes the leftmost, and says so in a notice. stdout is unchanged."""
    # One regex that stops at the 4th residue (any character but a gap or a frameshift digit), so the common case
    # does not scan the whole row
    if _MORE_THAN_3_RES.search(row):
        return  # the common case, decided without building tokens
    toks = toks if toks is not None else _tokens(row, START)
    nres = sum(1 for t in toks if t[0] == "r")
    if nres > AMBIGUOUS_MAX_RESIDUES:
        return
    n = sum(1 for st in range(len(dna) - _length(toks) + 1) if _matches_at(table, dna, st, toks))
    if n > 1:
        res.notices.append(f"{res.seq_id}: placement is ambiguous: {n} possible positions in the DNA for this "
                           f"{nres}-residue protein; used the leftmost (as PAL2NAL)")


def thread_one(seq_id: str, row: str, raw_dna: str, table: int, strict: bool = False) -> Threaded:
    bad = _VALID.search(row)
    if bad:
        chars = sorted(set(_VALID.findall(row)))
        raise CoalError(f"{seq_id}: invalid protein alignment characters {chars} (first at column {bad.start() + 1}): "
                        f"a protein alignment may hold only amino-acid letters, '-' or '.' gaps, '*' stops and "
                        f"frameshift digits", code="invalid_protein", seq_id=seq_id)

    res = Threaded(seq_id, row, "")
    for m in _AS_X.finditer(row):
        # PAL2NAL says "U unknown AA type. Taken as 'X'", but U is a known amino acid (selenocysteine); same
        # behaviour and stdout, clearer wording (also for B/Z/J/O).
        res.warnings.append(f"pepAlnPos {m.start() + 1}: {m.group()} ({AS_X_NAMES[m.group()]}) is not checked; "
                            f"any codon is accepted")
        # -nomismatch removes it, as PAL2NAL v14 removes every column it printed a message for
        res.mismatch_cols.add(m.start() + 1)

    dna = raw_dna.upper().replace("U", "T")
    residues = row.replace(GAP, "")
    if not residues:
        raise CoalError(f"{seq_id}: protein sequence is empty (its alignment row has only gaps): remove it from both "
                        f"files", code="empty", seq_id=seq_id)

    # 1. exact fit (quick check at offset 0 first: the common, clean case)
    if not _DIGIT.search(residues) and _quick_offset0(table, dna, residues):
        res.dna = raw_dna[:3 * len(residues)]
        _note_ambiguous(res, table, dna, row)
        return res
    target = _Target(table, dna)
    toks = _tokens(row, START)
    s = _find(target, toks)
    if s >= 0:
        res.dna = raw_dna[s:s + _length(toks)]
        _note_ambiguous(res, table, dna, row, toks)
        return res

    # 2. anchor fallback
    anchors = _anchors(row)
    if not anchors:
        raise CoalError(inconsistent_message(seq_id, "no match found", row, raw_dna), code="inconsistent",
                        seq_id=seq_id)
    pattern = []
    found_any = False
    expect = -1  # where this anchor starts if the previous anchors' matches are the real ones
    for at in _anchor_tokens(toks, anchors):
        # Each anchor only needs to occur somewhere in the DNA. A chunk of ANCHOR ordinary residues (the usual one)
        # is checked against the translated frames: at the expected position first, else looked up in the index of
        # translated windows (O(1) after a one-off O(len(dna)) build). Any other chunk (with a frameshift digit, an
        # X, the first Met, or the longer merged last chunk) is checked at the expected position and otherwise
        # searched for by _find. Scanning the DNA for every chunk would be quadratic in the gene length; the result
        # (found or not) is the same as _find's for each chunk.
        word = _plain_word(at)
        if word is not None:
            if expect >= 0 and target.frames()[expect % 3].startswith(word, expect // 3):
                hit = expect
            else:
                hit = target.find_plain(word)
        elif expect >= 0 and _matches_at(table, dna, expect, at):
            hit = expect
        else:
            hit = _find(target, at)
        n = _length(at)
        if hit < 0:
            at = [("w", 3 if t[0] == "r" else t[1], t[2]) for t in at]
            expect = expect + n if expect >= 0 else -1
        else:
            found_any = True
            expect = hit + n
        pattern += at
    s = _find(target, [("n", t[1], t[2]) if t[0] == "w" else t for t in pattern])
    if s < 0:
        raise CoalError(inconsistent_message(seq_id, "no match found", row, raw_dna), code="inconsistent",
                        seq_id=seq_id)
    total = _length(toks)
    res.dna = raw_dna[s:s + total]

    # mismatch report: first residue Met -> ATG/CTG/GTG, others -> genetic code
    nchk = nbad_chk = 0  # codons with only A/C/G/T(U) bases, and those of them that mismatch
    nres = 0
    if not _DIGIT.search(residues):
        # No frameshift digit (the usual case): codon i is dna[s+3i:s+3i+3]. A codon matches an ordinary residue
        # exactly when it translates to it (see _Target.frames), so translate the matched DNA once and compare it
        # with the residues (the first Met and X apart) instead of testing every codon one by one.
        nres = len(toks)
        nd = dna[s:s + total]
        if target.plain or _NON_ACGT.search(nd) is None:
            tr = "".join(map(fast_code(table).__getitem__, _CODON.findall(nd)))
            amb = frozenset()
        else:
            get = _iupac_code(table).get
            tr = "".join([get(c, "?") for c in _CODON.findall(nd)])
            amb = frozenset(m.start() // 3 for m in _NON_ACGT.finditer(nd))  # codons with an ambiguous base
        nchk = nres - len(amb)
        aas = residues.translate(_AS_X_TABLE)
        bad = [i for i, (x, y) in enumerate(zip(tr, aas)) if x != y and y != "X"]
        if aas[0] == "M":  # the first Met: start codons too
            ok0 = nd[:3] in accepted(table, "M", ANCHOR_START)
            if ok0 and bad and bad[0] == 0:
                del bad[0]
            elif not ok0 and (not bad or bad[0] != 0):
                bad.insert(0, 0)
        for i in bad:
            c = toks[i][2]
            res.warnings.append(f"pepAlnPos {c}: {row[c - 1]}{MISMATCH_WARNING}{res.dna[3 * i:3 * i + 3]}")
            res.mismatch_cols.add(c)
            nbad_chk += i not in amb
    else:
        pos, first = 0, True
        for t in toks:
            if t[0] == "n":
                pos += t[1]
                continue
            aa, c = t[1], t[2]
            nres += 1
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
    # The wildcard pattern fitting does not prove the protein and DNA are related: > half the codons mismatching
    # is probably-wrong data, and no anchor matching gives a notice. Codons containing any base other than
    # A/C/G/T(U) count in neither the numerator nor the denominator (no rule if none is checkable), so a correctly
    # placed gene with many NNN (masked / low coverage) codons is not flagged.
    namb = nres - nchk
    # An error only with --strict (strict=True). By default COAL writes the same row as PAL2NAL (which never refuses
    # on mismatches) and adds a notice saying the data is probably wrong.
    if nchk and 2 * nbad_chk > nchk:
        if strict:
            raise CoalError(inconsistent_message(seq_id, f"{nbad_chk} of {nchk} checkable codons do not match the "
                                                         f"protein", row, raw_dna),
                            code="inconsistent", seq_id=seq_id)
        res.notices.append(f"{seq_id}: protein and DNA are probably unrelated: {nbad_chk} of {nchk} checkable codons "
                           f"do not match the protein. The row was written as PAL2NAL writes it, but it is "
                           f"probably the wrong DNA or the wrong genetic code (--strict makes this an error)")
        return res
    # When no anchor matched and almost every codon is ambiguous, the all-wildcard pattern sits at offset 0, so a
    # 5' UTR can end up inside the codons (output shifted). So when no anchor matched, the DNA is longer than the
    # gene (the offset is not forced) and fewer than MIN_CHECKABLE codons can be checked, the placement is a guess.
    if not found_any and len(dna) > _length(toks) and nchk < MIN_CHECKABLE:
        # "the gene has only N codons" when none is ambiguous (not "only 7 of 7 codons have no ambiguous bases")
        few = (f"only {nchk} of {nres} codons have no ambiguous bases" if namb else f"the gene has only {nres} codons")
        if strict:
            raise CoalError(inconsistent_message(seq_id, f"no 10-residue segment of the protein was found in the "
                                                         f"DNA and {few}: too few to tell where "
                                                         f"{'the gene' if namb else 'it'} starts in the longer DNA",
                                                 row, raw_dna),
                            code="inconsistent", seq_id=seq_id)
        res.notices.append(f"{seq_id}: position of the gene in the DNA is a guess: no 10-residue segment of the "
                           f"protein was found in the DNA and {few}, too few to tell where "
                           f"{'the gene' if namb else 'it'} starts in the longer DNA. The codons were taken from the "
                           f"start of the DNA as PAL2NAL does, so the row may be shifted (--strict makes this an error)")
        return res
    if not found_any:
        # The same checkable-codon counts as the rule above (ambiguous codons are not mismatches), plus the number
        # of ambiguous codons
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
    # Deliberate difference from PAL2NAL, which uses the standard code's stops TAA/TAG/TGA (and TAR), upper case
    # only, whatever -codontable says (table 2 would remove TGA = Trp and keep AGA = stop; lower-case 'taa' would
    # never be removed).
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
                    table: int = 1) -> list[str]:  # table: the -nogap stops
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
    # tfastx/tfasty frameshift notation, as PAL2NAL v14 reads it: '\' = 1 base (digit 1), and '/' + gaps + residue =
    # a gap, the gaps, then 2 bases (digit 2).
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
    check_table(table)  # CoalError `unknown_table`, not tables.get_code's ValueError
    proteins = [(pid, prepare_row(pid, aln)) for pid, aln in proteins]
    lengths = {len(s) for _, s in proteins}
    if len(lengths) > 1:
        # the message names the likely cause (files swapped, or protein not aligned)
        first = proteins[0][1]
        k = next(i for i, (_, r) in enumerate(proteins) if len(r) != len(first))
        raise CoalError(row_length_message("the protein alignment", len(first), len(proteins[k][1]), proteins[k][0],
                                           first + proteins[k][1]), code="row_lengths")
    pair_warn: list[str] = []  # partial-ID-match warning added to warnings
    idx = pair_sequences([p[0] for p in proteins], [d[0] for d in dnas], warn=pair_warn.append, strict=strict)
    items = [thread_one(pid, row, dnas[j][1], table, strict) for (pid, row), j in zip(proteins, idx)]
    rows = build_alignment(items, nogap=nogap, nomismatch=nomismatch, table=table)
    warnings = [f"{it.seq_id} {w}" for it in items for w in it.warnings
                if not (nomismatch and MISMATCH_WARNING in w)]  # as PAL2NAL
    warnings += [n for it in items for n in it.notices]  # notices included
    warnings = pair_warn + warnings
    return [(it.seq_id, r) for it, r in zip(items, rows)], warnings
