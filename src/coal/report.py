"""Run report (schema "coal.report/1", read by the COAL GUI), the plain-language end-of-run summary
("coal: summary:" lines) and the genetic-code hint.

The report is built while the pipeline runs and written at the end, also when the run failed, so a person can
see why.

Per-sequence messages, statuses, hints and methods text follow the rules of the app report
(apps/coal-app/src/coal_app/core.py: build_fallback_report, mismatch_probe, id_hints, genetic_code_hint, methods_text,
viewer.locate_breaks), computed here from the engine's own results instead of parsed from stderr: the same
message codes and texts, suspect genes "failed" but written, mismatches counted only for A/C/G/T codons,
methods text counted by message code, genetic-code hint from real mismatches only and tested on all genes
(start-codon-only ties prefer code 11), sequences an aborted run never reached "not_checked" while the ones it
did check keep their status in the summary's counters (their entries say "not_checked" with a "not_written" message,
as the app's: no output holds them). The C++ engine mirrors this module (cpp/src/report.cpp).
"""
from __future__ import annotations

import os
import re
import stat
import sys
import time
from itertools import product

from . import __version__
from .checkpoint import ENGINE  # the report's "engine": the same string as the --resume checkpoint's
from .codon import IUPAC
from .tables import SUPPORTED, get_code

SCHEMA = "coal.report/1"
REPO_URL = "https://github.com/plumat-coal/COAL"
# The genetic-code hint re-tests at most HINT_SAMPLE sequences, problem genes first (50 keeps every bundled example
# whole, so the engine's report equals the app's), and drops a code as soon as it can no longer fit as many as the
# best code so far. Break locations are found for up to BREAK_SAMPLE failed genes.
HINT_SAMPLE = 50
# Genes with more residues than this are not re-tested for the hint: threading one 1,000,000-residue gene under the
# 17 other codes took 28 s of a 30 s run (the hint's value is in the other genes; the C++ engine does the same)
HINT_MAX_RESIDUES = 100_000
BREAK_SAMPLE = 200

# Names inside a sentence: common words lower case, genus names capitalised as NCBI writes them
CODE_NAMES = {
    1: "standard", 2: "vertebrate mitochondrial", 3: "yeast mitochondrial",
    4: "mold, protozoan and coelenterate mitochondrial and Mycoplasma/Spiroplasma", 5: "invertebrate mitochondrial",
    6: "ciliate, dasycladacean and Hexamita nuclear", 9: "echinoderm and flatworm mitochondrial",
    10: "euplotid nuclear", 11: "bacterial, archaeal and plant plastid", 12: "alternative yeast nuclear",
    13: "ascidian mitochondrial", 14: "alternative flatworm mitochondrial", 15: "Blepharisma nuclear",
    16: "chlorophycean mitochondrial", 21: "trematode mitochondrial",
    22: "Scenedesmus obliquus mitochondrial", 23: "Thraustochytrium mitochondrial",
}
TABLE_NAMES = CODE_NAMES  # alias

SUSPECT_CODES = ("suspect_unrelated", "suspect_position")
# The start and the end of the note before the methods text of a run with mostly problem genes
NOT_READY = "Fix the warnings first:"
AS_IT_STANDS = "Methods text for this run as it stands:"
# one character per sequence in RunReport.kinds
KIND_CHAR = {"ok": "o", "mismatch": "m", "unchecked": "u", "other": "x", "flagged": "f", "excluded": "e",
             "not_checked": "n"}
STATUSES = ("ok", "warning", "failed", "not_checked")
AS_X_NAMES = {"U": "selenocysteine", "O": "pyrrolysine", "B": "D or N", "Z": "E or Q", "J": "I or L"}
# Per-sequence flags kept next to `kinds` (one byte per sequence, 0 for a clean one): everything the summary lines,
# the methods text and the genetic-code hint read from a non-ok entry, so the entry itself (with every per-codon
# message) need not stay in memory (RunReport.flags; the C++ engine mirrors them).
F_REAL_MM = 1          # a definite (A/C/G/T) mismatch
F_FRAMESHIFT = 2       # a possible_frameshift message
F_UNRELATED = 4        # a suspect_unrelated message
F_ALL_START = 8        # every definite mismatch is a first-Met one (start_codon)
F_UNKNOWN_AA = 16      # an unknown_aa message
F_PAIRING = 32         # a suspect_pairing message
F_OTHER_CODE = 64      # a message whose code is neither unknown_aa nor suspect_pairing
F_FAILED_NOCOL = 128   # status failed with no column message: a break_event candidate
NO_CODON = b"   "      # RunReport.codons for a sequence without a definite mismatch

_MISMATCH = re.compile(r"^pepAlnPos (\d+): (\S+) does not correspond to (\S+)")
_UNKNOWN = re.compile(r"^pepAlnPos (\d+): (\S+) (?:unknown AA type|\([^)]*\) is not checked)")
_UNRELATED = re.compile(r"^protein and DNA are probably unrelated: (\d+) of (\d+) checkable codons do not match")
_GUESS = re.compile(r"^position of the gene in the DNA is a guess: .*?(?:only (\d+) of (\d+) codons have no "
                    r"ambiguous bases|the gene has only (\d+) codons)")
_NOANCHOR = re.compile(r"^no 10-residue segment of the protein was found in the DNA; codons were assigned by position "
                       r"\((\d+) of (\d+) (?:checkable codons )?do not match(?:; (\d+) contain N or other ambiguous "
                       r"bases)?\)")
_PARTIAL = re.compile(r"^only (\d+) of (\d+) protein IDs are DNA IDs; (.*)$")


def is_ambiguous_codon(codon) -> bool:
    """True if a codon has a base other than A/C/G/T(U) (N, IUPAC codes, ...): COAL can't check it."""
    return bool(codon) and bool(str(codon).upper().replace("U", "T").strip("ACGT"))


# ------------------------------------------------------------------ messages (same codes and texts as the app)

def warning_message(w: str, removed: bool = False) -> dict:
    """A per-codon warning of thread_one ("pepAlnPos N: ...") as a report message."""
    m = _MISMATCH.match(w)
    if m:
        col, aa, codon = m.groups()
        msg = {"code": "mismatch", "column": int(col), "text": f"Column {col}: amino acid {aa} does not match codon "
               f"{codon}.", "residue": aa, "codon": codon}
        if removed:
            msg["removed"] = True
        return msg
    m = _UNKNOWN.match(w)
    if m:
        col, aa = m.groups()
        what = AS_X_NAMES.get(aa)
        return {"code": "unknown_aa", "column": int(col), "residue": aa,
                "text": (f"Column {col}: '{aa}'" + (f" ({what})" if what else "") + " is not one of the 20 standard "
                         "amino acids, so COAL did not check its codon (any codon accepted).")}
    return {"code": "warning", "text": f"COAL warning: {w.strip()}"}


def notice_message(pid: str, n: str) -> dict:
    """A per-sequence notice of thread_one ("<id>: ...") as a report message."""
    rest = n[len(pid) + 2:] if n.startswith(pid + ": ") else n
    m = _UNRELATED.match(rest)
    if m:
        nbad, nchk = m.groups()
        return {"code": "suspect_unrelated",
                "text": (f"{nbad} of the {nchk} codons COAL could check do not match this protein, so this is "
                         f"almost certainly the wrong DNA for it (or the wrong genetic code). COAL still wrote "
                         f"its row, codon by position, as PAL2NAL does: do not trust it. Run COAL with --strict "
                         f"(in the app: 'Refuse probably-wrong genes') to refuse such genes.")}
    m = _GUESS.match(rest)
    if m:
        nchk, nres, short = m.groups()
        why = (f"Only {nchk} of this gene's {nres} codons have no unknown bases (N), and no "
               if short is None else f"The gene has only {short} codons, and no ")
        return {"code": "suspect_position",
                "text": (f"{why}10-residue stretch of the protein was found in the DNA, so COAL cannot tell where "
                         f"the gene starts in the longer DNA. Its row was taken from the start of the DNA (as "
                         f"PAL2NAL does) and may be shifted: do not trust it. Run COAL with --strict (in the app: "
                         f"'Refuse probably-wrong genes') to refuse such genes.")}
    m = _NOANCHOR.match(rest)
    if m:
        nbad, nres, namb = m.groups()
        amb = f"; {namb} contain N or other ambiguous bases" if namb else ""
        ask = ("Every codon COAL could check matches, so the placement looks right." if namb and int(nbad) == 0
               else "Is this the right DNA for this protein, and the right genetic code?")
        return {"code": "no_anchor",
                "text": (f"No 10-residue stretch of the protein was found in the DNA, so codons were assigned "
                         f"by position ({nbad} of {nres} checkable codons do not match{amb}). {ask}")}
    return {"code": "warning", "text": f"COAL warning: {rest.strip()}"}


def failure_message(pid: str, code: str, message: str) -> dict:
    rest = message[len(pid) + 2:] if message.startswith(pid + ": ") else message
    return {"code": code, "text": rest}


def entry_flags(entry: dict) -> tuple[int, bytes]:
    """The F_* flags of an entry and the codon of its first definite mismatch (NO_CODON if none): what the report
    keeps per sequence instead of the entry (RunReport.flags / codons)."""
    f, codon, all_start = 0, NO_CODON, True
    has_col = False
    for m in entry.get("messages") or []:
        code = m.get("code")
        if isinstance(m.get("column"), int):
            has_col = True
        if code == "mismatch":
            if not is_ambiguous_codon(m.get("codon")):
                if not f & F_REAL_MM:
                    codon = (str(m.get("codon")).encode("utf-8", "replace") + NO_CODON)[:3]
                f |= F_REAL_MM
                all_start = all_start and bool(m.get("start_codon"))
        elif code == "possible_frameshift":
            f |= F_FRAMESHIFT
        elif code == "suspect_unrelated":
            f |= F_UNRELATED
        elif code == "unknown_aa":
            f |= F_UNKNOWN_AA
        elif code == "suspect_pairing":
            f |= F_PAIRING
        if code not in ("unknown_aa", "suspect_pairing"):
            f |= F_OTHER_CODE
    if f & F_REAL_MM and all_start:
        f |= F_ALL_START
    if entry.get("status") == "failed" and not has_col:
        f |= F_FAILED_NOCOL
    return f, codon


def gene_kind(entry: dict) -> str:
    """What happened to one gene (as the app's gene_kind): ok, mismatch, unchecked, other, flagged (written but
    probably wrong), excluded (failed, no row) or not_checked."""
    st = entry["status"]
    msgs = entry.get("messages") or []
    if st == "failed":
        return "flagged" if any(m.get("code") in SUSPECT_CODES + ("suspect_pairing",) for m in msgs) else "excluded"
    if st == "warning":
        if any(m.get("code") == "mismatch" and not is_ambiguous_codon(m.get("codon")) for m in msgs):
            return "mismatch"
        if any(m.get("code") == "mismatch" for m in msgs):
            return "unchecked"
        return "other"
    return st


# The app's rules (apps/coal-app/src/coal_app/core.py real_mismatches, _fixable, _fits_now,
# _code_distance, annotate_mismatch_patterns, GENETIC_CODES), so the engine's report says what the app's would.
APP_CODE_NAMES = {
    1: "Standard (nuclear genes of most organisms)", 2: "Vertebrate mitochondrial", 3: "Yeast mitochondrial",
    4: "Mold, protozoan and coelenterate mitochondrial; Mycoplasma", 5: "Invertebrate mitochondrial",
    6: "Ciliate, dasycladacean and Hexamita nuclear", 9: "Echinoderm and flatworm mitochondrial",
    10: "Euplotid nuclear", 11: "Bacterial, archaeal and plant plastid", 12: "Alternative yeast nuclear",
    13: "Ascidian mitochondrial", 14: "Alternative flatworm mitochondrial", 15: "Blepharisma nuclear",
    16: "Chlorophycean mitochondrial", 21: "Trematode mitochondrial", 22: "Scenedesmus obliquus mitochondrial",
    23: "Thraustochytrium mitochondrial",
}


def real_mismatches(e: dict) -> list[dict]:
    return [m for m in e.get("messages") or [] if m.get("code") == "mismatch" and not is_ambiguous_codon(m.get("codon"))]


def fixable(e: dict) -> bool:
    """A gene another genetic code might fix (the app's _fixable)."""
    k = gene_kind(e)
    if k == "mismatch":
        return not any(m.get("code") == "possible_frameshift" for m in e.get("messages") or [])
    if k == "flagged":
        return any(m.get("code") == "suspect_unrelated" for m in e.get("messages") or []) or bool(real_mismatches(e))
    if k == "not_checked":
        written = any(m.get("code") == "not_written" for m in e.get("messages") or [])
        return bool(real_mismatches(e)) if written else True
    return k == "excluded"


def code_distance(a: int, b: int) -> tuple[int, int]:
    ca, cb = get_code(a), get_code(b)
    return (sum(ca.code[x] != cb.code[x] for x in ca.code), len(ca.starts ^ cb.starts))


def annotate_mismatch_pattern(row: str, e: dict) -> None:
    """The app's annotate_mismatch_patterns for one gene: a first-Met-only mismatch gets start_codon = True; a dense
    run of definite mismatches gets a possible_frameshift message."""
    index, k = {}, 0
    for c, ch in enumerate(row, 1):
        if ch not in "-.":
            index[c] = k
            k += 1
    mm = real_mismatches(e)
    idx = sorted({index[m["column"]] for m in mm if m.get("column") in index})
    if not idx:
        return
    if idx == [0] and all(m.get("residue", "").upper() == "M" for m in mm):
        for m in mm:
            m["start_codon"] = True
        return
    runs, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 2:
            cur.append(i)
        else:
            runs.append(cur)
            cur = [i]
    runs.append(cur)
    best = max(runs, key=len)
    span = best[-1] - best[0] + 1
    to_end = best[-1] >= k - 2
    if len(best) >= 3 and len(best) / span >= 0.6 and 2 * len(best) >= len(idx) and (to_end or len(best) >= 4):
        col_of = {v: c for c, v in index.items()}
        c0, c1 = col_of[best[0]], col_of[best[-1]]
        where = "to the end of the gene" if to_end else f"up to column {c1}"
        e.setdefault("messages", []).append({
            "code": "possible_frameshift", "column": c0, "end_column": c1,
            "text": (f"The mismatches form one run, from column {c0} {where}: this looks like a missing or extra "
                     f"base (a frameshift) or a small insertion/deletion in the DNA near column {c0}, not a "
                     f"genetic-code problem. Check the DNA there (a sequencing or assembly error?).")})


# ------------------------------------------------------------------ where a failed gene breaks (app viewer rule)

def _expand(codon: str) -> list[str]:
    try:
        return ["".join(c) for c in product(*(IUPAC[b] for b in codon))]
    except KeyError:
        return []


def _ok(aa, codon, code, starts, first) -> bool:
    if len(codon) < 3:
        return False
    if aa in "XBZJOU":
        return True
    exp = _expand(codon)
    return bool(exp) and all(code.get(c) == aa or (first and aa == "M" and c in starts) for c in exp)


SHIFTS = (-1, 1, -2, 2, -4, 4, -5, 5)  # +-4 / +-5 = a frameshift together with a lost or extra codon
SHIFT_WORDS = {-1: "one base is missing", -2: "two bases are missing", 1: "there is one extra base",
               2: "there are two extra bases", -4: "four bases are missing", -5: "five bases are missing",
               4: "there are four extra bases", 5: "there are five extra bases"}


def break_event(row: str, dna: str, table: int) -> dict | None:
    """The first event of the app viewer's diagnose(): where a failed gene stops matching (frameshift, missing
    DNA, or the column from which the DNA no longer spells the protein). Same algorithm, so the GUI says the same.
    The earliest credible start (not the best-scoring one, which can be a duplicated stretch further on); a shift
    must also hold for the next 20 residues and may include a whole lost or extra codon (+-4, +-5); "at column N"
    only when the old frame holds up to N and the new one from the next residue, else "near column N"; "the DNA
    ends" only after a walk where most codons matched (else no location).
    The C++ engine mirrors it (cpp/src/report.cpp)."""
    gc = get_code(table)
    code, starts = gc.code, gc.starts
    dna = re.sub(r"[^A-Za-z]", "", dna).upper().replace("U", "T")
    res = [(i + 1, a.upper()) for i, a in enumerate(row) if a.isalpha()]
    if not res or len(dna) < 3:
        return None
    want = min(len(res), 12)
    allowed = max(0, want // 4)
    need = max(min(3, len(res)), want - allowed)
    best = (0, -1)
    for off in range(0, len(dna) - 2):
        score = 0
        for k in range(want):
            if _ok(res[k][1], dna[off + 3 * k: off + 3 * k + 3], code, starts, k == 0):
                score += 1
            elif k + 1 - score > allowed:
                break
        if score >= need:  # the earliest credible start
            best = (off, score)
            break
        if score > best[1]:
            best = (off, score)
    if best[1] < need:
        return {"code": "inconsistent", "column": res[0][0],
                "text": "The DNA does not spell the start of this protein in any reading frame."}

    def hits(start, k0, n):
        return sum(_ok(a, dna[start + 3 * j: start + 3 * j + 3], code, starts, False)
                   for j, (_, a) in enumerate(res[k0:k0 + n]))
    pos, k, flags = best[0], 0, {}
    matched, last_good = 0, -1
    while k < len(res):
        col, aa = res[k]
        codon = dna[pos:pos + 3]
        if len(codon) < 3:
            if k and 5 * matched >= 4 * k:  # only after a credible walk
                return {"code": "missing_dna", "column": col, "text": f"The DNA ends before the protein does (at "
                                                                      f"column {col})."}
            return None
        if _ok(aa, codon, code, starts, k == 0):
            flags[col] = ""
            pos += 3
            k += 1
            matched += 1
            last_good = k - 1
            continue
        ahead = sum(_ok(a, dna[pos + 3 * j: pos + 3 * j + 3], code, starts, False)
                    for j, (_, a) in enumerate(res[k + 1:k + 6], 1))
        window = min(5, len(res) - k - 1)
        if window and ahead >= window - 1:
            flags[col] = "mismatch"
            pos += 3
            k += 1
            continue
        shifted = None
        n = min(10, len(res) - k - 1)
        n2 = min(20, len(res) - k - 1)
        for s in SHIFTS:
            if n <= 0:
                break
            if hits(pos + 3 + s, k + 1, n) >= max(1, int(0.8 * n)) and hits(pos + 3 + s, k + 1, n2) >= 0.8 * n2:
                shifted = s
                break
        if shifted is not None:
            exact = last_good == k - 1 and _ok(res[k + 1][1], dna[pos + 3 + shifted: pos + 6 + shifted], code,
                                               starts, False)
            return {"code": "frameshift", "column": col, "shift": shifted,
                    "text": f"Frameshift {'at' if exact else 'near'} column {col}: {SHIFT_WORDS[shifted]} in the DNA, "
                            f"so from here on the codons are read out of step."}
        flags[col] = "mismatch"
        pos += 3
        k += 1
    fl = [flags.get(col, "missing") for col, _ in res]
    bad_after, first = 0, None
    for i in range(len(fl) - 1, -1, -1):
        bad_after += fl[i] in ("mismatch", "missing")
        n = len(fl) - i
        if fl[i] == "mismatch" and n >= 10 and 2 * bad_after > n:
            first = i
    if first is not None:
        col = res[first][0]
        return {"code": "inconsistent", "column": col, "text": f"From column {col} on, the DNA no longer spells "
                                                               f"this protein."}
    return None


# ------------------------------------------------------------------ helpers

def sha256_file(path) -> str | None:
    import hashlib
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


LONG_ID = 200  # IDs longer than this are shortened on stderr (short_id); the report keeps them whole


def short_id(rid: str) -> str:
    """A very long ID as stderr shows it: the first 100 and the last 20 characters (the same in coal-cpp)."""
    return rid if len(rid) <= LONG_ID else rid[:100] + "..." + rid[-20:]


class IdShortener:
    """Replaces the IDs longer than LONG_ID in a message by short_id(): stderr lines stay readable when IDs are
    whole header lines of thousands of characters. Nothing to do (and no cost) when no ID is that long."""

    def __init__(self, rec):
        self.rec, self._seen, self.long = rec, (None, None), []

    def __call__(self, text: str) -> str:
        key = (id(self.rec.ids), id(getattr(self.rec, "dna_file_ids", None)))
        if key != self._seen:
            self._seen = key
            self.long = sorted({i for lst in (self.rec.ids, getattr(self.rec, "dna_file_ids", None)) if lst
                                for i in lst if len(i) > LONG_ID}, key=len, reverse=True)
        for rid in self.long:
            if rid in text:
                text = text.replace(rid, short_id(rid))
        return text


def _n(x: int) -> str:
    return f"{x:,}"


def _seqs(n: int) -> str:
    return f"{_n(n)} sequence" + ("" if n == 1 else "s")


def strip_version(seq_id: str) -> str:
    return re.sub(r"\.\d+$", "", seq_id)


def _fits(pid, row, dna, table) -> bool:
    """A gene fits a code: it threads, with no definite (A/C/G/T) mismatch and no suspect notice (the hint's test)."""
    from .core import CoalError, thread_one
    try:
        it = thread_one(pid, row, dna, table)
    except CoalError:
        return False  # (the app drops the whole code when one gene fails: --keep-going exits 2)
    return not (any(_MISMATCH.match(w) and not is_ambiguous_codon(w.rsplit(" ", 1)[-1]) for w in it.warnings)
                or any(notice_message(pid, x)["code"] in SUSPECT_CODES for x in it.notices))


def search_codes(tested, problems, data, table, memo=None):
    """The genetic-code hint's search: the codes under which most of the `tested` genes fit (data[k] = (id, row,
    DNA)), if better than `table` and fixing at least one of `problems`. Returns (tied codes, best first; the
    standard code is among them and first; how many fit the best; how many fit `table`) or None.
    A code is dropped as soon as it can no longer fit as many genes as the best so far (same result as testing every
    gene under every code). Shared by the end-of-run hint and the early notice; memo[(k, code)] keeps each answer, so
    the end-of-run hint does not redo the early notice's work. The C++ engine mirrors it (cpp/src/report.cpp
    search_codes)."""
    nt = len(tested)

    def bad_under(c, at_least=0):
        bad = set()
        for k in tested:
            key = (k, c)
            ok = memo.get(key) if memo is not None else None
            if ok is None:
                ok = _fits(*data[k], c)
                if memo is not None:
                    memo[key] = ok
            if not ok:
                bad.add(k)
                if len(bad) > nt - at_least:
                    return None
        return bad
    # fits under the code used: the same test (not the statuses; in a stopped run the genes never reached
    # would otherwise count as fixed by any code)
    bad_now = bad_under(table)
    fits_now = nt - len(bad_now)
    results = []
    best = fits_now + 1  # a code must fit more genes than the code used, and as many as the best one so far
    for c in SUPPORTED:
        if c == table:
            continue
        bad = bad_under(c, best)
        if bad is None:
            continue
        fit = nt - len(bad)
        fixed = sum(1 for k in problems if k in bad_now and k not in bad)
        if fixed and fit > fits_now:
            results.append((fit, c))
            best = max(best, fit)
    if not results:
        return None
    top = max(f for f, _ in results)
    tied = sorted((c for f, c in results if f == top), key=lambda c: (code_distance(table, c)[0], c))
    # The standard code 1 first when it is among the tied codes (by far the most common right answer for nuclear
    # genes), else the tied code that differs least from the one used
    standard = 1 in tied and len(tied) > 1
    if standard:
        tied = [1] + [x for x in tied if x != 1]
    return tied, standard, top, fits_now


def _do(k: int) -> str:
    return "none do" if k == 0 else ("only 1 does" if k == 1 else f"{k} do")


def early_code_notice(rows, dna_of, table: int, total: int, memo=None) -> str | None:
    """The genetic-code suggestion said before the per-codon warnings: the first HINT_SAMPLE
    genes (rows: iterable of (id, row) in input order; dna_of(k): the DNA paired with gene k) are tested under
    `table`; if any does not fit, under every other code, as the end-of-run hint does. None if no code does better.
    The C++ engine mirrors it (cpp/src/report.cpp early_code_notice)."""
    data, problems = {}, []
    for k, (pid, row) in enumerate(rows):
        if k >= HINT_SAMPLE:
            break
        if len(row) - row.count("-") > HINT_MAX_RESIDUES:
            continue
        data[k] = (pid, row, dna_of(k))
        ok = _fits(*data[k], table)
        if memo is not None:
            memo[(k, table)] = ok
        if not ok:
            problems.append(k)
    if not problems:
        return None
    found = search_codes(list(data), problems, data, table, memo)
    if found is None:
        return None
    tied, _standard, top, fits_now = found
    c, nt = tied[0], len(data)
    who = (f"all {total} sequences match" if top == total else f"{top} of {total} sequences match") if nt == total \
        else f"{top} of the first {nt} sequences match"
    return (f"wrong genetic code? Under genetic code {c} ({APP_CODE_NAMES[c]}), {who} their protein; under code "
            f"{table}, {_do(fits_now)}. Re-run with -codontable {c}.")


class RunReport:
    """Collects one run's results. keep_entries=False (no --report/--json): ok sequences are only counted.

    Per sequence the report keeps one `kinds` byte, one `flags` byte (F_*) and the first definite mismatch's codon
    (`codons`), plus the first flagged and the first excluded entry whole and the break events of up to
    BREAK_SAMPLE failed genes: everything the `coal: summary:` lines, the methods text and the genetic-code hint
    read. The per-sequence entries themselves (every per-codon message) are written to the entries file as they
    come (`attach_entries` / `open_entries_temp`: the pipeline attaches the `--resume` entries file, or a temporary
    file when a report is wanted without `-o`) and read back in order when the report is written, so memory does
    not grow with the number of warnings (a 400-codon gene of unrelated DNA has about 370). Without an entries
    file, the entries wanted for a report stay in memory (the Python API with a RunReport and no pipeline)."""

    def __init__(self, argv, pep: str, nucs, out: str | None, fmt: str, table: int, nogap: bool, nomismatch: bool,
                 keep_going: bool = False, keep_entries: bool = True, strict: bool = False, blockonly: bool = False,
                 hint: bool = True):
        self.t0 = time.perf_counter()
        self.want_hint = hint  # False (-q without --report/--json): the hint is never shown
        self.fit_memo: dict = {}  # (sequence, code) -> fits, shared with early_code_notice
        self.keep_entries = keep_entries
        self.flags = bytearray()                   # F_* per processed sequence (parallel to kinds)
        self.codons = bytearray()                  # 3 bytes per processed sequence: the first definite mismatch's codon
        self.first_flagged: dict | None = None     # the first entry whose gene_kind is "flagged" (summary lines)
        self.first_excluded: dict | None = None    # ... "excluded" (summary lines, the aborted entry views)
        self.nocol: list[int] = []                 # the first BREAK_SAMPLE failed sequences without a column message
        self.break_events: dict[int, dict] = {}    # sequence -> break_event message (finish), applied when read back
        self._pending_marks: set[int] = set()      # entries kept in memory that finish() still has to annotate
        self._store: dict[int, dict] = {}          # entries kept in memory (no entries file attached)
        self._sink = None                          # the entries file (text, append) while the run goes
        self._sink_path: str | None = None
        self._sink_owned = False                   # this report removes the entries file in close()
        self._sink_lines = 0                       # entries written to the file by this run
        self.argv = list(argv)
        self.pep, self.nucs, self.out, self.fmt, self.table = pep, list(nucs), out, fmt, table
        self.nogap, self.nomismatch, self.keep_going = nogap, nomismatch, keep_going
        self.strict, self.blockonly = strict, blockonly
        self.ids: list[str] | None = None          # protein IDs, input order
        self.dna_counts: list[int] | None = None
        self.dna_file_ids: list[str] | None = None
        self._dset: set | None = None
        self.pairing: str | None = None
        self.seq_pairing: list[str] | None = None
        self.dna_ids: list[str] | None = None      # DNA ID paired with each protein
        self.pair_hints: list[dict] = []
        self.pair_note: str | None = None
        self.partial = False
        self.kinds = bytearray()                   # gene_kind initial per processed sequence
        self._ff_cache: tuple | None = None        # (key, first failed index): _first_fail is called per entry
        self.cols: set[int] = set()                # -nomismatch: columns with a mismatch warning
        self.real_cols: set[int] = set()           # ... of which with a definite (A/C/G/T) mismatch
        self.error: dict | None = None
        self.errors: list[dict] = []
        self.output: dict | None = None
        self.gc_hint: dict | None = None
        self.gc_hint_cli: str | None = None        # the hint as the command line says it
        self.exit_code = 0
        self.resumed_from: int | None = None
        self._dict = None
        self.finished = False
        self.input_hashes: dict[int, str | None] = {}  # input number (0 = protein) -> sha256 of what was read

    # -------------------------------------------------- filled in by the pipeline
    def set_input_hashes(self, hashes: dict):
        """sha256 of inputs read through a temp copy (`-` = standard input, pipes), by input number (0 = protein):
        the argument itself cannot be read again when the report is written."""
        self.input_hashes = dict(hashes)

    def set_inputs(self, ids, dna_counts=None, dna_ids=None):
        """Protein IDs (input order); DNA records per file and DNA IDs in file order, once the DNA is read."""
        self.ids = list(ids)
        if dna_counts is not None:
            self.dna_counts = list(dna_counts)
        if dna_ids is not None:
            self.dna_file_ids = list(dna_ids)

    def set_pairing(self, mode: str, per_seq, dna_ids, notes=()):
        self.pairing, self.seq_pairing, self.dna_ids = mode, per_seq, list(dna_ids)
        for n in notes:
            m = _PARTIAL.match(n)
            if m:
                k, total, rest = m.groups()
                self.partial = True
                ex = re.search(r"protein '([^']*)' gets DNA '([^']*)'", rest)
                eg = f" (for example protein '{ex.group(1)}' got DNA '{ex.group(2)}')" if ex else ""
                self.pair_hints.append({"code": "partial_id_match", "text": (
                    f"Only {k} of {total} protein IDs were found in the DNA file, so COAL paired the sequences by "
                    f"their order in the files instead{eg}. Some proteins are probably paired with the wrong DNA: "
                    f"make the names match (or put both files in the same order) and run again.")})
            elif mode == "normalised":
                self.pair_note = n
                self.pair_hints.append({"code": "normalised_ids", "text": (
                    "Some protein IDs match their DNA only after tidying them (a '.1'-style version suffix, upper/"
                    "lower case or database fields such as 'sp|'), so COAL paired those sequences by their tidied "
                    "IDs. PAL2NAL pairs by file order here, which would have given some proteins the DNA of another "
                    "gene. " + n[0].upper() + n[1:] + ".")})

    def _pairing_of(self, k):
        return self.seq_pairing[k] if self.seq_pairing else self.pairing

    def _dna_of(self, k):
        return self.dna_ids[k] if self.dna_ids is not None and k < len(self.dna_ids) else None

    # -------------------------------------------------- the entries file (attached by the pipeline)
    def attach_entries(self, fh, path: str, owned: bool = False) -> None:
        """Write the non-ok entries to `fh` (open for appending at `path`) as they come, instead of keeping them: the
        `--resume` entries file next to the output (`owned=False`: the checkpoint removes it unless this report
        takes it over at the end of the run) or a file of this report's own. The report reads it back in order."""
        self._sink, self._sink_path, self._sink_owned, self._sink_lines = fh, path, owned, 0

    def open_entries_temp(self):
        """Attach a temporary entries file of this report's own (a report wanted for a run without -o, whose rows
        go through the temp folder too); removed in close(). Returns the open file."""
        import tempfile
        fd, path = tempfile.mkstemp(prefix=".coal_entries_", suffix=".jsonl")  # (not NamedTemporaryFile: its
        fh = os.fdopen(fd, "w", newline="\n", encoding="utf-8")                # wrapper closes the file when freed)
        self.attach_entries(fh, path, owned=True)
        return fh

    @property
    def entries_file(self) -> str | None:
        """The path of the attached entries file, if any."""
        return self._sink_path

    def finish_entries(self, keep: bool) -> None:
        """End of the pipeline: close the entries file. keep=True (the run ended, normally or with an error, and a
        report is wanted): the file stays for the report and this report removes it in close(); keep=False: a
        file of this report's own is removed now, the --resume one is left to the checkpoint."""
        if self._sink is not None:
            try:
                self._sink.close()
            except OSError:
                pass
            self._sink = None
        if self._sink_path is None:
            return
        if keep:
            self._sink_owned = True
        elif self._sink_owned:
            try:
                os.unlink(self._sink_path)
            except OSError:
                pass
            self._sink_path, self._sink_owned = None, False

    def owns_entries(self) -> bool:
        """True while this report will remove the entries file itself (close())."""
        return self._sink_path is not None and self._sink_owned

    def close(self) -> None:
        """Remove the entries file this report owns. Idempotent; the command line and the batch driver call it
        once the report is written (the report can be written several times before)."""
        self.finish_entries(keep=False)
        if self._sink_path is not None and self._sink_owned:
            try:
                os.unlink(self._sink_path)
            except OSError:
                pass
        self._sink_path, self._sink_owned = None, False

    def __del__(self):
        try:
            self.close()
        except Exception:  # noqa: BLE001 - never raise from a finalizer
            pass

    def _record(self, k: int, entry: dict) -> None:
        """Keep an entry for the report: in the entries file when one is attached (when a report or the hint wants
        it: the --resume one also feeds a resumed run's report), else in memory when a report is wanted."""
        if self._sink is not None:
            if self.keep_entries or self.want_hint:
                import json
                self._sink.write(json.dumps(entry) + "\n")
                self._sink_lines += 1
        elif self.keep_entries:
            self._store[k] = entry

    def _entries(self):
        """The non-ok entries in index order: (k, entry), from the entries file or the memory store, with the break
        events of finish() applied."""
        if self._sink_path is not None:
            import json
            if self._sink is not None:
                self._sink.flush()
            with open(self._sink_path, encoding="utf-8") as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    e = json.loads(line)
                    k = e["index"]
                    ev = self.break_events.get(k)
                    if ev is not None:
                        e["messages"].insert(0, ev)
                    yield k, e
        else:
            for k in sorted(self._store):
                yield k, self._store[k]

    def add_sequence(self, k: int, pid: str, warnings, notices, failure=None, *, mismatches=None, row=None):
        """Record sequence number k. warnings: thread_one's warning strings as printed (or, with -nomismatch, as
        they would be); with `mismatches` given, only the X/B/Z/J/O/U ones, and the per-codon mismatches come as
        (column, residue, codon) in column order (the strings are never built). notices: its notices; failure:
        (code, message) if it failed. row: the prepared protein row, so the mismatch pattern is annotated now
        (else finish() does it from the rows, for entries kept in memory). Returns the entry if not ok."""
        self._dict = None
        while len(self.kinds) < k:
            self._set_seq(len(self.kinds), "n", 0, NO_CODON)
        msgs = []
        if failure:
            status, msgs = "failed", [failure_message(pid, *failure)]
        else:
            if mismatches is None:
                mm = [warning_message(w, removed=self.nomismatch) for w in warnings if _MISMATCH.match(w)]
                others = [warning_message(w) for w in warnings if not _MISMATCH.match(w)]
                in_order = None if self.nomismatch else [warning_message(w) for w in warnings]
            else:
                mm = []
                for col, aa, codon in mismatches:
                    m = {"code": "mismatch", "column": col,
                         "text": f"Column {col}: amino acid {aa} does not match codon {codon}.", "residue": aa,
                         "codon": codon}
                    if self.nomismatch:
                        m["removed"] = True
                    mm.append(m)
                others = [warning_message(w) for w in warnings]
                in_order = None if self.nomismatch else others + mm  # thread_one's order: X/B/Z/J/O/U first
            if self.nomismatch:
                for m in mm:
                    self.cols.add(m["column"])
                    if not is_ambiguous_codon(m["codon"]):
                        self.real_cols.add(m["column"])
                msgs = mm + others
            else:
                msgs = in_order
            msgs += [notice_message(pid, n) for n in notices]
            d = self._dna_of(k)
            mispaired = self.pairing == "order" and self.partial and d is not None and d != pid
            if mispaired:
                if self._dset is None:
                    self._dset = set(self.dna_ids or ())
                in_dna = pid in self._dset
                msgs.append({"code": "suspect_pairing", "dna_id": d,
                             "text": (f"Paired with DNA '{d}' only by its position in the files"
                                      + (" (this protein's own name is also in the DNA file)" if in_dna else "")
                                      + ": check that this is the right DNA.")})
            if any(m["code"] in SUSPECT_CODES for m in msgs) or (mispaired and any(m["code"] == "mismatch"
                                                                                 for m in msgs)):
                status = "failed"
                msgs.sort(key=lambda m: (m["code"] not in SUSPECT_CODES, m["code"] != "suspect_pairing"))
            else:
                status = "warning" if msgs else "ok"
        if status == "ok":
            self._set_seq(k, "o", 0, NO_CODON)
            return None
        entry = {"id": pid, "index": k, "dna_id": self._dna_of(k), "pairing": self._pairing_of(k),
                 "status": status, "messages": msgs}
        # The mismatch pattern (start-codon-only, a run of mismatches) is marked now, with the row at hand, so the
        # entry written to the entries file is final; without the row, finish() marks the entries kept in memory
        if (self.keep_entries or self.want_hint) and any(m["code"] == "mismatch" and not is_ambiguous_codon(m["codon"])
                                                        for m in msgs):  # (with -q and no report nothing reads it)
            if row is not None:
                annotate_mismatch_pattern(row, entry)
            elif self._sink is None and self.keep_entries:
                self._pending_marks.add(k)
        self._digest(k, entry)
        self._record(k, entry)
        return entry

    def _digest(self, k: int, entry: dict) -> None:
        """Keep of a non-ok entry what the summary, the methods text and the hint need (kinds, flags, codons, the
        first flagged / excluded entry, the break_event candidates). Also used when a run is resumed."""
        kind = KIND_CHAR[gene_kind(entry)]
        f, codon = entry_flags(entry)
        self._set_seq(k, kind, f, codon)
        if kind == "f" and self.first_flagged is None:
            self.first_flagged = entry
        elif kind == "e" and self.first_excluded is None:
            self.first_excluded = entry
        if f & F_FAILED_NOCOL and len(self.nocol) < BREAK_SAMPLE:
            self.nocol.append(k)

    def _set_seq(self, k, ch, flags=0, codon=NO_CODON):
        self._ff_cache = None
        if k < len(self.kinds):
            self.kinds[k] = ord(ch)
            self.flags[k] = flags
            self.codons[3 * k:3 * k + 3] = codon
        else:
            self.kinds.append(ord(ch))
            self.flags.append(flags)
            self.codons += codon

    def _set_kind(self, k, ch):
        self._set_seq(k, ch, self.flags[k] if k < len(self.flags) else 0,
                      bytes(self.codons[3 * k:3 * k + 3]) if k < len(self.kinds) else NO_CODON)

    def fail(self, code: str, message: str, seq_id: str | None = None):
        self._dict = None
        self.error = {"code": code, "text": message, "detail": message}
        if seq_id is None:
            self.errors.append({"code": code, "text": message})

    def stopped(self):
        """The run was interrupted (signal / Ctrl-C)."""
        self._dict = None
        self.errors.insert(0, {"code": "stopped", "text": (
            "The run was stopped before it finished (COAL was interrupted, killed or ran out of time), so the "
            "sequences were not all checked and no alignment was written.")})

    def set_output(self, path, sha, n):
        self.output = {"path": path, "format": self.fmt, "sha256": sha, "sequences": n}

    # -------------------------------------------------- --resume support
    def checkpoint_state(self) -> dict:
        return {"kinds": self.kinds.decode("ascii"), "cols": sorted(self.cols), "real_cols": sorted(self.real_cols)}

    def restore(self, state, entries_path: str | None, entries_bytes: int, completed: int):
        """Continue from a checkpoint: `completed` sequences are done; the first `entries_bytes` bytes of the entries
        file hold their non-ok entries (checked by Checkpoint.load_for_resume), read once, one at a time, to rebuild
        what add_sequence kept of each (nothing is read with -q and no report, as nothing is read then)."""
        self._dict = None
        state = state or {}
        self.kinds = bytearray(state.get("kinds", "o" * completed).encode("ascii"))[:completed]
        self.flags = bytearray(len(self.kinds))
        self.codons = bytearray(NO_CODON * len(self.kinds))
        self.first_flagged = self.first_excluded = None
        self.nocol, self.break_events, self._pending_marks, self._store = [], {}, set(), {}
        self._ff_cache = None
        self.cols, self.real_cols = set(state.get("cols", ())), set(state.get("real_cols", ()))
        if (self.keep_entries or self.want_hint) and entries_path and entries_bytes > 0:
            import json
            left = entries_bytes
            with open(entries_path, "rb") as fh:
                for line in fh:
                    if left <= 0:
                        break
                    take = line[:left]
                    left -= len(take)
                    if not take.strip():
                        continue
                    e = json.loads(take)
                    if e["index"] < completed:
                        self._digest(e["index"], e)
        self.resumed_from = completed

    # -------------------------------------------------- end of run (called inside the pipeline, inputs open)
    def finish(self, rows, dna, pairing, aborted: bool = False):
        """Locate where failed genes break, mark mismatch patterns and compute the genetic-code hint, while the
        inputs are open. rows: callable giving (id, prepared row) in input order; dna: the DnaIndex; pairing: DNA
        index per protein; aborted: the run stops with an error (no output)."""
        self.finished = True
        try:
            self._finish(rows, dna, pairing, aborted)
        except Exception:  # noqa: BLE001 - extras must never hide the real result
            # Not a silent `pass`, which would hide bugs in the hint or the report: the run goes on, and COAL_DEBUG=1
            # prints the traceback.
            if os.environ.get("COAL_DEBUG"):
                import traceback
                print("coal: debug: the report's finish step failed (the run goes on):", file=sys.stderr)
                traceback.print_exc()

    def _finish(self, rows, dna, pairing, aborted):
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        # Where the first BREAK_SAMPLE failed genes without a column message break (the app viewer's rule), from
        # their rows and DNA: kept apart and put in front of the entry's messages when it is read back. Nothing to
        # locate with -q and no report (no entry is kept then, as before).
        failed = set(self.nocol) if (self.keep_entries or self.want_hint) else set()
        need = failed | self._pending_marks
        if need:
            last = max(need)
            for k, (pid, row) in enumerate(rows()):
                if k in failed:
                    ev = break_event(row, dna.get(pairing[k]), self.table)
                    if ev and self._sink_path is None and k in self._store:
                        self._store[k]["messages"].insert(0, ev)  # (the first flagged / excluded entry is this object)
                    elif ev:
                        self.break_events[k] = ev
                        for first in (self.first_flagged, self.first_excluded):
                            if first is not None and first["index"] == k:
                                first["messages"].insert(0, ev)
                if k in self._pending_marks:  # entries kept in memory and added without their row
                    annotate_mismatch_pattern(row, self._store[k])
                    f, codon = entry_flags(self._store[k])
                    self.flags[k], self.codons[3 * k:3 * k + 3] = f, codon
                if k >= last:
                    break
            self._pending_marks = set()
        # The app's hint (apps/coal-app/src/coal_app/core.py genetic_code_hint), computed here
        if not self.want_hint:
            return
        # From the per-sequence flags, one sequence at a time: which sequences another genetic code might fix
        # (fixable() of their entry views), and of those with definite mismatches whether every mismatch is a
        # first-Met one and the first such codon (the text below). A clean sequence of a run that was not aborted
        # is never fixable; in an aborted run the sequences before the failed gene count as written (not_written,
        # fixable only with a definite mismatch) and those after it as not reached (fixable).
        problems = []
        any_start, all_start, first_codon = False, True, None
        kinds, flags = self.kinds, self.flags
        nk = len(kinds)
        ff = self._first_fail() if aborted else None
        for k in range(n):
            kind = chr(kinds[k]) if k < nk else "n"
            f = flags[k] if k < nk else 0
            if not aborted or kind == "e":
                if kind == "o":
                    continue
                fix = (kind == "m" and not f & F_FRAMESHIFT) or (kind == "f" and f & (F_UNRELATED | F_REAL_MM)) \
                    or kind in ("e", "n")
                has_mm = f & F_REAL_MM
            elif ff is None or k > ff:
                fix, has_mm = True, False
            else:
                fix = has_mm = f & F_REAL_MM
            if fix:
                problems.append(k)
                if has_mm:
                    if not any_start:
                        first_codon = bytes(self.codons[3 * k:3 * k + 3]).decode("utf-8", "replace")
                    any_start = True
                    all_start = all_start and bool(f & F_ALL_START)
        if not problems:
            return
        pset = set(problems)
        tested = (problems + [k for k in range(n) if k not in pset])[:HINT_SAMPLE]
        data = {}
        want = set(tested)
        if not want:
            return
        for k, (pid, row) in enumerate(rows()):
            if k in want:
                if len(row) - row.count("-") > HINT_MAX_RESIDUES:
                    want.discard(k)
                else:
                    data[k] = (pid, row, dna.get(pairing[k]))
            if len(data) == len(want):
                break
        tested = [k for k in tested if k in data]
        if not tested:
            return
        found = search_codes(tested, problems, data, self.table, self.fit_memo)
        if found is None:
            return
        tied, standard, top, fits_now = found
        nt = len(tested)
        c = tied[0]
        total, table = n, self.table
        name = APP_CODE_NAMES[c]

        do = _do
        if nt == total:
            who = f"all {total} sequences match" if top == total else f"{top} of {total} sequences match"
            text = f"Under genetic code {c} ({name}), {who} their protein; under code {table}, {do(fits_now)}."
        else:
            which = (f"a sample of the {len(problems)} with problems" if len(problems) >= nt
                     else f"all {len(problems)} with problems, then others")
            text = (f"COAL tested {nt} of the {total} sequences ({which}): under "
                    f"genetic code {c} ({name}), {top} of them match their protein; under code {table}, "
                    f"{do(fits_now)}.")
        if any_start and all_start:
            text = (f"The only mismatches are at the first codon of {len(pset)} sequence{'s' if len(pset) != 1 else ''}"
                    f" (for example {first_codon} where the protein starts with M): an alternative start codon, "
                    f"which code {c} accepts. " + text)
        if len(tied) > 1:
            others = " and ".join([", ".join(str(x) for x in tied[1:-1]), str(tied[-1])] if len(tied) > 2
                                  else [str(tied[1])])
            closer = code_distance(table, c)[0] < code_distance(table, tied[1])[0]
            text += (f" Code{'s' if len(tied) > 2 else ''} {others} fit{'' if len(tied) > 2 else 's'} "
                     + ("as well, but code 1 is the standard code, used by most nuclear genes." if standard else
                        f"as well, but code {c} differs least from code {table}." if closer else
                        "these sequences too: choose the code of the organism and organelle the genes come from."))
        # On the command line the summary line names the option; the report keeps the app's text ("Re-run with code
        # 2?", the app's button wording).
        self.gc_hint_cli = text + f" Re-run with -codontable {c}."
        text += f" Re-run with code {c}?"
        self.gc_hint = {"code": "genetic_code", "text": text, "suggested_codontable": c,
                        "checked": nt, "fit": top, "fit_current": fits_now}

    def entry_view(self, k: int, aborted: bool | None = None) -> dict:
        """Sequence k as the report shows it (views() gives them all, in order, reading the entries once)."""
        e = next((x for i, x in self._entries() if i == k), None)
        return self._view(k, e, aborted)

    def _view(self, k: int, e: dict | None, aborted: bool | None = None) -> dict:
        """Sequence k as the report shows it, from its kept entry (None for a clean or unreached sequence). In a run
        that stopped at an error (no output), as the app's: genes before the first failed gene are 'not_checked'
        with a 'not_written' message (and their own messages), genes after it 'not_checked' with 'not_reached';
        failed genes stay failed."""
        if aborted is None:
            aborted = self.output is None
        if e is None:
            pid = self.ids[k] if self.ids is not None else None
            st = "ok" if k < len(self.kinds) and chr(self.kinds[k]) == "o" else "not_checked"
            e = {"id": pid, "index": k, "dna_id": self._dna_of(k) if self.pairing else None,
                 "pairing": self._pairing_of(k) if self.pairing else None, "status": st, "messages": []}
        if not aborted or (k < len(self.kinds) and chr(self.kinds[k]) == "e"):
            return e
        ff = self._first_fail()
        own = [m for m in e["messages"] if m.get("code") != "suspect_pairing"]
        if ff is None:
            msgs = []
        elif k < ff:
            msgs = [{"code": "not_written", "failed_id": self.first_excluded["id"]}] + own
        else:
            msgs = [{"code": "not_reached", "failed_id": self.first_excluded["id"]}]
        return dict(e, status="not_checked", messages=msgs)

    def views(self, aborted: bool | None = None):
        """Every sequence's entry view, in index order, reading the kept entries once (the report's "sequences")."""
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        it = self._entries()
        nxt = next(it, None)
        for k in range(n):
            e = None
            while nxt is not None and nxt[0] < k:  # (an entry beyond n is never written)
                nxt = next(it, None)
            if nxt is not None and nxt[0] == k:
                e = nxt[1]
                nxt = next(it, None)
            yield self._view(k, e, aborted)

    def _first_fail(self):
        # Cached: entry_view calls this for every sequence, and a full scan each time made a run that stopped on a
        # late gene error quadratic in the number of sequences. kinds only changes through add_sequence/_set_kind
        # (append or overwrite) and restore (new object), which the key covers together with _set_kind's reset.
        key = (id(self.kinds), len(self.kinds))
        if self._ff_cache is None or self._ff_cache[0] != key:
            i = self.kinds.find(b"e")
            self._ff_cache = (key, None if i < 0 else i)
        return self._ff_cache[1]

    # -------------------------------------------------- views
    def _all_entries(self) -> list[dict]:
        return list(self.views())

    def counts(self) -> dict:
        """The summary's counters, from what COAL found for each sequence it processed (its kind: ok, a warning,
        failed) plus `not_checked` for the sequences it never reached; they add up to the total. In a run that stopped
        at an error these are the numbers the `coal: summary:` lines give ("N sequences were not checked"), not the
        entry views' statuses, which also mark the sequences checked before the error `not_checked` (`not_written`)
        because no output holds them."""
        c = dict.fromkeys(STATUSES, 0)
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        st = {"o": "ok", "m": "warning", "u": "warning", "x": "warning", "f": "failed", "e": "failed",
              "n": "not_checked"}
        for ch in self.kinds:
            c[st[chr(ch)]] += 1
        c["not_checked"] += max(0, n - len(self.kinds))
        return c

    def kind_counts(self) -> dict:
        out: dict[str, int] = {}
        for ch in self.kinds:
            out[chr(ch)] = out.get(chr(ch), 0) + 1
        return out

    def _id_hints(self) -> list[dict]:
        """The app's id_hints: protein IDs that are not DNA IDs, when pairing by order (or not at all)."""
        if self.partial or self.pairing not in ("order", None) or self.ids is None or self.dna_file_ids is None:
            return []
        dids = self.dna_file_ids
        dset = set(dids)
        missing = [p for p in self.ids if p not in dset]
        if not missing:
            return []
        c = self.counts()
        total = len(self.ids)
        clean = self.exit_code == 0 and total and c["ok"] == total
        if self.pairing is None:
            verdict = "COAL stopped before pairing them."
        elif clean:
            verdict = ("Sequences were paired by their order in the files instead, and every pair matched codon "
                       "for codon, so the pairing looks right; renaming makes it safe.")
        else:
            verdict = ("Sequences were therefore paired by their order in the files, which is probably wrong "
                       "for the sequences that failed or have warnings.")
        dnorm: dict[str, str] = {}
        for d in dids:
            dnorm.setdefault(strip_version(d), d)
        fixable = [p for p in missing if strip_version(p) in dnorm]
        if fixable and len(fixable) == len(missing):
            ex = fixable[0]
            return [{"code": "id_mismatch",
                     "text": (f"The protein IDs don't exactly match the DNA IDs: they differ only by a version "
                              f"suffix (for example protein '{ex}' vs DNA '{dnorm[strip_version(ex)]}'). {verdict} "
                              f"Rename the sequences so the IDs match exactly, for example by removing the '.1'-style "
                              f"version suffixes from the protein alignment ('{ex}' -> '{strip_version(ex)}'), then "
                              f"run again.")}]
        return [{"code": "id_mismatch",
                 "text": (f"{len(missing)} protein ID(s) were not found in the DNA file (for example "
                          f"'{missing[0]}'). {verdict} Check that both files describe the same genes, in the same "
                          f"order.")}]

    def hints(self) -> list[dict]:
        return self.pair_hints + self._id_hints() + ([self.gc_hint] if self.gc_hint else [])

    def completed(self) -> bool:
        return (self.exit_code == 0 or (self.exit_code == 2 and self.keep_going)) and not self.errors \
            and self.output is not None

    def command(self) -> str:
        """The command that reproduces the output (options that change it, then -o), as the app writes it."""
        import shlex
        parts = ["coal", self.pep, *self.nucs]
        if self.fmt != "fasta":
            parts += ["-output", self.fmt]
        if self.table != 1:
            parts += ["-codontable", str(self.table)]
        for flag, on in (("-nogap", self.nogap), ("-nomismatch", self.nomismatch), ("-blockonly", self.blockonly),
                         ("--strict", self.strict), ("--keep-going", self.keep_going)):
            if on:
                parts.append(flag)
        if self.out:
            parts += ["-o", self.out]
        return " ".join(shlex.quote(str(p)) for p in parts)

    def summary_lines(self) -> list[str]:
        """Plain-language lines for the end of a run with problems ([] if none).
        They only sum up what stderr already reported: nothing about mismatches with -nomismatch (as PAL2NAL,
        -nomismatch prints no mismatch warnings) and nothing about pairing by order when no ID matched (silent as
        PAL2NAL); the report holds both."""
        kc = self.kind_counts()
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        lines = []
        hint = self.gc_hint_cli if self.gc_hint else None  # the command-line wording
        k = 0 if self.nomismatch else kc.get("m", 0)
        if k:
            fate = "removed by -nomismatch" if self.nomismatch else "kept in the alignment"
            s = (f"{_n(k)} of {_n(n)} sequences " + ("has" if k == 1 else "have")
                 + f" codons that do not match their amino acid ({fate}).")
            if hint:
                s += " " + hint
            lines.append(s)
        k = 0 if self.nomismatch else kc.get("u", 0)
        if k:
            lines.append(f"{_n(k)} of {_n(n)} sequences " + ("has" if k == 1 else "have") + " codons with unknown or "
                         "ambiguous bases (e.g. N) that could not be checked: missing sequence data, not a "
                         "genetic-code problem.")
        # (a gene's kind is the one recorded at add_sequence: the end-of-run break events and pattern marks never
        # change it)
        fl = kc.get("f", 0)
        if fl:
            first = self.first_flagged
            lines.append(f"{_seqs(fl)} of {_n(n)} {'was' if fl == 1 else 'were'} written to the output but "
                         f"{'is' if fl == 1 else 'are'} probably wrong (first: {first['id']}: "
                         f"{first['messages'][0]['text']}) Remove {'it' if fl == 1 else 'them'} before "
                         f"analysis, or use --strict to refuse such genes.")
        ex = kc.get("e", 0)
        c = self.counts()
        if ex:
            first = self.first_excluded
            why = next((m["text"] for m in first["messages"] if m["code"] not in ("frameshift", "missing_dna")
                        and not (m["code"] == "inconsistent" and "column" in m)), first["messages"][-1]["text"])
            where = next((m["text"] for m in first["messages"] if "column" in m), "")
            gc = " " + hint if hint and (self.nomismatch or not kc.get("m")) else ""
            if self.keep_going:
                left = (" and " + ("was" if ex == 1 else "were") + " left out") if not self.errors else ""
                where = " " + where if where else ""
                lines.append(f"{_seqs(ex)} of {_n(n)} could not be converted{left} (first: {first['id']}: "
                             f"{why}.{where}){gc}")
            else:
                # Where the gene breaks, without repeating the error just printed, and what to do (--keep-going).
                if where:
                    lines.append(f"{first['id']}: {where}{gc}")
                elif gc:
                    lines.append(gc.lstrip())
                left = c["not_checked"]
                s = f"The run stopped at {first['id']}"
                s += (f": {_seqs(left)} {'was' if left == 1 else 'were'} not checked and no output was written."
                      if left else ", so no output was written.")
                if n > 1:
                    s += " Add --keep-going to write the other sequences and leave this one out."
                lines.append(s)
        # A run-level error (e.g. the counts differ) is printed once, as the error, not repeated here.
        if self.pairing == "normalised":
            k = sum(1 for m in (self.seq_pairing or []) if m == "normalised")
            j = next((i for i, m in enumerate(self.seq_pairing or []) if m == "normalised"), None)
            eg = f", e.g. protein {self.ids[j]} = DNA {self.dna_ids[j]}" if j is not None and self.ids else ""
            lines.append(f"{_n(k)} of {_n(n)} IDs matched only after normalising (version suffix, case or database "
                         f"fields){eg}; COAL paired by those IDs where PAL2NAL's file order would mispair them.")
        unk = sum(1 for f in self.flags if f & F_UNKNOWN_AA)
        if unk:
            lines.append(f"{_n(unk)} of {_n(n)} sequences contain non-standard amino acids (U, O, B, Z or J) whose "
                         f"codons were not checked.")
        return lines

    def methods_text(self, rep: dict | None = None) -> str:
        """A ready-to-paste methods sentence: the app's rules (coal_app.core.methods_text) on this report."""
        import shlex
        rep = rep or self.to_dict(methods=False)
        if not self.completed():
            return ""
        s = rep["summary"]
        n = self.table
        total, ok = s.get("total", 0), s.get("ok", 0)
        phrase = f"{CODE_NAMES[n]} genetic code"
        parts = [f"Codon alignments were built with COAL v{__version__} ({REPO_URL}), which places each coding DNA "
                 f"sequence onto the protein multiple sequence alignment codon by codon (the method of PAL2NAL; "
                 f"Suyama et al. 2006), using the {phrase} (NCBI translation table {n})."]
        # The kind of every sequence's entry view: the run completed (its output was written), so a view is its
        # entry, whose kind is the one recorded in `kinds`; the message codes of the "other" entries are in the
        # flags (the C++ engine counts them the same way), so the entries are not read again here
        kinds: dict[str, int] = {}
        names = {"o": "ok", "m": "mismatch", "u": "unchecked", "x": "other", "f": "flagged", "e": "excluded",
                 "n": "not_checked"}
        for ch in self.kinds:
            g = names[chr(ch)]
            kinds[g] = kinds.get(g, 0) + 1
        codes: set = set()  # the message codes of the "other" entries
        for ch, f in zip(self.kinds, self.flags):
            if ch == 120:  # ord("x")
                if f & F_UNKNOWN_AA:
                    codes.add("unknown_aa")
                if f & F_PAIRING:
                    codes.add("suspect_pairing")
                if f & F_OTHER_CODE:
                    codes.add("other")
        nomm = self.nomismatch
        clauses = []

        def were(k):
            return "was" if k == 1 else "were"
        if kinds.get("mismatch"):
            k = kinds["mismatch"]
            clauses.append(f"{k} contained codons that did not match the aligned amino acid"
                           + ("" if nomm else " (kept in the alignment)"))
        if kinds.get("unchecked"):
            k = kinds["unchecked"]
            clauses.append(f"{k} contained codons with unknown or ambiguous bases (e.g. N) that could not be checked"
                           + ("" if nomm else " (kept in the alignment)"))
        other = kinds.get("other", 0)
        if other:
            why = []
            if "unknown_aa" in codes:
                why.append("contained non-standard amino acids (such as U, B or Z) whose codons were not checked")
            if "suspect_pairing" in codes:
                why.append("were paired with their DNA only by their order in the files")
            if codes - {"unknown_aa", "suspect_pairing"}:
                why.append("had other warnings (see the report)")
            k = other
            clauses.append(f"{k} " + (why[0] if len(why) == 1 else "had warnings: they " + " or ".join(why)))
        if kinds.get("flagged"):
            k = kinds["flagged"]
            clauses.append(f"{k} {were(k)} written to the output but flagged as probably wrong (most likely the wrong "
                           f"DNA, a shifted position or the wrong genetic code); remove "
                           f"{'it' if k == 1 else 'them'} before analysis")
        if kinds.get("excluded"):
            k = kinds["excluded"]
            clauses.append(f"{k} could not be aligned and {were(k)} excluded")
        counts = f"Of {total} sequences, {ok} matched their protein exactly"
        if clauses:
            counts += ("; " + "; ".join(clauses[:-1]) + "; and " + clauses[-1] if len(clauses) > 1
                       else ("; " if kinds.get("flagged") else ", and ") + clauses[0])
        parts.append(counts + ".")
        removed = s.get("columns_removed")
        if nomm and (kinds.get("mismatch") or kinds.get("unchecked")):
            parts.append(f"The {removed} codon columns holding such codons were removed (-nomismatch)." if removed
                         else "The alignment columns holding such codons were removed (-nomismatch).")
        if self.nogap:
            parts.append("Codon columns containing gaps or stop codons were removed.")
        cmd = rep.get("command") or ""
        for item in list(rep.get("inputs") or []) + [rep.get("output") or {}]:
            pth = item.get("path") if isinstance(item, dict) else None
            if pth:
                # The quoted full path is replaced first: on Windows every full path is quoted (backslashes), and
                # replacing only the folder would leave 'protein.aln.fasta' in quotes
                if shlex.quote(pth) != pth:
                    cmd = cmd.replace(shlex.quote(pth), shlex.quote(os.path.basename(pth)))
                cmd = cmd.replace(pth, os.path.basename(pth))
        if cmd:
            parts.append(f"Command: {cmd}")
        text = " ".join(parts)
        if total and 2 * ok < total:
            # A "ready-to-paste" methods text for a run where most genes have problems would read as if the
            # alignment were fine, so it starts with a note to fix the warnings first (the app shows the report's
            # methods_text as it is).
            text = (f"{NOT_READY} only {ok} of {total} sequences matched their protein exactly, so this alignment is "
                    f"probably not ready to use (see the hints and the messages for each sequence). {AS_IT_STANDS} "
                    + text)
        return text

    def to_dict(self, methods: bool = True) -> dict:
        return self._dict_with(self._all_entries(), methods)

    def _dict_with(self, entries, methods: bool = True) -> dict:
        """to_dict with the given entry views (a list, or in _write_json a callable giving a fresh iterator)."""
        import shlex
        views = entries if callable(entries) else (lambda: entries)
        c = self.counts()
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        # "path" is the argument as given (`-` for standard input); "sha256" is of the data read
        def _hash(k, p):
            if k in self.input_hashes:
                return self.input_hashes[k]
            return None if p == "-" else sha256_file(p)  # `-` never names a file in the current folder
        inputs = [{"role": "protein", "path": self.pep, "sha256": _hash(0, self.pep),
                   "sequences": len(self.ids) if self.ids is not None else None}]
        for i, p in enumerate(self.nucs):
            inputs.append({"role": "dna", "path": p, "sha256": _hash(i + 1, p),
                           "sequences": self.dna_counts[i] if self.dna_counts else None})
        p2n = ["pal2nal.pl", self.pep, *self.nucs, "-output", self.fmt]
        if self.table != 1:
            p2n += ["-codontable", str(self.table)]
        for flag, on in (("-nogap", self.nogap), ("-nomismatch", self.nomismatch), ("-blockonly", self.blockonly)):
            if on:
                p2n.append(flag)
        summary = {"total": n, "ok": c["ok"], "warning": c["warning"], "failed": c["failed"],
                   "pairing": self.pairing}
        if c["not_checked"]:
            summary["not_checked"] = c["not_checked"]
        if self.nomismatch and self.output is not None:  # as the app: counted for a run that wrote its output
            summary["columns_removed"] = len(self.cols)
            summary["columns_removed_ambiguous"] = len(self.cols - self.real_cols)
        d = {
            "schema": SCHEMA,
            "coal_version": __version__,
            "source": "coal",
            "engine": ENGINE,
            "command": self.command(),
            "command_line": " ".join(shlex.quote(a) for a in ["coal", *self.argv]),
            "pal2nal_equivalent": " ".join(shlex.quote(str(a)) for a in p2n),
            "options": {"output_format": self.fmt, "codontable": self.table, "nogap": self.nogap,
                        "nomismatch": self.nomismatch, "keep_going": self.keep_going, "strict": self.strict,
                        "blockonly": self.blockonly},
            "inputs": inputs,
            "output": self.output,
            "summary": summary,
            "hints": self.hints(),
            "errors": list(self.errors),
            "sequences": views(),
            "runtime_seconds": round(time.perf_counter() - self.t0, 3),
            "exit_code": self.exit_code,
        }
        d["summary_text"] = self.summary_lines()
        if methods:
            d["methods_text"] = self.methods_text(d)
        if self.error:
            d["error"] = self.error
        if self.resumed_from is not None:
            d["resumed_from"] = self.resumed_from
        return d

    def write_json(self, fh) -> None:
        """Write the report (json.dump(self.to_dict(), fh, indent=1), byte for byte) one sequence entry at a time, so a
        run of a million sequences never holds a million entry dicts at once."""
        import json
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        d = self._dict_with(self.views)
        head = {k: v for k, v in d.items() if k != "sequences"}
        keys = list(d)
        at = keys.index("sequences")
        # json.dumps(indent=1) of the keys before "sequences", without its closing "\n}"; then the array, each entry
        # as json.dumps writes it two levels deep; then the keys after it, without their opening "{"
        fh.write(json.dumps({k: head[k] for k in keys[:at]}, indent=1)[:-2])
        fh.write(',\n "sequences": ')
        if n == 0:
            fh.write("[]")
        else:
            fh.write("[")
            for k, view in enumerate(self.views()):
                fh.write(",\n  " if k else "\n  ")
                fh.write(json.dumps(view, indent=1).replace("\n", "\n  "))
            fh.write("\n ]")
        fh.write("," + json.dumps({k: head[k] for k in keys[at + 1:]}, indent=1)[1:])

    def write(self, path: str):
        """Write the report to `path`, streamed one sequence at a time (write_json_with: safe target handling)."""
        def writer(fh):
            self.write_json(fh)
            fh.write("\n")
        write_json_with(path, writer)


def _same_file(a, b) -> bool:
    """Whether two paths name the same file: by identity when both exist, else by resolved path."""
    try:
        return os.path.samefile(a, b)
    except OSError:
        pass
    try:
        return os.path.normcase(os.path.realpath(a)) == os.path.normcase(os.path.realpath(b))
    except OSError:
        return False


def report_path_conflict(report, inputs, out) -> str | None:
    """Why `--report REPORT` cannot be written: it names an input file or the -o output (the report replaces what is
    at its path, which would destroy them). A usage error (exit 2), checked before anything runs; also for the
    reports of a batch's job lines. None when the path is fine. The C++ engine mirrors this (main.cpp)."""
    if not report:
        return None
    for p in inputs:
        if str(p) != "-" and _same_file(report, p):
            return f"the report path {report} is an input file"
    if out and _same_file(report, out):
        return f"the report path {report} is the output file"
    return None


def report_folder_missing(report) -> str | None:
    """Why `--report REPORT` could not be written at the end: its folder does not exist (or is not a folder). Checked
    before anything runs (a usage error, exit 2), so a full run never ends by losing its report. None when the path
    is fine, or already exists (a FIFO, a device, an existing file). The C++ engine mirrors this (sysio.cpp)."""
    if not report or os.path.exists(report):
        return None
    folder = os.path.dirname(os.path.abspath(report))
    if os.path.isdir(folder):
        return None
    if os.path.exists(folder):
        return f"the report path {report}: {folder} is not a folder"
    return f"the report path {report}: folder {folder} does not exist"


def _test_pause_in_report() -> None:
    """Test hook: COAL_TEST_HOOKS=1 and COAL_TEST_PAUSE_IN_REPORT=1 stop here while the report's temp file exists (so
    a test can send SIGTERM at that point). Short sleeps, as pipeline._test_pause."""
    if os.environ.get("COAL_TEST_HOOKS") == "1" and os.environ.get("COAL_TEST_PAUSE_IN_REPORT") == "1":
        print("coal: test: paused in report", file=sys.stderr, flush=True)
        deadline = time.monotonic() + 3600
        while time.monotonic() < deadline:
            time.sleep(0.05)


def write_json(path: str, data: dict) -> None:
    """Write `data` as JSON to `path` (write_json_with)."""
    import json
    text = json.dumps(data, indent=1) + "\n"
    write_json_with(path, lambda fh: fh.write(text))


def write_json_with(path: str, writer) -> None:
    """Write a JSON file as the alignment is written (pipeline._run): a temp file renamed onto the real target,
    through a symlink, so a failure never leaves a half-written file; a path that exists and is not a regular file
    (a FIFO, /dev/stdout, /dev/null) or a file with other hard links is written in place instead of being replaced
    by a regular file. Never leaves the temp file behind. `writer(fh)` writes the text."""
    # An existing non-regular target is written through as it is, before any path resolution (as -o does): the real
    # path of /dev/fd/1 on a pipe is "pipe:[N]", which no folder holds
    try:
        if not os.path.isdir(path) and not stat.S_ISREG(os.stat(path).st_mode):
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                writer(fh)
            return
    except OSError:
        pass
    real = os.path.realpath(path)
    try:
        st = os.stat(real)
    except OSError:
        st = None
    if st is not None and (not os.path.isfile(real) or st.st_nlink > 1):
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            writer(fh)
        return
    tmp = f"{real}.tmp{os.getpid()}"
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
            writer(fh)
            _test_pause_in_report()
        if st is not None:
            import shutil
            shutil.copymode(real, tmp)  # keep the existing file's mode
        os.replace(tmp, real)
    except BaseException:  # an OSError, and SIGTERM / SIGHUP (SystemExit from cli._on_sigterm) or Ctrl-C too
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
