"""Run report (schema "coal.report/1", read by the COAL GUI), the plain-language end-of-run summary
("coal: summary:" lines) and the genetic-code hint.

The report is built while the pipeline runs and written at the end, also when the run failed, so a person can
see why.

CHANGED(MRG-2): this module comes from the teammate's branch (--report/--json, summary lines, hints, methods
text). Its per-sequence messages, statuses, hints and methods text now follow the rules of the reviewed app
report (src/coal_app/core.py: build_fallback_report, mismatch_probe, id_hints, genetic_code_hint, methods_text,
viewer.locate_breaks), computed here from the engine's own results instead of parsed from stderr: the same
message codes and texts, DEC-1 suspect genes "failed" but written, mismatches counted only for A/C/G/T codons,
methods text counted by message code, genetic-code hint from real mismatches only and tested on all genes
(start-codon-only ties prefer code 11), sequences an aborted run never reached "not_checked" while the ones it
did check keep their status. See overnight/integration/CHANGES.md. teammate's C++ must mirror.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import time
from itertools import product

from . import __version__
from .codon import IUPAC
from .tables import SUPPORTED, get_code

SCHEMA = "coal.report/1"
REPO_URL = "https://github.com/plumat-coal/COAL"  # CHANGED(MRG-2): was github.com/chloe-newt/COAL (WEB-7)
HINT_SAMPLE = 200          # the genetic-code hint re-tests at most this many sequences, problem genes first (app)

# Names inside a sentence: common words lower case, genus names capitalised as NCBI writes them (FIN-C11)
CODE_NAMES = {
    1: "standard", 2: "vertebrate mitochondrial", 3: "yeast mitochondrial",
    4: "mold, protozoan and coelenterate mitochondrial and Mycoplasma/Spiroplasma", 5: "invertebrate mitochondrial",
    6: "ciliate, dasycladacean and Hexamita nuclear", 9: "echinoderm and flatworm mitochondrial",
    10: "euplotid nuclear", 11: "bacterial, archaeal and plant plastid", 12: "alternative yeast nuclear",
    13: "ascidian mitochondrial", 14: "alternative flatworm mitochondrial", 15: "Blepharisma nuclear",
    16: "chlorophycean mitochondrial", 21: "trematode mitochondrial",
    22: "Scenedesmus obliquus mitochondrial", 23: "Thraustochytrium mitochondrial",
}
TABLE_NAMES = CODE_NAMES  # the teammate's name for it

SUSPECT_CODES = ("suspect_unrelated", "suspect_position")
# CHANGED(ENG-UX-9): the start and the end of the note before the methods text of a run with mostly problem genes
NOT_READY = "Fix the warnings first:"
AS_IT_STANDS = "Methods text for this run as it stands:"
# one character per sequence in RunReport.kinds
KIND_CHAR = {"ok": "o", "mismatch": "m", "unchecked": "u", "other": "x", "flagged": "f", "excluded": "e",
             "not_checked": "n"}
STATUSES = ("ok", "warning", "failed", "not_checked")
AS_X_NAMES = {"U": "selenocysteine", "O": "pyrrolysine", "B": "D or N", "Z": "E or Q", "J": "I or L"}

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
                         f"amino acids, so COAL did not check its codon (any codon accepted).")}
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


# CHANGED(MRG-2b): the app's FIN-C4/C5/C6/C7/C9 rules (src/coal_app/core.py real_mismatches, _fixable, _fits_now,
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


SHIFTS = (-1, 1, -2, 2, -4, 4, -5, 5)  # ENG-UX-19: +-4 / +-5 = a frameshift together with a lost or extra codon
SHIFT_WORDS = {-1: "one base is missing", -2: "two bases are missing", 1: "there is one extra base",
               2: "there are two extra bases", -4: "four bases are missing", -5: "five bases are missing",
               4: "there are four extra bases", 5: "there are five extra bases"}


def break_event(row: str, dna: str, table: int) -> dict | None:
    """The first event of the app viewer's diagnose(): where a failed gene stops matching (frameshift, missing
    DNA, or the column from which the DNA no longer spells the protein). Same algorithm, so the GUI says the same.
    CHANGED(ENG-UX-19): on real data it named wrong columns (a camelid gene with partial codons at columns 21 and 24:
    "Frameshift at column 94"; a gene missing its first 2 nt, with a duplicated stretch: "The DNA ends before the
    protein does (at column 3344)"). Now: the earliest credible start (not the best-scoring one, which can be a
    duplicated stretch further on); a shift must also hold for the next 20 residues and may include a whole lost or
    extra codon (+-4, +-5); "at column N" only when the old frame holds up to N and the new one from the next
    residue, else "near column N"; "the DNA ends" only after a walk where most codons matched (else no location).
    teammate's C++ mirrors it (cpp/src/report.cpp). See overnight/integration/CHANGES.md"""
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
        if score >= need:  # ENG-UX-19: the earliest credible start
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
            if k and 5 * matched >= 4 * k:  # ENG-UX-19: only after a credible walk
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
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def _n(x: int) -> str:
    return f"{x:,}"


def _seqs(n: int) -> str:
    return f"{_n(n)} sequence" + ("" if n == 1 else "s")


def strip_version(seq_id: str) -> str:
    return re.sub(r"\.\d+$", "", seq_id)


class RunReport:
    """Collects one run's results. keep_entries=False (no --report/--json): ok sequences are only counted."""

    def __init__(self, argv, pep: str, nucs, out: str | None, fmt: str, table: int, nogap: bool, nomismatch: bool,
                 keep_going: bool = False, keep_entries: bool = True, strict: bool = False, blockonly: bool = False):
        self.t0 = time.perf_counter()
        self.keep_entries = keep_entries
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
        self.entries: dict[int, dict] = {}         # non-ok entries (all entries are rebuilt from these)
        self.cols: set[int] = set()                # -nomismatch: columns with a mismatch warning
        self.real_cols: set[int] = set()           # ... of which with a definite (A/C/G/T) mismatch
        self.error: dict | None = None
        self.errors: list[dict] = []
        self.output: dict | None = None
        self.gc_hint: dict | None = None
        self.gc_hint_cli: str | None = None        # the hint as the command line says it (ENG-UX-7)
        self.exit_code = 0
        self.resumed_from: int | None = None
        self._dict = None
        self.finished = False

    # -------------------------------------------------- filled in by the pipeline
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

    def add_sequence(self, k: int, pid: str, warnings, notices, failure=None):
        """Record sequence number k. warnings: thread_one's warnings as printed (or, with -nomismatch, as they
        would be); notices: its notices; failure: (code, message) if it failed. Returns the entry if not ok."""
        self._dict = None
        while len(self.kinds) < k:
            self.kinds.append(ord("n"))
        msgs = []
        if failure:
            status, msgs = "failed", [failure_message(pid, *failure)]
        else:
            if self.nomismatch:
                mm = [warning_message(w, removed=True) for w in warnings if _MISMATCH.match(w)]
                for m in mm:
                    self.cols.add(m["column"])
                    if not is_ambiguous_codon(m["codon"]):
                        self.real_cols.add(m["column"])
                msgs = mm + [warning_message(w) for w in warnings if not _MISMATCH.match(w)]
            else:
                msgs = [warning_message(w) for w in warnings]
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
            self._set_kind(k, "o")
            return None
        entry = {"id": pid, "index": k, "dna_id": self._dna_of(k), "pairing": self._pairing_of(k),
                 "status": status, "messages": msgs}
        self.entries[k] = entry
        self._set_kind(k, KIND_CHAR[gene_kind(entry)])
        return entry

    def _set_kind(self, k, ch):
        if k < len(self.kinds):
            self.kinds[k] = ord(ch)
        else:
            self.kinds.append(ord(ch))

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

    def restore(self, state, nonok: list[dict], completed: int):
        self._dict = None
        state = state or {}
        self.kinds = bytearray(state.get("kinds", "o" * completed).encode("ascii"))[:completed]
        self.cols, self.real_cols = set(state.get("cols", ())), set(state.get("real_cols", ()))
        self.entries = {e["index"]: e for e in nonok if e["index"] < completed}
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
            pass

    def _finish(self, rows, dna, pairing, aborted):
        from .core import CoalError, thread_one
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        failed = {k for k, e in self.entries.items() if e["status"] == "failed"
                  and not any(isinstance(m.get("column"), int) for m in e["messages"])}
        marks = {k for k, e in self.entries.items() if real_mismatches(e)}
        need = set(sorted(failed)[:HINT_SAMPLE]) | marks
        rowof = {}
        if need:
            for k, (pid, row) in enumerate(rows()):
                if k in need:
                    rowof[k] = row
                    if k in failed:
                        ev = break_event(row, dna.get(pairing[k]), self.table)
                        if ev:
                            self.entries[k]["messages"].insert(0, ev)
                    if k in marks:
                        annotate_mismatch_pattern(row, self.entries[k])
        # CHANGED(MRG-2b): the app's FIN-A1/FIN-C4/FIN-C5 hint (src/coal_app/core.py genetic_code_hint), computed here
        view = [self.entry_view(k, aborted) for k in range(n)]
        problems = [k for k, e in enumerate(view) if fixable(e)]
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
                data[k] = (pid, row, dna.get(pairing[k]))
            if len(data) == len(want):
                break
        def bad_under(c):
            bad = set()
            for k in tested:
                pid, row, s = data[k]
                try:
                    it = thread_one(pid, row, s, c)
                except CoalError:
                    bad.add(k)  # (the app drops the whole code when one gene fails: --keep-going exits 2)
                    continue
                if any(_MISMATCH.match(w) and not is_ambiguous_codon(w.rsplit(" ", 1)[-1]) for w in it.warnings) or \
                        any(notice_message(pid, x)["code"] in SUSPECT_CODES for x in it.notices):
                    bad.add(k)
            return bad
        # fits under the code used: the same test (not the statuses; in a stopped run the genes never reached
        # would otherwise count as fixed by any code)
        bad_now = bad_under(self.table)
        fits_now = sum(1 for k in tested if k not in bad_now)
        results = []
        for c in SUPPORTED:
            if c == self.table:
                continue
            bad = bad_under(c)
            fit = sum(1 for k in tested if k not in bad)
            fixed = sum(1 for k in problems if k in bad_now and k not in bad)
            if fixed and fit > fits_now:
                results.append((fit, c))
        if not results:
            return
        top = max(f for f, _ in results)
        tied = sorted((c for f, c in results if f == top), key=lambda c: (code_distance(self.table, c)[0], c))
        # CHANGED(PRE-4): was always the tied code that differs least from the one used (one_broken_gene run with
        # code 2: "code 4 ... Codes 1, 10, 11, 15, 16 and 6 fit as well"); now the standard code 1 when it is among
        # them (by far the most common right answer for nuclear genes). See overnight/integration/CHANGES.md
        standard = 1 in tied and len(tied) > 1
        if standard:
            tied = [1] + [x for x in tied if x != 1]
        c = tied[0]
        nt, total, table = len(tested), n, self.table
        name = APP_CODE_NAMES[c]

        def do(k):
            return "none do" if k == 0 else ("only 1 does" if k == 1 else f"{k} do")
        if nt == total:
            who = f"all {total} sequences match" if top == total else f"{top} of {total} sequences match"
            text = f"Under genetic code {c} ({name}), {who} their protein; under code {table}, {do(fits_now)}."
        else:
            which = (f"a sample of the {len(problems)} with problems" if len(problems) >= nt
                     else f"all {len(problems)} with problems, then others")
            text = (f"COAL tested {nt} of the {total} sequences ({which}): under "
                    f"genetic code {c} ({name}), {top} of them match their protein; under code {table}, "
                    f"{do(fits_now)}.")
        starts = [m for k in problems for m in real_mismatches(view[k])]
        if starts and all(m.get("start_codon") for m in starts):
            ex = starts[0]
            text = (f"The only mismatches are at the first codon of {len(pset)} sequence{'s' if len(pset) != 1 else ''}"
                    f" (for example {ex.get('codon')} where the protein starts with M): an alternative start codon, "
                    f"which code {c} accepts. " + text)
        if len(tied) > 1:
            others = " and ".join([", ".join(str(x) for x in tied[1:-1]), str(tied[-1])] if len(tied) > 2
                                  else [str(tied[1])])
            closer = code_distance(table, c)[0] < code_distance(table, tied[1])[0]
            text += (f" Code{'s' if len(tied) > 2 else ''} {others} fit{'' if len(tied) > 2 else 's'} "
                     + ("as well, but code 1 is the standard code, used by most nuclear genes." if standard else
                        f"as well, but code {c} differs least from code {table}." if closer else
                        "these sequences too: choose the code of the organism and organelle the genes come from."))
        # CHANGED(ENG-UX-7): the summary line said "Re-run with code 2?" (the app's button wording); on the command line
        # it now names the option. The report keeps the app's text (the app shows its own button; contract
        # unchanged). See overnight/integration/CHANGES.md
        self.gc_hint_cli = text + f" Re-run with -codontable {c}."
        text += f" Re-run with code {c}?"
        self.gc_hint = {"code": "genetic_code", "text": text, "suggested_codontable": c,
                        "checked": nt, "fit": top, "fit_current": fits_now}

    def entry_view(self, k: int, aborted: bool | None = None) -> dict:
        """Sequence k as the report shows it. In a run that stopped at an error (no output), as the app's FIN-C7:
        genes before the first failed gene are 'not_checked' with a 'not_written' message (and their own messages),
        genes after it 'not_checked' with 'not_reached'; failed genes stay failed."""
        if aborted is None:
            aborted = self.output is None
        e = self.entries.get(k)
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
            msgs = [{"code": "not_written", "failed_id": self.entries[ff]["id"]}] + own
        else:
            msgs = [{"code": "not_reached", "failed_id": self.entries[ff]["id"]}]
        return dict(e, status="not_checked", messages=msgs)

    def _first_fail(self):
        return next((k for k, ch in enumerate(self.kinds) if chr(ch) == "e"), None)

    # -------------------------------------------------- views
    def _all_entries(self) -> list[dict]:
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        return [self.entry_view(k) for k in range(n)]

    def counts(self, entries=None) -> dict:
        c = dict.fromkeys(STATUSES, 0)
        if entries is not None:
            for e in entries:
                c[e["status"]] += 1
            return c
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
        CHANGED(MRG-2): they only sum up what stderr already reported: nothing about mismatches with -nomismatch
        (ENG-i2: as PAL2NAL, -nomismatch prints no mismatch warnings) and nothing about pairing by order when no
        ID matched (ENG-d: silent as PAL2NAL); the report holds both."""
        kc = self.kind_counts()
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        ent = self.entries
        lines = []
        hint = self.gc_hint_cli if self.gc_hint else None  # CHANGED(ENG-UX-7): the command-line wording
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
        fl = sorted(i for i, e in ent.items() if gene_kind(e) == "flagged")
        if fl:
            first = ent[fl[0]]
            lines.append(f"{_seqs(len(fl))} of {_n(n)} {'was' if len(fl) == 1 else 'were'} written to the output but "
                         f"{'is' if len(fl) == 1 else 'are'} probably wrong (first: {first['id']}: "
                         f"{first['messages'][0]['text']}) Remove {'it' if len(fl) == 1 else 'them'} before "
                         f"analysis, or use --strict to refuse such genes.")
        ex = sorted(i for i, e in ent.items() if gene_kind(e) == "excluded")
        c = self.counts()
        if ex:
            first = ent[ex[0]]
            why = next((m["text"] for m in first["messages"] if m["code"] not in ("frameshift", "missing_dna")
                        and not (m["code"] == "inconsistent" and "column" in m)), first["messages"][-1]["text"])
            where = next((m["text"] for m in first["messages"] if "column" in m), "")
            gc = " " + hint if hint and (self.nomismatch or not kc.get("m")) else ""
            if self.keep_going:
                left = (" and " + ("was" if len(ex) == 1 else "were") + " left out") if not self.errors else ""
                where = " " + where if where else ""
                lines.append(f"{_seqs(len(ex))} of {_n(n)} could not be converted{left} (first: {first['id']}: "
                             f"{why}.{where}){gc}")
            else:
                # CHANGED(ENG-UX-2, ENG-UX-5): was "1 sequence of 20 could not be converted (first: <id>: <the error
                # just printed>. <where>)" and "The run stopped there: ... no output was written." (no way out named);
                # now where the gene breaks, without repeating the error, and what to do (--keep-going).
                # See overnight/integration/CHANGES.md
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
        # CHANGED(ENG-UX-5): a run-level error (e.g. the counts differ) was printed twice: "coal: error: X" and
        # "coal: summary: The run failed: X" (+ "The run stopped there: ... not checked"); now once, as the error.
        if self.pairing == "normalised":
            k = sum(1 for m in (self.seq_pairing or []) if m == "normalised")
            j = next((i for i, m in enumerate(self.seq_pairing or []) if m == "normalised"), None)
            eg = f", e.g. protein {self.ids[j]} = DNA {self.dna_ids[j]}" if j is not None and self.ids else ""
            lines.append(f"{_n(k)} of {_n(n)} IDs matched only after normalising (version suffix, case or database "
                         f"fields){eg}; COAL paired by those IDs where PAL2NAL's file order would mispair them.")
        unk = sum(1 for e in ent.values() if any(m.get("code") == "unknown_aa" for m in e["messages"]))
        if unk:
            lines.append(f"{_n(unk)} of {_n(n)} sequences contain non-standard amino acids (U, O, B, Z or J) whose "
                         f"codons were not checked.")
        return lines

    def methods_text(self, rep: dict | None = None) -> str:
        """A ready-to-paste methods sentence: the app's rules (src/coal_app/core.py methods_text) on this report."""
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
        kinds: dict[str, list] = {}
        for sq in rep["sequences"]:
            kinds.setdefault(gene_kind(sq), []).append(sq)
        nomm = self.nomismatch
        clauses = []

        def were(k):
            return "was" if k == 1 else "were"
        if kinds.get("mismatch"):
            k = len(kinds["mismatch"])
            clauses.append(f"{k} contained codons that did not match the aligned amino acid"
                           + ("" if nomm else " (kept in the alignment)"))
        if kinds.get("unchecked"):
            k = len(kinds["unchecked"])
            clauses.append(f"{k} contained codons with unknown or ambiguous bases (e.g. N) that could not be checked"
                           + ("" if nomm else " (kept in the alignment)"))
        other = kinds.get("other") or []
        if other:
            codes = {m.get("code") for sq in other for m in sq.get("messages") or []}
            why = []
            if "unknown_aa" in codes:
                why.append("contained non-standard amino acids (such as U, B or Z) whose codons were not checked")
            if "suspect_pairing" in codes:
                why.append("were paired with their DNA only by their order in the files")
            if codes - {"unknown_aa", "suspect_pairing"}:
                why.append("had other warnings (see the report)")
            k = len(other)
            clauses.append(f"{k} " + (why[0] if len(why) == 1 else "had warnings: they " + " or ".join(why)))
        if kinds.get("flagged"):
            k = len(kinds["flagged"])
            clauses.append(f"{k} {were(k)} written to the output but flagged as probably wrong (most likely the wrong "
                           f"DNA, a shifted position or the wrong genetic code); remove "
                           f"{'it' if k == 1 else 'them'} before analysis")
        if kinds.get("excluded"):
            k = len(kinds["excluded"])
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
                # CHANGED(SHIP-2): the quoted full path is replaced first; on Windows every full path is quoted
                # (backslashes), and replacing only the folder left 'protein.aln.fasta' in quotes
                if shlex.quote(pth) != pth:
                    cmd = cmd.replace(shlex.quote(pth), shlex.quote(os.path.basename(pth)))
                cmd = cmd.replace(pth, os.path.basename(pth))
        if cmd:
            parts.append(f"Command: {cmd}")
        text = " ".join(parts)
        if total and 2 * ok < total:
            # CHANGED(ENG-UX-9): new. A "ready-to-paste" methods text for a run where most genes have problems read as
            # if the alignment were fine; it now starts with a note to fix the warnings first (the app shows the
            # report's methods_text as it is). See overnight/integration/CHANGES.md
            text = (f"{NOT_READY} only {ok} of {total} sequences matched their protein exactly, so this alignment is "
                    f"probably not ready to use (see the hints and the messages for each sequence). {AS_IT_STANDS} "
                    + text)
        return text

    def to_dict(self, methods: bool = True) -> dict:
        entries = self._all_entries()
        c = self.counts(entries)
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        inputs = [{"role": "protein", "path": self.pep, "sha256": sha256_file(self.pep),
                   "sequences": len(self.ids) if self.ids is not None else None}]
        for i, p in enumerate(self.nucs):
            inputs.append({"role": "dna", "path": p, "sha256": sha256_file(p),
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
            "sequences": entries,
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

    def write(self, path: str):
        """Write atomically (temp file + rename); never leaves the temp file behind."""
        tmp = f"{path}.tmp{os.getpid()}"
        try:
            with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
                json.dump(self.to_dict(), fh, indent=1)
                fh.write("\n")
            os.replace(tmp, path)
        except OSError:
            try:
                os.unlink(tmp)
            except OSError:
                pass
            raise
