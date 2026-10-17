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
HINT_ALL_MAX = 2000        # up to this many sequences, the genetic-code hint re-tests every sequence
HINT_SAMPLE = 200          # above it: this many problem sequences + this many others

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
                         f"to refuse such genes.")}
    m = _GUESS.match(rest)
    if m:
        nchk, nres, short = m.groups()
        why = (f"Only {nchk} of this gene's {nres} codons have no unknown bases (N), and no "
               if short is None else f"The gene has only {short} codons, and no ")
        return {"code": "suspect_position",
                "text": (f"{why}10-residue stretch of the protein was found in the DNA, so COAL cannot tell where "
                         f"the gene starts in the longer DNA. Its row was taken from the start of the DNA (as "
                         f"PAL2NAL does) and may be shifted: do not trust it. Run COAL with --strict to refuse "
                         f"such genes.")}
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


def break_event(row: str, dna: str, table: int) -> dict | None:
    """The first event of the app viewer's diagnose(): where a failed gene stops matching (frameshift, missing
    DNA, or the column from which the DNA no longer spells the protein). Same algorithm, so the GUI says the same."""
    gc = get_code(table)
    code, starts = gc.code, gc.starts
    dna = re.sub(r"[^A-Za-z]", "", dna).upper().replace("U", "T")
    res = [(i + 1, a.upper()) for i, a in enumerate(row) if a.isalpha()]
    if not res or len(dna) < 3:
        return None
    want = min(len(res), 12)
    allowed = max(0, want // 4)
    best = (0, -1)
    for off in range(0, len(dna) - 2):
        score = 0
        for k in range(want):
            if _ok(res[k][1], dna[off + 3 * k: off + 3 * k + 3], code, starts, k == 0):
                score += 1
            elif k + 1 - score > allowed:
                break
        if score > best[1]:
            best = (off, score)
            if score == want:
                break
    if best[1] < max(min(3, len(res)), want - allowed):
        return {"code": "inconsistent", "column": res[0][0],
                "text": "The DNA does not spell the start of this protein in any reading frame."}
    pos, k, flags = best[0], 0, {}
    while k < len(res):
        col, aa = res[k]
        codon = dna[pos:pos + 3]
        if len(codon) < 3:
            return {"code": "missing_dna", "column": col, "text": f"The DNA ends before the protein does (at column "
                                                                  f"{col})."}
        if _ok(aa, codon, code, starts, k == 0):
            flags[col] = ""
            pos += 3
            k += 1
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
        for s in (-1, 1, -2, 2):
            n = min(10, len(res) - k - 1)
            if n <= 0:
                break
            hits = sum(_ok(a, dna[pos + 3 + s + 3 * j: pos + 3 + s + 3 * j + 3], code, starts, False)
                       for j, (_, a) in enumerate(res[k + 1:k + 1 + n]))
            if hits >= max(1, int(0.8 * n)):
                shifted = s
                break
        if shifted is not None:
            what = {-1: "one base is missing", -2: "two bases are missing",
                    1: "there is one extra base", 2: "there are two extra bases"}[shifted]
            return {"code": "frameshift", "column": col, "shift": shifted,
                    "text": f"Frameshift at column {col}: {what} in the DNA, so from here on the codons are read out "
                            f"of step."}
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
    def finish(self, rows, dna, pairing):
        """Locate where failed genes break and compute the genetic-code hint, while the inputs are open.
        rows: callable giving (id, prepared row) in input order; dna: the DnaIndex; pairing: DNA index per protein."""
        self.finished = True
        try:
            self._finish(rows, dna, pairing)
        except Exception:  # noqa: BLE001 - extras must never hide the real result
            pass

    def _finish(self, rows, dna, pairing):
        from .core import CoalError, MISMATCH_WARNING, thread_one
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        kinds = self.kinds
        failed = [k for k, e in self.entries.items() if e["status"] == "failed"
                  and not any(isinstance(m.get("column"), int) for m in e["messages"])]
        problems = sorted(k for k, ch in enumerate(kinds) if chr(ch) in "mfe"
                          and not (chr(ch) == "e" and self.entries[k]["messages"][0]["code"] in ("invalid_protein",
                                                                                                  "empty")))
        if not failed and not problems:
            return
        if n <= HINT_ALL_MAX:
            test = list(range(n))
        else:
            others = [k for k, ch in enumerate(kinds) if chr(ch) == "o"][:HINT_SAMPLE]
            test = sorted(set(problems[:HINT_SAMPLE]) | set(others))
        want = set(test) | set(failed)
        data = {}
        for k, (pid, row) in enumerate(rows()):
            if k in want:
                data[k] = (pid, row, dna.get(pairing[k]))
            if len(data) == len(want):
                break
        for k in failed[:HINT_SAMPLE]:
            ev = break_event(data[k][1], data[k][2], self.table)
            if ev:
                self.entries[k]["messages"].insert(0, ev)
        if not problems:
            return

        def fits(t):
            good = set()
            for k in test:
                pid, row, s = data[k]
                try:
                    it = thread_one(pid, row, s, t)
                except CoalError:
                    continue
                if any(MISMATCH_WARNING in w and not is_ambiguous_codon(w.rsplit(" ", 1)[-1]) for w in it.warnings):
                    continue
                if any(": protein and DNA are probably unrelated" in x or ": position of the gene" in x
                       for x in it.notices):
                    continue
                good.add(k)
            return good
        cur = fits(self.table)
        res = {t: fits(t) for t in SUPPORTED if t != self.table}
        if not res:
            return
        top = max(len(v) for v in res.values())
        if top <= len(cur):
            return
        tested_problems = [k for k in test if k in set(problems)]
        tied = sorted(t for t, v in res.items() if len(v) == top)
        best = tied[0]
        if 11 in tied and len(tied) > 1 and self._start_only(tested_problems):
            best = 11
        fixed = len(res[best] & set(tested_problems))
        if 2 * fixed < len(tested_problems):
            return
        k, c, total = top, len(cur), len(test)
        allof = len(test) == n
        if allof:
            got = f"all {_n(total)} sequences match" if k == total else f"{_n(k)} of the {_n(total)} sequences match"
            now = "none do" if c == 0 else f"{'only ' if c * 2 < total else ''}{_n(c)} {'does' if c == 1 else 'do'}"
        else:
            got = f"{_n(k)} of the {_n(total)} sequences tested (a sample of the {_n(n)}) match"
            now = "none do" if c == 0 else f"{_n(c)} {'does' if c == 1 else 'do'}"
        name = CODE_NAMES[best]
        others = [t for t in tied if t != best]
        also = ""
        if others:
            ts = [str(t) for t in others]
            also = (f" Code{'s' if len(ts) > 1 else ''} {', '.join(ts[:-1]) + ' and ' + ts[-1] if len(ts) > 1 else ts[0]}"
                    f" fit{'' if len(ts) > 1 else 's'} these sequences equally well.")
        text = (f"With genetic code {best} ({name}), {got} their protein codon for codon; with code {self.table} "
                f"({CODE_NAMES[self.table]}), the one used, {now}. Re-run with code {best} (-codontable {best})."
                + also)
        self.gc_hint = {"code": "genetic_code", "text": text, "suggested_codontable": best,
                        "checked": total, "fit": k, "fit_current": c}

    def _start_only(self, problems) -> bool:
        """True if every definite mismatch of the problem genes is at the gene's first residue (a start codon)."""
        seen = False
        for k in problems:
            e = self.entries.get(k)
            if e is None:
                return False
            mm = [m for m in e["messages"] if m.get("code") == "mismatch" and not is_ambiguous_codon(m.get("codon"))]
            if not mm or e["status"] == "failed":
                return False
            firstcol = min(m["column"] for m in mm)
            if any(m["column"] != firstcol for m in mm) or mm[0].get("residue", "").upper() != "M":
                return False
            seen = True
        return seen

    # -------------------------------------------------- views
    def _all_entries(self) -> list[dict]:
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        out = []
        for k in range(n):
            e = self.entries.get(k)
            if e is None:
                pid = self.ids[k] if self.ids is not None else None
                st = "ok" if k < len(self.kinds) and chr(self.kinds[k]) == "o" else "not_checked"
                e = {"id": pid, "index": k, "dna_id": self._dna_of(k) if self.pairing else None,
                     "pairing": self._pairing_of(k) if self.pairing else None, "status": st, "messages": []}
            out.append(e)
        return out

    def counts(self) -> dict:
        n = len(self.ids) if self.ids is not None else len(self.kinds)
        c = dict.fromkeys(STATUSES, 0)
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
                              f"Remove the '.1'-style suffixes so the names match, or use a COAL version that "
                              f"matches IDs after normalising them.")}]
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
        hint = self.gc_hint
        k = 0 if self.nomismatch else kc.get("m", 0)
        if k:
            fate = "removed by -nomismatch" if self.nomismatch else "kept in the alignment"
            s = (f"{_n(k)} of {_n(n)} sequences " + ("has" if k == 1 else "have")
                 + f" codons that do not match their amino acid ({fate}).")
            if hint:
                s += " " + hint["text"]
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
        if ex:
            first = ent[ex[0]]
            left = (" and " + ("was" if len(ex) == 1 else "were") + " left out") if self.keep_going and \
                not self.errors else ""
            why = next((m["text"] for m in first["messages"] if m["code"] not in ("frameshift", "missing_dna")
                        and not (m["code"] == "inconsistent" and "column" in m)), first["messages"][-1]["text"])
            where = next((" " + m["text"] for m in first["messages"] if "column" in m), "")
            s = f"{_seqs(len(ex))} of {_n(n)} could not be converted{left} (first: {first['id']}: {why}.{where})"
            if hint and (self.nomismatch or not kc.get("m")):
                s += " " + hint["text"]
            lines.append(s)
        if self.errors and not ex:
            lines.append(f"The run failed: {self.errors[0]['text']}")
        c = self.counts()
        if (self.errors or ex) and self.exit_code not in (0, 2) and c["not_checked"]:
            left = c["not_checked"]
            lines.append(f"The run stopped there: {_seqs(left)} {'was' if left == 1 else 'were'} not checked and no "
                         f"output was written.")
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
                cmd = cmd.replace(pth, os.path.basename(pth))
        if cmd:
            parts.append(f"Command: {cmd}")
        return " ".join(parts)

    def to_dict(self, methods: bool = True) -> dict:
        c = self.counts()
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
            "sequences": self._all_entries(),
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
