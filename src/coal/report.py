"""Run report (schema "coal.report/1", read by the COAL GUI) and the
plain-language end-of-run summary with genetic-code hints.

The report is built while the pipeline runs (one small entry per sequence) and
written at the end, also when the run failed, so a person can see why.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import time

from . import __version__
from .tables import SUPPORTED

SCHEMA = "coal.report/1"
HINT_SAMPLE = 200          # problem sequences re-tested with every genetic code
HINT_THRESHOLD = 0.20      # hint only when more than this share mismatches or fails

TABLE_NAMES = {
    1: "standard", 2: "vertebrate mitochondrial", 3: "yeast mitochondrial",
    4: "mold, protozoan and coelenterate mitochondrial", 5: "invertebrate mitochondrial",
    6: "ciliate, dasycladacean and Hexamita nuclear", 9: "echinoderm and flatworm mitochondrial",
    10: "euplotid nuclear", 11: "bacterial, archaeal and plant plastid", 12: "alternative yeast nuclear",
    13: "ascidian mitochondrial", 14: "alternative flatworm mitochondrial", 15: "Blepharisma nuclear",
    16: "chlorophycean mitochondrial", 21: "trematode mitochondrial",
    22: "Scenedesmus obliquus mitochondrial", 23: "Thraustochytrium mitochondrial",
}

_MISMATCH = re.compile(r"pepAlnPos (\d+): (.) does not correspond to (\S+)")
_UNKNOWN = re.compile(r"pepAlnPos (\d+): (.) unknown AA type")

# Plain-English reasons for errors (code -> text); {detail} is the raw message.
_REASONS = {
    "no_match": "The protein could not be found in its DNA sequence in any reading frame, even allowing for "
                "a few wrong codons. Check that this is the right DNA for this protein and the genetic code "
                "(-codontable).",
    "invalid_chars": "The protein alignment contains characters that are not amino acids ({detail}).",
    "empty_protein": "This protein row contains only gaps.",
    "count_mismatch": "The protein alignment and the DNA file(s) contain different numbers of sequences "
                      "({detail}); every protein needs exactly one DNA sequence.",
    "no_sequences": "No sequences were found in the input files.",
    "row_lengths": "The protein alignment rows have different lengths ({detail}), so it is not an alignment.",
    "input_format": "An input file could not be read: {detail}.",
    "file_error": "A file could not be opened or written: {detail}.",
    "all_failed": "None of the sequences could be converted, so no output was written.",
}


REPO_URL = "https://github.com/chloe-newt/COAL"
FORMAT_NAMES = {"fasta": "FASTA", "clustal": "CLUSTAL", "paml": "PAML", "codon": "codon-table (PAL2NAL '-output codon')"}


def plain_reason(code: str, message: str) -> str:
    detail = message.split(": ", 1)[1] if code in ("invalid_chars", "no_match") and ": " in message else message
    m = re.search(r"\(([^()]*)\)\s*$", message)
    if code in ("count_mismatch", "row_lengths") and m:
        detail = m.group(1)
    if code == "invalid_chars":
        detail = message.rsplit("characters ", 1)[-1]
    return _REASONS.get(code, "{detail}").format(detail=detail)


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


class RunReport:
    def __init__(self, argv: list[str], pep: str, nucs: list[str], out: str | None, fmt: str, table: int,
                 nogap: bool, nomismatch: bool, keep_going: bool = False, keep_entries: bool = True):
        self.t0 = time.perf_counter()
        # keep_entries=False: only counts (enough for the summary lines); used when no
        # report was asked for, so big runs don't pay for 100,000 per-sequence entries.
        self.keep_entries = keep_entries
        self.n = {"ok": 0, "warning": 0, "failed": 0, "mismatch": 0, "unknown_aa": 0}
        self.failed_entries: list[dict] = []
        self.pep, self.nucs, self.out, self.fmt, self.table = pep, list(nucs), out, fmt, table
        self.nogap, self.nomismatch, self.keep_going = nogap, nomismatch, keep_going
        self.command = shlex.join(["coal", *argv])
        p2n = ["pal2nal.pl", pep, *nucs, "-output", fmt]
        if table != 1:
            p2n += ["-codontable", str(table)]
        if nogap:
            p2n.append("-nogap")
        if nomismatch:
            p2n.append("-nomismatch")
        self.pal2nal_equivalent = shlex.join(p2n)
        self.n_protein: int | None = None
        self.dna_counts: list[int] | None = None
        self.pairing = None
        self.pair_example = None
        self.dna_ids: list[str] | None = None     # DNA ID paired with each protein, input order
        self.dna_index: list[int] | None = None   # DNA record number paired with each protein
        self.seq_pairing: list[str] | None = None  # per-sequence pairing mode (when it varies)
        self.entries: list[dict] = []
        self.problems: list[int] = []             # input positions of mismatching/failed sequences
        self.error: dict | None = None
        self.output: dict | None = None
        self.hints: list[dict] = []
        self.exit_code = 0
        self.extra_summary: list[str] = []        # extra "coal: summary:" lines from the pipeline

    # -------------------------------------------------- filled in by the pipeline
    def set_inputs(self, n_protein: int, dna_counts: list[int]):
        self.n_protein, self.dna_counts = n_protein, dna_counts

    def set_pairing(self, mode: str, dna_ids, dna_index, per_seq: list[str] | None = None, example=None):
        """dna_ids: callable k -> DNA ID or a list; dna_index: DNA record number per protein;
        example: (protein ID, DNA ID) of one normalised match."""
        self.pairing, self.dna_ids, self.dna_index, self.seq_pairing = mode, dna_ids, dna_index, per_seq
        self.pair_example = example

    def add_sequence(self, pid: str, warnings: list[str], failure: tuple[str, str] | None = None,
                     k: int | None = None):
        """Record one sequence; k = its position in the input (default: next)."""
        k = sum(self.n[x] for x in ("ok", "warning", "failed")) if k is None else k
        if not warnings and failure is None and not self.keep_entries:
            self.n["ok"] += 1
            return
        msgs, mismatch = [], False
        for w in warnings:
            m = _MISMATCH.search(w)
            if m:
                mismatch = True
                msgs.append({"code": "mismatch", "column": int(m.group(1)),
                             "text": f"Column {m.group(1)}: {m.group(2)} does not match codon {m.group(3)}"})
                continue
            m = _UNKNOWN.search(w)
            if m:
                msgs.append({"code": "unknown_aa", "column": int(m.group(1)),
                             "text": f"Column {m.group(1)}: {m.group(2)} is not a standard amino acid; "
                                     f"treated as X (any codon)"})
                continue
            msgs.append({"code": "warning", "text": w})
        if failure:
            code, message = failure
            msgs.append({"code": code, "text": plain_reason(code, message)})
            status = "failed"
        else:
            status = "warning" if msgs else "ok"
        if (failure or mismatch) and len(self.problems) < HINT_SAMPLE:
            self.problems.append(k)
        self.n[status] += 1
        self.n["mismatch"] += mismatch
        self.n["unknown_aa"] += any(m["code"] == "unknown_aa" for m in msgs)
        dna_id = None
        if self.dna_ids is not None:
            dna_id = self.dna_ids(k) if callable(self.dna_ids) else self.dna_ids[k]
        entry = {"id": pid, "index": k, "dna_id": dna_id,
                 "pairing": (self.seq_pairing[k] if self.seq_pairing else self.pairing),
                 "status": status, "messages": msgs}
        if failure:
            self.failed_entries.append(entry)
        if self.keep_entries:
            self.entries.append(entry)

    def fail(self, code: str, message: str):
        self.error = {"code": code, "text": plain_reason(code, message), "detail": message}

    def set_output(self, path: str | None, sha: str | None, n: int):
        self.output = {"path": path, "format": self.fmt, "sha256": sha, "sequences": n}

    # -------------------------------------------------- end of run
    def counts(self):
        return {k: self.n[k] for k in ("ok", "warning", "failed")}

    def compute_hints(self):
        """Genetic-code hint: re-test up to 200 problem sequences with every code."""
        c = self.counts()
        done = sum(c.values())
        mism = self.n["mismatch"] + c["failed"]
        if not done or mism / done <= HINT_THRESHOLD or not self.problems or self.dna_index is None:
            return
        fits = genetic_code_fits(self.pep, self.nucs, self.problems, self.dna_index)
        if not fits:
            return
        n = len(self.problems)
        cur = fits.get(self.table, 0)
        best = max(fits.items(), key=lambda kv: (kv[1], -kv[0]))
        t, k = best
        if t == self.table or k < 0.8 * n or k - cur < 0.2 * n:
            return
        name = TABLE_NAMES[t]
        quant = "They all fit" if k == n else f"{_n(k)} of the {_n(n)} checked fit"
        self.hints.append({
            "code": "genetic_code",
            "text": f"Most sequences fit genetic code {t} ({name}). Try -codontable {t}.",
            "suggested_codontable": t,
            "summary": f"{quant} genetic code {t} ({name}): try -codontable {t}.",
            "checked": n, "fit": k, "fit_current": cur,
        })

    def summary_lines(self) -> list[str]:
        """Plain-language lines for the end of a run with problems ([] if none)."""
        c = self.counts()
        done = sum(c.values())
        total = self.n_protein if self.n_protein is not None else done
        lines = list(self.extra_summary)
        mism = self.n["mismatch"]
        hint = next((h for h in self.hints if h["code"] == "genetic_code"), None)
        if mism:
            s = f"{_n(mism)} of {_n(total)} sequences " + ("has" if mism == 1 else "have") + " codon mismatches."
            if hint:
                s += " " + hint["summary"]
            lines.append(s)
        failed = self.failed_entries
        if failed:
            first = failed[0]
            left_out = " and " + ("was" if len(failed) == 1 else "were") + " left out" \
                if self.keep_going and not self.error else ""
            s = (f"{_seqs(len(failed))} of {_n(total)} could not be converted{left_out} (first: {first['id']}: "
                 f"{first['messages'][-1]['text']})")
            if hint and not mism:
                s += " " + hint["summary"]
            lines.append(s)
        if self.error and not failed:
            lines.append(f"The run failed: {self.error['text']}")
        if self.pairing == "order" and done:
            lines.append("Protein and DNA IDs don't match, so sequences were paired by file order. "
                         "Check that the files list the same genes in the same order.")
        elif self.pairing == "normalised" and done:
            k = sum(1 for m in (self.seq_pairing or []) if m == "normalised")
            ex = f", e.g. protein {self.pair_example[0]} = DNA {self.pair_example[1]}" if self.pair_example else ""
            lines.append(f"{_n(k)} of {_n(total)} IDs matched only after normalising (version suffix, "
                         f"case or pipe fields){ex}.")
        unk = self.n["unknown_aa"]
        if unk:
            lines.append(f"{_n(unk)} of {_n(total)} sequences contain non-standard amino acids (treated as X).")
        if self.error and self.n_protein and done < self.n_protein and failed:
            left = self.n_protein - done
            lines.append(f"The run stopped there: {_seqs(left)} {'was' if left == 1 else 'were'} not processed "
                         f"and no output was written.")
        return lines

    def methods_text(self) -> str:
        """Ready-to-paste methods sentence describing this run."""
        c = self.counts()
        n = self.output["sequences"] if self.output else sum(c.values())
        t = self.table
        code = "standard genetic code" if t == 1 else f"{TABLE_NAMES[t]} genetic code"
        opts = []
        if self.nogap:
            opts.append("removing codon columns containing gaps or stop codons (-nogap)")
        if self.nomismatch:
            opts.append("removing columns with codons that do not match their amino acid (-nomismatch)")
        opt_txt = (", " + " and ".join(opts)) if opts else ""
        skipped = ""
        if c["failed"] and self.keep_going and not self.error:
            skipped = (f" ({_seqs(c['failed'])} that could not be matched to "
                       f"{'its' if c['failed'] == 1 else 'their'} DNA {'was' if c['failed'] == 1 else 'were'} excluded)")
        verified = "outputs are verified byte-identical to PAL2NAL v14 on the COAL test suite"
        if t in (10, 11, 13):
            verified += f" (COAL follows the current NCBI definition of table {t}, a documented difference)"
        return (f"Codon alignments were generated with COAL v{__version__} ({REPO_URL}) from {_n(n)} aligned "
                f"protein sequence{'' if n == 1 else 's'} and their coding sequences{skipped}, using the {code} "
                f"(NCBI table {t}) and {FORMAT_NAMES[self.fmt]} output{opt_txt}; {verified}.")

    def to_dict(self) -> dict:
        c = self.counts()
        done = sum(c.values())
        total = self.n_protein if self.n_protein is not None else done
        inputs = [{"role": "protein", "path": self.pep, "sha256": sha256_file(self.pep),
                   "sequences": self.n_protein}]
        for i, p in enumerate(self.nucs):
            inputs.append({"role": "dna", "path": p, "sha256": sha256_file(p),
                           "sequences": self.dna_counts[i] if self.dna_counts else None})
        d = {
            "schema": SCHEMA,
            "coal_version": __version__,
            "command": self.command,
            "pal2nal_equivalent": self.pal2nal_equivalent,
            "options": {"output_format": self.fmt, "codontable": self.table, "nogap": self.nogap,
                        "nomismatch": self.nomismatch, "keep_going": self.keep_going},
            "inputs": inputs,
            "output": self.output or {"path": self.out, "format": self.fmt, "sha256": None, "sequences": 0},
            "summary": {"total": total, **c, "pairing": self.pairing},
            "hints": [{k: v for k, v in h.items() if k != "summary"} for h in self.hints],
            "sequences": self.entries,
            "runtime_seconds": round(time.perf_counter() - self.t0, 3),
            "exit_code": self.exit_code,
        }
        d["summary"]["not_processed"] = max(0, total - done) if total else 0
        d["summary_text"] = self.summary_lines()
        d["methods_text"] = self.methods_text()
        if self.error:
            d["error"] = self.error
        return d

    def write(self, path: str):
        tmp = f"{path}.tmp{os.getpid()}"
        with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(self.to_dict(), fh, indent=1)
            fh.write("\n")
        os.replace(tmp, path)


def genetic_code_fits(pep, nucs, positions, dna_index) -> dict[int, int]:
    """For each supported code: how many of the sequences at `positions` thread
    with no codon mismatch at all."""
    from .core import CoalError, thread_one
    from .pipeline import _protein_rows
    from .seqio import DnaIndex

    want = set(positions)
    rows = {}
    for k, (pid, row) in enumerate(_protein_rows(pep)()):
        if k in want:
            rows[k] = (pid, row)
        if len(rows) == len(want):
            break
    dna = DnaIndex(nucs)
    try:
        seqs = [(rows[k][0], rows[k][1], dna.get(dna_index[k])) for k in positions if k in rows]
    finally:
        dna.close()
    fits = {}
    for t in SUPPORTED:
        n = 0
        for pid, row, s in seqs:
            try:
                if not thread_one(pid, row, s, t).mismatch_cols:
                    n += 1
            except CoalError:
                pass
        fits[t] = n
    return fits
