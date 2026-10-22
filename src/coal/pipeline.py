"""Streaming, optionally parallel conversion for large inputs.

Memory use does not grow with the number of sequences beyond their IDs:
  pass 1  read the protein alignment once: IDs, row length, frameshift columns
          index the DNA file(s): IDs and byte offsets only
  pass 2  read protein rows again one at a time, fetch the matching DNA from
          disk, thread codons (in worker processes if jobs > 1, in order), and
          write each output row as soon as it is ready.
With -nogap / -nomismatch / -blockonly, which depend on every row, pass 2 writes rows to a
temporary file and a final pass removes the flagged columns. The clustal and codon formats
interleave all sequences block by block, so they also take that path: the final rows go to a
fixed-width temporary file, which is then read back a group of blocks at a time, or, when one block of
all sequences is larger than BLOCK_MEMORY, one block a chunk of sequences at a time; the rows and
blocks held in memory stay bounded by BLOCK_MEMORY, not by the number or length of the sequences (the IDs
and padded names are kept, as in every format).

Output goes to a temporary file first and is only moved/copied to the real
destination when everything succeeded, so an error never leaves partial output
(PAL2NAL also prints nothing on error).
"""
from __future__ import annotations

import errno
import json
import os
import re
import time
from array import array
import shutil
import signal
import sys
import tempfile
from collections import deque
from itertools import islice

from . import checkpoint as _ck
from .checkpoint import Checkpoint, durable_size, open_append
from .core import (
    MISMATCH_WARNING, CoalError, apply_slices, check_table, column_widths, expand, keep_slices, nogap_bad_triplets,
    pair_ids, prepare_row, row_length_message, thread_one,
)
from .seqio import DnaIndex, InputFormatError, clustal_mask, is_clustal, iter_fasta, open_error_message, read_clustal

WRAP = 60
_DIGIT = re.compile(r"[0-9]")
# PAL2NAL v14's four output formats
FORMATS = ("fasta", "clustal", "paml", "codon")
BLOCK_FORMATS = ("clustal", "codon")  # interleaved blocks of 60 nt
BLOCK_MEMORY = 16_000_000  # bytes of rows held at once while writing clustal / codon blocks
_ALNPOS = re.compile(r" pepAlnPos (\d+):")
BATCH_NT = int(os.environ.get("COAL_BATCH_NT", 2_000_000))  # target DNA per work batch (keeps IPC overhead low)
# Starting worker processes costs ~0.1-0.2 s (more on Windows), about the time one
# core needs for a few million nucleotides. So each extra worker must have at
# least this much DNA to do, otherwise COAL runs in a single process.
MIN_NT_PER_JOB = 4_000_000
_WINDOWS = os.name == "nt"  # see the broken-pipe handling in _run
PROGRESS_MIN_SEQS = int(os.environ.get("COAL_PROGRESS_MIN_SEQS", 10_000))  # --progress: inputs larger than this
PROGRESS_EVERY = float(os.environ.get("COAL_PROGRESS_EVERY", 2.0))        # seconds between progress reports


class Progress:
    """--progress: percent done + estimated time left on stderr: one updating line on a
    terminal, plain lines otherwise. Never touches stdout."""

    def __init__(self, total: int, stream=None):
        self.total, self.stream = total, stream or sys.stderr
        self.t0 = self.last = time.monotonic()
        self.tty = bool(getattr(self.stream, "isatty", lambda: False)())
        self.printed = False

    def update(self, done: int, final: bool = False):
        now = time.monotonic()
        if not final and now - self.last < PROGRESS_EVERY:
            return
        self.last = now
        pct = 100 * done // self.total if self.total else 100
        el = now - self.t0
        if final:
            tail = f"done in {el:.1f} s"
        elif done:
            left = el / done * (self.total - done)
            tail = f"about {left:.0f} s left" if left < 120 else f"about {left / 60:.0f} min left"
        else:
            tail = "estimating time left"
        msg = f"coal: progress: {pct}% ({done:,} of {self.total:,} sequences), {tail}"
        self.stream.write("\r" + msg.ljust(78) + ("\n" if final else "") if self.tty else msg + "\n")
        self.stream.flush()
        self.printed = True


def _test_hooks() -> bool:
    """The test hooks below act only when COAL_TEST_HOOKS=1 is set as well, so a stray COAL_TEST_PAUSE_* variable
    can never make a real run sleep."""
    return os.environ.get("COAL_TEST_HOOKS") == "1"


def _test_pause(count: int) -> None:
    """Test hook: COAL_TEST_HOOKS=1 and COAL_TEST_PAUSE_AFTER_CHECKPOINT=k stop here after the k-th checkpoint (so a
    test can kill the run at a known point)."""
    if _test_hooks() and os.environ.get("COAL_TEST_PAUSE_AFTER_CHECKPOINT") == str(count):
        print(f"coal: test: paused after checkpoint {count}", file=sys.stderr, flush=True)
        time.sleep(3600)


def effective_jobs(requested: int, total_nt: int, min_nt_per_job: int | None = None) -> int:
    """Number of processes to actually use. requested: 0 = all cores, N = at most N."""
    # process_cpu_count (Python 3.13+) respects CPU affinity (taskset, Slurm, containers); cpu_count counts every
    # core of the machine, which over-subscribes a job limited to a few cores of a large node.
    cores = getattr(os, "process_cpu_count", os.cpu_count)() or 1
    cap = cores if requested <= 0 else requested
    per = MIN_NT_PER_JOB if min_nt_per_job is None else min_nt_per_job
    by_size = total_nt // per if per > 0 else cap
    return max(1, min(cap, by_size))


def _protein_rows(path):
    """Re-iterable source of (id, prepared row)."""
    if is_clustal(path):
        recs = [(pid, prepare_row(pid, r)) for pid, r in read_clustal(path)]
        return lambda: iter(recs)
    return lambda: ((pid, prepare_row(pid, r)) for pid, r in iter_fasta(path))


# ------------------------------------------------------------------ worker side
_W: dict = {}


def mp_child() -> bool:
    import multiprocessing
    return multiprocessing.parent_process() is not None


def _init_worker(table, widths, two_pass, nogap, strict=False, out_format="fasta", keep_going=False):
    # Workers ignore Ctrl-C (the main process reports it once, no worker tracebacks).
    if mp_child():
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    _W.update(table=table, widths=widths, wide=any(w != 3 for w in widths),
              two_pass=two_pass, nogap=nogap, strict=strict,
              record=paml_record if out_format == "paml" else fasta_record,
              keep_going=keep_going)


def _work(batch):
    """Thread one batch. Returns (text, items, mismatch columns, nogap triplets): text is finished FASTA/PAML, or
    'id<TAB>row' lines for the staged pass; items = [(id, warnings, notices, failure)] for every sequence in input
    order (failure = (code, message) of a sequence skipped by --keep-going, else None).
    Without --keep-going a CoalError is raised, with .index = position in the batch."""
    # Per-sequence items, so --keep-going can skip a sequence and the run report can keep each sequence's messages
    parts, items, mcols, bad = [], [], set(), set()
    two_pass = _W["two_pass"]
    for i, (pid, row, dna) in enumerate(batch):
        try:
            it = thread_one(pid, row, dna, _W["table"], _W["strict"])
        except CoalError as e:
            if _W["keep_going"] and e.seq_id is not None:
                items.append((pid, [], [], (e.code, str(e))))
                continue
            e.index = i
            e.done_items = items  # the batch's sequences before it (for the run report only)
            raise
        exp = expand(it, _W["widths"], _W["wide"])
        items.append((pid, it.warnings, it.notices, None))
        if two_pass:
            parts.append(f"{pid}\t{exp}\n")
            mcols |= it.mismatch_cols
            if _W["nogap"]:
                bad |= nogap_bad_triplets(exp, _W["table"])  # the selected code's stops
        else:
            parts.append(_W["record"](pid, exp))
    return "".join(parts), items, mcols, bad


# ------------------------------------------------------------------ main side

def fasta_record(rid: str, seq: str) -> str:
    """One FASTA record wrapped at 60 (an empty sequence gives an empty line, as PAL2NAL)."""
    if len(seq) <= WRAP:
        return f">{rid}\n{seq}\n"
    return f">{rid}\n" + "\n".join([seq[i:i + WRAP] for i in range(0, len(seq), WRAP)]) + "\n"


def paml_record(rid: str, seq: str) -> str:
    """One PAML (sequential PHYLIP-like) record as PAL2NAL writes it: the ID on its own line, then the row in
    lines of 60 (an empty row gives an empty line)."""
    return fasta_record(rid, seq)[1:]


def paml_header(nseq: int, alilen: int) -> str:
    return " %3d %6d\n" % (nseq, alilen)


def codon_aa_widths(ncol: int, digit_rows, last_row: str) -> list[int] | None:
    """Characters per column of the amino-acid lines of -output codon, as PAL2NAL v14 computes them: 1, except in a
    column whose largest frameshift digit is 4 or more, where the width comes from the *last* sequence's character
    in that column (int((d - 1) / 3) + 1 for a digit d, 1 for anything else). None if no row has a digit."""
    if not digit_rows:
        return None
    big = [0] * ncol
    for row in digit_rows:
        for m in _DIGIT.finditer(row):
            big[m.start()] = max(big[m.start()], int(m.group()))
    out = []
    for i in range(ncol):
        c = last_row[i]
        out.append((int(c) - 1) // 3 + 1 if big[i] >= 4 and c.isdigit() and c != "0" else 1)
    return out


def codon_aa_line(row: str, aa_widths) -> str:
    """The amino-acid line of -output codon for one row (each residue followed by '-' padding where widened)."""
    if aa_widths is None:
        return row
    return "".join(ch + "-" * (w - 1) for ch, w in zip(row, aa_widths))


def write_blocks(out, fmt: str, ids, rows_fh, alilen: int, aa_fh=None, aalen: int = 0, mask: str | None = None,
                 memory: int = BLOCK_MEMORY) -> None:
    """Write clustal or codon blocks of 60 nt. rows_fh / aa_fh: binary files holding one fixed-width row per
    sequence (alilen / aalen bytes each, in `ids` order). Reads a group of blocks for every sequence at a time.
    mask: the mark line ('#' / ' ', alilen characters) to print under each block, or None."""
    maxn = max([10] + [len(i) for i in ids])
    pad = " " * (maxn + 4)
    names = [i.ljust(maxn) + "    " for i in ids]
    if fmt == "clustal":
        out.write("CLUSTAL W multiple sequence alignment\n\n")
    nblk = -(-alilen // WRAP)
    n = len(ids)
    per_seq = WRAP + (20 if aa_fh else 0)
    if n * per_seq > memory:
        # One block of all sequences exceeds `memory` (~0.4 KB per sequence for clustal, ~0.8 KB for codon): each
        # block is written a chunk of sequences at a time, so memory stays bounded. Same bytes. The C++ engine
        # mirrors it.
        _write_blocks_chunked(out, fmt, names, pad, rows_fh, alilen, aa_fh, aalen, mask, nblk,
                              max(1, memory // (4 * per_seq)))
        return
    group = max(1, memory // max(1, n * per_seq))
    for b0 in range(0, nblk, group):
        b1 = min(nblk, b0 + group)
        segs, aas = [], []
        for k in range(n):
            rows_fh.seek(k * alilen + WRAP * b0)
            segs.append(rows_fh.read(min(WRAP * b1, alilen) - WRAP * b0).decode("ascii"))
            if aa_fh is not None:
                a0, a1 = min(20 * b0, aalen), min(20 * b1, aalen)
                aa_fh.seek(k * aalen + a0)
                aas.append(aa_fh.read(a1 - a0).decode("ascii"))
        buf = []
        for b in range(b0, b1):
            o = WRAP * (b - b0)
            for k in range(n):
                chunk = segs[k][o:o + WRAP]
                if fmt == "clustal":
                    buf.append(names[k] + chunk + "\n")
                else:
                    a = aas[k][20 * (b - b0):20 * (b - b0) + 20]
                    buf.append(pad + "   ".join(a) + "\n")
                    buf.append(names[k] + " ".join([chunk[i:i + 3] for i in range(0, len(chunk), 3)]) + "\n")
            if mask is not None:
                m = mask[WRAP * b:WRAP * b + WRAP]
                if fmt == "codon":
                    m = " ".join([m[i:i + 3] for i in range(0, len(m), 3)])
                buf.append(pad + m + "\n")
            buf.append("\n")
        out.write("".join(buf))


def _block_mark(fmt, mask, b, pad):
    m = mask[WRAP * b:WRAP * b + WRAP]
    if fmt == "codon":
        m = " ".join([m[i:i + 3] for i in range(0, len(m), 3)])
    return pad + m + "\n"


def _write_blocks_chunked(out, fmt, names, pad, rows_fh, alilen, aa_fh, aalen, mask, nblk, chunk) -> None:
    """write_blocks for many sequences: block by block, `chunk` sequences at a time."""
    n = len(names)
    for b in range(nblk):
        lo, hi = WRAP * b, min(WRAP * b + WRAP, alilen)
        a0, a1 = min(20 * b, aalen), min(20 * b + 20, aalen)
        for k0 in range(0, n, chunk):
            buf = []
            for k in range(k0, min(n, k0 + chunk)):
                rows_fh.seek(k * alilen + lo)
                seg = rows_fh.read(hi - lo).decode("ascii")
                if fmt == "clustal":
                    buf.append(names[k] + seg + "\n")
                else:
                    aa_fh.seek(k * aalen + a0)
                    buf.append(pad + "   ".join(aa_fh.read(a1 - a0).decode("ascii")) + "\n")
                    buf.append(names[k] + " ".join([seg[i:i + 3] for i in range(0, len(seg), 3)]) + "\n")
            out.write("".join(buf))
        out.write((_block_mark(fmt, mask, b, pad) if mask is not None else "") + "\n")


def _batches(rows, dna: DnaIndex, pairing, start: int = 0):
    """Batches of (id, row, DNA) of about BATCH_NT nucleotides, from sequence `start` on (--resume)."""
    batch, nt = [], 0
    for k, (pid, row) in enumerate(islice(rows, start, None), start):
        seq = dna.get(pairing[k])
        batch.append((pid, row, seq))
        nt += len(seq)
        if nt >= BATCH_NT:
            yield batch
            batch, nt = [], 0
    if batch:
        yield batch


def _results(batches, jobs: int, init_args):
    """Yield per-batch results in input order, with at most ~2*jobs batches in flight."""
    if jobs <= 1:
        _init_worker(*init_args)
        for b in batches:
            yield _work(b)
        return
    from concurrent.futures import ProcessPoolExecutor  # imported only when needed

    with ProcessPoolExecutor(max_workers=jobs, initializer=_init_worker, initargs=init_args) as ex:
        pending: deque = deque()
        for b in batches:
            pending.append(ex.submit(_work, b))
            if len(pending) >= 2 * jobs:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()


def _expand_mask(mask: str, widths, slices) -> str:
    """The '#' mark line under clustal / codon blocks, one character per output nucleotide (as PAL2NAL)."""
    full = "".join(("#" if mask[i:i + 1] == "#" else " ") * w for i, w in enumerate(widths))
    return apply_slices(full, slices)


def _write_block_output(out, fmt, ids, stage, slices, alilen, rows, aa_widths, mask, scratch, skip=()) -> None:
    """Final pass for clustal / codon: copy the kept columns of every staged row into a fixed-width binary temp
    file (and, for codon, the amino-acid lines into another), then write the blocks from those."""
    def temp():
        f = tempfile.NamedTemporaryFile("w+b", delete=False, prefix=".coal_rows_")
        scratch.append(f)
        return f
    fixed = temp()
    for line in stage:
        fixed.write(apply_slices(line.rstrip("\n").split("\t", 1)[1], slices).encode("ascii"))
    aa_fh, aalen = None, 0
    if fmt == "codon":
        aa_fh = temp()
        for k, (_pid, row) in enumerate(rows()):
            if k in skip:  # sequences skipped by --keep-going have no row
                continue
            line = codon_aa_line(row, aa_widths)
            aalen = len(line)
            aa_fh.write(line.encode("ascii", "replace"))
    fixed.flush()
    if aa_fh is not None:
        aa_fh.flush()
    write_blocks(out, fmt, ids, fixed, alilen, aa_fh, aalen, mask)


def _spool(path, spooled: list) -> str:
    """A path that can be opened several times and seeked: regular files as they are,
    anything else (pipe, /dev/stdin, <(...)) copied to a temp file first."""
    if os.path.isfile(path):
        return path
    base = re.sub(r"[^A-Za-z0-9._-]", "_", os.path.basename(str(path)))[-40:]
    with open(path, "rb") as src, tempfile.NamedTemporaryFile(
            "wb", delete=False, prefix=".coal_in_", suffix="_" + base) as dst:
        spooled.append(dst.name)
        shutil.copyfileobj(src, dst, 1 << 20)
    return dst.name


def run(pep_path, nuc_paths, *args, **kw) -> int:
    """Convert; write the codon alignment to out_path (or stdout). Arguments as for _run. Returns the number of
    sequences skipped by keep_going (0 otherwise)."""
    # Every input path is opened several times, so non-regular inputs (pipes, FIFOs) are copied to temp files first
    # and removed afterwards.
    if "fmt" in kw:  # fmt: accepted as an alias for out_format
        kw["out_format"] = kw.pop("fmt")
    spooled: list = []
    names = {}  # temp copy -> the argument it came from, for error messages
    inputs = {str(pep_path), *map(str, nuc_paths)}
    try:
        try:
            sp = _spool(pep_path, spooled)
            names[sp] = str(pep_path)
            nps = []
            for p in nuc_paths:
                nps.append(_spool(p, spooled))
                names[nps[-1]] = str(p)
            # A run whose input came through a pipe cannot be resumed (its data is gone), so it keeps no checkpoint.
            kw.setdefault("_resumable", not spooled)
            return _run(sp, nps, *args, **kw)
        except CoalError as e:  # InputFormatError included
            msg = str(e)
            for tmp_name in spooled:
                msg = msg.replace(tmp_name, names[tmp_name])
            if msg != str(e):
                if isinstance(e, InputFormatError):
                    raise InputFormatError(msg) from None
                raise CoalError(msg, e.code, e.seq_id, e.index) from None
            raise
        except OSError as e:
            # "cannot open 'x.fa': no such file" (also: permission denied, it is a folder), not Python's wording.
            fn = getattr(e, "filename", None)
            user = names.get(str(fn), str(fn)) if fn is not None else None
            if user in inputs:
                raise CoalError(open_error_message(user, e), code="file_error") from None
            raise
    finally:
        for name in spooled:
            try:
                os.unlink(name)
            except OSError:
                pass


def _output_error(out_path, e: OSError) -> CoalError:
    """The error for an output (or its partial file) that cannot be created, naming the user's -o path (not COAL's
    hidden temp file) in plain words. ENOENT in a folder that exists (e.g. /proc) is not "folder ... does not
    exist"; permission denied says what to do; with no -o (a temporary file for stdout) the temp folder is named."""
    folder = os.path.dirname(os.path.abspath(out_path)) if out_path else tempfile.gettempdir()
    if e.errno == errno.ENOENT and not os.path.isdir(folder):
        why = f"folder {folder} does not exist"
    elif e.errno == errno.ELOOP:
        why = f"folder {folder} is a symbolic link loop"
    elif e.errno == errno.ENOTDIR:
        why = f"{folder} is not a folder"
    elif e.errno in (errno.EACCES, errno.EPERM):
        why = "permission denied (choose a folder you can write to)"
    else:
        why = e.strerror or str(e)
    what = out_path if out_path else f"a temporary file in {folder}"
    return CoalError(f"cannot write {what}: {why}", code="file_error")


def _run(pep_path, nuc_paths, out_path=None, table: int = 1, nogap: bool = False,
         nomismatch: bool = False, jobs: int = 0, warn=None, min_nt_per_job: int | None = None,
         notice=None, strict: bool = False, out_format: str = "fasta", blockonly: bool = False,
         rec=None, keep_going: bool = False, progress: bool = False, resume: bool = False, note=None,
         _resumable: bool = True) -> int:
    """Convert; write the codon alignment to out_path (or stdout) in out_format (FORMATS).

    jobs: 0 = automatic (up to all cores), N = at most N processes. Small inputs
    always run in one process (see MIN_NT_PER_JOB).
    warn / notice: callables for per-codon warnings / important notices.
    blockonly: keep only the columns marked '#' in a CLUSTAL input's mark line (PAL2NAL -blockonly).
    rec: optional report.RunReport collecting per-sequence results (run report, summary lines).
    keep_going: skip sequences that cannot be converted (the others are written in order); returns how many.
    progress: percent done on stderr for inputs over PROGRESS_MIN_SEQS sequences.
    resume: continue an interrupted run with the same out_path from its checkpoint (checkpoint.py); note(text)
    receives informational messages."""
    if out_format not in FORMATS:
        raise CoalError(f"unknown output format '{out_format}' (one of {', '.join(FORMATS)})")
    if out_format == "codon" and (nogap or nomismatch or blockonly):
        raise CoalError("-output codon cannot be combined with -nogap, -nomismatch or -blockonly (as in PAL2NAL)")
    check_table(table)  # CoalError `unknown_table` (the command line never gets here: argparse checks the choices)
    rows = _protein_rows(pep_path)
    # PAL2NAL reads '#' block marks from a CLUSTAL input's conservation lines: -blockonly keeps only those columns,
    # and clustal / codon output show the marks under each block (when -blockonly is not given)
    mask = clustal_mask(pep_path) if is_clustal(pep_path) else ""
    has_marks = "#" in mask
    block_cols = blockonly and has_marks

    # pass 1: protein IDs, row length, frameshift column widths
    ids, ncol = [], None
    digit_rows = []
    last_row = ""
    for pid, row in rows():
        ids.append(pid)
        last_row = row
        if ncol is None:
            ncol = len(row)
            first_row = row
        elif len(row) != ncol:
            # The files-swapped case (the first file is DNA) and the not-aligned case are told apart, with what to do.
            raise CoalError(row_length_message(str(pep_path), ncol, len(row), pid, first_row + row),
                            code="row_lengths")
        if _DIGIT.search(row):
            digit_rows.append(row)
    if not ids:
        raise CoalError(f"{pep_path}: no sequences found: the file is empty or has no '>' header lines (expected a "
                        f"protein alignment in FASTA or CLUSTAL format)", code="no_sequences")
    if rec is not None:
        rec.set_inputs(ids)
    dna = DnaIndex(nuc_paths)
    try:
        if not len(dna):
            raise CoalError(f"{', '.join(map(str, nuc_paths))}: no sequences found: the file"
                            + ("s are" if len(nuc_paths) > 1 else " is")
                            + " empty or ha" + ("ve" if len(nuc_paths) > 1 else "s")
                            + " no '>' header lines (expected the coding DNA in FASTA format)", code="no_sequences")
        if rec is not None:
            per_file = [0] * len(dna.paths)
            for f in dna.fno:
                per_file[f] += 1
            rec.set_inputs(ids, per_file, dna.ids)
        # DNA record number for each protein, as a compact array; then free the ID tables
        # pair_ids reports a partial ID match and uses normalised IDs where file order would mispair
        pair_notes: list = []
        idx, mode, per_seq = pair_ids(ids, dna.ids, dna.pos, warn=pair_notes.append, strict=strict)
        pairing = array("q", idx)
        del idx
        if notice:
            for n in pair_notes:
                notice(n)
        nseq = len(ids)
        # Sequences paired by file order whose IDs differ: a failure names the DNA record used (a run where every
        # gene fits stays silent).
        order_why = None
        if mode == "order":
            hits = sum(1 for p in ids if p in dna.pos)
            order_why = ("because no protein ID matches a DNA ID" if not hits else
                         f"because only {hits} of {nseq} protein IDs are DNA IDs")
            order_dna = [dna.ids[j] if dna.ids[j] != ids[k] else None for k, j in enumerate(pairing)]

        def _paired(k: int, code: str, msg: str) -> str:
            if order_why and code == "inconsistent" and order_dna[k] is not None:
                return f"{msg} (paired by file order with DNA '{order_dna[k]}', {order_why})"
            return msg
        if rec is not None:
            rec.set_pairing(mode, per_seq, [dna.ids[j] for j in pairing], pair_notes)
        block_ids = ids if out_format in BLOCK_FORMATS else None  # interleaved formats need every ID at the end
        del ids
        dna.forget_ids()
        jobs = effective_jobs(jobs, dna.total_bytes(), min_nt_per_job)
        widths = column_widths(digit_rows, ncol)
        # --keep-going + paml stages the rows too (the header needs the final count)
        two_pass = (nogap or nomismatch or block_cols or out_format in BLOCK_FORMATS
                    or (keep_going and out_format == "paml"))
        keep_col = (lambda c: mask[c - 1:c] == "#") if block_cols else None

        # Temp-file-then-rename only onto a regular/new file (symlinks resolved; a symlink loop is a clear error);
        # FIFOs, devices and files with other hard links (renaming would break the links) are written directly.
        # An existing file keeps its permission bits.
        if out_path and os.path.islink(out_path):
            try:
                os.stat(out_path)
            except OSError as e:
                if e.errno == errno.ELOOP or getattr(e, "winerror", None) == 1921:
                    raise CoalError(f"{out_path}: symbolic link loop; cannot write the output there",
                                    code="file_error") from None
        real_out = os.path.realpath(out_path) if out_path else None
        # -o naming a folder or an input file is refused before any work.
        if real_out and os.path.isdir(real_out):
            raise CoalError(f"cannot write {out_path}: it is a folder, not a file (give a file name with -o)",
                            code="file_error")
        if real_out and os.path.isfile(real_out):
            for p in (pep_path, *nuc_paths):
                try:
                    same = os.path.samefile(p, real_out)
                except OSError:
                    same = False
                if same:
                    raise CoalError(f"cannot write {out_path}: it is one of the input files; choose another name for "
                                    f"the output", code="file_error")
        existing = bool(real_out) and os.path.isfile(real_out)
        direct = bool(real_out) and os.path.exists(real_out) and (
            not existing or os.stat(real_out).st_nlink > 1)

        # Checkpoints for --resume (checkpoint.py). Kept only for an output that is (or will be) a regular file and
        # inputs that are regular files: the partial output <real out>.coal-partial (rows: .coal-stage) lives next
        # to the real target, so symlinks, hard links and modes are handled as above; FIFOs/devices (e.g.
        # /dev/null) and piped inputs use hidden temp files and cannot be resumed.
        ck = state = None
        if out_path and _resumable and (not os.path.exists(real_out) or existing):
            ck = Checkpoint(real_out, {"fmt": out_format, "table": table, "nogap": nogap,
                                       "nomismatch": nomismatch, "blockonly": blockonly, "strict": strict,
                                       "keep_going": keep_going, "staged": two_pass},
                            [pep_path, *nuc_paths])
            if resume:
                state = ck.load_for_resume(ck.stage if two_pass else ck.partial)
                if note:
                    if state:
                        note(f"resuming {out_path}: {state['completed']:,} of {nseq:,} sequences were already done"
                             + (" (this format/option needs a second pass, which is redone in full)"
                                if two_pass else ""))
                    else:
                        note(("no checkpoint to resume" if ck.why == "no checkpoint" else
                              f"found {ck.why}, so nothing can be resumed") + f" for {out_path}; starting a normal run")
            if state is None:
                ck.remove()  # a stale checkpoint must never describe this run's files
        elif resume and note:
            note(f"{out_path} cannot be resumed (" + ("an input is not a regular file" if not _resumable else
                 "the output is not a regular file") + "); starting a normal run")

        dest_dir = os.path.dirname(real_out) if out_path and not direct else None
        scratch: list = []  # fixed-width temp files of the clustal / codon writer
        tmp = stage = None
        entries_fh = None
        finished = failed = False
        try:
            if ck is not None:
                try:
                    tmp_name = ck.partial
                    stage = open_append(ck.stage, state["sink_bytes"] if state else None) if two_pass else None
                    tmp = (open(ck.partial, "w", newline="\n", encoding="ascii") if two_pass else
                           open_append(ck.partial, state["sink_bytes"] if state else None))
                    entries_fh = open_append(ck.entries, state["entries_bytes"] if state else None)
                except OSError as e:
                    raise _output_error(out_path, e) from None
            else:
                try:
                    t = tempfile.NamedTemporaryFile("w", delete=False, newline="\n", dir=dest_dir,
                                                    prefix=".coal_", suffix=".tmp", encoding="ascii")
                except OSError as e:
                    raise _output_error(out_path, e) from None
                tmp, tmp_name = t.file, t.name
                scratch.append(t)
                if two_pass:
                    try:
                        st = tempfile.NamedTemporaryFile("w+", delete=False, newline="\n", prefix=".coal_rows_",
                                                         encoding="ascii")
                    except OSError as e:
                        raise _output_error(None, e) from None
                    scratch.append(st)
                    stage = st.file
            sink = stage if two_pass else tmp
            start = state["completed"] if state else 0
            if state:
                drop, bad, skipped = set(state["drop"]), set(state["bad"]), state["skipped"]
                skip_idx = set(state.get("skip_idx", ()))
                if rec is not None:
                    with open(ck.entries, encoding="utf-8") as fh:
                        rec.restore(state["rec"], [json.loads(x) for x in fh], start)
            else:
                drop, bad, skipped, skip_idx = set(), set(), 0, set()
                if not two_pass and out_format == "paml":
                    tmp.write(paml_header(nseq, sum(widths)))
            init = (table, widths, two_pass, nogap, strict, out_format, keep_going)

            def _shown(warns):  # the warnings -blockonly keeps (as on stderr)
                return [w for w in warns if not (keep_col and (m := _ALNPOS.search(" " + w))
                                                 and not keep_col(int(m.group(1))))]
            if note and warn and not nomismatch:
                # The genetic-code suggestion comes first, before the per-codon warnings it explains (not only in the
                # last summary line, after e.g. 834 warnings on coxI).
                from .report import early_code_notice
                hint = early_code_notice(rows(), lambda k: dna.get(pairing[k]), table, nseq,
                                         rec.fit_memo if rec is not None else None)
                if hint:
                    note(hint)
            done = start
            prog = Progress(nseq) if progress and nseq > PROGRESS_MIN_SEQS else None
            last_ck = time.monotonic()
            results = _results(_batches(rows(), dna, pairing, start), jobs, init)
            while True:
                try:
                    text, items, mcols, badt = next(results)
                except StopIteration:
                    break
                except CoalError as e:
                    if e.seq_id is not None:
                        e.index = done + (e.index or 0)
                        msg = _paired(e.index, e.code, str(e))
                        if msg != str(e):
                            e2 = CoalError(msg, e.code, e.seq_id, e.index)
                            e2.done_items = getattr(e, "done_items", ())
                            e = e2
                        if rec is not None:
                            # the sequences before it in its batch were checked (not written): report only, stderr
                            # stays as it was (they were never printed)
                            for k, (pid, warns, notes, _f) in enumerate(getattr(e, "done_items", ())):
                                rec.add_sequence(done + k, pid, _shown(warns), notes, None)
                            rec.add_sequence(e.index, e.seq_id, [], [], (e.code, str(e)))
                    raise e from None
                for pid, warns, notes, failure in items:
                    if failure:
                        failure = (failure[0], _paired(done, *failure))
                    for w in warns:
                        printed = not (nomismatch and MISMATCH_WARNING in w)
                        if keep_col:  # -blockonly: as PAL2NAL, only warnings in kept columns
                            m = _ALNPOS.search(" " + w)
                            if m and not keep_col(int(m.group(1))):
                                continue
                        # As PAL2NAL, no mismatch warnings with -nomismatch.
                        if warn and printed:
                            warn(f"{pid} {w}")
                    if notice:
                        for n in notes:
                            notice(n)
                        if failure:  # --keep-going names each skipped sequence
                            notice(f"skipped {failure[1]}")
                    if rec is not None:
                        entry = rec.add_sequence(done, pid, _shown(warns), notes, failure)
                        if entries_fh is not None and entry is not None:
                            entries_fh.write(json.dumps(entry) + "\n")
                    if failure:
                        skipped += 1
                        if block_ids is not None:
                            skip_idx.add(done)
                    done += 1
                if prog:
                    prog.update(done)
                sink.write(text)
                if nomismatch:
                    drop.update(c - 1 for c in mcols)
                bad |= badt
                if ck is not None and time.monotonic() - last_ck >= _ck.CHECKPOINT_EVERY and done < nseq:
                    ck.save({"completed": done, "nseq": nseq, "sink_bytes": durable_size(sink),
                             "entries_bytes": durable_size(entries_fh), "skipped": skipped,
                             "skip_idx": sorted(skip_idx), "drop": sorted(drop), "bad": sorted(bad),
                             "rec": rec.checkpoint_state() if rec is not None else None})
                    last_ck = time.monotonic()
                    _test_pause(ck.count)
            if prog and prog.printed:
                prog.update(done, final=True)
            if keep_going and skipped == nseq:
                raise CoalError(f"none of the {nseq} sequences could be converted: check that the DNA file holds the "
                                f"coding sequences of these proteins and that the genetic code (-codontable) is right",
                                code="all_failed")
            if rec is not None:
                rec.finish(rows, dna, pairing)
            if two_pass:
                if block_cols:  # -blockonly drops the columns not marked '#'
                    drop.update(i for i in range(ncol) if mask[i:i + 1] != "#")
                sl = keep_slices(widths, drop, bad)
                alilen = sum(b - a for a, b in sl)
                stage.flush()
                with open(stage.name, "r", newline="\n", encoding="ascii") as rows_in:
                    if out_format in BLOCK_FORMATS:
                        kept_ids = [i for k, i in enumerate(block_ids) if k not in skip_idx]
                        _write_block_output(tmp, out_format, kept_ids, rows_in, sl, alilen, rows,
                                            codon_aa_widths(ncol, digit_rows, last_row),
                                            _expand_mask(mask, widths, sl) if has_marks and not blockonly else None,
                                            scratch, skip_idx)
                    else:
                        if out_format == "paml":
                            tmp.write(paml_header(nseq - skipped, alilen))
                        record = paml_record if out_format == "paml" else fasta_record
                        for line in rows_in:
                            pid, exp = line.rstrip("\n").split("\t", 1)
                            tmp.write(record(pid, apply_slices(exp, sl)))
            tmp.close()
            if rec is not None:
                from .report import sha256_file
                # hashing a big output takes time: only when a report was asked for
                rec.set_output(out_path, sha256_file(tmp_name) if rec.keep_entries else None, nseq - skipped)
            if out_path and not direct:
                if existing:
                    shutil.copymode(real_out, tmp_name)  # keep the existing file's mode
                else:
                    umask = os.umask(0)
                    os.umask(umask)
                    os.chmod(tmp_name, 0o666 & ~umask)  # not mkstemp's 0600
                try:
                    os.replace(tmp_name, real_out)
                except OSError as e:
                    raise _output_error(out_path, e) from None
            elif out_path:
                with open(tmp_name, "r", encoding="ascii", newline="\n") as fh:
                    try:
                        dst = open(out_path, "w", encoding="ascii", newline="\n")
                    except OSError as e:
                        raise _output_error(out_path, e) from None
                    with dst:
                        shutil.copyfileobj(fh, dst, 1 << 20)
            else:
                with open(tmp_name, "r", encoding="ascii", newline="\n") as fh:
                    out = sys.stdout
                    if hasattr(out, "reconfigure"):
                        out.reconfigure(newline="\n")
                    try:
                        shutil.copyfileobj(fh, out, 1 << 20)
                        out.flush()
                    except OSError as e:
                        # On Windows, EINVAL from writing to a closed stdout pipe (how Windows often reports it) is a
                        # quiet exit 141 too, as BrokenPipeError.
                        if _WINDOWS and not isinstance(e, BrokenPipeError) and e.errno == errno.EINVAL:
                            raise BrokenPipeError(errno.EPIPE, "stdout closed") from None
                        raise
            finished = True
            return skipped
        except CoalError:
            failed = True  # a real error: nothing to resume
            if rec is not None and not rec.finished:
                rec.finish(rows, dna, pairing, aborted=True)
            raise
        finally:
            for f in (tmp, stage, entries_fh):
                if f is not None:
                    f.close()
            for f in scratch:
                f.close()
                try:
                    os.unlink(f.name)
                except OSError:
                    pass
            if ck is not None:
                # After a crash, kill or signal (cleanup) the partial files are kept only when a checkpoint describes
                # them, so --resume can continue; otherwise they are removed
                if finished or failed or not os.path.exists(ck.path):
                    ck.remove()
    finally:
        dna.close()
