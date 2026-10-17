"""Streaming, optionally parallel conversion for large inputs.

Memory use does not grow with the number of sequences beyond their IDs:
  pass 1  read the protein alignment once: IDs, row length, frameshift columns
          index the DNA file(s): IDs and byte offsets only
  pass 2  read protein rows again one at a time, fetch the matching DNA from
          disk, thread codons (in worker processes if jobs > 1, in order), and
          write each output row as soon as it is ready.
With -nogap / -nomismatch, which depend on every row, pass 2 writes rows to a
temporary file and a final pass removes the flagged columns.

Output goes to a temporary file first and is only moved/copied to the real
destination when everything succeeded, so an error never leaves partial output
(PAL2NAL also prints nothing on error).
"""
from __future__ import annotations

import json
import mmap
import os
import time
import re
from array import array
import shutil
import sys
import tempfile
from collections import deque
from itertools import islice

from .checkpoint import CHECKPOINT_EVERY, Checkpoint, durable_size, open_append
from .core import (
    CoalError, apply_slices, column_widths, expand, keep_slices, nogap_bad_triplets,
    order_pairing_warning, pair_ids, prepare_row, thread_one,
)
from .seqio import DnaIndex, is_clustal, iter_fasta, read_clustal
from .tables import get_code

WRAP = 60
FORMATS = ("fasta", "clustal", "paml", "codon")
BLOCK_FORMATS = ("clustal", "codon")   # interleaved blocks: need every row before writing
CLUSTAL_BLOCK = 60      # nucleotides per CLUSTAL block
CODON_BLOCK = 20        # residues / codons per "-output codon" block
_DIGIT = re.compile(r"[0-9]")
BATCH_NT = int(os.environ.get("COAL_BATCH_NT", 2_000_000))  # target DNA per work batch (keeps IPC overhead low)
# Starting worker processes costs ~0.1-0.2 s (more on Windows), about the time one
# core needs for a few million nucleotides. So each extra worker must have at
# least this much DNA to do, otherwise COAL runs in a single process.
MIN_NT_PER_JOB = 4_000_000
PROGRESS_MIN_SEQS = 10_000   # --progress only reports on inputs larger than this
PROGRESS_EVERY = 2.0         # seconds between progress reports


class Progress:
    """Percent done + estimated time left on stderr: one updating line on a
    terminal, plain lines otherwise. Never touches stdout."""

    def __init__(self, total: int, stream=None):
        import time
        self.time = time.monotonic
        self.total, self.stream = total, stream or sys.stderr
        self.t0 = self.last = self.time()
        self.tty = bool(getattr(self.stream, "isatty", lambda: False)())
        self.printed = False

    def update(self, done: int, final: bool = False):
        now = self.time()
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
        if self.tty:
            self.stream.write("\r" + msg.ljust(78) + ("\n" if final else ""))
        else:
            self.stream.write(msg + "\n")
        self.stream.flush()
        self.printed = True


def _test_pause(count: int) -> None:
    """Test hook: COAL_TEST_PAUSE_AFTER_CHECKPOINT=k stops here after the k-th checkpoint
    (so a test can kill the run at a known point)."""
    if os.environ.get("COAL_TEST_PAUSE_AFTER_CHECKPOINT") == str(count):
        print(f"coal: test: paused after checkpoint {count}", file=sys.stderr, flush=True)
        time.sleep(3600)


def effective_jobs(requested: int, total_nt: int, min_nt_per_job: int | None = None) -> int:
    """Number of processes to actually use. requested: 0 = all cores, N = at most N."""
    cores = os.cpu_count() or 1
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


def _init_worker(table, widths, two_pass, nogap, fmt="fasta", keep_going=False):
    _W.update(table=table, widths=widths, wide=any(w != 3 for w in widths),
              two_pass=two_pass, nogap=nogap, fmt=fmt, keep_going=keep_going)


def _work(batch):
    """Thread one batch. Returns (text, warnings, mismatch columns, nogap triplets, items):
    text is finished FASTA, or 'id<TAB>row' lines for the column-removal pass;
    items = [(id, None | (error code, message))] for every sequence, in order.
    Without keep_going a CoalError is raised, with .index = position in the batch."""
    parts, warns, mcols, bad, items = [], [], set(), set(), []
    keep_going = _W.get("keep_going", False)
    two_pass = _W["two_pass"]
    fmt = _W.get("fmt", "fasta")
    record = paml_record if fmt == "paml" else fasta_record
    for i, (pid, row, dna) in enumerate(batch):
        try:
            it = thread_one(pid, row, dna, _W["table"])
        except CoalError as e:
            if keep_going and e.seq_id is not None:
                items.append((pid, (e.code, str(e))))
                continue
            e.index = i
            raise
        items.append((pid, None))
        exp = expand(it, _W["widths"], _W["wide"])
        warns += [f"{pid} {w}" for w in it.warnings]
        if two_pass:
            # staged for a later pass; "-output codon" also needs the protein row
            parts.append(f"{pid}\t{row}\t{exp}\n" if fmt == "codon" else f"{pid}\t{exp}\n")
            mcols |= it.mismatch_cols
            if _W["nogap"]:
                bad |= nogap_bad_triplets(exp)
        else:
            parts.append(record(pid, exp))
    return "".join(parts), warns, mcols, bad, items


# ------------------------------------------------------------------ main side

def fasta_record(rid: str, seq: str) -> str:
    """One FASTA record wrapped at 60 (an empty sequence gives an empty line, as PAL2NAL)."""
    if len(seq) <= WRAP:
        return f">{rid}\n{seq}\n"
    return f">{rid}\n" + "\n".join([seq[i:i + WRAP] for i in range(0, len(seq), WRAP)]) + "\n"


def _wrap(seq: str) -> str:
    if len(seq) <= WRAP:
        return seq + "\n"
    return "\n".join([seq[i:i + WRAP] for i in range(0, len(seq), WRAP)]) + "\n"


def paml_record(rid: str, seq: str) -> str:
    """PAML (PHYLIP-like) record: ID line, then the sequence wrapped at 60."""
    return f"{rid}\n" + _wrap(seq)


def paml_header(nseq: int, length: int) -> str:
    return f"{nseq:4d}{length:7d}\n"


def name_width(ids) -> int:
    """Name column of CLUSTAL / codon output: longest ID (at least 10) + 4 spaces."""
    return max(10, max(map(len, ids), default=0)) + 4


def _view(f, size):
    f.flush()
    return mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) if size else b""


def write_blocks(stage, out, fmt: str, slices, length: int, ncol: int) -> None:
    """Write CLUSTAL or codon output from the staged rows.

    Every output row has the same length, so the rows are copied into a
    fixed-stride scratch file and read back block by block through mmap:
    memory holds the IDs only, never the whole alignment."""
    ids = []
    grid = tempfile.TemporaryFile("w+b", prefix=".coal_grid_")
    prot = tempfile.TemporaryFile("w+b", prefix=".coal_prot_") if fmt == "codon" else None
    try:
        for line in stage:
            parts = line.rstrip("\n").split("\t")
            ids.append(parts[0])
            grid.write(apply_slices(parts[-1], slices).encode("ascii"))
            if prot is not None:
                prot.write(parts[1].encode("ascii"))
        n = len(ids)
        gm = _view(grid, n * length)
        width = name_width(ids)
        names = [i.ljust(width) for i in ids]
        if fmt == "clustal":
            out.write("CLUSTAL W multiple sequence alignment\n\n")
            for s in range(0, length, CLUSTAL_BLOCK):
                e = min(s + CLUSTAL_BLOCK, length)
                out.write("".join([f"{names[k]}{gm[k * length + s:k * length + e].decode('ascii')}\n"
                                   for k in range(n)]))
                out.write("\n")
        else:
            pm = _view(prot, n * ncol)
            pad = " " * width
            ntrip = length // 3
            # The residue line and the codon line are cut into blocks of 20 separately
            # (as PAL2NAL does), so a widened frameshift column shifts the residues.
            nblk = max(-(-ncol // CODON_BLOCK), -(-ntrip // CODON_BLOCK))
            for b in range(nblk):
                a0, a1 = min(b * CODON_BLOCK, ncol), min((b + 1) * CODON_BLOCK, ncol)
                d0, d1 = min(3 * b * CODON_BLOCK, length), min(3 * (b + 1) * CODON_BLOCK, length)
                lines = []
                for k in range(n):
                    aa = pm[k * ncol + a0:k * ncol + a1].decode("ascii")
                    dna = gm[k * length + d0:k * length + d1].decode("ascii")
                    lines.append(pad + "   ".join(aa) + "\n")
                    lines.append(names[k] + " ".join([dna[i:i + 3] for i in range(0, len(dna), 3)]) + "\n")
                out.write("".join(lines))
                out.write("\n")
            if isinstance(pm, mmap.mmap):
                pm.close()
        if isinstance(gm, mmap.mmap):
            gm.close()
    finally:
        grid.close()
        if prot is not None:
            prot.close()


def _batches(rows, dna: DnaIndex, pairing, start: int = 0):
    """Batches of (id, row, DNA) of about BATCH_NT nucleotides, from sequence `start` on."""
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


def run(pep_path, nuc_paths, out_path=None, table: int = 1, nogap: bool = False,
        nomismatch: bool = False, jobs: int = 0, warn=None, min_nt_per_job: int | None = None,
        fmt: str = "fasta", rec=None, keep_going: bool = False, progress: bool = False,
        resume: bool = False, note=None) -> int:
    """Convert; write the alignment (fmt: fasta, clustal, paml or codon) to out_path (or stdout).

    jobs: 0 = automatic (up to all cores), N = at most N processes. Small inputs
    always run in one process (see MIN_NT_PER_JOB). rec: optional report.RunReport
    that collects per-sequence results for the run report.
    keep_going: skip sequences that cannot be threaded (written: all the others, in
    order). Returns the number of skipped sequences (0 without keep_going).
    progress: report percent done on stderr for inputs over PROGRESS_MIN_SEQS sequences.
    resume: continue an interrupted run with the same out_path from its checkpoint
    (see checkpoint.py); note(text) receives informational messages."""
    get_code(table)
    if fmt not in FORMATS:
        raise CoalError(f"unknown output format {fmt!r}; one of {', '.join(FORMATS)}")
    if fmt == "codon" and (nogap or nomismatch):
        raise CoalError('"-output codon" cannot be combined with -nogap or -nomismatch (as in PAL2NAL)')
    rows = _protein_rows(pep_path)

    # pass 1: protein IDs, row length, frameshift column widths
    ids, ncol = [], None
    digit_rows = []
    for pid, row in rows():
        ids.append(pid)
        if ncol is None:
            ncol = len(row)
        elif len(row) != ncol:
            raise CoalError(f"protein alignment rows have different lengths ({ncol} and {len(row)}, at {pid})",
                            code="row_lengths")
        if _DIGIT.search(row):
            digit_rows.append(row)
    dna = DnaIndex(nuc_paths)
    try:
        if rec is not None:
            per_file = [0] * len(dna.paths)
            for f in dna.fno:
                per_file[f] += 1
            rec.set_inputs(len(ids), per_file)
        # DNA record number for each protein, as a compact array; then free the ID tables
        idx, mode, per_seq = pair_ids(ids, dna.ids, dna.pos)
        pairing = array("q", idx)
        del idx
        nseq = len(pairing)
        if mode == "order" and warn:
            warn(order_pairing_warning(ids, dna.ids))
        if rec is not None:
            example = None
            if per_seq:
                k = per_seq.index("normalised")
                example = (ids[k], dna.ids[pairing[k]])
            rec.set_pairing(mode, [dna.ids[j] for j in pairing] if rec.keep_entries else None, pairing, per_seq,
                            example)
        ids_for_resume = ids if resume else None
        del ids
        dna.forget_ids()
        jobs = effective_jobs(jobs, dna.total_bytes(), min_nt_per_job)
        widths = column_widths(digit_rows, ncol)
        two_pass = nogap or nomismatch
        # keep_going + paml: the header needs the final count, so stage the rows
        staged = two_pass or fmt in BLOCK_FORMATS or (keep_going and fmt == "paml")

        ck, state = None, None
        if out_path:
            ck = Checkpoint(out_path, {"fmt": fmt, "table": table, "nogap": nogap, "nomismatch": nomismatch,
                                       "keep_going": keep_going, "staged": staged}, [pep_path, *nuc_paths])
            out_tmp, stage_path = ck.partial, (ck.stage if staged else None)
            if resume:
                state = ck.load_for_resume(stage_path or out_tmp)
                if note:
                    if state:
                        note(f"resuming {out_path}: {state['completed']:,} of {nseq:,} sequences were already done"
                             + (" (this format/option needs a second pass, which is redone in full)" if staged else ""))
                    else:
                        note(f"no checkpoint to resume for {out_path}; starting a normal run")
            if state is None:
                ck.remove(keep_partial=True)   # a stale checkpoint must never describe this run's files
        else:
            fd, out_tmp = tempfile.mkstemp(prefix=".coal_", suffix=".tmp")
            os.close(fd)
            stage_path = None
            if staged:
                fd, stage_path = tempfile.mkstemp(prefix=".coal_rows_", suffix=".tmp")
                os.close(fd)
        sink_path = stage_path if staged else out_tmp
        start = state["completed"] if state else 0
        sink = entries_fh = None
        finished = failed = False
        try:
            sink = open_append(sink_path, state["sink_bytes"] if state else None)
            if ck:
                entries_fh = open_append(ck.entries, state["entries_bytes"] if state else None)
            if state:
                drop, bad, skipped = set(state["drop"]), set(state["bad"]), state["skipped"]
                if rec is not None:
                    with open(ck.entries, encoding="ascii") as fh:
                        nonok = [json.loads(l) for l in fh]
                    rec.restore(state["rec"], nonok, start, ids_for_resume)
            else:
                drop, bad, skipped = set(), set(), 0
                if fmt == "paml" and not staged:
                    sink.write(paml_header(nseq, sum(widths)))
            ids_for_resume = None
            init = (table, widths, staged, nogap, fmt, keep_going)
            done = start
            prog = Progress(nseq) if progress and nseq > PROGRESS_MIN_SEQS else None
            last_ck = time.monotonic()
            results = _results(_batches(rows(), dna, pairing, start), jobs, init)
            while True:
                try:
                    text, warns, mcols, badt, items = next(results)
                except StopIteration:
                    break
                except CoalError as e:
                    if rec is not None and e.seq_id is not None:
                        e.index = done + (e.index or 0)
                        rec.add_sequence(e.seq_id, [], failure=(e.code, str(e)), k=e.index)
                    raise
                wi = 0
                for pid, failure in items:
                    mine, pre = [], pid + " "
                    while wi < len(warns) and warns[wi].startswith(pre):
                        mine.append(warns[wi])
                        wi += 1
                    if warn:
                        for w in mine:
                            warn(w)
                        if failure:
                            warn(f"skipped {failure[1]}")
                    if rec is not None:
                        entry = rec.add_sequence(pid, [w[len(pre):] for w in mine], failure=failure, k=done)
                        if entries_fh is not None and entry is not None and entry["status"] != "ok":
                            entries_fh.write(json.dumps(entry) + "\n")
                    done += 1
                    skipped += failure is not None
                if prog:
                    prog.update(done)
                sink.write(text)
                if nomismatch:
                    drop.update(c - 1 for c in mcols)
                bad |= badt
                if ck and time.monotonic() - last_ck >= CHECKPOINT_EVERY and done < nseq:
                    ck.save({"completed": done, "nseq": nseq, "sink_bytes": durable_size(sink),
                             "entries_bytes": durable_size(entries_fh), "skipped": skipped,
                             "drop": sorted(drop), "bad": sorted(bad),
                             "rec": rec.checkpoint_state() if rec is not None else None})
                    last_ck = time.monotonic()
                    _test_pause(ck.count)
            if prog and prog.printed:
                prog.update(done, final=True)
            if keep_going and skipped == nseq:
                raise CoalError(f"none of the {nseq} sequences could be converted", code="all_failed")
            sink.close()
            if staged:
                if two_pass:
                    sl = keep_slices(widths, drop if nomismatch else (), bad)
                else:
                    sl = [(0, sum(widths))]
                length = sum(b - a for a, b in sl)
                with open(out_tmp, "w", newline="\n", encoding="ascii") as out, \
                        open(stage_path, "r", newline="\n", encoding="ascii") as st:
                    if fmt in BLOCK_FORMATS:
                        write_blocks(st, out, fmt, sl, length, ncol)
                    else:
                        if fmt == "paml":
                            out.write(paml_header(nseq - skipped, length))
                        record = paml_record if fmt == "paml" else fasta_record
                        for line in st:
                            pid, exp = line.rstrip("\n").split("\t", 1)
                            out.write(record(pid, apply_slices(exp, sl)))
            if rec is not None:
                from .report import sha256_file
                rec.set_output(out_path, sha256_file(out_tmp), nseq - skipped)
            if out_path:
                os.replace(out_tmp, out_path)
            else:
                with open(out_tmp, "r", encoding="ascii", newline="\n") as fh:
                    out = sys.stdout
                    if hasattr(out, "reconfigure"):
                        out.reconfigure(newline="\n")
                    shutil.copyfileobj(fh, out, 1 << 20)
                    out.flush()
            finished = True
            return skipped
        except CoalError:
            failed = True   # a real error: nothing to resume
            raise
        finally:
            for f in (sink, entries_fh):
                if f is not None:
                    f.close()
            if ck is not None:
                if finished or failed:
                    ck.remove()
                # otherwise (crash, Ctrl-C, disk full ...) keep partial + checkpoint for --resume
            else:
                for path in (out_tmp, stage_path):
                    if path:
                        try:
                            os.unlink(path)
                        except OSError:
                            pass
    finally:
        dna.close()
