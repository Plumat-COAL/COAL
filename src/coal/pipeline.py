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

import errno
import os
import re
from array import array
import shutil
import signal
import sys
import tempfile
from collections import deque

from .core import (
    MISMATCH_WARNING, CoalError, apply_slices, column_widths, expand, keep_slices, nogap_bad_triplets,
    pair_sequences, prepare_row, thread_one,
)
from .seqio import DnaIndex, is_clustal, iter_fasta, read_clustal
from .tables import get_code

WRAP = 60
_DIGIT = re.compile(r"[0-9]")
BATCH_NT = 2_000_000      # target DNA per work batch (keeps IPC overhead low)
# Starting worker processes costs ~0.1-0.2 s (more on Windows), about the time one
# core needs for a few million nucleotides. So each extra worker must have at
# least this much DNA to do, otherwise COAL runs in a single process.
MIN_NT_PER_JOB = 4_000_000
_WINDOWS = os.name == "nt"  # CHANGED(R2-D1): see the broken-pipe handling in _run


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


def mp_child() -> bool:
    import multiprocessing
    return multiprocessing.parent_process() is not None


def _init_worker(table, widths, two_pass, nogap):
    # CHANGED(ENG-f): workers ignore Ctrl-C (the main process reports it once, no worker
    # tracebacks). See overnight/engine-fixes/CHANGES.md
    if mp_child():
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    _W.update(table=table, widths=widths, wide=any(w != 3 for w in widths),
              two_pass=two_pass, nogap=nogap)


def _work(batch):
    """Thread one batch. Returns (text, warnings, mismatch columns, nogap triplets, notices):
    text is finished FASTA, or 'id<TAB>row' lines for the column-removal pass."""
    # CHANGED(ENG-c): was a 4-tuple; now also returns notices. See overnight/engine-fixes/CHANGES.md
    parts, warns, mcols, bad, notes = [], [], set(), set(), []
    two_pass = _W["two_pass"]
    for pid, row, dna in batch:
        it = thread_one(pid, row, dna, _W["table"])
        exp = expand(it, _W["widths"], _W["wide"])
        warns += [f"{pid} {w}" for w in it.warnings]
        notes += it.notices
        if two_pass:
            parts.append(f"{pid}\t{exp}\n")
            mcols |= it.mismatch_cols
            if _W["nogap"]:
                bad |= nogap_bad_triplets(exp)
        else:
            parts.append(fasta_record(pid, exp))
    return "".join(parts), warns, mcols, bad, notes


# ------------------------------------------------------------------ main side

def fasta_record(rid: str, seq: str) -> str:
    """One FASTA record wrapped at 60 (an empty sequence gives an empty line, as PAL2NAL)."""
    if len(seq) <= WRAP:
        return f">{rid}\n{seq}\n"
    return f">{rid}\n" + "\n".join([seq[i:i + WRAP] for i in range(0, len(seq), WRAP)]) + "\n"


def _batches(rows, dna: DnaIndex, pairing):
    batch, nt = [], 0
    for k, (pid, row) in enumerate(rows):
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
    """Convert; write FASTA to out_path (or stdout). Arguments as for _run."""
    # CHANGED(ENG-b): was opening every input path several times (pipes/FIFOs gave
    # "protein: 0" or "not seekable"); now non-regular inputs are copied to temp files
    # first and removed afterwards. See overnight/engine-fixes/CHANGES.md
    spooled: list = []
    try:
        pep_path = _spool(pep_path, spooled)
        nuc_paths = [_spool(p, spooled) for p in nuc_paths]
        return _run(pep_path, nuc_paths, *args, **kw)
    finally:
        for name in spooled:
            try:
                os.unlink(name)
            except OSError:
                pass


def _run(pep_path, nuc_paths, out_path=None, table: int = 1, nogap: bool = False,
         nomismatch: bool = False, jobs: int = 0, warn=None, min_nt_per_job: int | None = None,
         notice=None) -> int:
    """Convert; write FASTA to out_path (or stdout).

    jobs: 0 = automatic (up to all cores), N = at most N processes. Small inputs
    always run in one process (see MIN_NT_PER_JOB).
    warn / notice: callables for per-codon warnings / important notices (the CLI
    passes notice even with -q)."""
    get_code(table)
    rows = _protein_rows(pep_path)

    # pass 1: protein IDs, row length, frameshift column widths
    ids, ncol = [], None
    digit_rows = []
    for pid, row in rows():
        ids.append(pid)
        if ncol is None:
            ncol = len(row)
        elif len(row) != ncol:
            raise CoalError(f"protein alignment rows have different lengths ({ncol} and {len(row)}, at {pid})")
        if _DIGIT.search(row):
            digit_rows.append(row)
    dna = DnaIndex(nuc_paths)
    try:
        # DNA record number for each protein, as a compact array; then free the ID tables
        # CHANGED(ENG-d): passes notice, so a partial ID match is reported. See overnight/engine-fixes/CHANGES.md
        pairing = array("q", pair_sequences(ids, dna.ids, dna.pos, warn=notice))
        del ids
        dna.forget_ids()
        jobs = effective_jobs(jobs, dna.total_bytes(), min_nt_per_job)
        widths = column_widths(digit_rows, ncol)
        two_pass = nogap or nomismatch

        # CHANGED(ENG-a): was always temp-file-then-rename next to out_path (replaced symlinks,
        # FIFOs and devices; mode 0600); now renames only onto a regular/new file (symlinks
        # resolved), writes other paths directly. See overnight/engine-fixes/CHANGES.md
        real_out = os.path.realpath(out_path) if out_path else None
        direct = bool(real_out) and os.path.exists(real_out) and not os.path.isfile(real_out)
        dest_dir = os.path.dirname(real_out) if out_path and not direct else None
        tmp = tempfile.NamedTemporaryFile("w", delete=False, newline="\n", dir=dest_dir,
                                          prefix=".coal_", suffix=".tmp", encoding="ascii")
        stage = tempfile.NamedTemporaryFile("w+", delete=False, newline="\n", prefix=".coal_rows_",
                                            encoding="ascii") if two_pass else None
        try:
            drop, bad = set(), set()
            init = (table, widths, two_pass, nogap)
            sink = stage.file if two_pass else tmp.file   # skip tempfile's per-call wrapper
            for text, warns, mcols, badt, notes in _results(_batches(rows(), dna, pairing), jobs, init):
                if warn:
                    for w in warns:
                        # CHANGED(ENG-i2): was printed with -nomismatch too; now, as PAL2NAL, no
                        # mismatch warnings with -nomismatch. See overnight/engine-fixes/CHANGES.md
                        if not (nomismatch and MISMATCH_WARNING in w):
                            warn(w)
                if notice:  # CHANGED(ENG-c): new; see overnight/engine-fixes/CHANGES.md
                    for n in notes:
                        notice(n)
                sink.write(text)
                if nomismatch:
                    drop.update(c - 1 for c in mcols)
                bad |= badt
            if two_pass:
                sl = keep_slices(widths, drop if nomismatch else (), bad)
                stage.seek(0)
                for line in stage.file:
                    pid, exp = line.rstrip("\n").split("\t", 1)
                    tmp.file.write(fasta_record(pid, apply_slices(exp, sl)))
            tmp.close()
            if out_path and not direct:
                umask = os.umask(0)
                os.umask(umask)
                os.chmod(tmp.name, 0o666 & ~umask)  # CHANGED(ENG-a): was 0600; see CHANGES.md
                os.replace(tmp.name, real_out)
            elif out_path:
                with open(tmp.name, "r", encoding="ascii", newline="\n") as fh, \
                        open(out_path, "w", encoding="ascii", newline="\n") as dst:
                    shutil.copyfileobj(fh, dst, 1 << 20)
            else:
                with open(tmp.name, "r", encoding="ascii", newline="\n") as fh:
                    out = sys.stdout
                    if hasattr(out, "reconfigure"):
                        out.reconfigure(newline="\n")
                    try:
                        shutil.copyfileobj(fh, out, 1 << 20)
                        out.flush()
                    except OSError as e:
                        # CHANGED(R2-D1): was only BrokenPipeError quiet; now, on Windows, EINVAL from writing
                        # to a closed stdout pipe (how Windows often reports it) is a quiet exit 141 too.
                        # See overnight/integration/CHANGES.md
                        if _WINDOWS and not isinstance(e, BrokenPipeError) and e.errno == errno.EINVAL:
                            raise BrokenPipeError(errno.EPIPE, "stdout closed") from None
                        raise
            return 0
        finally:
            for f in (tmp, stage):
                if f is not None:
                    f.close()
                    try:
                        os.unlink(f.name)
                    except OSError:
                        pass
    finally:
        dna.close()
