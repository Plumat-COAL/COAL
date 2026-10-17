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

import mmap
import os
import re
from array import array
import shutil
import sys
import tempfile
from collections import deque

from .core import (
    CoalError, apply_slices, column_widths, expand, keep_slices, nogap_bad_triplets,
    pair_sequences, prepare_row, thread_one,
)
from .seqio import DnaIndex, is_clustal, iter_fasta, read_clustal
from .tables import get_code

WRAP = 60
FORMATS = ("fasta", "clustal", "paml", "codon")
BLOCK_FORMATS = ("clustal", "codon")   # interleaved blocks: need every row before writing
CLUSTAL_BLOCK = 60      # nucleotides per CLUSTAL block
CODON_BLOCK = 20        # residues / codons per "-output codon" block
_DIGIT = re.compile(r"[0-9]")
BATCH_NT = 2_000_000      # target DNA per work batch (keeps IPC overhead low)
# Starting worker processes costs ~0.1-0.2 s (more on Windows), about the time one
# core needs for a few million nucleotides. So each extra worker must have at
# least this much DNA to do, otherwise COAL runs in a single process.
MIN_NT_PER_JOB = 4_000_000


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


def _init_worker(table, widths, two_pass, nogap, fmt="fasta"):
    _W.update(table=table, widths=widths, wide=any(w != 3 for w in widths),
              two_pass=two_pass, nogap=nogap, fmt=fmt)


def _work(batch):
    """Thread one batch. Returns (text, warnings, mismatch columns, nogap triplets, IDs):
    text is finished FASTA, or 'id<TAB>row' lines for the column-removal pass.
    A CoalError gets .index = position of the failing sequence in the batch."""
    parts, warns, mcols, bad, pids = [], [], set(), set(), []
    two_pass = _W["two_pass"]
    fmt = _W.get("fmt", "fasta")
    record = paml_record if fmt == "paml" else fasta_record
    for i, (pid, row, dna) in enumerate(batch):
        try:
            it = thread_one(pid, row, dna, _W["table"])
        except CoalError as e:
            e.index = i
            raise
        pids.append(pid)
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
    return "".join(parts), warns, mcols, bad, pids


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


def run(pep_path, nuc_paths, out_path=None, table: int = 1, nogap: bool = False,
        nomismatch: bool = False, jobs: int = 0, warn=None, min_nt_per_job: int | None = None,
        fmt: str = "fasta", rec=None) -> int:
    """Convert; write the alignment (fmt: fasta, clustal, paml or codon) to out_path (or stdout).

    jobs: 0 = automatic (up to all cores), N = at most N processes. Small inputs
    always run in one process (see MIN_NT_PER_JOB). rec: optional report.RunReport
    that collects per-sequence results for the run report."""
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
        pairing = array("q", pair_sequences(ids, dna.ids, dna.pos))
        nseq = len(pairing)
        if rec is not None:
            mode = "id" if all(p in dna.pos for p in ids) else "order"
            rec.set_pairing(mode, [dna.ids[j] for j in pairing] if rec.keep_entries else None, pairing)
        del ids
        dna.forget_ids()
        jobs = effective_jobs(jobs, dna.total_bytes(), min_nt_per_job)
        widths = column_widths(digit_rows, ncol)
        two_pass = nogap or nomismatch
        staged = two_pass or fmt in BLOCK_FORMATS

        dest_dir = os.path.dirname(os.path.abspath(out_path)) if out_path else None
        tmp = tempfile.NamedTemporaryFile("w", delete=False, newline="\n", dir=dest_dir,
                                          prefix=".coal_", suffix=".tmp", encoding="ascii")
        stage = tempfile.NamedTemporaryFile("w+", delete=False, newline="\n", prefix=".coal_rows_",
                                            encoding="ascii") if staged else None
        try:
            drop, bad = set(), set()
            init = (table, widths, staged, nogap, fmt)
            sink = stage.file if staged else tmp.file   # skip tempfile's per-call wrapper
            if fmt == "paml" and not staged:
                tmp.file.write(paml_header(nseq, sum(widths)))
            done = 0
            results = _results(_batches(rows(), dna, pairing), jobs, init)
            while True:
                try:
                    text, warns, mcols, badt, pids = next(results)
                except StopIteration:
                    break
                except CoalError as e:
                    if rec is not None and e.seq_id is not None:
                        e.index = done + (e.index or 0)
                        rec.add_sequence(e.seq_id, [], failure=(e.code, str(e)), k=e.index)
                    raise
                if warn:
                    for w in warns:
                        warn(w)
                if rec is not None:
                    wi = 0
                    for pid in pids:
                        mine, pre = [], pid + " "
                        while wi < len(warns) and warns[wi].startswith(pre):
                            mine.append(warns[wi][len(pre):])
                            wi += 1
                        rec.add_sequence(pid, mine, k=done)
                        done += 1
                else:
                    done += len(pids)
                sink.write(text)
                if nomismatch:
                    drop.update(c - 1 for c in mcols)
                bad |= badt
            if staged:
                if two_pass:
                    sl = keep_slices(widths, drop if nomismatch else (), bad)
                else:
                    sl = [(0, sum(widths))]
                length = sum(b - a for a, b in sl)
                stage.seek(0)
                if fmt in BLOCK_FORMATS:
                    write_blocks(stage.file, tmp.file, fmt, sl, length, ncol)
                else:
                    if fmt == "paml":
                        tmp.file.write(paml_header(nseq, length))
                    record = paml_record if fmt == "paml" else fasta_record
                    for line in stage.file:
                        pid, exp = line.rstrip("\n").split("\t", 1)
                        tmp.file.write(record(pid, apply_slices(exp, sl)))
            tmp.close()
            if rec is not None:
                from .report import sha256_file
                rec.set_output(out_path, sha256_file(tmp.name), nseq)
            if out_path:
                os.replace(tmp.name, out_path)
            else:
                with open(tmp.name, "r", encoding="ascii", newline="\n") as fh:
                    out = sys.stdout
                    if hasattr(out, "reconfigure"):
                        out.reconfigure(newline="\n")
                    shutil.copyfileobj(fh, out, 1 << 20)
                    out.flush()
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
