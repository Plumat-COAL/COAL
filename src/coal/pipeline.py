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
fixed-width temporary file, which is then read back a group of blocks at a time (memory stays
bounded by BLOCK_MEMORY, not by the number or length of the sequences).

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
from .seqio import DnaIndex, InputFormatError, clustal_mask, is_clustal, iter_fasta, read_clustal
from .tables import get_code

WRAP = 60
_DIGIT = re.compile(r"[0-9]")
# CHANGED(DEC-4): new. Was FASTA only; now PAL2NAL v14's four output formats. See overnight/integration/CHANGES.md
FORMATS = ("fasta", "clustal", "paml", "codon")
BLOCK_FORMATS = ("clustal", "codon")  # interleaved blocks of 60 nt
BLOCK_MEMORY = 16_000_000  # bytes of rows held at once while writing clustal / codon blocks
_ALNPOS = re.compile(r" pepAlnPos (\d+):")
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


def _init_worker(table, widths, two_pass, nogap, strict=False, out_format="fasta"):
    # CHANGED(ENG-f): workers ignore Ctrl-C (the main process reports it once, no worker
    # tracebacks). See overnight/engine-fixes/CHANGES.md
    if mp_child():
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    _W.update(table=table, widths=widths, wide=any(w != 3 for w in widths),
              two_pass=two_pass, nogap=nogap, strict=strict,  # CHANGED(DEC-1): strict is new
              record=paml_record if out_format == "paml" else fasta_record)  # CHANGED(DEC-4)


def _work(batch):
    """Thread one batch. Returns (text, warnings, mismatch columns, nogap triplets, notices):
    text is finished FASTA, or 'id<TAB>row' lines for the column-removal pass."""
    # CHANGED(ENG-c): was a 4-tuple; now also returns notices. See overnight/engine-fixes/CHANGES.md
    parts, warns, mcols, bad, notes = [], [], set(), set(), []
    two_pass = _W["two_pass"]
    for pid, row, dna in batch:
        it = thread_one(pid, row, dna, _W["table"], _W["strict"])
        exp = expand(it, _W["widths"], _W["wide"])
        warns += [f"{pid} {w}" for w in it.warnings]
        notes += it.notices
        if two_pass:
            parts.append(f"{pid}\t{exp}\n")
            mcols |= it.mismatch_cols
            if _W["nogap"]:
                bad |= nogap_bad_triplets(exp)
        else:
            parts.append(_W["record"](pid, exp))
    return "".join(parts), warns, mcols, bad, notes


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
    group = max(1, memory // max(1, n * (WRAP + (20 if aa_fh else 0))))
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


def _expand_mask(mask: str, widths, slices) -> str:
    """The '#' mark line under clustal / codon blocks, one character per output nucleotide (as PAL2NAL)."""
    full = "".join(("#" if mask[i:i + 1] == "#" else " ") * w for i, w in enumerate(widths))
    return apply_slices(full, slices)


def _write_block_output(out, fmt, ids, stage, slices, alilen, rows, aa_widths, mask, scratch) -> None:
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
        for _pid, row in rows():
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
    """Convert; write FASTA to out_path (or stdout). Arguments as for _run."""
    # CHANGED(ENG-b): was opening every input path several times (pipes/FIFOs gave
    # "protein: 0" or "not seekable"); now non-regular inputs are copied to temp files
    # first and removed afterwards. See overnight/engine-fixes/CHANGES.md
    spooled: list = []
    names = {}  # CHANGED(R2-E7): temp copy -> the argument it came from, for error messages
    try:
        sp = _spool(pep_path, spooled)
        names[sp] = str(pep_path)
        nps = []
        for p in nuc_paths:
            nps.append(_spool(p, spooled))
            names[nps[-1]] = str(p)
        try:
            return _run(sp, nps, *args, **kw)
        except (CoalError, InputFormatError) as e:
            msg = str(e)
            for tmp_name in spooled:
                msg = msg.replace(tmp_name, names[tmp_name])
            if msg != str(e):
                raise type(e)(msg) from None
            raise
    finally:
        for name in spooled:
            try:
                os.unlink(name)
            except OSError:
                pass


def _run(pep_path, nuc_paths, out_path=None, table: int = 1, nogap: bool = False,
         nomismatch: bool = False, jobs: int = 0, warn=None, min_nt_per_job: int | None = None,
         notice=None, strict: bool = False, out_format: str = "fasta", blockonly: bool = False) -> int:
    """Convert; write the codon alignment to out_path (or stdout) in out_format (FORMATS).

    jobs: 0 = automatic (up to all cores), N = at most N processes. Small inputs
    always run in one process (see MIN_NT_PER_JOB).
    warn / notice: callables for per-codon warnings / important notices.
    blockonly: keep only the columns marked '#' in a CLUSTAL input's mark line (PAL2NAL -blockonly)."""
    # CHANGED(DEC-4): out_format and blockonly are new (was FASTA only, -blockonly rejected).
    # See overnight/integration/CHANGES.md
    if out_format not in FORMATS:
        raise CoalError(f"unknown output format '{out_format}' (one of {', '.join(FORMATS)})")
    if out_format == "codon" and (nogap or nomismatch or blockonly):
        raise CoalError("-output codon cannot be combined with -nogap, -nomismatch or -blockonly (as in PAL2NAL)")
    get_code(table)
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
        elif len(row) != ncol:
            raise CoalError(f"protein alignment rows have different lengths ({ncol} and {len(row)}, at {pid})")
        if _DIGIT.search(row):
            digit_rows.append(row)
    dna = DnaIndex(nuc_paths)
    try:
        # DNA record number for each protein, as a compact array; then free the ID tables
        # CHANGED(ENG-d): passes notice, so a partial ID match is reported. See overnight/engine-fixes/CHANGES.md
        pairing = array("q", pair_sequences(ids, dna.ids, dna.pos, warn=notice, strict=strict))  # CHANGED(DEC-1)
        nseq = len(ids)
        block_ids = ids if out_format in BLOCK_FORMATS else None  # interleaved formats need every ID at the end
        del ids
        dna.forget_ids()
        jobs = effective_jobs(jobs, dna.total_bytes(), min_nt_per_job)
        widths = column_widths(digit_rows, ncol)
        two_pass = nogap or nomismatch or block_cols or out_format in BLOCK_FORMATS
        keep_col = (lambda c: mask[c - 1:c] == "#") if block_cols else None

        # CHANGED(ENG-a): was always temp-file-then-rename next to out_path (replaced symlinks,
        # FIFOs and devices; mode 0600); now renames only onto a regular/new file (symlinks
        # resolved), writes other paths directly. See overnight/engine-fixes/CHANGES.md
        # CHANGED(R2-E5): a symlink loop was silently replaced by a regular file; now a clear error. An existing
        # file with other hard links is written in place (renaming broke the links), and an existing file keeps
        # its permission bits (was widened to 0666 & ~umask). See overnight/integration/CHANGES.md
        if out_path and os.path.islink(out_path):
            try:
                os.stat(out_path)
            except OSError as e:
                if e.errno == errno.ELOOP or getattr(e, "winerror", None) == 1921:
                    raise CoalError(f"{out_path}: symbolic link loop; cannot write the output there") from None
        real_out = os.path.realpath(out_path) if out_path else None
        existing = bool(real_out) and os.path.isfile(real_out)
        direct = bool(real_out) and os.path.exists(real_out) and (
            not existing or os.stat(real_out).st_nlink > 1)
        dest_dir = os.path.dirname(real_out) if out_path and not direct else None
        # CHANGED(R3-E2): was the raw OSError naming COAL's hidden temp file ("No such file or directory:
        # '.../nodir/.coal_x.tmp'"); now the error names the user's -o path in plain words.
        # See overnight/integration/CHANGES.md
        try:
            tmp = tempfile.NamedTemporaryFile("w", delete=False, newline="\n", dir=dest_dir,
                                              prefix=".coal_", suffix=".tmp", encoding="ascii")
        except OSError as e:
            folder = os.path.dirname(os.path.abspath(out_path)) if out_path else ""
            if e.errno == errno.ENOENT:
                why = f"folder {folder} does not exist"
            elif e.errno == errno.ELOOP:
                why = f"folder {folder} is a symbolic link loop"
            elif e.errno == errno.ENOTDIR:
                why = f"{folder} is not a folder"
            else:
                why = e.strerror or str(e)
            raise CoalError(f"cannot write {out_path}: {why}") from None
        stage = tempfile.NamedTemporaryFile("w+", delete=False, newline="\n", prefix=".coal_rows_",
                                            encoding="ascii") if two_pass else None
        scratch: list = []  # CHANGED(DEC-4): fixed-width temp files of the clustal / codon writer
        try:
            drop, bad = set(), set()
            init = (table, widths, two_pass, nogap, strict, out_format)
            sink = stage.file if two_pass else tmp.file   # skip tempfile's per-call wrapper
            if not two_pass and out_format == "paml":
                tmp.file.write(paml_header(nseq, sum(widths)))
            for text, warns, mcols, badt, notes in _results(_batches(rows(), dna, pairing), jobs, init):
                if warn:
                    for w in warns:
                        # CHANGED(ENG-i2): was printed with -nomismatch too; now, as PAL2NAL, no
                        # mismatch warnings with -nomismatch. See overnight/engine-fixes/CHANGES.md
                        if not (nomismatch and MISMATCH_WARNING in w):
                            if keep_col:  # CHANGED(DEC-4): -blockonly: as PAL2NAL, only warnings in kept columns
                                m = _ALNPOS.search(w)
                                if m and not keep_col(int(m.group(1))):
                                    continue
                            warn(w)
                if notice:  # CHANGED(ENG-c): new; see overnight/engine-fixes/CHANGES.md
                    for n in notes:
                        notice(n)
                sink.write(text)
                if nomismatch:
                    drop.update(c - 1 for c in mcols)
                bad |= badt
            if two_pass:
                if block_cols:  # CHANGED(DEC-4): -blockonly drops the columns not marked '#'
                    drop.update(i for i in range(ncol) if mask[i:i + 1] != "#")
                sl = keep_slices(widths, drop, bad)
                alilen = sum(b - a for a, b in sl)
                stage.seek(0)
                if out_format in BLOCK_FORMATS:
                    _write_block_output(tmp.file, out_format, block_ids, stage.file, sl, alilen, rows,
                                        codon_aa_widths(ncol, digit_rows, last_row),
                                        _expand_mask(mask, widths, sl) if has_marks and not blockonly else None,
                                        scratch)
                else:
                    if out_format == "paml":
                        tmp.file.write(paml_header(nseq, alilen))
                    record = paml_record if out_format == "paml" else fasta_record
                    for line in stage.file:
                        pid, exp = line.rstrip("\n").split("\t", 1)
                        tmp.file.write(record(pid, apply_slices(exp, sl)))
            tmp.close()
            if out_path and not direct:
                if existing:
                    shutil.copymode(real_out, tmp.name)  # CHANGED(R2-E5): keep the existing file's mode
                else:
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
            for f in (tmp, stage, *scratch):
                if f is not None:
                    f.close()
                    try:
                        os.unlink(f.name)
                    except OSError:
                        pass
    finally:
        dna.close()
