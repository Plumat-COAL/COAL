"""Gradio-free logic behind the app's Ka/Ks and dN/dS pipeline features (the GUI's "Ka/Ks" and "dN/dS pipeline"
tabs and the MCP tools ``dnds_pairs`` / ``pipeline_pairs``).

As for codon alignment (core.py), the engine runs as ``python -m coal`` in a subprocess: ``coal dnds`` for Ka/Ks on
alignments the person already has, ``coal pipeline pairs`` for two files of coding sequences. Both need COAL's C++
engine coal-cpp (``coal --cpp-path``); the pipeline also needs DIAMOND, found by COAL's own lookup
(``coal.dnds_pipeline.find_diamond``: ``$COAL_DIAMOND``, ``PATH``, COAL's cache). Neither the GUI nor the tools
install anything: they say how (``coal pipeline install diamond``).

Limits: every input has a size limit (smaller on the hosted Space), checked before the engine starts, so a run
that would take hours or fill the disk is refused with a sentence that names the limit and the command line to use
instead. Results go to a fresh run folder (core.new_run_dir, deleted after core.RUNS_MAX_AGE_HOURS); the pipeline's
work folder (intermediate files) is removed when the run ends.
"""
from __future__ import annotations

import csv
import html
import os
import shutil
import subprocess
import sys
import threading
import time
import zipfile
from dataclasses import dataclass, field
from typing import Callable, Iterator

from . import core

# The Ka/Ks methods, in the engine's row order (coal.dnds_pipeline.METHODS; a test keeps them equal). YN is the
# default, as in `coal dnds`.
METHODS = ("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN", "Comeron", "Li")
DEFAULT_METHODS = ("YN",)
CODE1_ONLY = ("Comeron", "Li")   # implemented for the standard genetic code only (docs/kaks.md)
METHOD_NAMES = {
    "YN": "YN: Yang & Nielsen 2000 (yn00)",
    "NG": "NG: Nei & Gojobori 1986",
    "LWL": "LWL: Li, Wu & Luo 1985",
    "MLWL": "MLWL: modified LWL (Tzeng et al. 2004)",
    "LPB": "LPB: Li 1993; Pamilo & Bianchi 1993",
    "MLPB": "MLPB: modified LPB (Tzeng et al. 2004)",
    "MYN": "MYN: modified YN (Zhang et al. 2006)",
    "Comeron": "Comeron 1995 (orthologr's default)",
    "Li": "Li 1993 as seqinr::kaks",
}

# The statuses of a result row, most serious last (docs/kaks.md "Policy"), with what each means for the person.
STATUS_ORDER = ("ok", "warning", "not_converged", "saturated", "no_data", "error")
STATUS_TEXT = {
    "ok": "Estimated. Gaps, ambiguous bases and a stop codon at the very end were skipped, as usual.",
    "warning": "Estimated, but read the note: internal stop codons were skipped, or a filter you set (fewer than "
               "the minimum codons compared, dS above the maximum, low protein identity) flags the pair.",
    "not_converged": "YN or MYN only: the iteration did not settle (short or nearly saturated pairs), so there is "
                     "no estimate; dN and dS are NA. Try NG or LWL for this pair, or leave it out.",
    "saturated": "Too many differences for the correction: dN or dS (see the note) is NA because the "
                 "multiple-hit correction is undefined. The sequences are too divergent for this method.",
    "no_data": "Nothing to compare: no codon column was left, or the pair has no synonymous (or no nonsynonymous) "
               "site, so that estimate and dN/dS are NA.",
    "error": "The pair could not be read: the two rows differ in length or are not a multiple of 3 long (the "
             "input is probably not a codon alignment).",
}
STATUS_ICON = {"ok": ("✔", "ok"), "warning": ("⚠", "warning"), "not_converged": ("–", "unchecked"),
               "saturated": ("⚠", "warning"), "no_data": ("–", "unchecked"), "error": ("✖", "failed")}

# Shown columns of the result tables (the files have all of them)
SHOWN_KAKS = ("seq1", "seq2", "method", "dN", "dS", "dNdS", "codons_compared", "status", "note")
SHOWN_PIPELINE = ("seq1", "seq2", "method", "dN", "dS", "dNdS", "prot_identity", "status", "note")

# ---------------------------------------------------------------- limits

MB = 1_000_000
KAKS_MAX_BYTES = 100 * MB            # one alignment file (FASTA or AXT)
KAKS_MAX_BYTES_HOSTED = 5 * MB
KAKS_MAX_PASTE = 5 * MB              # characters of pasted text
KAKS_MAX_RESULTS = 500_000           # pairs x methods
KAKS_MAX_RESULTS_HOSTED = 20_000
KAKS_TIMEOUT = 3600                  # seconds
PIPELINE_MAX_BYTES = 1000 * MB       # each coding-sequence file (uploads; gzip counts compressed)
PIPELINE_MAX_BYTES_MCP = 500 * MB    # the MCP tool: an agent waits for the answer
PIPELINE_TIMEOUT = 12 * 3600         # GUI (it has a Stop button)
PIPELINE_TIMEOUT_MCP = 3600
TABLE_ROWS_SHOWN = 300               # rows drawn in the GUI's table (the download has all)


def kaks_limits() -> tuple[int, int]:
    """(largest file in bytes, most results) for Ka/Ks here: smaller on the hosted Space."""
    if core.hosted():
        return KAKS_MAX_BYTES_HOSTED, KAKS_MAX_RESULTS_HOSTED
    return KAKS_MAX_BYTES, KAKS_MAX_RESULTS


# ---------------------------------------------------------------- tools on this machine

def cpp_path(timeout: float = 60) -> str | None:
    """The coal-cpp that `coal dnds` / `coal pipeline` would run (`coal --cpp-path`: the lookup and the version
    check of the engine itself), or None when this install has the Python engine only."""
    try:
        r = subprocess.run(core._coal_cmd() + ["--cpp-path"], capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=timeout, stdin=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError):
        return None
    out = (r.stdout or "").strip().splitlines()
    return out[-1].strip() if r.returncode == 0 and out else None


def diamond_path() -> str | None:
    """The DIAMOND the pipeline would use (COAL's lookup: $COAL_DIAMOND, PATH, COAL's cache), or None."""
    try:
        from coal.dnds_pipeline import find_diamond
    except ImportError:  # pragma: no cover - an engine without the pipeline
        return None
    try:
        return find_diamond()
    except OSError:  # pragma: no cover
        return None


def diamond_version(exe: str) -> str:
    """'2.1.10' from `diamond version`; '' when it does not run."""
    try:
        r = subprocess.run([exe, "version"], capture_output=True, text=True, encoding="utf-8", errors="replace",
                           timeout=60, stdin=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError):
        return ""
    words = (r.stdout or r.stderr or "").split()
    return words[-1] if r.returncode == 0 and words else ""


INSTALL_DIAMOND = "coal pipeline install diamond"


def install_command() -> str:
    """How to install DIAMOND into COAL's cache from a terminal: `coal pipeline install diamond`; in the Windows app
    (a frozen COAL.exe, no `coal` command on PATH) the same through COAL.exe, which runs the engine with `-m coal`."""
    if getattr(sys, "frozen", False):
        return f'"{sys.executable}" -m coal pipeline install diamond'
    return INSTALL_DIAMOND


NEEDS_CPP = ("Ka/Ks and the dN/dS pipeline need COAL's C++ engine (coal-cpp), which this install does not have "
             "(it has the Python engine only, which makes codon alignments). The platform wheels and the Windows "
             "app include it: install coal_align-<version>-py3-none-<platform>.whl from "
             "https://github.com/plumat-coal/COAL/releases/latest, or build it (docs/cpp-engine.md) and set COAL_CPP. "
             "`coal --cpp-path` shows where COAL looks.")
NEEDS_DIAMOND = ("The pipeline searches for orthologs with DIAMOND, which was not found (COAL looks at "
                 "$COAL_DIAMOND, then PATH, then its own cache). Install it into COAL's cache with "
                 f"`{install_command()}` (or put DIAMOND 2.x on PATH), then check again.")


@dataclass
class ToolStatus:
    cpp: str | None
    diamond: str | None
    diamond_version: str = ""

    @property
    def kaks_ready(self) -> bool:
        return self.cpp is not None

    @property
    def pipeline_ready(self) -> bool:
        return self.cpp is not None and self.diamond is not None


def tool_status(need_diamond: bool = True) -> ToolStatus:
    cpp = cpp_path()
    dia = diamond_path() if need_diamond else None
    return ToolStatus(cpp, dia, diamond_version(dia) if dia else "")


# ---------------------------------------------------------------- checking inputs

def validate_methods(methods, codontable: int) -> list[str]:
    """The chosen methods in the engine's order; ValueError with one sentence when none, an unknown one, or
    Comeron / Li with a genetic code other than 1."""
    if isinstance(methods, str):
        methods = [methods]
    chosen = []
    by_lower = {m.lower(): m for m in METHODS}
    by_lower["c95"] = "Comeron"
    for m in methods or []:
        name = by_lower.get(str(m).strip().lower())
        if name is None:
            raise ValueError(f"unknown Ka/Ks method {m!r}: choose from {', '.join(METHODS)}")
        if name not in chosen:
            chosen.append(name)
    if not chosen:
        raise ValueError("choose at least one Ka/Ks method (YN is the usual choice)")
    if int(codontable) not in core.GENETIC_CODES:
        raise ValueError(f"unknown genetic code {codontable} (supported: "
                         f"{', '.join(str(c) for c in core.GENETIC_CODES)})")
    bad = [m for m in chosen if m in CODE1_ONLY]
    if bad and int(codontable) != 1:
        raise ValueError(f"{' and '.join(bad)} {'is' if len(bad) == 1 else 'are'} defined for the standard genetic "
                         f"code (1) only; choose another method for code {int(codontable)}")
    return [m for m in METHODS if m in chosen]


def sniff_format(path: str) -> str:
    """'fasta' when the first non-blank character is '>', else 'axt' (the engine's own rule)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s:
                return "fasta" if s.startswith(">") else "axt"
    return "empty"


def count_records(path: str) -> tuple[str, int]:
    """(format, records): FASTA sequences, or AXT records (blocks of a name line and two sequence lines)."""
    fmt = sniff_format(path)
    n = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        if fmt == "fasta":
            n = sum(1 for line in fh if line.startswith(">"))
        elif fmt == "axt":
            in_block = False
            for line in fh:
                s = line.strip()
                if not s:
                    in_block = False
                elif not in_block and not s.startswith("#"):
                    in_block = True
                    n += 1
    return fmt, n


def n_pairs(fmt: str, records: int, pairs: str) -> int:
    if fmt == "axt":
        return records
    if pairs == "consecutive":
        return records // 2
    return records * (records - 1) // 2


def check_kaks_input(path: str, pairs: str, n_methods: int) -> tuple[str, int]:
    """(format, number of pairs) of a Ka/Ks input, or ValueError with a plain sentence: too large, empty, too many
    pairs, an odd number of sequences with consecutive pairs."""
    max_bytes, max_results = kaks_limits()
    size = os.path.getsize(path)
    if size > max_bytes:
        raise ValueError(f"the file is {core.human_size(size)}; this app takes Ka/Ks inputs up to "
                         f"{core.human_size(max_bytes)}. For larger ones use the command line: "
                         "coal dnds IN -o dnds.tsv")
    fmt, recs = count_records(path)
    if fmt == "empty" or recs == 0:
        raise ValueError("the input is empty: give an aligned codon FASTA or an AXT file")
    if fmt == "fasta" and recs < 2:
        raise ValueError("the FASTA has one sequence; Ka/Ks compares pairs, so it needs at least two")
    np_ = n_pairs(fmt, recs, pairs)
    if np_ * max(1, n_methods) > max_results:
        how = ("every pair of " + f"{recs:,} sequences") if (fmt == "fasta" and pairs == "all") else f"{np_:,} pairs"
        raise ValueError(f"{how} with {n_methods} method(s) is {np_ * n_methods:,} results; this app computes up to "
                         f"{max_results:,} at a time. Use consecutive pairs, fewer methods, or the command line: "
                         "coal dnds IN -o dnds.tsv")
    return fmt, np_


def check_cds_input(path: str, what: str, max_bytes: int = PIPELINE_MAX_BYTES) -> None:
    size = os.path.getsize(path)
    if size > max_bytes:
        raise ValueError(f"the {what} file is {core.human_size(size)}; this app takes files up to "
                         f"{core.human_size(max_bytes)}. For larger ones use the command line: coal pipeline pairs")
    if size == 0:
        raise ValueError(f"the {what} file is empty")


def write_pasted(text: str, folder: str) -> str:
    """Pasted alignment text into folder/pasted.fasta or pasted.axt (by its first character)."""
    if len(text) > KAKS_MAX_PASTE:
        raise ValueError(f"the pasted text has {len(text):,} characters; paste up to {KAKS_MAX_PASTE:,} or upload "
                         "the file instead")
    first = text.lstrip()[:1]
    name = "pasted.fasta" if first == ">" else "pasted.axt"
    path = os.path.join(folder, name)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")
    return path


# ---------------------------------------------------------------- reading results

def read_table(path: str | None) -> list[dict]:
    if not path or not os.path.isfile(path):
        return []
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def status_counts(rows: list[dict]) -> dict:
    counts = {k: 0 for k in STATUS_ORDER}
    for r in rows:
        st = r.get("status") or "error"
        counts[st] = counts.get(st, 0) + 1
    counts["total"] = len(rows)
    return counts


def status_help(counts: dict, all_statuses: bool = False) -> list[dict]:
    """[{status, count, meaning}] for the statuses that occur (or all of them)."""
    return [{"status": s, "count": int(counts.get(s, 0)), "meaning": STATUS_TEXT[s]}
            for s in STATUS_ORDER if all_statuses or counts.get(s)]


def engine_lines(stderr: str, limit: int = 8) -> list[str]:
    return [ln.strip() for ln in (stderr or "").splitlines() if ln.strip()][:limit]


def cite(methods, codon_alignment: bool, diamond: bool = False) -> tuple[list[str], str]:
    """(keys, references text) for a run: coal.citations, as `coal cite` prints them."""
    try:
        from coal import citations
    except ImportError:  # pragma: no cover
        return [], ""
    keys = citations.keys_for(methods=tuple(methods), codon_alignment=codon_alignment, diamond=diamond)
    try:
        return keys, citations.format_references(keys)
    except KeyError:  # pragma: no cover
        return keys, ""


# ---------------------------------------------------------------- Ka/Ks on alignments (`coal dnds`)

@dataclass
class KaksResult:
    ok: bool                        # the engine finished and wrote the table
    exit_code: int
    tsv_path: str | None
    rows: list[dict]
    counts: dict
    fmt: str
    pairs: str
    methods: list[str]
    codontable: int
    stderr: str
    command: str
    error: str = ""                 # one sentence when the run did not finish


def kaks_command(in_path: str, out_path: str, methods, codontable: int = 1, pairs: str = "all",
                 min_codons: int = 0, max_ds: float | None = None, fmt: str | None = None) -> list[str]:
    cmd = core._coal_cmd() + ["dnds", in_path, "-o", out_path]
    for m in methods:
        cmd += ["-m", m]
    cmd += ["-c", str(int(codontable))]
    if fmt != "axt":
        cmd += ["--pairs", pairs]
    if min_codons:
        cmd += ["--min-codons", str(int(min_codons))]
    if max_ds:
        cmd += ["--max-ds", f"{float(max_ds):g}"]
    if core.hosted():
        cmd += ["-j", "2"]  # a shared machine
    return cmd


def shown_command(cmd: list[str], paths=()) -> str:
    """The command as the person would type it: `coal ...`, file names only (no server folders)."""
    parts = list(cmd)
    base = core._coal_cmd()
    if parts[:len(base)] == base:
        parts = ["coal"] + parts[len(base):]
    out = " ".join(core._q(p) for p in parts)  # POSIX shell quoting, or Windows double quotes
    return core.strip_folders(out, [p for p in paths if p])


def run_kaks(in_path: str, out_dir: str, methods=DEFAULT_METHODS, codontable: int = 1, pairs: str = "all",
             min_codons: int = 0, max_ds: float | None = None, timeout: float = KAKS_TIMEOUT,
             hide_paths: bool = True) -> KaksResult:
    """`coal dnds` on one aligned codon FASTA or AXT file, into out_dir/<name>.dnds.tsv. ValueError (one sentence)
    for an input refused before the engine starts: unknown method, too large, empty, coal-cpp missing.
    hide_paths: the command and messages name files without their folders (the GUI: no server folders on the
    page); the MCP tool keeps the full paths, so its command can be run as it is."""
    methods = validate_methods(methods, codontable)
    if pairs not in ("all", "consecutive"):
        raise ValueError(f"pairs must be 'all' or 'consecutive', not {pairs!r}")
    if min_codons is not None and int(min_codons) < 0:
        raise ValueError("the minimum number of codons cannot be negative")
    if max_ds is not None and float(max_ds) < 0:
        raise ValueError("the maximum dS cannot be negative")
    fmt, _np = check_kaks_input(in_path, pairs, len(methods))
    if cpp_path() is None:
        raise ValueError(NEEDS_CPP)
    stem = os.path.splitext(os.path.basename(in_path))[0] or "pairs"
    out = os.path.join(out_dir, f"{stem}.dnds.tsv")
    cmd = kaks_command(os.path.abspath(in_path), out, methods, codontable, pairs, int(min_codons or 0),
                       max_ds or None, fmt)
    hide = [os.path.abspath(in_path), out, out_dir] if hide_paths else []
    shown = shown_command(cmd, hide)
    try:
        code, _out, err = core.run_engine_command(cmd, timeout, temp_dir=out_dir)
    except subprocess.TimeoutExpired:
        return KaksResult(False, 1, None, [], status_counts([]), fmt, pairs, methods, int(codontable), "", shown,
                          f"Ka/Ks did not finish within {timeout:g} seconds; use the command line for this input")
    err = core.public_error(err, hide)
    rows = read_table(out) if os.path.isfile(out) else []
    finished = code in (0, 2) and os.path.isfile(out)
    error = ""
    if not finished:
        error = next((ln.split("error:", 1)[-1].strip() for ln in err.splitlines() if "error" in ln.lower()),
                     "") or (engine_lines(err, 1) or [f"coal dnds stopped with exit status {code}"])[0]
    elif code == 2:
        # finished, but no pair has an estimate (every result error or no_data): the engine's last line says why
        error = (engine_lines(err)[-1:] or ["no pair has an estimate"])[0]
    return KaksResult(finished, code, out if finished else None, rows, status_counts(rows), fmt, pairs, methods,
                      int(codontable), err, shown, error)


# ---------------------------------------------------------------- the pipeline (`coal pipeline pairs`)

# Stages of a pairs run as the person sees them, and the stderr line that starts each one
STAGES = ("Translating the coding sequences", "Searching for orthologs with DIAMOND",
          "Aligning the pairs and estimating Ka/Ks", "Writing the results")


def stage_of(line: str, current: int) -> int:
    """The stage a `coal pipeline:` progress line moves the run to (indexes into STAGES)."""
    t = line.split("coal pipeline:", 1)[-1].strip()
    if t.startswith("translated "):
        return max(current, 1)
    if t[:1].isdigit() and " pairs (" in t:
        return max(current, 2)
    if t.startswith("summary:"):
        return max(current, 3)
    return current


@dataclass
class PipelineRun:
    exit_code: int
    status: str                     # report.json's status (completed, failed, interrupted), or 'failed'
    report: dict
    outputs: dict                   # name -> path of a file that was written
    report_path: str | None
    zip_path: str | None
    rows: list[dict]
    counts: dict                    # the pair statuses of dnds.tsv
    stderr: str
    command: str
    error: str = ""
    log: list[str] = field(default_factory=list)


RESULT_FILES = ("dnds.tsv", "report.json", "failed.tsv", "genes.tsv", "orthologs.tsv", "unpaired.tsv",
                "dnds.kaks", "codon_alignments.fasta", "protein_alignments.fasta")


def pipeline_command(query: str, subject: str, out_dir: str, methods=DEFAULT_METHODS, codontable: int = 1,
                     longest_isoform: bool = False, ortho: str = "rbh", min_codons: int = 0,
                     max_ds: float | None = None, keep_alignments: bool = False, pairs_path: str | None = None,
                     diamond: str | None = None, threads: int = 0) -> list[str]:
    cmd = core._coal_cmd() + ["pipeline", "pairs", "--query", query, "--subject", subject, "-o", out_dir]
    for m in methods:
        cmd += ["-m", m]
    cmd += ["-c", str(int(codontable))]
    if longest_isoform:
        cmd.append("--longest-isoform")
    if pairs_path:
        cmd += ["--pairs", pairs_path]
    else:
        cmd += ["--ortho", ortho]
        if diamond:
            cmd += ["--diamond", diamond]
    if min_codons:
        cmd += ["--min-codons", str(int(min_codons))]
    if max_ds:
        cmd += ["--max-ds", f"{float(max_ds):g}"]
    if keep_alignments:
        cmd.append("--keep-alignments")
    if threads:
        cmd += ["-j", str(int(threads))]
    return cmd


def check_pipeline_args(methods, codontable, ortho, min_codons, max_ds) -> list[str]:
    methods = validate_methods(methods, codontable)
    if ortho not in ("rbh", "bh"):
        raise ValueError(f"orthologs must be 'rbh' (reciprocal best hits) or 'bh' (best hits), not {ortho!r}")
    if min_codons is not None and int(min_codons) < 0:
        raise ValueError("the minimum number of codons cannot be negative")
    if max_ds is not None and float(max_ds) < 0:
        raise ValueError("the maximum dS cannot be negative")
    return methods


def _start(cmd: list[str]) -> subprocess.Popen:
    """The pipeline with its progress lines (stderr) piped; stdout is not used. Its own process group on POSIX so
    that stopping it also stops DIAMOND and the engine (see core._start_engine)."""
    return subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                            text=True, encoding="utf-8", errors="replace", start_new_session=os.name != "nt")


def _stop(p: subprocess.Popen, grace: float = 10.0) -> None:
    """SIGTERM (the pipeline stops its tools and writes a report with status 'interrupted'), then the whole
    process group after `grace` seconds."""
    if p.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(p.pid)], capture_output=True, timeout=30,
                           stdin=subprocess.DEVNULL)
        except (OSError, subprocess.SubprocessError):
            pass
    else:
        p.terminate()
    try:
        p.wait(timeout=grace)
        return
    except subprocess.TimeoutExpired:
        pass
    core._kill_engine_tree(p)
    try:
        p.wait(timeout=grace)
    except subprocess.TimeoutExpired:  # pragma: no cover
        pass


_RUNNING: dict[str, subprocess.Popen] = {}   # key (the GUI's browser session) -> the running pipeline
_STOPPED: set[int] = set()                    # PIDs stopped on request (stop_run) or by the time limit
_LOCK = threading.Lock()


def stop_run(key: str) -> bool:
    """Stop the pipeline started under `key` (the GUI's Stop button). True if one was running."""
    with _LOCK:
        p = _RUNNING.get(key)
        if p is None or p.poll() is not None:
            return False
        _STOPPED.add(p.pid)
    _stop(p)
    return True


def _timed_out(p: subprocess.Popen) -> None:
    with _LOCK:
        _STOPPED.add(p.pid)
    _stop(p)


def iter_pipeline(cmd: list[str], out_dir: str, timeout: float = PIPELINE_TIMEOUT,
                  hide=(), key: str | None = None) -> Iterator[tuple[str, object]]:
    """Run the pipeline, yielding ("line", text) for each progress line as it comes and finally ("done",
    PipelineRun). stop_run(key) (the GUI's Stop button) or the time limit stops the pipeline and its tools; the
    run then ends as 'interrupted'. Closing the generator early (a cancelled request) also stops them. The work
    folder (intermediate files) is removed in every case. `hide`: paths shown as file names only."""
    hide = [h for h in hide if h]
    p = _start(cmd)
    if key is not None:
        with _LOCK:
            old = _RUNNING.get(key)
            _RUNNING[key] = p
        if old is not None and old.poll() is None:  # pragma: no cover - one run per session (concurrency_limit)
            _stop(old)
    timer = threading.Timer(timeout, _timed_out, (p,))
    timer.daemon = True
    timer.start()
    t0 = time.monotonic()
    log: list[str] = []
    try:
        assert p.stderr is not None
        for line in p.stderr:
            line = core.public_error(line.rstrip("\n"), hide)
            if line.strip():
                log.append(line)
                yield "line", line
        p.wait()
    finally:
        timer.cancel()
        _stop(p)  # nothing if it has ended; otherwise (generator closed) the pipeline and its tools are stopped
        if p.stderr is not None:
            p.stderr.close()
        if key is not None:
            with _LOCK:
                if _RUNNING.get(key) is p:
                    del _RUNNING[key]
        shutil.rmtree(os.path.join(out_dir, "work"), ignore_errors=True)  # intermediate files only
    with _LOCK:
        stopped = p.pid in _STOPPED
        _STOPPED.discard(p.pid)
    reason = ""
    if stopped:
        waited = time.monotonic() - t0
        reason = (f"the pipeline did not finish within {timeout:g} seconds and was stopped" if waited >= timeout
                  else "stopped: you pressed Stop")
    yield "done", finish_pipeline(out_dir, p.returncode, log, shown_command(cmd, hide), reason, hide)


def finish_pipeline(out_dir: str, exit_code: int, log: list[str], command: str, stopped: str = "",
                    hide=()) -> PipelineRun:
    """After the pipeline: read report.json and dnds.tsv, remove the work folder, zip the result files.
    `stopped`: why the run was stopped (Stop, the time limit), or ""."""
    import json
    shutil.rmtree(os.path.join(out_dir, "work"), ignore_errors=True)  # intermediate files only
    rpath = os.path.join(out_dir, "report.json")
    report: dict = {}
    try:
        with open(rpath, encoding="utf-8") as fh:
            report = json.load(fh)
    except (OSError, ValueError):
        rpath = None
    outputs = {}
    for k, v in (report.get("outputs") or {}).items():  # names relative to the output folder
        if isinstance(v, str) and os.path.isfile(os.path.join(out_dir, v)):
            outputs[k] = os.path.join(out_dir, v)
    if not outputs:  # no report (the run stopped early): what is there
        outputs = {os.path.splitext(n)[0]: os.path.join(out_dir, n) for n in RESULT_FILES
                   if os.path.isfile(os.path.join(out_dir, n)) and n != "report.json"}
    rows = read_table(outputs.get("dnds"))
    status = report.get("status") or ""
    error = ""
    if stopped:
        status, error = "interrupted", stopped
    elif report.get("error"):
        error = core.public_error(str(report["error"]), hide)
    elif exit_code not in (0, 2):
        error = next((ln.split("error:", 1)[-1].strip() for ln in log if "error:" in ln), "") or \
            f"the pipeline stopped with exit status {exit_code}"
    if status not in ("completed", "failed", "interrupted"):  # no report, or one the run could not finish
        status = "completed" if exit_code in (0, 2) and rows else "failed"
    zpath = None
    files = [os.path.join(out_dir, n) for n in RESULT_FILES if os.path.isfile(os.path.join(out_dir, n))]
    if files:
        zpath = os.path.join(out_dir, "coal_dnds_results.zip")
        try:
            with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
                for f in files:
                    z.write(f, os.path.basename(f))
        except OSError:
            zpath = None
    return PipelineRun(exit_code, status, report, outputs, rpath, zpath, rows, status_counts(rows),
                       "\n".join(log), command, error, log)


def run_pipeline(cmd: list[str], out_dir: str, timeout: float = PIPELINE_TIMEOUT,
                 on_line: Callable[[str], None] | None = None, hide=(), key: str | None = None) -> PipelineRun:
    result = None
    for kind, value in iter_pipeline(cmd, out_dir, timeout, hide, key):
        if kind == "line" and on_line:
            on_line(value)
        elif kind == "done":
            result = value
    assert result is not None
    return result


# ---------------------------------------------------------------- HTML for the GUI

def status_badge(status: str) -> str:
    icon, cls = STATUS_ICON.get(status, ("?", "unchecked"))
    return (f'<span class="coal-badge coal-{cls}" aria-hidden="true">{icon}</span> '
            f'{html.escape(status or "")}')


def table_html(rows: list[dict], columns=SHOWN_KAKS, limit: int = TABLE_ROWS_SHOWN, caption: str = "") -> str:
    """The first `limit` result rows as an HTML table (status as icon + word; values escaped)."""
    if not rows:
        return '<p class="coal-small">No result rows.</p>'
    head = "".join(f"<th scope=\"col\">{html.escape(c)}</th>" for c in columns)
    body = []
    for r in rows[:limit]:
        st = r.get("status") or ""
        cls = STATUS_ICON.get(st, ("", "unchecked"))[1]
        cells = []
        for c in columns:
            v = r.get(c)
            v = "" if v is None else str(v)
            cells.append(f"<td>{status_badge(v)}</td>" if c == "status" else f"<td>{html.escape(v)}</td>")
        body.append(f'<tr class="coal-{cls}">{"".join(cells)}</tr>')
    more = (f'<p class="coal-small">The first {limit:,} of {len(rows):,} rows; the download has them all.</p>'
            if len(rows) > limit else "")
    cap = f"<caption>{html.escape(caption)}</caption>" if caption else ""
    return (f'<div class="coal-dnds-scroll"><table class="coal-health coal-dnds-table">{cap}<thead><tr>{head}</tr>'
            f'</thead><tbody>{"".join(body)}</tbody></table></div>{more}')


def counts_chips(counts: dict) -> str:
    out = []
    for s in STATUS_ORDER:
        n = int(counts.get(s) or 0)
        if n:
            icon, cls = STATUS_ICON[s]
            out.append(f'<span class="coal-chip coal-{cls}"><span aria-hidden="true">{icon}</span> {n:,} '
                       f'{html.escape(s)}</span>')
    return f'<div class="coal-chips">{"".join(out)}</div>'


def status_help_md(counts: dict | None = None) -> str:
    """Markdown: what each status means (only those in `counts` when given, else all)."""
    items = status_help(counts or {}, all_statuses=counts is None)
    return "\n".join(f"- **{i['status']}**" + (f" ({i['count']:,})" if counts is not None else "")
                     + f": {i['meaning']}" for i in items)


def kaks_banner(r: KaksResult, source: str = "", notes=()) -> str:
    """The Ka/Ks result banner: icon + headline, a chip per status, the error or the engine's warnings (`notes`),
    what was computed and from which input."""
    c = r.counts
    est = sum(int(c.get(s) or 0) for s in ("ok", "warning"))
    if not r.ok:
        kind, icon, title = "failed", "✖", "Ka/Ks did not run"
    elif r.exit_code == 2 or est == 0:
        kind, icon, title = "failed", "✖", "No pair has an estimate"
    elif est == c.get("total"):
        kind, icon, title = ("ok", "✔", f"All {c['total']:,} results estimated" if c.get("total") != 1
                             else "The result is estimated")
    else:
        kind, icon, title = "warning", "⚠", f"{est:,} of {c['total']:,} results estimated"
    pairs = len({(row.get("seq1"), row.get("seq2")) for row in r.rows})
    detail = (f"{pairs:,} pair(s) × {len(r.methods)} method(s) ({', '.join(r.methods)}), genetic code "
              f"{r.codontable}; input read as {'AXT' if r.fmt == 'axt' else 'aligned FASTA'}"
              + (f", {r.pairs} pairs" if r.fmt == "fasta" else "") + ".")
    err = f'<p class="coal-detail">{html.escape(r.error)}</p>' if r.error else ""
    return (f'<div class="coal-banner coal-{kind}" role="status"><p class="coal-headline"><span class="coal-icon" '
            f'aria-hidden="true">{icon}</span> {html.escape(title)}</p>{counts_chips(c)}{err}'
            f'<p class="coal-detail">{html.escape(detail)}</p>'
            + "".join(f'<p class="coal-small">{html.escape(n)}</p>' for n in notes if n != r.error)
            + (f'<p class="coal-files">{html.escape(source)}</p>' if source else "") + "</div>")


def pipeline_banner(r: PipelineRun) -> str:
    c = r.counts
    rep_counts = r.report.get("counts") or {}
    pairs = (rep_counts.get("pairs") or {})
    genes = rep_counts.get("genes") or rep_counts.get("records") or {}
    est = sum(int(c.get(s) or 0) for s in ("ok", "warning"))
    if r.status == "interrupted":
        kind, icon, title = "failed", "✖", "The pipeline was stopped"
    elif r.status != "completed" or r.exit_code not in (0, 2):
        kind, icon, title = "failed", "✖", "The pipeline did not finish"
    elif not c.get("total"):
        kind, icon, title = "warning", "⚠", "The pipeline finished without any pair"
    elif est == c["total"] and r.exit_code == 0:
        kind, icon, title = "ok", "✔", f"dN/dS for {int(pairs.get('total') or 0) or c['total']:,} ortholog pairs"
    else:
        kind, icon, title = "warning", "⚠", f"{est:,} of {c['total']:,} results estimated"
    parts = []
    if genes.get("total") is not None:
        parts.append(f"{int(genes.get('total') or 0):,} records read ({int(genes.get('failed') or 0):,} not usable)")
    if pairs.get("total") is not None:
        parts.append(f"{int(pairs['total']):,} pairs")
    if r.report.get("failures_total"):
        parts.append(f"{int(r.report['failures_total']):,} failure(s), listed in failed.tsv")
    detail = "; ".join(parts) + ("." if parts else "")
    err = f'<p class="coal-detail">{html.escape(r.error)}</p>' if r.error else ""
    return (f'<div class="coal-banner coal-{kind}" role="status"><p class="coal-headline"><span class="coal-icon" '
            f'aria-hidden="true">{icon}</span> {html.escape(title)}</p>{counts_chips(c)}{err}'
            + (f'<p class="coal-detail">{html.escape(detail)}</p>' if detail else "") + "</div>")


def progress_html(stage: int, done: bool = False, elapsed: float | None = None, ended: bool = False) -> str:
    """The stages of a pipeline run as a short list: done ✔, current ▶ (✖ where a run that did not finish
    ended: `ended`), to come ○ (icon + words)."""
    items = []
    for i, name in enumerate(STAGES):
        if done or i < stage:
            mark, state = "✔", "done"
        elif i == stage:
            mark, state = ("✖", "not finished") if ended else ("▶", "running")
        else:
            mark, state = "○", "to come"
        items.append(f'<li class="coal-stage coal-stage-{state.replace(" ", "-")}"><span aria-hidden="true">{mark}'
                     f'</span> {html.escape(name)} <span class="coal-sr">({state})</span></li>')
    t = f'<p class="coal-small">{elapsed:,.0f} s</p>' if elapsed is not None else ""
    return f'<ol class="coal-stages" aria-live="polite">{"".join(items)}</ol>{t}'


def tools_md(st: ToolStatus) -> str:
    """The pipeline tab's line about the tools it needs."""
    lines = []
    if st.cpp:
        lines.append("✔ C++ engine (coal-cpp) found.")
    else:
        lines.append("✖ " + NEEDS_CPP)
    if st.diamond:
        v = f" {st.diamond_version}" if st.diamond_version else ""
        lines.append(f"✔ DIAMOND{v} found" + (f": `{st.diamond}`" if not core.hosted() else "") + ".")
    else:
        lines.append("✖ " + NEEDS_DIAMOND)
    return "\n\n".join(lines)
