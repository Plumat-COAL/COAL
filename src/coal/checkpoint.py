"""Checkpoints for resumable runs (--resume).

Every run with -o writes its output to <out>.coal-partial (or, for formats
that need a second pass, its rows to <out>.coal-stage) and, at most every
CHECKPOINT_EVERY seconds, a small JSON checkpoint <out>.coal-checkpoint.json:
engine + version, a hash of the options, the size and mtime of every input,
how many sequences (in input order) are finished and how many bytes of the
partial file they occupy. Checkpoints are written atomically (temp file +
rename), after the data they describe has been flushed to disk, so a crash
never leaves a checkpoint that points past the end of the data.

--resume truncates the partial file to the recorded length, skips the finished
sequences and carries on; the result is byte-identical to an uninterrupted run.
"""
from __future__ import annotations

import os

from . import __version__
from .core import CoalError

ENGINE = "coal-python"
ENGINE_NAMES = ("coal-python", "coal-cpp")  # what the two engines write as "engine"
CHECKPOINT_EVERY = float(os.environ.get("COAL_CHECKPOINT_EVERY", "2.0"))  # seconds (tests: 0)
CHANGED = ("The input files or options changed since the interrupted run, so it can't be resumed safely. "
           "Start a new run without --resume.")
ENGINE_WORDS = {"coal-python": ("Python", "python"), "coal-cpp": ("C++", "cpp")}


def other_engine_message(out: str, engine: str) -> str:
    """The error when a checkpoint was written by the other engine (the same text in coal-cpp)."""
    word, env = ENGINE_WORDS[engine]
    return (f"{out} was interrupted by the {word} engine ({engine}), which must resume it: run again with "
            f"COAL_ENGINE={env}, or start a new run without --resume")


def checkpoint_engine(out_path) -> str | None:
    """Which engine wrote the checkpoint of an interrupted run to `out_path` ("coal-python" / "coal-cpp"), or None
    (no checkpoint, or not one of ours). Read only, cheap: for the command line to pick that engine on --resume."""
    import json
    path = os.path.abspath(out_path) + ".coal-checkpoint.json"
    try:
        if os.path.getsize(path) > 50_000_000:
            return None
        with open(path, encoding="utf-8") as fh:
            st = json.load(fh)
    except Exception:  # noqa: BLE001 - unreadable, whatever the reason: the engine will say so
        return None
    eng = st.get("engine") if isinstance(st, dict) else None
    return eng if eng in ENGINE_NAMES else None


class Checkpoint:
    def __init__(self, out_path: str, options: dict, inputs: list[str]):
        import hashlib, json
        base = os.path.abspath(out_path)
        self.partial = base + ".coal-partial"
        self.stage = base + ".coal-stage"
        self.entries = base + ".coal-entries.jsonl"
        self.path = base + ".coal-checkpoint.json"
        opts = dict(options, inputs=[os.path.abspath(p) for p in inputs], out=base)
        self.options_hash = hashlib.sha256(json.dumps(opts, sort_keys=True).encode()).hexdigest()
        self.inputs = inputs
        self.count = 0
        self.why = "no checkpoint"   # why load_for_resume() returned None
        self.other_engine = None     # the engine named by a checkpoint this engine cannot resume (resume_mismatch)
        self.name = out_path         # how messages name the output (the -o argument)

    def _signature(self):
        sig = []
        for p in self.inputs:
            st = os.stat(p)
            sig.append({"path": os.path.abspath(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns})
        return sig

    def load_for_resume(self, sink_path: str, nseq: int | None = None, ncol: int | None = None,
                        ntriplets: int | None = None) -> dict | None:
        """The saved state if this run can resume, None if there is nothing to resume (self.why says why).
        Raises CoalError(code="resume_mismatch") if inputs or options changed.
        Every field is checked (types, ranges against `nseq`, `ncol` (protein columns, bounding `drop`), `ntriplets`
        (output codon columns, bounding `bad`) and the partial files' sizes, the entries file's lines) and any
        problem counts as an unreadable checkpoint, so an edited or corrupted checkpoint can never make a resumed
        run write a truncated or wrong output; the entries are checked line by line and not kept.
        The C++ engine applies the same rules (main.cpp Checkpoint::load_for_resume)."""
        import json
        if not os.path.exists(self.path):
            self.why = "no checkpoint"
            return None
        self.why = "an unreadable checkpoint"
        try:
            with open(self.path, encoding="utf-8") as fh:
                st = json.load(fh)
        except Exception:  # noqa: BLE001 - OSError, ValueError, UnicodeDecodeError, RecursionError (deep nesting) ...
            return None
        if not isinstance(st, dict):
            return None
        if st.get("engine") != ENGINE and st.get("engine") in ENGINE_NAMES:
            # Each engine resumes only its own runs (their partial files differ): the error says which one must
            # (the command line normally picks it by itself: dispatch.py checkpoint_engine)
            self.other_engine = st["engine"]
            raise CoalError(other_engine_message(self.name, st["engine"]), code="resume_mismatch")
        if (st.get("engine") != ENGINE or st.get("coal_version") != __version__
                or st.get("options_hash") != self.options_hash or st.get("inputs") != self._signature()):
            raise CoalError(CHANGED, code="resume_mismatch")
        if not _well_formed(st, nseq, ncol, ntriplets):
            return None
        self.why = "a checkpoint whose partial output is missing or incomplete"
        try:
            if os.path.getsize(sink_path) < st["sink_bytes"]:
                return None
            if os.path.getsize(self.entries) < st["entries_bytes"]:
                return None
        except OSError:
            return None
        self.why = "an unreadable checkpoint"
        # The entries are checked one line at a time and not kept (RunReport.restore reads them again the same
        # way): a run of a million genes of unrelated DNA has an entries file of tens of gigabytes
        try:
            left = st["entries_bytes"]
            with open(self.entries, "rb") as fh:
                for line in fh:
                    if left <= 0:
                        break
                    take = line[:left]
                    left -= len(take)
                    if not take.strip():
                        continue
                    e = json.loads(take.decode("utf-8"))
                    if not isinstance(e, dict) or not _is_int(e.get("index")) or not 0 <= e["index"] < st["completed"]:
                        return None
        except (OSError, UnicodeDecodeError, ValueError):
            return None
        return st

    def save(self, state: dict):
        import json
        self.count += 1
        st = dict(state, engine=ENGINE, coal_version=__version__, options_hash=self.options_hash,
                  inputs=self._signature())
        tmp = self.path + ".tmp"
        try:
            with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
                json.dump(st, fh)
                fh.flush()
                os.fsync(fh.fileno())
            os.replace(tmp, self.path)
        except OSError:  # the caller reports it; no temp file is left behind
            try:
                os.unlink(tmp)
            except OSError:
                pass
            raise

    def remove(self, *, keep_partial: bool = False, keep_entries: bool = False):
        """Remove the checkpoint and its files. keep_entries: the entries file stays (the run report reads it back
        for --report and removes it itself, RunReport.close)."""
        for p in ((self.path, self.path + ".tmp", self.stage) + (() if keep_entries else (self.entries,))
                  + (() if keep_partial else (self.partial,))):
            try:
                os.unlink(p)
            except OSError:
                pass


def _is_int(x) -> bool:
    return isinstance(x, int) and not isinstance(x, bool)


def _int_list(x, below: int | None = None) -> bool:
    return isinstance(x, list) and all(_is_int(v) and v >= 0 and (below is None or v < below) for v in x)


def _well_formed(st: dict, nseq: int | None, ncol: int | None = None, ntriplets: int | None = None) -> bool:
    """Whether a checkpoint's fields have the types and ranges a run of this engine writes (nseq: the number of
    sequences of this run, ncol: its protein columns, ntriplets: its output codon columns; None skips that check)."""
    for k in ("completed", "nseq", "sink_bytes", "entries_bytes", "skipped"):
        if not _is_int(st.get(k)) or st[k] < 0:
            return False
    if nseq is not None and st["nseq"] != nseq:
        return False
    if st["completed"] > st["nseq"] or st["skipped"] > st["completed"]:
        return False
    if not (_int_list(st.get("skip_idx"), st["completed"]) and _int_list(st.get("drop"), ncol)
            and _int_list(st.get("bad"), ntriplets)):
        return False
    if len(st["skip_idx"]) != st["skipped"]:
        return False
    rec = st.get("rec")
    if rec is not None:
        if not isinstance(rec, dict) or not isinstance(rec.get("kinds"), str) or len(rec["kinds"]) != st["completed"]:
            return False
        if not (_int_list(rec.get("cols")) and _int_list(rec.get("real_cols"))):
            return False
    return True


def open_append(path: str, resume_bytes: int | None):
    """Text file for writing: fresh, or truncated to resume_bytes and appended to."""
    if resume_bytes is None:
        return open(path, "w", newline="\n", encoding="ascii")
    with open(path, "r+b") as fh:
        fh.truncate(resume_bytes)
    return open(path, "a", newline="\n", encoding="ascii")


def durable_size(fh) -> int:
    """Flush a file to disk and return its size in bytes."""
    fh.flush()
    os.fsync(fh.fileno())
    return os.fstat(fh.fileno()).st_size
