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
CHECKPOINT_EVERY = float(os.environ.get("COAL_CHECKPOINT_EVERY", "2.0"))  # seconds (tests: 0)
CHANGED = ("The input files or options changed since the interrupted run, so it can't be resumed safely. "
           "Start a new run without --resume.")


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

    def _signature(self):
        sig = []
        for p in self.inputs:
            st = os.stat(p)
            sig.append({"path": os.path.abspath(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns})
        return sig

    def load_for_resume(self, sink_path: str) -> dict | None:
        """The saved state if this run can resume, None if there is nothing to resume.
        Raises CoalError(code="resume_mismatch") if inputs or options changed."""
        import json
        if not os.path.exists(self.path):
            self.why = "no checkpoint"
            return None
        try:
            with open(self.path, encoding="utf-8") as fh:
                st = json.load(fh)
        except (OSError, ValueError):
            self.why = "an unreadable checkpoint"
            return None
        if (st.get("engine") != ENGINE or st.get("coal_version") != __version__
                or st.get("options_hash") != self.options_hash or st.get("inputs") != self._signature()):
            raise CoalError(CHANGED, code="resume_mismatch")
        self.why = "a checkpoint whose partial output is missing or incomplete"
        try:
            if os.path.getsize(sink_path) < st["sink_bytes"]:
                return None
            if st.get("entries_bytes") and os.path.getsize(self.entries) < st["entries_bytes"]:
                return None
        except (OSError, KeyError):
            return None
        return st

    def save(self, state: dict):
        import json
        self.count += 1
        st = dict(state, engine=ENGINE, coal_version=__version__, options_hash=self.options_hash,
                  inputs=self._signature())
        tmp = self.path + ".tmp"
        with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(st, fh)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, self.path)

    def remove(self, *, keep_partial: bool = False):
        for p in (self.path, self.path + ".tmp", self.entries, self.stage) + (() if keep_partial else (self.partial,)):
            try:
                os.unlink(p)
            except OSError:
                pass


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
