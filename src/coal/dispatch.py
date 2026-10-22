"""Automatic use of the C++ engine (coal-cpp) for large inputs.

`coal` (the command, `python -m coal`, `COAL.exe -m coal`) hands a run to coal-cpp when
  * COAL_ENGINE is unset or "auto", every option given is one coal-cpp has (all of them today), every input is a
    regular file and the inputs hold at least CPP_MIN_BYTES bytes together (the measured crossover, with a margin
    for slower process start-up on Windows), and
  * a coal-cpp is found (find_cpp) whose `--version` is this package's version.
COAL_ENGINE=python never uses it; COAL_ENGINE=cpp always does (error if it is not found or does not run).
Output, stderr, the run report and the exit code are the same either way (tests/test_engine_cli.py). If coal-cpp
cannot be started, the run goes on in Python after one `coal: notice:` line. Only the command line dispatches:
`coal.cli.main(argv)` called from Python (the `pal2nal` drop-in, tests) always runs the Python engine, and so does
`coal --batch` (coal/batch.py): a batch and its jobs are never handed to coal-cpp.
"""
from __future__ import annotations

import os
import sys

from . import __version__

CPP_MIN_BYTES = 1_000_000  # protein alignment + DNA files, bytes
# argparse dests coal-cpp supports (cpp/src/main.cpp registers the same options)
CPP_OPTIONS = frozenset({"pep_aln", "nuc_fasta", "codontable", "nogap", "nomismatch", "output_format", "blockonly",
                         "out", "jobs", "quiet", "strict", "keep_going", "max_warnings", "progress", "resume",
                         "report", "json"})
NAMES = ("coal-cpp.exe", "coal-cpp") if os.name == "nt" else ("coal-cpp",)


class EngineError(Exception):
    """COAL_ENGINE asks for something that cannot be done (exit 1)."""


def _dirs() -> list[str]:
    """Folders next to the Python running COAL: a frozen app's folder and its _internal, a venv's bin/Scripts."""
    import sysconfig
    out = []
    exe = os.path.dirname(os.path.abspath(sys.executable)) if sys.executable else ""
    if exe:
        out += [exe, os.path.join(exe, "_internal")]
    try:
        scripts = sysconfig.get_path("scripts")
    except (KeyError, ValueError):
        scripts = None
    if scripts:
        out.append(scripts)
    seen, res = set(), []
    for d in out:
        key = os.path.normcase(os.path.abspath(d))
        if key not in seen:
            seen.add(key)
            res.append(d)
    return res


def find_cpp() -> str | None:
    """coal-cpp's path: COAL_CPP (that path only), else next to the Python running COAL, else on PATH."""
    import shutil
    env = os.environ.get("COAL_CPP")
    if env:
        return env if os.path.isfile(env) else None
    for d in _dirs():
        for n in NAMES:
            p = os.path.join(d, n)
            if os.path.isfile(p) and os.access(p, os.X_OK):
                return p
    for n in NAMES:
        p = shutil.which(n)
        if p:
            return p
    return None


def _input_bytes(args) -> int | None:
    """Total size of the input files, or None if one is not a regular file (a pipe: size unknown)."""
    total = 0
    for p in [args.pep_aln, *args.nuc_fasta]:
        try:
            st = os.stat(p)
        except OSError:
            return None  # let the Python engine give its usual error
        if not os.path.isfile(p):
            return None
        total += st.st_size
    return total


def _supported(args, defaults) -> bool:
    return all(k in CPP_OPTIONS or v == defaults.get(k) for k, v in vars(args).items())


def _version_ok(path: str) -> tuple[bool, str]:
    """(runs and is this version, why not)."""
    import subprocess
    try:
        r = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=30,
                           stdin=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError) as e:
        return False, f"could not start {path} ({getattr(e, 'strerror', None) or e})"
    first = (r.stdout or "").strip()
    if r.returncode != 0 or not first.startswith("coal "):
        return False, f"{path} --version failed (exit {r.returncode})"
    if first.split()[1:2] != [__version__]:
        return False, f"{path} is {first}, not coal {__version__}"
    return True, ""


def choose(args, argv, defaults, note=None) -> str | None:
    """The coal-cpp to run this command line with, or None for the Python engine. note(text): notices.
    Raises EngineError for an impossible COAL_ENGINE request."""
    mode = (os.environ.get("COAL_ENGINE") or "auto").strip().lower()
    if mode not in ("auto", "python", "cpp"):
        raise EngineError(f"COAL_ENGINE={os.environ.get('COAL_ENGINE')!r}: use auto (the default), python or cpp")
    if mode == "python":
        return None
    if not _supported(args, defaults):
        if mode == "cpp":
            raise EngineError("COAL_ENGINE=cpp: coal-cpp does not have every option of this command line")
        return None
    if mode == "auto":
        size = _input_bytes(args)
        by_size = size is not None and size >= CPP_MIN_BYTES
        if os.name == "nt" and not all(a.isascii() for a in argv):
            by_size = False  # coal-cpp reads its arguments in the ANSI code page on Windows
        # --resume: the engine that wrote the checkpoint resumes the run (each engine resumes only its own), so a
        # coal-cpp installed or removed between the run and the resume does not get in the way
        engine = None
        if getattr(args, "resume", False) and getattr(args, "out", None):
            from .checkpoint import checkpoint_engine
            engine = checkpoint_engine(args.out)
        if engine == "coal-python":
            if by_size and note:
                note(f"{args.out} was interrupted by the Python engine; resuming with it")
            return None
        if engine == "coal-cpp":
            path = find_cpp()
            if path is not None and _version_ok(path)[0]:
                if not by_size and note:
                    note(f"{args.out} was interrupted by the C++ engine (coal-cpp); resuming with it")
                return path
            return None  # the Python engine explains that coal-cpp must resume this run
        if not by_size:
            return None
    path = find_cpp()
    if path is None:
        if mode == "cpp":
            raise EngineError("COAL_ENGINE=cpp: coal-cpp was not found (set COAL_CPP to its path, or put it next to "
                              "Python or on PATH)")
        return None
    ok, why = _version_ok(path)
    if not ok:
        if mode == "cpp":
            raise EngineError(f"COAL_ENGINE=cpp: {why}")
        if note:
            note(f"{why}; using the Python engine")
        return None
    return path


def run_cpp(path: str, argv: list[str]) -> int | str:
    """Run coal-cpp with this command line in this process's place; its exit code, or why it could not start.
    POSIX: exec (signals, pipes and the exit code are coal-cpp's own); Windows: a child with the same handles."""
    import subprocess
    for s in (sys.stdout, sys.stderr):
        try:
            s.flush()
        except (OSError, ValueError, AttributeError):
            pass
    cmd = [path, *argv]
    if os.name != "nt":
        try:
            os.execv(path, cmd)
        except OSError as e:
            return f"could not start {path} ({e.strerror or e})"
    try:
        p = subprocess.Popen(cmd)
    except OSError as e:
        return f"could not start {path} ({e.strerror or e})"
    while True:
        try:
            return p.wait()
        except KeyboardInterrupt:  # Ctrl-C reaches coal-cpp too; it cleans up and exits 130
            continue
