"""The environment for a test's subprocess that runs this checkout's Python package: PYTHONPATH with the absolute
path of src/ first and every other entry made absolute.

A subprocess often runs in another folder (cwd=examples/clean, a temporary folder), where a relative PYTHONPATH
(`PYTHONPATH=src python -m unittest ...`, as docs/development.md runs the tests from the root) names nothing, and
`python -m coal` fails with "No module named coal" (REVIEW-COAL-3 L4). Use python_env() instead of os.environ.
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def absolute_pythonpath(value: str | None, first: Path = SRC, base: str | None = None) -> str:
    """`first`, then the entries of `value` (os.pathsep-separated) made absolute against `base` (default: the
    current folder), without empty entries or a second copy of `first`."""
    base = base or os.getcwd()
    out = [str(first)]
    for entry in (value or "").split(os.pathsep):
        if not entry:
            continue
        p = os.path.normpath(os.path.join(base, entry))
        if p not in out:
            out.append(p)
    return os.pathsep.join(out)


def python_env(extra: dict | None = None, base: dict | None = None) -> dict:
    """A copy of `base` (default os.environ) with PYTHONPATH = absolute_pythonpath(its PYTHONPATH), then `extra`."""
    env = dict(os.environ if base is None else base)
    env["PYTHONPATH"] = absolute_pythonpath(env.get("PYTHONPATH"))
    env.update(extra or {})
    return env
