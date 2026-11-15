"""A stand-in for coal-cpp that does one thing first, then runs the real binary with the same arguments, stdin and
exit status: how a test sees whether (and how) coal handed a run to coal-cpp. A POSIX sh script, or on Windows a
.cmd batch file (CreateProcess runs those through cmd.exe), so the tests that use it run on all three systems.
The arguments the tests pass are plain (file names, options): cmd.exe's own metacharacters are not escaped."""
from __future__ import annotations

import os
from pathlib import Path


def write_wrapper(path: Path, target, *, touch: Path | None = None, unless_version: bool = False,
                  stderr_line: str | None = None) -> Path:
    """Write the stand-in at `path` (".cmd" is added on Windows) and return its path. `touch`: a file it creates
    first (not for `--version` when `unless_version`); `stderr_line`: a line it writes to stderr first."""
    path, target = Path(path), str(target)
    if os.name == "nt":
        path = path.with_name(path.name + ".cmd")
        lines = ["@echo off"]
        if stderr_line is not None:
            lines.append(f"echo {stderr_line} 1>&2")
        if touch is not None:
            make = f'type nul > "{touch}"'
            lines.append(f'if not "%~1"=="--version" {make}' if unless_version else make)
        lines += [f'"{target}" %*', "exit /b %ERRORLEVEL%"]
        path.write_text("\r\n".join(lines) + "\r\n", encoding="utf-8")
        return path
    lines = ["#!/bin/sh"]
    if stderr_line is not None:
        lines.append(f'echo "{stderr_line}" >&2')
    if touch is not None:
        lines.append(f'case "$1" in --version) ;; *) touch "{touch}" ;; esac' if unless_version else f'touch "{touch}"')
    lines.append(f'exec "{target}" "$@"')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    path.chmod(0o755)
    return path
