"""A stand-in for coal-cpp that does one thing first, then runs the real binary with the same arguments, stdin and
exit status: how a test sees whether (and how) coal handed a run to coal-cpp. A POSIX sh script, or on Windows a
.cmd batch file (CreateProcess runs those through cmd.exe), so the tests that use it run on all three systems.

Limits (it is a test helper, not a launcher to ship):

* No escaping. The .cmd file passes its arguments on with `%*` and puts `target`, `touch` and `stderr_line` into
  the batch text as they are: cmd.exe's metacharacters (`&`, `|`, `<`, `>`, `^`, `%`, `!`, and `"` inside a
  word) in any of them change what the batch file does. The sh script quotes `"$@"`, but `target`, `touch` and
  `stderr_line` go inside double quotes, so `"`, `$`, `` ` `` and `\\` in them are the shell's. The tests pass plain
  file names and options only.
* On Windows, killing it does not kill the real coal-cpp. (On POSIX the script `exec`s coal-cpp, so the process is
  replaced, same PID, and a signal reaches coal-cpp.) On Windows cmd.exe stays the parent: TerminateProcess on the wrapper (Popen.kill,
  a subprocess timeout) ends cmd.exe only and leaves the real coal-cpp running (and holding its files) until it
  finishes; a test that kills a run must kill the process tree (taskkill /T) or not use the wrapper.
* On Windows the exit status is the one `exit /b %ERRORLEVEL%` passes on, and Ctrl+C prompts "Terminate batch
  job (Y/N)?" when a console is attached."""
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
