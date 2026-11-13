#!/usr/bin/env python3
"""Vendor COAL's C++ core into the R package (packaging/r/coalr/src/coal).

The R package `coalr` compiles the same sources as `coal-cpp` (cpp/src), copied into the package so that
`R CMD INSTALL` and `remotes::install_github(subdir = "packaging/r/coalr")` need nothing outside the package
folder. This script makes the copy; the copies are committed.

    python scripts/sync_r_sources.py           # copy, write the manifest and coal_version.h
    python scripts/sync_r_sources.py --check   # exit 1 when a copy differs from cpp/src (CI, tests)

The copies have LF line ends on every system (.gitattributes: packaging/r/** eol=lf; an R source package has LF line
ends, and VENDORED.txt holds the sha256 of exactly those bytes). cpp/src is read with either line end, so a Windows
checkout with core.autocrlf compares equal. A copy that differs from what is wanted only by its line ends is
reported as that (--check), not as "out of date".

Only the in-memory core is vendored: tables, core, kaks (with kaks_comeron and kaks_li), align, dnds and pool.hpp. The file readers (seqio) and
the command-line front ends are not: the R package takes sequences from R objects and writes nothing.
"""
from __future__ import annotations

import hashlib
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "cpp" / "src"
DST = ROOT / "packaging" / "r" / "coalr" / "src" / "coal"
FILES = ["tables.hpp", "tables.cpp", "core.hpp", "core.cpp", "kaks.hpp", "kaks.cpp", "kaks_comeron.cpp", "kaks_li.cpp",
         "align.hpp", "align.cpp", "dnds.hpp", "dnds.cpp", "pool.hpp"]
HEADER = ("// Vendored from COAL cpp/src/{name} by scripts/sync_r_sources.py. Do not edit here: edit cpp/src and sync.\n"
          "// Licence: GPL-3.0-only, as COAL (https://github.com/Plumat-COAL/COAL); copyright and third-party notices\n"
          "// (KaKs_Calculator 2.0 and PAML for kaks.cpp / kaks.hpp, libsequence for kaks_comeron.cpp, seqinr for\n"
          "// kaks_li.cpp): inst/COPYRIGHTS of this package.\n")


def coal_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    return re.search(r'^version\s*=\s*"([^"]+)"', text, re.M).group(1)


def wanted() -> dict[str, bytes]:
    out = {}
    for name in FILES:
        body = (SRC / name).read_text(encoding="utf-8").replace("\r\n", "\n")
        out[name] = (HEADER.format(name=name) + body).encode("utf-8")
    version = coal_version()
    out["coal_version.h"] = (
        "// Written by scripts/sync_r_sources.py: the COAL version the vendored sources come from.\n"
        "#pragma once\n"
        f'#define COAL_VERSION_STRING "{version}"\n').encode()
    lines = [f"COAL {version}", "Licence: GPL-3.0-only (see inst/COPYRIGHTS for the third-party notices)",
             "sha256 of the vendored files (scripts/sync_r_sources.py --check verifies them)"]
    lines += [f"{hashlib.sha256(out[n]).hexdigest()}  {n}" for n in sorted(out)]
    out["VENDORED.txt"] = ("\n".join(lines) + "\n").encode()
    return out


def main(argv: list[str]) -> int:
    check = "--check" in argv
    files = wanted()
    stale = []
    crlf = []
    for name, data in files.items():
        path = DST / name
        have = path.read_bytes() if path.exists() else None
        if have != data:
            (crlf if have is not None and have.replace(b"\r\n", b"\n") == data else stale).append(name)
            if not check:
                DST.mkdir(parents=True, exist_ok=True)
                path.write_bytes(data)
    extra = sorted(p.name for p in DST.glob("*") if p.is_file() and p.name not in files and p.suffix in
                   (".cpp", ".hpp", ".h", ".txt")) if DST.exists() else []
    if check:
        if crlf and not (stale or extra):
            print("packaging/r/coalr/src/coal: the content is in sync with cpp/src, but these files have CRLF line "
                  "ends: " + ", ".join(crlf) + "\nAn R source package has LF line ends (.gitattributes: "
                  "packaging/r/** eol=lf; this checkout was made without that rule).\n"
                  "run: python scripts/sync_r_sources.py   (writes them with LF)", file=sys.stderr)
            return 1
        stale += crlf
        if stale or extra:
            print("packaging/r/coalr/src/coal is out of date: " + ", ".join(stale + [f"{e} (not vendored)" for e in extra])
                  + "\nrun: python scripts/sync_r_sources.py", file=sys.stderr)
            return 1
        print(f"packaging/r/coalr/src/coal: {len(files)} files in sync with cpp/src")
        return 0
    stale += crlf
    for e in extra:
        (DST / e).unlink()
    print(f"packaging/r/coalr/src/coal: {len(stale)} file(s) written, {len(files) - len(stale)} unchanged"
          + (f", {len(extra)} removed" if extra else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
