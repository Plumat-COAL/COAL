#!/usr/bin/env python3
"""Vendor COAL's C++ core into the R package (packaging/r/coalr/src/coal).

The R package `coalr` compiles the same sources as `coal-cpp` (cpp/src), copied into the package so that
`R CMD INSTALL` and `remotes::install_github(subdir = "packaging/r/coalr")` need nothing outside the package
folder. This script makes the copy; the copies are committed.

    python scripts/sync_r_sources.py           # copy, write the manifest and coal_version.h
    python scripts/sync_r_sources.py --check   # exit 1 when a copy differs from cpp/src (CI, tests)

Only the in-memory core is vendored: tables, core, kaks, align, dnds and pool.hpp. The file readers (seqio) and
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
FILES = ["tables.hpp", "tables.cpp", "core.hpp", "core.cpp", "kaks.hpp", "kaks.cpp", "align.hpp", "align.cpp",
         "dnds.hpp", "dnds.cpp", "pool.hpp"]
HEADER = "// Vendored from COAL cpp/src/{name} by scripts/sync_r_sources.py. Do not edit here: edit cpp/src and sync.\n"


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
    lines = [f"COAL {version}", "sha256 of the vendored files (scripts/sync_r_sources.py --check verifies them)"]
    lines += [f"{hashlib.sha256(out[n]).hexdigest()}  {n}" for n in sorted(out)]
    out["VENDORED.txt"] = ("\n".join(lines) + "\n").encode()
    return out


def main(argv: list[str]) -> int:
    check = "--check" in argv
    files = wanted()
    stale = []
    for name, data in files.items():
        path = DST / name
        if not path.exists() or path.read_bytes() != data:
            stale.append(name)
            if not check:
                DST.mkdir(parents=True, exist_ok=True)
                path.write_bytes(data)
    extra = sorted(p.name for p in DST.glob("*") if p.is_file() and p.name not in files and p.suffix in
                   (".cpp", ".hpp", ".h", ".txt")) if DST.exists() else []
    if check:
        if stale or extra:
            print("packaging/r/coalr/src/coal is out of date: " + ", ".join(stale + [f"{e} (not vendored)" for e in extra])
                  + "\nrun: python scripts/sync_r_sources.py", file=sys.stderr)
            return 1
        print(f"packaging/r/coalr/src/coal: {len(files)} files in sync with cpp/src")
        return 0
    for e in extra:
        (DST / e).unlink()
    print(f"packaging/r/coalr/src/coal: {len(stale)} file(s) written, {len(files) - len(stale)} unchanged"
          + (f", {len(extra)} removed" if extra else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
