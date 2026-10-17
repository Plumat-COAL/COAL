#!/usr/bin/env python3
"""Fetch the PAL2NAL v14 reference distribution into PAL2NAL/ (git-ignored).

PAL2NAL is licensed under GPL-2.0 and COAL contains none of its code, so it is
never committed to or redistributed from this repository.
Every developer fetches it locally with this script.

Usage:
    python scripts/fetch_pal2nal.py                 # official site; if that fails, the orthologr GitHub copy
    python scripts/fetch_pal2nal.py --archive FILE  # use an already-downloaded tarball

Sources, in order (CHANGED(DEC-4): the GitHub fallback is new; see overnight/integration/CHANGES.md):
  1. the official tarball (checked against ARCHIVE_SHA256, then pal2nal.pl against SCRIPT_SHA256);
  2. the copy of pal2nal.pl in the public orthologr R package on GitHub (checked against SCRIPT_SHA256, the
     checksum of the script in the official tarball, so it is byte for byte the official v14 script).
Nothing is kept unless its checksum matches. Prints the script's path on success (for PAL2NAL=... in CI).
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

URL = "https://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz"
GITHUB_URL = "https://raw.githubusercontent.com/HajkD/orthologr/master/inst/pal2nal/pal2nal.v14/pal2nal.pl"
# SHA-256 of the tarball and of pal2nal.pl as originally downloaded (16 Oct 2026).
ARCHIVE_SHA256 = "0e458d38571ab6d61078047a9e80d843bb2e733414eeff8b1144d1c6bc13848f"
SCRIPT_SHA256 = "2b419417e304e9a63129056d61fa513a4c5834257a7b162cc4389c3ac53e23f2"

REPO_ROOT = Path(__file__).resolve().parent.parent
DEST = REPO_ROOT / "PAL2NAL"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def safe_extract(tar: tarfile.TarFile, dest: Path) -> None:
    """Extract, refusing absolute paths, '..' and links (archive is third-party)."""
    dest = dest.resolve()
    for m in tar.getmembers():
        target = (dest / m.name).resolve()
        if not str(target).startswith(str(dest)) or m.issym() or m.islnk():
            sys.exit(f"Refusing unsafe archive member: {m.name}")
    tar.extractall(dest)


def _download(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()


def _from_archive(data: bytes, script: Path) -> None:
    got = sha256(data)
    if got != ARCHIVE_SHA256:
        raise ValueError(f"checksum mismatch for the archive: expected {ARCHIVE_SHA256}, got {got}")
    DEST.mkdir(exist_ok=True)
    (DEST / "pal2nal.v14.tar.gz").write_bytes(data)
    with tarfile.open(DEST / "pal2nal.v14.tar.gz", "r:gz") as tar:
        safe_extract(tar, DEST)
    if sha256(script.read_bytes()) != SCRIPT_SHA256:
        script.unlink()
        raise ValueError(f"checksum mismatch for {script}")


def _from_github(script: Path) -> None:
    data = _download(GITHUB_URL)
    got = sha256(data)
    if got != SCRIPT_SHA256:
        raise ValueError(f"checksum mismatch for {GITHUB_URL}: expected {SCRIPT_SHA256}, got {got}")
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_bytes(data)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--archive", type=Path, help="local pal2nal.v14.tar.gz to use instead of downloading")
    ap.add_argument("--force", action="store_true", help="re-extract even if PAL2NAL/ already exists")
    args = ap.parse_args(argv)

    script = DEST / "pal2nal.v14" / "pal2nal.pl"
    if script.exists() and not args.force:
        ok = sha256(script.read_bytes()) == SCRIPT_SHA256
        print(f"{script} already present (checksum {'OK' if ok else 'MISMATCH'}).")
        return 0 if ok else 1

    if args.archive:
        try:
            _from_archive(args.archive.read_bytes(), script)
        except ValueError as e:
            sys.exit(str(e))
    else:
        problems = []
        for name, get in (("official site", lambda: _from_archive(_download(URL), script)),
                          ("orthologr copy on GitHub", lambda: _from_github(script))):
            print(f"Fetching PAL2NAL v14 from the {name}")
            try:
                get()
                break
            except (OSError, ValueError, tarfile.TarError) as e:
                problems.append(f"{name}: {e}")
                print(f"  failed: {e}", file=sys.stderr)
        else:
            sys.exit("Could not fetch PAL2NAL v14:\n  " + "\n  ".join(problems))
    print(f"PAL2NAL v14 ready: {script}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
