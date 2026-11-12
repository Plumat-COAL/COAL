"""Multiple protein alignments of clusters for the codon step: Clustal Omega (default) or MAFFT (docs/alignment.md).

    from coal import msa
    exe = msa.find_clustalo() or msa.install_clustalo()      # PATH, COAL's cache, else install into the cache
    rows = msa.run_msa([("g1", "MKV..."), ("g2", "MRV...")])   # [(id, gapped row)], input order
    results = msa.align_clusters({"c1": [...], "c2": [...]}, jobs=8)   # many clusters, one process each

Neither tool is bundled with COAL (docs/alignment.md "Licences": COAL is GPL-3.0-only and Clustal Omega
GPL-2.0-or-later, so bundling would be licence-compatible, but it is not done).
install_clustalo() puts one into COAL's cache folder (cache_dir()), like orthologr's diamond_install(): through
conda / mamba / micromamba from Bioconda (micromamba itself is one pinned release fetched from GitHub and kept only
if its SHA-256 is the pinned one, if no conda is installed), or the official static binary from clustal.org
(method="official"; the only way on Windows), which has no published checksum: that route needs sha256=... or
allow_unverified=True. No downloaded file is run before its SHA-256 was checked.

Stop codons: '*' at the end of a sequence is removed before aligning and given back as a last column of its own
('*' in the rows that had one, '-' in the others); an internal '*' is aligned as X and put back afterwards. (Clustal
Omega would otherwise align '*' as a residue.) The returned rows hold exactly the input residues.

Not part of the stable `import coal` API yet (python-api.md): it may change before 1.0.
"""
from __future__ import annotations

import hashlib
import os
import platform
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

__all__ = ["cache_dir", "find_clustalo", "find_mafft", "install_clustalo", "install_mafft", "tool_version",
           "run_msa", "align_clusters", "MsaResult", "MsaError", "CLUSTALO_VERSION", "MICROMAMBA_VERSION"]

CLUSTALO_VERSION = "1.2.4"
_EXE = ".exe" if os.name == "nt" else ""

# ---- Downloads that COAL runs afterwards: each is pinned to one file by its SHA-256 ------------------------------
# A downloaded file is moved into the cache only after its SHA-256 matched the pinned value (or the one the caller
# gave); a file with no known hash is refused unless the caller passes allow_unverified=True. Nothing fetched is
# ever executed before this check.
#
# Official Clustal Omega binaries (www.clustal.org/omega/). clustal.org publishes no checksums and serves plain
# HTTP (https is tried first, in case the host offers it), and the host could not be reached from the machine
# these tables were made on (2026-11-12), so no hash is pinned for these files yet: this route refuses to run
# unless the caller gives a hash (sha256=...) or accepts an unverified download (allow_unverified=True). To pin
# one: download the file from a network you trust, check it, and put its `sha256sum` here.
# (system, machine) -> (path under www.clustal.org, member of the zip or None, SHA-256 or None)
_OFFICIAL_HOST = "www.clustal.org"
_OFFICIAL = {
    ("Linux", "x86_64"): ("/omega/clustalo-1.2.4-Ubuntu-x86_64", None, None),
    ("Darwin", "x86_64"): ("/omega/clustal-omega-1.2.3-macosx", None, None),
    ("Darwin", "arm64"): ("/omega/clustal-omega-1.2.3-macosx", None, None),  # x86_64, Rosetta 2
    ("Windows", "AMD64"): ("/omega/clustal-omega-1.2.2-win64.zip", "clustalo.exe", None),
}
# micromamba (static, BSD-3-Clause) from https://github.com/mamba-org/micromamba-releases: one fixed release, not
# "latest", with the SHA-256 of each asset (computed from the downloads on 2026-11-12; each equals the
# <asset>.sha256 file the release publishes).
MICROMAMBA_VERSION = "2.10.0-0"
_MICROMAMBA = {
    ("Linux", "x86_64"): ("micromamba-linux-64",
                          "701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd"),
    ("Linux", "aarch64"): ("micromamba-linux-aarch64",
                           "3bb56786e21b9b24024f599fe9185b4a59dc3f9f59937a42311ee583adbbb23c"),
    ("Darwin", "x86_64"): ("micromamba-osx-64",
                           "8421462a1ef1df9f7acd114f9e2a24c1332ce4f4e29f44465b4acb44e2755172"),
    ("Darwin", "arm64"): ("micromamba-osx-arm64",
                          "c547639f5caf742236641c362746abb4a3505774753ecd08fe341643670c4c01"),
}
_MICROMAMBA_URL = "https://github.com/mamba-org/micromamba-releases/releases/download/{version}/{asset}"


class MsaError(RuntimeError):
    """The aligner is missing, failed, or the input cannot be aligned."""


def _platform() -> tuple[str, str]:
    """(system, machine) as the tables above spell them: Linux x86_64 / aarch64, Darwin x86_64 / arm64, Windows
    AMD64."""
    sysname, m = platform.system(), platform.machine()
    if sysname == "Windows":
        return sysname, "AMD64" if m.lower() in ("amd64", "x86_64") else m
    return sysname, {"amd64": "x86_64", "arm64": ("arm64" if sysname == "Darwin" else "aarch64")}.get(m.lower(), m)


def cache_dir() -> str:
    """COAL's cache folder for tools: $COAL_CACHE_DIR, else the user cache (~/.cache/coal or $XDG_CACHE_HOME/coal,
    ~/Library/Caches/coal on macOS, %LOCALAPPDATA%\\coal\\Cache on Windows)."""
    env = os.environ.get("COAL_CACHE_DIR")
    if env:
        return env
    home = os.path.expanduser("~")
    if sys.platform == "win32":
        return os.path.join(os.environ.get("LOCALAPPDATA", os.path.join(home, "AppData", "Local")), "coal", "Cache")
    if sys.platform == "darwin":
        return os.path.join(home, "Library", "Caches", "coal")
    return os.path.join(os.environ.get("XDG_CACHE_HOME", os.path.join(home, ".cache")), "coal")


def _cached(tool: str) -> list[str]:
    root = cache_dir()
    sub = "Library/bin" if os.name == "nt" else "bin"
    return [os.path.join(root, "tools", tool, tool + _EXE), os.path.join(root, "envs", tool, *sub.split("/"),
                                                                          tool + _EXE)]


def _find(tool: str, env_var: str) -> str | None:
    env = os.environ.get(env_var)
    if env:
        return env if os.path.isfile(env) else None
    p = shutil.which(tool)
    if p:
        return p
    for c in _cached(tool):
        if os.path.isfile(c) and (os.name == "nt" or os.access(c, os.X_OK)):
            return c
    return None


def find_clustalo() -> str | None:
    """Clustal Omega: $COAL_CLUSTALO, else `clustalo` on PATH, else the one install_clustalo() put in cache_dir()."""
    return _find("clustalo", "COAL_CLUSTALO")


def find_mafft() -> str | None:
    """MAFFT: $COAL_MAFFT, else `mafft` on PATH, else the one install_mafft() put in cache_dir()."""
    return _find("mafft", "COAL_MAFFT")


def tool_version(exe: str, timeout: float = 60) -> str:
    """The first line `exe --version` prints ('' if it does not run)."""
    try:
        r = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=timeout,
                           encoding="utf-8", errors="replace")
    except (OSError, subprocess.SubprocessError):
        return ""
    out = (r.stdout or r.stderr).strip().splitlines()
    return out[0].strip() if out else ""


def _download(url: str, dest: str, quiet: bool) -> str:
    """Fetch url into dest + ".part" (dest itself is not touched) and return the SHA-256 of what was fetched."""
    if not quiet:
        print(f"coal: downloading {url}", file=sys.stderr)
    tmp = dest + ".part"
    h = hashlib.sha256()
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fh:  # noqa: S310 (fixed https/http URLs)
        for chunk in iter(lambda: r.read(1 << 20), b""):
            h.update(chunk)
            fh.write(chunk)
    return h.hexdigest()


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _check_hash_argument(sha256: str | None) -> str | None:
    if sha256 is None:
        return None
    v = sha256.strip().lower()
    if len(v) != 64 or any(c not in "0123456789abcdef" for c in v):
        raise ValueError("sha256: expected the 64 hexadecimal digits of a SHA-256 sum")
    return v


def _fetch_verified(urls: Sequence[str], dest: str, quiet: bool, pinned: str | None, sha256: str | None = None,
                    allow_unverified: bool = False, what: str = "the file", hint: str = "") -> str:
    """Download one file (the first of `urls` that answers; https before http) to `dest`, and keep it only if its
    SHA-256 is the expected one: `pinned` (COAL's own table) or else `sha256` (the caller's). With neither, the
    download is refused before anything is fetched, unless allow_unverified is True; then the file is kept, a
    warning names its hash, and the caller has accepted that nobody checked it. Returns the SHA-256.

    Raises MsaError: no hash and no permission, every URL failed, or the hash does not match (the file is deleted,
    never moved to `dest`, never run)."""
    given = _check_hash_argument(sha256)
    if pinned and given and given != pinned.lower():
        raise MsaError(f"the SHA-256 given for {what} ({given}) is not the one COAL pins for it ({pinned}); "
                       f"leave sha256 out to use the pinned value")
    want = (pinned or given or "").lower() or None
    if want is None and not allow_unverified:
        plain = any(u.startswith("http://") for u in urls)
        raise MsaError(
            f"refusing to download and run {what}: no SHA-256 is known for {urls[-1]} (its publisher gives none"
            + (", and the host serves plain HTTP" if plain else "") + "). Download and check the file yourself and "
            "pass its hash (sha256=\"...\"; `coal pipeline install ... --sha256 HEX`), or accept an unverified "
            "download (allow_unverified=True; `--allow-unverified`)" + (". " + hint if hint else ""))
    tmp = dest + ".part"
    errors, got = [], None
    for url in urls:
        try:
            got = _download(url, dest, quiet)
            break
        except (OSError, ValueError) as e:   # URLError and HTTPError are OSErrors
            errors.append(f"{url}: {e}")
            if os.path.exists(tmp):
                os.remove(tmp)
    if got is None:
        raise MsaError("download failed (" + "; ".join(errors) + ")")
    if want is not None and got != want:
        os.remove(tmp)
        raise MsaError(f"{what} failed its SHA-256 check and was deleted: expected {want}, got {got} ({url})")
    if want is None:
        print(f"coal: warning: {what} from {url} was NOT verified (no SHA-256 is known for it); its SHA-256 is "
              f"{got}. Pass sha256=\"{got}\" next time to accept exactly this file.", file=sys.stderr)
    os.replace(tmp, dest)
    return got


def _conda_exe() -> str | None:
    """micromamba, mamba or conda on PATH, else $CONDA_EXE, else the micromamba of an earlier install."""
    for n in ("micromamba", "mamba", "conda"):
        p = shutil.which(n)
        if p:
            return p
    if os.environ.get("CONDA_EXE") and os.path.isfile(os.environ["CONDA_EXE"]):
        return os.environ["CONDA_EXE"]
    p = os.path.join(cache_dir(), "tools", "micromamba", "micromamba" + _EXE)
    return p if os.path.isfile(p) else None


def _get_micromamba(quiet: bool) -> str:
    """micromamba MICROMAMBA_VERSION into the cache, kept only if its SHA-256 is the pinned one."""
    entry = _MICROMAMBA.get(_platform())
    if entry is None:
        raise MsaError(f"no micromamba binary for {_platform()}; install conda, or use method='official'")
    asset, pinned = entry
    d = os.path.join(cache_dir(), "tools", "micromamba")
    os.makedirs(d, exist_ok=True)
    exe = os.path.join(d, "micromamba")
    _fetch_verified([_MICROMAMBA_URL.format(version=MICROMAMBA_VERSION, asset=asset)], exe, quiet, pinned,
                    what=f"micromamba {MICROMAMBA_VERSION}")
    os.chmod(exe, os.stat(exe).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return exe


def _conda_install(tool: str, spec: str, conda: str, quiet: bool, timeout: float) -> str:
    prefix = os.path.join(cache_dir(), "envs", tool)
    env = dict(os.environ)
    if os.path.basename(conda).startswith("micromamba"):
        env.setdefault("MAMBA_ROOT_PREFIX", os.path.join(cache_dir(), "mamba"))
    cmd = [conda, "create", "-y", "-p", prefix, "-c", "conda-forge", "-c", "bioconda", spec]
    if not quiet:
        print("coal: " + " ".join(cmd), file=sys.stderr)
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout, encoding="utf-8",
                       errors="replace")
    exe = _cached(tool)[1]
    if r.returncode != 0 or not os.path.isfile(exe):
        raise MsaError(f"installing {spec} with {conda} failed: {(r.stderr or r.stdout).strip()[-800:]}")
    return exe


def install_clustalo(method: str = "auto", force: bool = False, quiet: bool = False,
                     timeout: float = 1800, sha256: str | None = None, allow_unverified: bool = False) -> str:
    """Install Clustal Omega into cache_dir() and return its path (an existing one when force is False).

    method: "auto" (conda/mamba/micromamba if one is installed, else micromamba fetched into the cache, else the
    official binary; Windows: the official binary), "conda", "micromamba", "official".

    Verification: conda fetches from conda-forge / Bioconda over HTTPS and checks each package against the
    channel's index; micromamba itself is one pinned release with a pinned SHA-256. The official binary of
    clustal.org has no published checksum and no hash is pinned for it here, so that route needs `sha256` (the
    SHA-256 of the file you expect: the download is deleted if it differs) or allow_unverified=True (keep and run
    whatever arrives; the hash is printed). Without either it refuses before downloading."""
    _check_hash_argument(sha256)
    if not force:
        for c in _cached("clustalo"):
            if os.path.isfile(c) and tool_version(c).startswith("1.2"):
                return c
    sysname, machine = _platform()
    if method == "auto":
        order = ["official"] if sysname == "Windows" else (["conda", "official"] if _conda_exe() else
                                                           ["micromamba", "official"])
    elif method in ("conda", "micromamba", "official"):
        order = [method]
    else:
        raise ValueError("method: 'auto', 'conda', 'micromamba' or 'official'")
    errors = []
    for m in order:
        try:
            if m in ("conda", "micromamba"):
                conda = _conda_exe() if m == "conda" else _get_micromamba(quiet)
                if conda is None:
                    raise MsaError("no conda, mamba or micromamba found")
                exe = _conda_install("clustalo", f"clustalo={CLUSTALO_VERSION}", conda, quiet, timeout)
            else:
                exe = _install_official(sysname, machine, quiet, sha256=sha256, allow_unverified=allow_unverified)
            v = tool_version(exe)
            if not v.startswith("1.2"):
                raise MsaError(f"{exe} does not run or is not Clustal Omega 1.2 ({v or 'no output'})")
            if not quiet:
                print(f"coal: Clustal Omega {v} installed: {exe}", file=sys.stderr)
            return exe
        except (MsaError, OSError, subprocess.SubprocessError) as e:
            errors.append(f"{m}: {e}")
    raise MsaError("could not install Clustal Omega (" + "; ".join(errors) + "). Install it yourself (`conda install "
                   "-c bioconda clustalo`, `apt install clustalo`, `brew install clustal-omega`) or set COAL_CLUSTALO.")


def _install_official(sysname: str, machine: str, quiet: bool, sha256: str | None = None,
                      allow_unverified: bool = False) -> str:
    entry = _OFFICIAL.get((sysname, machine))
    if entry is None:
        raise MsaError(f"no official Clustal Omega binary for {sysname} {machine}")
    path, member, pinned = entry
    urls = [f"https://{_OFFICIAL_HOST}{path}", f"http://{_OFFICIAL_HOST}{path}"]
    hint = ""  # (install_clustalo's own message names the other ways to install it)
    d = os.path.join(cache_dir(), "tools", "clustalo")
    os.makedirs(d, exist_ok=True)
    exe = os.path.join(d, "clustalo" + _EXE)
    what = "the official Clustal Omega binary"
    if member is None:
        _fetch_verified(urls, exe, quiet, pinned, sha256, allow_unverified, what, hint)
    else:
        z = exe + ".zip"
        _fetch_verified(urls, z, quiet, pinned, sha256, allow_unverified, what, hint)
        with zipfile.ZipFile(z) as zf:
            names = [n for n in zf.namelist() if n.replace("\\", "/").rsplit("/", 1)[-1].lower() == member]
            if not names:
                raise MsaError(f"{urls[-1]} has no {member}")
            top = names[0].rsplit("/", 1)[0] + "/" if "/" in names[0] else ""
            # the Windows build needs the DLLs next to it: extract the whole folder, flattened
            for n in zf.namelist():
                if n.startswith(top) and not n.endswith("/"):
                    target = os.path.join(d, n[len(top):].replace("/", os.sep))
                    if os.path.commonpath([os.path.abspath(target), os.path.abspath(d)]) != os.path.abspath(d):
                        continue  # no path escapes
                    os.makedirs(os.path.dirname(target), exist_ok=True)
                    with zf.open(n) as src, open(target, "wb") as dst:
                        shutil.copyfileobj(src, dst)
        os.remove(z)
    if os.name != "nt":
        os.chmod(exe, os.stat(exe).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return exe


def install_mafft(method: str = "auto", force: bool = False, quiet: bool = False, timeout: float = 1800) -> str:
    """Install MAFFT (BSD licence) into cache_dir() through conda / mamba / micromamba (Bioconda). Not on Windows
    (use the official Windows package from https://mafft.cbrc.jp and set COAL_MAFFT)."""
    if not force:
        c = _cached("mafft")[1]
        if os.path.isfile(c):
            return c
    if platform.system() == "Windows":
        raise MsaError("install MAFFT from https://mafft.cbrc.jp/alignment/software/windows.html and set COAL_MAFFT")
    conda = _conda_exe() if method in ("auto", "conda") else None
    if conda is None:
        if method == "conda":
            raise MsaError("no conda, mamba or micromamba found")
        conda = _get_micromamba(quiet)
    return _conda_install("mafft", "mafft", conda, quiet, timeout)


# ---------------------------------------------------------------- running

@dataclass
class MsaResult:
    rows: list[tuple[str, str]] | None
    error: str | None = None


def _prepare(records: Sequence[tuple[str, str]]):
    """Sequences without stop marks for the aligner: trailing '*' removed, internal '*' -> X (positions kept)."""
    seqs, tails, internal = [], [], []
    for _, s in records:
        s = "".join(s.split())
        tail = s.endswith("*")
        if tail:
            s = s.rstrip("*")
        pos = [i for i, c in enumerate(s) if c == "*"]
        seqs.append(s.replace("*", "X"))
        tails.append(tail)
        internal.append(pos)
    return seqs, tails, internal


def _restore(rows: list[str], seqs, tails, internal) -> list[str]:
    out = []
    any_tail = any(tails)
    for row, s, tail, pos in zip(rows, seqs, tails, internal):
        if pos:
            r, k, want = list(row), -1, set(pos)
            for c, ch in enumerate(r):
                if ch != "-":
                    k += 1
                    if k in want:
                        r[c] = "*"
            row = "".join(r)
        out.append(row + ("*" if tail else "-") if any_tail else row)
    return out


def _read_fasta(path: str) -> dict[str, str]:
    recs, name, buf = {}, None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    recs[name] = "".join(buf)
                name, buf = line[1:].split()[0] if line[1:].split() else "", []
            elif name is not None:
                buf.append(line)
    if name is not None:
        recs[name] = "".join(buf)
    return recs


# Aligner processes that are running now (run_msa starts them from align_clusters' worker threads): kept so that an
# interrupted batch can end them (kill_running) and not leave them to finish on their own.
_running: set = set()
_running_lock = threading.Lock()


def _run_tracked(cmd, timeout, **kw) -> subprocess.CompletedProcess:
    """subprocess.run(cmd, timeout=..., **kw) with the process registered in _running while it runs."""
    proc = subprocess.Popen(cmd, **kw)
    with _running_lock:
        _running.add(proc)
    try:
        try:
            out, err = proc.communicate(timeout=timeout)
        except BaseException:   # the timeout, or an interrupt: never leave the aligner behind
            proc.kill()
            proc.communicate()
            raise
    finally:
        with _running_lock:
            _running.discard(proc)
    return subprocess.CompletedProcess(cmd, proc.returncode, out, err)


def kill_running() -> int:
    """End every Clustal Omega / MAFFT process this module started that is still running (an interrupted run).
    Returns how many were ended."""
    with _running_lock:
        procs = list(_running)
    n = 0
    for proc in procs:
        if proc.poll() is None:
            try:
                proc.kill()
                n += 1
            except OSError:
                pass
    return n


def run_msa(records: Sequence[tuple[str, str]], tool: str = "clustalo", threads: int = 1, exe: str | None = None,
            extra_args: Sequence[str] = (), timeout: float | None = None) -> list[tuple[str, str]]:
    """Align the protein records [(id, sequence)] with Clustal Omega ("clustalo") or MAFFT ("mafft"); returns
    [(id, gapped row)] in input order. One record comes back as it is (no alignment needed). Raises MsaError."""
    if tool not in ("clustalo", "mafft"):
        raise ValueError("tool: 'clustalo' or 'mafft'")
    records = list(records)
    if not records:
        return []
    seqs, tails, internal = _prepare(records)
    for (rid, _), s in zip(records, seqs):
        if not s:
            raise MsaError(f"sequence {rid!r} is empty")
    if len(records) == 1:
        return list(zip([records[0][0]], _restore([seqs[0]], seqs, tails, internal)))
    exe = exe or (find_clustalo() if tool == "clustalo" else find_mafft())
    if exe is None:
        raise MsaError(f"{tool} was not found: run coal.msa.install_{tool}() or put it on PATH")
    with tempfile.TemporaryDirectory(prefix="coal_msa_") as tmp:
        inp, outp = os.path.join(tmp, "in.fasta"), os.path.join(tmp, "out.fasta")
        with open(inp, "w", encoding="utf-8", newline="\n") as fh:
            for k, s in enumerate(seqs):  # positional names: any ID works, duplicates included
                fh.write(f">s{k}\n{s}\n")
        if tool == "clustalo":
            cmd = [exe, "-i", inp, "-o", outp, "--outfmt=fa", "--seqtype=Protein", "--output-order=input-order",
                   f"--threads={max(1, threads)}", "--force", *extra_args]
            r = _run_tracked(cmd, timeout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                             encoding="utf-8", errors="replace")
        else:
            cmd = [exe, "--auto", "--amino", "--anysymbol", "--quiet", "--thread", str(max(1, threads)), *extra_args,
                   inp]
            with open(outp, "w", encoding="utf-8") as fh:
                r = _run_tracked(cmd, timeout, stdout=fh, stderr=subprocess.PIPE, text=True, encoding="utf-8",
                                 errors="replace")
        if r.returncode != 0:
            raise MsaError(f"{tool} failed (exit {r.returncode}): {(r.stderr or '').strip()[-800:]}")
        got = _read_fasta(outp)
    rows = []
    for k, s in enumerate(seqs):
        row = got.get(f"s{k}")
        if row is None:
            raise MsaError(f"{tool} output lacks sequence {records[k][0]!r}")
        if row.replace("-", "").upper() != s.upper():
            raise MsaError(f"{tool} changed the residues of {records[k][0]!r}")
        rows.append(row)
    if len({len(r) for r in rows}) != 1:
        raise MsaError(f"{tool} returned rows of different lengths")
    return list(zip([rid for rid, _ in records], _restore(rows, seqs, tails, internal)))


def align_clusters(clusters: Mapping[str, Sequence[tuple[str, str]]] | Iterable[Sequence[tuple[str, str]]],
                   tool: str = "clustalo", jobs: int = 0, threads: int = 1, pairs_with: str = "nw",
                   exe: str | None = None, timeout: float | None = None, align_params=None) -> dict:
    """Align many clusters: {name: [(id, sequence)]} (or a list; names are then 0, 1, ...). Returns
    {name: MsaResult}; a cluster that fails has .error set and the others are unaffected.

    jobs aligner processes run at once (0 = one per CPU), each with `threads` threads: many small clusters are
    fastest with jobs = CPUs and threads = 1 (Clustal Omega's own threading helps only large clusters), a few big
    ones with fewer jobs and more threads. Clusters are started biggest first. pairs_with="nw": two-sequence
    clusters are aligned with COAL's global aligner (coal.align, BLOSUM62 needle defaults or align_params) instead
    of starting the MSA tool; "msa": with the tool too."""
    items = list(clusters.items()) if isinstance(clusters, Mapping) else list(enumerate(clusters))
    if pairs_with not in ("nw", "msa"):
        raise ValueError("pairs_with: 'nw' or 'msa'")
    results: dict = {}
    pairs = [(k, list(recs)) for k, recs in items if pairs_with == "nw" and len(recs) == 2]
    if pairs:
        from . import align as A
        params = align_params or A.AlignParams()
        prepared = [_prepare(recs) for _, recs in pairs]
        done = A.align_pairs([(p[0][0], p[0][1]) for p in prepared], params, jobs=jobs)
        for (k, recs), prep, res in zip(pairs, prepared, done):
            empty = [rid for (rid, _), s in zip(recs, prep[0]) if not s]
            if empty:
                results[k] = MsaResult(None, f"sequence {empty[0]!r} is empty")
            elif res.error:
                results[k] = MsaResult(None, res.error)
            else:
                seqs, tails, internal = prep
                rows = _restore([res.a, res.b], seqs, tails, internal)
                results[k] = MsaResult([(recs[0][0], rows[0]), (recs[1][0], rows[1])])
    rest = [(k, list(recs)) for k, recs in items if k not in results]
    if rest and any(len(r) > 1 for _, r in rest):
        exe = exe or (find_clustalo() if tool == "clustalo" else find_mafft())
    rest.sort(key=lambda kr: -sum(len(s) for _, s in kr[1]))
    workers = jobs if jobs > 0 else (os.cpu_count() or 1)

    def one(kr):
        k, recs = kr
        try:
            return k, MsaResult(run_msa(recs, tool, threads, exe, timeout=timeout))
        except (MsaError, OSError, subprocess.SubprocessError, ValueError) as e:
            return k, MsaResult(None, str(e))
    ex = ThreadPoolExecutor(max_workers=max(1, workers))
    try:
        for k, res in ex.map(one, rest):
            results[k] = res
    except BaseException:
        # an interrupt in the calling thread: drop the clusters not started yet and end the running aligners,
        # where the executor would otherwise align every remaining cluster before the interrupt is seen
        ex.shutdown(wait=False, cancel_futures=True)
        kill_running()
        raise
    ex.shutdown(wait=True)
    return {k: results[k] for k, _ in items}
