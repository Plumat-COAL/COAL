"""Multiple protein alignments of clusters for the codon step: Clustal Omega (default) or MAFFT (docs/alignment.md).

    from coal import msa
    exe = msa.find_clustalo() or msa.install_clustalo()      # PATH, COAL's cache, else install into the cache
    rows = msa.run_msa([("g1", "MKV..."), ("g2", "MRV...")])   # [(id, gapped row)], input order
    results = msa.align_clusters({"c1": [...], "c2": [...]}, jobs=8)   # many clusters, one process each

Neither tool is bundled with COAL: Clustal Omega is GPL-2.0 and COAL is MIT (docs/alignment.md "Licences").
install_clustalo() puts one into COAL's cache folder (cache_dir()), like orthologr's diamond_install(): through
conda / mamba / micromamba from Bioconda (micromamba itself is fetched from its GitHub releases, SHA-256 checked, if no
conda is installed), or the official static binary from clustal.org (method="official"; the only way on Windows).

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
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

__all__ = ["cache_dir", "find_clustalo", "find_mafft", "install_clustalo", "install_mafft", "tool_version",
           "run_msa", "align_clusters", "MsaResult", "MsaError", "CLUSTALO_VERSION"]

CLUSTALO_VERSION = "1.2.4"
_EXE = ".exe" if os.name == "nt" else ""

# Official Clustal Omega binaries (http://www.clustal.org/omega/). Their SHA-256 sums are not published; the binary
# is accepted when it runs and reports a 1.2 version.
_OFFICIAL = {
    ("Linux", "x86_64"): ("http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64", None),
    ("Darwin", "x86_64"): ("http://www.clustal.org/omega/clustal-omega-1.2.3-macosx", None),
    ("Darwin", "arm64"): ("http://www.clustal.org/omega/clustal-omega-1.2.3-macosx", None),  # x86_64, Rosetta 2
    ("Windows", "AMD64"): ("http://www.clustal.org/omega/clustal-omega-1.2.2-win64.zip", "clustalo.exe"),
}
# micromamba (static, MIT-licensed) from https://github.com/mamba-org/micromamba-releases; the release publishes a
# <asset>.sha256 next to every asset, which is checked.
_MICROMAMBA = {
    ("Linux", "x86_64"): "micromamba-linux-64",
    ("Linux", "aarch64"): "micromamba-linux-aarch64",
    ("Darwin", "x86_64"): "micromamba-osx-64",
    ("Darwin", "arm64"): "micromamba-osx-arm64",
}
_MICROMAMBA_URL = "https://github.com/mamba-org/micromamba-releases/releases/latest/download/"


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


def _download(url: str, dest: str, quiet: bool) -> None:
    if not quiet:
        print(f"coal: downloading {url}", file=sys.stderr)
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fh:  # noqa: S310 (fixed https/http URLs)
        shutil.copyfileobj(r, fh)
    os.replace(tmp, dest)


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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
    asset = _MICROMAMBA.get(_platform())
    if asset is None:
        raise MsaError(f"no micromamba binary for {_platform()}; install conda, or use method='official'")
    d = os.path.join(cache_dir(), "tools", "micromamba")
    os.makedirs(d, exist_ok=True)
    exe = os.path.join(d, "micromamba")
    _download(_MICROMAMBA_URL + asset, exe, quiet)
    sums = exe + ".sha256"
    _download(_MICROMAMBA_URL + asset + ".sha256", sums, quiet)
    with open(sums, encoding="utf-8") as fh:
        want = fh.read().split()[0].strip().lower()
    got = _sha256(exe)
    if got != want:
        os.remove(exe)
        raise MsaError(f"micromamba download failed its SHA-256 check (expected {want}, got {got})")
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
                     timeout: float = 1800) -> str:
    """Install Clustal Omega into cache_dir() and return its path (an existing one when force is False).

    method: "auto" (conda/mamba/micromamba if one is installed, else micromamba fetched into the cache, else the
    official binary; Windows: the official binary), "conda", "micromamba", "official"."""
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
                exe = _install_official(sysname, machine, quiet)
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


def _install_official(sysname: str, machine: str, quiet: bool) -> str:
    entry = _OFFICIAL.get((sysname, machine))
    if entry is None:
        raise MsaError(f"no official Clustal Omega binary for {sysname} {machine}")
    url, member = entry
    d = os.path.join(cache_dir(), "tools", "clustalo")
    os.makedirs(d, exist_ok=True)
    exe = os.path.join(d, "clustalo" + _EXE)
    if member is None:
        _download(url, exe, quiet)
    else:
        z = exe + ".zip"
        _download(url, z, quiet)
        with zipfile.ZipFile(z) as zf:
            names = [n for n in zf.namelist() if n.replace("\\", "/").rsplit("/", 1)[-1].lower() == member]
            if not names:
                raise MsaError(f"{url} has no {member}")
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
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, encoding="utf-8",
                               errors="replace")
        else:
            cmd = [exe, "--auto", "--amino", "--anysymbol", "--quiet", "--thread", str(max(1, threads)), *extra_args,
                   inp]
            with open(outp, "w", encoding="utf-8") as fh:
                r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True, timeout=timeout,
                                   encoding="utf-8", errors="replace")
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
    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        for k, res in ex.map(one, rest):
            results[k] = res
    return {k: results[k] for k, _ in items}
