"""`coal pipeline`: coding sequences in, dN/dS table out, in one step (docs/pipeline.md).

    coal pipeline pairs    --query q_cds.fa --subject s_cds.fa -o out_dir     # orthologs between two species
    coal pipeline clusters --cds all_cds.fa [--clusters clusters.tsv] -o out_dir   # gene families / clusters

Pair mode:    translate -> DIAMOND blastp (reciprocal best hit, or best hit) -> pairwise protein alignment (COAL's
              global aligner by default) -> codon alignment -> Ka/Ks -> out_dir/dnds.tsv, one row per pair and method.
Cluster mode: translate -> DIAMOND DeepClust (or a given clusters file) -> one protein MSA per cluster (Clustal Omega;
              two members: the global aligner) -> codon MSA -> Ka/Ks for every pair of the cluster ->
              out_dir/dnds.tsv and out_dir/cluster_summary.tsv.

This module only orchestrates. The work is done by `coal-cpp translate` and `coal-cpp codon-dnds` (cpp/src/dnds.cpp:
alignment of pairs, codon threading and Ka/Ks in one process, on threads, fed through a pipe: no process and no
temporary file per pair), by DIAMOND and by Clustal Omega / MAFFT (one process per cluster of three or more,
`-j` at a time). Sequences are fetched from disk by offset and streamed, so memory holds the IDs, the pair or
cluster lists and one batch of sequences, not the proteome.

One bad gene, pair or cluster never stops a run: it is written with status `error` (pairs) or `failed`, listed in
out_dir/failed.tsv and in the JSON report, and the exit code is 2 (0: nothing failed; 1: the run itself could not
be done, e.g. an unreadable input or a missing tool).

From Python:

    from coal import dnds_pipeline as P
    res = P.run_pairs("q_cds.fa", "s_cds.fa", "out", method=["YN"], jobs=8)       # -> PipelineResult
    res = P.run_clusters(["all_cds.fa"], "out", clusters="clusters.tsv", max_cluster_size=200)
    res.exit_code, res.report["counts"], res.outputs["dnds"]

Needs coal-cpp (the platform wheels have it; `coal --cpp-path`). Imports nothing outside the engine package.
The supported programmatic surface is what __all__ names, documented in docs/python-api.md ("dN/dS and the
pipeline") with its versioning promise; everything with a leading underscore, and the classes and helpers not in
__all__, are internal.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import math
import os
import platform
import random
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tarfile
import tempfile
import time
import zipfile
import zlib
from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence, Union

from . import __version__
from . import msa as _msa
from .core import CoalError
from .pipeline import release_output_lock, take_output_lock
from .tables import SUPPORTED

__all__ = ["run_pairs", "run_clusters", "main", "PipelineResult", "PipelineError", "find_diamond",
           "install_diamond", "codon_dnds", "summarise", "PAIR_COLUMNS", "SUMMARY_COLUMNS", "FAILED_COLUMNS",
           "METHODS", "SCHEMA"]

_Path = Union[str, "os.PathLike[str]"]
SCHEMA = "coal.pipeline/1"
DIAMOND_VERSION = "2.2.8"   # the version install_diamond() fetches and the one COAL is tested with
# The oldest DIAMOND the pipeline runs: 2.1.0 is the first release with `diamond deepclust` (cluster mode without
# --clusters stops with an older one; pair mode only warns: its blastp options are older than that, but no older
# version has been tested). Ubuntu 22.04's package is 2.0.14.
DIAMOND_MIN_VERSION = (2, 1, 0)
# KaKs_Calculator's order, then the two methods it does not have (docs/kaks.md): the order of the rows
METHODS = ("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN", "Comeron", "Li")
_METHOD_NAMES = {**{m.upper(): m for m in METHODS}, "C95": "Comeron"}


def _method_name(name: str) -> str:
    """The engine's name of a method, whatever the case ("comeron", "C95" -> "Comeron"); unknown names are
    returned as given, for the caller's error message."""
    return _METHOD_NAMES.get(str(name).upper(), str(name))
SENSITIVITIES = ("faster", "fast", "default", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive",
                 "ultra-sensitive")
# out_dir/dnds.tsv (cluster mode has `cluster_id` first); docs/pipeline.md and docs/pipeline.schema.json
PAIR_COLUMNS = ("seq1", "seq2", "method", "dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "Nd", "Sd",
                "codons_compared", "prot_aln_length", "prot_identity", "status", "note")
_STATS = ("median", "mean", "q25", "q75", "min", "max")
# out_dir/cluster_summary.tsv
SUMMARY_COLUMNS = ("cluster_id", "representative", "n_members", "n_used", "method", "pair_selection", "aligner",
                   "n_pairs", "n_ok", "n_not_converged", "n_saturated", "n_no_data", "n_failed",
                   *(f"{v}_{s}" for v in ("dN", "dS", "dNdS") for s in _STATS), "n_dNdS", "status", "note")
FAILED_COLUMNS = ("level", "cluster_id", "seq1", "seq2", "code", "message")
HASH_LIMIT = 2 * 1024 ** 3   # --hash-inputs auto: files up to this size get a SHA-256 in the report (2 GiB)
REPORT_MAX_FAILURES = 1000   # the JSON report lists this many failures; failed.tsv has all of them
_CLUSTER_HEADERS = {("centroid", "member"), ("representative", "member"), ("representative_id", "member_id"),
                    ("cluster", "member"), ("cluster_id", "member_id"), ("rep", "member")}
_EXE = ".exe" if os.name == "nt" else ""


# Files and folders a run writes into the output folder (besides work/): whatever an earlier run left of them is
# removed when a run starts, so that the folder never mixes two runs (a dnds.kaks or codon_alignments.fasta of a
# run with other options, the dnds.tsv.tmp of a killed run).
# Every table and alignment file is written as NAME.tmp and renamed when it is complete (_atomic, and the engine's
# --new-files), so the .tmp names are a run's files too: a run that was killed leaves them, the next one removes them.
_RUN_TABLES = ("dnds.tsv", "dnds.kaks", "codon_alignments.fasta", "protein_alignments.fasta", "cluster_summary.tsv",
               "genes.tsv", "orthologs.tsv", "unpaired.tsv", "clusters.tsv", "realign.tsv", "failed.tsv",
               "report.json")
_RUN_OUTPUTS = _RUN_TABLES + tuple(n + ".tmp" for n in _RUN_TABLES) + ("realign.tsv.body",)


def _open_new(path: str):
    """A text file created anew for writing: whatever has this name (a leftover, or a symbolic link, which is
    removed and never followed) goes first, then the file is created with O_CREAT | O_EXCL | O_NOFOLLOW (where
    the system has it). COAL writes its tables to NAME.tmp this way and renames them when they are complete."""
    try:
        os.remove(path)
    except FileNotFoundError:
        pass
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_BINARY", 0)
    return open(path, "w", encoding="utf-8", errors="backslashreplace", newline="\n",
                opener=lambda p, _flags: os.open(p, flags, 0o666))


@contextlib.contextmanager
def _atomic(path: str):
    """Write the table `path` through `path`.tmp: it appears under its name when it is complete (os.replace: a
    symbolic link at `path` is replaced, not followed), and nothing is left under either name when writing fails."""
    tmp = path + ".tmp"
    fh = _open_new(tmp)
    try:
        yield fh
        fh.close()
        os.replace(tmp, path)
    except BaseException:
        try:
            fh.close()
        except OSError:
            pass
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


class _Terminated(SystemExit):
    """A termination signal (SIGTERM, SIGHUP, Ctrl-Break) reached the Python process during a run: raised in the
    main thread so that the run ends as an interrupted one (children killed, report written, lock released). A
    SystemExit with the shell's code for the signal (128 + number), so an uncaught one ends the process quietly."""

    def __init__(self, signum: int):
        super().__init__(128 + signum)
        self.signum = signum
        try:
            self.name = signal.Signals(signum).name
        except ValueError:
            self.name = f"signal {signum}"


class PipelineError(RuntimeError):
    """The run cannot be done (exit 1): an unreadable input, a missing or failing tool, an unwritable folder."""


@dataclass
class PipelineResult:
    """What run_pairs() / run_clusters() return, for a failed run too (they do not raise for one).

    exit_code: 0 nothing failed; 2 the run finished and some genes, pairs or clusters failed (they are in
    outputs["failed"]); 1 the run could not be done (report["error"] says why).
    report: what out_dir/report.json holds (schema "coal.pipeline/1", docs/pipeline.schema.json); for a run that
    could not write into its folder, the few fields schema / exit_code / status / error.
    outputs: name -> path of the files that were written: "dnds", "genes", "failed", "report", and by mode and
    options "orthologs", "unpaired", "clusters", "realign", "cluster_summary", "kaks", "codon_alignments",
    "protein_alignments"."""
    exit_code: int
    report: dict[str, Any]
    outputs: dict[str, str] = field(default_factory=dict)


# ------------------------------------------------------------------------------------------------ tools

def usable_cpus() -> int:
    """The CPUs this process may run on: its affinity mask where the system has one (Linux: taskset, SLURM,
    cgroup cpusets), else os.process_cpu_count() (Python 3.13+), else the CPUs of the machine. At least 1."""
    try:
        n = len(os.sched_getaffinity(0))    # type: ignore[attr-defined]
    except (AttributeError, OSError):
        n = 0
    if n <= 0:
        try:    # (process_cpu_count asks for the affinity itself, and can fail the same way)
            n = getattr(os, "process_cpu_count", os.cpu_count)() or 0
        except OSError:
            n = 0
    return max(1, n or os.cpu_count() or 1)


def find_diamond() -> str | None:
    """DIAMOND: $COAL_DIAMOND, else `diamond` on PATH, else the one install_diamond() put in COAL's cache
    (coal.msa.cache_dir())."""
    return _msa._find("diamond", "COAL_DIAMOND")


_VERSION_RE = None


def _diamond_version(exe: str) -> str:
    """'2.2.8' from `diamond version` ("diamond version 2.2.8"): the first dotted number of the output, whatever
    follows it; when there is none, the last word (reported as it is; _version_tuple() then gives None and the
    caller warns). '' if the program does not run as DIAMOND."""
    global _VERSION_RE
    import re
    if _VERSION_RE is None:
        _VERSION_RE = re.compile(r"(?<![\w.])v?(\d+(?:\.\d+)+)")
    try:
        r = subprocess.run([exe, "version"], capture_output=True, text=True, timeout=60, stdin=subprocess.DEVNULL,
                           encoding="utf-8", errors="backslashreplace")
    except (OSError, subprocess.SubprocessError):
        return ""
    text = r.stdout or r.stderr
    words = text.split()
    if r.returncode != 0 or "diamond" not in text.lower() or not words:
        return ""
    m = _VERSION_RE.search(text)
    return m.group(1) if m else words[-1]


def _version_tuple(text: str) -> tuple[int, ...] | None:
    """'2.1.0' -> (2, 1, 0); 'v2.0.14' -> (2, 0, 14); '2.1' -> (2, 1, 0) (padded to three numbers, so that 2.1
    is not older than 2.1.0); None when the text does not start with numbers."""
    nums = []
    for part in text.strip().lstrip("vV").split("."):
        digits = ""
        for ch in part:
            if not ch.isdigit():
                break
            digits += ch
        if not digits:
            break
        nums.append(int(digits))
        if len(digits) != len(part):
            break       # "9-beta": the number ends here
    if not nums:
        return None
    return tuple(nums + [0] * (3 - len(nums)))


_DIAMOND_ASSETS = {("Linux", "x86_64"): "diamond-linux64.tar.gz", ("Darwin", "x86_64"): "diamond-macos.tar.gz",
                   ("Darwin", "arm64"): "diamond-macos.tar.gz", ("Windows", "AMD64"): "diamond-windows.zip"}
_DIAMOND_URL = "https://github.com/bbuchfink/diamond/releases/download/v{version}/{asset}"
# SHA-256 of the release assets COAL installs, by (version, asset). DIAMOND's releases publish no checksums: these
# were computed from the files downloaded from the URL above on 2026-11-12 (the Linux binary reports "diamond
# version 2.2.8"). An archive is unpacked only after its hash matched; a version that is not listed here needs
# sha256=... or allow_unverified=True (install_diamond).
_DIAMOND_SHA256 = {
    ("2.2.8", "diamond-linux64.tar.gz"): "33a765d192ae4f5b388f15fe2e7c23d3aaea0f9e3f6872efe8c324aec261fbc7",
    ("2.2.8", "diamond-macos.tar.gz"): "7d17c9c72b927d74ae1e72b4b7ca44f6100af02e9f753ad1ecded9abb30a3b22",
    ("2.2.8", "diamond-windows.zip"): "8035ece14562e13576da6e9583c74ab1aaff9d0c96086d93ab337e95cd288ce1",
}


def install_diamond(method: str = "auto", version: str = DIAMOND_VERSION, force: bool = False,
                    quiet: bool = False, timeout: float = 1800, sha256: str | None = None,
                    allow_unverified: bool = False) -> str:
    """Install DIAMOND into COAL's cache folder (coal.msa.cache_dir()) and return its path (an existing one when
    force is False). Works like coal.msa.install_clustalo() and orthologr's diamond_install().

    method: "official" (the pre-compiled binary of that version from DIAMOND's GitHub releases, over HTTPS: Linux
    x86_64, macOS, Windows x64), "conda" / "micromamba" (Bioconda, into the cache), or "auto" (official, then
    conda or micromamba).

    Verification: the archive of DIAMOND_VERSION is kept and unpacked only if its SHA-256 is the one pinned in
    this module (_DIAMOND_SHA256). For another `version` there is no pinned hash: give `sha256` (the download is
    deleted if it differs) or allow_unverified=True (unpack and run whatever arrives; the hash is printed);
    without either the official route refuses before downloading. Nothing downloaded runs before this check."""
    _msa._check_hash_argument(sha256)

    def have():
        for c in _msa._cached("diamond"):   # only one COAL installed and that is unchanged since
            if _msa._trusted_cached("diamond", c) and _diamond_version(c):
                return c
        return None
    if not force and have():
        return have()
    if method == "auto":
        order = ["official"] + ([] if platform.system() == "Windows" else
                                ["conda" if _msa._conda_exe() else "micromamba"])
    elif method in ("official", "conda", "micromamba"):
        order = [method]
    else:
        raise ValueError("method: 'auto', 'official', 'conda' or 'micromamba'")
    errors = []
    with _msa._install_lock("diamond", timeout):
        if not force and have():    # another install finished while this one waited
            return have()
        for m in order:
            try:
                if m == "official":
                    exe = _install_diamond_official(version, quiet, sha256, allow_unverified)
                else:
                    conda = _msa._conda_exe() if m == "conda" else _msa._get_micromamba(quiet)
                    if conda is None:
                        raise _msa.MsaError("no conda, mamba or micromamba found")
                    exe = _msa._conda_install("diamond", f"diamond={version}", conda, quiet, timeout)
                v = _diamond_version(exe)
                if not v:
                    raise _msa.MsaError(f"{exe} does not run")
                if not quiet:
                    print(f"coal: DIAMOND {v} installed: {exe}", file=sys.stderr)
                return exe
            except (_msa.MsaError, OSError, subprocess.SubprocessError, tarfile.TarError, zipfile.BadZipFile,
                    EOFError, zlib.error) as e:
                errors.append(f"{m}: {e}")
    raise PipelineError("could not install DIAMOND (" + "; ".join(errors) + "). Install it yourself (`conda install "
                        "-c bioconda diamond`, or https://github.com/bbuchfink/diamond/releases) or set COAL_DIAMOND.")


def _install_diamond_official(version: str, quiet: bool, sha256: str | None = None,
                              allow_unverified: bool = False) -> str:
    asset = _DIAMOND_ASSETS.get(_msa._platform())
    if asset is None:
        raise _msa.MsaError(f"no pre-compiled DIAMOND for {_msa._platform()}")
    version = version.lstrip("v")
    d = os.path.join(_msa.cache_dir(), "tools", "diamond")
    os.makedirs(d, exist_ok=True)
    archive = os.path.join(d, asset)
    _msa._fetch_verified([_DIAMOND_URL.format(version=version, asset=asset)], archive, quiet,
                         _DIAMOND_SHA256.get((version, asset)), sha256, allow_unverified,
                         what=f"DIAMOND {version} ({asset})",
                         hint=f"COAL pins the hashes of DIAMOND {DIAMOND_VERSION} only. Or install it another way: "
                              "`conda install -c bioconda diamond`; or set COAL_DIAMOND")
    url = _DIAMOND_URL.format(version=version, asset=asset)
    pinned = _DIAMOND_SHA256.get((version, asset))
    exe = os.path.join(d, "diamond" + _EXE)
    written: list[str] = []

    def put(src, name: str) -> None:
        # through a new file of its own, then renamed: never written through a link at the final name, never a
        # half-written binary under the final name, and a binary that is running can be replaced
        fd, tmp = tempfile.mkstemp(dir=d, prefix=name + ".", suffix=".part")
        try:
            with os.fdopen(fd, "wb") as dst:
                shutil.copyfileobj(src, dst)
            if os.name != "nt":
                os.chmod(tmp, 0o755)
            os.replace(tmp, os.path.join(d, name))
        except BaseException:
            try:
                os.remove(tmp)
            except OSError:
                pass
            raise
        written.append(os.path.join(d, name))
    try:
        try:
            if asset.endswith(".zip"):
                with zipfile.ZipFile(archive) as zf:
                    members = [(n, _safe_member_name(n, asset)) for n in zf.namelist() if not n.endswith(("/", "\\"))]
                    for n, base in members:   # diamond.exe and its DLLs, flattened; no name leaves the folder
                        with zf.open(n) as src:
                            put(src, base)
            else:
                with tarfile.open(archive) as tf:
                    member = next((m for m in tf.getmembers()
                                   if m.isfile() and os.path.basename(m.name) == "diamond"), None)
                    if member is None:
                        raise _msa.MsaError(f"{asset} holds no diamond binary")
                    with tf.extractfile(member) as src:
                        put(src, "diamond")
        except (EOFError, zlib.error, tarfile.TarError, zipfile.BadZipFile) as e:
            # (only an archive nobody verified can be broken here: allow_unverified, or a wrong sha256 given)
            raise _msa.MsaError(f"{asset} from {url} is not a complete archive ({type(e).__name__}: {e}); "
                                "the download was probably cut short. Nothing was installed") from None
    finally:
        os.remove(archive)
    if exe not in written:
        raise _msa.MsaError(f"{asset} holds no diamond binary")
    _msa._write_manifest(exe, url, "pinned SHA-256 of the archive" if pinned else
                         ("SHA-256 given by the caller" if sha256 else "NOT verified (allow_unverified)"),
                         also=[w for w in written if w != exe])
    return exe


_WINDOWS_RESERVED = {"CON", "PRN", "AUX", "NUL", *(f"COM{k}" for k in range(1, 10)), *(f"LPT{k}" for k in range(1, 10))}


def _safe_member_name(name: str, asset: str) -> str:
    """The file name under which a zip member is unpacked into the tool folder, or MsaError: members are
    flattened to their last path component, and a name that could mean another place on any system is refused
    (an absolute path, a `..` component, a drive letter or an alternate data stream (`:`), a reserved Windows
    device name, a trailing dot or space, anything outside letters, digits, `.`, `_`, `+`, `-`)."""
    parts = name.replace("\\", "/").split("/")
    base = parts[-1]
    ok = (bool(base) and not name.startswith(("/", "\\")) and ".." not in parts and ":" not in name
          and all(c.isascii() and (c.isalnum() or c in "._+-") for c in base)
          and not base.startswith(".") and not base.endswith(".")
          and base.split(".")[0].upper() not in _WINDOWS_RESERVED)
    if not ok:
        raise _msa.MsaError(f"{asset} holds a member with an unsafe name ({name!r}); nothing was installed")
    return base


def _find_cpp() -> str:
    from . import dispatch
    path = dispatch.find_cpp()
    if path is None:
        from .cli import NEEDS_CPP
        raise PipelineError(f"coal pipeline needs coal-cpp, which was not found (looked at {dispatch.LOOKUP_ORDER}): "
                            + NEEDS_CPP)
    ok, why = dispatch._version_ok(path)
    if not ok:
        raise PipelineError(why)
    return path


# ------------------------------------------------------------------------------------------------ small helpers

# Text and bytes (the same on every system; docs/pipeline.md "Text encoding"). Sequence IDs are bytes that COAL never
# re-encodes: the FASTA files are read as bytes, the engine is fed bytes through its pipe and writes its tables in
# binary mode. Every table and log this module reads or writes is opened as UTF-8 with the "backslashreplace" error
# handler and (when written) "\n" line ends, never with the system's default encoding (cp1252 on Windows, ASCII in
# the C locale) and never with its line ends. Bytes that are not valid UTF-8 (an ID in Latin-1, say) are therefore
# shown as \xhh, exactly as `coal-cpp translate` writes them into its genes table: nothing is dropped or replaced
# by a look-alike, two different IDs never become the same text, and every output file is valid UTF-8.
# tests/test_dnds_pipeline.py (Encoding) runs the pipeline under an ASCII locale and with Python's
# warn_default_encoding switched on.

class FastaIndex:
    """IDs and byte offsets of the records of FASTA files (the first record of a repeated ID); sequences are read
    from disk when asked for. Memory: the IDs."""

    def __init__(self, paths):
        self.paths = [os.fspath(p) for p in ([paths] if isinstance(paths, (str, os.PathLike)) else paths)]
        self.pos: dict[str, tuple[int, int, int]] = {}
        self._fh: dict[int, object] = {}
        for fno, path in enumerate(self.paths):
            self._scan(fno, path)

    def _scan(self, fno: int, path: str) -> None:
        try:
            fh = open(path, "rb")
        except OSError as e:
            raise PipelineError(f"cannot open {path!r}: {e.strerror or e}") from None
        with fh:
            off, cur, start = 0, None, 0
            pos = self.pos
            for line in fh:
                if line[:1] == b">":
                    if cur is not None:
                        pos.setdefault(cur, (fno, start, off))
                    words = line[1:].split()
                    cur = words[0].decode("utf-8", "backslashreplace") if words else ""
                    start = off + len(line)
                off += len(line)
            if cur is not None:
                pos.setdefault(cur, (fno, start, off))

    def __contains__(self, rid: str) -> bool:
        return rid in self.pos

    def __len__(self) -> int:
        return len(self.pos)

    def get(self, rid: str) -> bytes:
        """The sequence without whitespace (bytes); KeyError if there is no such record."""
        fno, a, b = self.pos[rid]
        fh = self._fh.get(fno)
        if fh is None:
            if len(self._fh) >= 32:   # many input files: keep few open
                self._fh.pop(next(iter(self._fh))).close()
            fh = self._fh[fno] = open(self.paths[fno], "rb")
        fh.seek(a)
        return fh.read(b - a).translate(None, b" \t\r\n\x0b\x0c")

    def close(self) -> None:
        for fh in self._fh.values():
            fh.close()
        self._fh.clear()


_BOM = b"\xef\xbb\xbf"


def fasta_lines(fh):
    """(offset, line) of the lines of a file opened in binary mode, split as coal-cpp splits them (cpp/src/seqio.cpp):
    a UTF-8 byte order mark at the start is skipped (Windows Notepad writes one), and a line ends at \\n, \\r\\n or a
    lone \\r (classic Mac OS). Each line keeps its line end; offset is where it starts in the file."""
    off = 0
    first = True
    for line in fh:
        if first:
            first = False
            if line.startswith(_BOM):
                off, line = 3, line[3:]
        k = line.find(b"\r")
        if k < 0 or k >= len(line) - 2 and line[k:] in (b"\r\n", b"\r"):
            yield off, line
            off += len(line)
            continue
        start, n = 0, len(line)
        while start < n:   # (one \n-line that holds lone \r line ends: a classic Mac file is one such "line")
            k = line.find(b"\r", start)
            if k < 0:
                end = n
            else:
                end = k + 2 if line[k + 1:k + 2] == b"\n" else k + 1
            yield off + start, line[start:end]
            start = end
        off += n


def fasta_ids(path: str) -> list[bytes]:
    """The IDs (first word of each header, bytes) of a FASTA file, in file order, repeats included."""
    ids = []
    with open(path, "rb") as fh:
        for _off, line in fasta_lines(fh):
            if line[:1] == b">":
                words = line[1:].split()
                ids.append(words[0] if words else b"")
    return ids


_FASTA_SUFFIXES = (".gz", ".fasta", ".fas", ".fa", ".fna", ".ffn", ".cds", ".fsa", ".seq", ".txt")


def file_stem(path: str) -> str:
    """The name of a file without its folder and its FASTA / gzip suffixes: dir/S_cer.cds.fa.gz -> S_cer."""
    name = os.path.basename(os.fspath(path))
    while True:
        low = name.lower()
        cut = next((x for x in _FASTA_SUFFIXES if low.endswith(x) and len(name) > len(x)), None)
        if cut is None:
            return name
        name = name[:-len(cut)]


ID_PREFIX_SEP = "|"


def write_prefixed_fasta(src: str, dst: str, prefix: bytes) -> int:
    """Copy the FASTA file src to dst with prefix + ID_PREFIX_SEP before the ID of every record (the rest of the
    header and the sequence lines as they are; line ends as \\n). Returns the number of records."""
    n = 0
    tag = b">" + prefix + ID_PREFIX_SEP.encode()
    with open(src, "rb") as fh, open(dst + ".tmp", "wb") as out:
        for _off, line in fasta_lines(fh):
            body = line.rstrip(b"\r\n")
            if body[:1] == b">":
                n += 1
                rest = body[1:]
                # (a header without an ID stays without one: coal-cpp reports it as such)
                out.write((tag + rest if rest[:1].strip() else body) + b"\n")
            else:
                out.write(body + b"\n")
    os.replace(dst + ".tmp", dst)
    return n


def sort_fasta_by_id(src: str, dst: str) -> int:
    """Write the records of the FASTA file src to dst in the order of their IDs (the first word of the header,
    compared as bytes; records with the same ID keep their order). Returns the number of records. Memory: the
    IDs and offsets, not the sequences."""
    recs: list[tuple[bytes, int, int, int]] = []
    with open(src, "rb") as fh:
        off, cur = 0, None
        for line in fh:
            if line[:1] == b">":
                if cur is not None:
                    recs.append((*cur, off))
                words = line[1:].split()
                cur = (words[0] if words else b"", len(recs), off)
            off += len(line)
        if cur is not None:
            recs.append((*cur, off))
        recs.sort(key=lambda r: (r[0], r[1]))
        with open(dst + ".tmp", "wb") as out:
            for _rid, _k, a, b in recs:
                fh.seek(a)
                block = fh.read(b - a)
                out.write(block if block.endswith(b"\n") else block + b"\n")
    os.replace(dst + ".tmp", dst)
    return len(recs)


def _utc_iso(seconds: float) -> str:
    """Seconds since the epoch as UTC, ISO 8601 with seconds: 2026-11-13T09:15:02Z."""
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(seconds))


def _utc_now() -> str:
    return _utc_iso(time.time())


def _file_sig(paths) -> list:
    out = []
    for p in paths:
        st = os.stat(p)
        out.append([os.path.abspath(p), st.st_size, st.st_mtime_ns])
    return out


def _sig(*parts) -> str:
    return hashlib.sha256(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:24]


def _num(x: float | None) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if x == 0:
        return "0"
    return format(x, ".8g")


def _quantile(v: list[float], q: float) -> float:
    """Linear interpolation between order statistics (R type 7, numpy's default); v sorted, not empty."""
    h = (len(v) - 1) * q
    lo = math.floor(h)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (h - lo)


def _stats(values: list[float]) -> list[str]:
    """median, mean, q25, q75, min, max as text (NA when there is no value)."""
    if not values:
        return ["NA"] * 6
    v = sorted(values)
    return [_num(_quantile(v, 0.5)), _num(math.fsum(v) / len(v)), _num(_quantile(v, 0.25)), _num(_quantile(v, 0.75)),
            _num(v[0]), _num(v[-1])]


def _tail(path: str, n: int = 800) -> str:
    try:
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            fh.seek(max(0, size - n))
            return fh.read().decode("utf-8", "backslashreplace").strip()
    except OSError:
        return ""


def _clean(text: str) -> str:
    return " ".join(str(text).split())


# ------------------------------------------------------------------------------------------------ the run

# ------------------------------------------------------------------------------------------------ genetic-code hint

GENETIC_CODE_NAMES = {1: "standard", 2: "vertebrate mitochondrial", 3: "yeast mitochondrial",
                      4: "mold / protozoan mitochondrial, Mycoplasma", 5: "invertebrate mitochondrial",
                      6: "ciliate nuclear", 9: "echinoderm / flatworm mitochondrial", 10: "euplotid nuclear",
                      11: "bacterial, archaeal and plant plastid", 12: "alternative yeast nuclear",
                      13: "ascidian mitochondrial", 14: "alternative flatworm mitochondrial",
                      15: "Blepharisma nuclear", 16: "chlorophycean mitochondrial", 21: "trematode mitochondrial",
                      22: "Scenedesmus obliquus mitochondrial", 23: "Thraustochytrium mitochondrial"}
# Codons a code does not use although NCBI's table gives them an amino acid ("CGA, CGC absent" in code 3): that a
# gene avoids them is evidence for the code, as avoiding a code's extra stop codons is (AGA / AGG for code 2).
_ABSENT_CODONS = {3: frozenset({b"CGA", b"CGC"})}
_CODON_FLOOR = 0.001        # the lowest usage per codon assumed when the set has no other genes to measure it on
_PSEUDO = re.compile(r"\[pseudo=true\]|biotype:\S*pseudogene")
_REPLICON = re.compile(r"^lcl\|(.+?)_cds_|\s(?:chromosome|primary_assembly|scaffold|supercontig|"
                                     r"plasmid|contig):[^:\s]*:([^:\s]+):")
HINT_BACKGROUND_CODONS = 1_000_000   # codons of the other genes read to measure the set's codon usage


def _and(items) -> str:
    items = [str(x) for x in items]
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _hint_inputs(paths, flagged: dict, usable: set):
    """For the genetic-code hint, one pass over the FASTA files: the flagged genes ({id, seq, header}), the codon
    usage of the other usable genes (internal codons, at most HINT_BACKGROUND_CODONS), and the replicon of every
    record when every header names one (NCBI lcl|NC_..._cds_, Ensembl chromosome:GRCh38:MT:), else None."""
    from collections import Counter
    genes: dict[str, dict] = {}
    background: Counter = Counter()
    n_bg = 0
    rep_of: dict[str, str] = {}
    rep_n: Counter = Counter()
    all_named = True
    seen: set[str] = set()

    def done(rid, header, chunks):
        nonlocal n_bg, all_named
        if rid is None or rid in seen:
            return
        seen.add(rid)
        m = _REPLICON.search(header)
        if m:
            rep = m.group(1) or m.group(2)
            rep_of[rid] = rep
            rep_n[rep] += 1
        else:
            all_named = False
        if rid in flagged:
            genes[rid] = {"id": rid, "header": header,
                          "seq": b"".join(chunks).upper().replace(b"U", b"T").translate(None, b" \t\r\n\x0b\x0c")}
        elif rid in usable and n_bg < HINT_BACKGROUND_CODONS:
            seq = b"".join(chunks).upper().replace(b"U", b"T").translate(None, b" \t\r\n\x0b\x0c")
            n = len(seq) - len(seq) % 3 - 3
            background.update(seq[k:k + 3] for k in range(0, n, 3))
            n_bg += max(0, n // 3)

    for path in paths:
        rid, header, chunks = None, "", []
        with open(path, "rb") as fh:
            for _off, line in fasta_lines(fh):
                if line[:1] == b">":
                    done(rid, header, chunks)
                    header = line[1:].decode("utf-8", "backslashreplace").strip()
                    words = header.split()
                    rid, chunks = (words[0] if words else ""), []
                    if rid not in flagged and (rid not in usable or n_bg >= HINT_BACKGROUND_CODONS):
                        chunks = None   # (nothing of its sequence is needed)
                elif chunks is not None:
                    chunks.append(line)
            done(rid, header, chunks or [])
    order = [genes[i] for i in flagged if i in genes]
    return order, background, ({"of": rep_of, "records": rep_n} if all_named and rep_of else None)


def suggest_genetic_code(chosen: int, genes: list[dict], background, n_records: int, with_stops: int,
                         replicons: dict | None = None, min_genes: int = 2) -> dict | None:
    """Which genetic code explains genes that have internal stop codons under `chosen`, or None (no hint).

    genes: {id, seq (upper case DNA), header} of the genes with two or more internal stops. background: codon
    counts of the set's other genes. The rules (deep test of 2026-11-14, finding D02: the old hint named the code
    that made the most genes stop-free, and that was 14 for yeast mitochondria and frameshifted pseudogenes, 4
    for vertebrate mitochondria):

    1. Evidence comes only from genes that look complete: not marked as pseudogenes ([pseudo=true],
       *_pseudogene biotype), and with a length that is a multiple of 3 or ends in the T / TA of an incomplete
       stop codon (vertebrate mitochondria); a frameshifted CDS has random stop codons that "fit" some code.
    2. When every header names its replicon, only genes on a replicon of which at least half the records are
       flagged count (the mitochondrial chromosome; not scattered nuclear genes).
    3. Genes are grouped by the codes under which they have no internal stop; the largest group decides. Of its
       codes, those win whose extra stop codons or unused codons (AGA / AGG for code 2, CGA / CGC for code 3) the
       genes avoid where the set's other genes use them (and lose when they do not), then those that reassign
       no stop codon the genes do not need (code 14 also reads TAA, which TGA-only genes do not need).
    4. Codes that tie are all named (their sequences cannot tell 4, 5 and 13 apart), less those whose start
       codons do not fit at least two and a quarter of the genes when another candidate's do.
    5. No hint unless at least min_genes genes and half of the evidence genes fit, and a code that reads TAA or
       TAG as an amino acid is suggested only when half the set has internal stops (the whole set's code).

    Returns {"codes": [...], "genes": [ids], "also": [codes that also fit these genes], "evidence": text}."""
    from collections import Counter
    from .tables import get_code
    std = get_code(chosen)
    stops_c = frozenset(x.encode() for x in std.codons.get("*", ()))
    tables = {}
    for t in SUPPORTED:
        if t == chosen:
            continue
        code = get_code(t)
        st = frozenset(x.encode() for x in code.codons.get("*", ()))
        if st != stops_c:   # (code 11 for code 1: the same stops, nothing to explain)
            tables[t] = (st, frozenset(x.encode() for x in code.starts))
    ev = []
    for g in genes:
        if _PSEUDO.search(g["header"]):
            continue
        seq = g["seq"]
        tr = len(seq) % 3
        if tr and seq[-tr:] != (b"T", b"TA")[tr - 1]:
            continue
        cod = Counter(seq[k:k + 3] for k in range(0, len(seq) - tr - (0 if tr else 3), 3))
        types = frozenset(cod) & stops_c
        if not types:
            continue
        clean = frozenset(t for t, (st, _) in tables.items() if not (cod.keys() & st))
        ev.append({"id": g["id"], "cod": cod, "types": types, "clean": clean, "first": seq[:3],
                   "n": sum(cod.values())})
    if replicons and all(g["id"] in replicons["of"] for g in ev):
        per = Counter(replicons["of"][g["id"]] for g in ev)
        keep = {r for r, n in per.items() if 2 * n >= replicons["records"].get(r, 0)}
        ev = [g for g in ev if replicons["of"][g["id"]] in keep]
    groups: dict[frozenset, list] = {}
    for g in ev:
        if g["clean"]:
            groups.setdefault(g["clean"], []).append(g)
    if not groups:
        return None
    G = max(groups.values(), key=len)
    sig = G[0]["clean"]
    needed = frozenset().union(*(g["types"] for g in G))
    n_cod = sum(g["n"] for g in G)
    seen = Counter()
    for g in G:
        seen.update(g["cod"])
    bg_total = sum(background.values())
    keys, notes = {}, {}
    for t in sig:
        st, _starts = tables[t]
        unneeded = len((stops_c - st) - needed)
        extra = (st - stops_c) | (_ABSENT_CODONS.get(t, frozenset()) - stops_c)
        e = 0
        if extra:
            obs = sum(seen[x] for x in extra)
            bg = sum(background[x] for x in extra)
            rate = bg / bg_total if (bg >= 20 or bg_total >= 20_000) and bg_total else 0.0
            source = "the set's other genes"
            if rate <= 0:
                rate, source = _CODON_FLOOR * len(extra), "the rarest usage of nuclear genes"
            exp = rate * n_cod
            if exp >= 5:
                e = 1 if obs <= 0.2 * exp else -1 if obs > 0.5 * exp else 0
                what = "reads as stops" if extra <= st else "does not use"
                notes[t] = (f"{obs:,} {' / '.join(sorted(x.decode() for x in extra))} codons, which code {t} {what}, "
                            f"where about {exp:,.0f} were expected from {source}")
        keys[t] = (-e, unneeded)
    best = min(keys.values())
    W = sorted(t for t in sig if keys[t] == best)
    if len(W) > 1:
        drop = set()
        for t in W:
            bad = sum(1 for g in G if g["first"] not in tables[t][1]
                      and any(g["first"] in tables[u][1] for u in W if u != t))
            if bad >= 2 and 4 * bad >= len(G):
                drop.add(t)
        if len(drop) < len(W):
            W = [t for t in W if t not in drop]
    support = [g for g in ev if all(t in g["clean"] for t in W)]
    if len(support) < min_genes or 2 * len(support) < len(ev):
        return None
    if needed - {b"TGA"} and 2 * with_stops < n_records:
        return None
    also = sorted(t for t in tables if t not in W and all(t in g["clean"] for g in support))
    kinds = sorted(x.decode() for x in needed)
    text = (f"their only internal stop codon is {kinds[0]}" if len(kinds) == 1
            else "their internal stop codons are " + " / ".join(kinds))
    if len(W) == 1 and best[0] < 0 and W[0] in notes:
        text += "; they have " + notes[W[0]]
    return {"codes": W, "genes": [g["id"] for g in support], "also": also, "evidence": text}


class _Run:
    """State shared by the stages of one run: folders, tools, the report, the failures file."""

    def __init__(self, mode: str, a: argparse.Namespace, argv):
        self.mode, self.a = mode, a
        self.out = os.path.abspath(a.out_dir)
        self.work = os.path.join(self.out, "work")
        self.t0 = time.time()
        self.quiet = a.quiet
        # Threads: -j 0 (the default) is one per CPU this process may use (its affinity: taskset, a SLURM or
        # cgroup allocation), not one per CPU of the machine; and no program is given more than that, coal-cpp
        # included (`-j 5000` made DIAMOND reserve memory for 5,000 threads and be killed, and each of the
        # engine's alignment threads may hold --max-cells of memory). The report has what was asked and used.
        self.cpus = usable_cpus()
        self.jobs_requested = a.jobs if a.jobs and a.jobs > 0 else 0
        self.threads = max(1, min(self.jobs_requested or self.cpus, self.cpus))
        self.jobs = self.threads           # coal-cpp's -j, and the size of the batches fed to it
        self.jobs_clamped = bool(self.jobs_requested) and self.threads != self.jobs_requested
        self.stages: list[dict] = []
        self.failures: list[dict] = []
        self.n_failures = 0
        self.counts: dict = {}
        self.tools: dict = {}
        self.outputs: dict = {}
        self.argv = list(argv) if argv is not None else None
        self._failed_fh = None
        self.state: dict = {}
        self.cpp = None
        self.gene_failures: dict = {}
        self.lock = None        # <out_dir>.coal-lock while this run owns the output folder
        self._inputs: list[tuple[str, str]] = []
        self._earlier: tuple[set, set] = (set(), set())
        self.pairs_in_place = False   # --pairs is DIR/orthologs.tsv itself: read, not written again
        self._engine = None     # the running coal-cpp codon-dnds (killed if the run is interrupted)
        self.started = _utc_now()
        self.input_records: list[dict] = []   # the report's `inputs`
        self._commands: list[list[str]] = []  # the external commands of the stage that is running
        self._plain: dict[str, str] = {}   # a gzip-compressed input -> its decompressed copy in work/inputs
        self._prefixed: dict[str, str] = {}   # --prefix-ids: a --cds file -> its copy with prefixed IDs
        self.id_prefixes: dict[str, str] = {}  # --prefix-ids: a --cds file (as given) -> its prefix
        self.notes: list[str] = []         # warnings of the run, also in the report
        self.isoforms: dict = {}           # the report's `isoforms`: per set what was found and done
        self._gene_map: dict | None = None
        self._first_ids: dict = {}
        self._plain_sha256: dict = {}
        self.code_hints: list = []         # sets whose genes with internal stops are clean under another code
        self._engine_files: list = []
        self.set_counts: dict = {}         # per input set: records, not usable, dropped (as genes.tsv lists them)

    # -- set-up
    def start(self) -> None:
        try:
            os.makedirs(self.out, exist_ok=True)
        except OSError as e:
            raise PipelineError(f"cannot create the output folder {self.out!r}: {e.strerror or e}") from None
        # One run per output folder: COAL's lock convention (<out>.coal-lock next to the output, holding the
        # process id; a lock whose process is gone is replaced). Two runs into one folder shared work/ and
        # dnds.tsv.tmp and one of them failed or, worse, finished with the other's files.
        try:
            self.lock = take_output_lock(self.out, self.out)
        except CoalError as e:
            raise PipelineError(str(e).replace("write to another file", "use another output folder")) from None
        try:
            self._check_inputs()
            self._earlier = self._earlier_outputs()
        except PipelineError:
            # nothing is written into the folder (no report either: report.json could be the input itself)
            self.release()
            raise
        self._clear_earlier_outputs()
        try:
            os.makedirs(os.path.join(self.work, "logs"), exist_ok=True)
            self._record_outputs()
        except OSError as e:
            raise PipelineError(f"cannot create the output folder {self.out!r}: {e.strerror or e}") from None
        self.cpp = _find_cpp()
        self.tools["coal-cpp"] = {"path": self.cpp, "version": __version__}
        sp = os.path.join(self.work, "state.json")
        if self.a.resume and os.path.isfile(sp):
            try:
                with open(sp, encoding="utf-8", errors="backslashreplace") as fh:
                    self.state = json.load(fh)
            except (OSError, ValueError):
                self.state = {}
        self.outputs["failed"] = os.path.join(self.out, "failed.tsv")
        # failed.tsv.tmp while the run goes on; finish() gives it its name, for a failed run too (it is complete)
        self._failed_fh = _open_new(self.outputs["failed"] + ".tmp")
        self._failed_fh.write("\t".join(FAILED_COLUMNS) + "\n")
        if self.jobs_clamped:
            self.warn(f"-j {self.jobs_requested} is more than the {self.cpus} CPU(s) this process may use: "
                      f"coal-cpp, DIAMOND and the aligners are run with {self.threads} thread(s)")
        self._same_file_warning()
        self._gunzip_inputs()
        self._describe_inputs()

    def _describe_inputs(self) -> None:
        """The report's `inputs`: every file the run reads, with its size, modification time and SHA-256 (of the
        file as given, compressed or not). --hash-inputs auto (default) skips the hash of a file above HASH_LIMIT
        and says so; always hashes everything; never hashes nothing."""
        mode = getattr(self.a, "hash_inputs", "auto")
        t = time.perf_counter()
        for opt, p in self._inputs:
            if opt == "--diamond":
                continue
            rec: dict = {"option": opt, "path": os.path.abspath(p)}
            try:
                st = os.stat(p)
            except OSError as e:
                rec["error"] = e.strerror or str(e)
                self.input_records.append(rec)
                continue
            rec["bytes"] = st.st_size
            rec["modified"] = _utc_iso(st.st_mtime)
            if os.path.abspath(p) in self._plain:
                # sha256 below is of the file as given (the compressed bytes: the same FASTA compressed twice has
                # two of them); sha256_uncompressed is of its content, computed while it was decompressed
                rec["compression"] = "gzip"
                if mode != "never":
                    digest = self._plain_sha256.get(os.path.abspath(p))
                    if digest is None:     # the copy of a run made before this field existed: hash it now
                        h = hashlib.sha256()
                        with open(self._plain[os.path.abspath(p)], "rb") as fh:
                            for block in iter(lambda: fh.read(1 << 20), b""):
                                h.update(block)
                        digest = h.hexdigest()
                    rec["sha256_uncompressed"] = digest
            if mode == "never":
                rec["sha256"] = None
                rec["sha256_skipped"] = "--hash-inputs never"
            elif mode == "auto" and st.st_size > HASH_LIMIT:
                rec["sha256"] = None
                rec["sha256_skipped"] = (f"larger than {HASH_LIMIT:,} bytes (--hash-inputs auto); run with "
                                         "--hash-inputs always to hash it")
            else:
                h = hashlib.sha256()
                with open(p, "rb") as fh:
                    for block in iter(lambda: fh.read(1 << 20), b""):
                        h.update(block)
                rec["sha256"] = h.hexdigest()
            self.input_records.append(rec)
        self.hash_seconds = round(time.perf_counter() - t, 3)

    def warn(self, text: str) -> None:
        """One warning line on stderr (not with -q), kept for the report's `warnings`."""
        text = self.shown(text)
        self.notes.append(text)
        self.say(f"warning: {text}")

    def _same_file_warning(self) -> None:
        a = self.a
        if self.mode != "pairs":
            return
        same = False
        try:
            same = os.path.samefile(a.query, a.subject)
        except OSError:
            pass
        if same:
            self.warn(f"--query and --subject are the same file ({a.query}): every gene's best hit is itself, so "
                      "the pairs are self-pairs with dN = dS = 0. For paralogs within one genome use `coal pipeline "
                      "clusters`")

    def _gunzip_inputs(self) -> None:
        """A gzip-compressed input (found by its first two bytes, whatever its name) is decompressed into
        work/inputs/ and read from there: coal-cpp, DIAMOND's makedb aside, reads plain text. The copy of an
        earlier run is reused when the compressed file has not changed (so --resume can reuse the stages)."""
        from .cli import gunzip_to, is_gzip
        todo = [(opt, p) for opt, p in self._inputs if opt != "--diamond" and os.path.isfile(p) and is_gzip(p)]
        if not todo:
            return
        d = self.wpath("inputs")
        os.makedirs(d, exist_ok=True)
        record_path = os.path.join(d, "sources.json")
        try:
            with open(record_path, encoding="utf-8", errors="backslashreplace") as fh:
                record = json.load(fh)
        except (OSError, ValueError):
            record = {}
        t = time.perf_counter()
        n_bytes = 0
        for k, (opt, p) in enumerate(todo):
            name = os.path.basename(p)
            name = name[:-3] if name.lower().endswith(".gz") and len(name) > 3 else name + ".txt"
            dst = os.path.join(d, f"{k}_{name}")
            st = os.stat(p)
            sig = [os.path.abspath(p), st.st_size, st.st_mtime_ns]
            old = record.get(dst)
            if isinstance(old, list) and old[:3] == sig and os.path.isfile(dst):
                digest = old[3] if len(old) > 3 else None    # (a copy made by an earlier version has none)
            else:
                try:
                    digest = gunzip_to(p, dst)
                except ValueError as e:
                    raise PipelineError(f"{opt}: {e}") from None
                record[dst] = [*sig, digest]
                # the copy keeps the source's times, so an unchanged input gives the same --resume signature
                os.utime(dst, ns=(st.st_atime_ns, st.st_mtime_ns))
            self._plain[os.path.abspath(p)] = dst
            self._plain_sha256[os.path.abspath(p)] = digest
            n_bytes += os.path.getsize(dst)
        tmp = record_path + ".tmp"
        with open(tmp, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as fh:
            json.dump(record, fh)
        os.replace(tmp, record_path)
        self.stages.append({"name": "decompress_inputs", "seconds": round(time.perf_counter() - t, 3),
                            "cached": False, "files": len(todo), "bytes": n_bytes})
        self.say(f"{len(todo)} gzip-compressed input file(s) decompressed into {d}")

    def plain(self, path):
        """The file to read for an input path: itself, or its decompressed copy."""
        key = os.path.abspath(os.fspath(path))
        return self._prefixed.get(key) or self._plain.get(key, os.fspath(path))

    # -- the folder: this run's inputs, what it will write, what an earlier run wrote
    def _input_paths(self) -> list[tuple[str, str]]:
        """(option, path) of every file this run reads."""
        a = self.a
        got = []
        for opt in ("query", "subject", "pairs", "clusters", "gene_map", "diamond"):
            v = getattr(a, opt, None)
            if v:
                got.append(("--" + opt.replace("_", "-"), os.fspath(v)))
        for v in (getattr(a, "cds", None) or ()):
            got.append(("--cds", os.fspath(v)))
        return got

    def _is_input(self, path: str) -> str | None:
        """The option whose file `path` is (the same file by any name: a link, another spelling), or None."""
        for opt, inp in self._inputs:
            try:
                if os.path.samefile(path, inp):
                    return opt
            except OSError:
                pass
            if os.path.realpath(path) == os.path.realpath(inp):
                return opt
        return None

    def _planned_outputs(self) -> list[str]:
        """Names (in the output folder) of the files this run may write, temporary ones included."""
        a = self.a
        names = ["failed.tsv", "genes.tsv", "dnds.tsv"]
        if not a.report:
            names += ["report.json"]
        if a.kaks_format:
            names.append("dnds.kaks")
        if self.mode == "pairs":
            names += ["orthologs.tsv"]
            if not a.pairs:
                names.append("unpaired.tsv")
            if a.keep_alignments:
                names.append("codon_alignments.fasta")
            if a.keep_protein_alignments:
                names.append("protein_alignments.fasta")
        else:
            names += ["cluster_summary.tsv"]
            if not a.clusters:
                names.append("clusters.tsv")
                if a.realign:
                    names += ["realign.tsv", "realign.tsv.body"]
        return names + [n + ".tmp" for n in names if not n.endswith(".body")]

    def _check_inputs(self) -> None:
        """Refuse, before anything in the folder is removed or written, a run that would write onto one of its
        own inputs. The one exception is the documented `--pairs DIR/orthologs.tsv` (the pairs of an earlier
        run): that file is read and left as it is instead of being written again."""
        self._inputs = self._input_paths()
        self.pairs_in_place = False
        planned = [self.path(n) for n in self._planned_outputs()]
        if self.a.report:
            planned.append(os.path.abspath(self.a.report))
        for p in planned:
            opt = self._is_input(p)
            if opt is None:
                continue
            if opt == "--pairs" and os.path.basename(p) == "orthologs.tsv" and self.mode == "pairs":
                self.pairs_in_place = True
                continue
            raise PipelineError(f"{opt} {p} is a file this run writes: nothing was changed. Use another output "
                                "folder, or move the input out of this one")
        inside = (os.path.realpath(self.work), os.path.realpath(self.path("codon_alignments")),
                  os.path.realpath(self.path("protein_alignments")))
        for opt, inp in self._inputs:
            real = os.path.realpath(inp)
            for d in inside:
                if real.startswith(d + os.sep):
                    raise PipelineError(f"{opt} {inp} is inside {d}, a folder this run writes: nothing was "
                                        "changed. Use another output folder, or move the input out of it")

    # the folders of per-cluster alignment files (cluster mode): name -> the option that writes it
    _ALIGNMENT_DIRS = {"codon_alignments": "keep_alignments", "protein_alignments": "keep_protein_alignments"}

    def _earlier_outputs(self) -> tuple[set, set]:
        """What an earlier run is known to have written here: the names in work/outputs.json (written when a
        run starts, so a killed run counts), or, for a run of a version without that file, the outputs its
        report.json lists. Second value: the folders (codon_alignments, protein_alignments) whose .fasta files
        are that run's (a cluster run with --keep-alignments / --keep-protein-alignments)."""
        names: set = set()
        dirs: set = set()
        codon_dir = False
        try:
            with open(self.wpath("outputs.json"), encoding="utf-8", errors="backslashreplace") as fh:
                rec = json.load(fh)
            names = {n for n in rec["files"] if isinstance(n, str)}
            codon_dir = rec.get("codon_alignments_dir") is True
            if rec.get("protein_alignments_dir") is True:
                dirs.add("protein_alignments")
        except (OSError, ValueError, KeyError, TypeError):
            try:
                with open(self.path("report.json"), encoding="utf-8", errors="backslashreplace") as fh:
                    rep = json.load(fh)
                if isinstance(rep, dict) and rep.get("schema") == SCHEMA and isinstance(rep.get("outputs"), dict):
                    names = {v for v in rep["outputs"].values() if isinstance(v, str)}
                    names |= {"report.json", "failed.tsv"}
                    codon_dir = rep.get("mode") == "clusters" and "codon_alignments" in names
                    names |= {n + ".tmp" for n in list(names)}
            except (OSError, ValueError):
                pass
        if codon_dir:
            dirs.add("codon_alignments")
        for name, opt in self._ALIGNMENT_DIRS.items():
            d = self.path(name)
            if not (os.path.isdir(d) and not os.path.islink(d)):
                dirs.discard(name)
            if getattr(self.a, opt) and self.mode == "clusters" and os.path.lexists(d) and name not in dirs:
                # not COAL's: the engine would write its files among (or over) what the folder holds
                if not os.path.isdir(d) or os.path.islink(d) or os.listdir(d):
                    raise PipelineError(f"{d} exists and was not written by a COAL run: --{opt.replace('_', '-')} "
                                        "writes one file per cluster there and would mix with (or replace) what "
                                        "it holds. Nothing was changed: move it away or use another output folder")
        return {n for n in names if n in _RUN_OUTPUTS}, dirs

    def _clear_earlier_outputs(self) -> None:
        """Remove the files an earlier COAL run wrote into this folder, so the folder holds this run's results
        only. Only files a run is known to have written are removed (_earlier_outputs), never one of this run's
        inputs, and never a file COAL has no record of: a `clusters.tsv` or a codon_alignments/ (or
        protein_alignments/) folder that was the user's stays. --resume keeps the stage outputs it can reuse (clusters.tsv, realign.tsv; everything
        under work/ is kept anyway)."""
        a = self.a
        keep = set()
        if a.resume and self.mode == "clusters" and not a.clusters:
            keep.add("clusters.tsv")
            if a.realign:
                keep.add("realign.tsv")
        known, dirs = self._earlier
        try:
            for name in sorted(known):
                p = self.path(name)
                if name not in keep and (os.path.isfile(p) or os.path.islink(p)) and self._is_input(p) is None:
                    os.remove(p)
            for d in (self.path(name) for name in sorted(dirs)):
                if not os.path.isdir(d) or os.path.islink(d):
                    continue
                for n in os.listdir(d):
                    p = os.path.join(d, n)
                    if n.endswith((".fasta", ".fasta.tmp")) and self._is_input(p) is None:
                        os.remove(p)
                if not os.listdir(d):
                    os.rmdir(d)
        except OSError as e:
            raise PipelineError(f"cannot remove the earlier run's {e.filename or 'outputs'} in {self.out!r}: "
                                f"{e.strerror or e}") from None

    def _record_outputs(self) -> None:
        """work/outputs.json: the names this run may write, for the next run's clean-up. Written before the first
        output, so that what a killed run left is known too."""
        planned = [n for n in self._planned_outputs() if not (n == "orthologs.tsv" and self.pairs_in_place)]
        tmp = self.wpath("outputs.json.tmp")
        with open(tmp, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as fh:
            json.dump({"coal_version": __version__, "files": sorted(planned),
                       "codon_alignments_dir": bool(self.a.keep_alignments and self.mode == "clusters"),
                       "protein_alignments_dir": bool(self.a.keep_protein_alignments and self.mode == "clusters")},
                      fh)
            fh.write("\n")
        os.replace(tmp, self.wpath("outputs.json"))

    def remove_temporaries(self) -> None:
        """Remove what is left of this run's NAME.tmp files (a stage that failed or was interrupted while it was
        writing): no cut-off table or alignment stays in the folder, under a final name or any other. The names
        are also in work/outputs.json, so the next run removes what a killed run could not."""
        names = [n for n in self._planned_outputs() if n.endswith((".tmp", ".body"))]
        if self._failed_fh is not None:
            names.remove("failed.tsv.tmp")     # still being written: finish() gives it its name
        for n in names:
            try:
                os.remove(self.path(n))
            except OSError:
                pass
        if self.mode == "clusters":
            for name, opt in self._ALIGNMENT_DIRS.items():
                d = self.path(name)
                if getattr(self.a, opt) and os.path.isdir(d) and not os.path.islink(d):
                    try:
                        for n in os.listdir(d):
                            if n.endswith(".fasta.tmp"):
                                os.remove(os.path.join(d, n))
                    except OSError:
                        pass

    def release(self) -> None:
        release_output_lock(self.lock)
        self.lock = None

    def kill_children(self) -> None:
        """End the programs this run started (interrupt): the engine, and any aligner process of coal.msa."""
        proc = self._engine
        if proc is not None and proc.poll() is None:
            try:
                proc.kill()
                proc.wait(timeout=30)
            except (OSError, subprocess.SubprocessError):
                pass
        _msa.kill_running()

    def shown(self, text: str) -> str:
        """text with the paths of decompressed copies (work/inputs/...) replaced by the files as the user named
        them; lines that quote a command line (`$ ...`, from a log) are left as they ran."""
        copies = [*self._prefixed.items(), *self._plain.items()]
        if not copies or not any(copy in text for _, copy in copies):
            return text
        lines = text.split("\n")
        for k, line in enumerate(lines):
            if line.startswith("$ "):
                continue
            for orig, copy in copies:
                line = line.replace(copy, orig)
            lines[k] = line
        return "\n".join(lines)

    def say(self, text: str) -> None:
        text = self.shown(text)
        if not self.quiet:
            print(f"coal pipeline: {text}", file=sys.stderr, flush=True)

    def path(self, *parts) -> str:
        return os.path.join(self.out, *parts)

    def wpath(self, *parts) -> str:
        return os.path.join(self.work, *parts)

    # -- failures
    def fail(self, level: str, code: str, message: str, cluster: str = "", seq1: str = "", seq2: str = "") -> None:
        self.n_failures += 1
        message = _clean(message)
        self._failed_fh.write("\t".join((level, cluster, seq1, seq2, code, message)) + "\n")
        if len(self.failures) < REPORT_MAX_FAILURES:
            e = {"level": level, "code": code, "message": message}
            if cluster:
                e["cluster_id"] = cluster
            if seq1:
                e["seq1"] = seq1
            if seq2:
                e["seq2"] = seq2
            self.failures.append(e)

    # -- stages
    def stage(self, name: str, sig: str, outputs: list[str], fn) -> dict:
        """Run fn() -> detail dict, or reuse the outputs of an earlier run (--resume) made with the same inputs
        and options (sig)."""
        old = self.state.get(name)
        if self.a.resume and old and old.get("sig") == sig and all(os.path.isfile(p) for p in outputs):
            rec = {"name": name, "seconds": 0.0, "cached": True, **(old.get("detail") or {})}
            self.say(f"{name}: reused from the earlier run (--resume)")
        else:
            t = time.perf_counter()
            self._commands = []
            detail = fn() or {}
            if self._commands:   # the argv of every program the stage ran (kept for --resume: a reused stage
                detail = {**detail, "commands": self._commands}   # still says how its files were made)
                self._commands = []
            rec = {"name": name, "seconds": round(time.perf_counter() - t, 3), "cached": False, **detail}
            self.state[name] = {"sig": sig, "detail": detail}
            self._save_state()
        self.stages.append(rec)
        return rec

    def _save_state(self) -> None:
        tmp = self.wpath("state.json.tmp")
        with open(tmp, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as fh:
            json.dump(self.state, fh)
        os.replace(tmp, self.wpath("state.json"))

    def run_tool(self, cmd: list[str], log: str, what: str) -> None:
        """One external program, no shell, its output appended to work/logs/<log>.log; PipelineError if it fails."""
        logp = self.wpath("logs", log + ".log")
        self._commands.append([str(c) for c in cmd])
        with open(logp, "ab") as fh:
            fh.write(("$ " + " ".join(shlex.quote(c) for c in cmd) + "\n").encode())
            fh.flush()
            try:
                r = subprocess.run(cmd, stdin=subprocess.DEVNULL, stdout=fh, stderr=subprocess.STDOUT)
            except OSError as e:
                raise PipelineError(f"{what}: could not start {cmd[0]} ({e.strerror or e})") from None
        if r.returncode != 0:
            raise PipelineError(f"{what} failed ({_msa.exit_text(r.returncode)}); the end of {logp}:\n{_tail(logp)}")

    def diamond(self) -> str:
        if "diamond" in self.tools:
            return self.tools["diamond"]["path"]
        exe = self.a.diamond or find_diamond()
        if exe is None and self.a.install_tools:
            exe = install_diamond(quiet=self.quiet)
        if exe is None or not (os.path.isfile(exe) or shutil.which(exe)):
            raise PipelineError("DIAMOND was not found (PATH, $COAL_DIAMOND, COAL's cache). Run again with "
                                "--install-tools, run `coal pipeline install diamond`, or give --diamond PATH"
                                + (" or --pairs TSV" if self.mode == "pairs" else " or --clusters TSV"))
        v = _diamond_version(exe)
        if not v:
            raise PipelineError(f"{exe} does not run as DIAMOND (`{exe} version` failed)")
        vt = _version_tuple(v)
        need = ".".join(map(str, DIAMOND_MIN_VERSION))
        if vt is not None and vt < DIAMOND_MIN_VERSION:
            newer = ("Install a newer one (`coal pipeline install diamond` puts DIAMOND " + DIAMOND_VERSION +
                     " into COAL's cache; `conda install -c bioconda diamond`) and give it with --diamond PATH or "
                     "COAL_DIAMOND, since an older `diamond` on PATH is found first")
            if self.mode == "clusters":
                raise PipelineError(f"DIAMOND {v} ({exe}) is too old: `diamond deepclust` needs DIAMOND {need} or "
                                    f"later (COAL is tested with {DIAMOND_VERSION}). {newer}; or give the clusters "
                                    "with --clusters TSV")
            self.warn(f"DIAMOND {v} ({exe}) is older than {need}, the oldest version COAL supports (it is tested "
                      f"with {DIAMOND_VERSION}): the search may fail or differ. {newer}")
        if vt is None:
            self.warn(f"the version of DIAMOND could not be read from `{exe} version` (it says {v!r}): COAL "
                      f"supports DIAMOND {need} and later and is tested with {DIAMOND_VERSION}; the check is skipped")
        self.tools["diamond"] = {"path": exe, "version": v}
        return exe

    def msa_tool(self, tool: str) -> str:
        if tool in self.tools:
            return self.tools[tool]["path"]
        exe = getattr(self.a, tool, None) or (_msa.find_clustalo() if tool == "clustalo" else _msa.find_mafft())
        if exe is None and self.a.install_tools:
            try:
                exe = _msa.install_clustalo(quiet=self.quiet) if tool == "clustalo" else \
                    _msa.install_mafft(quiet=self.quiet)
            except _msa.MsaError as e:
                raise PipelineError(str(e)) from None
        if exe is None:
            name = "Clustal Omega" if tool == "clustalo" else "MAFFT"
            if platform.system() == "Windows":   # nothing to install it from there: say what works instead
                raise PipelineError(
                    f"{name} was not found (PATH, $COAL_{tool.upper()}), and " + (_msa.NO_WINDOWS_CLUSTALO
                    if tool == "clustalo" else
                    "there is no automatic MAFFT install on Windows: install MAFFT yourself and set COAL_MAFFT "
                    "(untested with COAL), or run COAL under WSL, where `coal pipeline install mafft` works"))
            raise PipelineError(f"{name} was not found (PATH, $COAL_{tool.upper()}, COAL's cache). Run again with "
                                f"--install-tools, or run `coal pipeline install {tool}`")
        self.tools[tool] = {"path": exe, "version": _msa.tool_version(exe)}
        return exe

    # -- translate
    def translate(self, name: str, paths: list[str]) -> tuple[str, str, dict]:
        """coal-cpp translate -> (protein FASTA, genes table, counts)."""
        for p in paths:
            if not os.path.isfile(p):
                raise PipelineError(f"cannot open {p!r}: no such file")
        faa, genes = self.wpath(f"{name}.faa"), self.wpath(f"{name}.genes.tsv")

        plain = [self.plain(p) for p in paths]

        def go():
            self.run_tool([self.cpp, "translate", *plain, "-o", faa, "--genes", genes, "-c",
                           str(self.a.codontable)], "translate", f"translating {paths[0]}")
        self.stage(f"translate_{name}", _sig("translate", _file_sig(plain), self.a.codontable, __version__),
                   [faa, genes], go)
        return faa, genes, {}

    # -- clusters: the IDs of several --cds files
    def check_cds_ids(self) -> None:
        """Clusters mode: the --cds files are one set, so an ID must name one record across all of them. IDs that
        occur in two files (BRAKER / AUGUSTUS g1.t1, gene1 in every species) are refused before anything runs:
        coal-cpp translate keeps the first record of an ID, so every later file's gene of that ID would be left
        out, a whole species when all IDs clash. With --prefix-ids every file's IDs get the file's name as a
        prefix (S_cer|g1.t1) instead, in a copy in work/inputs/ that the run reads."""
        a = self.a
        paths = list(a.cds)
        prefix = bool(getattr(a, "prefix_ids", False))
        if len(paths) < 2 and not prefix:
            return
        if prefix:
            stems = [file_stem(p) for p in paths]
            bad = next((x for x in stems if not x or x != x.strip() or any(c in x for c in "|> \t")
                        or not x.isascii()), None)
            if bad is not None:
                raise PipelineError(f"--prefix-ids: the file name {bad!r} cannot be an ID prefix (it must be ASCII, "
                                    "without blanks, '|' or '>'): rename the file")
            seen: dict[str, str] = {}
            for p, x in zip(paths, stems):
                if x in seen and os.path.abspath(seen[x]) != os.path.abspath(p):
                    raise PipelineError(f"--prefix-ids: {seen[x]} and {p} would get the same prefix {x!r} (the "
                                        "file name without folder and suffix): rename one of them")
                seen.setdefault(x, p)
            d = self.wpath("inputs")
            os.makedirs(d, exist_ok=True)
            t = time.perf_counter()
            done: dict[str, str] = {}
            for k, (p, x) in enumerate(zip(paths, stems)):
                key = os.path.abspath(os.fspath(p))
                if key in done:
                    continue
                src = self.plain(p)
                dst = os.path.join(d, f"prefixed_{k}_{x}.fa")
                write_prefixed_fasta(src, dst, x.encode())
                st = os.stat(src)   # the copy keeps the source's time: an unchanged input -> the same --resume
                os.utime(dst, ns=(st.st_atime_ns, st.st_mtime_ns))   # signature (as for decompressed copies)
                self._prefixed[key] = done[key] = dst
                self.id_prefixes[os.fspath(p)] = x
            for rec in self.input_records:   # the report's `inputs`
                if rec.get("option") == "--cds":
                    x = next((x for p, x in zip(paths, stems) if os.path.abspath(p) == rec.get("path")), None)
                    if x is not None:
                        rec["id_prefix"] = x + ID_PREFIX_SEP
            self.stages.append({"name": "prefix_ids", "seconds": round(time.perf_counter() - t, 3), "cached": False,
                                "files": len(done)})
            self.say("--prefix-ids: the IDs of " + ", ".join(f"{p} start with {x}{ID_PREFIX_SEP}"
                                                              for p, x in list(self.id_prefixes.items())[:3])
                     + (" ..." if len(self.id_prefixes) > 3 else ""))
            return
        first: dict[bytes, int] = {}
        clash: dict[tuple[int, int], list[bytes]] = {}
        files: list[str] = []
        for k, p in enumerate(paths):
            key = os.path.abspath(os.fspath(p))
            if any(os.path.abspath(os.fspath(q)) == key for q in files):
                files.append(os.fspath(p))   # the same file named twice: translate reports its records as repeats
                continue
            files.append(os.fspath(p))
            try:
                ids = fasta_ids(self.plain(p))
            except OSError:
                return   # (translate says what is wrong with the file)
            for rid in dict.fromkeys(ids):   # (a repeat within one file is that file's own duplicate_id)
                j = first.setdefault(rid, k)
                if j != k and rid:
                    clash.setdefault((j, k), []).append(rid)
        if not clash:
            return
        total = sum(len(v) for v in clash.values())
        parts = []
        for (j, k), ids in list(clash.items())[:3]:
            ex = ", ".join(i.decode("utf-8", "backslashreplace") for i in ids[:3])
            parts.append(f"{len(ids):,} in both {paths[j]} and {paths[k]} (e.g. {ex})")
        raise PipelineError(
            f"--cds: {total:,} sequence ID(s) occur in more than one of the --cds files: " + "; ".join(parts)
            + ("; ..." if len(clash) > 3 else "") + ". The files are one set of genes, so every ID must name one "
            "record across them: only the first file's record of each such ID would be used and the other files' "
            "genes left out. Nothing was run. Run again with --prefix-ids, which puts each file's name before its "
            f"IDs ({file_stem(paths[next(iter(clash))[1]])}{ID_PREFIX_SEP}"
            f"{next(iter(clash.values()))[0].decode('utf-8', 'backslashreplace')}), or make the IDs unique")

    def read_genes(self, genes_path: str, set_name: str, out_fh, order: dict | None = None,
                   dropped: dict | None = None) -> dict:
        """Copy a genes table into out_dir/genes.tsv (with the `set` column), report the failed genes, and return
        {id: reason} of the failed ones. order: filled with {id: rank} of the usable genes (file order).
        dropped: {id: (gene, kept id)} of the records --longest-isoform leaves out: written with status `dropped`
        (code not_longest_isoform; the note names the gene and the record kept), counted apart, not in `order`."""
        failed: dict[str, str] = {}
        dropped = dropped or {}
        c = self.counts.setdefault("genes", {"total": 0, "ok": 0, "warning": 0, "failed": 0})
        if dropped:
            c.setdefault("dropped", 0)
        codes = self.gene_failures.setdefault(set_name, {"usable": 0, "codes": {}, "first": None})
        mine = self.set_counts.setdefault(set_name, {"total": 0, "failed": 0, "dropped": 0})
        with open(genes_path, encoding="utf-8", errors="backslashreplace") as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] in dropped and f[7] != "failed":
                    gene, kept = dropped[f[0]]
                    f[7], f[8] = "dropped", "not_longest_isoform"
                    f[9] = "; ".join(x for x in (f[9], f"gene {gene}: {kept} is kept (--longest-isoform)") if x)
                    out_fh.write(set_name + "\t" + "\t".join(f) + "\n")
                    c["total"] += 1
                    c["dropped"] += 1
                    mine["total"] += 1
                    mine["dropped"] += 1
                    continue
                out_fh.write(set_name + "\t" + line)
                c["total"] += 1
                mine["total"] += 1
                c[f[7]] = c.get(f[7], 0) + 1
                if f[7] == "failed":
                    mine["failed"] += 1
                    if f[8] != "duplicate_id":   # a repeated ID: its first record stands
                        failed.setdefault(f[0], f[9])
                    self.fail("gene", f[8], f"{set_name}: {f[9]}", seq1=f[0])
                    codes["codes"][f[8]] = codes["codes"].get(f[8], 0) + 1
                    if codes["first"] is None:
                        codes["first"] = (f[0], f[9])
                else:
                    codes["usable"] += 1
                    if order is not None:
                        order[f[0]] = len(order)
        return failed

    # -- a set that looks like another genetic code
    _CODE_NAMES = GENETIC_CODE_NAMES
    CODE_HINT_MIN_STOPS = 2      # a gene counts with this many internal stops under the chosen code ...
    CODE_HINT_MIN_GENES = 2      # ... and the hint needs this many such genes that are clean under another code
    CODE_HINT_MAX_GENES = 5000   # genes looked at per set (the first ones)

    def code_hint(self, set_name: str, paths: list[str], genes_path: str, shown: list[str] | None = None) -> None:
        """One warning when genes of this set have internal stop codons under the chosen genetic code (-c) and
        none under another code COAL supports: mitochondrial or plastid genes in a nuclear set, or the wrong -c
        for the whole set. The code (or the codes, when the sequences cannot tell them apart) is chosen by
        suggest_genetic_code(). The files are read again only when two or more genes have two or more internal
        stops, so a set without such genes costs nothing. (A gene with one internal TGA is not counted: a
        selenoprotein reads "clean" under every code where TGA is Trp.)"""
        chosen = int(self.a.codontable)
        flagged: dict[str, int] = {}     # id -> trailing nucleotides
        usable: set[str] = set()
        with_stops = 0
        n_records = 0
        with open(genes_path, encoding="utf-8", errors="backslashreplace") as fh:
            col = {n: k for k, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
            ki, si, ci, ti = col.get("internal_stops"), col.get("status"), col.get("code"), col.get("trailing_nt")
            if ki is None or si is None:
                return
            for line in fh:
                f = line.rstrip("\n").split("\t")
                try:
                    n = int(f[ki])
                except (ValueError, IndexError):
                    continue
                if ci is not None and f[ci] == "duplicate_id":
                    continue
                n_records += 1
                if f[si] != "failed":
                    usable.add(f[0])
                if n < 1:
                    continue
                with_stops += 1
                if n >= self.CODE_HINT_MIN_STOPS and len(flagged) < self.CODE_HINT_MAX_GENES:
                    try:
                        flagged[f[0]] = int(f[ti]) if ti is not None else 0
                    except (ValueError, IndexError):
                        flagged[f[0]] = 0
        if len(flagged) < self.CODE_HINT_MIN_GENES:
            return
        genes, background, replicons = _hint_inputs(paths, flagged, usable)
        hint = suggest_genetic_code(chosen, genes, background, n_records, with_stops, replicons,
                                    min_genes=self.CODE_HINT_MIN_GENES)
        if hint is None:
            return
        codes = hint["codes"]
        best = len(hint["genes"])
        examples = hint["genes"][:3]
        self.code_hints.append({"set": set_name, "genetic_code": chosen, "genes_with_internal_stops": with_stops,
                                "genes_looked_at": len(flagged), "clean_under": codes[0] if len(codes) == 1 else None,
                                "candidates": codes, "genes": best, "also_clean_under": hint["also"],
                                "evidence": hint["evidence"], "examples": examples})
        shown = [str(p) for p in (shown or paths)]      # the files as the user named them (not work/ copies)
        files = " ".join(shown[:2]) + (" ..." if len(shown) > 2 else "")
        head = (f"{set_name}: {best:,} of the {with_stops:,} genes with internal stop codons under genetic code "
                f"{chosen}")
        names = "; ".join(f"{t} {GENETIC_CODE_NAMES.get(t, 'NCBI')}" for t in codes)
        why = f" ({hint['evidence']})" if hint["evidence"] else ""
        if len(codes) == 1:
            t = codes[0]
            self.warn(head + f" have none under NCBI genetic code {t} ({GENETIC_CODE_NAMES.get(t, 'NCBI')}), the "
                      f"code that explains them best{why}, e.g. {', '.join(examples)}. They look like genes of "
                      "another genetic code (mitochondrial or plastid genes in a nuclear set?): their rows are "
                      f"estimates from mistranslated codons. Run them on their own: take them out of {files} into "
                      f"a file of their own and run `coal pipeline` twice, the rest with -c {chosen} and these with "
                      f"-c {t} (docs/pipeline.md \"Genes of another genetic code\")")
        else:
            self.warn(head + f" have none under any of the NCBI genetic codes {_and(codes)}{why}, e.g. "
                      f"{', '.join(examples)}. They look like genes of another genetic code (mitochondrial genes "
                      "in a nuclear set?), and their sequences cannot tell which of these it is: " + names + ". "
                      "Their rows are estimates from mistranslated codons. Run them on their own with the code of "
                      f"their organism and organelle: take them out of {files} into a file of their own and run "
                      f"`coal pipeline` twice, the rest with -c {chosen} (docs/pipeline.md \"Genes of another "
                      "genetic code\")")

    def n_records(self, set_name: str, what: str = "total") -> int:
        """The records of a set as genes.tsv lists them (a repeated ID is a record, and not usable): the same
        counts as the summary line."""
        return self.set_counts.get(set_name, {}).get(what, 0)

    # -- several transcripts per gene (coal/isoforms.py)
    def find_isoforms(self, set_name: str, paths: list[str], genes_path: str) -> dict:
        """What the records of one set say about their genes; with --longest-isoform the records to drop,
        {id: (gene, kept id)}. Fills self.isoforms[set_name]; the warning is printed by isoform_warning()."""
        from . import isoforms as I
        a = self.a
        lengths: dict[str, int] = {}
        ids: list[str] = []
        with open(genes_path, encoding="utf-8", errors="backslashreplace") as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[8] == "duplicate_id":
                    continue
                ids.append(f[0])
                if f[7] != "failed":
                    lengths[f[0]] = int(f[1])
        headers = I.read_headers([self.plain(p) for p in paths])
        rule = getattr(a, "isoform_rule", "auto") or "auto"
        if a.gene_map:
            if self._gene_map is None:
                try:
                    self._gene_map = I.read_gene_map(self.plain(a.gene_map), a.gene_map)
                except ValueError as e:
                    raise PipelineError(str(e)) from None
                except OSError as e:
                    raise PipelineError(f"cannot open {a.gene_map!r}: {e.strerror or e}") from None
                if not self._gene_map:
                    raise PipelineError(f"--gene-map {a.gene_map} names no record (it has no line with a record ID "
                                        "and a gene ID)")
                if self.id_prefixes:   # --prefix-ids: the map names the IDs as in the files (or prefixed ones)
                    xs = set(self.id_prefixes.values())
                    gm: dict[str, str] = {}
                    for r, g in self._gene_map.items():
                        x = r.split(ID_PREFIX_SEP, 1)[0] if ID_PREFIX_SEP in r else None
                        if x in xs:
                            gm.setdefault(r, g)
                        else:   # the same gene ID in two species is two genes
                            for x in xs:
                                gm.setdefault(f"{x}{ID_PREFIX_SEP}{r}", f"{x}{ID_PREFIX_SEP}{g}")
                    self._gene_map = gm
            # the map wins for the records it names; the others fall back to what the headers and IDs say
            det = I.from_gene_map(self._gene_map, headers, detection=rule != "map")
            unknown = det.records - (det.in_map or 0)
        else:
            det = I.detect(headers) if rule != "map" else I.from_gene_map({}, headers, detection=False)
            unknown = None
        dropped = I.choose_longest(det, lengths) if a.longest_isoform else {}
        by_suffix = sum(1 for r in dropped if r in det.by_suffix)
        if by_suffix * 2 > det.records and rule != "suffix":
            # a suffix alone is a guess: gene families (HSP70.1, HSP70.2) and genes numbered per scaffold
            # (scaffold7.1, scaffold7.2) look the same as transcripts. Never drop most of a file on a guess.
            ex = det.example
            suffix = det.fallback if det.kind == "map" else det.source
            raise PipelineError(
                f"--longest-isoform would drop {by_suffix:,} of the {det.records:,} {set_name} records on the "
                f"evidence of their IDs alone ({suffix}"
                + (f"; e.g. {', '.join(ex[1])} would be one gene, {ex[0]}" if ex else "") + "): more than half. "
                "IDs that differ in a number are transcripts in some annotations and different genes in others "
                "(paralogs HSP70.1 / HSP70.2, genes numbered per scaffold), and nothing in the headers names a "
                "gene. Nothing was dropped. If these are transcripts of one gene, run again with --isoform-rule "
                "suffix; otherwise give the gene of each record with --gene-map TSV (record ID, gene ID)")
        rec = {"records": det.records, "genes": det.genes, "genes_with_several_records": det.multi,
               "records_beyond_one_per_gene": det.extra, "gene_source": det.source, "dropped": len(dropped)}
        if unknown is not None:
            rec["records_not_in_gene_map"] = unknown
            rec["records_in_gene_map"] = det.in_map or 0
            if det.kind == "map" and det.fallback:
                rec["gene_source_of_records_not_in_gene_map"] = det.fallback
        if by_suffix:
            rec["dropped_by_id_suffix"] = by_suffix
        if det.conflicts:
            rec["gene_symbols_at_different_loci"] = det.conflicts
            sym, members, why = det.conflict_example
            rec["gene_symbol_example"] = {"symbol": sym, "records": members, "evidence": why}
        if det.refseq_ids:
            rec["refseq_transcript_ids"] = det.refseq_ids
        if det.example:
            rec["example"] = {"gene": det.example[0], "records": det.example[1]}
        self._first_ids[set_name] = [rid for rid, _ in headers[:3]]
        self.isoforms[set_name] = rec
        return dropped

    def check_gene_map(self) -> None:
        """--gene-map that names no record of any input set is an error (the usual cause: IDs of another release,
        or protein IDs for a file of transcripts). Called when every set has been read."""
        a = self.a
        sets = self.isoforms
        ids = self._first_ids
        if not a.gene_map or not sets or any(r.get("records_in_gene_map", 0) for r in sets.values()):
            return
        mapped = list(self._gene_map or {})[:3]
        raise PipelineError(
            f"--gene-map {a.gene_map} names none of the records of the input: its record IDs (first column) are "
            f"e.g. {', '.join(mapped)}; the input's are e.g. "
            + "; ".join(f"{', '.join(v)} ({n})" for n, v in ids.items() if v)
            + ". The first column must hold the IDs of the FASTA records (the header up to the first blank)")

    @staticmethod
    def _source_text(r: dict) -> str:
        if r["gene_source"] != "gene map":
            return str(r["gene_source"])
        fb = r.get("gene_source_of_records_not_in_gene_map")
        return "--gene-map" + (f", and {fb} for the records it does not name" if fb else "")

    def isoform_warning(self) -> None:
        """One line for the run: input that looks like several transcripts per gene (and nothing done about it),
        or what --longest-isoform did."""
        a = self.a
        self.check_gene_map()
        sets = self.isoforms
        if a.gene_map:
            for n, r in sets.items():
                if not r.get("records_in_gene_map") and not (a.longest_isoform and r["gene_source"] is None):
                    self.warn(f"--gene-map {a.gene_map} names no record of the {n} file(s): "
                              + (f"their genes are read from the {r['gene_source']}" if r["gene_source"]
                                 else "nothing else names their genes, so every record there is a gene of its own"))
        for n, r in sets.items():
            if r.get("gene_symbols_at_different_loci"):
                ex = r["gene_symbol_example"]
                self.warn(f"{n}: {r['gene_symbols_at_different_loci']:,} gene symbol(s) are shared by records at "
                          f"different loci (e.g. {ex['symbol']}: {', '.join(ex['records'])}; {ex['evidence']}): "
                          "these are different genes with one name (copies of a gene, several strains or genomes "
                          "in one file), so they are not merged; each record is a gene of its own. Give the genes "
                          "with --gene-map TSV if some of them are transcripts of one gene")
        refseq = [n for n, r in sets.items() if r["gene_source"] is None and r.get("refseq_transcript_ids", 0) * 2
                  >= max(1, r["records"])]
        refseq_text = ("RefSeq transcript accessions (NM_, XM_) hold no gene in the ID, and the headers of the "
                       f"{' and '.join(refseq)} file(s) have no gene field: several transcripts of one gene cannot "
                       "be seen. Give the genes with --gene-map TSV (record ID, gene ID; e.g. from the GFF or "
                       "gene2refseq), or use NCBI's CDS file (cds_from_genomic), whose headers name the gene")
        if a.longest_isoform:
            total = sum(r["dropped"] for r in sets.values())
            none = [n for n, r in sets.items() if r["gene_source"] is None]
            parts = [f"{n}: {r['dropped']:,} of {r['records']:,} records dropped, {r['genes']:,} genes "
                     f"(from {self._source_text(r)})"
                     for n, r in sets.items() if r["gene_source"] is not None]
            self.say("--longest-isoform: " + ("; ".join(parts) if parts else "no record dropped")
                     + (f"; {total:,} in all, listed in genes.tsv with status dropped" if total else ""))
            plain = [n for n in none if n not in refseq]
            if plain:
                self.warn(f"--longest-isoform: the headers of the {' and '.join(plain)} file(s) name no gene (no "
                          "locus_tag= / GeneID: / gene: / [gene= field, no shared ID stem with a transcript "
                          "suffix)" + (f", and --gene-map {a.gene_map} names none of their records" if a.gene_map
                                       else "") + ", so every record was kept there; give the genes with "
                          "--gene-map TSV (record ID, gene ID)")
            if refseq:
                self.warn("--longest-isoform: " + refseq_text + ". Every record was kept there")
            if a.gene_map:
                missing = sum(r.get("records_not_in_gene_map", 0) for r in sets.values()
                              if r.get("records_in_gene_map"))
                named = sum(1 for r in sets.values() if r.get("records_in_gene_map")
                            and r.get("gene_source_of_records_not_in_gene_map"))
                if missing:
                    self.warn(f"{missing:,} record(s) are not in --gene-map {a.gene_map}: "
                              + ("their genes are read from the headers or IDs where those name one, and each "
                                 "other record is kept as a gene of its own" if named
                                 else "each is kept as a gene of its own"))
            return
        if refseq:
            self.warn(refseq_text)
        multi = {n: r for n, r in sets.items() if r["genes_with_several_records"]}
        if not multi:
            return
        parts = []
        for n, r in multi.items():
            ex = r.get("example")
            parts.append(f"{n}: {r['records']:,} records of {r['genes']:,} genes, "
                         f"{r['genes_with_several_records']:,} genes with several records"
                         + (f" (e.g. {ex['gene']}: {', '.join(ex['records'])})" if ex else "")
                         + f", read from the {self._source_text(r)}")
        self.warn("the input looks like several transcripts per gene (" + "; ".join(parts) + "). The pipeline "
                  "treats each record as a gene: such a gene is compared once per transcript, and reciprocal best "
                  "hits are split between its transcripts. Use --longest-isoform to keep the longest CDS of each "
                  "gene (with --gene-map TSV if the genes are not in the headers), or give one sequence per gene")

    def kept_proteins(self, faa: str, dropped: dict, name: str) -> str:
        """The protein file without the dropped records (what DIAMOND gets); faa itself when nothing is dropped.
        Rewritten only when its content changes, so --resume can reuse the search."""
        if not dropped:
            return faa
        out = self.wpath(f"{name}.kept.faa")
        with open(faa, "rb") as src, open(out + ".tmp", "wb") as dst:
            keep = True
            for line in src:
                if line[:1] == b">":
                    words = line[1:].split()
                    keep = (words[0].decode("utf-8", "backslashreplace") if words else "") not in dropped
                if keep:
                    dst.write(line)
        import filecmp
        if os.path.isfile(out) and filecmp.cmp(out, out + ".tmp", shallow=False):
            os.remove(out + ".tmp")
        else:
            os.replace(out + ".tmp", out)
        return out

    def require_usable_genes(self, set_name: str, paths: list[str]) -> None:
        """Stop here, before DIAMOND or an aligner is started, when a CDS file has no usable gene at all: the
        usual cause is a protein file given where coding DNA is expected. (DIAMOND would otherwise fail on the
        empty protein file with "Error detecting input file format".) Every gene is in genes.tsv and failed.tsv."""
        g = self.gene_failures.get(set_name)
        if g is None or g["usable"] > 0:
            return
        total = sum(g["codes"].values())
        files = ", ".join(paths)
        first = f"; the first: {g['first'][0]}: {g['first'][1]}" if g["first"] else ""
        bad_letters = g["codes"].get("cds_invalid_character", 0)
        if total and bad_letters * 2 >= total:
            raise PipelineError(
                f"{files} does not hold coding DNA: {bad_letters:,} of its {total:,} sequences have letters that "
                f"are no nucleotide codes{first}. Is this a protein file? The pipeline needs the coding sequences "
                "(CDS, nucleotides) of the genes and translates them itself")
        what = ", ".join(f"{n:,} {code}" for code, n in sorted(g["codes"].items(), key=lambda kv: (-kv[1], kv[0])))
        raise PipelineError(f"none of the {total:,} sequences of {files} can be used as a coding sequence "
                            f"({what}){first}. Every gene and its reason is in {self.path('failed.tsv')}")

    # -- the engine step
    def engine_cmd(self, group_column: str | None) -> tuple[list[str], str, str]:
        a = self.a
        out = self.path("dnds.tsv")
        fail = self.wpath("engine_failures.tsv")
        # Every file the engine writes into the output folder is NAME.tmp, created anew (--new-files: a leftover or
        # a symbolic link of that name is removed, never written through), and gets its name when the engine has
        # ended with success (run_engine): a failed or interrupted stage leaves no cut-off file under a final name.
        self._engine_files = [("dnds", out)]
        cmd = [self.cpp, "codon-dnds", "-i", "-", "-o", out + ".tmp", "--new-files", "--failures", fail,
               "-c", str(a.codontable),
               "-j", str(self.jobs), "--digits", str(a.digits), "--matrix", a.matrix, "--gap-open", str(a.gap_open),
               "--gap-extend", str(a.gap_extend), "--end-gaps", a.end_gaps, "--max-cells", str(a.max_cells),
               "--flanks", a.flanks, "--yn-flavour", a.yn_flavour, "--yn-weighting", str(a.yn_weighting),
               "--comeron-flavour", a.comeron_flavour, "--li-flavour", a.li_flavour,
               "--terminal-stops", a.terminal_stops]
        for m in a.method:
            cmd += ["-m", m]
        if a.min_codons:
            cmd += ["--min-codons", str(a.min_codons)]
        if a.max_ds is not None:
            cmd += ["--max-ds", repr(float(a.max_ds))]
        if getattr(a, "min_identity", None) is not None:
            cmd += ["--min-identity", repr(float(a.min_identity))]
        if group_column:
            cmd += ["--group-column", group_column]
        if a.kaks_format:
            self._engine_files.append(("kaks", self.path("dnds.kaks")))
            cmd += ["--kaks-out", self.path("dnds.kaks") + ".tmp"]
        if a.keep_alignments:
            if self.mode == "pairs":
                self._engine_files.append(("codon_alignments", self.path("codon_alignments.fasta")))
                cmd += ["--codon-out", self.path("codon_alignments.fasta") + ".tmp"]
            else:
                d = self.outputs["codon_alignments"] = self.path("codon_alignments")
                os.makedirs(d, exist_ok=True)   # (an earlier run's files were removed when this run started)
                cmd += ["--codon-dir", d]
        if a.keep_protein_alignments:   # the protein alignments the codons were threaded on, in the same layout
            if self.mode == "pairs":
                self._engine_files.append(("protein_alignments", self.path("protein_alignments.fasta")))
                cmd += ["--protein-out", self.path("protein_alignments.fasta") + ".tmp"]
            else:
                d = self.outputs["protein_alignments"] = self.path("protein_alignments")
                os.makedirs(d, exist_ok=True)
                cmd += ["--protein-dir", d]
        return cmd, out, fail

    def run_engine(self, group_column: str | None, feed) -> dict:
        """Start `coal-cpp codon-dnds`, let feed(write) send it the groups through its standard input, and collect
        its failures. feed may return a dict of details for the report."""
        cmd, out, fail = self.engine_cmd(group_column)
        logp = self.wpath("logs", "codon-dnds.log")
        with open(logp, "wb") as log:
            try:
                proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=log)
            except OSError as e:
                raise PipelineError(f"could not start {cmd[0]} ({e.strerror or e})") from None
            self._engine = proc
            detail, broken = {}, False
            try:
                try:
                    detail = feed(proc.stdin.write) or {}
                    proc.stdin.close()
                except (BrokenPipeError, OSError) as e:
                    if isinstance(e, BrokenPipeError) or getattr(e, "errno", None) in (22, 32):
                        broken = True   # coal-cpp ended early: its message is in the log
                    else:
                        raise
                if proc.poll() is None and not broken:
                    try:
                        proc.stdin.close()
                    except OSError:
                        pass
                rc = proc.wait()
            except BaseException:
                # An interrupt (Ctrl-C, SIGTERM) or an error while feeding: the engine must not be left running.
                # Closing its input would only make it finish the pairs it has and write its files, as an orphan
                # of a run that no longer exists; it is killed, and what it wrote so far is removed.
                self.kill_children()
                self.remove_temporaries()
                raise
            finally:
                self._engine = None
        if rc != 0 or broken:
            self.remove_temporaries()
            raise PipelineError(f"coal-cpp codon-dnds failed ({_msa.exit_text(rc)}); the end of {logp}:\n"
                                f"{_tail(logp)}")
        for key, final in reversed(self._engine_files):   # dnds.tsv last: with it the stage's files are all there
            os.replace(final + ".tmp", final)
            self.outputs[key] = final
        with open(fail, encoding="utf-8", errors="backslashreplace") as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4 or (not f[1] and self.mode == "pairs"):
                    continue   # a failed pair is reported from its row of dnds.tsv (_count_pairs)
                self.fail("gene" if f[1] else "cluster", f[2], f[3],
                          cluster=f[0] if self.mode == "clusters" else "", seq1=f[1])
        for line in _tail(logp, 2000).splitlines():
            if "timing:" in line:
                w = line.replace(",", " ").split()
                try:
                    detail["align_codon_seconds"] = round(float(w[w.index("align+codon") + 1]), 3)
                    detail["kaks_seconds"] = round(float(w[w.index("kaks") + 1]), 3)
                except (ValueError, IndexError):
                    pass
        detail["commands"] = [[str(c) for c in cmd]]   # (the groups arrive on its standard input)
        return detail

    # -- the end
    def finish(self, exit_code: int, error: str | None = None, status: str | None = None) -> dict:
        """Close failed.tsv and write the report. status: "completed" (exit 0), "completed_with_failures" (2),
        "failed" (1: the run could not be done) or "interrupted" (a signal or Ctrl-C ended it: the tables in the
        folder are incomplete or missing)."""
        if self._failed_fh:
            fh, self._failed_fh = self._failed_fh, None
            try:
                fh.close()
                os.replace(self.outputs["failed"] + ".tmp", self.outputs["failed"])
            except OSError as e:   # a full disk shows when the buffer is written out
                why = f"cannot write {self.outputs.get('failed', 'failed.tsv')}: {e.strerror or e}"
                self.outputs.pop("failed", None)   # not complete: it does not get its name (removed below)
                if exit_code in (0, 2):   # the run would have looked complete with an incomplete failures table
                    print(f"coal pipeline: error: {why}", file=sys.stderr)
                    exit_code, error, status = 1, why, "failed"
        self.remove_temporaries()
        a = self.a
        skip = {"out_dir", "quiet", "report", "resume", "install_tools"}
        if status is None:
            status = {0: "completed", 2: "completed_with_failures"}.get(exit_code, "failed")
        # every tool's command lines, from the stages (also the reused ones); the aligners run once per cluster
        # or pair, so they have the template of that command
        tools = {name: dict(t) for name, t in self.tools.items()}
        for st in self.stages:
            for cmd in st.get("commands", ()):
                for t in tools.values():
                    if cmd and cmd[0] == t.get("path"):
                        t.setdefault("commands", []).append(cmd)
        for name in ("clustalo", "mafft"):
            if name in tools:
                threads = 1 if self.mode == "pairs" else getattr(a, "msa_threads", 1)
                tools[name]["command_template"] = _msa.command_line(name, tools[name]["path"], threads)
        report = {
            "schema": SCHEMA, "coal_version": __version__, "mode": self.mode, "command": self.argv,
            "out_dir": self.out, "exit_code": exit_code, "status": status,
            "started": self.started, "finished": _utc_now(),
            "platform": {"system": platform.system(), "release": platform.release(), "machine": platform.machine(),
                         "python": platform.python_version(), "python_implementation": platform.python_implementation(),
                         "cpus": os.cpu_count()},
            "threads": {"requested": a.jobs, "engine": self.jobs, "tools": self.threads, "usable_cpus": self.cpus},
            "inputs": self.input_records,
            "parameters": {k: v for k, v in sorted(vars(a).items()) if k not in skip},
            "tools": tools, "stages": self.stages,
            "seconds": round(time.time() - self.t0, 3), "counts": self.counts,
            "failures_total": self.n_failures, "failures": self.failures,
            "failures_truncated": self.n_failures > len(self.failures),
            "outputs": {k: os.path.relpath(v, self.out) for k, v in sorted(self.outputs.items())},
        }
        from . import citations
        report["cite"] = citations.keys_for(
            methods=a.method or (), diamond="diamond" in self.tools,
            deepclust=self.mode == "clusters" and not a.clusters and "diamond" in self.tools,
            clustalo="clustalo" in self.tools, mafft="mafft" in self.tools)
        if error:
            report["error"] = error
        if self.isoforms:
            report["isoforms"] = {"longest_isoform": bool(a.longest_isoform), "gene_map": a.gene_map or None,
                                  "isoform_rule": getattr(a, "isoform_rule", "auto"),
                                  "sets": self.isoforms}
        if self.code_hints:
            report["genetic_code_hints"] = self.code_hints
        if self.notes:
            report["warnings"] = self.notes
        rp = a.report or self.path("report.json")

        def write(obj) -> None:
            # through a temporary file: a report is complete or absent, never cut short (a full disk)
            with _atomic(rp) as fh:
                json.dump(obj, fh, indent=1, default=str)
                fh.write("\n")
        try:
            os.makedirs(os.path.dirname(os.path.abspath(rp)), exist_ok=True)
            write(report)
            self.outputs["report"] = rp
        except OSError as e:
            print(f"coal pipeline: error: could not write the report {rp}: {e.strerror or e}", file=sys.stderr)
            if report["exit_code"] == 0:   # (as `coal --report`: a run whose report is missing is not a clean run)
                report["exit_code"] = 1
                report["status"] = "failed"
                report["error"] = f"could not write the report {rp}: {e.strerror or e}"
            try:   # the few lines that say what happened may still fit where the whole report did not
                write({k: report[k] for k in ("schema", "coal_version", "mode", "out_dir", "exit_code", "status",
                                               "error") if k in report} | {"report_truncated": True})
                self.outputs["report"] = rp
            except OSError:
                pass
        return report


def _count_pairs(run: _Run, dnds_path: str, group_col: bool) -> None:
    """Pair counts by status from dnds.tsv (the first method's rows), and the error pairs into failed.tsv."""
    c = {"total": 0, "ok": 0, "warning": 0, "not_converged": 0, "saturated": 0, "no_data": 0, "error": 0}
    k = 1 if group_col else 0
    first = None
    # the notes of --min-codons / --max-ds as the engine writes them (cpp/src/kaks_cli.cpp PairFlags)
    a = getattr(run, "a", None)
    min_codons, max_ds = getattr(a, "min_codons", 0), getattr(a, "max_ds", None)
    short = f"fewer than {min_codons} codons compared" if min_codons else None
    high = f"dS above {float(max_ds):.6g}" if max_ds is not None else None
    min_identity = getattr(a, "min_identity", None)
    low = f"protein identity below {float(min_identity):.6g}%" if min_identity is not None else None
    flagged = {"min_codons": 0, "max_ds": 0, "min_identity": 0}
    with open(dnds_path, encoding="utf-8", errors="backslashreplace") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if first is None:
                first = f[k + 2]
            if f[k + 2] != first:
                continue
            c["total"] += 1
            c[f[-2]] = c.get(f[-2], 0) + 1
            if short or high or low:
                notes = f[-1].split("; ")
                flagged["min_codons"] += short in notes
                flagged["max_ds"] += high in notes
                flagged["min_identity"] += low in notes
            if f[-2] == "error":
                run.fail("pair", "pair_failed", f[-1], cluster=f[0] if group_col else "", seq1=f[k], seq2=f[k + 1])
    run.counts["pairs"] = c
    if short or high or low:
        run.counts["flagged"] = {key: n for key, n, on in (("min_codons", flagged["min_codons"], short),
                                                           ("max_ds", flagged["max_ds"], high),
                                                           ("min_identity", flagged["min_identity"], low)) if on}


# ------------------------------------------------------------------------------------------------ pair mode

_HIT_COLS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qstart", "qend", "sstart", "send"]


def _best_hits(path: str, keep: int) -> dict[str, list[str]]:
    """{query: the fields after qseqid of its first line} (DIAMOND writes a query's best hit first)."""
    best: dict[str, list[str]] = {}
    with open(path, encoding="utf-8", errors="backslashreplace") as fh:
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            if len(f) >= 2 and f[0] not in best:
                best[f[0]] = f[1:keep]
    return best


def _pairs_mode(run: _Run) -> None:
    a = run.a
    for p in (a.query, a.subject):
        if not os.path.isfile(p):
            raise PipelineError(f"cannot open {p!r}: no such file")
    q_faa, q_genes, _ = run.translate("query", [a.query])
    s_faa, s_genes, _ = run.translate("subject", [a.subject])
    run.outputs["genes"] = run.path("genes.tsv")
    q_order: dict[str, int] = {}
    s_order: dict[str, int] = {}
    q_drop = run.find_isoforms("query", [a.query], q_genes)
    s_drop = run.find_isoforms("subject", [a.subject], s_genes)
    with _atomic(run.outputs["genes"]) as gf:
        with open(q_genes, encoding="utf-8", errors="backslashreplace") as fh:
            gf.write("set\t" + fh.readline())
        q_failed = run.read_genes(q_genes, "query", gf, q_order, q_drop)
        s_failed = run.read_genes(s_genes, "subject", gf, s_order, s_drop)
    run.say(f"translated {run.n_records('query'):,} query and {run.n_records('subject'):,} subject records "
            f"({run.n_records('query', 'failed') + run.n_records('subject', 'failed'):,} not usable)")
    run.require_usable_genes("query", [a.query])
    run.require_usable_genes("subject", [a.subject])
    run.code_hint("query", [run.plain(a.query)], q_genes, [a.query])
    run.code_hint("subject", [run.plain(a.subject)], s_genes, [a.subject])
    run.isoform_warning()
    if q_drop or s_drop:   # DIAMOND searches the kept records only; a pair that names a dropped one says why
        q_faa, s_faa = run.kept_proteins(q_faa, q_drop, "query"), run.kept_proteins(s_faa, s_drop, "subject")
        for failed, drop in ((q_failed, q_drop), (s_failed, s_drop)):
            for rid, (gene, kept) in drop.items():
                failed.setdefault(rid, f"dropped by --longest-isoform: {kept} is the longest record of gene {gene}")

    orth = run.outputs["orthologs"] = run.path("orthologs.tsv")
    hits: dict[str, list[str]] = {}
    if a.pairs:
        if a.aligner == "diamond":
            raise PipelineError("--aligner diamond uses the alignments of the DIAMOND search, so it cannot be "
                                "combined with --pairs")
        pairs = _read_pairs_file(run.plain(a.pairs), a.pairs)
        run.stages.append({"name": "orthologs", "seconds": 0.0, "cached": False, "source": "given",
                           "pairs": len(pairs)})
        if not run.pairs_in_place:   # (--pairs DIR/orthologs.tsv: the input is the table; it is not rewritten)
            with _atomic(orth) as fh:
                fh.write("query\tsubject\n")
                for q, s in pairs:
                    fh.write(f"{q}\t{s}\n")
    else:
        pairs, hits = _search(run, q_faa, s_faa, q_order, s_order, orth)
    run.say(f"{len(pairs):,} pairs ({'given' if a.pairs else a.ortho})")
    run.counts["orthologs"] = len(pairs)

    q_idx, s_idx = FastaIndex(run.plain(a.query)), FastaIndex(run.plain(a.subject))
    if a.pairs:
        every = bool(pairs) and all((q not in q_idx and q not in q_failed) or (s not in s_idx and s not in s_failed)
                                    for q, s in pairs)
        _unmatched_ids_warning(run, "--pairs", a.pairs,
                               [("query", [q for q, _ in pairs], q_idx, q_failed),
                                ("subject", [s for _, s in pairs], s_idx, s_failed)], every)
    need_prot = a.aligner in ("clustalo", "mafft")
    exe = run.msa_tool(a.aligner) if need_prot else None
    qp_idx = FastaIndex(q_faa) if need_prot else None
    sp_idx = FastaIndex(s_faa) if need_prot else None

    def member(idx: FastaIndex, failed: dict, rid: str, which: str, row: bytes = b"") -> bytes:
        if rid in failed:
            return b"F\t" + rid.encode() + b"\t" + _clean(failed[rid]).encode() + b"\n"
        if rid not in idx:
            return b"F\t" + rid.encode() + f"\tno CDS record with this ID in the {which} file\n".encode()
        return b"M\t" + rid.encode() + b"\t" + idx.get(rid) + (b"\t" + row if row else b"") + b"\n"

    def feed(write) -> dict:
        t_msa = 0.0
        if not need_prot:
            for q, s in pairs:
                buf = [b"G\t", q.encode(), b"&", s.encode(), b"\tall\n", member(q_idx, q_failed, q, "query"),
                       member(s_idx, s_failed, s, "subject")]
                h = hits.get(q) if a.aligner == "diamond" else None
                if h is not None and len(h) >= 11:
                    buf.append(("H\t" + "\t".join(h[5:11]) + "\n").encode())
                buf.append(b"E\n")
                write(b"".join(buf))
            return {}
        step = 64 * run.jobs
        for at in range(0, len(pairs), step):
            chunk = pairs[at:at + step]
            todo = {}
            for k, (q, s) in enumerate(chunk):
                if q in q_failed or s in s_failed or q not in qp_idx or s not in sp_idx:
                    continue
                todo[k] = [(q, qp_idx.get(q).decode()), (s, sp_idx.get(s).decode())]
            t = time.perf_counter()
            res = _msa.align_clusters(todo, tool=a.aligner, jobs=run.threads, threads=1, pairs_with="msa", exe=exe,
                                      timeout=a.msa_timeout)
            t_msa += time.perf_counter() - t
            for k, (q, s) in enumerate(chunk):
                r = res.get(k)
                buf = [b"G\t", q.encode(), b"&", s.encode(), b"\tall\n"]
                if r is not None and r.error:
                    buf.append(("X\tmsa_failed\t" + _clean(f"{a.aligner} failed: {r.error}") + "\n").encode())
                rows = [x[1].encode() for x in r.rows] if r is not None and r.rows else [b"", b""]
                buf += [member(q_idx, q_failed, q, "query", rows[0]), member(s_idx, s_failed, s, "subject", rows[1]),
                        b"E\n"]
                write(b"".join(buf))
        return {"msa_seconds": round(t_msa, 3)}

    t = time.perf_counter()
    detail = run.run_engine(None, feed)
    run.stages.append({"name": "align_codon_kaks", "seconds": round(time.perf_counter() - t, 3), "cached": False,
                       "aligner": a.aligner, **detail})
    for i in (q_idx, s_idx, qp_idx, sp_idx):
        if i is not None:
            i.close()
    _count_pairs(run, run.outputs["dnds"], False)


def _unmatched_ids_warning(run: _Run, option: str, path: str, sides, every_group_fails: bool = False) -> None:
    """One line when half or more of the IDs of a --pairs / --clusters file on one side are no records of the
    sequence files, or when every pair names at least one unknown ID (every_group_fails):
    the usual causes are a version suffix on one side (AT1G01010 against AT1G01010.1), another ID type (gene
    against transcript IDs), or swapped columns. Every pair still gets its own row and its failure; this line
    says what is going on. sides: [(name, the IDs asked for, the index of the file, the failed genes)]."""
    for name, ids, idx, failed in sides:
        if not ids:
            continue
        missing = [i for i in ids if i not in idx and i not in failed]
        # half of a side unknown, or no pair at all that can be made: something systematic, worth one line
        # (below that, single missing genes are reported one by one in failed.tsv)
        if not missing or (len(missing) * 2 < len(ids) and not every_group_fails):
            continue
        have = next(iter(idx.pos), None)
        ex = missing[0]
        hint = ""
        if have is not None:
            stem = ex.split(".")[0]
            if any(h == stem or h.startswith(ex + ".") or h.split(".")[0] == stem for h in list(idx.pos)[:5000]):
                hint = (" The IDs differ by a version or transcript suffix (.1, .2): use the IDs exactly as "
                        "they are after '>' in the sequence file (the first word of the header)")
            else:
                hint = (" Check that the IDs are the first word of the FASTA headers, that the columns are in "
                        "the right order, and whether one side has a version suffix (.1) the other lacks")
        run.warn(f"{len(missing):,} of the {len(ids):,} {name} IDs of {option} {path} are not in the {name} "
                 f"sequence file(s): e.g. {ex!r} in {os.path.basename(path)}"
                 + (f", while the sequences are named like {have!r}." if have is not None else ".") + hint)


def _read_pairs_file(path: str, shown: str | None = None) -> list[tuple[str, str]]:
    """The pairs of a --pairs file. shown: the name for messages (the file as given, when path is its
    decompressed copy)."""
    pairs = []
    try:
        fh = open(path, encoding="utf-8", errors="backslashreplace")
    except OSError as e:
        raise PipelineError(f"cannot open {shown or path!r}: {e.strerror or e}") from None
    path = shown or path
    with fh:
        for n, line in enumerate(fh):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t") if "\t" in line else line.split()
            if len(f) < 2:
                raise PipelineError(f"{path}: line {n + 1}: expected two columns (query ID, subject ID)")
            if n == 0 and (f[0].lower(), f[1].lower()) in (("query", "subject"), ("query_id", "subject_id"),
                                                             ("seq1", "seq2"), ("qseqid", "sseqid")):
                continue
            pairs.append((f[0].strip(), f[1].strip()))
    return pairs


def _search(run: _Run, q_faa: str, s_faa: str, q_order: dict, s_order: dict, orth: str):
    a = run.a
    dmd = run.diamond()
    gapped = a.aligner == "diamond"
    extra = shlex.split(a.diamond_args) if a.diamond_args else []
    common = ["--max-target-seqs", "1", "--evalue", str(a.evalue), "--threads", str(run.threads), "--quiet"]
    if a.sensitivity != "default":
        common.append("--" + a.sensitivity)
    if gapped:
        common += ["--masking", "0"]   # DIAMOND's alignment inside a masked stretch is arbitrary (docs/alignment.md)
    fwd, rev = run.wpath("query_vs_subject.tsv"), run.wpath("subject_vs_query.tsv")
    both = a.ortho == "rbh"

    def go():
        run.run_tool([dmd, "makedb", "--in", s_faa, "-d", run.wpath("subject"), "--threads", str(run.threads),
                      "--quiet"], "diamond", "diamond makedb")
        cols = _HIT_COLS + (["qseq_gapped", "sseq_gapped"] if gapped else [])
        run.run_tool([dmd, "blastp", "-q", q_faa, "-d", run.wpath("subject"), "-o", fwd, "--outfmt", "6", *cols,
                      *common, *extra], "diamond", "diamond blastp (query against subject)")
        if both:
            run.run_tool([dmd, "makedb", "--in", q_faa, "-d", run.wpath("query"), "--threads", str(run.threads),
                          "--quiet"], "diamond", "diamond makedb")
            run.run_tool([dmd, "blastp", "-q", s_faa, "-d", run.wpath("query"), "-o", rev, "--outfmt", "6",
                          *_HIT_COLS, *common, *extra], "diamond", "diamond blastp (subject against query)")
        for d in ("subject.dmnd", "query.dmnd"):
            if os.path.isfile(run.wpath(d)):
                os.remove(run.wpath(d))
    run.stage("diamond_search", _sig("search", _file_sig([q_faa, s_faa]), run.tools["diamond"]["version"], a.ortho,
                                     a.evalue, a.sensitivity, gapped, extra), [fwd] + ([rev] if both else []), go)
    t = time.perf_counter()
    hits = _best_hits(fwd, 12 if gapped else 10)
    back = _best_hits(rev, 2) if both else {}
    pairs = [(q, h[0]) for q, h in hits.items() if not both or back.get(h[0], [None])[0] == q]
    pairs.sort(key=lambda p: (q_order.get(p[0], len(q_order)), p[0]))
    with _atomic(orth) as fh:
        fh.write("query\tsubject\tpident\tlength\tevalue\tbitscore\n")
        for q, s in pairs:
            fh.write("\t".join([q, s, *hits[q][1:5]]) + "\n")
    # genes without a partner
    unp = run.outputs["unpaired"] = run.path("unpaired.tsv")
    paired_q = {q for q, _ in pairs}
    paired_s = {s for _, s in pairs}
    n = {"no_hit": 0, "not_reciprocal": 0}
    with _atomic(unp) as fh:
        fh.write("set\tid\treason\n")
        for q in q_order:
            if q not in paired_q:
                why = "no_hit" if q not in hits else "not_reciprocal"
                n[why] += 1
                fh.write(f"query\t{q}\t{why}\n")
        if both:
            for s in s_order:
                if s not in paired_s:
                    why = "no_hit" if s not in back else "not_reciprocal"
                    n[why] += 1
                    fh.write(f"subject\t{s}\t{why}\n")
    run.counts["unpaired"] = n
    run.stages.append({"name": "orthologs", "seconds": round(time.perf_counter() - t, 3), "cached": False,
                       "source": a.ortho, "pairs": len(pairs)})
    if not gapped:
        hits = {}
    return pairs, hits


# ------------------------------------------------------------------------------------------------ cluster mode

def read_clusters(path: str) -> tuple[dict[str, list[str]], int]:
    """DIAMOND's clusters table (representative, member; tab-separated) -> ({representative: [members]}, repeated
    member lines). Clusters in the order their representative first appears, the representative first, then the
    members in file order. A member named twice stays in its first cluster."""
    clusters: dict[str, list[str]] = {}
    seen: set[str] = set()
    dup = 0
    try:
        fh = open(path, encoding="utf-8", errors="backslashreplace")
    except OSError as e:
        raise PipelineError(f"cannot open {path!r}: {e.strerror or e}") from None
    with fh:
        for n, line in enumerate(fh):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 2:
                raise PipelineError(f"{path}: line {n + 1}: expected two tab-separated columns (representative, "
                                    "member), as `diamond deepclust` writes them")
            rep, mem = f[0].strip(), f[1].strip()
            if not clusters and (rep.lower(), mem.lower()) in _CLUSTER_HEADERS:
                continue
            members = clusters.get(rep)
            if members is None:
                members = clusters[rep] = []
                if rep not in seen:
                    seen.add(rep)
                    members.append(rep)
            if mem in seen:
                if mem != rep:
                    dup += 1
                continue
            seen.add(mem)
            members.append(mem)
    return clusters, dup


def sample_pairs(n: int, k: int, seed, cluster_id: str) -> list[tuple[int, int]]:
    """k distinct pairs (i < j) of n members, the same for the same (seed, cluster_id) whatever else is in the run;
    sorted. All pairs when k is at least n(n-1)/2."""
    total = n * (n - 1) // 2
    if k >= total:
        return [(i, j) for i in range(n) for j in range(i + 1, n)]
    rng = random.Random(f"{seed}:{cluster_id}")
    out = []
    for t in sorted(rng.sample(range(total), k)):
        # t-th pair in the order (0,1), (0,2), ..., (1,2), ...: row i starts at i*n - i*(i+1)/2
        i = int((2 * n - 1 - math.sqrt((2 * n - 1) ** 2 - 8 * t)) / 2)
        while i * n - i * (i + 1) // 2 > t:
            i -= 1
        while (i + 1) * n - (i + 1) * (i + 2) // 2 <= t:
            i += 1
        out.append((i, i + 1 + t - (i * n - i * (i + 1) // 2)))
    return out


def _clusters_mode(run: _Run) -> None:
    a = run.a
    run.check_cds_ids()
    faa, genes, _ = run.translate("proteins", list(a.cds))
    run.outputs["genes"] = run.path("genes.tsv")
    order: dict[str, int] = {}
    drop = run.find_isoforms("cds", list(a.cds), genes)
    with _atomic(run.outputs["genes"]) as gf:
        with open(genes, encoding="utf-8", errors="backslashreplace") as fh:
            gf.write("set\t" + fh.readline())
        failed = run.read_genes(genes, "cds", gf, order, drop)
    run.say(f"translated {run.n_records('cds'):,} records ({run.n_records('cds', 'failed'):,} not usable)")
    run.require_usable_genes("cds", list(a.cds))
    run.code_hint("cds", [run.plain(p) for p in a.cds], genes, list(a.cds))
    run.isoform_warning()
    all_faa = faa
    if drop:   # DeepClust clusters the kept records only; a given cluster that names a dropped one says why
        faa = run.kept_proteins(faa, drop, "proteins")
        for rid, (gene, kept) in drop.items():
            failed.setdefault(rid, f"dropped by --longest-isoform: {kept} is the longest record of gene {gene}")

    if a.clusters:
        cfile = run.plain(a.clusters)
        run.stages.append({"name": "clusters", "seconds": 0.0, "cached": False, "source": "given"})
    else:
        cfile = run.outputs["clusters"] = run.path("clusters.tsv")
        dmd = run.diamond()
        extra = shlex.split(a.deepclust_args) if a.deepclust_args else []
        opts = []
        for flag, v in (("--approx-id", a.approx_id), ("--member-cover", a.member_cover),
                        ("--mutual-cover", a.mutual_cover), ("--evalue", a.cluster_evalue)):
            if v is not None:
                opts += [flag, str(v)]

        # DeepClust's result depends on the order of its input (which sequence becomes a representative, and
        # with that the clusters): it is given the proteins sorted by ID (bytes), so the clusters are the same
        # whatever the order of the --cds files and of the records in them (docs/pipeline.md "Reproducibility").
        sfaa = run.wpath("proteins.sorted.faa")

        def go():
            n = sort_fasta_by_id(faa, sfaa)
            try:   # (a leftover or a symbolic link of that name is not written through)
                os.remove(cfile + ".tmp")
            except OSError:
                pass
            run.run_tool([dmd, "deepclust", "--db", sfaa, "--out", cfile + ".tmp", "--threads", str(run.threads),
                          "--quiet", *opts, *extra], "diamond", "diamond deepclust")
            os.replace(cfile + ".tmp", cfile)
            return {"input_order": "sorted by ID", "sequences": n}
        run.stage("deepclust", _sig("deepclust", "sorted-by-id", _file_sig([faa]), run.tools["diamond"]["version"],
                                    opts, extra), [cfile, sfaa], go)
        if a.realign:
            rfile = run.outputs["realign"] = run.path("realign.tsv")
            rextra = shlex.split(a.realign_args) if a.realign_args else []

            def go2():
                try:
                    os.remove(rfile + ".body")
                except OSError:
                    pass
                run.run_tool([dmd, "realign", "--db", sfaa, "--clusters", cfile, "--out", rfile + ".body", "--threads",
                              str(run.threads), "--quiet", "--outfmt", "6", "qseqid", "sseqid", "approx_pident",
                              "evalue", "bitscore", "qcovhsp", "scovhsp", *rextra], "diamond", "diamond realign")
                with _atomic(rfile) as out, \
                        open(rfile + ".body", encoding="utf-8", errors="backslashreplace") as body:
                    # (the columns are fixed: --realign-args cannot change --outfmt)
                    out.write("member\trepresentative\tapprox_pident\tevalue\tbitscore\tqcovhsp\tscovhsp\n")
                    shutil.copyfileobj(body, out)
                os.remove(rfile + ".body")
            run.stage("realign", _sig("realign", _file_sig([cfile]), run.tools["diamond"]["version"], rextra),
                      [rfile], go2)
    t0 = time.perf_counter()
    clusters, dup = read_clusters(cfile)
    if dup:
        run.say(f"warning: {dup:,} member line(s) of {a.clusters or cfile} name a sequence that is already in a cluster; the "
                "first one is used")
    run.say(f"{len(clusters):,} clusters, {sum(len(m) for m in clusters.values()):,} sequences")

    cds_idx, prot_idx = FastaIndex([run.plain(p) for p in a.cds]), FastaIndex(all_faa)
    if a.clusters:
        _unmatched_ids_warning(run, "--clusters", a.clusters,
                               [("CDS", [m for mem in clusters.values() for m in mem], cds_idx, failed)])
    tool = a.aligner
    cap = a.max_cluster_size
    plans: list[dict] = []
    cc = run.counts["clusters"] = {"total": len(clusters), "singletons": 0, "all_pairs": 0, "above_max_size": 0,
                                   "skipped": 0, "msa_fallback": 0, "members_dropped": 0}
    msa_exe = [None]
    t_msa = [0.0]

    def usable(cid: str, members: list[str]) -> list[str]:
        ok = []
        for m in members:
            if m in failed:
                run.fail("gene", "gene_failed", f"left out of the cluster: {failed[m]}", cluster=cid, seq1=m)
            elif m not in cds_idx:
                run.fail("gene", "no_cds", "left out of the cluster: no CDS record with this ID", cluster=cid, seq1=m)
            else:
                ok.append(m)
        return ok

    def emit_pairwise(write, cid: str, ids: list[str], pairs) -> None:
        """Each pair as a two-member group of its own (aligned pairwise by the engine)."""
        cache: dict[int, bytes] = {}
        head = b"G\t" + cid.encode() + b"\tall\tnocodon\n"
        for i, j in pairs:
            if len(cache) > 256:   # a bounded cache of member lines (member 0, the representative, stays)
                keep = cache.get(0)
                cache.clear()
                if keep is not None:
                    cache[0] = keep
            for k in (i, j):
                if k not in cache:
                    cache[k] = b"M\t" + ids[k].encode() + b"\t" + cds_idx.get(ids[k]) + b"\n"
            write(b"".join((head, cache[i], cache[j], b"E\n")))

    def feed(write) -> dict:
        pending: list[tuple[dict, list[str]]] = []   # clusters waiting for their MSA, in order
        weight = 0   # members waiting

        def flush():
            nonlocal weight
            todo = {k: [(m, prot_idx.get(m).decode()) for m in ids] for k, (p, ids) in enumerate(pending)
                    if p["aligner"] == tool}
            res = {}
            if todo:
                if msa_exe[0] is None:
                    msa_exe[0] = run.msa_tool(tool)
                t = time.perf_counter()
                res = _msa.align_clusters(todo, tool=tool, jobs=run.threads, threads=a.msa_threads, pairs_with="msa",
                                          exe=msa_exe[0], timeout=a.msa_timeout)
                t_msa[0] += time.perf_counter() - t
            for k, (p, ids) in enumerate(pending):
                cid = p["cluster_id"]
                r = res.get(k)
                if r is not None and r.error:
                    # one cluster's MSA failed: its pairs are aligned pairwise instead, and the summary says so
                    p["aligner"] = "nw_pairwise"
                    p["note"] = _clean(f"{tool} failed ({r.error}); pairs aligned pairwise")
                    cc["msa_fallback"] += 1
                    run.fail("cluster", "msa_failed", p["note"], cluster=cid)
                    emit_pairwise(write, cid, ids, ((i, j) for i in range(len(ids)) for j in range(i + 1, len(ids))))
                    continue
                if p["aligner"] == "nw_pairwise":
                    emit_pairwise(write, cid, ids, p.pop("_pairs"))
                    continue
                rows = [x[1].encode() for x in r.rows] if r is not None else [b""] * len(ids)
                buf = [b"G\t", cid.encode(), b"\tall\n"]
                for m, row in zip(ids, rows):
                    buf += [b"M\t", m.encode(), b"\t", cds_idx.get(m), b"\t" if row else b"", row, b"\n"]
                buf.append(b"E\n")
                write(b"".join(buf))
            pending.clear()
            weight = 0

        for cid, members in clusters.items():
            ids = usable(cid, members)
            n = len(ids)
            p = {"cluster_id": cid, "representative": cid, "n_members": len(members), "n_used": n,
                 "pair_selection": "all", "aligner": "none", "n_planned": 0, "note": ""}
            cc["members_dropped"] += len(members) - n
            plans.append(p)
            notes = []
            if len(members) != n:
                notes.append(f"{len(members) - n} member(s) left out (see failed.tsv)")
            if n < 2:
                p["pair_selection"] = "singleton" if len(members) < 2 else "none"
                cc["singletons"] += 1
            elif n > cap:
                cc["above_max_size"] += 1
                if a.large_clusters == "skip":
                    p["pair_selection"] = "skipped"
                    cc["skipped"] += 1
                    notes.append(f"{n} members: above --max-cluster-size {cap}, not computed")
                else:
                    p["aligner"] = "nw_pairwise"
                    if a.large_clusters == "representative":
                        p["pair_selection"] = "representative"
                        p["_pairs"] = [(0, j) for j in range(1, n)]
                        if ids[0] != cid:
                            notes.append(f"the representative is not usable; {ids[0]} stands in for it")
                    else:
                        p["_pairs"] = sample_pairs(n, a.pair_sample, a.seed, cid)
                        p["pair_selection"] = f"sample:{len(p['_pairs'])}"
                    p["n_planned"] = len(p["_pairs"])
                    notes.append(f"{n} members: above --max-cluster-size {cap}")
            else:
                cc["all_pairs"] += 1
                p["n_planned"] = n * (n - 1) // 2
                p["aligner"] = "nw" if n == 2 else tool
            p["note"] = "; ".join(notes)
            if p["n_planned"]:
                pending.append((p, ids))
                weight += n
                if len(pending) >= 64 * run.jobs or weight >= 2000 * run.jobs:
                    flush()
        flush()
        return {}

    detail = run.run_engine("cluster_id", feed)
    cds_idx.close()
    prot_idx.close()
    run.stages.append({"name": "align_codon_kaks", "seconds": round(time.perf_counter() - t0, 3), "cached": False,
                       "aligner": tool, "msa_seconds": round(t_msa[0], 3), **detail})
    t = time.perf_counter()
    _count_pairs(run, run.outputs["dnds"], True)
    run.outputs["cluster_summary"] = run.path("cluster_summary.tsv")
    summarise(run.outputs["dnds"], plans, run.outputs["cluster_summary"])
    run.stages.append({"name": "cluster_summary", "seconds": round(time.perf_counter() - t, 3), "cached": False})


def summarise(dnds_path: str, plans: list[dict], out_path: str) -> None:
    """out_dir/cluster_summary.tsv from out_dir/dnds.tsv: one row per cluster and method (SUMMARY_COLUMNS).

    plans: one dict per cluster, in the order of dnds.tsv: cluster_id, representative, n_members, n_used,
    pair_selection, aligner, note. n_ok counts the pairs with status ok or warning (both have estimates); a pair
    whose YN / MYN iteration did not converge has none and is counted in n_not_converged, not in n_ok; the
    statistics are over those pairs' values that are numbers (dN/dS is NA when dS is 0: n_dNdS says how many
    values the dN/dS statistics are over)."""
    with open(dnds_path, encoding="utf-8", errors="backslashreplace") as fh, _atomic(out_path) as out:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: k for k, name in enumerate(header)}
        ci, mi, si = col["cluster_id"], col["method"], col["status"]
        vi = (col["dN"], col["dS"], col["dNdS"])
        out.write("\t".join(SUMMARY_COLUMNS) + "\n")
        methods: list[str] = []
        line = fh.readline()
        for p in plans:
            cid = p["cluster_id"]
            acc: dict[str, dict] = {}
            while line:
                f = line.rstrip("\n").split("\t")
                if f[ci] != cid:
                    break
                d = acc.get(f[mi])
                if d is None:
                    d = acc[f[mi]] = {"n": 0, "ok": 0, "not_converged": 0, "saturated": 0, "no_data": 0, "error": 0,
                                     "v": ([], [], [])}
                    if f[mi] not in methods:
                        methods.append(f[mi])
                d["n"] += 1
                st = f[si]
                if st in ("ok", "warning"):
                    d["ok"] += 1
                    for k in range(3):
                        x = f[vi[k]]
                        if x != "NA":
                            d["v"][k].append(float(x))
                else:
                    d[st] = d.get(st, 0) + 1
                line = fh.readline()
            lead = [cid, p["representative"], str(p["n_members"]), str(p["n_used"])]
            if not acc:
                status = "skipped" if p["pair_selection"] == "skipped" else "ok"
                out.write("\t".join(lead + ["NA", p["pair_selection"], p["aligner"], "0", "0", "0", "0", "0", "0"]
                                    + ["NA"] * 18 + ["0", status, p.get("note", "")]) + "\n")
                continue
            for m in (x for x in METHODS if x in acc):
                d = acc[m]
                status = "ok" if d["ok"] else "failed"
                out.write("\t".join(lead + [m, p["pair_selection"], p["aligner"], str(d["n"]), str(d["ok"]),
                                            str(d["not_converged"]), str(d["saturated"]), str(d["no_data"]),
                                            str(d["error"])]
                                    + _stats(d["v"][0]) + _stats(d["v"][1]) + _stats(d["v"][2])
                                    + [str(len(d["v"][2])), status, p.get("note", "")]) + "\n")


# ------------------------------------------------------------------------------------------------ low level

def codon_dnds(groups: Iterable[Sequence[Any]], out_tsv: _Path, method: Sequence[str] = ("YN",),
               codontable: int = 1, jobs: int = 0, group_column: str | None = None,
               codon_out: _Path | None = None, kaks_out: _Path | None = None, failures: _Path | None = None,
               extra_args: Sequence[str] = (), cpp: str | None = None) -> int:
    """Protein alignment -> codon alignment -> Ka/Ks for groups held in Python, through one `coal-cpp codon-dnds`
    process. groups: an iterable of (group_id, members[, pairs]) with members = [(id, cds)] or [(id, cds,
    protein_row)]: two members without rows are aligned here; three or more need rows. pairs: "all" (default),
    "rep" or a list of (i, j). Writes out_tsv (PAIR_COLUMNS, with group_column first when given); returns the
    number of groups sent. Raises PipelineError when coal-cpp is missing or fails."""
    cpp = cpp or _find_cpp()
    if isinstance(method, str):
        method = [method]
    cmd = [cpp, "codon-dnds", "-i", "-", "-o", os.fspath(out_tsv), "-c", str(codontable), "-j", str(jobs), "-q"]
    for m in method:
        cmd += ["-m", m]
    for flag, v in (("--group-column", group_column), ("--codon-out", codon_out), ("--kaks-out", kaks_out),
                    ("--failures", failures)):
        if v:
            cmd += [flag, os.fspath(v)]
    cmd += [str(x) for x in extra_args]
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    except OSError as e:
        raise PipelineError(f"could not start {cpp} ({e.strerror or e})") from None
    import threading
    err: list[bytes] = []
    reader = threading.Thread(target=lambda: err.append(proc.stderr.read()), daemon=True)
    reader.start()
    n = 0
    try:
        for g in groups:
            gid, members = g[0], g[1]
            pairs = g[2] if len(g) > 2 else "all"
            mode = pairs if isinstance(pairs, str) else "list"
            buf = [f"G\t{gid}\t{mode}\n"]
            for m in members:
                buf.append("M\t" + "\t".join("".join(str(x).split()) for x in m) + "\n")
            if mode == "list":
                buf += [f"P\t{i}\t{j}\n" for i, j in pairs]
            buf.append("E\n")
            proc.stdin.write("".join(buf).encode())
            n += 1
        proc.stdin.close()
    except (BrokenPipeError, OSError):
        pass
    rc = proc.wait()
    reader.join()
    proc.stderr.close()
    if rc != 0:
        raise PipelineError(f"coal-cpp codon-dnds failed ({_msa.exit_text(rc)}): "
                            f"{b''.join(err).decode('utf-8', 'backslashreplace').strip()[-800:]}")
    return n


# ------------------------------------------------------------------------------------------------ command line

def _positive(kind):
    def conv(text):
        v = kind(text)
        if v <= 0:
            raise argparse.ArgumentTypeError("must be greater than 0")
        return v
    return conv


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="coal pipeline",
        description="Coding sequences in, dN/dS table out: translation, DIAMOND (orthologs or clusters), protein "
                    "alignment, codon alignment and Ka/Ks in one step.",
        epilog="Details, defaults and the output columns: docs/pipeline.md.")
    sub = ap.add_subparsers(dest="mode", metavar="{pairs,clusters,install}")

    common = argparse.ArgumentParser(add_help=False)
    g = common.add_argument_group("Output")
    g.add_argument("-o", "--out-dir", required=True, metavar="DIR", help="output folder (created if missing)")
    g.add_argument("--kaks-format", action="store_true",
                   help="also write DIR/dnds.kaks in KaKs_Calculator's format")
    g.add_argument("--keep-alignments", action="store_true",
                   help="also write the codon alignments (pairs: DIR/codon_alignments.fasta; clusters: one file per "
                        "cluster in DIR/codon_alignments/)")
    g.add_argument("--keep-protein-alignments", action="store_true",
                   help="also write the protein alignments the codon alignments are made from (pairs: "
                        "DIR/protein_alignments.fasta; clusters: one file per cluster in DIR/protein_alignments/), "
                        "with the same records in the same order as --keep-alignments writes")
    g.add_argument("--report", metavar="FILE.json", help="the JSON run report (default: DIR/report.json)")
    g.add_argument("--hash-inputs", choices=("auto", "always", "never"), default="auto",
                   help="SHA-256 of the input files in the report: auto (default) hashes files up to 2 GiB and "
                        "records why a larger one was skipped; always; never")
    g.add_argument("--digits", type=int, default=8, metavar="N", help="significant digits (default 8)")
    g.add_argument("-q", "--quiet", action="store_true", help="no progress or summary on stderr")
    g = common.add_argument_group("Genes with several transcripts (every FASTA record is treated as a gene)")
    g.add_argument("--longest-isoform", action="store_true",
                   help="keep one record per gene, the longest CDS (ties: the ID that sorts first); the others are "
                        "listed in DIR/genes.tsv with status dropped. The gene of a record comes from --gene-map, "
                        "else from its header (locus_tag=, GeneID:, gene:, parent=FBgn, locus=, then [gene= / gene= "
                        "/ gene_symbol:) or a shared ID stem with a transcript suffix (AT1G01010.1 / .2, g1.t1, "
                        "_i1, -RA, -mRNA-1, _T001; docs/pipeline.md lists the styles)")
    g.add_argument("--gene-map", metavar="TSV",
                   help="record ID and gene ID, tab-separated, for --longest-isoform (and for the warning about "
                        "several transcripts per gene). It wins for the records it names; the genes of the others "
                        "are read from their headers or IDs as without it. A map that names no record of the "
                        "input, or gives a record two genes, is an error")
    g.add_argument("--isoform-rule", choices=("auto", "suffix", "map"), default="auto",
                   help="auto (default): header fields, else ID suffixes, and --longest-isoform refuses to drop "
                        "more than half of a file's records on the evidence of ID suffixes alone (HSP70.1 / "
                        "HSP70.2 may be paralogs, not transcripts); suffix: the ID suffixes are transcripts, drop "
                        "them whatever their share; map: use --gene-map only, read nothing from headers or IDs")
    g = common.add_argument_group("Ka/Ks")
    g.add_argument("-m", "--method", action="append", metavar="METHOD", type=_method_name, choices=METHODS,
                   help="YN (default), NG, LWL, MLWL, LPB, MLPB, MYN, Comeron (alias C95), Li; repeat for several "
                        "(Comeron and Li: genetic code 1 only)")
    g.add_argument("-c", "--codontable", type=int, default=1, choices=SUPPORTED, metavar="N",
                   help="NCBI genetic code for translation, the codon alignment and Ka/Ks (default 1)")
    g.add_argument("--yn-flavour", "--yn-flavor", default="kaks", choices=("kaks", "paml"),
                   help="YN00 as KaKs_Calculator (default) or as PAML yn00")
    g.add_argument("--yn-weighting", type=int, default=1, choices=(0, 1), help="YN00 pathway weighting (default 1)")
    g.add_argument("--comeron-flavour", "--comeron-flavor", default="coal", choices=("coal", "gestimator"),
                   help="Comeron with COAL's stop-codon and saturation policy (default) or with the numbers of "
                        "orthologr's gestimator")
    g.add_argument("--li-flavour", "--li-flavor", default="coal", choices=("coal", "seqinr"),
                   help="Li with COAL's stop-codon policy (default) or reading stop codons as seqinr::kaks does")
    g.add_argument("--terminal-stops", default="drop", choices=("drop", "compare"),
                   help="with --comeron-flavour gestimator / --li-flavour seqinr: compare takes the two terminal "
                        "stop codons of a pair as the last compared column (as orthologr 0.4.3); default: drop")
    g.add_argument("--min-codons", type=int, default=0, metavar="N",
                   help="flag estimates from fewer than N compared codons: status warning and the note \"fewer than "
                        "N codons compared\"; no row is dropped and no number changes (default 0: off)")
    g.add_argument("--max-ds", type=float, default=None, metavar="X",
                   help="flag estimates with dS above X (saturation): status warning and the note \"dS above X\" "
                        "(default: off; docs/pipeline.md \"Filtering before interpretation\")")
    g.add_argument("--min-identity", type=float, default=None, metavar="P",
                   help="flag pairs whose protein identity (the column prot_identity, percent of the alignment's "
                        "columns) is below P: status warning and the note \"protein identity below P%%\"; no row is "
                        "dropped and no number changes (default: off). A cluster can hold pairs that are not "
                        "related at all")
    g = common.add_argument_group("Pairwise protein alignment (COAL's global aligner)")
    g.add_argument("--preset", choices=("needle", "pwalign"),
                   help="needle: BLOSUM62, open 10, extend 0.5, end gaps free (the default); pwalign: open 14, "
                        "extend 4, end gaps penalised")
    g.add_argument("--matrix", default=None, help="BLOSUM62 (default), BLOSUM45, BLOSUM80, PAM250")
    g.add_argument("--gap-open", type=float, default=None, metavar="X", help="cost of a gap's first position")
    g.add_argument("--gap-extend", type=float, default=None, metavar="X", help="cost of each further position")
    g.add_argument("--end-gaps", choices=("free", "penalized"), default=None, help="default free")
    g.add_argument("--max-cells", type=int, default=4_000_000_000, metavar="N",
                   help="refuse a pair whose alignment table has more cells (default 4e9, about 2 GB)")
    g = common.add_argument_group("Running")
    g.add_argument("-j", "--jobs", type=int, default=0, metavar="N",
                   help="threads, and aligner processes at a time (default 0 = one per CPU this process may use; "
                        "never more than that)")
    g.add_argument("--resume", action="store_true",
                   help="reuse the translation and DIAMOND results of an earlier run in DIR when the inputs and "
                        "options are the same")
    g.add_argument("--diamond", metavar="PATH", help="the DIAMOND to run (default: PATH, $COAL_DIAMOND, COAL's cache)")
    g.add_argument("--install-tools", action="store_true",
                   help="install DIAMOND / Clustal Omega / MAFFT into COAL's cache when they are not found")
    g.add_argument("--msa-timeout", type=_positive(float), default=None, metavar="SECONDS",
                   help="give up on one Clustal Omega / MAFFT run after this long")

    p = sub.add_parser("pairs", parents=[common], help="orthologs between two species",
                       description="Orthologous pairs between two sets of coding sequences: translate, DIAMOND "
                                   "reciprocal best hits, pairwise protein alignment, codon alignment, Ka/Ks.")
    g = p.add_argument_group("Pairs")
    g.add_argument("--query", required=True, metavar="CDS.fasta", help="coding sequences of the first species")
    g.add_argument("--subject", required=True, metavar="CDS.fasta", help="coding sequences of the second species")
    g.add_argument("--ortho", choices=("rbh", "bh"), default="rbh",
                   help="rbh: reciprocal best hits (default); bh: every query's best hit")
    g.add_argument("--pairs", metavar="TSV", help="take the pairs from this file (query ID, subject ID) and skip "
                                                  "DIAMOND")
    g.add_argument("--aligner", choices=("nw", "clustalo", "mafft", "diamond"), default="nw",
                   help="nw: COAL's global aligner (default); clustalo / mafft: one process per pair (slow); "
                        "diamond: DIAMOND's own alignment, completed to full length (no realignment)")
    g.add_argument("--flanks", choices=("align", "gap"), default="align",
                   help="--aligner diamond: align the ends DIAMOND left out (default), or put them against gaps")
    g.add_argument("--evalue", type=float, default=1e-5, help="DIAMOND e-value cutoff (default 1e-5)")
    g.add_argument("--sensitivity", choices=SENSITIVITIES, default="default",
                   help="DIAMOND sensitivity mode (default: DIAMOND's default)")
    g.add_argument("--diamond-args", metavar="'ARGS'",
                   help="further arguments for `diamond blastp`, in one string (e.g. '--block-size 4 --tmpdir /scratch'); "
                        "the options COAL sets itself (--query, --db, --out, --outfmt, --threads, --max-target-seqs, "
                        "--evalue, a sensitivity) are refused")

    p = sub.add_parser("clusters", parents=[common], help="gene families: every pair within each cluster",
                       description="Clusters of coding sequences: translate, DIAMOND DeepClust (or a given clusters "
                                   "file), a protein MSA per cluster, codon MSA, Ka/Ks for every pair of a cluster, "
                                   "and a per-cluster summary.")
    g = p.add_argument_group("Clusters")
    g.add_argument("--cds", required=True, nargs="+", metavar="CDS.fasta",
                   help="coding sequences (one or more files; an ID must not occur in two of them)")
    g.add_argument("--prefix-ids", action="store_true",
                   help="put each --cds file's name (without folder and suffix) and '|' before its IDs: S_cer.fa's "
                        "g1.t1 becomes S_cer|g1.t1 in every output. For files whose IDs repeat across species "
                        "(BRAKER / AUGUSTUS g1.t1), which are refused without it")
    g.add_argument("--clusters", metavar="TSV", help="take the clusters from this file (representative, member; as "
                                                     "`diamond deepclust` writes them) and skip DIAMOND")
    g.add_argument("--aligner", choices=("clustalo", "mafft"), default="clustalo",
                   help="MSA program for clusters of 3 or more (default clustalo); two members are aligned by COAL")
    g.add_argument("--msa-threads", type=_positive(int), default=1, metavar="N",
                   help="threads of each MSA process (default 1; -j processes run at a time)")
    g.add_argument("--max-cluster-size", type=_positive(int), default=200, metavar="N",
                   help="clusters with more members are handled as --large-clusters says (default 200)")
    g.add_argument("--large-clusters", choices=("representative", "sample", "skip"), default="representative",
                   help="above --max-cluster-size: the representative against every member (default), --pair-sample "
                        "random pairs, or nothing; these pairs are aligned pairwise, without an MSA")
    g.add_argument("--pair-sample", type=_positive(int), default=1000, metavar="K",
                   help="--large-clusters sample: pairs per cluster (default 1000)")
    g.add_argument("--seed", type=int, default=1, help="seed of the pair sample (default 1)")
    g.add_argument("--approx-id", type=float, default=None, metavar="PCT", help="DeepClust --approx-id")
    g.add_argument("--member-cover", type=float, default=None, metavar="PCT",
                   help="DeepClust --member-cover (DIAMOND's default: 80)")
    g.add_argument("--mutual-cover", type=float, default=None, metavar="PCT", help="DeepClust --mutual-cover")
    g.add_argument("--cluster-evalue", type=float, default=None, metavar="E", help="DeepClust --evalue")
    g.add_argument("--deepclust-args", metavar="'ARGS'",
                   help="further arguments for `diamond deepclust`, in one string (e.g. '--memory-limit 12G'); the "
                        "options COAL sets itself (--db, --out, --threads) are refused")
    g.add_argument("--realign", action="store_true",
                   help="also run `diamond realign` (members against their representatives) into DIR/realign.tsv")
    g.add_argument("--realign-args", metavar="'ARGS'",
                   help="with --realign: further arguments for `diamond realign`, in one string; the options COAL "
                        "sets itself (--db, --clusters, --out, --outfmt, --threads) are refused")

    p = sub.add_parser("install", help="install DIAMOND, Clustal Omega or MAFFT into COAL's cache",
                       description="Install the external tools into COAL's cache folder (they are not bundled: "
                                   "docs/alignment.md \"Licences\").")
    p.add_argument("tools", nargs="+", choices=("diamond", "clustalo", "mafft"))
    p.add_argument("--force", action="store_true", help="install again even if one is there")
    p.add_argument("--sha256", metavar="HEX",
                   help="the SHA-256 the downloaded file must have (one tool at a time), for a download COAL has "
                        "no pinned hash for; the file is deleted, not run, if it differs. No download needs it "
                        "at present (DIAMOND and micromamba are pinned); it does not apply to clustalo or mafft")
    p.add_argument("--allow-unverified", action="store_true",
                   help="accept a download that nobody can verify (no pinned hash, no --sha256): it is kept and "
                        "run as it arrives, and its SHA-256 is printed. Without this such a download is refused. "
                        "It does not apply to clustalo or mafft (Bioconda only)")
    p.add_argument("-q", "--quiet", action="store_true")
    return ap


def _method_code_error(a: argparse.Namespace) -> str | None:
    """Comeron and Li exist for the standard genetic code only (their references have no other)."""
    ms = [_method_name(m) for m in (a.method or [])]
    for m in ("Comeron", "Li"):
        if m in ms and a.codontable != 1:
            return f"method {m} is implemented for the standard genetic code only (codontable 1)"
    return None


def _flag_error(a: argparse.Namespace) -> str | None:
    if a.min_codons is None or isinstance(a.min_codons, bool) or not isinstance(a.min_codons, int) \
            or a.min_codons < 0:
        return "--min-codons: expected a whole number of codons, 0 or more (0 = off)"
    if a.max_ds is not None and not (isinstance(a.max_ds, (int, float)) and not isinstance(a.max_ds, bool)
                                     and math.isfinite(a.max_ds) and a.max_ds > 0):
        return "--max-ds: expected a number greater than 0"
    mi = getattr(a, "min_identity", None)
    if mi is not None and not (isinstance(mi, (int, float)) and not isinstance(mi, bool) and math.isfinite(mi)
                               and 0 < mi <= 100):
        return "--min-identity: expected a percentage, above 0 and at most 100"
    return None


# --diamond-args / --deepclust-args / --realign-args: the DIAMOND options COAL sets itself. A string that
# repeats one of them is refused before anything runs: DIAMOND would take the user's value (its output would land
# elsewhere, or in another format, and the run would end with an unrelated error) or the report would describe a
# search that was not the one run. The lists follow the command lines of _search() and _clusters_mode();
# tests/test_pipeline_robustness.py compares them with the argv a run records in its report.
_DIAMOND_ALIASES = {"-d": "--db", "-o": "--out", "-q": "--query", "-p": "--threads", "-f": "--outfmt",
                    "-k": "--max-target-seqs", "-e": "--evalue", "-v": "--verbose"}
_OWN_DIAMOND_OPTIONS: dict[str, dict[str, str]] = {
    "diamond_args": {      # diamond blastp (and makedb)
        "--query": "the translated query set", "--db": "the database COAL builds from the other set",
        "--in": "the proteins COAL translated", "--out": "COAL reads the hits from its work folder",
        "--outfmt": "COAL reads fixed columns", "--threads": "use -j", "--quiet": "",
        "--max-target-seqs": "1: the best hit of each query is the rule of --ortho rbh / bh",
        "--evalue": "use --evalue", "--masking": "set by --aligner diamond",
    },
    "deepclust_args": {    # diamond deepclust
        "--db": "the translated proteins, sorted by ID", "--out": "DIR/clusters.tsv", "--threads": "use -j",
        "--quiet": "", "--approx-id": "use --approx-id", "--member-cover": "use --member-cover",
        "--mutual-cover": "use --mutual-cover", "--evalue": "use --cluster-evalue",
    },
    "realign_args": {      # diamond realign
        "--db": "the translated proteins", "--clusters": "the clusters of this run", "--out": "DIR/realign.tsv",
        "--outfmt": "the columns of DIR/realign.tsv are fixed", "--threads": "use -j", "--quiet": "",
    },
}
# options COAL sets only when its own option is given: the same one in the string is then a second value
_OWN_ONLY_IF = {"--approx-id": "approx_id", "--member-cover": "member_cover", "--mutual-cover": "mutual_cover",
                ("deepclust_args", "--evalue"): "cluster_evalue", ("diamond_args", "--masking"): "aligner"}


def _own_option(flag: str) -> str:
    """The long form of a DIAMOND option as written in an ARGS string (-o, --out, --out=FILE -> --out)."""
    flag = flag.split("=", 1)[0]
    return _DIAMOND_ALIASES.get(flag, flag)


def _args_error(a: argparse.Namespace) -> str | None:
    """Why --diamond-args / --deepclust-args / --realign-args cannot be used as given, or None."""
    for dest, own in _OWN_DIAMOND_OPTIONS.items():
        text = getattr(a, dest, None)
        if not text:
            continue
        opt = "--" + dest.replace("_", "-")
        if not isinstance(text, str):
            return f"{opt}: expected one string of DIAMOND options"
        try:
            words = shlex.split(text)
        except ValueError as e:
            return f"{opt}: {e} (in {text!r})"
        for w in words:
            if not w.startswith("-") or len(w) < 2 or w[1:2].isdigit() or w[1:2] == ".":
                continue   # a value (also a negative number)
            flag = _own_option(w)
            if dest == "diamond_args" and flag[2:] in SENSITIVITIES:
                return (f"{opt}: {w} is DIAMOND's sensitivity, which COAL sets itself: use --sensitivity "
                        f"{flag[2:]} (the report and --resume then know it)")
            if flag not in own:
                continue
            cond = _OWN_ONLY_IF.get((dest, flag), _OWN_ONLY_IF.get(flag) if dest == "deepclust_args" else None)
            if cond == "aligner":
                if getattr(a, "aligner", None) != "diamond":
                    continue
            elif cond is not None and getattr(a, cond, None) is None:
                continue
            why = own[flag]
            return (f"{opt}: {w} is an option COAL sets itself for this DIAMOND run"
                    + (f" ({why})" if why else "") + ": remove it from the string. COAL sets "
                    + ", ".join(sorted(own)) + f"; everything else in {opt} is passed on")
    return None


def _finish_args(a: argparse.Namespace) -> None:
    """Defaults that depend on other options, and the method list in the order of the rows."""
    a.method = [m for m in METHODS if m in set(a.method or ["YN"])]
    pw = a.preset == "pwalign"
    a.matrix = (a.matrix or "BLOSUM62").upper()
    a.gap_open = a.gap_open if a.gap_open is not None else (14.0 if pw else 10.0)
    a.gap_extend = a.gap_extend if a.gap_extend is not None else (4.0 if pw else 0.5)
    a.end_gaps = a.end_gaps or ("penalized" if pw else "free")
    if not hasattr(a, "flanks"):
        a.flanks = "align"


def _on_signals(handler):
    """Install handler for the termination signals of this platform (main thread only: elsewhere, and where a
    signal cannot be set, nothing changes). Returns {signal: previous handler} for _restore_signals."""
    old = {}
    for name in ("SIGTERM", "SIGHUP", "SIGBREAK"):
        sig = getattr(signal, name, None)
        if sig is None:
            continue
        try:
            old[sig] = signal.signal(sig, handler)
        except (ValueError, OSError):   # not the main thread, or not settable here
            pass
    return old


def _restore_signals(old) -> None:
    for sig, h in old.items():
        try:
            signal.signal(sig, h)
        except (ValueError, OSError):
            pass


def _execute(mode: str, a: argparse.Namespace, argv) -> PipelineResult:
    _finish_args(a)
    run = _Run(mode, a, argv)

    def terminated(signum, frame):
        raise _Terminated(signum)
    old_handlers = _on_signals(terminated)
    try:
        try:
            run.start()
            (_pairs_mode if mode == "pairs" else _clusters_mode)(run)
        except (PipelineError, OSError) as e:
            if isinstance(e, OSError):
                # a full disk, a folder that went away, a quota: one line with the file and the system's reason
                e = PipelineError(_os_error_text(e, run.out))
            e = PipelineError(run.shown(str(e)))   # the user's file names, not those of decompressed copies
            print(f"coal pipeline: error: {e}", file=sys.stderr)
            if run.lock is not None and os.path.isdir(run.out):
                report = run.finish(1, str(e))
            else:   # the folder could not be made, or it belongs to another run: nothing is written into it
                report = {"schema": SCHEMA, "exit_code": 1, "status": "failed", "error": str(e)}
            return PipelineResult(1, report, run.outputs)
        except (KeyboardInterrupt, _Terminated) as e:
            # Ctrl-C or a termination signal that reached this process (also when only the Python process got it,
            # not its process group): end the programs the run started, say in the report that the run is not
            # complete, release the folder, then let the interrupt go on (the command exits 130 or 128 + signal).
            _restore_signals(old_handlers)
            old_handlers = {}
            name = "SIGINT" if isinstance(e, KeyboardInterrupt) else e.name
            code = 130 if isinstance(e, KeyboardInterrupt) else int(e.code)
            run.kill_children()
            if run.lock is not None and os.path.isdir(run.out):
                run.finish(code, f"interrupted by {name}: the run did not finish and the tables in this folder are "
                                 "incomplete; run it again (--resume reuses the finished stages)", "interrupted")
            raise
        pairs = run.counts.get("pairs", {})
        code = 2 if run.n_failures else 0
        try:
            report = run.finish(code)
        except OSError as e:
            text = _os_error_text(e, run.out)
            print(f"coal pipeline: error: {text}", file=sys.stderr)
            return PipelineResult(1, {"schema": SCHEMA, "exit_code": 1, "status": "failed", "error": text},
                                  run.outputs)
        if report["exit_code"] != code:   # failed.tsv or the report could not be written
            return PipelineResult(report["exit_code"], report, run.outputs)
        genes = run.counts.get("genes", {})
        n_drop = genes.get("dropped", 0)
        kept = genes.get("total", 0) - genes.get("failed", 0) - n_drop
        run.say(f"summary: {genes.get('total', 0):,} records ({genes.get('failed', 0):,} not usable"
                + (f", {n_drop:,} dropped by --longest-isoform: {kept:,} genes" if n_drop else "") + "); "
                f"{pairs.get('total', 0):,} pairs: {pairs.get('ok', 0):,} ok, {pairs.get('warning', 0):,} warning, "
                f"{pairs.get('not_converged', 0):,} not_converged, {pairs.get('saturated', 0):,} saturated, "
                f"{pairs.get('no_data', 0):,} no_data, "
                f"{pairs.get('error', 0):,} error")
        if mode == "clusters":
            c = run.counts["clusters"]
            run.say(f"summary: {c['total']:,} clusters: {c['all_pairs']:,} with every pair, {c['singletons']:,} without "
                    f"a pair, {c['above_max_size']:,} above --max-cluster-size ({c['skipped']:,} skipped)")
        run.say("summary: " + ", ".join(f"{s['name']} {s['seconds']:.1f} s" + (" (reused)" if s.get("cached") else "")
                                        for s in run.stages))
        fl = run.counts.get("flagged")
        if fl is not None:
            parts = ([f"{fl['min_codons']:,} pair(s) with fewer than {a.min_codons} codons compared"]
                     if "min_codons" in fl else []) + \
                    ([f"{fl['max_ds']:,} with dS above {float(a.max_ds):.6g}"] if "max_ds" in fl else []) + \
                    ([f"{fl['min_identity']:,} with protein identity below {float(a.min_identity):.6g}%"]
                     if "min_identity" in fl else [])
            run.say("summary: flagged (status warning, numbers unchanged): " + ", ".join(parts))
        if run.n_failures:
            run.say(f"summary: {run.n_failures:,} failure(s) listed in {run.outputs['failed']}")
        run.say(f"summary: results in {run.outputs.get('dnds', run.out)}")
        run.say(f"summary: what to cite for this run: coal cite {run.outputs.get('report', '')}")
        return PipelineResult(code, report, run.outputs)
    finally:
        _restore_signals(old_handlers)
        run.release()


def _os_error_text(e: OSError, out_dir: str) -> str:
    """An operating-system error anywhere in a run as one line: what could not be done, to which file, and why."""
    import errno
    name = e.filename if isinstance(e.filename, (str, bytes)) else None
    if isinstance(name, bytes):
        name = name.decode("utf-8", "backslashreplace")
    why = e.strerror or str(e)
    what = f"cannot write {name}" if name else f"cannot write into the output folder {out_dir}"
    if e.errno in (errno.ENOENT, errno.EACCES, errno.EPERM, errno.EISDIR, errno.ENOTDIR) and name:
        what = f"cannot read or write {name}"
    text = f"{what}: {why}"
    if e.errno in (errno.ENOSPC, getattr(errno, "EDQUOT", -1)):
        text += (". The disk (or quota) of the output folder is full: free space or use another -o folder; "
                 "--resume reuses the finished stages")
    return text


def _namespace(mode: str, required: list[str], options: dict) -> argparse.Namespace:
    a = build_parser().parse_args([mode, *required])
    for k, v in options.items():
        if not hasattr(a, k):
            raise TypeError(f"unknown option {k!r} (the options are those of `coal pipeline {mode} --help`, with "
                            "underscores)")
        setattr(a, k, os.fspath(v) if isinstance(v, os.PathLike) else v)
    if isinstance(a.method, str):
        a.method = [a.method]
    if a.method:
        a.method = [_method_name(m) for m in a.method]
        bad = [m for m in a.method if m not in METHODS]
        if bad:
            raise ValueError(f"method {bad[0]!r}: choose from {', '.join(METHODS)}")
    why = _method_code_error(a) or _flag_error(a) or _args_error(a)
    if why:
        raise ValueError(why)
    for opt, allowed in (("hash_inputs", ("auto", "always", "never")), ("comeron_flavour", ("coal", "gestimator")), ("li_flavour", ("coal", "seqinr")),
                         ("terminal_stops", ("drop", "compare"))):
        if getattr(a, opt) not in allowed:
            raise ValueError(f"{opt} {getattr(a, opt)!r}: choose from {', '.join(allowed)}")
    return a


def run_pairs(query: _Path, subject: _Path, out_dir: _Path, **options: Any) -> PipelineResult:
    """Pair mode from Python: the options of `coal pipeline pairs` as keyword arguments (method=["YN", "NG"],
    ortho="rbh", aligner="nw", pairs=None, jobs=0, codontable=1, kaks_format=False, keep_alignments=False,
    keep_protein_alignments=False, resume=False, quiet=False, ...).

    Never raises for a failed run: an unreadable input, a missing tool or a full disk give exit_code 1 and
    report["error"] (the message is also printed on stderr). Raises TypeError for an option that does not exist
    and ValueError for a value that is not allowed (an unknown method, Comeron with another genetic code), before
    anything is written. KeyboardInterrupt and a termination signal pass through, after the run has ended its
    programs and written an "interrupted" report."""
    a = _namespace("pairs", ["--query", os.fspath(query), "--subject", os.fspath(subject), "-o", os.fspath(out_dir)],
                   options)
    return _execute("pairs", a, None)


def run_clusters(cds: _Path | Iterable[_Path], out_dir: _Path, clusters: _Path | None = None,
                 **options: Any) -> PipelineResult:
    """Cluster mode from Python: the options of `coal pipeline clusters` as keyword arguments (clusters=None runs
    DIAMOND DeepClust; aligner="clustalo", max_cluster_size=200, large_clusters="representative", pair_sample=1000,
    seed=1, method=["YN"], jobs=0, ...). cds: one path or a list. Returns and raises as run_pairs()."""
    files = [os.fspath(cds)] if isinstance(cds, (str, os.PathLike)) else [os.fspath(c) for c in cds]
    if clusters is not None:
        options = {**options, "clusters": os.fspath(clusters)}
    a = _namespace("clusters", ["--cds", *files, "-o", os.fspath(out_dir)], options)
    return _execute("clusters", a, None)


def main(argv: Sequence[str] | None = None) -> int:
    """`coal pipeline ...`. Exit codes: 0 done and nothing failed; 2 done, but some genes, pairs or clusters failed
    (failed.tsv, report.json), or a command-line error; 1 the run could not be done."""
    argv = sys.argv[1:] if argv is None else list(argv)
    ap = build_parser()
    a = ap.parse_args(argv)
    if a.mode is None:
        ap.print_help(sys.stderr)
        return 2
    if a.mode == "install":
        if a.sha256 and len(a.tools) != 1:
            ap.error("--sha256 names one file: install one tool at a time with it")
        try:
            _msa._check_hash_argument(a.sha256)
        except ValueError as e:
            ap.error(f"--{e}")
        try:
            for t in a.tools:
                if t == "diamond":
                    path = install_diamond(force=a.force, quiet=a.quiet, sha256=a.sha256,
                                           allow_unverified=a.allow_unverified)
                elif t == "clustalo":
                    path = _msa.install_clustalo(force=a.force, quiet=a.quiet, sha256=a.sha256,
                                                 allow_unverified=a.allow_unverified)
                else:
                    path = _msa.install_mafft(force=a.force, quiet=a.quiet)
                print(path)
        except (PipelineError, _msa.MsaError) as e:
            print(f"coal pipeline: error: {e}", file=sys.stderr)
            return 1
        return 0
    if a.digits < 1 or a.digits > 17:
        ap.error("--digits: expected 1 to 17")
    if a.jobs < 0:
        ap.error("-j: expected 0 (one per CPU) or more")
    why = _flag_error(a) or _args_error(a)
    if why:
        ap.error(why)
    if _method_code_error(a):
        ap.error(_method_code_error(a))
    try:
        return _execute(a.mode, a, ["pipeline", *argv]).exit_code
    except KeyboardInterrupt:
        print("coal pipeline: interrupted (the report says so; nothing is left running)", file=sys.stderr)
        return 130
    except _Terminated as e:
        print(f"coal pipeline: interrupted by {e.name} (the report says so; nothing is left running)",
              file=sys.stderr)
        return int(e.code)
