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
Not part of the stable `import coal` API yet (python-api.md): it may change before 1.0.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import random
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

from . import __version__
from . import msa as _msa
from .core import CoalError
from .pipeline import release_output_lock, take_output_lock
from .tables import SUPPORTED

__all__ = ["run_pairs", "run_clusters", "main", "PipelineResult", "PipelineError", "find_diamond",
           "install_diamond", "codon_dnds", "summarise", "PAIR_COLUMNS", "SUMMARY_COLUMNS", "SCHEMA"]

SCHEMA = "coal.pipeline/1"
DIAMOND_VERSION = "2.2.8"   # the version install_diamond() fetches; any DIAMOND >= 2.1 on PATH is used as it is
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
REPORT_MAX_FAILURES = 1000   # the JSON report lists this many failures; failed.tsv has all of them
_CLUSTER_HEADERS = {("centroid", "member"), ("representative", "member"), ("representative_id", "member_id"),
                    ("cluster", "member"), ("cluster_id", "member_id"), ("rep", "member")}
_EXE = ".exe" if os.name == "nt" else ""


# Files and folders a run writes into the output folder (besides work/): whatever an earlier run left of them is
# removed when a run starts, so that the folder never mixes two runs (a dnds.kaks or codon_alignments.fasta of a
# run with other options, the dnds.tsv.tmp of a killed run).
_RUN_OUTPUTS = ("dnds.tsv", "dnds.tsv.tmp", "dnds.kaks", "codon_alignments.fasta", "cluster_summary.tsv",
                "cluster_summary.tsv.tmp", "genes.tsv", "orthologs.tsv", "unpaired.tsv", "clusters.tsv",
                "realign.tsv", "realign.tsv.body", "failed.tsv", "report.json")


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
    exit_code: int                 # 0 nothing failed; 2 finished, some genes / pairs / clusters failed
    report: dict                   # what out_dir/report.json holds
    outputs: dict = field(default_factory=dict)   # name -> path (dnds, cluster_summary, genes, failed, report, ...)


# ------------------------------------------------------------------------------------------------ tools

def find_diamond() -> str | None:
    """DIAMOND: $COAL_DIAMOND, else `diamond` on PATH, else the one install_diamond() put in COAL's cache
    (coal.msa.cache_dir())."""
    return _msa._find("diamond", "COAL_DIAMOND")


def _diamond_version(exe: str) -> str:
    """'2.2.8' ('' if it does not run)."""
    try:
        r = subprocess.run([exe, "version"], capture_output=True, text=True, timeout=60, stdin=subprocess.DEVNULL,
                           encoding="utf-8", errors="backslashreplace")
    except (OSError, subprocess.SubprocessError):
        return ""
    words = (r.stdout or r.stderr).split()
    return words[-1] if r.returncode == 0 and "diamond" in (r.stdout or r.stderr).lower() and words else ""


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
        raise PipelineError(f"coal pipeline needs coal-cpp, which was not found (looked at {dispatch.LOOKUP_ORDER}); "
                            "install a platform wheel of coal-align or build it (docs/cpp-engine.md)")
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

class _Run:
    """State shared by the stages of one run: folders, tools, the report, the failures file."""

    def __init__(self, mode: str, a: argparse.Namespace, argv):
        self.mode, self.a = mode, a
        self.out = os.path.abspath(a.out_dir)
        self.work = os.path.join(self.out, "work")
        self.t0 = time.time()
        self.quiet = a.quiet
        self.jobs = a.jobs if a.jobs and a.jobs > 0 else (os.cpu_count() or 1)
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
        self._earlier: tuple[set, bool] = (set(), False)
        self.pairs_in_place = False   # --pairs is DIR/orthologs.tsv itself: read, not written again
        self._engine = None     # the running coal-cpp codon-dnds (killed if the run is interrupted)

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
        self._failed_fh = open(self.outputs["failed"], "w", encoding="utf-8", errors="backslashreplace", newline="\n")
        self._failed_fh.write("\t".join(FAILED_COLUMNS) + "\n")

    # -- the folder: this run's inputs, what it will write, what an earlier run wrote
    def _input_paths(self) -> list[tuple[str, str]]:
        """(option, path) of every file this run reads."""
        a = self.a
        got = []
        for opt in ("query", "subject", "pairs", "clusters", "diamond"):
            v = getattr(a, opt, None)
            if v:
                got.append(("--" + opt, os.fspath(v)))
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
        names = ["failed.tsv", "genes.tsv", "dnds.tsv", "dnds.tsv.tmp"]
        if not a.report:
            names.append("report.json")
        if a.kaks_format:
            names.append("dnds.kaks")
        if self.mode == "pairs":
            names += ["orthologs.tsv"]
            if not a.pairs:
                names.append("unpaired.tsv")
            if a.keep_alignments:
                names.append("codon_alignments.fasta")
        else:
            names += ["cluster_summary.tsv", "cluster_summary.tsv.tmp"]
            if not a.clusters:
                names.append("clusters.tsv")
                if a.realign:
                    names += ["realign.tsv", "realign.tsv.body"]
        return names

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
        inside = (os.path.realpath(self.work), os.path.realpath(self.path("codon_alignments")))
        for opt, inp in self._inputs:
            real = os.path.realpath(inp)
            for d in inside:
                if real.startswith(d + os.sep):
                    raise PipelineError(f"{opt} {inp} is inside {d}, a folder this run writes: nothing was "
                                        "changed. Use another output folder, or move the input out of it")

    def _earlier_outputs(self) -> tuple[set, bool]:
        """What an earlier run is known to have written here: the names in work/outputs.json (written when a
        run starts, so a killed run counts), or, for a run of a version without that file, the outputs its
        report.json lists. Second value: the .fasta files of codon_alignments/ are that run's (a cluster run with
        --keep-alignments)."""
        names: set = set()
        codon_dir = False
        try:
            with open(self.wpath("outputs.json"), encoding="utf-8", errors="backslashreplace") as fh:
                rec = json.load(fh)
            names = {n for n in rec["files"] if isinstance(n, str)}
            codon_dir = rec.get("codon_alignments_dir") is True
        except (OSError, ValueError, KeyError, TypeError):
            try:
                with open(self.path("report.json"), encoding="utf-8", errors="backslashreplace") as fh:
                    rep = json.load(fh)
                if isinstance(rep, dict) and rep.get("schema") == SCHEMA and isinstance(rep.get("outputs"), dict):
                    names = {v for v in rep["outputs"].values() if isinstance(v, str)}
                    names |= {"report.json", "failed.tsv"}
                    codon_dir = rep.get("mode") == "clusters" and "codon_alignments" in names
                    names |= {n + ".tmp" for n in ("dnds.tsv", "cluster_summary.tsv") if n in names}
            except (OSError, ValueError):
                pass
        d = self.path("codon_alignments")
        codon_dir = codon_dir and os.path.isdir(d) and not os.path.islink(d)
        if self.a.keep_alignments and self.mode == "clusters" and os.path.lexists(d) and not codon_dir:
            # not COAL's: the engine would write its files among (or over) what the folder holds
            if not os.path.isdir(d) or os.path.islink(d) or os.listdir(d):
                raise PipelineError(f"{d} exists and was not written by a COAL run: --keep-alignments writes one "
                                    "file per cluster there and would mix with (or replace) what it holds. "
                                    "Nothing was changed: move it away or use another output folder")
        return {n for n in names if n in _RUN_OUTPUTS}, codon_dir

    def _clear_earlier_outputs(self) -> None:
        """Remove the files an earlier COAL run wrote into this folder, so the folder holds this run's results
        only. Only files a run is known to have written are removed (_earlier_outputs), never one of this run's
        inputs, and never a file COAL has no record of: a `clusters.tsv` or a codon_alignments/ folder that was
        the user's stays. --resume keeps the stage outputs it can reuse (clusters.tsv, realign.tsv; everything
        under work/ is kept anyway)."""
        a = self.a
        keep = set()
        if a.resume and self.mode == "clusters" and not a.clusters:
            keep.add("clusters.tsv")
            if a.realign:
                keep.add("realign.tsv")
        known, codon_dir = self._earlier
        d = self.path("codon_alignments")
        try:
            for name in sorted(known):
                p = self.path(name)
                if name not in keep and (os.path.isfile(p) or os.path.islink(p)) and self._is_input(p) is None:
                    os.remove(p)
            if codon_dir and os.path.isdir(d) and not os.path.islink(d):
                for n in os.listdir(d):
                    p = os.path.join(d, n)
                    if n.endswith(".fasta") and self._is_input(p) is None:
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
                       "codon_alignments_dir": bool(self.a.keep_alignments and self.mode == "clusters")}, fh)
            fh.write("\n")
        os.replace(tmp, self.wpath("outputs.json"))

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

    def say(self, text: str) -> None:
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
            detail = fn() or {}
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
        with open(logp, "ab") as fh:
            fh.write(("$ " + " ".join(shlex.quote(c) for c in cmd) + "\n").encode())
            fh.flush()
            try:
                r = subprocess.run(cmd, stdin=subprocess.DEVNULL, stdout=fh, stderr=subprocess.STDOUT)
            except OSError as e:
                raise PipelineError(f"{what}: could not start {cmd[0]} ({e.strerror or e})") from None
        if r.returncode != 0:
            raise PipelineError(f"{what} failed (exit {r.returncode}); the end of {logp}:\n{_tail(logp)}")

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

        def go():
            self.run_tool([self.cpp, "translate", *paths, "-o", faa, "--genes", genes, "-c",
                           str(self.a.codontable)], "translate", f"translating {paths[0]}")
        self.stage(f"translate_{name}", _sig("translate", _file_sig(paths), self.a.codontable, __version__),
                   [faa, genes], go)
        return faa, genes, {}

    def read_genes(self, genes_path: str, set_name: str, out_fh, order: dict | None = None) -> dict:
        """Copy a genes table into out_dir/genes.tsv (with the `set` column), report the failed genes, and return
        {id: reason} of the failed ones. order: filled with {id: rank} of the usable genes (file order)."""
        failed: dict[str, str] = {}
        c = self.counts.setdefault("genes", {"total": 0, "ok": 0, "warning": 0, "failed": 0})
        codes = self.gene_failures.setdefault(set_name, {"usable": 0, "codes": {}, "first": None})
        with open(genes_path, encoding="utf-8", errors="backslashreplace") as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                out_fh.write(set_name + "\t" + line)
                c["total"] += 1
                c[f[7]] = c.get(f[7], 0) + 1
                if f[7] == "failed":
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
        cmd = [self.cpp, "codon-dnds", "-i", "-", "-o", out + ".tmp", "--failures", fail, "-c", str(a.codontable),
               "-j", str(self.jobs), "--digits", str(a.digits), "--matrix", a.matrix, "--gap-open", str(a.gap_open),
               "--gap-extend", str(a.gap_extend), "--end-gaps", a.end_gaps, "--max-cells", str(a.max_cells),
               "--flanks", a.flanks, "--yn-flavour", a.yn_flavour, "--yn-weighting", str(a.yn_weighting),
               "--comeron-flavour", a.comeron_flavour, "--li-flavour", a.li_flavour,
               "--terminal-stops", a.terminal_stops]
        for m in a.method:
            cmd += ["-m", m]
        if group_column:
            cmd += ["--group-column", group_column]
        if a.kaks_format:
            self.outputs["kaks"] = self.path("dnds.kaks")
            cmd += ["--kaks-out", self.outputs["kaks"]]
        if a.keep_alignments:
            if self.mode == "pairs":
                self.outputs["codon_alignments"] = self.path("codon_alignments.fasta")
                cmd += ["--codon-out", self.outputs["codon_alignments"]]
            else:
                d = self.outputs["codon_alignments"] = self.path("codon_alignments")
                os.makedirs(d, exist_ok=True)   # (an earlier run's files were removed when this run started)
                cmd += ["--codon-dir", d]
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
                for tmp in (out + ".tmp",):
                    try:
                        os.remove(tmp)
                    except OSError:
                        pass
                raise
            finally:
                self._engine = None
        if rc != 0 or broken:
            raise PipelineError(f"coal-cpp codon-dnds failed (exit {rc}); the end of {logp}:\n{_tail(logp)}")
        os.replace(out + ".tmp", out)
        self.outputs["dnds"] = out
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
        return detail

    # -- the end
    def finish(self, exit_code: int, error: str | None = None, status: str | None = None) -> dict:
        """Close failed.tsv and write the report. status: "completed" (exit 0), "completed_with_failures" (2),
        "failed" (1: the run could not be done) or "interrupted" (a signal or Ctrl-C ended it: the tables in the
        folder are incomplete or missing)."""
        if self._failed_fh:
            self._failed_fh.close()
            self._failed_fh = None
        a = self.a
        skip = {"out_dir", "quiet", "report", "resume", "install_tools"}
        if status is None:
            status = {0: "completed", 2: "completed_with_failures"}.get(exit_code, "failed")
        report = {
            "schema": SCHEMA, "coal_version": __version__, "mode": self.mode, "command": self.argv,
            "out_dir": self.out, "exit_code": exit_code, "status": status,
            "parameters": {k: v for k, v in sorted(vars(a).items()) if k not in skip},
            "tools": self.tools, "stages": self.stages,
            "seconds": round(time.time() - self.t0, 3), "counts": self.counts,
            "failures_total": self.n_failures, "failures": self.failures,
            "failures_truncated": self.n_failures > len(self.failures),
            "outputs": {k: os.path.relpath(v, self.out) for k, v in sorted(self.outputs.items())},
        }
        if error:
            report["error"] = error
        rp = a.report or self.path("report.json")
        try:
            os.makedirs(os.path.dirname(os.path.abspath(rp)), exist_ok=True)
            with open(rp, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as fh:
                json.dump(report, fh, indent=1, default=str)
                fh.write("\n")
            self.outputs["report"] = rp
        except OSError as e:
            print(f"coal pipeline: error: could not write the report {rp}: {e.strerror or e}", file=sys.stderr)
        return report


def _count_pairs(run: _Run, dnds_path: str, group_col: bool) -> None:
    """Pair counts by status from dnds.tsv (the first method's rows), and the error pairs into failed.tsv."""
    c = {"total": 0, "ok": 0, "warning": 0, "not_converged": 0, "saturated": 0, "no_data": 0, "error": 0}
    k = 1 if group_col else 0
    first = None
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
            if f[-2] == "error":
                run.fail("pair", "pair_failed", f[-1], cluster=f[0] if group_col else "", seq1=f[k], seq2=f[k + 1])
    run.counts["pairs"] = c


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
    with open(run.outputs["genes"], "w", encoding="utf-8", errors="backslashreplace", newline="\n") as gf:
        with open(q_genes, encoding="utf-8", errors="backslashreplace") as fh:
            gf.write("set\t" + fh.readline())
        q_failed = run.read_genes(q_genes, "query", gf, q_order)
        s_failed = run.read_genes(s_genes, "subject", gf, s_order)
    run.say(f"translated {len(q_order) + len(q_failed):,} query and {len(s_order) + len(s_failed):,} subject CDS "
            f"({len(q_failed) + len(s_failed):,} not usable)")
    run.require_usable_genes("query", [a.query])
    run.require_usable_genes("subject", [a.subject])

    orth = run.outputs["orthologs"] = run.path("orthologs.tsv")
    hits: dict[str, list[str]] = {}
    if a.pairs:
        if a.aligner == "diamond":
            raise PipelineError("--aligner diamond uses the alignments of the DIAMOND search, so it cannot be "
                                "combined with --pairs")
        pairs = _read_pairs_file(a.pairs)
        run.stages.append({"name": "orthologs", "seconds": 0.0, "cached": False, "source": "given",
                           "pairs": len(pairs)})
        if not run.pairs_in_place:   # (--pairs DIR/orthologs.tsv: the input is the table; it is not rewritten)
            with open(orth, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as fh:
                fh.write("query\tsubject\n")
                for q, s in pairs:
                    fh.write(f"{q}\t{s}\n")
    else:
        pairs, hits = _search(run, q_faa, s_faa, q_order, s_order, orth)
    run.say(f"{len(pairs):,} pairs ({'given' if a.pairs else a.ortho})")
    run.counts["orthologs"] = len(pairs)

    q_idx, s_idx = FastaIndex(a.query), FastaIndex(a.subject)
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
            res = _msa.align_clusters(todo, tool=a.aligner, jobs=run.jobs, threads=1, pairs_with="msa", exe=exe,
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


def _read_pairs_file(path: str) -> list[tuple[str, str]]:
    pairs = []
    try:
        fh = open(path, encoding="utf-8", errors="backslashreplace")
    except OSError as e:
        raise PipelineError(f"cannot open {path!r}: {e.strerror or e}") from None
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
    common = ["--max-target-seqs", "1", "--evalue", str(a.evalue), "--threads", str(run.jobs), "--quiet"]
    if a.sensitivity != "default":
        common.append("--" + a.sensitivity)
    if gapped:
        common += ["--masking", "0"]   # DIAMOND's alignment inside a masked stretch is arbitrary (docs/alignment.md)
    fwd, rev = run.wpath("query_vs_subject.tsv"), run.wpath("subject_vs_query.tsv")
    both = a.ortho == "rbh"

    def go():
        run.run_tool([dmd, "makedb", "--in", s_faa, "-d", run.wpath("subject"), "--threads", str(run.jobs),
                      "--quiet"], "diamond", "diamond makedb")
        cols = _HIT_COLS + (["qseq_gapped", "sseq_gapped"] if gapped else [])
        run.run_tool([dmd, "blastp", "-q", q_faa, "-d", run.wpath("subject"), "-o", fwd, "--outfmt", "6", *cols,
                      *common, *extra], "diamond", "diamond blastp (query against subject)")
        if both:
            run.run_tool([dmd, "makedb", "--in", q_faa, "-d", run.wpath("query"), "--threads", str(run.jobs),
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
    with open(orth, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as fh:
        fh.write("query\tsubject\tpident\tlength\tevalue\tbitscore\n")
        for q, s in pairs:
            fh.write("\t".join([q, s, *hits[q][1:5]]) + "\n")
    # genes without a partner
    unp = run.outputs["unpaired"] = run.path("unpaired.tsv")
    paired_q = {q for q, _ in pairs}
    paired_s = {s for _, s in pairs}
    n = {"no_hit": 0, "not_reciprocal": 0}
    with open(unp, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as fh:
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
    faa, genes, _ = run.translate("proteins", list(a.cds))
    run.outputs["genes"] = run.path("genes.tsv")
    order: dict[str, int] = {}
    with open(run.outputs["genes"], "w", encoding="utf-8", errors="backslashreplace", newline="\n") as gf:
        with open(genes, encoding="utf-8", errors="backslashreplace") as fh:
            gf.write("set\t" + fh.readline())
        failed = run.read_genes(genes, "cds", gf, order)
    run.say(f"translated {len(order) + len(failed):,} CDS ({len(failed):,} not usable)")
    run.require_usable_genes("cds", list(a.cds))

    if a.clusters:
        cfile = a.clusters
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

        def go():
            run.run_tool([dmd, "deepclust", "--db", faa, "--out", cfile, "--threads", str(run.jobs), "--quiet",
                          *opts, *extra], "diamond", "diamond deepclust")
        run.stage("deepclust", _sig("deepclust", _file_sig([faa]), run.tools["diamond"]["version"], opts, extra),
                  [cfile], go)
        if a.realign:
            rfile = run.outputs["realign"] = run.path("realign.tsv")
            rextra = shlex.split(a.realign_args) if a.realign_args else []

            def go2():
                run.run_tool([dmd, "realign", "--db", faa, "--clusters", cfile, "--out", rfile + ".body", "--threads",
                              str(run.jobs), "--quiet", "--outfmt", "6", "qseqid", "sseqid", "approx_pident",
                              "evalue", "bitscore", "qcovhsp", "scovhsp", *rextra], "diamond", "diamond realign")
                with open(rfile, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as out, \
                        open(rfile + ".body", encoding="utf-8", errors="backslashreplace") as body:
                    if not rextra:
                        out.write("member\trepresentative\tapprox_pident\tevalue\tbitscore\tqcovhsp\tscovhsp\n")
                    shutil.copyfileobj(body, out)
                os.remove(rfile + ".body")
            run.stage("realign", _sig("realign", _file_sig([cfile]), run.tools["diamond"]["version"], rextra),
                      [rfile], go2)
    t0 = time.perf_counter()
    clusters, dup = read_clusters(cfile)
    if dup:
        run.say(f"warning: {dup:,} member line(s) of {cfile} name a sequence that is already in a cluster; the "
                "first one is used")
    run.say(f"{len(clusters):,} clusters, {sum(len(m) for m in clusters.values()):,} sequences")

    cds_idx, prot_idx = FastaIndex(list(a.cds)), FastaIndex(faa)
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
                res = _msa.align_clusters(todo, tool=tool, jobs=run.jobs, threads=a.msa_threads, pairs_with="msa",
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
    with open(dnds_path, encoding="utf-8", errors="backslashreplace") as fh, \
            open(out_path, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as out:
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

def codon_dnds(groups, out_tsv: str, method=("YN",), codontable: int = 1, jobs: int = 0,
               group_column: str | None = None, codon_out: str | None = None, kaks_out: str | None = None,
               failures: str | None = None, extra_args=(), cpp: str | None = None) -> int:
    """Protein alignment -> codon alignment -> Ka/Ks for groups held in Python, through one `coal-cpp codon-dnds`
    process. groups: an iterable of (group_id, members[, pairs]) with members = [(id, cds)] or [(id, cds,
    protein_row)]: two members without rows are aligned here; three or more need rows. pairs: "all" (default),
    "rep" or a list of (i, j). Writes out_tsv (PAIR_COLUMNS, with group_column first when given); returns the
    number of groups sent. Raises PipelineError when coal-cpp is missing or fails."""
    cpp = cpp or _find_cpp()
    cmd = [cpp, "codon-dnds", "-i", "-", "-o", out_tsv, "-c", str(codontable), "-j", str(jobs), "-q"]
    for m in method:
        cmd += ["-m", m]
    for flag, v in (("--group-column", group_column), ("--codon-out", codon_out), ("--kaks-out", kaks_out),
                    ("--failures", failures)):
        if v:
            cmd += [flag, v]
    cmd += list(extra_args)
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
        raise PipelineError(f"coal-cpp codon-dnds failed (exit {rc}): "
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
    g.add_argument("--report", metavar="FILE.json", help="the JSON run report (default: DIR/report.json)")
    g.add_argument("--digits", type=int, default=8, metavar="N", help="significant digits (default 8)")
    g.add_argument("-q", "--quiet", action="store_true", help="no progress or summary on stderr")
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
                   help="threads, and aligner processes at a time (default 0 = one per CPU)")
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
    g.add_argument("--diamond-args", metavar="'ARGS'", help="further arguments for `diamond blastp`, in one string")

    p = sub.add_parser("clusters", parents=[common], help="gene families: every pair within each cluster",
                       description="Clusters of coding sequences: translate, DIAMOND DeepClust (or a given clusters "
                                   "file), a protein MSA per cluster, codon MSA, Ka/Ks for every pair of a cluster, "
                                   "and a per-cluster summary.")
    g = p.add_argument_group("Clusters")
    g.add_argument("--cds", required=True, nargs="+", metavar="CDS.fasta", help="coding sequences (one or more files)")
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
    g.add_argument("--deepclust-args", metavar="'ARGS'", help="further arguments for `diamond deepclust`")
    g.add_argument("--realign", action="store_true",
                   help="also run `diamond realign` (members against their representatives) into DIR/realign.tsv")
    g.add_argument("--realign-args", metavar="'ARGS'", help="further arguments for `diamond realign`")

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
                e = PipelineError(f"{e.strerror or e}" + (f": {e.filename}" if e.filename else ""))
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
                for tmp in ("dnds.tsv.tmp", "cluster_summary.tsv.tmp", "realign.tsv.body"):
                    try:
                        os.remove(run.path(tmp))
                    except OSError:
                        pass
                run.finish(code, f"interrupted by {name}: the run did not finish and the tables in this folder are "
                                 "incomplete; run it again (--resume reuses the finished stages)", "interrupted")
            raise
        pairs = run.counts.get("pairs", {})
        code = 2 if run.n_failures else 0
        report = run.finish(code)
        genes = run.counts.get("genes", {})
        run.say(f"summary: {genes.get('total', 0):,} genes ({genes.get('failed', 0):,} not usable); "
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
        if run.n_failures:
            run.say(f"summary: {run.n_failures:,} failure(s) listed in {run.outputs['failed']}")
        run.say(f"summary: results in {run.outputs.get('dnds', run.out)}")
        return PipelineResult(code, report, run.outputs)
    finally:
        _restore_signals(old_handlers)
        run.release()


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
    why = _method_code_error(a)
    if why:
        raise ValueError(why)
    for opt, allowed in (("comeron_flavour", ("coal", "gestimator")), ("li_flavour", ("coal", "seqinr")),
                         ("terminal_stops", ("drop", "compare"))):
        if getattr(a, opt) not in allowed:
            raise ValueError(f"{opt} {getattr(a, opt)!r}: choose from {', '.join(allowed)}")
    return a


def run_pairs(query, subject, out_dir, **options) -> PipelineResult:
    """Pair mode from Python: the options of `coal pipeline pairs` as keyword arguments (method=["YN", "NG"],
    ortho="rbh", aligner="nw", pairs=None, jobs=0, codontable=1, kaks_format=False, keep_alignments=False,
    resume=False, quiet=False, ...). Never raises for a failed run: see PipelineResult.exit_code and .report."""
    a = _namespace("pairs", ["--query", os.fspath(query), "--subject", os.fspath(subject), "-o", os.fspath(out_dir)],
                   options)
    return _execute("pairs", a, None)


def run_clusters(cds, out_dir, clusters=None, **options) -> PipelineResult:
    """Cluster mode from Python: the options of `coal pipeline clusters` as keyword arguments (clusters=None runs
    DIAMOND DeepClust; aligner="clustalo", max_cluster_size=200, large_clusters="representative", pair_sample=1000,
    seed=1, method=["YN"], jobs=0, ...). cds: one path or a list."""
    files = [os.fspath(cds)] if isinstance(cds, (str, os.PathLike)) else [os.fspath(c) for c in cds]
    if clusters is not None:
        options = {**options, "clusters": os.fspath(clusters)}
    a = _namespace("clusters", ["--cds", *files, "-o", os.fspath(out_dir)], options)
    return _execute("clusters", a, None)


def main(argv=None) -> int:
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
