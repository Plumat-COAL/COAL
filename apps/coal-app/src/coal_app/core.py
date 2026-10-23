"""Gradio-free logic behind the COAL GUI and MCP tools.

Everything here talks to the engine only through ``python -m coal`` in a
subprocess. Newer engine features (``-output clustal|paml|codon``,
``--report FILE.json``, ``--keep-going``) are feature-detected from
``python -m coal --help``; when the engine lacks ``--report`` we build a report
with the same schema (``coal.report/1``) from COAL's stderr instead.
"""
from __future__ import annotations

import errno
import hashlib
import html
import json
import os
import re
import shlex
import signal
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE_REPORT = os.path.join(HERE, "sample_report.json")


def _find_examples() -> str:
    """The examples folder: COAL_EXAMPLES_DIR, then a repo checkout, then the copy the coal-align wheel installs under
    <prefix>/share/coal/examples."""
    import sysconfig
    env = os.environ.get("COAL_EXAMPLES_DIR")
    if env:
        return os.path.normpath(env)
    # repo checkout: <root>/apps/coal-app/src/coal_app -> <root>/examples
    repo = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", "examples"))
    # a folder made by space/build_space.py: <space>/src/coal_app -> <space>/examples
    bundle = os.path.normpath(os.path.join(HERE, "..", "..", "examples"))
    cands = [repo, bundle]
    for base in (sys.prefix, sysconfig.get_path("data"), getattr(sys, "base_prefix", sys.prefix)):
        if base:
            cands.append(os.path.join(base, "share", "coal", "examples"))
    for c in cands:
        if os.path.isdir(os.path.join(c, "clean")):
            return c
    return repo


EXAMPLES_DIR = _find_examples()
EXAMPLES = {
    "clean": ("Clean nuclear genes", "Eight genes that fit the standard code. Everything should be OK."),
    # runs them with code 1 on purpose (the demo of the genetic-code tip)
    "mito": ("Mitochondrial (wrong code)", "Vertebrate mitochondrial genes, run with the standard code 1 on purpose."),
    "messy_ids": ("Messy IDs", "Protein IDs like AT1G01010.1, DNA IDs like AT1G01010."),
    "one_broken_gene": ("One broken gene", "19 good genes and one with a frameshift."),
}

# The settings each example runs with; an example button puts them into the form (so the dropdown and checkboxes
# show what ran) before it runs, so an example never inherits a code left in the form by an earlier run. Every
# example uses the standard code 1: the Mitochondrial genes deliberately so (the demo shows the genetic-code tip).
# The output format stays the user's choice (it does not change what is found).
EXAMPLE_SETTINGS: dict[str, dict] = {
    name: {"codontable": 1, "nogap": False, "nomismatch": False, "strict": False, "blockonly": False}
    for name in EXAMPLES
}

# NCBI genetic codes supported by COAL/PAL2NAL, with plain names.
GENETIC_CODES: dict[int, str] = {
    1: "Standard (nuclear genes of most organisms)",
    2: "Vertebrate mitochondrial",
    3: "Yeast mitochondrial",
    4: "Mold, protozoan and coelenterate mitochondrial; Mycoplasma",
    5: "Invertebrate mitochondrial",
    6: "Ciliate, dasycladacean and Hexamita nuclear",
    9: "Echinoderm and flatworm mitochondrial",
    10: "Euplotid nuclear",
    11: "Bacterial, archaeal and plant plastid",
    12: "Alternative yeast nuclear",
    13: "Ascidian mitochondrial",
    14: "Alternative flatworm mitochondrial",
    15: "Blepharisma nuclear",
    16: "Chlorophycean mitochondrial",
    21: "Trematode mitochondrial",
    22: "Scenedesmus obliquus mitochondrial",
    23: "Thraustochytrium mitochondrial",
}

OUTPUT_FORMATS = {
    "fasta": "FASTA (one sequence per record; works everywhere)",
    "clustal": "CLUSTAL (aligned blocks, easy to read)",
    "paml": "PAML (for codeml / dN/dS analysis)",
    "codon": "Codon (codons separated by spaces)",
}

STATUS = {
    # status: (icon, word, css class)
    "ok": ("✔", "OK", "ok"),
    "warning": ("⚠", "Warning", "warning"),
    "failed": ("✖", "Failed", "failed"),
    "not_checked": ("–", "Not checked", "unchecked"),
}

ENGINE_TIMEOUT = 600  # seconds for one engine run
PROBE_TIMEOUT = 30    # seconds for one genetic-code probe
PROBE_BUDGET = 90     # seconds for all probes of one run together
PROBE_SAMPLE = 200    # at most this many problem sequences are re-run per probe

# Per-run output folders (GUI and MCP tools). Pruned by age; see prune_runs().
RUNS_DIR = os.environ.get("COAL_RUNS_DIR") or os.path.join(tempfile.gettempdir(), "coal_app_runs")
RUNS_MAX_AGE_HOURS = 24


RUN_MARKER = ".coal_run"  # every folder new_run_dir() makes holds this file


def prune_runs(max_age_hours: float = RUNS_MAX_AGE_HOURS) -> int:
    """Delete run folders in RUNS_DIR not modified for max_age_hours. Returns how many were removed.
    Only folders named run_* that hold COAL's marker file are touched, so COAL_RUNS_DIR pointing at a folder of
    the user's own deletes nothing of theirs."""
    import shutil
    try:
        names = os.listdir(RUNS_DIR)
    except OSError:
        return 0
    cutoff = time.time() - max_age_hours * 3600
    removed = 0
    for n in names:
        p = os.path.join(RUNS_DIR, n)
        try:
            if (n.startswith("run_") and os.path.isdir(p) and not os.path.islink(p)
                    and os.path.isfile(os.path.join(p, RUN_MARKER)) and os.path.getmtime(p) < cutoff):
                shutil.rmtree(p, ignore_errors=True)
                removed += not os.path.exists(p)
        except OSError:
            continue
    return removed


def ensure_runs_dir() -> str:
    """Create RUNS_DIR if needed and return it. An unusable folder (COAL_RUNS_DIR under a path that cannot hold
    one, an existing file, no permission) is a ValueError with one plain sentence, so the GUI shows a banner, the
    MCP tool a tool error and the launcher a message instead of a traceback."""
    try:
        os.makedirs(RUNS_DIR, exist_ok=True)
        if not os.access(RUNS_DIR, os.W_OK):
            raise PermissionError(errno.EACCES, os.strerror(errno.EACCES))
    except OSError as e:
        raise ValueError(f"Could not create the results folder {RUNS_DIR}: {e.strerror or e}. "
                         "Set COAL_RUNS_DIR to a folder COAL can write to.") from e
    return RUNS_DIR


def new_run_dir() -> str:
    """A fresh folder for one run's outputs inside RUNS_DIR (old runs are pruned first). ValueError when RUNS_DIR
    cannot be used (see ensure_runs_dir)."""
    ensure_runs_dir()
    prune_runs()
    try:
        d = tempfile.mkdtemp(prefix="run_", dir=RUNS_DIR)
        with open(os.path.join(d, RUN_MARKER), "w", encoding="utf-8") as fh:
            fh.write("Made by COAL; this folder is deleted automatically after "
                     f"{RUNS_MAX_AGE_HOURS} hours.\n")
    except OSError as e:
        raise ValueError(f"Could not create a run folder in {RUNS_DIR}: {e.strerror or e}. "
                         "Set COAL_RUNS_DIR to a folder COAL can write to.") from e
    return d


# ---------------------------------------------------------------- local vs hosted

def hosted() -> bool:
    """True when running as a shared web service (the Hugging Face Space sets COAL_HOSTED=1)."""
    return os.environ.get("COAL_HOSTED", "") == "1"


def sandbox_roots() -> list[str]:
    """Folders the tools may read when hosted: the bundled examples, run outputs and Gradio's uploads."""
    gradio_tmp = os.environ.get("GRADIO_TEMP_DIR") or os.path.join(tempfile.gettempdir(), "gradio")
    return [os.path.realpath(p) for p in (EXAMPLES_DIR, RUNS_DIR, gradio_tmp)]


def check_path_allowed(path: str) -> str | None:
    """None if the tools may read `path`, else a plain-English reason. Locally every path is allowed;
    when hosted only files inside sandbox_roots() (after resolving symlinks and '..')."""
    if not hosted():
        return None
    real = os.path.realpath(path)
    for root in sandbox_roots():
        if real == root or real.startswith(root.rstrip(os.sep) + os.sep):
            return None
    return (f"reading '{path}' is not allowed on this shared server: only uploaded files, the bundled "
            f"examples and this server's own results can be used. Install COAL locally to use your own paths.")


# ---------------------------------------------------------------- inputs given as paths (APP-L)

def path_inputs_enabled() -> bool:
    """Path boxes are offered only when running locally (never on the hosted Space)."""
    return not hosted()


def human_size(n: int) -> str:
    for unit in ("bytes", "KB", "MB", "GB"):
        if n < 1000 or unit == "GB":
            return f"{n:,} {unit}" if unit == "bytes" else f"{n:.1f} {unit}"
        n /= 1000
    return f"{n:.1f} GB"


_INVISIBLE = "\u202a\u202b\u202c\u202d\u202e\u200e\u200f\u200b\ufeff\u2066\u2067\u2068\u2069"


def clean_pasted_path(text: str | None) -> str:
    """A pasted path as the user meant it: surrounding whitespace, quotes and invisible
    direction marks removed (Windows "Properties > Security" copies a leading U+202A), and file:// URLs
    (as browsers and file managers copy them, with %20 etc.) turned into a path."""
    t = (text or "")
    prev = None
    while t != prev:
        prev = t
        t = t.strip().strip(_INVISIBLE).strip().strip('"').strip("'")
    if t.lower().startswith("file:"):
        from urllib.parse import urlparse
        from urllib.request import url2pathname
        u = urlparse(t)
        path = url2pathname(u.path)
        if u.netloc and u.netloc.lower() != "localhost":
            path = "\\\\" + u.netloc + url2pathname(u.path) if os.name == "nt" else "//" + u.netloc + path
        t = path
    return t


def resolve_input(upload: str | None, path_text: str | None) -> tuple[str | None, str]:
    """(path to use, short description) for one GUI input. A typed/pasted path wins over an upload and is
    read in place (never copied). Raises ValueError with a plain message if the path is unusable."""
    text = clean_pasted_path(path_text)
    if text:
        if not path_inputs_enabled():
            raise ValueError("this shared server cannot read paths; please upload the file instead")
        p = os.path.abspath(os.path.expanduser(text))
        if not os.path.exists(p):
            raise ValueError(f"'{text}' does not exist on this computer")
        if os.path.isdir(p):
            raise ValueError(f"'{text}' is a folder; give the file inside it")
        if not os.path.isfile(p):
            raise ValueError(f"'{text}' is not a regular file")
        try:
            with open(p, "rb") as fh:
                fh.read(1)
        except OSError as e:
            raise ValueError(f"'{text}' cannot be read ({e.strerror or e})") from None
        return p, f"{os.path.basename(p)}, {human_size(os.path.getsize(p))} (read in place)"
    if upload:
        try:
            return upload, f"{os.path.basename(upload)}, {human_size(os.path.getsize(upload))} (uploaded)"
        except OSError:
            return upload, os.path.basename(upload)
    return None, ""


def path_status(path_text: str | None) -> str:
    """One line for the GUI under a path box: the file size, or what is wrong."""
    if not (path_text or "").strip():
        return ""
    try:
        return "✔ " + resolve_input(None, path_text)[1]
    except ValueError as e:
        return f"✖ {e}"


def paths_from_query(params) -> tuple[str, str]:
    """Pre-fill the path boxes from the page URL: http://127.0.0.1:7860/?protein=PATH&dna=PATH (URL-encoded).
    This is the hook for desktop wrappers. Ignored on the hosted Space."""
    if not path_inputs_enabled() or not params:
        return "", ""
    get = params.get
    return str(get("protein") or ""), str(get("dna") or "")


# No note when running locally (that nothing leaves the computer is said in docs/install.md, desktop/README.md and the
# site). The hosted note matters: there, files do go to a shared server, and people need to know.
INTRO_LOCAL = ""
INTRO_HOSTED = ("This copy runs on a shared Hugging Face server: files you upload are processed there and "
                f"deleted within {RUNS_MAX_AGE_HOURS} hours (or when the Space restarts). For unpublished or "
                "sensitive data, install COAL on your own computer (pip install) instead.")


def intro_text() -> str:
    return INTRO_HOSTED if hosted() else INTRO_LOCAL


# ---------------------------------------------------------------- engine features

@dataclass(frozen=True)
class Features:
    output_formats: tuple[str, ...] = ("fasta",)
    report: bool = False
    json: bool = False
    keep_going: bool = False
    version: str = "unknown"


def _coal_cmd() -> list[str]:
    return [sys.executable, "-m", "coal"]


@lru_cache(maxsize=1)
def detect_features() -> Features:
    """Read ``python -m coal --help`` to find which options this engine has."""
    try:
        # stdin=DEVNULL: see _start_engine
        h = subprocess.run(_coal_cmd() + ["--help"], capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=60, stdin=subprocess.DEVNULL).stdout
        v = subprocess.run(_coal_cmd() + ["--version"], capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=60, stdin=subprocess.DEVNULL)
        version = (v.stdout or v.stderr).strip().replace("coal ", "") or "unknown"
    except (OSError, subprocess.SubprocessError):
        return Features()
    fmts = ("fasta",)
    m = re.search(r"-output(?:,\s*--[\w-]+)?\s+\{([\w,]+)\}", h)
    if m:
        fmts = tuple(f for f in m.group(1).split(",") if f)
    return Features(
        output_formats=fmts,
        report="--report" in h,
        json=re.search(r"(^|\s)--json\b", h) is not None,
        keep_going="--keep-going" in h,
        version=version,
    )


# ---------------------------------------------------------------- small file helpers

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_ids(path: str) -> list[str]:
    """Sequence IDs from a FASTA or CLUSTAL file (first word of each name), in file order.
    Read by the engine's rules, via seqformat."""
    from .seqformat import read_ids as _read_ids
    return _read_ids(path)


def strip_version(seq_id: str) -> str:
    """AT1G01010.1 -> AT1G01010 (removes a trailing .N version suffix)."""
    return re.sub(r"\.\d+$", "", seq_id)


# ---------------------------------------------------------------- commands

def _q(p: str) -> str:
    return shlex.quote(p) if os.name != "nt" else (f'"{p}"' if " " in p else p)


def strip_folders(text: str, paths) -> str:
    """Show each path in a command by its file name only, however it was quoted.
    The engine quotes with shlex.quote on every OS, so on Windows a full path (backslashes, drive colon) arrives as
    'C:\\...\\protein.aln.fasta' and removing only the folder would leave the quotes: 'protein.aln.fasta'. Each
    quoted form of the full path is replaced by the file name quoted the app's way."""
    for p in sorted({p for p in paths if p}, key=len, reverse=True):
        base = os.path.basename(p)
        for form in (shlex.quote(p), _q(p)):
            if form != p:
                text = text.replace(form, _q(base))
        text = text.replace(p, base)
    return text


def coal_command(protein: str, dna: str, output_format="fasta", codontable=1, nogap=False,
                 nomismatch=False, out="codon_aln.fasta", keep_going=False, report: str | None = None,
                 strict=False, blockonly=False) -> str:
    parts = ["coal", _q(protein), _q(dna)]
    if output_format != "fasta":
        parts += ["-output", output_format]
    if int(codontable) != 1:
        parts += ["-codontable", str(int(codontable))]
    if nogap:
        parts.append("-nogap")
    if nomismatch:
        parts.append("-nomismatch")
    if blockonly:
        parts.append("-blockonly")
    if strict:
        parts.append("--strict")
    if keep_going:
        parts.append("--keep-going")
    if report:
        parts += ["--report", _q(report)]
    parts += ["-o", _q(out)]
    return " ".join(parts)


def pal2nal_command(protein: str, dna: str, output_format="fasta", codontable=1, nogap=False,
                    nomismatch=False, out="codon_aln.fasta", blockonly=False) -> str:
    """The closest PAL2NAL v14 command line (see pal2nal_caveats for when it is not the same run)."""
    parts = ["pal2nal.pl", _q(protein), _q(dna), "-output", output_format]
    if int(codontable) != 1:
        parts += ["-codontable", str(int(codontable))]
    if nogap:
        parts.append("-nogap")
    if nomismatch:
        parts.append("-nomismatch")
    if blockonly:
        parts.append("-blockonly")
    return " ".join(parts) + f" > {_q(out)}"


def pal2nal_caveats(opts: dict) -> list[str]:
    """When the closest PAL2NAL v14 command does not do the same as this COAL run."""
    n = int(opts.get("codontable") or 1)
    out = []
    old = {3: "accepts only ATA/ATG as start codons (current NCBI table 3 also GTG)",
           10: "reads TGA as a stop codon (current NCBI table 10: Cys)",
           11: "lacks the start codons ATT, ATC and ATA of the current NCBI table 11",
           13: "lacks the start codon ATA of the current NCBI table 13"}
    if n in old:
        out.append(f"Genetic code {n}: PAL2NAL v14's table {old[n]}, so the two tools can differ (warnings, "
                   f"and the columns removed by -nomismatch).")
    if opts.get("nogap") and n != 1:
        out.append("-nogap: PAL2NAL v14 removes the standard code's stop codons (TAA/TAG/TGA) whatever "
                   "-codontable says, so its output differs from COAL's, which removes this code's stop codons.")
    if opts.get("strict"):
        out.append("--strict is COAL's own option: PAL2NAL v14 has no way to refuse probably-wrong genes.")
    return out


def commands_text(coal_cmd: str, pal2nal_cmd: str, opts: dict) -> str:
    """The command box of the GUI."""
    cav = pal2nal_caveats(opts or {})
    notes = "".join(f"# Note: {c}\n" for c in cav)
    return (f"# COAL\n{coal_cmd}\n\n# Closest PAL2NAL v14 command (COAL's documented deliberate differences from "
            f"PAL2NAL apply)\n{notes}{pal2nal_cmd}")


# ---------------------------------------------------------------- running

@dataclass
class RunResult:
    exit_code: int
    output_path: str | None
    report_path: str
    report: dict
    stderr: str
    coal_cmd: str
    pal2nal_cmd: str
    hints: list[dict] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        # exit 2 is also argparse's usage error: see run_completed()
        return run_completed(self.exit_code, self.stderr, (self.report.get("options") or {}).get("keep_going"))\
            and self.output_path is not None


def run_completed(code, stderr: str = "", keep_going=False) -> bool:
    """True if COAL finished: exit 0, or exit 2 only from a real --keep-going run (never argparse's usage error,
    which also exits 2)."""
    if code == 0:
        return True
    return (code == 2 and bool(keep_going) and detect_features().keep_going
            and "usage:" not in (stderr or ""))


def was_stopped(code, stderr: str = "") -> bool:
    """Killed by a signal / Ctrl-C (exit < 0 or >= 128, or 'coal: interrupted')."""
    return code is not None and (code < 0 or code >= 128 or "coal: interrupted" in (stderr or ""))


STOPPED_TEXT = ("The run was stopped before it finished (COAL was interrupted, killed or ran out of time), "
                "so the sequences were not all checked and no alignment was written.")


_WARN = re.compile(r"^coal: warning: (\S+) pepAlnPos (\d+): (\S+) does not correspond to (\S+)")
# PAL2NAL's "X unknown AA type" wording and the engine's "U (selenocysteine) is not checked; any codon is accepted"
_WARN_X = re.compile(r"^coal: warning: (\S+) pepAlnPos (\d+): (\S+) (?:unknown AA type|\([^)]*\) is not checked)")
# The colon after the ID is stripped ("s0:" -> "s0"); a line naming no known sequence becomes a run-level hint.
_WARN_ANY = re.compile(r"^coal: warning: (\S+?):? (.+)$")  # any other per-sequence warning
# both "(K of N do not match)" and "(K of N checkable codons do not match; A contain N or other ambiguous bases)"
_WARN_NOANCHOR = re.compile(r"^coal: warning: (\S+?): no 10-residue segment of the protein was found in the DNA; "
                            r"codons were assigned by position \((\d+) of (\d+) (?:checkable codons )?do not match"
                            r"(?:; (\d+) contain N or other ambiguous bases)?\)")
_WARN_PAIRING = re.compile(r"^coal: warning: only (\d+) of (\d+) protein IDs are DNA IDs; (.*)$")
# Without --strict the engine writes these genes (as PAL2NAL does) and prints these notices; the report marks such
# genes "failed" (probably wrong data).
_WARN_UNRELATED = re.compile(r"^coal: warning: (\S+?): protein and DNA are probably unrelated: (\d+) of (\d+) "
                             r"checkable codons do not match")
_WARN_GUESS = re.compile(r"^coal: warning: (\S+?): position of the gene in the DNA is a guess: .*?(?:only (\d+) of "
                         r"(\d+) codons have no ambiguous bases|the gene has only (\d+) codons)")
SUSPECT_CODES = ("suspect_unrelated", "suspect_position")  # a gene with one of these is reported as "failed"
PROBLEM_HINTS = ("partial_id_match", "engine_warning", "nomismatch_unknown")  # hints that mean "something is off"


def _parse_warning(line: str) -> tuple[str, dict] | None:
    """(sequence id, message) for any 'coal: warning: <id> ...' line, else None."""
    m = _WARN.match(line)
    if m:
        sid, col, aa, codon = m.groups()
        return sid, {"code": "mismatch", "column": int(col),
                     "text": f"Column {col}: amino acid {aa} does not match codon {codon}.",
                     "residue": aa, "codon": codon}
    m = _WARN_X.match(line)
    if m:
        sid, col, aa = m.groups()
        what = {"U": "selenocysteine", "O": "pyrrolysine", "B": "D or N", "Z": "E or Q", "J": "I or L"}.get(aa)
        return sid, {"code": "unknown_aa", "column": int(col), "residue": aa,
                     "text": (f"Column {col}: '{aa}'" + (f" ({what})" if what else "") + " is not one of the 20 "
                              f"standard amino acids, so COAL did not check its codon (any codon accepted).")}
    m = _WARN_UNRELATED.match(line)
    # the --strict advice also names the app's option (the GUI and the MCP tool have it)
    if m:
        sid, nbad, nchk = m.groups()
        return sid, {"code": "suspect_unrelated",
                     "text": (f"{nbad} of the {nchk} codons COAL could check do not match this protein, so this is "
                              f"almost certainly the wrong DNA for it (or the wrong genetic code). COAL still wrote "
                              f"its row, codon by position, as PAL2NAL does: do not trust it. Run COAL with --strict "
                              f"(in the app: 'Refuse probably-wrong genes') to refuse such genes.")}
    m = _WARN_GUESS.match(line)
    if m:
        sid, nchk, nres, short = m.groups()
        why = (f"Only {nchk} of this gene's {nres} codons have no unknown bases (N), and no "
               if short is None else f"The gene has only {short} codons, and no ")
        return sid, {"code": "suspect_position",
                     "text": (f"{why}"
                              f"10-residue stretch of the protein was found in the DNA, so COAL cannot tell where "
                              f"the gene starts in the longer DNA. Its row was taken from the start of the DNA (as "
                              f"PAL2NAL does) and may be shifted: do not trust it. Run COAL with --strict (in the app: "
                              f"'Refuse probably-wrong genes') to refuse such genes.")}
    m = _WARN_NOANCHOR.match(line)
    if m:
        sid, nbad, nres, namb = m.groups()
        amb = f"; {namb} contain N or other ambiguous bases" if namb else ""
        # when every checkable codon matches and the rest are ambiguous, the placement is fine and the code is not
        # the problem
        if namb and int(nbad) == 0:
            ask = "Every codon COAL could check matches, so the placement looks right."
        else:
            ask = "Is this the right DNA for this protein, and the right genetic code?"
        return sid, {"code": "no_anchor",
                     "text": (f"No 10-residue stretch of the protein was found in the DNA, so codons were assigned "
                              f"by position ({nbad} of {nres} checkable codons do not match{amb}). {ask}")}
    if _WARN_PAIRING.match(line):
        return None  # run-level; see pairing_warning()
    m = _WARN_ANY.match(line)
    if m:
        sid, rest = m.groups()
        return sid, {"code": "warning", "text": f"COAL warning: {rest.strip()}"}
    return None
_ERR = re.compile(r"^coal: error: (.*)$")


def pairing_warning(line: str) -> dict | None:
    """The engine's line 'only K of N protein IDs are DNA IDs; ...' as a hint."""
    m = _WARN_PAIRING.match(line)
    if not m:
        return None
    k, n, rest = m.groups()
    ex = re.search(r"protein '([^']*)' gets DNA '([^']*)'", rest)
    eg = f" (for example protein '{ex.group(1)}' got DNA '{ex.group(2)}')" if ex else ""
    return {"code": "partial_id_match",
            "text": (f"Only {k} of {n} protein IDs were found in the DNA file, so COAL paired the sequences by "
                     f"their order in the files instead{eg}. Some proteins are probably paired with the wrong "
                     f"DNA: make the names match (or put both files in the same order) and run again.")}

_EXT = {"fasta": "fasta", "clustal": "aln", "paml": "paml", "codon": "codon.txt"}


def _run_engine(protein, dna, out, fmt, table, nogap, nomismatch, keep_going, report_path, timeout=None,
                strict=False, blockonly=False):
    f = detect_features()
    timeout = timeout or ENGINE_TIMEOUT
    cmd = _coal_cmd() + ["-o", out, "-codontable", str(int(table))]
    if fmt != "fasta":
        cmd += ["-output", fmt]
    if nogap:
        cmd.append("-nogap")
    if nomismatch:
        cmd.append("-nomismatch")
    if blockonly:
        cmd.append("-blockonly")
    if strict:
        cmd.append("--strict")
    if keep_going and f.keep_going:
        cmd.append("--keep-going")
    if report_path and f.report:
        cmd += ["--report", report_path]
    cmd += ["--", protein, dna]  # after "--", so a file named '-p.fa' is not read as an option
    t0 = time.perf_counter()
    p = _start_engine(cmd)
    try:
        _, err = p.communicate(timeout=timeout)
        code = p.returncode
    except subprocess.TimeoutExpired:
        stop_engine(p)
        _remove_engine_temp(os.path.dirname(os.path.abspath(out)))
        code, err = 1, f"coal: error: COAL did not finish within {timeout:g} seconds\n"
    except BaseException:  # e.g. KeyboardInterrupt in the app: never leave the engine running
        stop_engine(p)
        raise
    return code, err, time.perf_counter() - t0


def _start_engine(cmd: list) -> subprocess.Popen:
    """The engine as a child with piped output.
    stdin=DEVNULL: in the desktop app the app's stdin is the wrapper's lifeline pipe, which serve() keeps a
    blocking read on; on Windows a child given that same (synchronous) pipe blocks at startup behind the pending
    read, so every engine run would hang. The engine never reads stdin.
    POSIX start_new_session: the engine leads its own process group, so stop_engine's SIGKILL can take its worker
    processes with it (see there). The app's Ctrl-C then no longer reaches the engine directly; _run_engine and
    run_engine_command stop it themselves on KeyboardInterrupt."""
    return subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            encoding="utf-8", errors="replace", start_new_session=os.name != "nt")


def _kill_engine_tree(p: subprocess.Popen) -> None:
    """SIGKILL the engine and, when it leads its own process group (started by _start_engine), its workers too."""
    if os.name != "nt":
        try:
            if os.getpgid(p.pid) == p.pid:
                os.killpg(p.pid, signal.SIGKILL)
                return
        except OSError:  # already gone (or not ours to signal): fall through to the plain kill
            pass
    p.kill()


def stop_engine(p: subprocess.Popen, grace: float = 5.0) -> None:
    """Stop the engine and its workers: POSIX SIGTERM (the engine cleans up, stops its workers, exits 143), wait
    `grace` seconds, then SIGKILL the engine's process group; Windows: taskkill /T /F (the whole process tree),
    then kill(). A plain SIGKILL of the engine alone would orphan its worker processes (they wait forever for
    work from a parent that is gone) and leave its .coal_*.tmp files behind. The engine's workers also hold its
    stdout/stderr pipes, so after the kill the pipes are read for at most `grace` more seconds and then closed:
    an orphan that survived the kill must never make the app wait forever."""
    if p.poll() is not None:
        p.communicate()
        return
    if os.name == "nt":
        try:
            # stdin=DEVNULL: see _start_engine
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(p.pid)], capture_output=True, timeout=30,
                           stdin=subprocess.DEVNULL)
        except (OSError, subprocess.SubprocessError):
            pass
    else:
        p.terminate()
    try:
        p.communicate(timeout=grace)
        return
    except subprocess.TimeoutExpired:
        pass
    _kill_engine_tree(p)
    try:
        p.communicate(timeout=grace)
    except subprocess.TimeoutExpired:  # something still holds the pipes: stop reading them
        for fh in (p.stdout, p.stderr):
            if fh is not None:
                fh.close()
        p.wait()


def _remove_engine_temp(folder: str) -> None:
    """After a forced stop: remove the engine's unfinished .coal_*.tmp files in our own run folder."""
    try:
        names = os.listdir(folder)
    except OSError:
        return
    for n in names:
        if n.startswith(".coal_") and n.endswith(".tmp"):
            try:
                os.remove(os.path.join(folder, n))
            except OSError:
                pass


def run_engine_command(cmd: list, timeout: float, temp_dir: str | None = None) -> tuple[int, str, str]:
    """For the demo scripts (subprocess.run(timeout=) would SIGKILL only the engine: workers orphaned, .coal_*.tmp
    left). Runs `cmd`, returns (returncode, stdout, stderr); on timeout or interrupt stops the engine as stop_engine
    does, removes its .coal_*.tmp files in `temp_dir` and re-raises (subprocess.TimeoutExpired on timeout)."""
    p = _start_engine(cmd)
    try:
        out, err = p.communicate(timeout=timeout)
    except BaseException:
        stop_engine(p)
        if temp_dir:
            _remove_engine_temp(temp_dir)
        raise
    return p.returncode, out, err


def public_error(e, paths=()) -> str:
    """An error message for the user with every given path shown as its file name only (no server folders, e.g.
    on the hosted Space)."""
    text = str(e)
    for p in sorted((p for p in paths if p), key=len, reverse=True):
        text = text.replace(p, os.path.basename(p))
    return text


def _classify_error(text: str) -> str:
    t = text.lower()
    if "inconsistent" in t or "no match" in t:
        return "inconsistent"
    if "number of sequences differs" in t:
        return "count_mismatch"
    if "invalid protein" in t:
        return "invalid_protein"
    if "empty" in t:
        return "empty"
    if "no such file" in t:
        return "file_not_found"
    return "error"


def build_fallback_report(protein, dna, out, fmt, table, nogap, nomismatch, keep_going,
                          code, stderr, runtime, strict=False, blockonly=False) -> dict:
    """Build a coal.report/1 report from COAL's stderr (for engines without --report)."""
    try:
        pids = read_ids(protein)
    except OSError:
        pids = []
    try:
        dids = read_ids(dna)
    except OSError:
        dids = []
    dset = set(dids)
    if pids and all(p in dset for p in pids):
        pairing = "id"
    else:
        pairing = "order"
    per: dict[str, list[dict]] = {p: [] for p in pids}
    failed: dict[str, list[dict]] = {}
    general: list[dict] = []
    run_hints: list[dict] = []
    for line in stderr.splitlines():
        pw = pairing_warning(line)
        if pw:
            run_hints.append(pw)
            continue
        w = _parse_warning(line)
        if w and w[0] in per:
            per[w[0]].append(w[1])
            continue
        if w or line.startswith("coal: warning:"):
            run_hints.append({"code": "engine_warning", "text": "COAL warning: " + line[len("coal: warning:"):].strip()})
            continue
        m = _ERR.match(line)
        if m:
            msg = m.group(1)
            sid, _, rest = msg.partition(": ")
            if rest and sid in per:
                failed[sid] = [{"code": _classify_error(rest), "text": rest}]
            else:
                general.append({"code": _classify_error(msg), "text": msg})
    aborted = not run_completed(code, stderr, keep_going)
    if was_stopped(code, stderr) and not any(g.get("code") == "stopped" for g in general):
        general.insert(0, {"code": "stopped", "text": STOPPED_TEXT})
    elif aborted and not general and not failed:
        general.append({"code": "error", "text": "COAL stopped without saying why (exit code "
                                                 f"{code}); no alignment was written."})
    if aborted and general and not failed:
        pairing = None  # COAL stopped before pairing anything (e.g. counts differ, unreadable input)
    keep_going = bool(keep_going and detect_features().keep_going)  # only what was really passed
    seqs = []
    # In a stopped run, genes before the failing one were processed: they say "Not written ..." (with their own
    # warnings); genes after it name the failing gene.
    first_fail = next((i for i, p in enumerate(pids) if p in failed), None) if aborted else None
    for i, p in enumerate(pids):
        if pairing is None:
            d = None
        elif pairing == "id":
            d = p
        else:
            d = dids[i] if i < len(dids) else None
        # In a run where only some protein IDs are DNA IDs, a gene paired by order with a DNA of another name says
        # so; it is "failed" if it also has mismatches.
        msgs_p = list(per.get(p) or [])
        mispaired = (pairing == "order" and d is not None and d != p
                     and any(h.get("code") == "partial_id_match" for h in run_hints))
        if mispaired and not aborted and p not in failed:
            msgs_p.append({"code": "suspect_pairing", "dna_id": d,
                           "text": (f"Paired with DNA '{d}' only by its position in the files"
                                    + (" (this protein's own name is also in the DNA file)" if p in dset else "")
                                    + ": check that this is the right DNA.")})
        if p in failed:
            status, msgs = "failed", failed[p]
        elif aborted:
            status, msgs = "not_checked", []
            if first_fail is not None:
                fid = pids[first_fail]
                if i < first_fail:
                    msgs = [{"code": "not_written", "failed_id": fid}] + list(per.get(p) or [])
                else:
                    msgs = [{"code": "not_reached", "failed_id": fid}]
        elif any(m.get("code") in SUSPECT_CODES for m in msgs_p) or (  # probably wrong data
                mispaired and any(m.get("code") == "mismatch" for m in msgs_p)):
            status = "failed"
            msgs = sorted(msgs_p, key=lambda m: (m.get("code") not in SUSPECT_CODES,
                                                 m.get("code") != "suspect_pairing"))
        elif msgs_p:
            status, msgs = "warning", msgs_p
        else:
            status, msgs = "ok", []
        seqs.append({"id": p, "dna_id": d, "pairing": pairing, "status": status, "messages": msgs})
    counts = {s: sum(1 for x in seqs if x["status"] == s) for s in STATUS}
    produced = (not aborted) and os.path.exists(out)
    rep = {
        "schema": "coal.report/1",
        "coal_version": detect_features().version,
        "source": "coal_app (built from COAL's messages; this engine has no --report yet)",
        "command": coal_command(protein, dna, fmt, table, nogap, nomismatch, out, keep_going,
                                strict=strict, blockonly=blockonly),
        "pal2nal_equivalent": pal2nal_command(protein, dna, fmt, table, nogap, nomismatch, out, blockonly),
        "options": {"output_format": fmt, "codontable": int(table), "nogap": bool(nogap),
                    "nomismatch": bool(nomismatch), "keep_going": bool(keep_going), "strict": bool(strict),
                    "blockonly": bool(blockonly)},
        "inputs": [
            {"role": "protein", "path": protein, "sha256": _safe_sha(protein), "sequences": len(pids)},
            {"role": "dna", "path": dna, "sha256": _safe_sha(dna), "sequences": len(dids)},
        ],
        "output": ({"path": out, "format": fmt, "sha256": sha256_file(out), "sequences": len(pids)}
                   if produced else None),
        "summary": {"total": len(pids), "ok": counts["ok"], "warning": counts["warning"],
                    "failed": counts["failed"], "pairing": pairing},
        "hints": run_hints,
        "errors": general,
        "sequences": seqs,
        "runtime_seconds": round(runtime, 3),
        "exit_code": code,
    }
    if counts["not_checked"]:
        rep["summary"]["not_checked"] = counts["not_checked"]
    return rep


def _safe_sha(p):
    try:
        return sha256_file(p)
    except OSError:
        return None


def id_hints(protein: str, dna: str, report: dict | None = None) -> list[dict]:
    """Spot IDs that only differ by a version suffix (AT1G01010.1 vs AT1G01010).

    With a report: says the order pairing is probably wrong only if the run had failures, warnings or
    stopped; if every pair matched codon for codon it says so and only suggests renaming."""
    try:
        pids, dids = read_ids(protein), read_ids(dna)
    except OSError:
        return []
    dset = set(dids)
    missing = [p for p in pids if p not in dset]
    if not missing:
        return []
    if any(h.get("code") == "partial_id_match" for h in (report or {}).get("hints") or []):
        return []  # the engine's own partial-match warning already says it (and that it's wrong)
    summ = (report or {}).get("summary") or {}
    if report is not None and summ.get("pairing") not in ("order", None):
        return []  # the engine paired them some other way (e.g. normalised IDs)
    clean = (report is not None and report.get("exit_code") == 0 and summ.get("total")
             and summ.get("ok") == summ.get("total"))
    if report is not None and summ.get("pairing") is None:
        verdict = "COAL stopped before pairing them."
    elif clean:
        verdict = ("Sequences were paired by their order in the files instead, and every pair matched codon "
                   "for codon, so the pairing looks right; renaming makes it safe.")
    else:
        verdict = ("Sequences were therefore paired by their order in the files, which is probably wrong "
                   "for the sequences that failed or have warnings.") if report is not None else (
                   "Sequences were therefore paired by their order in the files, which may be wrong.")
    dnorm: dict[str, str] = {}
    for d in dids:
        dnorm.setdefault(strip_version(d), d)
    fixable = [p for p in missing if strip_version(p) in dnorm]
    if fixable and len(fixable) == len(missing):
        ex = fixable[0]
        return [{"code": "id_mismatch",
                 "text": (f"The protein IDs don't exactly match the DNA IDs: they differ only by a version "
                          f"suffix (for example protein '{ex}' vs DNA '{dnorm[strip_version(ex)]}'). {verdict} "
                          f"Rename the sequences so the IDs match exactly, for example by removing the '.1'-style "
                          f"version suffixes from the protein alignment ('{ex}' -> '{strip_version(ex)}'), then "
                          f"run again.")}]
    return [{"code": "id_mismatch",
             "text": (f"{len(missing)} protein ID(s) were not found in the DNA file (for example '{missing[0]}'). "
                      f"{verdict} Check that both files describe the same genes, in the same order.")}]


CANDIDATE_CODES = (2, 3, 4, 5, 9, 11, 13, 14, 21)


def _write_sample(protein, dna, ids, report, work_dir):
    """Protein alignment + DNA FASTA with only the sequences in `ids` (kept in file order)."""
    from .viewer import read_alignment
    want = set(ids)
    prot = [(i, r) for i, r in read_alignment(protein) if i in want]
    dna_recs = read_alignment(dna, dna=True)  # DNA rules ('#' lines, letters only)
    pairs = {s["id"]: s.get("dna_id") for s in report.get("sequences", [])}
    dmap = dict(dna_recs)
    nuc = [(i, dmap.get(pairs.get(i) or i, "")) for i, _ in prot]
    if any(not seq for _, seq in nuc):
        return None
    ps = os.path.join(work_dir, "probe_sample.pep.fasta")
    ds = os.path.join(work_dir, "probe_sample.nuc.fasta")
    for path, recs in ((ps, prot), (ds, nuc)):
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.writelines(f">{i}\n{seq}\n" for i, seq in recs)
    return ps, ds


def annotate_mismatch_patterns(protein: str, report: dict) -> None:
    """Marks, in place, the mismatch patterns that are not a genetic-code question:
    a gene whose only definite mismatch is its first Met gets m["start_codon"] = True; a gene whose definite mismatches
    form a dense run (to the end of the gene, or a block of 4 or more) gets a "possible_frameshift" message (a missing
    or extra base, or an indel, then in-frame again or the gene ends)."""
    want = {sq["id"]: sq for sq in report.get("sequences") or [] if real_mismatches(sq)}
    if not want:
        return
    from .viewer import read_alignment
    try:
        rows = {i: r for i, r in read_alignment(protein) if i in want}
    except (OSError, ValueError):
        return
    for sid, sq in want.items():
        row = rows.get(sid)
        if row is None:
            continue
        index, k = {}, 0  # alignment column (1-based) -> residue index in the gene
        for c, ch in enumerate(row, 1):
            if ch not in "-.":
                index[c] = k
                k += 1
        mm = real_mismatches(sq)
        idx = sorted({index[m["column"]] for m in mm if m.get("column") in index})
        if not idx:
            continue
        if idx == [0] and all(m.get("residue", "").upper() == "M" for m in mm):
            for m in mm:
                m["start_codon"] = True
            continue
        runs, cur = [], [idx[0]]
        for i in idx[1:]:
            if i - cur[-1] <= 2:  # at most one matching codon in between
                cur.append(i)
            else:
                runs.append(cur)
                cur = [i]
        runs.append(cur)
        best = max(runs, key=len)
        span = best[-1] - best[0] + 1
        to_end = best[-1] >= k - 2
        if len(best) >= 3 and len(best) / span >= 0.6 and 2 * len(best) >= len(idx) and (to_end or len(best) >= 4):
            col_of = {v: c for c, v in index.items()}
            c0, c1 = col_of[best[0]], col_of[best[-1]]
            where = "to the end of the gene" if to_end else f"up to column {c1}"
            sq.setdefault("messages", []).append({
                "code": "possible_frameshift", "column": c0, "end_column": c1,
                "text": (f"The mismatches form one run, from column {c0} {where}: this looks like a missing or extra "
                         f"base (a frameshift) or a small insertion/deletion in the DNA near column {c0}, not a "
                         f"genetic-code problem. Check the DNA there (a sequencing or assembly error?).")})


def _code_distance(a: int, b: int) -> tuple[int, int]:
    """(codons with a different amino acid, codons with a different start flag) between two genetic codes."""
    from coal.tables import get_code
    ca, cb = get_code(a), get_code(b)
    return (sum(ca.code[x] != cb.code[x] for x in ca.code), len(ca.starts ^ cb.starts))


def _fits_now(sq: dict) -> bool:
    if sq.get("status") == "not_checked":  # processed before the error, without a definite mismatch
        return any(m.get("code") == "not_written" for m in sq.get("messages") or []) and not real_mismatches(sq)
    return gene_kind(sq) in ("ok", "unchecked", "other")


def _fixable(sq: dict) -> bool:
    """A gene another genetic code might fix: definite mismatches (not a run that looks like a frameshift), probably
    unrelated, excluded or never checked. Not: pairing-only, non-standard residues, ambiguous bases."""
    k = gene_kind(sq)
    if k == "mismatch":
        return not any(m.get("code") == "possible_frameshift" for m in sq.get("messages") or [])
    if k == "flagged":
        return any(m.get("code") == "suspect_unrelated" for m in sq.get("messages") or []) or bool(
            real_mismatches(sq))
    if k == "not_checked":  # a gene processed before the error is fixable only if it mismatched
        written = any(m.get("code") == "not_written" for m in sq.get("messages") or [])
        return bool(real_mismatches(sq)) if written else True
    return k == "excluded"


def genetic_code_hint(protein, dna, table, report, work_dir) -> dict | None:
    """If sequences have mismatches, see whether another genetic code fits more of them.
    Only genes with definite mismatches (not a frameshift-like run), probably-unrelated, excluded or unchecked genes
    ask for a hint (pairing-only, U/B/Z-only and ambiguous-placement genes would be "fixed" by any code). Every gene
    (up to PROBE_SAMPLE, problem genes first) is re-run under each candidate code, so a gene that fits the current
    code but not the candidate is counted too; a gene fits when it has no definite mismatch and no
    probably-unrelated notice; the text says what was tested; ties go to code 1 when it is tied, otherwise to the
    code closest to the current one, and other codes that fit as well are named; a start-codon-only problem says so.

    Each probe has a short time limit and all probes run within PROBE_BUDGET seconds; every probe file is deleted."""
    import shutil
    seqs = report.get("sequences", [])
    problems = [sq["id"] for sq in seqs if _fixable(sq)]
    if not seqs or not problems:
        return None
    try:
        if os.path.getsize(protein) + os.path.getsize(dna) > 50_000_000:
            return None  # keep the GUI responsive on very big inputs
    except OSError:
        return None
    pset = set(problems)
    tested = (problems + [sq["id"] for sq in seqs if sq["id"] not in pset])[:PROBE_SAMPLE]
    tset = set(tested)
    fits_now = sum(1 for sq in seqs if sq["id"] in tset and _fits_now(sq))
    table = int(table)
    probe_dir = tempfile.mkdtemp(prefix="probe_", dir=work_dir)
    results = []
    try:
        files = _write_sample(protein, dna, tested, report, probe_dir)
        if files is None:
            return None
        ps, ds = files
        t_end = time.perf_counter() + PROBE_BUDGET
        for c in CANDIDATE_CODES:
            left = t_end - time.perf_counter()
            if c == table or left <= 1:
                continue
            out = os.path.join(probe_dir, f"probe_{c}.fasta")
            code, err, _ = _run_engine(ps, ds, out, "fasta", c, False, False, True, None,
                                       timeout=min(PROBE_TIMEOUT, left))
            if os.path.exists(out):
                os.remove(out)
            if code != 0:
                continue
            bad = set()
            for w in map(_parse_warning, err.splitlines()):
                if w and ((w[1].get("code") == "mismatch" and not is_ambiguous_codon(w[1].get("codon")))
                          or w[1].get("code") in SUSPECT_CODES):
                    bad.add(w[0])
            fit = sum(1 for i in tested if i not in bad)
            fixed = sum(1 for i in problems if i in tset and i not in bad)
            if fixed and fit > fits_now:
                results.append((fit, c))
    finally:
        shutil.rmtree(probe_dir, ignore_errors=True)
    if not results:
        return None
    top = max(f for f, _ in results)
    tied = sorted((c for f, c in results if f == top),
                  key=lambda c: (_code_distance(table, c)[0], CANDIDATE_CODES.index(c)))
    # the standard code 1 first when it is among the tied codes (as the engine, coal/report.py)
    standard = 1 in tied and len(tied) > 1
    if standard:
        tied = [1] + [x for x in tied if x != 1]
    c = tied[0]
    n, total = len(tested), len(seqs)
    name = GENETIC_CODES[c]

    def do(k):
        return "none do" if k == 0 else ("only 1 does" if k == 1 else f"{k} do")
    if n == total:
        who = f"all {total} sequences match" if top == total else f"{top} of {total} sequences match"
        text = f"Under genetic code {c} ({name}), {who} their protein; under code {table}, {do(fits_now)}."
    else:
        which = (f"a sample of the {len(problems)} with problems" if len(problems) >= n
                 else f"all {len(problems)} with problems, then others")
        text = (f"COAL tested {n} of the {total} sequences ({which}): under "
                f"genetic code {c} ({name}), {top} of them match their protein; under code {table}, {do(fits_now)}.")
    starts = [m for sq in seqs if sq["id"] in pset for m in real_mismatches(sq)]
    if starts and all(m.get("start_codon") for m in starts):
        ex = starts[0]
        text = (f"The only mismatches are at the first codon of {len(pset)} sequence{'s' if len(pset) != 1 else ''} "
                f"(for example {ex.get('codon')} where the protein starts with M): an alternative start codon, which "
                f"code {c} accepts. " + text)
    if len(tied) > 1:
        others = " and ".join([", ".join(str(x) for x in tied[1:-1]), str(tied[-1])] if len(tied) > 2
                              else [str(tied[1])])
        closer = _code_distance(table, c)[0] < _code_distance(table, tied[1])[0]
        text += (f" Code{'s' if len(tied) > 2 else ''} {others} fit{'' if len(tied) > 2 else 's'} "
                 + ("as well, but code 1 is the standard code, used by most nuclear genes." if standard else
                    f"as well, but code {c} differs least from code {table}." if closer else
                    "these sequences too: choose the code of the organism and organelle the genes come from."))
    text += f" Re-run with code {c}?"
    # also the counts, under the names the engine's hint uses (coal/report.py), for code_fix
    return {"code": "genetic_code", "text": text, "suggested_codontable": c, "checked": n, "fit": top,
            "fit_current": fits_now}


def mismatch_probe(protein, dna, table, keep_going, report, work_dir, blockonly=False) -> dict:
    """With -nomismatch the engine (like PAL2NAL) prints no mismatch warnings, so the report would say every gene
    matched. This re-runs the same inputs without -nomismatch (output
    to a temporary file, deleted) to find which genes had mismatching codons, marks them "warning" (columns
    removed) and stores the number of removed codon columns in summary["columns_removed"]."""
    import shutil
    probe_dir = tempfile.mkdtemp(prefix="probe_nomismatch_", dir=work_dir)
    try:
        code, err, _ = _run_engine(protein, dna, os.path.join(probe_dir, "probe.fasta"), "fasta", table,
                                   False, False, keep_going, None, blockonly=blockonly)
    finally:
        shutil.rmtree(probe_dir, ignore_errors=True)
    if code != 0:
        report.setdefault("hints", []).append({
            "code": "nomismatch_unknown",
            "text": ("Columns were removed (-nomismatch), but COAL could not count which genes had mismatching "
                     "codons. Run once without -nomismatch to see them.")})
        return report
    per: dict[str, list[dict]] = {}
    cols: set[int] = set()
    real_cols: set[int] = set()  # columns with at least one definite (A/C/G/T) mismatch
    for line in err.splitlines():
        w = _parse_warning(line)
        if w and w[1].get("code") == "mismatch":
            m = dict(w[1], removed=True)
            per.setdefault(w[0], []).append(m)
            cols.add(m["column"])
            if not is_ambiguous_codon(m.get("codon")):
                real_cols.add(m["column"])
    for sq in report.get("sequences", []):
        mm = per.get(sq["id"])
        if mm and sq["status"] in ("ok", "warning"):
            sq["status"] = "warning"
            sq["messages"] = mm + list(sq.get("messages") or [])
    summ = report["summary"]
    for k in ("ok", "warning"):
        summ[k] = sum(1 for x in report.get("sequences", []) if x["status"] == k)
    summ["columns_removed"] = len(cols)
    summ["columns_removed_ambiguous"] = len(cols - real_cols)
    return report


def has_problem_hint(hints) -> bool:
    return any(isinstance(h, dict) and h.get("code") in PROBLEM_HINTS for h in hints or [])


def banner_head(ok: bool, summary: dict, hints=()) -> str:
    """"Every sequence matched" also needs no removed columns and no run-level problem (e.g. a partial ID
    match)."""
    if (ok and not summary.get("failed") and not summary.get("warning") and not summary.get("columns_removed")
            and not has_problem_hint(hints)):
        return "### ✔ Done: every sequence matched"
    if ok:
        return "### ⚠ Done, with some problems (see the table)"
    return "### ✖ COAL stopped: no alignment was written"


def run_coal(protein_path: str, dna_path: str, output_format: str = "fasta", codontable: int = 1,
             nogap: bool = False, nomismatch: bool = False, keep_going: bool = True,
             out_dir: str | None = None, probe_codes: bool = True, strict: bool = False,
             blockonly: bool = False) -> RunResult:
    """Run COAL once and return the output, a coal.report/1 report and hints.
    strict: COAL's --strict (refuse probably-wrong genes); blockonly: -blockonly."""
    f = detect_features()
    if output_format not in f.output_formats:
        raise ValueError(f"output format '{output_format}' is not supported by this COAL "
                         f"(supported: {', '.join(f.output_formats)})")
    if int(codontable) not in GENETIC_CODES:
        raise ValueError(f"unknown genetic code {codontable}")
    out_dir = out_dir or new_run_dir()  # pruned by prune_runs()
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(protein_path))[0].replace(".aln", "")
    out = os.path.join(out_dir, f"{stem}.codon.{_EXT.get(output_format, 'txt')}")
    report_path = os.path.join(out_dir, f"{stem}.report.json")
    if os.path.exists(report_path):
        os.remove(report_path)
    code, err, runtime = _run_engine(protein_path, dna_path, out, output_format, codontable,
                                     nogap, nomismatch, keep_going, report_path, strict=strict, blockonly=blockonly)
    report = None
    if f.report and os.path.exists(report_path):
        try:
            with open(report_path, encoding="utf-8") as fh:
                report = json.load(fh)
        except (OSError, ValueError):
            report = None
    if report is None:
        report = build_fallback_report(protein_path, dna_path, out, output_format, codontable, nogap,
                                       nomismatch, keep_going, code, err, runtime, strict=strict, blockonly=blockonly)
        if nomismatch and code == 0:
            mismatch_probe(protein_path, dna_path, codontable, keep_going, report, out_dir, blockonly=blockonly)
        from .viewer import locate_breaks
        locate_breaks(protein_path, dna_path, report)
        annotate_mismatch_patterns(protein_path, report)
        hints = list(report.get("hints") or []) + id_hints(protein_path, dna_path, report)
        if probe_codes:
            h = genetic_code_hint(protein_path, dna_path, codontable, report, out_dir)
            if h:
                hints.append(h)
        report["hints"] = hints
        with open(report_path, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2)
    produced = run_completed(code, err, keep_going) and os.path.exists(out)
    return RunResult(
        exit_code=code,
        output_path=out if produced else None,
        report_path=report_path,
        report=report,
        stderr=err,
        coal_cmd=coal_command(os.path.basename(protein_path), os.path.basename(dna_path), output_format,
                              codontable, nogap, nomismatch, os.path.basename(out),
                              keep_going and f.keep_going, strict=strict, blockonly=blockonly),
        pal2nal_cmd=pal2nal_command(os.path.basename(protein_path), os.path.basename(dna_path),
                                    output_format, codontable, nogap, nomismatch, os.path.basename(out), blockonly),
        hints=list(report.get("hints") or []),
    )


# ---------------------------------------------------------------- reading reports

def _validate_report(rep) -> None:
    """Raise ValueError unless `rep` has the shapes explain() and the GUI rely on."""
    def bad(what):
        raise ValueError(f"not a valid COAL report: {what}")
    if not isinstance(rep, dict):
        bad("the file does not hold a JSON object")
    if not str(rep.get("schema", "")).startswith("coal.report/"):
        bad("expected schema 'coal.report/1'")
    for key, typ, name in (("summary", dict, "an object"), ("options", dict, "an object"),
                           ("output", (dict, type(None)), "an object or null"), ("inputs", list, "a list"),
                           ("hints", list, "a list"), ("errors", list, "a list"), ("sequences", list, "a list")):
        if key in rep and not isinstance(rep[key], typ):
            bad(f"'{key}' should be {name}")
    for key in ("inputs", "hints", "errors", "sequences"):
        for i, item in enumerate(rep.get(key) or []):
            if not isinstance(item, dict):
                bad(f"'{key}' item {i + 1} should be an object")
    for i, sq in enumerate(rep.get("sequences") or []):
        if not isinstance(sq.get("messages", []), list) or not all(
                isinstance(m, dict) for m in sq.get("messages") or []):
            bad(f"'messages' of sequence {i + 1} should be a list of objects")
    for k, v in (rep.get("summary") or {}).items():
        if k != "pairing" and not isinstance(v, int):
            bad(f"summary '{k}' should be a whole number")


def load_report(path: str) -> dict:
    # utf-8-sig (a BOM is fine); every field explain() uses is type-checked (e.g. {"summary": []})
    with open(path, encoding="utf-8-sig") as fh:
        try:
            rep = json.load(fh)
        except ValueError as e:
            raise ValueError(f"{path} is not JSON ({e})") from None
    try:
        _validate_report(rep)
    except ValueError as e:
        raise ValueError(f"{path} is {e}") from None
    return rep


def load_sample_report() -> dict:
    return load_report(SAMPLE_REPORT)


def sample_banner(rep: dict) -> str:
    """GUI banner for the sample report: a clearly labelled demo sample (a real run of examples/mito with the wrong
    genetic code)."""
    return ("### Demo sample (not your data)\nThe 24 mitochondrial genes of the *Mitochondrial (wrong code)* example, run "
            "with genetic code 1, the wrong code for them: see the tip to re-run with code 2. "
            + summary_sentence(rep))


def _mismatch_sentence(msgs: list[dict]) -> str:
    mm = [m for m in msgs if m.get("code") == "mismatch"]
    others = " ".join(m.get("text", "") for m in msgs if m.get("code") != "mismatch")
    if not mm:
        return others or "COAL printed a warning for this sequence."
    # a frameshift-like run of mismatches gets no genetic-code advice
    fs = any(m.get("code") == "possible_frameshift" for m in msgs)
    return _mm_sentence(mm, code_advice=not fs) + (" " + others if others else "")


def is_ambiguous_codon(codon) -> bool:
    """True if a codon has a base other than A/C/G/T(U) (N, IUPAC codes, ...): COAL can't check it."""
    return bool(codon) and bool(str(codon).upper().replace("U", "T").strip("ACGT"))


def _mm_sentence(mm: list[dict], code_advice: bool = True) -> str:
    # codons with unknown/ambiguous bases get their own sentence and no genetic-code advice
    amb = [m for m in mm if is_ambiguous_codon(m.get("codon"))]
    real = [m for m in mm if not is_ambiguous_codon(m.get("codon"))]
    removed = any(m.get("removed") for m in mm)
    out = []
    if real:
        first, n = real[0], len(real)
        where = f" (first at column {first['column']})" if first.get("column") else ""
        detail = ""
        if first.get("residue") and first.get("codon"):
            detail = f": the protein says {first['residue']} but the DNA has {first['codon']}"
        # mismatches found by the -nomismatch probe say their columns were removed
        fate = ("Those columns were removed from the alignment (-nomismatch)" if first.get("removed")
                else "Kept in the alignment")
        if n == 1 and first.get("start_codon"):
            out.append(f"The first codon, {first.get('codon')}, is not a start codon of this genetic code, although "
                       f"the protein starts with M: probably an alternative start codon (some codes accept it). "
                       f"{fate}.")
        else:
            out.append(f"{n} codon{'s' if n != 1 else ''} {'do' if n != 1 else 'does'} not match "
                       f"{'their' if n != 1 else 'its'} amino acid{where}{detail}. {fate}"
                       + ("; check the genetic code." if code_advice else "."))
    if amb:
        first, n = amb[0], len(amb)
        where = f" (first at column {first['column']}: {first['codon']})" if first.get("column") else ""
        fate = ("Those columns were removed from the alignment (-nomismatch)." if removed
                else "Kept in the alignment.")
        out.append(f"{n} codon{'s' if n != 1 else ''} {'contain' if n != 1 else 'contains'} unknown or ambiguous "
                   f"bases (e.g. N){where}, so COAL could not check {'them' if n != 1 else 'it'} against the "
                   f"protein. {fate} This is missing sequence data (masked or low-coverage DNA), not a "
                   f"genetic-code problem.")
    return " ".join(out)


_FAIL_SENTENCES = {
    "inconsistent": ("The DNA does not spell this protein, even allowing for a few mismatches. "
                     "Usually a frameshift (a missing or extra base), the wrong DNA paired with this protein, "
                     "or the wrong genetic code."),
    "frameshift": "A base is missing or extra in the DNA (a frameshift), so later codons are out of step.",
    "missing_dna": "No DNA sequence was found for this protein.",
    "invalid_protein": "The protein alignment contains characters that are not amino acids.",
    "empty": "This protein sequence is empty.",
    "internal_stop": "The DNA has a stop codon in the middle of the gene.",
}


def sequence_sentence(seq: dict) -> str:
    st = seq.get("status")
    msgs = seq.get("messages") or []
    if st == "ok":
        return "Protein and DNA match codon for codon."
    if st == "not_checked":
        # genes before / after the gene that stopped the run
        mark = next((m for m in msgs if m.get("code") in ("not_written", "not_reached")), None)
        if mark is None:
            return "Not checked: COAL stopped at an earlier error."
        fid = mark.get("failed_id")
        if mark["code"] == "not_reached":
            return f"Not checked: COAL stopped at an error in {fid} before reaching this gene."
        own = [m for m in msgs if m is not mark]
        if not own:
            return (f"Not written: COAL stopped at an error in {fid}, so no alignment was written; this gene itself "
                    f"had no problem.")
        return f"Not written: COAL stopped at an error in {fid}, so no alignment was written. " + \
            _mismatch_sentence(own)
    if st == "warning":
        return _mismatch_sentence(msgs)
    if st == "failed":
        lead = [m for m in msgs if m.get("code") in SUSPECT_CODES + ("suspect_pairing",)]
        if lead:  # written but probably wrong; say why first
            return " ".join(m.get("text", "") for m in lead)
        for m in msgs:
            if isinstance(m.get("column"), int) and m.get("code") in ("frameshift", "inconsistent", "missing_dna"):
                tail = {"frameshift": " Check the DNA near there (a sequencing or assembly error?).",
                        "inconsistent": " Is this the right DNA for this protein, and the right genetic code?",
                        "missing_dna": " The DNA looks truncated."}[m["code"]]
                if seq.get("pairing") == "order" and seq.get("dna_id") and seq.get("dna_id") != seq.get("id"):
                    # paired by order: the wrong DNA is more likely than a truncated one
                    tail = (f" It was paired with DNA '{seq['dna_id']}' only by their order in the files (the IDs "
                            f"don't match), so this is most likely the wrong DNA for it: make the IDs match.")
                return m.get("text", "") + tail
        for m in msgs:
            if m.get("code") in _FAIL_SENTENCES:
                return _FAIL_SENTENCES[m["code"]]
        return " ".join(m.get("text", "") for m in msgs) or "COAL could not align this sequence."
    return " ".join(m.get("text", "") for m in msgs)


def health_rows(report: dict) -> list[list[str]]:
    """One row per sequence: [icon, word, id, plain-English sentence]."""
    rows = []
    for s in report.get("sequences", []):
        icon, word, _ = STATUS.get(s.get("status"), ("?", str(s.get("status")), ""))
        rows.append([icon, word, s.get("id", ""), sequence_sentence(s)])
    return rows


def _jump_links(seq: dict, limit: int = 12, preview: dict | None = None) -> str:
    """preview (viewer.view's info: ids, first, last): a column outside the drawn part is named
    without a link ("not in the preview"), so no link leads nowhere."""
    from .viewer import cell_id
    # only codes the viewer gives an anchor to (it outlines mismatched and frameshift cells)
    cols = sorted({m["column"] for m in seq.get("messages") or []
                   if isinstance(m.get("column"), int) and m.get("code") in ("mismatch", "frameshift")})
    if not cols:
        return ""
    sid = str(seq.get("id", ""))
    drawn = preview is None or not preview.get("preview") or sid in (preview.get("ids") or ())
    lo, hi = (preview or {}).get("first", 1), (preview or {}).get("last", 10 ** 9)

    # "Columns 1, 7, 36" (the number is the link, as in a page index)
    def one(c):
        if drawn and lo <= c <= hi:
            return f'<a href="#{cell_id(sid, c)}" aria-label="column {c:,}">{c:,}</a>'
        return f"{c:,}"
    links = ", ".join(one(c) for c in cols[:limit])
    more = f" and {len(cols) - limit} more" if len(cols) > limit else ""
    lead = "Column" if len(cols) == 1 else "Columns"
    if not drawn:
        return f'<br><span class="coal-small">{lead} {links}{more} (this gene is not in the alignment preview)</span>'
    inside = [c for c in cols[:limit] if lo <= c <= hi]
    if len(inside) < len(cols[:limit]):
        more += (f" ({'some are ' if inside else ''}outside the preview, which draws columns {lo:,}–{hi:,}; type a "
                 f"column in \u201cCentre the preview on column\u201d under the alignment to see it)")
    return f'<br><span class="coal-small">{lead} {links}{more}</span>'


def _closing_sentence(text: str) -> str:
    """The last sentence of a row's text when there is more than one ("" otherwise)."""
    text = (text or "").rstrip()
    i = text.rfind(". ", 0, len(text) - 1)
    return text[i + 2:] if i > 0 else ""


CODE_FIX_SHARE = 0.9  # share of the checked genes that must fit the suggested code


def code_fix(report: dict) -> dict | None:
    """When the genetic-code hint explains the run (all, or at least CODE_FIX_SHARE, of the
    checked sequences fit the suggested code, and more of them than fit the code used), what the result page
    says first: the code used, the code to use, the counts, and the other codes that fit too (the hint's own
    sentence). None otherwise. Reads the report's hints only; the report and the engine are unchanged."""
    cur = (report.get("options") or {}).get("codontable") if isinstance(report.get("options"), dict) else None
    for h in report.get("hints") or []:
        sug = h.get("suggested_codontable")
        try:
            checked, fit, fit_cur = int(h.get("checked") or 0), int(h.get("fit") or 0), int(h.get("fit_current") or 0)
        except (TypeError, ValueError):
            continue
        if h.get("code") != "genetic_code" or not sug or not checked or sug == cur:
            continue
        if fit < CODE_FIX_SHARE * checked or fit <= fit_cur:
            continue
        m = re.search(r"Codes? [\d, and]+ fits? [^.]*\.", h.get("text") or "")
        return {"current": cur, "code": int(sug), "checked": checked, "fit": fit, "fit_current": fit_cur,
                "others": m.group(0) if m else ""}
    return None


def _plural_genes(n: int) -> str:
    return "1 gene" if n == 1 else f"{n:,} genes"


def _code_fix_genes(seqs: list[dict]) -> list[dict]:
    """The warnings a wrong genetic code explains (a warning with a mismatched codon)."""
    return [s for s in seqs if s.get("status") == "warning"
            and any(m.get("code") == "mismatch" for m in (s.get("messages") or []))]


HEALTH_MAX_ROWS = 200
HEALTH_GROUP_SHOWN = 3  # rows shown of a group that shares its closing sentence; the rest on demand
HEALTH_OK_NAMES = 2000  # OK genes are named in one collapsed row, at most this many


def health_table_html(report: dict, links: bool = False, preview: dict | None = None) -> str:
    """Accessible HTML table (icon + word + sentence, never colour alone).
    With links=True, problem columns link to the alignment viewer's cells (same page). No links when
    -nogap/-nomismatch removed columns: report columns then don't exist in the output.
    preview: what the alignment preview drew (viewer.view's info); only cells in it are linked."""
    opts = report.get("options") if isinstance(report.get("options"), dict) else {}
    if opts.get("nogap") or opts.get("nomismatch"):  # no jump links in codon-only mode
        links = False
    order = {"failed": 0, "warning": 1, "not_checked": 2, "ok": 3}
    allseqs = sorted(report.get("sequences", []), key=lambda s: order.get(s.get("status"), 4))
    if not allseqs:  # no report yet, or an input error: no empty table
        return ""
    # OK genes without a note are one collapsed row listing their names, not one row each
    oks = [s for s in allseqs if s.get("status") == "ok" and not (s.get("dna_id") and s.get("dna_id") != s.get("id"))]
    ok_ids = {id(s) for s in oks}
    allseqs = [s for s in allseqs if id(s) not in ok_ids]
    # When the genetic-code hint explains the run (code_fix), the warnings it explains are one
    # closed row after the other problems ("24 genes with mismatched codons under code 1: expected with the wrong
    # code"), not 24 rows
    fix = code_fix(report)
    fixed = _code_fix_genes(allseqs) if fix else []
    if fixed:
        fixed_ids = {id(s) for s in fixed}
        allseqs = [s for s in allseqs if id(s) not in fixed_ids]
    # at most HEALTH_MAX_ROWS rows (one row per sequence: 22,000 rows = 4.2 MB of HTML), problems first, then one
    # summary row for the rest
    seqs = allseqs[:HEALTH_MAX_ROWS]
    rest = allseqs[HEALTH_MAX_ROWS:]
    caption = "Sequences"
    out = [f'<table class="coal-health"><caption>{caption}</caption>',
           '<thead><tr><th scope="col">Status</th><th scope="col">Sequence</th>'
           '<th scope="col">What this means</th></tr></thead><tbody>']
    # Three or more problem genes with the same status and the same sentence, and nothing of their
    # own (no jump link, no pairing note), are one row naming them (e.g. the 19 "Not checked: COAL stopped at an
    # error in X" rows after a fatal error).
    groups: dict = {}
    for s in seqs:
        paired = s.get("dna_id") and s.get("dna_id") != s.get("id")
        if not paired and not (links and _jump_links(s, preview=preview)):
            groups.setdefault((s.get("status"), sequence_sentence(s)), []).append(s)
    grouped = {id(x) for g in groups.values() if len(g) >= 3 for x in g}
    # A closing sentence that three or more rows of the same status share (e.g. "Kept in the
    # alignment; check the genetic code." under 24 warnings) is said once, in a row above the first of them
    tails: dict = {}
    for s in seqs:
        if id(s) not in grouped:
            t = _closing_sentence(sequence_sentence(s))
            if t:
                tails.setdefault((s.get("status"), t), []).append(s)
    shared = {id(x): k[1] for k, g in tails.items() if len(g) >= 3 for x in g}
    # Such a group shows its first HEALTH_GROUP_SHOWN rows, then the others under one "Show N more" disclosure (a
    # native <details>, as gr.HTML runs no scripts) right after them, each with its own jump links
    later = {}
    for g in tails.values():
        if len(g) > HEALTH_GROUP_SHOWN:
            later[id(g[HEALTH_GROUP_SHOWN - 1])] = g[HEALTH_GROUP_SHOWN:]
    hidden = {id(x) for g in later.values() for x in g}
    told = set()
    done = set()

    def row(s, sentence):
        icon, word, cls = STATUS.get(s.get("status"), ("?", str(s.get("status")), ""))
        pair = ""
        if s.get("dna_id") and s.get("dna_id") != s.get("id"):
            pair = f'<br><span class="coal-small">paired with DNA {html.escape(str(s["dna_id"]))}</span>'
        return (f'<tr class="coal-{cls}"><td><span class="coal-badge coal-{cls}" aria-hidden="true">{icon}</span> '
                f'<strong>{word}</strong></td><td><code>{html.escape(str(s.get("id", "")))}</code>{pair}</td>'
                f'<td>{html.escape(sentence)}{_jump_links(s, preview=preview) if links else ""}</td></tr>')

    def short(s):
        sentence = sequence_sentence(s)
        return sentence[:-len(shared[id(s)])].rstrip() if id(s) in shared else sentence
    for s in seqs:
        if id(s) in hidden:
            continue
        icon, word, cls = STATUS.get(s.get("status"), ("?", str(s.get("status")), ""))
        if id(s) in grouped:
            key = (s.get("status"), sequence_sentence(s))
            if key in done:
                continue
            done.add(key)
            g = groups[key]
            names = " ".join(f"<code>{html.escape(str(x.get('id', '')))}</code>" for x in g)
            out.append(
                f'<tr class="coal-{cls}"><td><span class="coal-badge coal-{cls}" aria-hidden="true">{icon}</span> '
                f'<strong>{word}</strong></td><td colspan="2">{html.escape(key[1])}<details><summary>{len(g):,} genes'
                f'</summary><p class="coal-oklist">{names}</p></details></td></tr>')
            continue
        if id(s) in shared:
            tail = shared[id(s)]
            if (s.get("status"), tail) not in told:
                told.add((s.get("status"), tail))
                n = len(tails[(s.get("status"), tail)])
                out.append(f'<tr class="coal-{cls} coal-shared"><td colspan="3"><span class="coal-badge coal-{cls}" '
                           f'aria-hidden="true">{icon}</span> <strong>{n:,} genes:</strong> {html.escape(tail)}</td></tr>')
        out.append(row(s, short(s)))
        if id(s) in later:  # the rest of its group, on demand
            more = later[id(s)]
            out.append(f'<tr class="coal-morerows"><td colspan="3"><details><summary>Show {len(more):,} '
                       f'more</summary><table class="coal-health coal-nested"><thead class="coal-sr"><tr>'
                       f'<th scope="col">Status</th><th scope="col">Sequence</th><th scope="col">What this means</th>'
                       f'</tr></thead><tbody>{"".join(row(x, short(x)) for x in more)}</tbody></table></details>'
                       f'</td></tr>')
    if fixed:
        icon, word, cls = STATUS["warning"]
        shown = fixed[:HEALTH_MAX_ROWS]
        more = (f'<p class="coal-small">… and {len(fixed) - len(shown):,} more; the report JSON lists every '
                f'sequence.</p>' if len(fixed) > len(shown) else "")
        what = (f"{_plural_genes(len(fixed))} with mismatched codons under code {fix['current']}: expected with "
                f"the wrong code" if fix.get("current") else f"{_plural_genes(len(fixed))} with mismatched codons")
        out.append(f'<tr class="coal-{cls} coal-codefix-row"><td><span class="coal-badge coal-{cls}" '
                   f'aria-hidden="true">{icon}</span> <strong>{word}</strong></td><td colspan="2"><details><summary>'
                   f'{html.escape(what)}</summary><table class="coal-health coal-nested"><thead class="coal-sr"><tr>'
                   f'<th scope="col">Status</th><th scope="col">Sequence</th><th scope="col">What this means</th>'
                   f'</tr></thead><tbody>{"".join(row(x, sequence_sentence(x)) for x in shown)}</tbody></table>'
                   f'{more}</details></td></tr>')
    if rest:
        kinds = {}
        for s in rest:
            kinds[s.get("status")] = kinds.get(s.get("status"), 0) + 1
        what = ", ".join(f"{n:,} {STATUS.get(k, ('', str(k), ''))[1]}" for k, n in kinds.items())
        out.append(f'<tr><td colspan="3">… and {len(rest):,} more sequences ({what}); the report JSON lists '
                   f'every sequence.</td></tr>')
    if oks:  # one row for every OK gene, names listed on demand
        icon, word, cls = STATUS["ok"]
        names = " ".join(f"<code>{html.escape(str(s.get('id', '')))}</code>" for s in oks[:HEALTH_OK_NAMES])
        more = (f" … and {len(oks) - HEALTH_OK_NAMES:,} more (the report JSON lists every sequence)"
                if len(oks) > HEALTH_OK_NAMES else "")
        n = len(oks)
        if n == 1:
            what = "The gene matches" if not allseqs else "1 gene matches"
        else:
            what = f"{'All ' if not allseqs else ''}{n:,} genes match"
        out.append(
            f'<tr class="coal-{cls} coal-oksum"><td><span class="coal-badge coal-{cls}" aria-hidden="true">{icon}'
            f'</span> <strong>{word}</strong></td><td colspan="2"><details><summary>{what} codon for codon'
            f'</summary><p class="coal-oklist">{names}{more}</p></details></td></tr>')
    out.append("</tbody></table>")
    return "".join(out)


def summary_sentence(report: dict) -> str:
    s = report.get("summary", {})
    total = s.get("total", 0)
    parts = [f"{s.get('ok', 0)} OK", f"{s.get('warning', 0)} with warnings", f"{s.get('failed', 0)} failed"]
    if s.get("not_checked"):
        parts.append(f"{s['not_checked']} not checked")
    text = f"{total} sequence{'s' if total != 1 else ''}: " + ", ".join(parts) + "."
    if s.get("columns_removed"):  # say how much -nomismatch removed
        n = s["columns_removed"]
        # columns removed only for codons with unknown or ambiguous bases are named as such
        a = int(s.get("columns_removed_ambiguous") or 0)
        if a >= n:
            why = (f"{'they' if n != 1 else 'it'} had codons with unknown or ambiguous bases (e.g. N) that COAL "
                   f"could not check")
        elif a:
            why = (f"{'they' if n != 1 else 'it'} didn't match ({a} of them only because of codons with unknown or "
                   f"ambiguous bases, e.g. N)")
        else:
            why = f"{'they' if n != 1 else 'it'} didn't match"
        text += f" {n} codon column{'s were' if n != 1 else ' was'} removed because {why} (-nomismatch)."
    return text


def _report_completed(report: dict) -> bool:
    return run_completed(report.get("exit_code"), "",
                         (report.get("options") or {}).get("keep_going")) and not any(
        (e or {}).get("code") == "stopped" for e in report.get("errors") or [])


# The code names as they read inside a sentence: common words lower case, genus names
# capitalised (NCBI: "The Scenedesmus obliquus Mitochondrial Code", "... and the Mycoplasma/Spiroplasma Code").
METHODS_CODE_NAMES: dict[int, str] = {
    1: "standard genetic code",
    2: "vertebrate mitochondrial genetic code",
    3: "yeast mitochondrial genetic code",
    4: "mold, protozoan and coelenterate mitochondrial and Mycoplasma/Spiroplasma genetic code",
    5: "invertebrate mitochondrial genetic code",
    6: "ciliate, dasycladacean and Hexamita nuclear genetic code",
    9: "echinoderm and flatworm mitochondrial genetic code",
    10: "euplotid nuclear genetic code",
    11: "bacterial, archaeal and plant plastid genetic code",
    12: "alternative yeast nuclear genetic code",
    13: "ascidian mitochondrial genetic code",
    14: "alternative flatworm mitochondrial genetic code",
    15: "Blepharisma nuclear genetic code",
    16: "chlorophycean mitochondrial genetic code",
    21: "trematode mitochondrial genetic code",
    22: "Scenedesmus obliquus mitochondrial genetic code",
    23: "Thraustochytrium mitochondrial genetic code",
}


def real_mismatches(seq: dict) -> list[dict]:
    """The gene's mismatch messages with a definite (A/C/G/T) codon."""
    return [m for m in seq.get("messages") or []
            if m.get("code") == "mismatch" and not is_ambiguous_codon(m.get("codon"))]


def gene_kind(seq: dict) -> str:
    """What happened to one gene, by its messages: 'ok', 'mismatch' (definite
    mismatching codons), 'unchecked' (only codons with unknown/ambiguous bases), 'other' (other warnings: non-standard
    residues, pairing by order, ambiguous placement), 'flagged' (written but probably wrong), 'excluded'
    (failed, no row in the output) or 'not_checked'."""
    st = seq.get("status")
    msgs = seq.get("messages") or []
    if st == "failed":
        return "flagged" if any(m.get("code") in SUSPECT_CODES + ("suspect_pairing",) for m in msgs) else "excluded"
    if st == "warning":
        if real_mismatches(seq):
            return "mismatch"
        if any(m.get("code") == "mismatch" for m in msgs):
            return "unchecked"
        return "other"
    return st if st in ("ok", "not_checked") else "other"


def methods_text(report: dict) -> str:
    """A ready-to-paste methods sentence for a run (report["methods_text"] if the engine wrote one)."""
    if report.get("methods_text"):
        return str(report["methods_text"])
    if not _report_completed(report) or not report.get("output"):
        return ""
    opts = report.get("options") or {}
    s = report.get("summary") or {}
    n = int(opts.get("codontable") or 1)
    ver = report.get("coal_version") or "unknown"
    total, ok = s.get("total", 0), s.get("ok", 0)
    # genus names capitalised as NCBI writes them (METHODS_CODE_NAMES)
    phrase = METHODS_CODE_NAMES.get(n, f"NCBI genetic code {n}")
    # a non-numeric version (the demo sample) is shown in brackets
    shown = f"COAL v{ver}" if ver[:1].isdigit() else f"COAL ({ver})"
    parts = [f"Codon alignments were built with {shown} (https://github.com/plumat-coal/COAL), which places "
             f"each coding DNA sequence onto the protein multiple sequence alignment codon by codon (the method "
             f"of PAL2NAL; Suyama et al. 2006), using the {phrase} (NCBI translation table {n})."]
    # counted by message code, not by status: a "warning" gene need not have mismatched codons (U/B/Z-only,
    # pairing-only) and a "failed" gene may still be written (probably wrong data)
    kinds = {}
    for sq in report.get("sequences") or []:
        k = gene_kind(sq)
        kinds.setdefault(k, []).append(sq)
    nomm = opts.get("nomismatch")
    clauses = []

    def were(k):
        return "was" if k == 1 else "were"
    if kinds.get("mismatch"):
        k = len(kinds["mismatch"])
        clauses.append(f"{k} contained codons that did not match the aligned amino acid"
                       + ("" if nomm else " (kept in the alignment)"))
    if kinds.get("unchecked"):
        k = len(kinds["unchecked"])
        clauses.append(f"{k} contained codons with unknown or ambiguous bases (e.g. N) that could not be checked"
                       + ("" if nomm else " (kept in the alignment)"))
    other = kinds.get("other") or []
    if other:
        codes = {m.get("code") for sq in other for m in sq.get("messages") or []}
        why = []
        if "unknown_aa" in codes:
            why.append("contained non-standard amino acids (such as U, B or Z) whose codons were not checked")
        if "suspect_pairing" in codes:
            why.append("were paired with their DNA only by their order in the files")
        if codes - {"unknown_aa", "suspect_pairing"}:
            why.append("had other warnings (see the report)")
        k = len(other)
        lead = f"{k} " + (why[0] if len(why) == 1 else "had warnings: they " + " or ".join(why))
        clauses.append(lead)
    if kinds.get("flagged"):
        k = len(kinds["flagged"])
        clauses.append(f"{k} {were(k)} written to the output but flagged as probably wrong (most likely the wrong "
                       f"DNA, a shifted position or the wrong genetic code); remove "
                       f"{'it' if k == 1 else 'them'} before analysis")
    if kinds.get("excluded"):
        k = len(kinds["excluded"])
        clauses.append(f"{k} could not be aligned and {were(k)} excluded")
    counts = f"Of {total} sequences, {ok} matched their protein exactly"
    if clauses:
        counts += ("; " + "; ".join(clauses[:-1]) + "; and " + clauses[-1] if len(clauses) > 1
                   else ("; " if kinds.get("flagged") else ", and ") + clauses[0])
    parts.append(counts + ".")
    removed = s.get("columns_removed")
    if nomm and (kinds.get("mismatch") or kinds.get("unchecked")):
        parts.append(f"The {removed} codon columns holding such codons were removed (-nomismatch)." if removed
                     else "The alignment columns holding such codons were removed (-nomismatch).")
    if opts.get("nogap"):
        parts.append("Codon columns containing gaps or stop codons were removed.")
    cmd = strip_folders(report.get("command") or "",
                        [item.get("path") for item in list(report.get("inputs") or []) + [report.get("output") or {}]
                         if isinstance(item, dict)])
    if cmd:
        parts.append(f"Command: {cmd}")
    return " ".join(parts)


_PAIRING_TEXT = {
    "id": "Proteins and DNA were paired by matching their IDs exactly.",
    "normalised": "Proteins and DNA were paired by ID after tidying the IDs (for example removing '.1' suffixes).",
    "order": "Proteins and DNA were paired by their order in the two files, because the IDs did not all match.",
}


def explain(report: dict) -> str:
    """Plain-English explanation of a COAL report, focusing on what went wrong and how to fix it."""
    s = report.get("summary", {})
    lines = [summary_sentence(report)]
    if s.get("pairing") in _PAIRING_TEXT:
        lines.append(_PAIRING_TEXT[s["pairing"]])
    code = report.get("exit_code")
    out = report.get("output")
    if out:
        lines.append(f"The codon alignment was written to {out.get('path')} ({out.get('format')} format).")
    elif not _report_completed(report):
        lines.append("No alignment was written, because COAL stopped at an error.")
    for e in report.get("errors") or []:
        lines.append(f"Error: {e.get('text')}")
    failed = [x for x in report.get("sequences", []) if x.get("status") == "failed"]
    warned = [x for x in report.get("sequences", []) if x.get("status") == "warning"]
    if failed:
        lines.append("")
        lines.append("Failed sequences:")
        for x in failed:
            lines.append(f"- {x.get('id')}: {sequence_sentence(x)}")
    if warned:
        lines.append("")
        lines.append("Sequences with warnings (still in the alignment):")
        for x in warned:
            lines.append(f"- {x.get('id')}: {sequence_sentence(x)}")
    hints = report.get("hints") or []
    if hints:
        lines.append("")
        lines.append("Suggestions:")
        for h in hints:
            lines.append(f"- {h.get('text')}")
    if (not failed and not warned and not report.get("errors") and not s.get("columns_removed")
            and not has_problem_hint(hints) and _report_completed(report)):  # never when aborted
        lines.append("Nothing went wrong: every protein matched its DNA.")
    elif failed and not any(h.get("code") == "genetic_code" for h in hints):
        lines.append("")
        lines.append("How to fix a failed sequence: check it is the right DNA for that protein (same gene, "
                     "same species), look for a missing or extra base near where it stops matching, and make "
                     "sure the genetic code is right for where the gene comes from (for example code 2 for "
                     "vertebrate mitochondria).")
    m = methods_text(report)
    if m:
        lines += ["", "Methods text (ready to paste):", m]
    return "\n".join(lines)


def example_paths(name: str) -> tuple[str, str]:
    d = os.path.join(EXAMPLES_DIR, name)
    return os.path.join(d, "protein.aln.fasta"), os.path.join(d, "dna.fasta")


def available_examples() -> list[str]:
    return [n for n in EXAMPLES if all(os.path.isfile(p) for p in example_paths(n))]
