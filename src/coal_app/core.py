"""Gradio-free logic behind the COAL GUI and MCP tools.

Everything here talks to the engine only through ``python -m coal`` in a
subprocess. Newer engine features (``-output clustal|paml|codon``,
``--report FILE.json``, ``--keep-going``) are feature-detected from
``python -m coal --help``; when the engine lacks ``--report`` we build a report
with the same schema (``coal.report/1``) from COAL's stderr instead.
"""
from __future__ import annotations

import hashlib
import html
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE_REPORT = os.path.join(HERE, "sample_report.json")


def _find_examples() -> str:
    """CHANGED(APP-C12): was only <repo>/examples; now COAL_EXAMPLES_DIR, then a repo checkout, then the copy a
    wheel installs under <prefix>/share/coal/examples. See overnight/app/CHANGES.md"""
    import sysconfig
    env = os.environ.get("COAL_EXAMPLES_DIR")
    if env:
        return os.path.normpath(env)
    repo = os.path.normpath(os.path.join(HERE, "..", "..", "examples"))  # <root>/src/coal_app -> <root>/examples
    cands = [repo]
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
    "mito": ("Mitochondrial genes", "Genes that only fit the vertebrate mitochondrial code (2)."),
    "messy_ids": ("Messy IDs", "Protein IDs like AT1G01010.1, DNA IDs like AT1G01010."),
    "one_broken_gene": ("One broken gene", "19 good genes and one with a frameshift."),
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


def prune_runs(max_age_hours: float = RUNS_MAX_AGE_HOURS) -> int:
    """Delete run folders in RUNS_DIR not modified for max_age_hours. Returns how many were removed."""
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
            if os.path.isdir(p) and not os.path.islink(p) and os.path.getmtime(p) < cutoff:
                shutil.rmtree(p, ignore_errors=True)
                removed += not os.path.exists(p)
        except OSError:
            continue
    return removed


def new_run_dir() -> str:
    """A fresh folder for one run's outputs inside RUNS_DIR (old runs are pruned first)."""
    os.makedirs(RUNS_DIR, exist_ok=True)
    prune_runs()
    return tempfile.mkdtemp(prefix="run_", dir=RUNS_DIR)


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


def resolve_input(upload: str | None, path_text: str | None) -> tuple[str | None, str]:
    """(path to use, short description) for one GUI input. A typed/pasted path wins over an upload and is
    read in place (never copied). Raises ValueError with a plain message if the path is unusable."""
    text = (path_text or "").strip().strip('"').strip("'").strip()
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


INTRO_LOCAL = "Nothing leaves your computer: COAL runs here, on this machine."
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
        h = subprocess.run(_coal_cmd() + ["--help"], capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=60).stdout
        v = subprocess.run(_coal_cmd() + ["--version"], capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=60)
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
    """Sequence IDs from a FASTA or CLUSTAL file (first word of each name), in file order."""
    ids: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        first = ""
        for line in fh:
            if line.strip():
                first = line
                break
        if first.startswith(">"):
            ids.append(first[1:].split()[0] if first[1:].split() else "")
            for line in fh:
                if line.startswith(">"):
                    parts = line[1:].split()
                    ids.append(parts[0] if parts else "")
            return ids
        if first.upper().startswith("CLUSTAL"):
            seen = set()
            for line in fh:
                if not line.strip() or line[0].isspace():
                    continue
                name = line.split()[0]
                if name not in seen:
                    seen.add(name)
                    ids.append(name)
    return ids


def strip_version(seq_id: str) -> str:
    """AT1G01010.1 -> AT1G01010 (removes a trailing .N version suffix)."""
    return re.sub(r"\.\d+$", "", seq_id)


# ---------------------------------------------------------------- commands

def _q(p: str) -> str:
    return shlex.quote(p) if os.name != "nt" else (f'"{p}"' if " " in p else p)


def coal_command(protein: str, dna: str, output_format="fasta", codontable=1, nogap=False,
                 nomismatch=False, out="codon_aln.fasta", keep_going=False, report: str | None = None) -> str:
    parts = ["coal", _q(protein), _q(dna)]
    if output_format != "fasta":
        parts += ["-output", output_format]
    if int(codontable) != 1:
        parts += ["-codontable", str(int(codontable))]
    if nogap:
        parts.append("-nogap")
    if nomismatch:
        parts.append("-nomismatch")
    if keep_going:
        parts.append("--keep-going")
    if report:
        parts += ["--report", _q(report)]
    parts += ["-o", _q(out)]
    return " ".join(parts)


def pal2nal_command(protein: str, dna: str, output_format="fasta", codontable=1, nogap=False,
                    nomismatch=False, out="codon_aln.fasta") -> str:
    parts = ["pal2nal.pl", _q(protein), _q(dna), "-output", output_format]
    if int(codontable) != 1:
        parts += ["-codontable", str(int(codontable))]
    if nogap:
        parts.append("-nogap")
    if nomismatch:
        parts.append("-nomismatch")
    return " ".join(parts) + f" > {_q(out)}"


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
        return self.exit_code in (0, 2) and self.output_path is not None


_WARN = re.compile(r"^coal: warning: (\S+) pepAlnPos (\d+): (\S+) does not correspond to (\S+)")
_WARN_X = re.compile(r"^coal: warning: (\S+) pepAlnPos (\d+): (\S+) unknown AA type")
# CHANGED(R2-A2/R2-E2): was (\S+) so "s0:" kept its colon (and the notice was dropped); now the colon is stripped,
# the engine's ENG-c / ENG-d notices are recognised, and a line naming no known sequence becomes a run-level hint.
# See overnight/integration/CHANGES.md
_WARN_ANY = re.compile(r"^coal: warning: (\S+?):? (.+)$")  # any other per-sequence warning
_WARN_NOANCHOR = re.compile(r"^coal: warning: (\S+?): no 10-residue segment of the protein was found in the DNA; "
                            r"codons were assigned by position \((\d+) of (\d+) do not match\)")
_WARN_PAIRING = re.compile(r"^coal: warning: only (\d+) of (\d+) protein IDs are DNA IDs; (.*)$")
PROBLEM_HINTS = ("partial_id_match", "engine_warning", "nomismatch_unknown")  # hints that mean "something is off"


# CHANGED(APP-3): was only "does not correspond to" lines; now every "coal: warning: <id> ..." line. See overnight/app/CHANGES.md
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
        return sid, {"code": "unknown_aa", "column": int(col), "residue": aa,
                     "text": (f"Column {col}: '{aa}' is not one of the 20 standard amino acids; COAL treated it "
                              f"as X (any codon accepted).")}
    m = _WARN_NOANCHOR.match(line)
    if m:
        sid, nbad, nres = m.groups()
        return sid, {"code": "no_anchor",
                     "text": (f"No 10-residue stretch of the protein was found in the DNA, so codons were assigned "
                              f"by position ({nbad} of {nres} do not match). Is this the right DNA for this protein, "
                              f"and the right genetic code?")}
    if _WARN_PAIRING.match(line):
        return None  # run-level; see pairing_warning()
    m = _WARN_ANY.match(line)
    if m:
        sid, rest = m.groups()
        return sid, {"code": "warning", "text": f"COAL warning: {rest.strip()}"}
    return None
_ERR = re.compile(r"^coal: error: (.*)$")


def pairing_warning(line: str) -> dict | None:
    """CHANGED(R2-A2): new. The engine's ENG-d line 'only K of N protein IDs are DNA IDs; ...' as a hint."""
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


def _run_engine(protein, dna, out, fmt, table, nogap, nomismatch, keep_going, report_path, timeout=None):
    f = detect_features()
    timeout = timeout or ENGINE_TIMEOUT
    cmd = _coal_cmd() + [protein, dna, "-o", out, "-codontable", str(int(table))]
    if fmt != "fasta":
        cmd += ["-output", fmt]
    if nogap:
        cmd.append("-nogap")
    if nomismatch:
        cmd.append("-nomismatch")
    if keep_going and f.keep_going:
        cmd.append("--keep-going")
    if report_path and f.report:
        cmd += ["--report", report_path]
    t0 = time.perf_counter()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace",
                           timeout=timeout)
        code, err = p.returncode, p.stderr
    except subprocess.TimeoutExpired:
        code, err = 1, f"coal: error: COAL did not finish within {timeout:g} seconds\n"
    return code, err, time.perf_counter() - t0


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
                          code, stderr, runtime) -> dict:
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
        if w or line.startswith("coal: warning:"):  # CHANGED(R2-A2): was kept under an unknown id and lost
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
    aborted = code not in (0, 2)
    # CHANGED(APP-6): was "order"/"id" even when COAL stopped before pairing; now None (no pairing sentence). See overnight/app/CHANGES.md
    if aborted and general and not failed:
        pairing = None  # COAL stopped before pairing anything (e.g. counts differ, unreadable input)
    keep_going = bool(keep_going and detect_features().keep_going)  # only what was really passed
    seqs = []
    for i, p in enumerate(pids):
        if pairing is None:
            d = None
        elif pairing == "id":
            d = p
        else:
            d = dids[i] if i < len(dids) else None
        if p in failed:
            status, msgs = "failed", failed[p]
        elif aborted:
            status, msgs = "not_checked", []
        elif per.get(p):
            status, msgs = "warning", per[p]
        else:
            status, msgs = "ok", []
        seqs.append({"id": p, "dna_id": d, "pairing": pairing, "status": status, "messages": msgs})
    counts = {s: sum(1 for x in seqs if x["status"] == s) for s in STATUS}
    produced = (not aborted) and os.path.exists(out)
    rep = {
        "schema": "coal.report/1",
        "coal_version": detect_features().version,
        "source": "coal_app (built from COAL's messages; this engine has no --report yet)",
        "command": coal_command(protein, dna, fmt, table, nogap, nomismatch, out, keep_going),
        "pal2nal_equivalent": pal2nal_command(protein, dna, fmt, table, nogap, nomismatch, out),
        "options": {"output_format": fmt, "codontable": int(table), "nogap": bool(nogap),
                    "nomismatch": bool(nomismatch), "keep_going": bool(keep_going)},
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


# CHANGED(APP-7): was the stripped protein ID shown as "DNA" and "probably wrong" even on clean runs; now the
# real DNA ID, and the verdict depends on the run's outcome. See overnight/app/CHANGES.md
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
        return []  # CHANGED(R2-A2): the engine's own partial-match warning already says it (and that it's wrong)
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
                          f"Remove the '.1'-style suffixes so the names match, or use a COAL version that "
                          f"matches IDs after normalising them.")}]
    return [{"code": "id_mismatch",
             "text": (f"{len(missing)} protein ID(s) were not found in the DNA file (for example '{missing[0]}'). "
                      f"{verdict} Check that both files describe the same genes, in the same order.")}]


CANDIDATE_CODES = (2, 3, 4, 5, 9, 11, 13, 14, 21)


def _write_sample(protein, dna, ids, report, work_dir):
    """Protein alignment + DNA FASTA with only the sequences in `ids` (kept in file order)."""
    from .viewer import read_alignment
    want = set(ids)
    prot = [(i, r) for i, r in read_alignment(protein) if i in want]
    dna_recs = read_alignment(dna)
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


def genetic_code_hint(protein, dna, table, report, work_dir) -> dict | None:
    """If sequences have mismatches, see whether another genetic code fits more of them.
    CHANGED(APP-13): was 9 full-input engine runs (600 s limit each, probe files left behind); now as below.
    See overnight/app/CHANGES.md.

    Each candidate code is tried on (at most PROBE_SAMPLE of) the problem sequences only, each probe with a
    short time limit and all probes within PROBE_BUDGET seconds; every probe file is deleted afterwards."""
    import shutil
    seqs = report.get("sequences", [])
    bad_ids = [s["id"] for s in seqs if s["status"] in ("warning", "failed")]
    if not seqs or not bad_ids:
        return None
    # sequences COAL never reached (it stopped at an error) may fit another code too
    bad_ids += [s["id"] for s in seqs if s["status"] == "not_checked"]
    try:
        if os.path.getsize(protein) + os.path.getsize(dna) > 50_000_000:
            return None  # keep the GUI responsive on very big inputs
    except OSError:
        return None
    sample = bad_ids[:PROBE_SAMPLE]
    current_ok = sum(1 for s in seqs if s["status"] == "ok")
    probe_dir = tempfile.mkdtemp(prefix="probe_", dir=work_dir)
    best = None
    try:
        files = _write_sample(protein, dna, sample, report, probe_dir)
        if files is None:
            return None
        ps, ds = files
        t_end = time.perf_counter() + PROBE_BUDGET
        for c in CANDIDATE_CODES:
            left = t_end - time.perf_counter()
            if c == int(table) or left <= 1:
                continue
            out = os.path.join(probe_dir, f"probe_{c}.fasta")
            code, err, _ = _run_engine(ps, ds, out, "fasta", c, False, False, True, None,
                                       timeout=min(PROBE_TIMEOUT, left))
            if os.path.exists(out):
                os.remove(out)
            if code not in (0, 2):
                continue
            warned = {w[0] for w in map(_parse_warning, err.splitlines()) if w}
            fixed = sum(1 for i in sample if i not in warned)
            if best is None or fixed > best[1]:
                best = (c, fixed)
    finally:
        shutil.rmtree(probe_dir, ignore_errors=True)
    if not best or best[1] == 0:
        return None
    c, fixed = best
    if len(sample) == len(bad_ids):
        n_ok = current_ok + fixed
        most = "All" if n_ok == len(seqs) else ("Most" if n_ok * 2 > len(seqs) else f"{n_ok} of {len(seqs)}")
        text = (f"{most} sequences fit genetic code {c} ({GENETIC_CODES[c]}) without any mismatches, "
                f"compared with {current_ok} under code {int(table)}. Re-run with code {c}?")
    else:
        text = (f"{fixed} of {len(sample)} sequences checked (a sample of the {len(bad_ids)} with problems or "
                f"not checked) fit genetic code {c} "
                f"({GENETIC_CODES[c]}) without any mismatches; {current_ok} of {len(seqs)} sequences fit code "
                f"{int(table)}. Re-run with code {c}?")
    return {"code": "genetic_code", "text": text, "suggested_codontable": c}


def mismatch_probe(protein, dna, table, keep_going, report, work_dir) -> dict:
    """CHANGED(R2-A1/R2-D5/R2-E1): new. With -nomismatch the engine (like PAL2NAL) prints no mismatch
    warnings, so the report said every gene matched. This re-runs the same inputs without -nomismatch (output
    to a temporary file, deleted) to find which genes had mismatching codons, marks them "warning" (columns
    removed) and stores the number of removed codon columns in summary["columns_removed"].
    See overnight/integration/CHANGES.md"""
    import shutil
    probe_dir = tempfile.mkdtemp(prefix="probe_nomismatch_", dir=work_dir)
    try:
        code, err, _ = _run_engine(protein, dna, os.path.join(probe_dir, "probe.fasta"), "fasta", table,
                                   False, False, keep_going, None)
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
    for line in err.splitlines():
        w = _parse_warning(line)
        if w and w[1].get("code") == "mismatch":
            m = dict(w[1], removed=True)
            per.setdefault(w[0], []).append(m)
            cols.add(m["column"])
    for sq in report.get("sequences", []):
        mm = per.get(sq["id"])
        if mm and sq["status"] in ("ok", "warning"):
            sq["status"] = "warning"
            sq["messages"] = mm + list(sq.get("messages") or [])
    summ = report["summary"]
    for k in ("ok", "warning"):
        summ[k] = sum(1 for x in report.get("sequences", []) if x["status"] == k)
    summ["columns_removed"] = len(cols)
    return report


def has_problem_hint(hints) -> bool:
    return any(isinstance(h, dict) and h.get("code") in PROBLEM_HINTS for h in hints or [])


def banner_head(ok: bool, summary: dict, hints=()) -> str:
    """CHANGED(R2-A1/R2-A2): was inline in ui.run_ui; "every sequence matched" now also needs no removed columns
    and no run-level problem (e.g. a partial ID match)."""
    if (ok and not summary.get("failed") and not summary.get("warning") and not summary.get("columns_removed")
            and not has_problem_hint(hints)):
        return "### ✔ Done: every sequence matched"
    if ok:
        return "### ⚠ Done, with some problems (see the table)"
    return "### ✖ COAL stopped: no alignment was written"


def run_coal(protein_path: str, dna_path: str, output_format: str = "fasta", codontable: int = 1,
             nogap: bool = False, nomismatch: bool = False, keep_going: bool = True,
             out_dir: str | None = None, probe_codes: bool = True) -> RunResult:
    """Run COAL once and return the output, a coal.report/1 report and hints."""
    f = detect_features()
    if output_format not in f.output_formats:
        raise ValueError(f"output format '{output_format}' is not supported by this COAL "
                         f"(supported: {', '.join(f.output_formats)})")
    if int(codontable) not in GENETIC_CODES:
        raise ValueError(f"unknown genetic code {codontable}")
    out_dir = out_dir or new_run_dir()  # CHANGED(APP-12): was mkdtemp("coal_") in /tmp, never removed. See overnight/app/CHANGES.md
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(protein_path))[0].replace(".aln", "")
    out = os.path.join(out_dir, f"{stem}.codon.{_EXT.get(output_format, 'txt')}")
    report_path = os.path.join(out_dir, f"{stem}.report.json")
    if os.path.exists(report_path):
        os.remove(report_path)
    code, err, runtime = _run_engine(protein_path, dna_path, out, output_format, codontable,
                                     nogap, nomismatch, keep_going, report_path)
    report = None
    if f.report and os.path.exists(report_path):
        try:
            with open(report_path, encoding="utf-8") as fh:
                report = json.load(fh)
        except (OSError, ValueError):
            report = None
    if report is None:
        report = build_fallback_report(protein_path, dna_path, out, output_format, codontable, nogap,
                                       nomismatch, keep_going, code, err, runtime)
        if nomismatch and code == 0:
            mismatch_probe(protein_path, dna_path, codontable, keep_going, report, out_dir)
        from .viewer import locate_breaks
        locate_breaks(protein_path, dna_path, report)
        hints = list(report.get("hints") or []) + id_hints(protein_path, dna_path, report)
        if probe_codes:
            h = genetic_code_hint(protein_path, dna_path, codontable, report, out_dir)
            if h:
                hints.append(h)
        report["hints"] = hints
        with open(report_path, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2)
    produced = code in (0, 2) and os.path.exists(out)
    return RunResult(
        exit_code=code,
        output_path=out if produced else None,
        report_path=report_path,
        report=report,
        stderr=err,
        coal_cmd=coal_command(os.path.basename(protein_path), os.path.basename(dna_path), output_format,
                              codontable, nogap, nomismatch, os.path.basename(out),
                              keep_going and f.keep_going),
        pal2nal_cmd=pal2nal_command(os.path.basename(protein_path), os.path.basename(dna_path),
                                    output_format, codontable, nogap, nomismatch, os.path.basename(out)),
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
    # CHANGED(APP-14): was plain utf-8 and no type checks (AttributeError on {"summary": []}); now utf-8-sig
    # (a BOM is fine) and every field explain() uses is type-checked. See overnight/app/CHANGES.md
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
    """GUI banner for the sample report. CHANGED(APP-15): was "what a run looks like with a newer COAL";
    now a clearly labelled demo sample (a real run of examples/mito with the wrong genetic code)."""
    return ("### Demo sample (not your data)\nThe 24 mitochondrial genes of the *Mitochondrial genes* example, run "
            "with genetic code 1, the wrong code for them: see the tip to re-run with code 2. "
            + summary_sentence(rep))


def _mismatch_sentence(msgs: list[dict]) -> str:
    mm = [m for m in msgs if m.get("code") == "mismatch"]
    others = " ".join(m.get("text", "") for m in msgs if m.get("code") != "mismatch")
    if not mm:
        return others or "COAL printed a warning for this sequence."
    return _mm_sentence(mm) + (" " + others if others else "")


def _mm_sentence(mm: list[dict]) -> str:
    first = mm[0]
    n = len(mm)
    where = f" (first at column {first['column']})" if first.get("column") else ""
    detail = ""
    if first.get("residue") and first.get("codon"):
        detail = f": the protein says {first['residue']} but the DNA has {first['codon']}"
    # CHANGED(R2-A1): mismatches found by the -nomismatch probe say their columns were removed
    fate = ("Those columns were removed from the alignment (-nomismatch)" if first.get("removed")
            else "Kept in the alignment")
    return (f"{n} codon{'s' if n != 1 else ''} {'do' if n != 1 else 'does'} not match "
            f"{'their' if n != 1 else 'its'} amino acid{where}{detail}. {fate}; check the genetic code.")


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
        return "Not checked: COAL stopped at an earlier error."
    if st == "warning":
        return _mismatch_sentence(msgs)
    if st == "failed":
        for m in msgs:
            if isinstance(m.get("column"), int) and m.get("code") in ("frameshift", "inconsistent", "missing_dna"):
                tail = {"frameshift": " Check the DNA near there (a sequencing or assembly error?).",
                        "inconsistent": " Is this the right DNA for this protein, and the right genetic code?",
                        "missing_dna": " The DNA looks truncated."}[m["code"]]
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


def _jump_links(seq: dict, limit: int = 12) -> str:
    from .viewer import cell_id
    # CHANGED(APP-3/APP-1): was every message column linked; now only cells the viewer anchors. See overnight/app/CHANGES.md
    # only codes the viewer gives an anchor to (it outlines mismatched and frameshift cells)
    cols = sorted({m["column"] for m in seq.get("messages") or []
                   if isinstance(m.get("column"), int) and m.get("code") in ("mismatch", "frameshift")})
    if not cols:
        return ""
    links = ", ".join(f'<a href="#{cell_id(str(seq.get("id", "")), c)}">column {c}</a>' for c in cols[:limit])
    more = f" and {len(cols) - limit} more" if len(cols) > limit else ""
    return f'<br><span class="coal-small">Show in alignment: {links}{more}</span>'


def health_table_html(report: dict, links: bool = False) -> str:
    """Accessible HTML table (icon + word + sentence, never colour alone).
    With links=True, problem columns link to the alignment viewer's cells (same page). No links when
    -nogap/-nomismatch removed columns: report columns then don't exist in the output."""
    opts = report.get("options") if isinstance(report.get("options"), dict) else {}
    if opts.get("nogap") or opts.get("nomismatch"):  # CHANGED(APP-1): no jump links in codon-only mode
        links = False
    order = {"failed": 0, "warning": 1, "not_checked": 2, "ok": 3}
    seqs = sorted(report.get("sequences", []), key=lambda s: order.get(s.get("status"), 4))
    out = ['<table class="coal-health"><caption>Sequence health, problems first</caption>',
           '<thead><tr><th scope="col">Status</th><th scope="col">Sequence</th>'
           '<th scope="col">What this means</th></tr></thead><tbody>']
    for s in seqs:
        icon, word, cls = STATUS.get(s.get("status"), ("?", str(s.get("status")), ""))
        pair = ""
        if s.get("dna_id") and s.get("dna_id") != s.get("id"):
            pair = f'<br><span class="coal-small">paired with DNA {html.escape(str(s["dna_id"]))}</span>'
        out.append(
            f'<tr class="coal-{cls}"><td><span class="coal-badge coal-{cls}" aria-hidden="true">{icon}</span> '
            f'<strong>{word}</strong></td><td><code>{html.escape(str(s.get("id", "")))}</code>{pair}</td>'
            f'<td>{html.escape(sequence_sentence(s))}{_jump_links(s) if links else ""}</td></tr>')
    out.append("</tbody></table>")
    return "".join(out)


def summary_sentence(report: dict) -> str:
    s = report.get("summary", {})
    total = s.get("total", 0)
    parts = [f"{s.get('ok', 0)} OK", f"{s.get('warning', 0)} with warnings", f"{s.get('failed', 0)} failed"]
    if s.get("not_checked"):
        parts.append(f"{s['not_checked']} not checked")
    text = f"{total} sequence{'s' if total != 1 else ''}: " + ", ".join(parts) + "."
    if s.get("columns_removed"):  # CHANGED(R2-A1): say how much -nomismatch removed
        n = s["columns_removed"]
        text += (f" {n} codon column{'s were' if n != 1 else ' was'} removed because "
                 f"{'they' if n != 1 else 'it'} didn't match (-nomismatch).")
    return text


def methods_text(report: dict) -> str:
    """A ready-to-paste methods sentence for a run (report["methods_text"] if the engine wrote one)."""
    if report.get("methods_text"):
        return str(report["methods_text"])
    if report.get("exit_code") not in (0, 2) or not report.get("output"):
        return ""
    opts = report.get("options") or {}
    s = report.get("summary") or {}
    n = int(opts.get("codontable") or 1)
    name = GENETIC_CODES.get(n, f"NCBI code {n}")
    ver = report.get("coal_version") or "unknown"
    total, ok, warn, failed = (s.get(k, 0) for k in ("total", "ok", "warning", "failed"))
    # CHANGED(APP-15): was always "COAL v<version>"; a non-numeric version (the demo sample) is shown in brackets
    shown = f"COAL v{ver}" if ver[:1].isdigit() else f"COAL ({ver})"
    parts = [f"Codon alignments were built with {shown} (https://github.com/chloe-newt/COAL), which places "
             f"each coding DNA sequence onto the protein multiple sequence alignment codon by codon (the method "
             f"of PAL2NAL; Suyama et al. 2006), using the {name.split(' (')[0].lower()} genetic code "
             f"(NCBI translation table {n})."]
    counts = f"Of {total} sequences, {ok} matched their protein exactly"
    if warn:
        removed = s.get("columns_removed")
        counts += (f", {warn} contained codons that did not match the aligned amino acid"
                   + ((f" (those {removed} codon columns were removed)" if removed else
                       " (those alignment columns were removed)") if opts.get("nomismatch") else
                      " (kept in the alignment)"))
    if failed:
        counts += f", and {failed} could not be aligned and {'was' if failed == 1 else 'were'} excluded"
    parts.append(counts + ".")
    if opts.get("nogap"):
        parts.append("Codon columns containing gaps or stop codons were removed.")
    cmd = report.get("command") or ""
    for item in list(report.get("inputs") or []) + [report.get("output") or {}]:
        pth = item.get("path") if isinstance(item, dict) else None
        if pth:
            cmd = cmd.replace(pth, os.path.basename(pth))
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
    elif code not in (0, 2):
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
            and not has_problem_hint(hints)):
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
