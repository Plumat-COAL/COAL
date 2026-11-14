"""The MCP tools AI agents use (codon_align, explain_report, dnds_pairs, pipeline_pairs). Plain functions: served by Gradio's MCP server
(``python -m coal_app``) or by the stand-alone stdio server (``python -m coal_app --mcp``).

The result says what happened truthfully (``completed`` = COAL finished; ``all_ok`` /
``needs_attention`` = whether every sequence matched), errors raise :class:`CoalToolError` (MCP clients see
``isError: true``), the arguments carry descriptions and enums, and stderr is summarised instead of copied.
"""
from __future__ import annotations

import os
from typing import Annotated, Literal

from . import core, dnds

# pydantic / typing_extensions come with both extras (Gradio and mcp); the core install (no extras) has neither,
# and only needs the plain functions, so the schema metadata degrades to nothing there.
try:
    from pydantic import Field
except ImportError:  # pragma: no cover - core install without extras
    def Field(**_kw):  # noqa: N802 - stands in for pydantic.Field
        return None
try:
    from typing_extensions import TypedDict  # pydantic needs this one on Python < 3.12
except ImportError:  # pragma: no cover
    from typing import TypedDict

try:  # mcp >= 2 (the `mcp` extra): a ToolError is an anticipated failure, its message reaches the model
    from mcp.server.mcpserver.exceptions import ToolError as _ToolErrorBase
except ImportError:
    try:  # mcp 1.x (Gradio's, and the Windows app's stdio server) has ToolError under fastmcp
        from mcp.server.fastmcp.exceptions import ToolError as _ToolErrorBase
    except ImportError:  # no MCP at all: any exception becomes isError
        _ToolErrorBase = None

# Without an MCP ToolError base, the bases are just (ValueError,): (Exception, ValueError) would be an invalid MRO
# (Exception before its subclass ValueError) and importing this module would fail.
_TOOL_ERROR_BASES = (_ToolErrorBase, ValueError) if _ToolErrorBase is not None else (ValueError,)


class CoalToolError(*_TOOL_ERROR_BASES):
    """A failure the caller can fix (missing file, unknown genetic code, unsupported format, path not allowed).
    Raised instead of returning {"ok": false}, so MCP clients get a real tool error (isError: true)."""


# The NCBI genetic codes COAL supports (same list as core.GENETIC_CODES; a test keeps them equal).
GeneticCode = Literal[1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23]
OutputFormat = Literal["fasta", "clustal", "paml", "codon"]

MAX_PROBLEMS = 10       # problems listed in full in the result; the rest are in the report
MAX_ENGINE_LINES = 5    # non-warning engine messages (errors, notices) kept in the result


class Counts(TypedDict):
    total: int
    ok: int
    warning: int
    failed: int
    not_checked: int


class Problem(TypedDict):
    id: str | None
    status: str | None
    explanation: str


class CodonAlignResult(TypedDict):
    """What codon_align returns (also published as the MCP outputSchema)."""
    completed: bool
    all_ok: bool
    needs_attention: bool
    counts: Counts
    suggested_codontable: int | None
    next_step: str
    output_path: str | None
    report_path: str
    exit_code: int
    summary: dict
    summary_text: str
    hints: list[dict]
    problems: list[Problem]
    problems_total: int
    engine_messages: list[str]
    warning_lines: int
    methods_text: str
    coal_command: str
    pal2nal_command: str
    pal2nal_caveats: list[str]


def _counts(summary: dict) -> Counts:
    return Counts(**{k: int(summary.get(k) or 0) for k in ("total", "ok", "warning", "failed", "not_checked")})


def _suggested_code(hints: list[dict]) -> int | None:
    for h in hints or []:
        if h.get("code") == "genetic_code" and h.get("suggested_codontable") is not None:
            return int(h["suggested_codontable"])
    return None


def _engine_messages(stderr: str) -> tuple[list[str], int]:
    """First few non-warning lines (errors, notices) and the number of per-codon warning lines."""
    lines = [ln.strip() for ln in (stderr or "").splitlines() if ln.strip()]
    warn = [ln for ln in lines if ln.startswith("coal: warning:")]
    other = [ln for ln in lines if not ln.startswith("coal: warning:")]
    return other[:MAX_ENGINE_LINES], len(warn)


def _next_step(completed: bool, all_ok: bool, counts: Counts, code: int | None) -> str:
    if code is not None:
        return f"Re-run codon_align with codontable={code} (see hints), then check again."
    if all_ok:
        return "All sequences match their protein; the alignment at output_path is ready to use."
    if not completed:
        return "COAL stopped before writing an alignment: call explain_report(report_path) for the cause and fix."
    return (f"{counts['warning']} sequence(s) with warnings and {counts['failed']} failed: call "
            "explain_report(report_path) and tell the user which genes need checking before using the alignment.")


def summarise(r: core.RunResult) -> CodonAlignResult:
    """Turn a core.RunResult into the compact, truthful tool result (pure function, tested without MCP)."""
    summary = r.report.get("summary", {}) or {}
    counts = _counts(summary)
    completed = bool(r.ok)  # core.RunResult.ok: the process finished and wrote an alignment
    clean = counts["total"] > 0 and counts["ok"] == counts["total"]
    risky_hint = any(h.get("code") in ("genetic_code", "partial_id_match") for h in r.hints or [])
    all_ok = completed and clean and not risky_hint
    code = _suggested_code(r.hints)
    seqs = [s for s in r.report.get("sequences", []) if s.get("status") in ("warning", "failed")]
    problems = [Problem(id=s.get("id"), status=s.get("status"), explanation=core.sequence_sentence(s))
                for s in seqs[:MAX_PROBLEMS]]
    engine_lines, n_warn = _engine_messages(r.stderr)
    return CodonAlignResult(
        completed=completed,
        all_ok=all_ok,
        needs_attention=not all_ok,
        counts=counts,
        suggested_codontable=code,
        next_step=_next_step(completed, all_ok, counts, code),
        output_path=r.output_path,
        report_path=r.report_path,
        exit_code=r.exit_code,
        summary=summary,
        summary_text=core.summary_sentence(r.report),
        hints=list(r.hints or []),
        problems=problems,
        problems_total=len(seqs),
        engine_messages=engine_lines,
        warning_lines=n_warn,
        methods_text=core.methods_text(r.report),
        coal_command=r.coal_cmd,
        pal2nal_command=r.pal2nal_cmd,
        pal2nal_caveats=core.pal2nal_caveats(r.report.get("options") or {}),
    )


def _check_input(path: str) -> None:
    # only sandbox paths are allowed when hosted
    why = core.check_path_allowed(path)
    if why:
        raise CoalToolError(why)
    if not os.path.isfile(path):
        raise CoalToolError(f"file not found: {path} (give an absolute path; relative paths resolve against "
                            f"the server's working directory {os.getcwd()})")


# This tool runs COAL with --keep-going (core.run_coal default), so with strict refused genes are left out and the
# run goes on.
def codon_align(
    protein_alignment_path: Annotated[str, Field(
        description="Absolute path to the protein multiple sequence alignment (FASTA or CLUSTAL).")],
    dna_fasta_path: Annotated[str, Field(
        description="Absolute path to the unaligned coding DNA (FASTA), one sequence per protein, same IDs.")],
    output_format: Annotated[OutputFormat, Field(
        description="Output format of the codon alignment: fasta (default), clustal, paml (codeml) or codon.")]
    = "fasta",
    codontable: Annotated[GeneticCode, Field(
        description="NCBI genetic code: 1 standard (nuclear), 2 vertebrate mitochondrial, 3 yeast mito, "
                    "4 mold/protozoan mito and Mycoplasma, 5 invertebrate mito, 9 echinoderm/flatworm mito, "
                    "11 bacteria, archaea and plastids; also 6, 10, 12-16, 21-23.")] = 1,
    nogap: Annotated[bool, Field(description="Remove codon columns that contain a gap or a stop codon.")] = False,
    nomismatch: Annotated[bool, Field(
        description="Remove columns where a codon doesn't match its amino acid.")] = False,
    strict: Annotated[bool, Field(
        description="Refuse probably-wrong genes (left out and reported as failed; the others are written).")]
    = False,
    blockonly: Annotated[bool, Field(
        description="Keep only the columns marked '#' under the blocks of a CLUSTAL protein alignment.")] = False,
    save_as: Annotated[str, Field(
        description="Optional absolute path to save the codon alignment to, e.g. next to the input files "
                    "(C:\\data\\genes.codon.fasta). An existing folder means: save there under COAL's file name. "
                    "Never overwrites an existing file. Default: COAL's own run folder (output_path).")] = "",
) -> dict:
    """Build a codon alignment from a protein alignment and the matching unaligned coding DNA (PAL2NAL replacement).

    Use this when you have aligned protein sequences and the DNA (CDS) of the same genes and need
    the DNA aligned codon by codon, for example before dN/dS analysis (PAML codeml, HyPhy).
    Sequences are paired by ID. Codons that don't match their amino acid are reported as warnings.
    Read all_ok / needs_attention, not completed, to know whether the alignment can be used.

    Args:
        protein_alignment_path: Absolute path to the protein multiple sequence alignment (FASTA or CLUSTAL).
        dna_fasta_path: Absolute path to the unaligned coding DNA sequences (FASTA), one per protein.
        output_format: "fasta" (default), "clustal", "paml" or "codon".
        codontable: NCBI genetic code number: 1 = standard (nuclear), 2 = vertebrate mitochondrial,
            5 = invertebrate mitochondrial, 11 = bacteria/plastids. Also 3, 4, 6, 9, 10, 12-16, 21-23.
        nogap: If true, remove codon columns that contain a gap or a stop codon.
        nomismatch: If true, remove columns where a codon doesn't match its amino acid.
        strict: If true, refuse probably-wrong genes (protein and DNA probably unrelated, or position of the gene a
            guess): they are left out of the alignment and reported as failed, and the other genes are still
            written; only a partial ID match (sequences that would be paired by order) stops the run. By default
            such genes are written, as PAL2NAL does, and reported as failed.
        blockonly: If true, keep only the columns marked '#' under the blocks of a CLUSTAL protein alignment.
        save_as: Optional absolute path (or existing folder) to save the codon alignment to, for example next to
            the input files, so the person finds it. Never overwrites a file. Empty: COAL's own run folder.

    Returns:
        A dict with: completed (COAL finished and wrote an alignment), all_ok (every sequence matched its
        protein and no hint says the data looks wrong), needs_attention (the opposite of all_ok), counts
        (total / ok / warning / failed / not_checked sequences), suggested_codontable (a better genetic code,
        or null), next_step (what to do now), output_path (the codon alignment file, or null), report_path
        (a coal.report/1 JSON file, pass it to explain_report), summary and summary_text, hints (suggestions,
        e.g. {"code": "genetic_code", "suggested_codontable": 2, "text": ...}), problems (the first few
        {id, status, explanation} with a warning or failure; problems_total counts them all), engine_messages
        (COAL's first errors and notices) and warning_lines (how many per-codon warnings it printed),
        methods_text (a ready-to-paste methods sentence), coal_command, pal2nal_command (the closest PAL2NAL v14
        command line) and pal2nal_caveats (when that command does not do the same, e.g. -nogap with another
        genetic code). A missing file, an unknown genetic code or an unsupported format is a tool error.
    """
    for p in (protein_alignment_path, dna_fasta_path):
        _check_input(p)
    save_as = (save_as or "").strip()
    if save_as:
        # lets agents (Claude Desktop) write the result where they choose, not only into COAL's temporary run folder
        if core.hosted():
            raise CoalToolError("save_as is not available on this shared server; use output_path")
        if not os.path.isabs(save_as):
            raise CoalToolError(f"save_as must be an absolute path, not {save_as!r}")
    if int(codontable) not in core.GENETIC_CODES:
        raise CoalToolError(f"unknown genetic code {codontable} (supported: "
                            f"{', '.join(str(c) for c in core.GENETIC_CODES)})")
    try:
        # a run folder in core.RUNS_DIR, pruned over time
        r = core.run_coal(protein_alignment_path, dna_fasta_path, output_format, int(codontable),
                          bool(nogap), bool(nomismatch), out_dir=core.new_run_dir(), strict=bool(strict),
                          blockonly=bool(blockonly))
    except ValueError as e:
        raise CoalToolError(str(e)) from e
    result = dict(summarise(r))
    if save_as and result.get("output_path"):
        dest = save_as
        if os.path.isdir(dest):
            dest = os.path.join(dest, os.path.basename(result["output_path"]))
        if os.path.exists(dest):
            raise CoalToolError(f"not saved: {dest} already exists (COAL never overwrites a file; choose another "
                                f"name). The alignment is at {result['output_path']}")
        parent = os.path.dirname(dest)
        if not os.path.isdir(parent):
            raise CoalToolError(f"not saved: the folder {parent} does not exist. The alignment is at "
                                f"{result['output_path']}")
        import shutil
        try:
            shutil.copyfile(result["output_path"], dest)
        except OSError as e:
            raise CoalToolError(f"not saved to {dest}: {e.strerror or e}. The alignment is at "
                                f"{result['output_path']}") from e
        result["output_path"] = dest
    return result


def explain_report(
    report_path: Annotated[str, Field(
        description="Absolute path to a coal.report/1 JSON file (report_path from codon_align, or COAL --report).")],
) -> str:
    """Explain a COAL run report (coal.report/1 JSON) in plain English: what failed, why, and how to fix it.

    Args:
        report_path: Absolute path to the report JSON written by COAL (--report) or returned by codon_align.

    Returns:
        A short plain-English explanation: counts of OK / warning / failed sequences, how proteins and
        DNA were paired, one line per problem sequence with its likely cause (wrong genetic code,
        mismatched IDs, frameshift), suggested fixes, and a ready-to-paste methods sentence.
        A report that cannot be read is a tool error.
    """
    why = core.check_path_allowed(report_path)  # sandboxed when hosted
    if why:
        raise CoalToolError(f"Could not read the report: {why}")
    try:
        rep = core.load_report(report_path)
    except (OSError, ValueError) as e:
        raise CoalToolError(f"Could not read the report: {e}") from e
    return core.explain(rep)


# ---------------------------------------------------------------- Ka/Ks and the dN/dS pipeline (dnds.py)

KaksMethod = Literal["YN", "NG", "LWL", "MLWL", "LPB", "MLPB", "MYN", "Comeron", "Li"]
MAX_ROWS = 20       # result rows returned in full; the table file has them all
MAX_FLAGGED = 10    # rows with a status other than ok, listed with their note


class StatusCounts(TypedDict):
    total: int
    ok: int
    warning: int
    not_converged: int
    saturated: int
    no_data: int
    error: int


class StatusMeaning(TypedDict):
    status: str
    count: int
    meaning: str


class DndsPairsResult(TypedDict):
    """What dnds_pairs returns (also published as the MCP outputSchema)."""
    completed: bool
    all_estimated: bool
    exit_code: int
    error: str
    counts: StatusCounts
    pairs: int
    methods: list[str]
    input_format: str
    output_path: str | None
    rows: list[dict]
    rows_total: int
    flagged: list[dict]
    flagged_total: int
    statuses: list[StatusMeaning]
    engine_messages: list[str]
    cite: list[str]
    references: str
    next_step: str
    coal_command: str
    limits: dict


class PipelinePairsResult(TypedDict):
    """What pipeline_pairs returns (also published as the MCP outputSchema)."""
    completed: bool
    status: str
    exit_code: int
    error: str
    counts: dict
    pair_counts: StatusCounts
    rows: list[dict]
    rows_total: int
    flagged: list[dict]
    flagged_total: int
    failures_total: int
    failures: list[dict]
    statuses: list[StatusMeaning]
    out_dir: str
    outputs: dict
    report_path: str | None
    zip_path: str | None
    diamond: str | None
    cite: list[str]
    references: str
    next_step: str
    coal_command: str
    progress: list[str]
    limits: dict


_ROW_KEYS = ("seq1", "seq2", "method", "dN", "dS", "dNdS", "SE_dN", "SE_dS", "codons_compared", "prot_identity",
             "status", "note")


def _row(r: dict) -> dict:
    return {k: r[k] for k in _ROW_KEYS if k in r}


def _flagged(rows: list[dict]) -> list[dict]:
    return [_row(r) for r in rows if (r.get("status") or "") != "ok"]


def _check_methods(methods, codontable) -> list[str]:
    try:
        return dnds.validate_methods(list(methods or []), int(codontable))
    except ValueError as e:
        raise CoalToolError(str(e)) from e


def _check_number(value, name: str) -> float:
    try:
        v = float(value or 0)
    except (TypeError, ValueError):
        raise CoalToolError(f"{name} must be a number, not {value!r}") from None
    if v < 0:
        raise CoalToolError(f"{name} cannot be negative")
    return v


def _kaks_next_step(r: dnds.KaksResult) -> str:
    c = r.counts
    if not r.ok:
        return f"Ka/Ks did not run: {r.error or 'see engine_messages'}. Fix the input and call dnds_pairs again."
    if r.exit_code == 2 or not (c.get("ok") or c.get("warning")):
        return ("No pair has an estimate: the input is probably not a codon alignment (protein sequences, or rows "
                "of different lengths). Make one with codon_align first, then call dnds_pairs on its output_path.")
    if c.get("ok") == c.get("total"):
        return "Every result is estimated; the table is at output_path. Cite the references when you publish."
    return ("Some results have no estimate or need care: read flagged (status and note) and statuses, and tell the "
            "person which pairs are affected before using the numbers.")


def dnds_pairs(
    codon_alignment_path: Annotated[str, Field(
        description="Absolute path to aligned coding DNA: an aligned codon FASTA (e.g. codon_align's output_path) "
                    "or a KaKs_Calculator AXT file. Up to 100 MB and 500,000 results (pairs x methods).")],
    pairs: Annotated[Literal["all", "consecutive"], Field(
        description="FASTA input: 'all' compares every pair of sequences (n*(n-1)/2 pairs), 'consecutive' "
                    "compares records 1+2, 3+4, ... AXT files are read pair by pair.")] = "all",
    methods: Annotated[list[KaksMethod], Field(
        description="Ka/Ks methods, one result row per pair and method: YN (Yang & Nielsen 2000, the default), "
                    "NG, LWL, MLWL, LPB, MLPB, MYN, Comeron, Li (Comeron and Li: genetic code 1 only).")]
    = ["YN"],  # noqa: B006 - never mutated
    codontable: Annotated[GeneticCode, Field(
        description="NCBI genetic code: 1 standard, 2 vertebrate mitochondrial, 11 bacteria and plastids, ...")]
    = 1,
    min_codons: Annotated[int, Field(
        description="Flag (status warning, nothing dropped) estimates from fewer compared codons than this. "
                    "0: off.", ge=0)] = 0,
    max_ds: Annotated[float, Field(
        description="Flag (status warning, nothing dropped) estimates with dS above this (saturation). 0: off.",
        ge=0)] = 0,
    save_as: Annotated[str, Field(
        description="Optional absolute path (or existing folder) to save the table (TSV) to. Never overwrites a "
                    "file. Default: COAL's own run folder (output_path).")] = "",
) -> dict:
    """Estimate Ka (dN), Ks (dS) and Ka/Ks (dN/dS) for pairs of aligned coding sequences (COAL's `coal dnds`).

    Use this when you already have a codon alignment (aligned coding DNA, for example the output_path of
    codon_align) or an AXT file, and need dN/dS per pair. For two files of unaligned coding sequences of two
    species use pipeline_pairs instead; for a protein alignment plus DNA call codon_align first.
    Every pair gets one row per method with a status: ok, warning (estimated; read the note), not_converged,
    saturated, no_data or error (no estimate, or only part of one); one bad pair never stops the run.

    Args:
        codon_alignment_path: Absolute path to an aligned codon FASTA or an AXT file (at most 100 MB).
        pairs: "all" (every pair of a FASTA) or "consecutive" (records 1+2, 3+4, ...). Ignored for AXT.
        methods: Ka/Ks methods, default ["YN"]. Also NG, LWL, MLWL, LPB, MLPB, MYN, Comeron, Li.
        codontable: NCBI genetic code number (1 standard, 2 vertebrate mitochondrial, 11 bacteria/plastids, ...).
        min_codons: Flag estimates from fewer compared codons than this (status warning). 0: off.
        max_ds: Flag estimates with dS above this (status warning). 0: off.
        save_as: Optional absolute path (or existing folder) for the table; never overwrites a file.

    Returns:
        A dict with: completed (the table was written), all_estimated (every result has status ok), exit_code,
        error (why nothing could be estimated, or ""), counts (results per status), pairs, methods,
        input_format (fasta or axt), output_path (the TSV: seq1 seq2 method dN dS dNdS SE_dN SE_dS N_sites
        S_sites Nd Sd codons_compared status note), rows (the first 20 rows; rows_total counts them all),
        flagged (the first rows whose status is not ok, with their note; flagged_total), statuses (what each
        status that occurs means), engine_messages, cite and references (what to cite), next_step, coal_command
        (the same run on the command line) and limits. Too large an input, an unknown method, Comeron/Li with a
        code other than 1, or an install without COAL's C++ engine is a tool error.
    """
    _check_input(codon_alignment_path)
    save_as = (save_as or "").strip()
    if save_as:
        if core.hosted():
            raise CoalToolError("save_as is not available on this shared server; use output_path")
        if not os.path.isabs(save_as):
            raise CoalToolError(f"save_as must be an absolute path, not {save_as!r}")
    if int(codontable) not in core.GENETIC_CODES:
        raise CoalToolError(f"unknown genetic code {codontable} (supported: "
                            f"{', '.join(str(c) for c in core.GENETIC_CODES)})")
    chosen = _check_methods(methods, codontable)
    min_c = int(_check_number(min_codons, "min_codons"))
    mds = _check_number(max_ds, "max_ds")
    if pairs not in ("all", "consecutive"):
        raise CoalToolError(f"pairs must be 'all' or 'consecutive', not {pairs!r}")
    max_bytes, max_results = dnds.kaks_limits()
    try:
        r = dnds.run_kaks(codon_alignment_path, core.new_run_dir(), chosen, int(codontable), pairs, min_c,
                          mds or None, hide_paths=False)
    except (ValueError, OSError) as e:
        raise CoalToolError(str(e)) from e
    flagged = _flagged(r.rows)
    keys, refs = dnds.cite(r.methods, codon_alignment=False)
    out_path = r.tsv_path
    if save_as and out_path:
        out_path = _save_copy(out_path, save_as, "table")
    return dict(DndsPairsResult(
        completed=r.ok,
        all_estimated=bool(r.ok and r.counts["total"] and r.counts["ok"] == r.counts["total"]),
        exit_code=r.exit_code,
        error=r.error,
        counts=StatusCounts(**{k: int(r.counts.get(k, 0)) for k in StatusCounts.__annotations__}),
        pairs=len({(x.get("seq1"), x.get("seq2")) for x in r.rows}),
        methods=r.methods,
        input_format=r.fmt,
        output_path=out_path,
        rows=[_row(x) for x in r.rows[:MAX_ROWS]],
        rows_total=len(r.rows),
        flagged=flagged[:MAX_FLAGGED],
        flagged_total=len(flagged),
        statuses=[StatusMeaning(**s) for s in dnds.status_help(r.counts)],
        engine_messages=dnds.engine_lines(r.stderr),
        cite=keys,
        references=refs,
        next_step=_kaks_next_step(r),
        coal_command=r.command,
        limits={"max_bytes": max_bytes, "max_results": max_results},
    ))


def _save_copy(path: str, save_as: str, what: str) -> str:
    """Copy `path` to `save_as` (a file path, or an existing folder to save in under the same name); never
    overwrites. Returns the new path."""
    import shutil
    dest = save_as
    if os.path.isdir(dest):
        dest = os.path.join(dest, os.path.basename(path))
    if os.path.exists(dest):
        raise CoalToolError(f"not saved: {dest} already exists (COAL never overwrites a file; choose another name). "
                            f"The {what} is at {path}")
    parent = os.path.dirname(dest)
    if not os.path.isdir(parent):
        raise CoalToolError(f"not saved: the folder {parent} does not exist. The {what} is at {path}")
    try:
        shutil.copyfile(path, dest)
    except OSError as e:
        raise CoalToolError(f"not saved to {dest}: {e.strerror or e}. The {what} is at {path}") from e
    return dest


def _pipeline_out_dir(out_dir: str) -> str:
    out_dir = (out_dir or "").strip()
    if not out_dir:
        try:
            return core.new_run_dir()
        except ValueError as e:
            raise CoalToolError(str(e)) from e
    if not os.path.isabs(out_dir):
        raise CoalToolError(f"out_dir must be an absolute path, not {out_dir!r}")
    if os.path.isdir(out_dir):
        if os.listdir(out_dir):
            raise CoalToolError(f"out_dir {out_dir} is not empty: give a new or empty folder (COAL never "
                                "overwrites results)")
        return out_dir
    if os.path.exists(out_dir):
        raise CoalToolError(f"out_dir {out_dir} is a file, not a folder")
    if not os.path.isdir(os.path.dirname(out_dir)):
        raise CoalToolError(f"the folder {os.path.dirname(out_dir)} (where out_dir would be made) does not exist")
    try:
        os.mkdir(out_dir)
    except OSError as e:
        raise CoalToolError(f"could not create out_dir {out_dir}: {e.strerror or e}") from e
    return out_dir


def _pipeline_next_step(r: dnds.PipelineRun) -> str:
    if r.status != "completed":
        return f"The pipeline did not finish: {r.error or 'see progress'}. Fix the cause and call pipeline_pairs again."
    c = r.counts
    if not c.get("total"):
        return ("No ortholog pair was found or estimated: check that both files hold coding sequences of the same "
                "genetic code, and read report_path and the failures.")
    if r.exit_code == 0 and c.get("ok") == c.get("total"):
        return "Every pair is estimated; the table is outputs['dnds'] (and everything in zip_path). Cite references."
    return ("Some pairs or genes need care: read flagged, failures (failed.tsv) and statuses, and tell the person "
            "which ones are affected before using the numbers.")


def pipeline_pairs(
    query_cds_path: Annotated[str, Field(
        description="Absolute path to the coding sequences (CDS, unaligned FASTA, may be .gz) of the first species. "
                    "Up to 500 MB.")],
    subject_cds_path: Annotated[str, Field(
        description="Absolute path to the coding sequences of the second species. Up to 500 MB.")],
    methods: Annotated[list[KaksMethod], Field(
        description="Ka/Ks methods: YN (default), NG, LWL, MLWL, LPB, MLPB, MYN, Comeron, Li.")]
    = ["YN"],  # noqa: B006 - never mutated
    codontable: Annotated[GeneticCode, Field(
        description="NCBI genetic code for translation, the codon alignment and Ka/Ks (1 standard, 11 bacteria).")]
    = 1,
    longest_isoform: Annotated[bool, Field(
        description="Keep one record per gene, the longest CDS (for files with several transcripts per gene).")]
    = False,
    ortho: Annotated[Literal["rbh", "bh"], Field(
        description="rbh: reciprocal best hits (default); bh: every query's best hit.")] = "rbh",
    pairs_path: Annotated[str, Field(
        description="Optional absolute path to a two-column TSV of ortholog pairs (query ID, subject ID): the given "
                    "pairs are used and DIAMOND is not needed.")] = "",
    min_codons: Annotated[int, Field(
        description="Flag (status warning) estimates from fewer compared codons than this. 0: off.", ge=0)] = 0,
    max_ds: Annotated[float, Field(
        description="Flag (status warning) estimates with dS above this (saturation). 0: off.", ge=0)] = 0,
    keep_alignments: Annotated[bool, Field(
        description="Also write the codon alignments (codon_alignments.fasta).")] = False,
    out_dir: Annotated[str, Field(
        description="Optional absolute path of a new or empty folder for the results (e.g. next to the inputs). "
                    "Default: a run folder of COAL's, deleted after 24 hours.")] = "",
    threads: Annotated[int, Field(description="Threads; 0: one per CPU.", ge=0)] = 0,
) -> dict:
    """Run COAL's one-step dN/dS pipeline for two species (`coal pipeline pairs`): translate both CDS files, find
    ortholog pairs with DIAMOND (reciprocal best hits), align each pair, build the codon alignment and estimate
    Ka/Ks.

    Use this when the person has the coding sequences of two species (unaligned) and wants dN/dS per ortholog pair.
    It needs DIAMOND on this machine (COAL looks at $COAL_DIAMOND, PATH and its cache); without it the tool error
    says how to install it (`coal pipeline install diamond`), unless pairs_path gives the pairs. A run takes
    seconds for hundreds of genes and minutes for whole genomes (limit: 500 MB per file and one hour; larger runs
    belong on the command line). For a codon alignment you already have, use dnds_pairs.

    Args:
        query_cds_path: Absolute path to the first species' coding sequences (FASTA, may be gzip-compressed).
        subject_cds_path: Absolute path to the second species' coding sequences.
        methods: Ka/Ks methods, default ["YN"].
        codontable: NCBI genetic code number.
        longest_isoform: Keep only the longest transcript per gene.
        ortho: "rbh" (reciprocal best hits, default) or "bh" (best hits).
        pairs_path: Optional TSV of given ortholog pairs; then DIAMOND is not run.
        min_codons: Flag estimates from fewer compared codons than this. 0: off.
        max_ds: Flag estimates with dS above this. 0: off.
        keep_alignments: Also write the codon alignments.
        out_dir: Optional absolute path of a new or empty folder for the results.
        threads: Threads (0: all CPUs).

    Returns:
        A dict with: completed (status is completed), status (completed / failed / interrupted), exit_code
        (0 nothing failed, 2 some genes or pairs could not be used, 1 the run failed), error, counts (the report's
        counts of genes, orthologs and pairs), pair_counts (results per status), rows (the first 20 rows of
        dnds.tsv; rows_total), flagged (rows with a status other than ok; flagged_total), failures (the first
        entries of failed.tsv; failures_total), statuses (what each status that occurs means), out_dir, outputs
        (name -> path: dnds, failed, genes, orthologs, unpaired, ...), report_path (report.json, schema
        coal.pipeline/1), zip_path (all result files), diamond (the DIAMOND used), cite and references,
        next_step, coal_command, progress (the pipeline's own progress lines) and limits. A missing tool, too
        large an input or an out_dir that is not empty is a tool error.
    """
    if core.hosted():
        raise CoalToolError("the dN/dS pipeline is not available on this shared server; install COAL locally")
    for p in (query_cds_path, subject_cds_path):
        _check_input(p)
    pairs_path = (pairs_path or "").strip()
    if pairs_path:
        _check_input(pairs_path)
    if int(codontable) not in core.GENETIC_CODES:
        raise CoalToolError(f"unknown genetic code {codontable} (supported: "
                            f"{', '.join(str(c) for c in core.GENETIC_CODES)})")
    chosen = _check_methods(methods, codontable)
    min_c = int(_check_number(min_codons, "min_codons"))
    mds = _check_number(max_ds, "max_ds")
    thr = int(_check_number(threads, "threads"))
    if ortho not in ("rbh", "bh"):
        raise CoalToolError(f"ortho must be 'rbh' or 'bh', not {ortho!r}")
    try:
        dnds.check_cds_input(query_cds_path, "query", dnds.PIPELINE_MAX_BYTES_MCP)
        dnds.check_cds_input(subject_cds_path, "subject", dnds.PIPELINE_MAX_BYTES_MCP)
    except (ValueError, OSError) as e:
        raise CoalToolError(str(e)) from e
    st = dnds.tool_status(need_diamond=not pairs_path)
    if not st.cpp:
        raise CoalToolError(dnds.NEEDS_CPP)
    if not pairs_path and not st.diamond:
        raise CoalToolError(dnds.NEEDS_DIAMOND + " Or pass pairs_path with the ortholog pairs.")
    out = _pipeline_out_dir(out_dir)
    cmd = dnds.pipeline_command(os.path.abspath(query_cds_path), os.path.abspath(subject_cds_path), out, chosen,
                                int(codontable), bool(longest_isoform), ortho, min_c, mds or None,
                                bool(keep_alignments), os.path.abspath(pairs_path) if pairs_path else None,
                                st.diamond, thr)
    r = dnds.run_pipeline(cmd, out, dnds.PIPELINE_TIMEOUT_MCP)
    flagged = _flagged(r.rows)
    keys = list(r.report.get("cite") or [])
    _k, refs = dnds.cite(chosen, codon_alignment=True, diamond=not pairs_path)
    if keys:
        try:
            from coal import citations
            refs = citations.format_references(keys)
        except (ImportError, KeyError):
            pass
    else:
        keys = _k
    return dict(PipelinePairsResult(
        completed=r.status == "completed",
        status=r.status,
        exit_code=r.exit_code,
        error=r.error,
        counts=dict(r.report.get("counts") or {}),
        pair_counts=StatusCounts(**{k: int(r.counts.get(k, 0)) for k in StatusCounts.__annotations__}),
        rows=[_row(x) for x in r.rows[:MAX_ROWS]],
        rows_total=len(r.rows),
        flagged=flagged[:MAX_FLAGGED],
        flagged_total=len(flagged),
        failures_total=int(r.report.get("failures_total") or 0),
        failures=list(r.report.get("failures") or [])[:MAX_FLAGGED],
        statuses=[StatusMeaning(**s) for s in dnds.status_help(r.counts)],
        out_dir=out,
        outputs=dict(r.outputs),
        report_path=r.report_path,
        zip_path=r.zip_path,
        diamond=None if pairs_path else st.diamond,
        cite=keys,
        references=refs,
        next_step=_pipeline_next_step(r),
        coal_command=r.command,
        progress=r.log[-30:],
        limits={"max_bytes_per_file": dnds.PIPELINE_MAX_BYTES_MCP, "timeout_seconds": dnds.PIPELINE_TIMEOUT_MCP},
    ))
