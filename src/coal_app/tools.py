"""The two MCP tools AI agents use. Plain functions: served by Gradio's MCP server
(``python -m coal_app``) or by the stand-alone stdio server (``python -m coal_app --mcp``).

CHANGED(IMP-M-1): the result says what happened truthfully (``completed`` = COAL finished; ``all_ok`` /
``needs_attention`` = whether every sequence matched), errors raise :class:`CoalToolError` (MCP clients see
``isError: true``), the arguments carry descriptions and enums, and stderr is summarised instead of copied.
See overnight/integration/CHANGES.md
"""
from __future__ import annotations

import os
from typing import Annotated, Literal

from . import core

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
except ImportError:  # Gradio's MCP server (mcp 1.x) or no MCP at all: any exception becomes isError
    _ToolErrorBase = Exception


class CoalToolError(_ToolErrorBase, ValueError):
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
        pal2nal_caveats=core.pal2nal_caveats(r.report.get("options") or {}),  # CHANGED(FIN-A5/FIN-C8)
    )


def _check_input(path: str) -> None:
    # CHANGED(APP-11): was any server path accepted; now only sandbox paths when hosted. See overnight/app/CHANGES.md
    why = core.check_path_allowed(path)
    if why:
        raise CoalToolError(why)
    if not os.path.isfile(path):
        raise CoalToolError(f"file not found: {path} (give an absolute path; relative paths resolve against "
                            f"the server's working directory {os.getcwd()})")


# CHANGED(GUI-5): docstring only (strict): was "the run stops with an error and no alignment is written", but this tool
# runs COAL with --keep-going (core.run_coal default, MRG-3), so refused genes are left out and the run goes on.
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
    if int(codontable) not in core.GENETIC_CODES:
        raise CoalToolError(f"unknown genetic code {codontable} (supported: "
                            f"{', '.join(str(c) for c in core.GENETIC_CODES)})")
    try:
        # CHANGED(APP-12): was a /tmp/coal_* folder never cleaned; now a pruned folder in core.RUNS_DIR
        # CHANGED(FIN-A6): strict and blockonly passed on (were not available to agents)
        r = core.run_coal(protein_alignment_path, dna_fasta_path, output_format, int(codontable),
                          bool(nogap), bool(nomismatch), out_dir=core.new_run_dir(), strict=bool(strict),
                          blockonly=bool(blockonly))
    except ValueError as e:
        raise CoalToolError(str(e)) from e
    return dict(summarise(r))


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
    why = core.check_path_allowed(report_path)  # CHANGED(APP-11): sandboxed when hosted
    if why:
        raise CoalToolError(f"Could not read the report: {why}")
    try:
        rep = core.load_report(report_path)
    except (OSError, ValueError) as e:
        raise CoalToolError(f"Could not read the report: {e}") from e
    return core.explain(rep)
