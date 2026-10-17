"""The two MCP tools AI agents use. Plain functions: served by Gradio's MCP server
(``python -m coal_app``) or by the stand-alone stdio server (``python -m coal_app --mcp``).
"""
from __future__ import annotations

import os

from . import core


def codon_align(protein_alignment_path: str, dna_fasta_path: str, output_format: str = "fasta",
                codontable: int = 1, nogap: bool = False, nomismatch: bool = False) -> dict:
    """Build a codon alignment from a protein alignment and the matching unaligned coding DNA (PAL2NAL replacement).

    Use this when you have aligned protein sequences and the DNA (CDS) of the same genes and need
    the DNA aligned codon by codon, for example before dN/dS analysis (PAML codeml, HyPhy).
    Sequences are paired by ID. Codons that don't match their amino acid are reported as warnings.

    Args:
        protein_alignment_path: Path to the protein multiple sequence alignment (FASTA or CLUSTAL).
        dna_fasta_path: Path to the unaligned coding DNA sequences (FASTA), one per protein.
        output_format: "fasta" (always available); "clustal", "paml" or "codon" if this COAL supports them.
        codontable: NCBI genetic code number: 1 = standard (nuclear), 2 = vertebrate mitochondrial,
            5 = invertebrate mitochondrial, 11 = bacteria/plastids. Also 3, 4, 6, 9, 10, 12-16, 21-23.
        nogap: If true, remove codon columns that contain a gap or a stop codon.
        nomismatch: If true, remove columns where a codon doesn't match its amino acid.

    Returns:
        A dict with: ok (bool), output_path (str or None; the codon alignment file), report_path
        (str; a coal.report/1 JSON file, pass it to explain_report), summary (counts of ok / warning /
        failed sequences and how they were paired), hints (list of suggestions, e.g. a better genetic
        code as {"code": "genetic_code", "suggested_codontable": 2, "text": ...}), problems (list of
        {id, status, explanation} for every sequence that is not OK), coal_command and pal2nal_command
        (the equivalent command lines), and stderr (COAL's raw messages).
    """
    for p in (protein_alignment_path, dna_fasta_path):
        if not os.path.isfile(p):
            return {"ok": False, "error": f"file not found: {p}", "output_path": None}
    try:
        r = core.run_coal(protein_alignment_path, dna_fasta_path, output_format, int(codontable),
                          bool(nogap), bool(nomismatch))
    except ValueError as e:
        return {"ok": False, "error": str(e), "output_path": None}
    problems = [{"id": s.get("id"), "status": s.get("status"), "explanation": core.sequence_sentence(s)}
                for s in r.report.get("sequences", []) if s.get("status") != "ok"]
    return {
        "ok": r.ok,
        "exit_code": r.exit_code,
        "output_path": r.output_path,
        "report_path": r.report_path,
        "summary": r.report.get("summary", {}),
        "summary_text": core.summary_sentence(r.report),
        "hints": r.hints,
        "problems": problems,
        "coal_command": r.coal_cmd,
        "pal2nal_command": r.pal2nal_cmd,
        "stderr": r.stderr[-4000:],
    }


def explain_report(report_path: str) -> str:
    """Explain a COAL run report (coal.report/1 JSON) in plain English: what failed, why, and how to fix it.

    Args:
        report_path: Path to the report JSON written by COAL (--report) or returned by codon_align.

    Returns:
        A short plain-English explanation: counts of OK / warning / failed sequences, how proteins and
        DNA were paired, one line per problem sequence with its likely cause (wrong genetic code,
        mismatched IDs, frameshift), and suggested fixes.
    """
    try:
        rep = core.load_report(report_path)
    except (OSError, ValueError) as e:
        return f"Could not read the report: {e}"
    return core.explain(rep)
