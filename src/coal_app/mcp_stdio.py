"""Stand-alone MCP server over stdio (no Gradio needed): ``python -m coal_app --mcp``.

Uses the official ``mcp`` Python SDK, version 2.x (``pip install -e ".[mcp]"``), or 1.x (the one Gradio's ``mcp`` extra
installs, so the one inside the Windows app: ``COAL.exe -m coal_app --mcp``).
Example client config (Claude Code / Claude Desktop), with the absolute path of the Python that has COAL:
    {"mcpServers": {"coal": {"command": "/abs/path/.venv-mcp/bin/python", "args": ["-m", "coal_app", "--mcp"]}}}
"""
from __future__ import annotations

import inspect

from . import tools

INSTRUCTIONS = (
    "COAL turns a protein multiple sequence alignment plus the unaligned coding DNA of the same genes "
    "into a codon alignment (a fast PAL2NAL replacement). Call codon_align with absolute file paths. "
    "Use the alignment only if all_ok is true; otherwise call explain_report on the returned report_path. "
    "If suggested_codontable is set, re-run codon_align with it."
)


def server_version() -> str:
    """CHANGED(IMP-M-1): new. serverInfo.version (was empty): the version of the COAL engine actually imported
    (not the installed distribution's metadata, which can be stale for an editable or PYTHONPATH checkout)."""
    try:
        from coal import __version__
        return str(__version__)
    except ImportError:  # pragma: no cover
        return "unknown"


def _structured_codon_align():
    """codon_align with its result TypedDict as the return annotation, so mcp publishes an outputSchema.
    CHANGED(IMP-M-1): new. tools.codon_align itself keeps `-> dict` for Gradio's gr.api."""
    fn = tools.codon_align

    def codon_align(*args, **kwargs):
        return fn(*args, **kwargs)

    codon_align.__doc__ = fn.__doc__
    codon_align.__module__ = fn.__module__
    codon_align.__signature__ = inspect.signature(fn, eval_str=True).replace(
        return_annotation=tools.CodonAlignResult)
    return codon_align


def build_server():
    try:
        from mcp.server.mcpserver import MCPServer  # mcp >= 2.0
    except ImportError:
        return _build_fastmcp()

    server = MCPServer(name="coal", title="COAL codon alignment", instructions=INSTRUCTIONS,
                       version=server_version())
    server.add_tool(_structured_codon_align(), name="codon_align")
    server.add_tool(tools.explain_report)
    return server


def _build_fastmcp():
    """CHANGED(SHIP-18): new. The same server on mcp 1.x (FastMCP). Gradio 5.49.1's `mcp` extra pins mcp==1.10.1, so
    that is the mcp inside the Windows app, which now runs this server for Claude Desktop (`COAL.exe -m coal_app
    --mcp`). Was: mcp 2 only (ImportError under 1.x). mcp 1.x has no server version field."""
    from mcp.server.fastmcp import FastMCP  # mcp 1.x

    server = FastMCP(name="coal", instructions=INSTRUCTIONS)
    server.add_tool(_structured_codon_align(), name="codon_align")
    server.add_tool(tools.explain_report)
    return server


def main() -> int:
    from . import core
    core.prune_runs()  # CHANGED(APP-12): old run folders are removed at startup
    build_server().run("stdio")
    return 0
