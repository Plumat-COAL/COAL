"""Stand-alone MCP server over stdio (no Gradio needed): ``python -m coal_app --mcp``.

Uses the official ``mcp`` Python SDK, version 2.x (``pip install -e "apps/coal-app[mcp]"``), or 1.x (the one Gradio's ``mcp`` extra
installs, so the one inside the Windows app: ``COAL.exe -m coal_app --mcp``).
Example client config (Claude Code / Claude Desktop), with the absolute path of the Python that has COAL:
    {"mcpServers": {"coal": {"command": "/abs/path/.venv-mcp/bin/python", "args": ["-m", "coal_app", "--mcp"]}}}
"""
from __future__ import annotations

import inspect

from . import tools

# What the server covers, said in the instructions, the GUI's About and the README (the GUI's tabs are the same four
# jobs): codon alignment, Ka/Ks on alignments, and the pairs pipeline; the pipeline's cluster mode and `coal kaks`
# (KaKs_Calculator's format) stay in the `coal` command line.
SCOPE = ("Tools: codon_align (protein alignment + coding DNA -> codon alignment), explain_report, dnds_pairs (Ka/Ks, "
         "i.e. dN/dS, for pairs of an aligned codon FASTA or an AXT file) and pipeline_pairs (two species' coding "
         "sequences -> ortholog pairs with DIAMOND -> dN/dS per pair; needs DIAMOND: if it is missing, tell the "
         "person to run `coal pipeline install diamond`). For dN/dS from a protein alignment and DNA, call "
         "codon_align first, then dnds_pairs on its output_path. Each tool states its input limits; larger inputs "
         "and the cluster mode of the pipeline belong on the `coal` command line (`coal dnds`, `coal pipeline`, "
         "`coal kaks`; `coal help`): tell the person to run the command, or run it yourself if you have a shell.")

INSTRUCTIONS = (
    "COAL turns a protein multiple sequence alignment plus the unaligned coding DNA of the same genes "
    "into a codon alignment (a fast PAL2NAL replacement), and estimates Ka/Ks (dN/dS). Call the tools with "
    "absolute file paths. "
    "Use an alignment only if all_ok is true; otherwise call explain_report on the returned report_path. "
    "If suggested_codontable is set, re-run codon_align with it. "
    "When the person names a folder, pass save_as (codon_align, dnds_pairs) or out_dir (pipeline_pairs) so they find "
    "the results there. Read the status of every Ka/Ks row (ok, warning, not_converged, saturated, no_data, error) "
    "before using its numbers. "
    + SCOPE
)


def server_version() -> str:
    """serverInfo.version: the version of the COAL engine actually imported
    (not the installed distribution's metadata, which can be stale for an editable or PYTHONPATH checkout)."""
    try:
        from coal import __version__
        return str(__version__)
    except ImportError:  # pragma: no cover
        return "unknown"


def _structured(fn, result_type):
    """`fn` with its result TypedDict as the return annotation, so mcp publishes an outputSchema.
    The functions in tools keep `-> dict` for Gradio's gr.api."""

    def tool(*args, **kwargs):
        return fn(*args, **kwargs)

    tool.__name__ = tool.__qualname__ = fn.__name__
    tool.__doc__ = fn.__doc__
    tool.__module__ = fn.__module__
    tool.__signature__ = inspect.signature(fn, eval_str=True).replace(return_annotation=result_type)
    return tool


def _structured_codon_align():
    return _structured(tools.codon_align, tools.CodonAlignResult)


def _add_tools(server) -> None:
    server.add_tool(_structured_codon_align(), name="codon_align")
    server.add_tool(tools.explain_report)
    server.add_tool(_structured(tools.dnds_pairs, tools.DndsPairsResult), name="dnds_pairs")
    server.add_tool(_structured(tools.pipeline_pairs, tools.PipelinePairsResult), name="pipeline_pairs")


def build_server():
    try:
        from mcp.server.mcpserver import MCPServer  # mcp >= 2.0
    except ImportError:
        return _build_fastmcp()

    server = MCPServer(name="coal", title="COAL codon alignment and Ka/Ks", instructions=INSTRUCTIONS,
                       version=server_version())
    _add_tools(server)
    return server


def _build_fastmcp():
    """The same server on mcp 1.x (FastMCP). Gradio 5.49.1's `mcp` extra pins mcp==1.10.1, so that is the mcp
    inside the Windows app, which runs this server for Claude Desktop (`COAL.exe -m coal_app --mcp`).
    FastMCP has no version parameter; it reports the mcp library's own version unless the low-level server's is
    set, so that is set to COAL's. log_level: FastMCP's default INFO logs tens of kilobytes to stderr per session
    (every request), which Claude Desktop keeps in its MCP log."""
    from mcp.server.fastmcp import FastMCP  # mcp 1.x

    server = FastMCP(name="coal", instructions=INSTRUCTIONS, log_level="WARNING")
    server._mcp_server.version = server_version()  # serverInfo.version
    _add_tools(server)
    return server


def main() -> int:
    from . import core
    core.prune_runs()  # old run folders are removed at startup
    build_server().run("stdio")
    return 0
