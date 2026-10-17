"""Stand-alone MCP server over stdio (no Gradio needed): ``python -m coal_app --mcp``.

Uses the official ``mcp`` Python SDK, version 2.x (``pip install -e .[mcp]``).
Example client config (Claude Code / Claude Desktop):
    {"mcpServers": {"coal": {"command": "python", "args": ["-m", "coal_app", "--mcp"]}}}
"""
from __future__ import annotations

from . import tools

INSTRUCTIONS = (
    "COAL turns a protein multiple sequence alignment plus the unaligned coding DNA of the same genes "
    "into a codon alignment (a fast PAL2NAL replacement). Call codon_align with file paths, then "
    "explain_report on the returned report_path if anything is not OK. If hints suggest another "
    "genetic code (suggested_codontable), re-run codon_align with it."
)


def build_server():
    from mcp.server.mcpserver import MCPServer  # mcp >= 2.0

    server = MCPServer(name="coal", title="COAL codon alignment", instructions=INSTRUCTIONS)
    server.add_tool(tools.codon_align)
    server.add_tool(tools.explain_report)
    return server


def main() -> int:
    from . import core
    core.prune_runs()  # CHANGED(APP-12): old run folders are removed at startup
    build_server().run("stdio")
    return 0
