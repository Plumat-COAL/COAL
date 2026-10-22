"""Print how to register COAL's stdio MCP server with Claude, for the Python running this script.

    .venv-mcp\\Scripts\\python desktop\\mcp_config.py

Used by desktop\\setup_mcp.bat (IMP-M-2). The JSON is made with json.dumps, so Windows paths come out
escaped (``C:\\\\Users\\\\...``) and can be pasted into claude_desktop_config.json as they are.
"""
from __future__ import annotations

import json
import subprocess
import sys


def server_entry(python: str) -> dict:
    return {"command": python, "args": ["-m", "coal_app", "--mcp"]}


def claude_desktop_config(python: str) -> str:
    return json.dumps({"mcpServers": {"coal": server_entry(python)}}, indent=2)


def claude_code_command(python: str) -> str:
    return "claude mcp add coal -- " + subprocess.list2cmdline([python]) + " -m coal_app --mcp"


def main(python: str | None = None) -> int:
    python = python or sys.executable
    print("Claude Code (run once in a terminal):")
    print()
    print("  " + claude_code_command(python))
    print()
    print("Claude Desktop: Settings > Developer > Edit Config opens claude_desktop_config.json")
    print("(%APPDATA%\\Claude\\claude_desktop_config.json). Put this in it, or add the \"coal\" entry")
    print("to the \"mcpServers\" you already have, save, and restart Claude Desktop:")
    print()
    print(claude_desktop_config(python))
    print()
    print("Give Claude absolute paths to your files (for example C:\\data\\protein.aln.fasta).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
