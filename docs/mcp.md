# AI agents (MCP)

COAL's stdio MCP server gives an agent two tools: `codon_align` (protein alignment + coding DNA -> codon alignment,
with a report; `save_as` puts the alignment where you want it) and `explain_report` (that report in plain English).

**Windows app:** keep **Connect COAL to Claude Desktop** ticked on the installer's last page (or run
`COAL.exe --connect-claude` later), then quit and reopen Claude Desktop. This adds a `coal` entry that starts
`COAL.exe -m coal_app --mcp` to `claude_desktop_config.json`, keeping everything else in it (the old file is kept as
`claude_desktop_config.json.bak-coal`); uninstalling COAL removes the entry again.

**From source:** the MCP server is part of the app package (`apps/coal-app`, `coal-align-app`); install it with
the engine in its **own** virtual environment (not the one with the `gui` extra), from the repository root:

```
python -m venv .venv-mcp
.venv-mcp/bin/python -m pip install -e . -e "apps/coal-app[mcp]"   # Windows: apps\coal-app\desktop\setup_mcp.bat
claude mcp add coal -- "$PWD/.venv-mcp/bin/python" -m coal_app --mcp
```

Without a checkout: `pip install "coal-align @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz" "coal-align-app[mcp] @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_app_src.tar.gz"`.

For Claude Desktop, add this to `claude_desktop_config.json` (Settings > Developer > Edit Config), with the full path
of that environment's Python (on Windows, `apps\coal-app\desktop\setup_mcp.bat` prints it ready to paste,
backslashes escaped):

```json
{"mcpServers": {"coal": {"command": "/full/path/to/COAL/.venv-mcp/bin/python", "args": ["-m", "coal_app", "--mcp"]}}}
```

Give the tools **absolute paths**. `codon_align` returns `completed` (COAL finished), `all_ok` / `needs_attention`
(whether every sequence matched its protein), `counts`, `suggested_codontable`, `next_step`, `output_path` and
`report_path`; a missing file, an unknown genetic code or an unsupported format is a tool error. The web GUI
(`python -m coal_app`) also serves the tools at `http://127.0.0.1:7860/gradio_api/mcp/sse`; the Windows app does
not by default (`COAL.exe --mcp`). Agent instructions: `.claude/skills/coal/SKILL.md`; a stage demo:
`demo/mcp_demo.md`.
