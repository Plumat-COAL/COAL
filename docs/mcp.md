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

## Scope: codon alignment only

The MCP server (and the GUI it belongs to) has two tools, `codon_align` and `explain_report`. **Ka/Ks (dN/dS) and
`coal pipeline` are not MCP tools**: an agent that is asked for dN/dS should run the `coal` command line
(`coal dnds`, `coal kaks`, `coal pipeline pairs` / `clusters`; [kaks.md](kaks.md), [pipeline.md](pipeline.md)) if
it has a shell, or tell the person the command. The server's instructions say so, so that an agent does not try
to get dN/dS out of `codon_align`.

## Planned: tools for Ka/Ks and the pipeline

Not built yet (a proposal, written down so that the names and shapes are settled before the code):

| Tool | Arguments | Returns |
|---|---|---|
| `dnds_pairs` | `codon_alignment` (path; pairs `all` or `consecutive`), or `cds_a`, `cds_b` and `pairs` (paths: two CDS files and a two-column table of IDs); `methods` (default `["YN"]`), `codontable`, `min_codons`, `max_ds`, `save_as` | `completed`, `output_path` (the pair table, `coal.dnds_pipeline.PAIR_COLUMNS`), `counts` by status, `flagged`, `cite` (keys and references), `next_step` |
| `pipeline_pairs` | `query`, `subject` (CDS files), `out_dir`, `longest_isoform`, `gene_map`, `methods`, `aligner`, `threads` | the summary of `report.json`: `status`, `exit_code`, `counts`, `warnings`, `outputs` (paths), `cite`, `report_path` |

Both would wrap `coal.dnds_pipeline` ([python-api.md](python-api.md#dnds-and-the-pipeline)): `codon_dnds` /
`run_pairs`, which return instead of raising for a failed run, so the tool result can carry the report's `error`.
Open points: `pipeline_pairs` needs DIAMOND on the machine (the tool would report `find_diamond()` and offer
`install_diamond()` only when asked), a run can take minutes (progress notifications, or a job handle), and the
tools write into a folder the person names (`out_dir` / `save_as`), never elsewhere. A cluster-mode tool would
follow the same shape.
