# coal-align-app: the COAL web GUI and MCP server

The app for [COAL](https://github.com/plumat-coal/COAL), a fast, MIT-licensed codon aligner (a PAL2NAL alternative):

- a web GUI on this computer (`coal-gui`, or `python -m coal_app`): a per-gene health check, an alignment view that
  marks where a gene breaks, genetic-code suggestions and a ready-to-paste methods sentence;
- an MCP server for AI agents (`python -m coal_app --mcp`), with the tools `codon_align` and `explain_report`;
- the Windows app (`desktop/`) and the Hugging Face Space (`space/`), built from this package.

The engine is a separate package, `coal-align` (the repository root), which this one depends on at exactly the same
version. Pipelines that only need the `coal` / `pal2nal` commands install `coal-align` alone.

## Install

Neither package is on PyPI yet. From the latest release:

```
pip install "coal-align @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz" \
    "coal-align-app[gui] @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_app_src.tar.gz"
```

From a checkout, at the repository root: `pip install -e . -e "apps/coal-app[gui]"`. Use `[mcp]` instead of `[gui]`
for the stand-alone MCP server, in a separate environment (the two extras pin conflicting `mcp` versions).

## Tests

From the repository root: `python -m unittest discover -s apps/coal-app/tests`, and in `apps/coal-app/desktop/`:
`python -m unittest test_wrapper test_smoke`. More: [docs/install.md](../../docs/install.md),
[docs/mcp.md](../../docs/mcp.md), [docs/development.md](../../docs/development.md).
