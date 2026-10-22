---
title: COAL codon alignment
colorFrom: blue
colorTo: green
sdk: gradio
sdk_version: 5.49.1
app_file: app.py
pinned: false
license: mit
short_description: Protein alignment + coding DNA -> codon alignment (PAL2NAL)
---

# COAL: codon alignment made easy

Upload a protein alignment and the coding DNA of the same genes; COAL lines the DNA up codon by
codon, shows a plain-English health check per gene and an alignment view that marks exactly where
a gene breaks. A fast, MIT-licensed alternative to PAL2NAL.

Files you upload are processed on this shared Hugging Face server and deleted within 24 hours (or when
the Space restarts). For unpublished or sensitive data, install COAL on your own computer instead.

AI agents: the MCP tools (`codon_align`, `explain_report`) work on files on the computer that runs them,
so use them locally: COAL is not on PyPI (a package called "coal" there is something else), so install it
from the repository (`git clone https://github.com/plumat-coal/COAL`, then `pip install -e ./COAL -e "./COAL/apps/coal-app[mcp]"`),
then run `python -m coal_app --mcp`, or `coal-gui`, which also serves MCP.
This public Space does not serve MCP.

Source: https://github.com/plumat-coal/COAL
