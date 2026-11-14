---
title: COAL codon alignment
colorFrom: blue
colorTo: green
sdk: gradio
sdk_version: 5.49.1
app_file: app.py
pinned: false
license: gpl-3.0
short_description: Protein alignment + coding DNA -> codon alignment (PAL2NAL)
---

# COAL: codon alignment made easy

Upload a protein alignment and the coding DNA of the same genes; COAL lines the DNA up codon by
codon, shows a plain-English health check per gene and an alignment view that marks exactly where
a gene breaks. A fast alternative to PAL2NAL; free software (GPL-3.0-only).

Files you upload are processed on this shared Hugging Face server and deleted within 24 hours (or when
the Space restarts). For unpublished or sensitive data, install COAL on your own computer instead.

This Space makes codon alignments. Its Ka/Ks tab needs COAL's C++ engine (`coal-cpp`), which the Space does not
have, and the dN/dS pipeline does not run on a shared server: for both, install COAL on your own computer (the
platform wheel or the Windows app; the pipeline also needs DIAMOND: `coal pipeline install diamond`).

AI agents: the MCP tools (`codon_align`, `explain_report`, `dnds_pairs`, `pipeline_pairs`) work on files on the
computer that runs them,
so use them locally: COAL is not on PyPI (a package called "coal" there is something else), so install it
from the repository (`git clone https://github.com/plumat-coal/COAL`, then `pip install -e ./COAL -e "./COAL/apps/coal-app[mcp]"`),
then run `python -m coal_app --mcp`, or `coal-gui`, which also serves MCP.
This public Space does not serve MCP.

Source: https://github.com/plumat-coal/COAL
