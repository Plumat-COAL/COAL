# Overnight B (app/GUI/MCP) log

(Morning summary will be written here at the end.)

## Progress log

- Sat 02:06 start. Engine tests: 11 OK. Python 3.13 on Linux. cpp/ does not exist on this branch.
- Sat 02:06 Task 1, problem "install gradio": attempt 1 `pip index versions gradio` -> DNS failure for pypi.org (pypi.org is in NO_PROXY but does not resolve).
- Sat 02:06 Task 1, problem "install gradio": attempt 2 curl pypi.org via the agent proxy -> 403 connect_rejected (host not allowed for this session). Not routing around it.
- Sat 02:06 Task 1, problem "install gradio": attempt 3 look for a cached gradio wheel in ~/.cache/pip and ~/.cache/uv -> none. **PARKED**: no PyPI access in this sandbox, so Gradio cannot be installed, the pin cannot be checked against PyPI, and uv.lock cannot be generated (uv lock needs the index). Diagnosis: egress policy, not a config error. Workaround: app logic is written in a Gradio-free core (fully tested), the Gradio UI is written but untested tonight, and the MCP tools are also served with the locally installed `mcp` 2.2.0 package (API read from its installed source).
- Sat 02:07 Task 1 done: pyproject gains pinned extras `gui` (gradio[mcp]==5.49.1, unverified pin) and `mcp` (mcp==2.2.0, the version installed here), `coal-gui` entry point, package data for coal_app; .github/workflows/ci.yml (3 OS x Py 3.10/3.13, both test suites; C++ job skips cleanly when cpp/CMakeLists.txt is absent) parses with PyYAML. uv.lock not generated (parked, no index). Local editable install needs --no-build-isolation here (no PyPI).
- Sat 02:09 Task 2 done: examples/make_examples.py (seed 20261017, self-contained tables for codes 1/2, no coal import) writes clean (8 genes), mito (24 genes, every gene has an internal TGA=Trp so code 1 warns, code 2 is clean), messy_ids (10 genes, AT1G.....1 vs AT1G....., DNA shuffled so order pairing fails too), one_broken_gene (20 genes, OG0001_sp14 has 1 base deleted). All < 17 KB. Checked: clean and mito/code 2 outputs equal the construction. On main, messy_ids and one_broken_gene make COAL exit 1 (expected until the engine gains normalised IDs / --keep-going).
