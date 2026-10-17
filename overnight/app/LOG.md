# Overnight B (app/GUI/MCP) log

(Morning summary will be written here at the end.)

## Progress log

- Sat 02:06 start. Engine tests: 11 OK. Python 3.13 on Linux. cpp/ does not exist on this branch.
- Sat 02:06 Task 1, problem "install gradio": attempt 1 `pip index versions gradio` -> DNS failure for pypi.org (pypi.org is in NO_PROXY but does not resolve).
- Sat 02:06 Task 1, problem "install gradio": attempt 2 curl pypi.org via the agent proxy -> 403 connect_rejected (host not allowed for this session). Not routing around it.
- Sat 02:06 Task 1, problem "install gradio": attempt 3 look for a cached gradio wheel in ~/.cache/pip and ~/.cache/uv -> none. **PARKED**: no PyPI access in this sandbox, so Gradio cannot be installed, the pin cannot be checked against PyPI, and uv.lock cannot be generated (uv lock needs the index). Diagnosis: egress policy, not a config error. Workaround: app logic is written in a Gradio-free core (fully tested), the Gradio UI is written but untested tonight, and the MCP tools are also served with the locally installed `mcp` 2.2.0 package (API read from its installed source).
- Sat 02:07 Task 1 done: pyproject gains pinned extras `gui` (gradio[mcp]==5.49.1, unverified pin) and `mcp` (mcp==2.2.0, the version installed here), `coal-gui` entry point, package data for coal_app; .github/workflows/ci.yml (3 OS x Py 3.10/3.13, both test suites; C++ job skips cleanly when cpp/CMakeLists.txt is absent) parses with PyYAML. uv.lock not generated (parked, no index). Local editable install needs --no-build-isolation here (no PyPI).
