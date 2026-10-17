# Overnight B (app/GUI/MCP): morning summary

Tasks 1-7 done, plus task 8 (alignment viewer) on branch `overnight/app` (local commits only, nothing pushed).

## Done
| Task | Commit |
|---|---|
| 1. Extras (`gui`, `mcp`), `coal-gui` entry point, CI workflow (3 OS x Py 3.10/3.13 + optional C++ job) | 987fbb3 |
| 2. Example datasets: clean, mito, messy_ids, one_broken_gene (+ generator) | 7a0ec15 |
| 3. `coal_app`: Gradio GUI, MCP tools `codon_align` / `explain_report`, stdio MCP server, sample report, tests | 4ce2b08 |
| 4. `.claude/skills/coal/SKILL.md` | dc1553e |
| 5. `pitch/explainer.html` (20 s, keyboard, reduced motion) | 4113091 |
| 6. `demo/race.py` (PAL2NAL per gene vs one batched COAL run) | e87f62f |
| 7. `coal_compat`: drop-in `pal2nal` command | 49d3a5c |
| Final log | 31d50b9 |
| 8. Alignment viewer: `coal_app/viewer.py` + `python -m coal_app.viewer`, GUI panel, jump links | (see log) |

## Test results (final run, Linux, Python 3.13)
- `python -m unittest discover -s tests`: 11 tests OK (engine untouched).
- `python -m unittest discover -s tests_app`: 34 tests OK, 2 skipped (Gradio GUI tests: Gradio not installable here).
- Windows/macOS not tested here; the CI workflow covers them on the first push.

## What switches on once the engine branch is merged (no app changes needed)
- `-output clustal|paml|codon`: detected from `coal --help`; the GUI dropdown, MCP tool and `pal2nal` command offer them automatically.
- `--report`: the app uses the engine's report instead of rebuilding one from stderr (and the engine's hints instead of the app's own genetic-code probing).
- `--keep-going` / exit code 2: passed automatically; exit 2 is treated as success with output, so one_broken_gene will give 19 aligned genes + 1 failed row.
- ID normalisation (`pairing: "normalised"`): messy_ids should then pass. The app's own ID-suffix hint stays as a fallback.
- C++ build: the CI `cpp` job runs the engine tests with `COAL_IMPL=cpp` once `cpp/CMakeLists.txt` exists; `demo/race.py` races a binary found in `cpp/build/`.

## Parked problems
- **No PyPI access in this sandbox** (pypi.org: DNS failure direct, 403 via the proxy, nothing cached). So: Gradio is not installed, **`src/coal_app/ui.py` has never run against real Gradio** (only checked with a mocked module), the pin `gradio[mcp]==5.49.1` was not checked against PyPI, and **`uv.lock` was not generated**.

## Needs a human decision / check
1. On a laptop: `pip install -e .[gui]`, then `python -m coal_app`; check the GUI and `python -m unittest discover -s tests_app` (runs the 2 GUI tests). Fix the Gradio pin to the current stable version if needed (written for the Gradio 5 API: `show_api=False`, `gr.api`, `launch(mcp_server=True)`; Gradio 6 renamed some of these). Then run `uv lock` and commit `uv.lock`.
2. `gradio[mcp]` 5.x probably needs `mcp` 1.x, while the `mcp` extra pins `mcp==2.2.0` (the stand-alone stdio server uses the 2.x `MCPServer` API). Don't install both extras in one environment, or drop one MCP route.
3. `pal2nal` command name (not `pal2nal.pl`, see log). PAL2NAL's default output is clustal, so `pal2nal pep nuc` without `-output fasta` errors until the engine has clustal: OK, or default to fasta with a warning?
4. Suggested README text (README is not mine to edit):
   > **GUI and AI agents.** `pip install -e .[gui]` then `python -m coal_app` (or `coal-gui`) opens a web page on your computer: upload files or try an example, see a plain-English health check per gene, download the result. The same server offers MCP tools (`codon_align`, `explain_report`) to AI agents; `python -m coal_app --mcp` is a stdio-only MCP server (`pip install -e .[mcp]`). **Old pipelines:** replace `pal2nal.pl` with `pal2nal` (installed with COAL). **Examples:** `examples/` has four small datasets with READMEs.
5. A real race against PAL2NAL (`python demo/race.py`) must be tried on a laptop with PAL2NAL fetched; here only the "not installed" path and a stand-in were tested.

## How to try it
Linux / macOS:
```
pip install -e ".[gui]"
python -m coal_app                      # GUI on http://127.0.0.1:7860 ; --demo shows the sample report
python -m coal_app --check              # which engine features were detected
python -m coal examples/mito/protein.aln.fasta examples/mito/dna.fasta -codontable 2
pal2nal examples/clean/protein.aln.fasta examples/clean/dna.fasta -output fasta
python demo/race.py                     # needs: python scripts/fetch_pal2nal.py
open pitch/explainer.html               # (Linux: xdg-open)
```
Windows (PowerShell):
```
py -m pip install -e ".[gui]"
py -m coal_app
py -m coal examples\mito\protein.aln.fasta examples\mito\dna.fasta -codontable 2
py -m coal_compat examples\clean\protein.aln.fasta examples\clean\dna.fasta -output fasta
py demo\race.py --perl C:\Strawberry\perl\bin\perl.exe
start pitch\explainer.html
```
MCP for Claude Code / Desktop (stdio): `{"mcpServers": {"coal": {"command": "python", "args": ["-m", "coal_app", "--mcp"]}}}`.
With the GUI running, Gradio serves MCP at `http://127.0.0.1:7860/gradio_api/mcp/sse`.

---

## Progress log

- Sat 02:06 start. Engine tests: 11 OK. Python 3.13 on Linux. cpp/ does not exist on this branch.
- Sat 02:06 Task 1, problem "install gradio": attempt 1 `pip index versions gradio` -> DNS failure for pypi.org (pypi.org is in NO_PROXY but does not resolve).
- Sat 02:06 Task 1, problem "install gradio": attempt 2 curl pypi.org via the agent proxy -> 403 connect_rejected (host not allowed for this session). Not routing around it.
- Sat 02:06 Task 1, problem "install gradio": attempt 3 look for a cached gradio wheel in ~/.cache/pip and ~/.cache/uv -> none. **PARKED**: no PyPI access in this sandbox, so Gradio cannot be installed, the pin cannot be checked against PyPI, and uv.lock cannot be generated (uv lock needs the index). Diagnosis: egress policy, not a config error. Workaround: app logic is written in a Gradio-free core (fully tested), the Gradio UI is written but untested tonight, and the MCP tools are also served with the locally installed `mcp` 2.2.0 package (API read from its installed source).
- Sat 02:07 Task 1 done: pyproject gains pinned extras `gui` (gradio[mcp]==5.49.1, unverified pin) and `mcp` (mcp==2.2.0, the version installed here), `coal-gui` entry point, package data for coal_app; .github/workflows/ci.yml (3 OS x Py 3.10/3.13, both test suites; C++ job skips cleanly when cpp/CMakeLists.txt is absent) parses with PyYAML. uv.lock not generated (parked, no index). Local editable install needs --no-build-isolation here (no PyPI).
- Sat 02:09 Task 2 done: examples/make_examples.py (seed 20261017, self-contained tables for codes 1/2, no coal import) writes clean (8 genes), mito (24 genes, every gene has an internal TGA=Trp so code 1 warns, code 2 is clean), messy_ids (10 genes, AT1G.....1 vs AT1G....., DNA shuffled so order pairing fails too), one_broken_gene (20 genes, OG0001_sp14 has 1 base deleted). All < 17 KB. Checked: clean and mito/code 2 outputs equal the construction. On main, messy_ids and one_broken_gene make COAL exit 1 (expected until the engine gains normalised IDs / --keep-going).
- Sat 02:14 Task 3 done (Gradio parts untested): src/coal_app/ = core.py (Gradio-free: feature detection from `coal --help`, run via subprocess, coal.report/1 from --report or rebuilt from stderr, genetic-code hint by probing codes 2,3,4,5,9,11,13,14,21 on inputs < 5 MB, ID-suffix hint, health table, plain-English explain), tools.py (MCP tools codon_align / explain_report), ui.py (Gradio 5 Blocks; system fonts, analytics off, Okabe-Ito palette + icon + word, glossary, example buttons, hint banner with one-click re-run, download buttons, equivalent commands), mcp_stdio.py (stdio MCP via installed mcp 2.2.0 MCPServer; handshake + tools/list verified by piping JSON-RPC with timeout 20), sample_report.json (20 seqs: 17 ok, 2 warning, 1 failed, code-2 hint). 25 app tests pass, 2 GUI tests skipped (no gradio). ui.py checked only with a mocked gradio module (wiring/NameErrors), not the real API.
- Sat 02:14 Task 4 done: .claude/skills/coal/SKILL.md (when to use, install, 5 commands, report fields, failure/fix table).
- Sat 02:17 Task 5 done: pitch/explainer.html (single file, no external resources, 5 steps x 4 s = 20 s autoplay; Space/arrows/PageUp/PageDown step, P pause; reduced-motion: no transitions or flashing). Rendered in headless Chromium at 1280x720, no JS errors; fixed one label-overlap bug after the first screenshot.
- Sat 02:18 Task 6 done: demo/race.py (PAL2NAL once per gene vs one batched COAL run, plus C++ binary if found in cpp/build or --cpp; live bars; times, speed-up and output match). Tested: missing-PAL2NAL path (unit test) and the race machinery with a stand-in .pl that forwards to COAL (40 genes, outputs MATCH). A real race against PAL2NAL needs testing on a teammate's laptop. C++ path untested (no cpp/ here; assumes the binary takes the same arguments as coal).
- Sat 02:19 Task 7 done: src/coal_compat (`pal2nal` entry point and `python -m coal_compat`): PAL2NAL v14 command line (options in any order, several DNA files) -> `python -m coal`, output to stdout; -blockonly, -html, unknown options and unsupported formats give clear errors. Named `pal2nal`, not `pal2nal.pl`: a console script called pal2nal.pl becomes pal2nal.pl.exe on Windows and may clash with the .pl -> Perl file association; can't verify Windows here. Since PAL2NAL's default output is clustal, calls without -output fail clearly until the engine has clustal (no silent format change).
- Sat 02:19 Final: all tests run once (engine 11 OK; app 34 OK, 2 skipped). Morning summary written. Stopping.
- Sat 02:21 Task 8 (alignment viewer) started. Note: on the current engine one_broken_gene produces no codon alignment (exit 1), so the viewer also locates breaks itself from protein + DNA (frameshift = codons stop matching and a +/-1/2 base shift makes them match again).
- Sat 02:25 Task 8, layout check 1 (headless Chromium 1280 px): every row of a failed run carried the per-row note "COAL could not align this gene" -> now only rows that break get a note (their break text).
- Sat 02:25 Task 8, layout check 2 (390 px phone): sequence-name column took most of the width -> narrower name column and hidden row notes under 600 px. Rechecked: OK, no page-level horizontal scroll.
- Sat 02:25 Task 8 done: src/coal_app/viewer.py (Gradio-free; self-contained HTML; Okabe-Ito tints by amino-acid property, letter always shown; mismatches outlined + the symbol for 'not equal'; frameshift double outline + hatching + lightning symbol; sticky names and column numbers; horizontal scroll; anchor links per problem cell, small optional JS highlights the whole column; capped at 200 seqs x 2,000 cols and 200,000 cells with a note). gcode.py holds NCBI tables (data copy, app still never imports coal). When COAL wrote no row for a gene, the viewer places its codons and finds the break: one_broken_gene -> frameshift at column 102 (1 base missing). Fallback reports now carry that column, so the health table, explain_report and the viewer agree. GUI: 'View alignment' panel under the health table (anchors work there without JS) + HTML download. 8 new tests (42 app tests OK, 2 skipped).
- Sat 02:26 Task 9 done (not deployed): src/coal_app/space_app.py (0.0.0.0:7860, MCP on), space/ (README front matter sdk gradio / sdk_version 5.49.1 = gui pin / app_file app.py, requirements.txt, app.py that puts src/ on sys.path and PYTHONPATH for the engine subprocess, build_space.py that assembles a Space folder, DEPLOY.md step by step). Tested: front matter matches the pin; build_space output contains everything and its copied engine runs. Space itself not run (no Gradio here).
- Sat 02:28 Task 10 done: demo/diamond_pipeline.py (longest isoform per gene from gene:/gene= tag or isoform suffix; reciprocal best hits with the user's diamond (makedb + blastp both ways, outfmt 6); pure-Python Needleman-Wunsch with BLOSUM62 + affine gaps (Gotoh, pointer traceback); one batch padded with trailing gaps so COAL runs once; split into pairs/*.codon.fasta; --keep-going/--report used if the engine has them). Clear error when diamond is missing. Tested with tests_app/diamond_stub.py (a 3-mer stand-in, NOT DIAMOND): 11 proteins -> 6 genes -> 5 RBH -> 5 pair alignments that translate back to the proteins. A real run needs DIAMOND on a laptop. Pure-Python alignment is O(n*m): about 0.1 s per 400-aa pair here; --max-length 5000 skips giants.
