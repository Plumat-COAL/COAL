# Review fixes (branch fix/app-review, from keanu_overnight/app; local commits only, nothing pushed)

Two independent reviews (app findings 1-19, config/CI/docs findings C1-C15) plus user decisions (APP-L large
files by path, APP-M MCP default). Every change has a row in `overnight/app/CHANGES.md` and a
`CHANGED(<ID>)` comment at the code site. Each finding first got a failing test in `tests_app/`
(mostly `tests_app/test_review.py`), then the fix.

| Finding | Status | Commit |
|---|---|---|
| 1 codon-only viewer marks/links (blocker) | fixed | a9636bc |
| 2 tolerant start search (blocker) | fixed | a9636bc |
| 3 every engine warning counted (blocker) | fixed | a9636bc |
| 4 viewer status without codon file | fixed | c7e76e9 |
| 5 IUPAC codons like the engine | fixed | c7e76e9 |
| 6 no pairing sentence when nothing was paired | fixed | c7e76e9 |
| 7 ID-suffix hint (real DNA ID, verdict by outcome) | fixed | c7e76e9 |
| 8 race: PAL2NAL once per gene family | fixed (coal_compat timeout left out on purpose: drop-in for long runs) | 1444e98 |
| 9 diamond pipeline: all pairs skipped | fixed | 8c288e4 |
| 10 NW memory (3000x3000: 15 s/1.25 GB -> 4.3 s/~28 MB), --max-length 2000 | fixed | 8c288e4 |
| 11 hosted Space: intro text, path sandbox, MCP off on Space | fixed | e884ce4 |
| 12 temp dirs (RUNS_DIR, 24 h pruning, probe files, tests clean up) | fixed | e884ce4 |
| 13 genetic-code probes (sample <= 200, 30 s each, 90 s total; 22 s -> 1.6 s) | fixed; probes are not run with -q (they count warnings) | e884ce4 |
| 14 malformed / BOM reports | fixed | 655a862 |
| 15 sample report = real run of examples/mito with code 1, labelled demo | fixed | 655a862 |
| 16 CI gui/mcp jobs, required-extras switch, stdio MCP handshake, unconditional asserts | fixed; gui/mcp jobs not yet run (nothing pushed); the gradio[mcp] 5.49.1 vs mcp 2.2.0 question stays **for a human / CI** | a879b52 |
| 17 "COAL's app estimates" wording; mcp_demo.md regenerated | fixed | 655a862 |
| 18 utf-8 everywhere | fixed | 8c288e4 |
| 19 gcode drift test; locate_breaks above main | fixed | 655a862 |
| C1 PAL2NAL also a web form | fixed | 3dbe45d |
| C2 drop-in claim limited to -output fasta; -html hint; -nostderr accepted | fixed (-verbose still rejected: PAL2NAL source not available to confirm it) | 95d73bd |
| C3 165-check split recomputed | fixed (numbers below) | 3dbe45d |
| C4 honest C++ CI job | fixed | a879b52 |
| C5 extras CI jobs | fixed (= 16) | a879b52 |
| C6 explainer columns at step 1 | fixed (Playwright test) | 3dbe45d |
| C7 "TGA = W (mito)" on the last step | fixed | 3dbe45d |
| C8 explainer at 390 px | fixed | 3dbe45d |
| C9 "10 s with 4 processes, 25 on one" | fixed | 3dbe45d |
| C10 SKILL.md: current engine only | fixed | 4351f4a, 4d89937 |
| C11 extras in separate environments | documented; real resolution **needs CI / a human** | 4351f4a |
| C12 examples shipped in the wheel | fixed (wheel built + installed offline in a test) | a3ff7bc |
| C13 public Space path sandbox | fixed (= 11) | e884ce4 |
| C14 this log | fixed (here) | (this commit) |
| C15 explainer "about 16 s" in PITCH.md | fixed | 3dbe45d |
| APP-L1/L2 local path inputs, `?protein=&dna=` URL hook | done | be69a5c |
| APP-M1 MCP on by default locally, off on the Space | done (decision recorded in CHANGES.md) | be69a5c |

Note: commit 4351f4a was made with one failing app test (my own SKILL.md marker comment named the upcoming
flags); fixed in the next commit 4d89937. All other commits had both suites passing.

**C3 numbers** (recomputed by instrumenting `tests/test_differential.py::test_all_cases`): 165 checks =
143 compare stdout byte for byte (94 of them also compare the mismatch warnings; 47 are -nomismatch variants
whose warnings are not compared; 2 are table-10 deliberate warning-only differences) + 12 where PAL2NAL failed
and COAL must exit non-zero with empty output (error text not compared) + 3 deliberate invalid-character
errors (exit code only) + 7 deliberate "differs" (tables 10/11/13, exit code 0 only; 4 of them happen to give
identical stdout). Golden .rc files are never read.

**Tests now:** engine `python -m unittest discover -s tests`: 11 OK. App `python -m unittest discover -s
tests_app`: 115 tests OK, 5 skipped here (3 need Gradio, 1 needs a non-root user, 1 needs network or
`COAL_WHEELHOUSE` to build a wheel). With `COAL_WHEELHOUSE=<folder with a setuptools wheel>` the wheel test
runs and passes. On CI the plain jobs will additionally skip the 2 MCP-server tests (no `mcp`) and the 2
Playwright explainer tests; the new gui/mcp jobs require their tests. CI did run before these fixes: green
runs exist for keanu_overnight/app (`gh run list`), each with 3 app tests skipped (2 GUI + 1 MCP), not 2.

**Needs a human**
1. Push the branch and read the new CI jobs: `gui` (first real run of ui.py and the Gradio 5.49.1 pin, incl.
   the new path boxes and URL hook, only mock-tested here), `mcp`, and "Info: gui + mcp extras in one
   environment" (does gradio[mcp] conflict with mcp==2.2.0?).
2. README.md corrections (not mine to edit): line 69 says every check requires byte-identical output and equal
   warnings; true only for 143 (output) and 94 (warnings) of the 165, see C3 numbers. Line 10 "byte-identical
   output on the benchmarked datasets" is fine. Suggested README text for the GUI (earlier summary item 4)
   should now also mention path boxes for big files and that the public Space has no MCP.
3. Race (`demo/race.py`) CLI changed (`--families/--seqs/--codons`, `--family PEP NUC`); measure it on a
   laptop with PAL2NAL before quoting any race number. No race number is quoted in the pitch today.
4. PAL2NAL `-verbose`: if PAL2NAL v14 has it, add it to coal_compat's accepted-and-ignored flags.
5. Engine-side (not app): unrelated protein and DNA currently exit 0 with per-codon warnings; SKILL.md now
   tells agents to count warnings per gene.

---

# Overnight B (app/GUI/MCP): morning summary

All 13 tasks done (1-7 from the brief, 8-13 added overnight) on branch `overnight/app` (local commits only, nothing pushed).

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
| 8. Alignment viewer: `coal_app/viewer.py` + `python -m coal_app.viewer`, GUI panel, jump links, frameshift location | ab2c8cd |
| 9. Hugging Face Space files in `space/` (prepared, NOT deployed; `space/DEPLOY.md`) | 6d69028 |
| 10. `demo/diamond_pipeline.py`: two species -> longest isoforms -> DIAMOND RBH -> pairwise NW -> one COAL run | aff4a4f |
| 11. Ready-to-paste methods paragraph (GUI box with copy button, `codon_align`, `explain_report`) | 1abff5b |
| 12. `demo/mcp_demo.md` 60-second MCP stage script with real tool outputs (`demo/make_mcp_demo.py` regenerates it) | a2351d0 |
| 13. `pitch/PITCH.md` (3-minute script with timings + backup plan), `pitch/QA.md` (10 judge questions), `pitch/backup/*.png` | (final commit) |

## Test results (final run, Linux, Python 3.13)
- `python -m unittest discover -s tests`: 11 tests OK (engine untouched).
- `python -m unittest discover -s tests_app`: 52 tests OK, 2 skipped (Gradio GUI tests: Gradio not installable here).
- Viewer checked visually in headless Chromium (1280 px and 390 px phone); 2 layout fixes, see log.
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
5. Real runs still to try on a laptop: the DIAMOND pipeline with real `diamond` (tested here only with a 3-mer stand-in), the Space deploy (`space/DEPLOY.md`), and the MCP demo in a real Claude client.
6. Pitch: "165 golden checks byte-identical" is not accurate. It is 165 checks: 143 byte-identical, 12 where both tools reject the input, 10 documented deliberate differences (PITCH.md and QA.md say this). CI has never run (nothing pushed), so the pitch says "configured" unless a green run exists. Fill in the [Team: ...] placeholders in QA.md (maintainer, who reviewed what) and record a demo video for pitch/backup/.
7. A real race against PAL2NAL (`python demo/race.py`) must be tried on a laptop with PAL2NAL fetched; here only the "not installed" path and a stand-in were tested.

## How to try it
Linux / macOS:
```
pip install -e ".[gui]"
python -m coal_app                      # GUI on http://127.0.0.1:7860 ; --demo shows the sample report
python -m coal_app --check              # which engine features were detected
python -m coal_app.viewer examples/one_broken_gene/protein.aln.fasta --dna examples/one_broken_gene/dna.fasta -o view.html
python -m coal examples/mito/protein.aln.fasta examples/mito/dna.fasta -codontable 2
pal2nal examples/clean/protein.aln.fasta examples/clean/dna.fasta -output fasta
python demo/race.py                     # needs: python scripts/fetch_pal2nal.py
open pitch/explainer.html               # (Linux: xdg-open)
```
Windows (PowerShell):
```
py -m pip install -e ".[gui]"
py -m coal_app
py -m coal_app.viewer examples\one_broken_gene\protein.aln.fasta --dna examples\one_broken_gene\dna.fasta -o view.html
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
- Sat 02:29 Task 11 done: core.methods_text(report) (report['methods_text'] if the engine writes one, else built from version, genetic-code name + number, counts, -nogap/-nomismatch and the command with file names only; empty when nothing was aligned). Shown in the GUI as a read-only textbox with a copy button (Gradio show_copy_button, unverified), returned by codon_align as methods_text and appended to explain_report. Also fixed: rebuilt reports claimed --keep-going in the command when the engine doesn't have it. 4 new tests.
- Sat 02:30 Task 12 done: demo/mcp_demo.md (setup with `claude mcp add --transport stdio coal -- python -m coal_app --mcp`, syntax checked against code.claude.com/docs/en/mcp; Claude Desktop config; 60-second timed script; real outputs of codon_align / explain_report on one_broken_gene and mito with code 1 then 2), regenerated by demo/make_mcp_demo.py. Making it showed `problems` listed 19 'not checked' rows: now only warnings/failures (count stays in summary.not_checked). On the current engine the script says COAL stopped at OG0001_sp14 (not '19 of 20 fine'); it notes how the answer changes with --keep-going.
- Sat 02:30 Final: engine tests 11 OK; app tests 52 OK, 2 skipped. Summary updated. Stopping.
- Sat 02:36 Task 13 done: pitch/PITCH.md (3 minutes: explainer 0:00-0:30, problems 0:30-0:55, GUI demo 0:55-2:05 with clicks and words to the second, proof 2:05-2:35, vision 2:35-3:00; [after merge] labels; backup plan), pitch/QA.md (10 questions incl. the 7 required + data privacy, AI agents, what's missing), pitch/backup/ (4 PNGs from the real viewer and explainer, headless Chromium). Facts checked: counted the 165 differential checks by kind (143 identical, 12 both-fail, 10 deliberate); speed numbers copied from README (measured, Windows); engine ~1,045 lines, app ~1,770. Hajk Drost mentioned only generically ("the kind of questions Hajk Drost asks"); team should confirm the wording. Small viewer fix found while making screenshots: the 19 'not checked' rows are now one summary line in the problem list. 
- Sat 02:36 Final: engine tests 11 OK; app tests 52 OK, 2 skipped. Stopping.
