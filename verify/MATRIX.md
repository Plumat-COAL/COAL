# Coverage matrix (phase 0)

Every option, tool parameter, control, installer option and numeric or factual claim the verification covers, the phase that tests it and its result: PASS, FINDING V-n, or SKIPPED (reason). Rows are grouped by area; `python verify/matrix.py` rewrites this file from `verify/matrix.py`.

**174 of 174 rows filled.**

| id | area | item | phase | result |
|---|---|---|---|---|
| C01 | coal CLI | pep_aln positional (FASTA or CLUSTAL) | 1,2 | PASS |
| C02 | coal CLI | nuc_fasta positional, one or more files | 1,2 | PASS |
| C03 | coal CLI | `-` reads standard input (protein or one DNA file; two is a usage error) | 2,4 | PASS |
| C04 | coal CLI | -codontable N for all 17 codes (1-6, 9-16, 21-23); an unknown code is a usage error | 2,4 | PASS |
| C05 | coal CLI | -o FILE; `-o -` means stdout; -o equal to an input path refused | 2,4 | PASS |
| C06 | coal CLI | -q silences warnings and notices, not errors; stdout unchanged | 4 | PASS |
| C07 | coal CLI | -strict / --strict: the three rules (unrelated, guessed position, partial IDs) -> exit 1 | 4 | PASS |
| C08 | coal CLI | -nogap: removes gap and stop columns of the chosen code | 2 | PASS |
| C09 | coal CLI | -nomismatch: removes mismatched columns, prints no mismatch warnings | 2,4 | PASS |
| C10 | coal CLI | -output fasta (default) | 2 | PASS |
| C11 | coal CLI | -output clustal | 2 | PASS |
| C12 | coal CLI | -output paml | 2 | PASS |
| C13 | coal CLI | -output codon; refused with -nogap/-nomismatch/-blockonly | 2,4 | PASS |
| C14 | coal CLI | -blockonly with '#' marks of a CLUSTAL input | 2 | PASS |
| C15 | coal CLI | -j N: 1 / 2 / all identical in every format | 3 | PASS |
| C16 | coal CLI | --keep-going: exit 2, skipped genes named, others in order, all formats, -j > 1 | 4 | PASS |
| C17 | coal CLI | --progress: only for inputs over 10,000 sequences, never on stdout | 4 | PASS |
| C18 | coal CLI | --resume: kill and resume byte-identical, every format; 'inputs changed' refusal | 4 | PASS |
| C19 | coal CLI | --max-warnings N (default 20, 0 = all) and the '... and N more warnings' line | 4 | PASS |
| C20 | coal CLI | --report FILE.json on success, warnings, failure and --keep-going; every schema field | 4 | PASS |
| C21 | coal CLI | --json to stdout (needs -o) | 4 | PASS |
| C22 | coal CLI | --batch JOBS: one coal command line per line, byte-identical to single runs, exit codes | 4 | PASS |
| C23 | coal CLI | --batch-out-dir DIR | 4 | PASS |
| C24 | coal CLI | --version prints `coal 0.5.0` | 9 | PASS |
| C25 | coal CLI | -h / --help text | 0 | PASS |
| C26 | coal CLI | exit codes: 0 / 1 / 2 / 130 Ctrl-C / 143 SIGTERM / 141 closed pipe silent / 129 SIGHUP; no temp files left | 4 | PASS (SIGKILL of a stdout run leaves .coal_*.tmp in the temp folder: V-3, docs) |
| C27 | coal CLI | safe ID pairing: normalised matches, order fallback warning, no false pairing when ambiguous | 4 | PASS |
| C28 | coal CLI | summary and hint lines: the mito example suggests code 2; nothing suggests a worse code | 4 | PASS |
| C29 | coal CLI | methods_text: version, plumat-coal URL (no 'chloe-newt' anywhere), code name, counts, true statements | 4 | PASS |
| C30 | coal CLI | automatic hand-over to coal-cpp for inputs >= 1 MB; fallback / message when coal-cpp is missing | 3 | PASS |
| C31 | coal CLI | --cpp-path | 3,9 | PASS |
| X01 | coal-cpp CLI | help text identical to Python's except the documented -j wording (threads vs processes) | 0 | PASS |
| X02 | coal-cpp CLI | every option of coal except --batch / --batch-out-dir (refused with a clear error) | 2,3,4 | PASS |
| X03 | coal-cpp CLI | -j N threads: output identical for 1 / 2 / all, every format | 3 | PASS |
| X04 | coal-cpp CLI | --version prints `coal 0.5.0 (C++)` | 9 | PASS |
| X05 | coal-cpp CLI | stdout, stderr, exit status and report identical to the Python engine on every phase-2 case | 3 | PASS |
| X06 | coal-cpp CLI | --resume of its own runs; refuses the Python engine's checkpoint with the documented error | 4 | PASS |
| X07 | coal-cpp CLI | signals and exit codes as the Python engine; no temp files left | 4 | PASS |
| P01 | pal2nal drop-in | no -output -> CLUSTAL by default | 5 | PASS |
| P02 | pal2nal drop-in | -output clustal / paml / fasta / codon | 5 | PASS |
| P03 | pal2nal drop-in | -codontable N; `-codontable 2x` means 2; a dangling -codontable or -output is ignored | 5 | PASS |
| P04 | pal2nal drop-in | -nogap | 5 | PASS |
| P05 | pal2nal drop-in | -nomismatch | 5 | PASS |
| P06 | pal2nal drop-in | -blockonly | 5 | PASS |
| P07 | pal2nal drop-in | -nostderr maps to -q (errors still printed) | 5 | PASS |
| P08 | pal2nal drop-in | -html: clear error, non-zero exit | 5 | PASS |
| P09 | pal2nal drop-in | -h / no arguments: usage text; `pal2nal | head` exits 141 | 5 | FINDING V-2 (no arguments: usage on stdout, exit 0; -h fine; pal2nal | head exits 141) |
| P10 | pal2nal drop-in | COAL-only options (--strict, -j, -o ...) refused with a clear error, non-zero exit | 5 | PASS |
| P11 | pal2nal drop-in | `-` standard input | 5 | PASS |
| P12 | pal2nal drop-in | exit codes: 1 for usage error / bad input (not 2), 130 / 141 / 143 | 5 | FINDING V-2 (a call without arguments exits 0; 1 / 130 / 141 / 143 otherwise as documented) |
| P13 | pal2nal drop-in | entry point works from a fresh shell (venv bin) | 5 | PASS |
| P14 | pal2nal drop-in | output equals PAL2NAL v14's for every option combination (fuzz, new seeds) | 2,5 | PASS (every raw difference is a documented one); FINDING V-1 in the fuzz harness (two seeds) |
| M01 | MCP | `python -m coal_app --mcp` starts, initialize / tools/list / tools/call over stdio, clean shutdown | 6 | PASS |
| M02 | MCP | server never writes to stdout outside the protocol | 6 | PASS |
| M03 | MCP | codon_align: protein_alignment_path (absolute; relative, spaces, missing, wrong type) | 6 | PASS |
| M04 | MCP | codon_align: dna_fasta_path (same cases) | 6 | PASS |
| M05 | MCP | codon_align: output_format fasta / clustal / paml / codon; an unsupported format is a tool error | 6 | PASS |
| M06 | MCP | codon_align: codontable (valid codes; an unknown code is a tool error) | 6 | PASS |
| M07 | MCP | codon_align: nogap | 6 | PASS |
| M08 | MCP | codon_align: nomismatch | 6 | PASS |
| M09 | MCP | codon_align: strict (refused genes left out and reported, run goes on) | 6 | PASS |
| M10 | MCP | codon_align: blockonly | 6 | PASS |
| M11 | MCP | codon_align: save_as (new file, existing file never overwritten, existing folder, missing folder) | 6 | PASS |
| M12 | MCP | codon_align result: completed / all_ok / needs_attention / counts / suggested_codontable / next_step / problems / output_path / report_path right on each example | 6 | PASS |
| M13 | MCP | codon_align on a large input (>= 1 MB: coal-cpp hand-over through the app layer) | 3,6 | PASS |
| M14 | MCP | explain_report(report_path): correct and plain for ok / warning / failed / aborted reports; unreadable report is a tool error | 6 | PASS |
| M15 | MCP | tool docstrings and parameter descriptions accurate and enough for an agent | 6 | PASS |
| M16 | MCP | each example dataset through codon_align | 6 | PASS |
| M17 | MCP | `COAL.exe -m coal_app --mcp` (installed Windows app) | 8 | PASS (initialize, tools/list = codon_align + explain_report, codon_align on the clean example and on 13.7 MB (coal-cpp), explain_report, clean stdout, no traceback) |
| G01 | GUI | Align tab: Protein alignment file box (upload) and path textbox | 7 | PASS |
| G02 | GUI | Align tab: Coding DNA file box and path textbox | 7 | PASS |
| G03 | GUI | Try an example buttons: clean, mito (wrong code), messy_ids, one_broken_gene | 7 | PASS |
| G04 | GUI | Sample report button | 7 | PASS |
| G05 | GUI | Genetic code dropdown (17 codes) | 7 | PASS |
| G06 | GUI | Output format dropdown (4 formats) | 7 | PASS |
| G07 | GUI | Remove gap columns checkbox | 7 | PASS |
| G08 | GUI | Remove mismatched columns checkbox | 7 | PASS |
| G09 | GUI | Refuse probably-wrong genes checkbox | 7 | PASS |
| G10 | GUI | Keep only '#'-marked columns checkbox (CLUSTAL input) | 7 | PASS |
| G11 | GUI | Start alignment button; banner says whether everything matched | 7 | PASS |
| G12 | GUI | Hint banner's one-click Re-run with the suggested code | 7 | PASS |
| G13 | GUI | Download codon alignment: byte-identical to the CLI in every format | 7 | PASS |
| G14 | GUI | Download report (JSON) | 7 | PASS |
| G15 | GUI | Health table: problem genes first, symbol + word, clicking a warning jumps to the right column | 7 | PASS |
| G16 | GUI | Alignment view: centre on column + Show these columns; Download alignment view (HTML) | 7 | PASS |
| G17 | GUI | Methods paragraph with copy button | 7 | PASS; V-11 (the label is not linked to the box) |
| G18 | GUI | Raw COAL messages box | 7 | PASS |
| G19 | GUI | Settings tab: code, format, three checkboxes persist across restarts; Reset to defaults; path shown | 7 | PASS (persist across page loads; Reset to defaults; settings.json path shown; a server restart was not exercised) |
| G20 | GUI | errors for bad uploads are plain and recoverable | 7 | FINDING V-12 (an unreadable file gives 'COAL stopped: no alignment was written' with the reason only in the raw messages box); recoverable |
| G21 | GUI | light and dark mode | 7 | PASS |
| G22 | GUI | scale: 100,000-sequence input (time to first result, responsive, honest preview) | 7 | PASS (10.7 s to the result, page responsive, honest preview); V-10 (tip wording) |
| G23 | GUI | offline: no external requests (fonts, CDNs, analytics) | 7 | FINDING V-8 (iframeResizer from cdnjs on every page load) |
| G24 | GUI | accessibility: keyboard-only main flow, visible focus, axe-core, contrast >= 4.5:1, status never by colour alone | 7 | FINDING V-9 (axe aria-hidden-focus, Gradio's hidden tab bar); keyboard flow, focus and contrast PASS; V-11 |
| G25 | GUI | screenshots at 1366x768, 1920x1080, 200 % zoom | 7 | PASS |
| G26 | GUI | usability review against docs/ui-principles.md: clean run, wrong code -> fix, one broken gene | 7 | PASS (review in section 7, phase 7) |
| G27 | GUI | the GUI also serves the MCP tools at /gradio_api/mcp/sse | 6,7 | PASS (GET /gradio_api/mcp/sse answers 200; /gradio_api/mcp/schema lists codon_align and explain_report) |
| W01 | Windows | download latest release; SHA256SUMS verified | 8 | PASS (gh release download v0.5.0; both files match SHA256SUMS.txt, re-hashed in the Windows pass) |
| W02 | Windows | per-user install without admin; SmartScreen text recorded; site/README tell users what to click | 8 | PASS (per-user, no UAC; SmartScreen: Windows protected your PC / Unknown publisher, recorded by the owner; desktop README and site/install.html say More info then Run anyway, as seen; Smart App Control is off on this PC: human checklist) |
| W03 | Windows | installer task: desktop shortcut (unchecked by default) | 8 | PASS (a default install makes no desktop shortcut; /TASKS=desktopicon makes COAL.lnk; uninstall removes it) |
| W04 | Windows | installer task: 'Open with' for .aln/.fasta/.fa/.fas/.faa/.fna (unchecked by default; defaults unchanged) | 8 | PASS (off by default; the task registers COAL.SequenceFile under OpenWithProgids for .aln .fasta .fa .fas .faa .fna with the command COAL.exe %1; .fasta default program and UserChoice unchanged; uninstall removes the entries and the ProgId) |
| W05 | Windows | installer run page: Connect COAL to Claude Desktop (claude_desktop_config.json entry, other servers kept, .bak-coal backup) | 8 | PASS (--connect-claude adds coal = COAL.exe -m coal_app --mcp, keeps another server and other keys, writes claude_desktop_config.json.bak-coal byte for byte, is idempotent; on the first-ever install there was no file to back up, as expected) |
| W06 | Windows | installer run page: Open COAL now; first launch | 8 | PASS (window started by the owner at install time; here COAL.exe window: local app ready 3.7 s after start, 2 COAL.exe processes, killing the window leaves none) |
| W07 | Windows | each example in the desktop app; drag-and-drop; path with spaces; OneDrive-style folder; long paths | 8 | PASS for the four examples, a path with spaces in a OneDrive-style folder, 13.7 MB through the path boxes; drag-and-drop: human checklist; long paths: not testable (LongPathsEnabled is off on this PC: the harness itself could not create a 270-character path) |
| W08 | Windows | high-DPI and dark mode | 8 | PASS for dark mode (the app page in dark theme); high-DPI: human checklist |
| W09 | Windows | uninstall removes app, shortcut, Open-with entries, Claude entry; restores nothing it should not | 8 | PASS (silent uninstall removes the app folder, Start Menu and desktop shortcuts, the Add/Remove entry, the Open-with entries and ProgId, LOCALAPPDATA/COAL and TEMP/coal_app_runs, and the coal entry in the Claude config while keeping other servers and keys; leaves TEMP/gradio as documented; never restores .bak-coal over later edits) |
| W10 | Windows | portable zip runs | 8 | PASS (zip extracted: COAL.exe --version, -m coal --cpp-path names its own coal-cpp.exe, 13.7 MB run in coal-cpp, --serve ready in 6.0 s and the four examples) |
| W11 | Windows | bundled coal-cpp.exe is the one used for >= 1 MB inputs | 8 | PASS (installed and portable: --cpp-path is the bundled _internal/coal-cpp.exe, 0.5.0 (C++); a 13.7 MB run through -m coal, through MCP and through the path boxes reports engine coal-cpp) |
| W12 | Windows | COAL.exe --connect-claude / --disconnect-claude / --mcp | 8 | PASS (--connect-claude / --disconnect-claude with --quiet exit 0 and edit only the coal entry; -m coal_app --mcp serves both tools) |
| D01 | Distribution | release wheel installs in clean venvs on Python 3.10 and 3.13 | 9 | PASS |
| D02 | Distribution | sdist builds and installs | 9 | PASS |
| D03 | Distribution | `pip install git+https://...` works | 9 | PASS |
| D04 | Distribution | version strings agree: pyproject (both), coal --version, coal-cpp --version, report, CHANGELOG, site, release tag, CITATION.cff, SKILL.md, desktop README | 9 | FINDING V-5 (desktop README names COAL-Setup-0.4.2.exe; every other string is 0.5.0) |
| D05 | Distribution | platform wheels with coal-cpp inside (release assets) install and work | 9 | PASS |
| R01 | README | every command example runs with the shown output and exit code | 10 | PASS |
| R02 | README | '300x faster on 1 million sequences (12.5 s vs 62 min, 172 MB vs 18.6 GB)' traces to benchmarks/results | 10,11 | PASS (traces to benchmarks/results/windows_47f57f1; this machine's numbers: phase 11) |
| R03 | README | 'over 1,000x faster on long genes (1,000 x 10,000 codons: 0.25 s vs 280 s)' traces to a committed result | 10,11 | PASS (traces to benchmarks/results/windows_47f57f1; this machine's numbers: phase 11) |
| R04 | README | 'PAL2NAL v14 cannot process proteins longer than about 22,000 residues; COAL can' | 10 | PASS (adversarial very_long_gene, 1,000,000 codons: PAL2NAL v14 did not finish in 600 s, both engines aligned it identically; the 22,000 figure itself was not re-measured) |
| R05 | README | '0 unexplained differences in random comparisons' reproducible | 2,10 | FINDING V-1 (true for the documented seeds; fresh seeds hit the harness gap, COAL right) |
| R06 | README | exit-code list (0/1/2/141/130/143/129) matches behaviour | 4,10 | PASS |
| R07 | README | --batch 'about 40x less wall time on tiny jobs, byte-identical' | 10 | PASS (byte-identical, exit codes; the 40x timing not re-measured) |
| R08 | README | install commands (release sdist URLs, source build, app install line) work | 9,10 | PASS |
| R09 | README | codon_align Python example gives the shown rows | 10 | PASS |
| R10 | README | 'Python 3.10+', 'no dependencies' for the engine | 9,10 | PASS |
| S01 | SKILL.md | '238 byte-for-byte comparisons', '285 checks', '20,000 random comparisons: 13,826 / 4,681 / 1,493 / 0' reproducible | 1,2,10 | PASS (285 / 238 / 24 / 23 / 39 / 4 reproduced in phase 1; default-seed fuzz clean, both engines) |
| S02 | SKILL.md | install lines and `python -m coal --version` check | 10 | PASS |
| S03 | SKILL.md | common commands 1-7 run as described | 10 | PASS |
| S04 | SKILL.md | notice texts quoted match the engine's | 4,10 | PASS |
| S05 | SKILL.md | exit-code list incl. 'On Windows (expected, not yet verified there)' | 4,8,10 | PASS (Windows, both engines: test_decisions and test_robustness exit-code tests pass: 130 / 141 / 143 as documented; so the SKILL.md wording not yet verified there is stale: V-6) |
| S06 | SKILL.md | GUI on 7860 'serves MCP there too'; `--mcp` gives two tools; absolute paths; result fields | 6,7,10 | PASS |
| S07 | SKILL.md | report field descriptions (summary, sequences, hints, output, command) match the schema | 4,10 | PASS |
| S08 | SKILL.md | 'After the engine update (only if --help lists these)' section: --keep-going, --report, normalised pairing are in the current engine | 10 | FINDING V-6 (stale 'engine update' section) |
| S09 | SKILL.md | failure table rows (messages quoted) match the engine's wording | 4,10 | PASS |
| H01 | CHANGELOG | numbers in 'Unreleased' and 0.5.0 trace to committed results or reproducible commands | 10,11 | PASS (commands and behaviours; the memory and timing numbers: B01 / B02) |
| H02 | CHANGELOG | exit codes, options and behaviours it describes exist | 4,10 | PASS |
| E01 | desktop/README | 'Download COAL-Setup-0.4.2.exe' vs the current release | 8,9,10 | FINDING V-5 |
| E02 | desktop/README | SmartScreen steps, 'Publisher: Unknown publisher', Smart App Control note | 8,10 | PASS (the SmartScreen steps and Unknown publisher match what the installer showed; the Smart App Control note could not be exercised, SAC is off here) |
| E03 | desktop/README | install folder %LOCALAPPDATA%\Programs\COAL, 127.0.0.1 only, no admin | 8 | PASS (LOCALAPPDATA/Programs/COAL; --serve binds 127.0.0.1 only; no UAC, install log: Administrative install mode: No) |
| E04 | desktop/README | use steps (drop files, genetic code, Start alignment, examples), results (banner, table, tip, view, downloads, methods, commands) | 7,8 | PASS for the GUI steps (phase 7); the desktop app itself: Windows pass |
| E05 | desktop/README | uninstall removes what it says and leaves %TEMP%\gradio | 8 | PASS (uninstall leaves TEMP/gradio and removes TEMP/coal_app_runs and LOCALAPPDATA/COAL as the README says) |
| E06 | desktop/README | run from source: run_from_source.bat (.venv-gui), setup_mcp.bat (.venv-mcp), Python launcher order | 8,10 | PASS by reading (run_from_source.bat makes .venv-gui with [gui] + pywebview, setup_mcp.bat makes .venv-mcp with [mcp], both try py -3.12, py -3, python); not executed (they create environments in the checkout) |
| E07 | desktop/README | privacy: Connect COAL to Claude Desktop 'ticked by default', --connect-claude, --mcp, logs folder | 8,10 | PASS (Connect COAL to Claude Desktop is ticked by default in installer.iss and was in the wizard; --connect-claude and --mcp work; logs in LOCALAPPDATA/COAL/logs as documented) |
| E08 | desktop/README | 'first start can take ten seconds or so' | 8 | PASS (first start here: window to local app ready in 3.7 s; the documented ten seconds is an upper bound) |
| T01 | site | index: version 0.5.0, download links resolve (installer, zip, sdists, SHA256SUMS), checksums match | 10 | PASS |
| T02 | site | index/differences: '285 checks', '238', '24', '23 (6/3/14)', '39 input cases', '4 output formats' | 1,10 | PASS |
| T03 | site | differences timing table (74x / 223x / 297x / 1,115x; 12.5 s; 172 MB; 62 min; 18.6 GB; 'Python engine under a minute at a million genes, 74x') | 10,11 | PASS (traces to benchmarks/results/windows_47f57f1; this machine's numbers: phase 11) |
| T04 | site | guide: clean example has 8 genes; mito 24 genes, 179 warnings, summary text; one_broken_gene column 102, 6 not checked, 19 written with --keep-going | 10 | PASS |
| T05 | site | guide: '4 million bases of DNA per extra process', '-j 2', '1 MB or more' hand-over, build instructions | 10 | PASS |
| T06 | site | install: pip commands for Windows / macOS / Linux, http://127.0.0.1:7860, coal-gui, checksum commands, 'first start about 10 seconds' | 10 | PASS |
| T07 | site | install: Smart App Control text; file table (installer, zip, sdists, SHA256SUMS) | 8,10 | PASS for the file table (installer, zip, sdists, SHA256SUMS all release assets with matching checksums); the Smart App Control text: human checklist (SAC off on this PC) |
| T08 | site | cite: version 0.5.0, CITATION.cff agrees | 9,10 | PASS |
| T09 | site | no broken links; 'Last updated 23 October 2026' | 10 | PASS |
| T10 | site | explainer page renders (Chromium), no external requests | 10 | PASS |
| B01 | docs/benchmarks | 47f57f1 table reproducible in shape (phase 11 numbers compared); raw data committed | 11 | PASS (shape reproduced on ext4; numbers in section 3) |
| B02 | docs | docs/usage.md memory figures (350 B / 190 B per sequence; 150 MB per worker; 100,000 unrelated genes 118 / 69 / 885 MB) against phase 11 | 11 | PASS (500,000 x 150: Python 190 / 191 / 190 / 197 MB, coal-cpp 112 / 137 / 113 / 149 MB for fasta / clustal / paml / codon, identical outputs; usage.md figures met) |
| B03 | docs | every command example in docs/*.md runs as shown | 10 | PASS |
| B04 | docs | pitch/QA.md: listed as a source by the prompt but absent from the repository | 10 | SKIPPED (pitch/QA.md is not in the repository) |
| N01 | claims absent | NVIDIA Brev, hosted web version, PyPI availability, code signing, C++ via pip, DIAMOND on real genomes | 10 | PASS |
| Q01 | suites | python -m unittest discover -s tests: counts vs claimed (e.g. 285 golden checks); skipped tests listed | 1 | PASS (WSL: 502 OK, 8 skipped; Windows: 502 run, 88 skipped, 2 environment failures: V-14 symlink privilege, V-15 venv launcher count) |
| Q02 | suites | the same with COAL_IMPL=cpp | 1 | PASS (WSL: 502 OK, 8 skipped; Windows through coal-cpp.exe (MSVC): 502 run, 88 skipped, the same 2 environment failures) |
| Q03 | suites | tests/test_cpp_vs_python.py | 1 | PASS |
| Q04 | suites | apps/coal-app/tests (gui venv and mcp venv) | 1 | PASS |
| Q05 | suites | desktop tests (test_wrapper, test_smoke) | 1 | PASS |
| Q06 | suites | PAL2NAL fuzz harness (default seed), both engines | 1 | PASS (default seed: 0 unexplained, both engines) |
| Q07 | suites | COAL_REF_SRC and similar env vars default to the repo's own src | 1 | PASS |
| Q08 | suites | CI status of the latest main commit | 1 | PASS |
| F01 | correctness | fuzz with 3 new seeds x FUZZ_N >= 2,000, all 4 formats, -nogap, -nomismatch, all 17 codes: every difference maps to a DELIBERATE row | 2 | PASS (0 wrong alignments in 12,000 cases); FINDING V-1 in the harness |
| F02 | correctness | adversarial generator: every case of the prompt's list through PAL2NAL, Python and C++ | 2 | PASS |
| F03 | correctness | real data: published alignments rebuilt in all 4 formats with both engines; '17 rebuilt exactly' claim | 2 | PASS (12 datasets rebuilt exactly, both engines); V-7: the '17' figure is local-only |
| K01 | perf | benchmarks/run_bench.py at 10,000x500, 100x10,000, 100,000x500 (Python vs C++, -j 1 and default), PAL2NAL at the two smaller | 11 | PASS (shape reproduced on ext4; numbers in section 3) |
| K02 | perf | peak memory at 500,000 x 150 codons in every format, both engines (streaming claim) | 11 | PASS (500,000 x 150: Python 190 / 191 / 190 / 197 MB, coal-cpp 112 / 137 / 113 / 149 MB for fasta / clustal / paml / codon, identical outputs; usage.md figures met) |
| J01 | demos | demo/race.py against the real PAL2NAL | 12 | SKIPPED (demo/race.py removed with the hackathon material in dfe0049) |
| J02 | demos | demo/diamond_pipeline.py (stand-in test only without diamond) | 12 | PASS |
| J03 | demos | demo/mcp_demo.md tool outputs still match the tools | 12 | FINDING V-4 (demo/mcp_demo.md shows the old stop-at-first-error output) |
| U01 | mutations | about 15 realistic bugs in a scratch copy; each not caught is a test gap | 13 | DONE: 15 bugs, 12 caught, 2 equivalent mutants (Python start-codon set is a fast path only; -blockonly keep_col), 1 test gap: V-16 (GUI download button file untested); the C++ reverse-order mutant caught by a hang, not an assertion |
