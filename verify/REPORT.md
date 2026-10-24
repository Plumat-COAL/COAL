# COAL full verification, 2026-10-24 (Linux pass under WSL 2, then a native Windows pass for phase 8)

Branch `verify/full-2026-10-24` of a fresh clone of `main` at `4159e10` (merge of PR #43, the memory improvements);
latest release tag `v0.5.0`; every version string in the repository is `0.5.0` (the one exception is a finding below).
Run unattended by Claude Code following `full-test-prompt.md`; the decisions it took are in `verify/LOG.md`.
**Update, later on 2026-10-24:** every finding V-1 to V-16 was fixed and merged into `main` as PR #44 (3d525ce);
section 4 records a second pass at the human checklist, which added the 1,000,000-sequence benchmark and V-17.

## 1. Verdict per area

| Area | Verdict | Notes |
|---|---|---|
| Engine, Python | PASS | every suite, 3 new fuzz seeds, 50 adversarial cases, flags, signals, resume: no wrong output found |
| Engine, C++ | PASS | byte-identical to the Python engine on everything above; the only help-text difference is the documented `-j` wording |
| `pal2nal` drop-in | PASS (1 minor) | V-2: no-argument usage goes to stdout with exit 0 |
| Report and flags | PASS | schema valid, checksums right, counts consistent, `--resume` identical after a kill in every format, both engines |
| MCP server | PASS | both tools, every parameter, clean stdout, `save_as` rules, large-input hand-over |
| GUI | PASS (6 usability / minor) | every control and flow right, downloads byte-identical to the CLI; V-8 CDN script vs the offline claim, V-9 axe (Gradio), V-10 tip wording, V-11 unlabelled box, V-12 error reason hidden (all fixed in PR #44); V-17 the residue palette merges for colour-blind readers (second pass, open) |
| Windows app and installer | PASS | installer, uninstaller, both optional tasks, Claude Desktop connect/disconnect, installed and portable apps, bundled `coal-cpp.exe`, MCP via `COAL.exe`: all as documented; long paths, drag-and-drop, high-DPI, Smart App Control and a real Claude Desktop chat are on the human checklist |
| Tests and harness | PASS (5 minor) | 502 tests pass on WSL both engines; on Windows 2 environment failures (V-14 symlink privilege, V-15 venv launcher count); V-1 fuzz neutralisation gap, V-13 explainer test skipped by its Linux-only lookup; 15 planted bugs: 12 caught, 2 equivalent mutants, 1 gap (V-16, the GUI download file) |
| Distribution | PASS (1 docs) | wheels, sdists, git install, platform wheel with `coal-cpp` inside; V-5 stale installer name in the desktop README |
| Docs and site | ISSUES (docs only) | V-4 stale MCP demo output, V-5 stale installer name, V-6 stale skill sections, V-7 unreproducible '17', V-3 SIGKILL note; every command example and number checked traces |
| Performance | PASS | benchmarks re-run on this machine keep the documented shape; memory flat at 500,000 sequences in every format, under the documented figures |

## 2. BLOCKER and MAJOR findings

None found. (V-1 is a gap in the PAL2NAL fuzz harness, not a wrong alignment; see the findings.)

## 3. Numbers

Everything in this section was measured on this machine during the run (WSL 2 Ubuntu 26.04, 16 CPUs, 30 GB, the engines built with `-O2`, Python 3.13); medians of 3 timed runs after a warm-up (`benchmarks/run_bench.py`). "DrvFs" rows had the data on the Windows drive mounted into WSL (slow file I/O; kept to show the effect), "ext4" rows had the data on the Linux file system, which is the comparable setting. The raw CSVs are in `verify/p11/`.

**Engines (`--impl both`).**

| Data | Workload | Python, 1 process | Python, all cores | coal-cpp, 1 thread | coal-cpp, all threads | identical |
|---|---|---|---|---|---|---|
| DrvFs | 10,000 × 500 | 2.63 s | 2.31 s | 0.812 s | 0.807 s | True |
| DrvFs | 100 × 10,000 | 0.639 s | 0.655 s | 0.0615 s | 0.0615 s | True |
| DrvFs | 100,000 × 500 | 23.6 s | 22.7 s | 8.43 s | 8.45 s | True |
| ext4 | 10,000 × 500 | 0.595 s | 0.44 s | 0.108 s | 0.0779 s | True |
| ext4 | 100 × 10,000 | 0.21 s | 0.231 s | 0.0311 s | 0.0262 s | True |
| ext4 | 100,000 × 500 | 4.18 s | 2.19 s | 0.898 s | 0.749 s | True |
| ext4 | 1,000,000 × 400 | 36.7 s | 17.8 s | 8.44 s | 7.53 s | True |

**Against PAL2NAL v14 (Perl 5, same inputs).**

| Data | Workload | PAL2NAL v14 (peak memory) | Python engine, 1 process (peak) | speed-up | Python, all cores | speed-up | identical |
|---|---|---|---|---|---|---|---|
| DrvFs | 10,000 × 500 | 11.7 s (0.187 GB) | 2.61 s (44.5 MB) | 4.48× | 2.63 s | 4.44× | True |
| DrvFs | 100 × 10,000 | 12.9 s (0.0417 GB) | 0.628 s (21.9 MB) | 20.6× | 0.598 s | 21.6× | True |
| ext4 | 10,000 × 500 | 11.3 s (0.187 GB) | 0.581 s (44.4 MB) | 19.5× | 0.431 s | 26.2× | True |
| ext4 | 100 × 10,000 | 12.8 s (0.0417 GB) | 0.215 s (21.4 MB) | 59.5× | 0.215 s | 59.4× | True |
| ext4 | 1,000,000 × 400 | 1.04e+03 s (14.6 GB) | 39.1 s (345 MB) | 26.5× | – s | –× | True |

Compared with the committed tables: `docs/benchmarks.md` (v0.3 table, commit 58a1219) has 10,000 × 500 at 24.7 s PAL2NAL / 0.864 s COAL and 100 × 10,000 at 24.5 s / 0.207 s; here (ext4) 11.3 s / 0.58 s and 12.8 s / 0.215 s, so PAL2NAL is about twice as fast on this CPU and the Python engine's long-gene case is the same; the speed-ups keep their shape (sequences: ~20×, long genes: ~60×). The 47f57f1 Windows table has 100,000 × 400 at 1.07 s (coal-cpp, 1 thread) and 5.04 s (Python, 1 process); here 100,000 × 500 took 0.90 s and 4.13 s. The 1,000,000 × 400 row was run in the second pass (one timed run each, data on the WSL ext4 home, no all-core PAL2NAL comparison): PAL2NAL v14 1,035 s and 14.9 GB here against 3,714 s and 18.6 GB in the committed Windows table, coal-cpp 8.4 s (1 thread) against 12.5 s, the Python engine 36 s (345 MB) against 50.3 s (1.2 GB), so the speed-ups on this Linux CPU are about 120× (coal-cpp) and 27× (Python) rather than 297× and 74×: PAL2NAL's Perl runs much faster here than on the Windows machine of the table, COAL's shape and memory hold.

**Peak memory at 500,000 sequences × 150 codons (`tests/tools/stress_run.py`, summed PSS of the process tree, `-j 1`, output to a file).**

| Output | Python engine | wall | coal-cpp | wall | outputs identical |
|---|---|---|---|---|---|
| fasta | 190 MB | 52 s | 112 MB | 38 s | True |
| clustal | 191 MB | 54 s | 137 MB | 40 s | True |
| paml | 190 MB | 49 s | 113 MB | 40 s | True |
| codon | 196 MB | 62 s | 149 MB | 45 s | True |

`docs/usage.md` says 203 / 203 / 209 MB (Python, fasta / clustal / codon) and 112 / 136 / 148 MB (coal-cpp) for this shape, measured before PR #43; the figures above are at or under those, so the "about 350 / 190 bytes per sequence" statement still holds (the Python engine is now nearer 380 bytes per sequence including its 15 MB base, coal-cpp 220, both flat across formats).

**Other numbers.** Engine suite: 502 tests OK on both engines (8 skipped, each explained in phase 1); fuzz: 3 new seeds × 2,000 cases × 4 formats × both engines, 2 harness-level unexplained differences (V-1), 0 wrong alignments; adversarial: 50 cases, engines identical on every one; GUI: 100,000 sequences through the path boxes in 10.7 s; mutation testing: see phase 13.

## 4. Human checklist (what an unattended run could not verify), and the second pass at it

Each item names the phase it belongs to and why it needs a person. After the first report the owner asked for this
list to be worked through; the second pass (2026-10-24, later the same day) settled what a script could and records
the rest.

- **Done in the second pass.** (1) The 1,000,000-sequence benchmark (phase 11): PAL2NAL v14 and both engines on this
  machine, numbers in section 3. (2) Colour (phase 7): the residue-class palette simulated for the four common
  deficiencies, with Chromium's own emulation and numerically: finding V-17 (the classes merge; the letters and
  outlines carry the meaning). (3) The exit codes on Windows (SKILL.md's "not yet verified there"): verified by the
  Windows pass, docs updated in PR #44.
- **Could not be done by the agent, still for a person.**
  - *A real Claude Desktop chat (phase 8).* The agent runs inside Claude Desktop and may not drive that application.
    Steps: restart Claude Desktop, open a new chat, check that the tools menu lists `coal`, and ask "Build a codon
    alignment from `D:\coal-handoff-kit\COAL\examples\clean\protein.aln.fasta` and `...\dna.fasta`"; the answer
    should name the output file and say all 8 genes match (`COAL.exe -m coal_app --mcp` answers the protocol, and the
    config entry is in place).
  - *Drag-and-drop onto the COAL window, the "Open with" entry in Explorer's menu, the wizard's and the uninstaller's
    own pages (phase 8).* The owner declined desktop control for Explorer and COAL in the second pass, so these stay by
    hand: drop two files from Explorer onto the window (it should fill both boxes); right-click a `.fasta` file and look
    for COAL under "Open with" (only after an install with that task ticked; the registry command is right); run
    `COAL-Setup-0.5.0.exe` once more and read each page (the silent pass checked what each page does).
  - *High-DPI (phase 8).* Changing the display scaling is a system setting the agent does not touch: set 150 % and
    200 % in Settings > System > Display and open COAL from the Start menu; the page itself was checked at 200 % zoom.
  - *Long paths (phase 8).* `LongPathsEnabled` is off on this PC (a registry key under HKLM the agent does not set):
    turn it on, then align `examples/clean` from a folder whose path is longer than 260 characters through the path
    boxes.
  - *Smart App Control (phase 8, E02 / T07).* Off on this PC and only switchable on by reinstalling Windows: the docs'
    text for it cannot be checked here.
  - *Screen reader (phase 7, G24).* axe-core and the keyboard flow were run; a read-through with NVDA or Narrator was
    not (after PR #44 the methods box is named; the hidden duplicate tabs are gone from the tab order).
  - *PAL2NAL v16 (phase 2).* The repository compares against v14 only; v16 is not committed and was not fetched.
  - *Root-only tests (phase 1).* Two tmpfs tests need root in WSL (no passwordless sudo on this account).
  - *`demo/race.py` (phase 12).* Removed from the repository; nothing to run.

## 5. Suggested fix order

V-1 to V-16 below are fixed in PR #44 (merged into `main`); V-17, found in the second pass, is open.

0. **V-17** (GUI): a residue-class palette that also varies lightness, or a switch to turn the residue colours off.

1. **V-12** (GUI): show the engine's reason in the banner when a run stops before pairing (one line from the raw
   messages); the message exists, it is just hidden.
2. **V-8** (GUI / docs): drop or serve locally Gradio's iframeResizer script, or soften "nothing leaves your computer"
   to "none of your data leaves your computer".
3. **V-2** (drop-in): `pal2nal` without arguments should print usage to stderr and exit 1 as its own help says.
4. **V-1** (fuzz harness): make the `-codontable 10 -nogap` neutralisation compare like for like, or exclude those
   cases from the "exactly the documented change" test; otherwise the nightly fuzz fails on roughly one new seed in
   three for a non-bug.
5. **V-4** (docs): regenerate `demo/mcp_demo.md`.
6. **V-5** (docs): `COAL-Setup-0.4.2.exe` in the desktop README.
7. **V-6** (docs): the stale "engine update" and "not yet verified on Windows" sections of SKILL.md (the Windows pass
   settles the second).
8. **V-10** (GUI): say "50 of the 100,000 genes were sampled" in the genetic-code tip.
9. **V-11** (GUI): link the "Methods paragraph" label to its textbox.
10. **V-9** (GUI): hide Gradio's duplicate tab bar from the accessibility tree (`display:none`).
11. **V-3** (docs): say that a killed (SIGKILL) run writing to stdout can leave `.coal_*.tmp` in the temp folder.
12. **V-7** (docs): either commit what makes "17 rebuilt exactly" reproducible or say "12 in the repository, 5 local".
13. **V-14, V-15** (tests, Windows): copy instead of symlink in `test_cpp_vs_reference`; count only multiprocessing
    workers in `test_resume_with_workers`.
14. **V-13** (tests): let `test_docs_explainer` find Playwright's Chromium on every platform.

## 6. Findings

Severity: BLOCKER (wrong scientific output, crash, data loss, overwrites user files, install/launch fails), MAJOR (wrong
exit code, misleading message, broken flag, broken flow), USABILITY, DOCS (a claim that does not match behaviour),
MINOR. Every reproducer is under `verify/repro/V-n/`.

### V-1 (MINOR, test harness; both engines) — the PAL2NAL fuzz reports an unexplained difference for `-codontable 10 -nogap` cases that are COAL being right

- **Area:** `tests/fuzz_vs_pal2nal.py` (the DELIBERATE neutralisation), docs claim "0 unexplained differences".
- **Steps:** `PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl FUZZ_N=2000 FUZZ_SEED=20261024 python tests/fuzz_vs_pal2nal.py`
  (case 1957) and `FUZZ_SEED=20261025` (case 1380); the same with `COAL_FUZZ_ENGINE=cpp`. Seed 31337 is clean.
- **Expected:** every difference maps to a DELIBERATE row (tests/README.md "Random differential testing": a
  difference is accepted only if it is exactly the documented one).
- **Actual:** one unexplained difference per seed, both with `-codontable 10 -nogap` (features `genetic_code` +
  `nogap_code`, case 1380 also `rna_first_met`). The checker's message: `without -nogap: s3 at nt 160: COAL GAN,
  expected gcN` / `3FhjC45zBMvXbk2u at nt 64: COAL GAG, expected gcG`. On the original input COAL's output equals
  PAL2NAL's for that gene; the checker compares COAL with PAL2NAL run on a *neutralised* copy (table-10 TGA codons
  rewritten to TGC), and in these genes the rewrite moves where PAL2NAL places the gene, so the "exactly the documented
  change" test fails. Where the raw outputs do differ (genes s4, s5 of case 1380) it is the documented table-10
  difference: COAL finds the gene behind its 5' UTR with TGA = Cys, PAL2NAL (TGA = stop) places it at offset 0 with
  many mismatches. Both engines agree byte for byte; untouched `main` before PR #43 fails identically.
- **Reproducer:** `verify/repro/V-1/c1380/` and `c1957/` (the fuzz's kept cases with `args.txt`).
- **Why it matters:** the nightly fuzz with a fresh seed will fail on roughly one seed in three for this reason, and
  the "0 unexplained differences" claim (README, SKILL.md, tests/README.md, site) holds for the documented seeds only.

### V-2 (MINOR; drop-in) — `pal2nal` without arguments prints its usage to stdout and exits 0

- **Area:** `src/coal_compat/cli.py`.
- **Steps:** `pal2nal > out.txt; echo $?` (and `perl pal2nal.pl > out.txt; echo $?`).
- **Expected:** the drop-in's own help says "1 for a usage error or bad input"; PAL2NAL v14 prints its usage to
  stderr (exit 0, as it always does), so a mis-called pipeline step gets an empty output file.
- **Actual:** the usage text (1,587 bytes) goes to stdout and the exit status is 0: a pipeline that calls `pal2nal`
  with no files gets the usage text as its "alignment" and no error status. `-h` behaves the same (fine for `-h`).
- **Affects:** the drop-in only (`coal` without arguments is an argparse usage error, exit 2, usage on stderr).

### V-3 (MINOR / DOCS; both engines) — a `kill -9` of a run writing to standard output leaves `.coal_<random>.tmp` in the temp folder

- **Area:** `docs/usage.md` ("an error never leaves a partial file"), the signal section.
- **Steps:** `coal big.pep big.nuc > out.fa &` then `kill -9 $!` (the test suite's `test_workers_die_with_the_main_process`
  does this); look in `$TMPDIR`.
- **Expected / actual:** SIGTERM, SIGINT, SIGHUP and a closed pipe leave nothing (checked, both engines); SIGKILL
  cannot run a handler, so the output's temporary file stays in the temp folder (two files of 13 and 15 MB after the
  engine suite here). Not a bug the program can fix; the docs should say that a killed (not terminated) run to stdout
  can leave `.coal_*.tmp` in the temp folder, as they do for the resumable files next to `-o`.

### V-4 (DOCS; app) — `demo/mcp_demo.md` shows `codon_align` stopping at the broken gene, but the tool now keeps going

- **Steps:** `python demo/make_mcp_demo.py` (in the `[mcp]` venv) and `git diff demo/mcp_demo.md`.
- **Expected:** the file says "Everything below the script is the real output of the tool functions, regenerated with
  `python demo/make_mcp_demo.py`".
- **Actual:** the committed output has `"completed": false`, `ok 13`, `not_checked 6`, `exit_code 1`, `output_path null`
  and the next step "COAL stopped before writing an alignment"; the tool now runs with `--keep-going` and returns
  `"completed": true`, `ok 19`, `failed 1`, `exit_code 2` and an `output_path` (27 lines differ). Regenerate it.

### V-5 (DOCS; Windows app) — `apps/coal-app/desktop/README.md` tells users to download `COAL-Setup-0.4.2.exe`

- **Steps:** `grep -n "COAL-Setup" apps/coal-app/desktop/README.md`; `gh release view` (latest v0.5.0, assets
  `COAL-Setup-0.5.0.exe` and `COAL-Setup-windows-x64.exe`).
- **Expected:** the version-free name the site uses (`COAL-Setup-windows-x64.exe`) or the current version. Every other
  version string in the repository is 0.5.0 (`verify/p9/versions.json`).

### V-6 (DOCS; agent skill) — `.claude/skills/coal/SKILL.md` carries sections written before the current engine

- "After the engine update (only if `python -m coal --help` lists these)": `--keep-going`, `--report` and normalised
  pairing are in the engine (`--help` lists them; phase 4 exercised them), so the "being built on the engine branch;
  do not use them" wording is stale and tells an agent to hesitate over released options.
- "On Windows (expected, not yet verified there): Ctrl-C gives 130 and a closed pipe 141, but a process killed from
  outside exits 1 without cleaning up": the Windows CI runs those tests (`test_decisions`, `test_robustness` subsets);
  the Windows pass of this verification should settle the wording.
- "`codon_align` returns a `report_path` even with the current engine (the app builds the report from COAL's
  messages)": the engine writes `coal.report/1` itself (`--report`), the app passes it through.
- The install block says `pip install -e . -e "apps/coal-app[gui]"` and, two lines later, "Install `gui` and `mcp` in
  separate environments": right, but the hand-off kit's setup scripts (not in the repository) still install both.

### V-7 (DOCS; tests) — "17 local codon alignments rebuilt exactly" is not reproducible from the repository

- **Area:** `tests/README.md`, "Real data".
- **Actual:** the repository holds 12 datasets (all rebuilt exactly in every format by `test_realdata.py`, both
  engines, phase 1); the 17-alignment figure covers PAML's GPL examples, PAL2NAL's own test set and GISAID-derived
  data that are "run only locally (not redistributable or not committed)", as the same paragraph says. The number is
  honest but cannot be checked by a reader; the prompt asked for it to be reproducible from the repo.

### V-8 (DOCS / MINOR; GUI) — the web page loads one script from a CDN, so "nothing leaves your computer" is not quite true

- **Area:** `apps/coal-app/desktop/README.md` ("Nothing leaves your computer. COAL works offline. It doesn't upload
  your files, check for updates or collect statistics"), `ui.py`'s own note "Offline: no CDN fonts (system fonts only)
  and analytics disabled".
- **Steps:** open the GUI with every request to a host other than 127.0.0.1 blocked and logged (`verify/phase7.py`).
- **Actual:** each page load requests `https://cdnjs.cloudflare.com/ajax/libs/iframe-resizer/4.3.1/iframeResizer.contentWindow.min.js`
  (Gradio's own `<head>`; three requests over a session). No file content leaves the computer and the page works with
  the request blocked, but a request to a third-party CDN is made on every start, which the privacy text denies. Either
  drop the script (Gradio's `head`/`analytics` settings, or serve it locally) or say "no data leaves your computer".

### V-9 (USABILITY / MINOR; GUI, upstream Gradio) — axe-core: hidden duplicate tab buttons are focusable

- **Steps:** `verify/phase7.py` injects axe-core 4.14 (WCAG 2 A/AA) after a run; `verify/p7/axe.json`.
- **Actual:** one serious violation, `aria-hidden-focus`: `<div class="tab-container visually-hidden" aria-hidden="true">`
  holds two focusable `<button>`s ("Align", "Settings"), Gradio's own second copy of the tab bar; a screen-reader user
  tabs onto buttons that are announced as hidden. No colour-contrast violation; the keyboard-only flow (Tab to an
  example, Enter, Tab to Start alignment, Enter) works and focus is visible. Upstream (Gradio); COAL could hide the
  duplicate with CSS `display:none` in its own stylesheet.

### V-10 (USABILITY / MINOR; GUI) — the genetic-code tip on a large input says "50 of 50 genes", not that 50 of 100,000 were sampled

- **Steps:** with the Settings default code set to 2, align 100,000 clean standard-code genes through the path boxes.
- **Actual:** the tip reads "With code 2 (Vertebrate mitochondrial) 50 of 50 genes have mismatched codons; with code 1
  (Standard) all 50 match", while the command line says "COAL tested 50 of the 100,000 sequences (a sample of the
  100,000 with problems): under genetic code 1 ... 50 of them match". A biologist reading the GUI cannot tell that the
  numbers are a sample of the file's 100,000 genes. Same data, phrasing only.

### V-11 (USABILITY / MINOR; GUI) — the "Methods paragraph" label is not linked to its text box

- **Steps:** `verify/phase7.py`: Playwright's `get_by_label("Methods paragraph")` finds nothing while
  `get_by_label("Raw COAL messages")` finds its textarea; the DOM shows the first textarea without an accessible name.
- **Actual:** the methods text box has no programmatic label (a Gradio `Textbox` rendered without the `for`/`aria`
  link, probably `show_label=False` with a Markdown heading above), so assistive technology announces an unnamed
  text field. axe reports it under "incomplete" rather than as a violation.

### V-12 (USABILITY / MINOR; GUI) — an unreadable input file gives "COAL stopped: no alignment was written" with the reason hidden

- **Steps:** `verify/gui_followup.py` and `verify/gui_probe2.py`: upload a binary file named `garbage.fasta` as the
  protein alignment (a tiny valid DNA file as the other), Start alignment.
- **Expected:** the banner (the one thing a user reads) says what was wrong and what to do, as the engine does:
  `coal: error: garbage.fasta:1: control character: sequence files must be plain text (is this a binary or
  word-processor file?)`.
- **Actual:** the banner reads "✖ COAL stopped: no alignment was written" with the "Files used" line; the health table
  is empty; the engine's sentence appears only in the "Raw COAL messages" box inside the collapsed "Methods, commands
  and glossary" accordion, and the alignment view says "Could not draw the alignment: garbage.fasta: not FASTA or
  CLUSTAL". A good run afterwards works (recoverable). The `.bin` extension case is fine: Gradio refuses the upload
  with a toast "Invalid file type only .fasta,.fa,.faa,.aln,.fas,.clustal,.txt allowed." (the banner then says "Please
  add both files", which is right, the file was never added).
- **Reproducer:** `verify/repro/V-12/` (the garbage file and the steps); screenshot `verify/screenshots/11_bad_upload.png`.

### V-13 (MINOR; tests) — `tests/test_docs_explainer.py` only ever runs on Linux with `PLAYWRIGHT_BROWSERS_PATH` or `/opt/pw-browsers`

- **Area:** `tests/test_docs_explainer.py`, `_chromium()`.
- **Actual:** the test looks for `chromium-*/chrome-linux/chrome` under `$PLAYWRIGHT_BROWSERS_PATH` or `/opt/pw-browsers`,
  so with Playwright's default browser folder (`~/.cache/ms-playwright`, `%LOCALAPPDATA%\ms-playwright`) and on
  Windows or macOS (`chrome-win`, `chrome-mac`) it is always skipped ("playwright with a Chromium build not
  available"): skipped in WSL in phase 1 and on Windows in the Windows pass, both with Playwright and Chromium
  installed. Run here by hand against the Windows Chromium (`verify/logs/p8-explainer.log`): see the Windows pass.
  `playwright.chromium.launch()` without `executable_path` would find the browser everywhere.

### V-14 (MINOR; tests, Windows) — `tests/test_cpp_vs_reference.py` errors on a Windows account without the symlink privilege

- **Steps:** `python -m unittest tests.test_cpp_vs_reference` on Windows outside Developer Mode (a normal user account).
- **Actual:** `build_shim` calls `os.symlink` and dies with `OSError: [WinError 1314] A required privilege is not held
  by the client` (both engines, `verify/logs/p8-tests-python.log`, `p8-tests-cpp.log`). GitHub's Windows runners
  have the privilege, so CI never sees it; a developer running the suite on a stock Windows account gets an ERROR
  that is not about COAL. A copy (`shutil.copytree`) or a `skipUnless` on the privilege would make the suite honest
  on Windows.

### V-15 (MINOR; tests, Windows) — `test_resume_with_workers` counts the venv launcher's child as a worker on Windows

- **Steps:** on Windows, from a venv (`.venv\Scripts\python.exe -m unittest tests.test_engine_extras.TestResume.test_resume_with_workers`).
- **Actual:** `AssertionError: Lists differ: [3] != [2]` every time (5 of 5 runs, Git Bash and PowerShell). The venv's
  `python.exe` on Windows is a launcher that runs the real interpreter as a child, so `worker_processes(pid)` (every
  descendant whose command line mentions python) sees the coal interpreter plus its two workers. CI runs the suite
  without a venv and passes. The run itself resumes byte-identically (the rest of `TestResume` passes on Windows,
  both engines). Counting only processes whose command line is a multiprocessing worker, or walking down through a
  same-executable child, would fix the count.

### V-16 (MINOR; tests, GUI) — no test checks what the "Download codon alignment" button hands out

- **Steps:** phase 13, `verify/mutations.py` recipe "app: GUI download hands out the report instead of the alignment"
  (`ui.py`: the download component given `report_path` instead of `output_path`), then `python -m unittest discover
  -s apps/coal-app/tests` from both app environments.
- **Actual:** 257 tests, OK. The GUI layout tests cover the banner, table and view texts but not the file behind the
  download button; a test that runs an example through `run_alignment`'s outputs and asserts the download value is
  the alignment path (and the report button's the report path) would catch it. The other app mutant (MCP `all_ok` true
  with a failed gene) is caught by four tests.

### V-17 (USABILITY / MINOR; GUI) — the residue-class colours of the alignment view are not distinguishable for colour-blind readers

- **Steps:** `verify/colour_vision.py` (added after the first report, for the human-checklist item "colour"): the
  one-broken-gene example rendered in Chromium with its vision-deficiency emulation (`verify/p7/colour/view_*.png`)
  and the seven class colours simulated numerically (Machado et al. 2009), smallest CIE76 distance between any two
  classes: normal vision 5.8 (Hydrophobic `rgb(212,235,249)` against Positive `rgb(201,220,239)`), protanopia 4,
  deuteranopia 4 (Aromatic against Polar), tritanopia 4, greyscale 0 (every class the same grey).
- **Actual:** the pastel palette tells the classes apart mainly by hue at nearly equal lightness, so two of them are
  already close for everyone and most pairs merge under red-green deficiency. Nothing is lost in meaning: every cell
  carries its residue letter and codon, mismatches and frameshifts are outlines and hatching, and the legend names
  each class, so colour is a secondary cue (docs/ui-principles.md: never colour alone). A palette that also varies
  lightness (or a "no residue colours" switch) would make the classes readable for colour-blind users.

## 7. What was run, phase by phase

Machine: Ubuntu 26.04 LTS under WSL 2 (kernel 6.6.114) on Windows 11 Pro 10.0.26200; Intel Core i5-13600KF, 16 threads
and 29 GB given to WSL (20 threads, 32 GB on the host); Python 3.13.13 (conda), g++ 15.2.0, CMake 4.2.3, Perl 5.40.1,
gh 2.46.0, Node 22.22.1 (axe-core 4.14.0 ran on the Windows side), Playwright 1.64 with Chromium 1248 (the Linux
Chromium needs system libraries that cannot be installed without sudo, so the headless browser ran on the Windows host
against the GUI served from WSL); `diamond` absent; Inno Setup absent; Claude Desktop installed on the host (this run
happened inside it); Windows available (its own pass of the prompt follows). PAL2NAL v14 fetched by
`scripts/fetch_pal2nal.py`, sha256 `2b419417e304e9a63129056d61fa513a4c5834257a7b162cc4389c3ac53e23f2`. The clone's
engine venv has `[fast]`, psutil, pyyaml, playwright, jsonschema and the app's `[gui]` extra; `.venv-mcp` has the engine
and the app's `[mcp]` extra (the two extras pin different `mcp` versions and cannot share a venv, as
`apps/coal-app/pyproject.toml` says). The generated inputs lived under `/mnt/d/coal-stress/verify/` (D:) because the
WSL disk image sits on the small C: drive; every log is under `/mnt/d/coal-stress/verify/logs/` (not committed).

**Phase 0 (matrix).** `verify/matrix.py` holds 174 rows (every option of `coal`, `coal-cpp` and the drop-in, every MCP
parameter, every GUI control, every installer option, the numeric and factual claims of README, SKILL.md, CHANGELOG,
the desktop README, the site and docs/, and the prompt's other phases); `verify/MATRIX.md` is rendered from it. The
Python and C++ help texts differ in one line only, `-j`: "worker processes (default 0 = automatic; 1 = one process)"
against "worker threads (default 0 = all hardware threads)", which `tests/test_cpp_cli.py` substitutes on purpose
(HELP_DRIVER). `pitch/QA.md` and `demo/race.py`, which the prompt names, were removed from the repository with the
rest of the hackathon material (commit `dfe0049`): their rows are SKIPPED.

**Phase 1 (existing suites, with timeouts).** `python -m unittest discover -s tests` with `COAL_CPP_BIN`: 502 tests,
OK, 8 skipped; the same with `COAL_IMPL=cpp`: 502, OK, 8 skipped; `tests/test_cpp_vs_python.py`: 9, OK;
`apps/coal-app/tests` in the gui venv: 257, OK, 1 skipped; in the mcp venv: 257, OK, 6 skipped;
`apps/coal-app/desktop` `test_wrapper test_smoke`: 77, OK, 2 skipped; `tests/fuzz_vs_pal2nal.py` with the default seed
(3,000 cases × 4 formats), Python engine and `COAL_FUZZ_ENGINE=cpp`: OK. The printed golden line is the one the docs
quote: `39 cases, 285 checks ... stdout byte-identical 238, warnings identical 234, with -nomismatch 48, both fail 24,
deliberate error 6, deliberate same_output 4, deliberate differs 3, deliberate nogap_code 14`. Every skipped test and
its reason: `test_robustness.TestDiskFull` (2: a tiny tmpfs needs root); `test_docs_explainer` (2: Chromium's shared
libraries are missing in WSL, see the machine note; the explainer page is covered by phase 10's render check on the
Windows side); `test_cpp_vs_reference.TestCppVsPal2nalFuzz.test_fuzz` (needs `PAL2NAL`: re-run with it, OK);
`test_stress_findings.QuadraticExactFitScan` and `UnrelatedTwentyThousand` (COAL_SLOW: re-run with it, OK);
`test_mcp_tools.TestFrozenAppMcp` (needs the frozen Windows exe: Windows pass); in the mcp venv the three `test_gui`
tests and `test_gradio_api_keeps_the_positional_order` (no gradio there, by design) and `test_packaging.TestSdists`
(needs `build`: re-run with it, OK); in `desktop/`, `test_exclusive_holder_is_detected` (Windows-only socket option)
and `WindowApi.test_documented_signatures` (pywebview is Windows-only here). `COAL_REF_SRC`, `COAL_CPP_BIN` and
`COAL_CPP_DUMP` default to this checkout's `src/` and `build/` (read from the test files). CI of `main` at `4159e10`:
green (one earlier run of the PR had the known intermittent worker-kill test on an Ubuntu 3.10 runner; it passed 16 of
16 local repetitions on 3.10 and on the rerun).

**Phase 2 (correctness against the real PAL2NAL).** Fuzz with three new seeds, 2,000 cases × 4 formats each, Python
engine and `coal-cpp` separately: seed 31337 0 unexplained; seeds 20261024 and 20261025 one unexplained difference
each, both `-codontable 10 -nogap`, which is V-1 (a neutralisation gap of the checker; COAL's output is right, both
engines agree). Per-format counts for the Python engine: 20261024 fasta/clustal/paml 1,590 identical, 223 both fail,
186 deliberate, 1 unexplained, codon 819 / 1,138 / 43 / 0; 20261025 1,580 / 220 / 199 / 1 and 836 / 1,109 / 55 / 0;
31337 1,558 / 246 / 196 / 0 and 817 / 1,139 / 44 / 0. The adversarial generator (`verify/adversarial.py`, 50 cases of
the prompt's list, every case through PAL2NAL, the Python engine and `coal-cpp` with `-output fasta`): the two engines
were identical (stdout, stderr, exit status) on all 50; against PAL2NAL 28 identical, 11 both fail, 3 "PAL2NAL output,
COAL error" (duplicate IDs, non-ASCII headers, a CLUSTAL file without its header line: COAL's errors name the file and
line and say what to do), 1 "COAL output, PAL2NAL failed" (a UTF-8 BOM), 2 different outputs that are the documented
deliberate differences (versioned IDs in another order: COAL pairs each protein with its own DNA, PAL2NAL pairs by
order and gets every gene wrong; `-nogap` with code 2: COAL keeps TGA = Trp columns), 3 not run through PAL2NAL
(standard input, which PAL2NAL has no mode for; the 1,000,000-codon gene, on which PAL2NAL ran past its 600 s
timeout while both engines took under 10 s). The reverse-complement gene is written with per-codon warnings and the
"probably unrelated" notice, as the docs say (the report marks it failed; `--strict` refuses it). `-o` naming an
input file is refused before anything runs and the input is intact; a read-only output folder (checked on the Linux
file system, where permissions apply) gives `cannot write ro/out.fa: permission denied (choose a folder you can write
to)`, exit 1, no leftovers, in both engines. Real data: `test_realdata.py` rebuilds the 12 committed datasets in every
format with both engines (phase 1); the "17 rebuilt exactly" figure also counts local-only sets (V-7).

**Phase 3 (engine parity and the automatic switch).** Everything in phase 2 was run through both engines and compared
(above). `-j 1`, `-j 2` and `-j 0` (all) on 20,000 genes of 120 codons with late mismatches gave one distinct output
and one distinct stderr per format in both engines (`verify/p3`). The app's hand-over: a 1.5 MB input through the MCP
tool ran in `coal-cpp` (the report's `engine` field says so, phase 6); a 12,000-gene input just under 1 MB ran in
Python; with `COAL_CPP` unset the app uses the Python engine (phase 1's app tests cover the missing-binary message).

**Phase 4 (every flag and contract).** `verify/phase4.py`, 199 checks, both engines: `--report` and `--json` on a
success, a warnings run, a failure and two `--keep-going` runs (fasta and paml, `-j 2`) validate against
`docs/report.schema.json`, the sha256 of every input and of the output match the files, `output.sequences` equals the
records written, the counts add up, the engines' reports are identical apart from `engine` and `runtime_seconds`;
`--strict`'s three rules (unrelated DNA, guessed position via the short-gene case in the suite, partial IDs) exit 1
with nothing written; `-q` leaves stdout and the exit status unchanged with an empty stderr; `--max-warnings`:
default 20 lines then "... and N more warnings", 0 = all, 5 = 5, -1 a usage error (exit 2), and the drop-in prints
them all (PAL2NAL compatibility; usage.md and CHANGELOG 0.4.2 agree: the "unlimited by default" of earlier specs is
gone from the docs); `--progress` prints nothing for 50 sequences and progress lines on stderr only (never stdout)
for 12,000, also with `-q`; `--resume`: a 20,000-sequence run with `--report` killed after its fourth checkpoint and
resumed gives a byte-identical output and an identical report in all four formats, both engines, with no checkpoint
or temp file left, and a changed input is refused with "changed since the interrupted run"; SIGINT 130, SIGTERM 143
(only the resumable files stay, the lock goes), a closed pipe 141 and silent, no `.coal_` files next to the output
(V-3 for `kill -9`); pairing: normalised IDs (`AT1G0k.1` vs `AT1G0k`, DNA reversed) give `pairing: normalised` with the
notice and 4 of 4 ok, two DNA IDs that tidy to one protein ID are not paired by tidied IDs, no matching ID gives
`pairing: order` with the `id_mismatch` hint; the mito example suggests code 2 (hint fit 24 > fit_current 0; the early
notice says `Re-run with -codontable 2`) and is clean under code 2, one wrong codon under the right code gives no
hint; methods_text names COAL v0.5.0, the plumat-coal URL, the code's name and the counts, and `chloe-newt` appears
nowhere in the repository.

**Phase 5 (the drop-in).** `verify/phase5.py` against PAL2NAL v14: 18 old-pipeline command lines (default CLUSTAL,
each `-output`, `-codontable`, `-nogap`, `-nomismatch`, `-blockonly`, `-nostderr`, options before files, dangling
`-codontable` / `-output`, `-codontable 2x`, two DNA files) all byte-identical or both failing, except
`-codontable 2 -output paml -nogap`, the documented deliberate `-nogap` difference; `-nostderr` is silent;
`-nomismatch` prints no mismatch warnings; `-html`, `--strict`, `-j`, `-o`, `--report` and `-foo` are refused with a
clear error and exit 1; `-` reads standard input; a closed pipe exits 141 silently on a 20,000-gene output; a count
mismatch is exit 1 with a message; the `pal2nal` entry point works from a fresh login shell. V-2 is the only note.

**Phase 6 (MCP server).** `verify/phase6.py`, a stdio JSON-RPC client with timeouts: initialize (serverInfo.version
0.5.0), tools/list (exactly `codon_align` and `explain_report`, nine described parameters, the 17 codes and the 4
formats as enums, a docstring that names all_ok / needs_attention), 41 checks, 0 findings: clean, warning, broken and
mito inputs give the right completed / all_ok / needs_attention / counts / suggested_codontable / next_step /
problems; every `output_format` and every boolean option gives the CLI's bytes; a missing file, a file of the wrong
kind, code 7 and format `nexus` are tool or validation errors; `save_as` writes a new file, never overwrites an
existing one, saves into an existing folder, and a missing folder is a clear error; a 1.5 MB input runs in
`coal-cpp`; `explain_report` is plain English and right for a clean, a broken and a mito report (texts in
`verify/p6/explain_*.txt`), and a non-report or missing file is a tool error; nothing but JSON-RPC appeared on stdout
and no traceback on stderr. A relative path is accepted (resolved against the server's folder) although the docs ask
for absolute paths: harmless, noted.

**Phase 7 (web GUI, Playwright).** The GUI ran in WSL (`python -m coal_app --port 7860` in the `[gui]` venv, 90-minute
timeout); Playwright's Chromium ran on the Windows host against it (WSL's Chromium lacks `libnspr4`, which needs root):
`verify/phase7.py` (39 checks), `verify/gui_followup.py` (18 checks on fresh pages for the steps the first driver got
wrong: it clicked the already-open "View alignment" accordion shut and re-used file inputs that disappear after an
upload) and two probes (`gui_probe2.py`, `gui_probe3.py`). Every control: two uploads align (2 genes, the download has
2 records); the path boxes take the 100,000-sequence input (result in 10.7 s, the page answered within 2 s throughout);
each example (clean: "All 8 genes match codon for codon"; mitochondrial: "Wrong genetic code?" with the counts and a
"Re-run with genetic code 2" button that gives a clean result; messy IDs; one broken gene: "19 of 20 genes aligned, 1
failed", the failed row first with a plain sentence and a "Column 102" link); the Sample report button; the 17 codes
and 4 formats (every download byte-identical to the CLI's output for the same example and options); the four
checkboxes; the report download (valid `coal.report/1`); the health table (problem rows first, symbol plus the word;
the column link scrolls the page to the cell `#coalv-OG0001_sp14-c102`, which is drawn as a frameshift cell); the
alignment view (the preview line says what it shows: "100 of 100,000 sequences, spread evenly through the file; all
100 columns", and on a 10,000-codon input "columns 1–600 of 10,000" with the centre-on-column box; centring on 15,000
or 99,999,999 clamps to the last window, no error; "Download alignment view (HTML)" gives a 168 KB page); the methods
paragraph and its copy button; the raw messages box; the Settings tab (code, format and the three checkboxes persist
across page loads, "Reset to defaults" restores code 1, the `settings.json` path is shown; a server restart was not
exercised). Bad uploads: a `.bin` file is refused by the file box itself ("Invalid file type ..."); a binary file named
`.fasta` gives "COAL stopped: no alignment was written" with the engine's reason only in the raw messages box (V-12);
a good run afterwards works. Light and dark themes render (dark body background checked), 200 % zoom (683 CSS px, scale
2) has no horizontal overflow, screenshots at 1366×768 and 1920×1080 are in `verify/screenshots/`. Offline: every
request went to 127.0.0.1 except Gradio's iframeResizer script from cdnjs (V-8). Accessibility: the keyboard-only
flow (Tab to an example, Enter, Tab to Start alignment, Enter) works with visible focus; axe-core 4.14 (WCAG 2 A/AA)
reports one serious violation, Gradio's hidden duplicate tab buttons (V-9), no contrast violation; the methods box has
no programmatic label (V-11); the genetic-code tip on the 100,000-sequence run says "50 of 50 genes" (V-10). The GUI
also serves the MCP tools (`/gradio_api/mcp/sse` answers, the schema lists both tools).

*Usability review (G26, against `docs/ui-principles.md`, from the screenshots and page texts).* The three flows read
as the principles ask: one sentence of verdict first ("All 8 genes match codon for codon"; "Wrong genetic code?";
"19 of 20 genes aligned, 1 failed"), then one next step (the re-run button; the failed gene's row with what to check),
never a code or a raw message above the fold; status is a symbol plus a word, never colour alone; the per-gene strip
under the banner shows the failure's position; the methods paragraph starts with "Fix the warnings first" when the run
is not ready to use; the glossary explains every term the page uses. Rough edges seen: the hint banner uses the blue ℹ
while the health row below calls the same thing a ⚠ Warning (two severities for one fact); the strip of squares has
no legend or tooltip; the 100,000-sequence health table has an empty "Sequence" cell on its one row; and the two
usability findings above (V-12: a file-level error with the reason hidden; V-10: the sampled-genes tip). None of these
would stop a biologist getting a right alignment or understanding why they did not.

**Phase 8 (Windows desktop and installer; the native Windows pass).** Run from the Windows host (Windows 11 Pro
10.0.26200, Python 3.13.7, the Windows clone `D:\coal-handoff-kit\COAL` at the same tree as `main` 4159e10, MSVC-built
`coal-cpp.exe`), with `verify/phase8.py` (45 checks) and `verify/phase8_install.ps1` (25 checks, logs in
`verify/logs/p8-*.log`, states in `verify/p8/state_*.json`). Release: `COAL-Setup-0.5.0.exe` and
`COAL-0.5.0-windows-x64.zip` downloaded with `gh`, both match `SHA256SUMS.txt`. SmartScreen: the owner launched the
installer with the mark of the web set before the run and recorded "Windows protected your PC", App
`COAL-Setup-0.5.0.exe`, Publisher "Unknown publisher", "More info" then "Run anyway", which is what the desktop README
and site/install.html tell users to click; Smart App Control is off on this PC, so that text could not be exercised.
The wizard (run by the owner) installed per-user without a UAC prompt with "Connect COAL to Claude Desktop" ticked
(the installer's default) and created `%APPDATA%\Claude\claude_desktop_config.json` with the `coal` entry; no
`.bak-coal` because there was no file before. Installed app: `COAL.exe --version` (COAL desktop 0.5.0), `-m coal
--version` (coal 0.5.0), `-m coal --cpp-path` names `_internal\coal-cpp.exe` (coal 0.5.0 (C++)); a 13.7 MB input
through `-m coal`, through the MCP server and through the window's path boxes reports engine `coal-cpp`; `-m other` is
refused with a plain message. `COAL.exe -m coal_app --mcp` (what Claude Desktop runs): initialize, both tools,
`codon_align` on the clean example and the 13.7 MB input, `explain_report`, nothing but JSON-RPC on stdout, no
traceback. `--connect-claude` against a config holding another server adds the `coal` entry, keeps the other server
and unrelated keys, writes `claude_desktop_config.json.bak-coal` byte for byte, is idempotent; `--disconnect-claude`
removes only `coal`. The local app (`COAL.exe --serve --port 7870`, the window's child) is ready in 5.4 s, listens on
127.0.0.1 only, runs the four examples, a path with spaces inside a `OneDrive - Plumat Lab\my sequences (2026)` folder,
the 13.7 MB input and the dark theme; the real window (`COAL.exe`) logs "local app ready" 3.7 s after start (two
`COAL.exe` processes), and killing it leaves none. Long paths could not be tested: `LongPathsEnabled` is off on this
PC and the harness itself could not create a 270-character path (human checklist, with drag-and-drop and high-DPI).
Uninstall (silent) removes the app folder, the Start Menu shortcut, the Add/Remove entry, `%LOCALAPPDATA%\COAL` and
`%TEMP%\coal_app_runs`, and the `coal` entry from Claude's config while keeping the other server and keys; it leaves
`%TEMP%\gradio` as the README says and does not put `.bak-coal` back over later edits. Reinstalled with
`/TASKS=desktopicon,openwith`: a desktop shortcut and `COAL.SequenceFile` under `OpenWithProgids` for `.aln .fasta
.fa .fas .faa .fna` with the command `"COAL.exe" "%1"`, the `.fasta` default program and UserChoice untouched; the
next uninstall removes all of it. The machine was left as it started: a default install with the `coal` entry in
Claude's config (the owner's choice). Portable zip: `COAL.exe --version`, `--cpp-path` names its own `coal-cpp.exe`,
the 13.7 MB run in `coal-cpp`, `--serve` ready in 6.0 s with the four examples. The engine suite on Windows
(`verify/logs/p8-tests-python.log`, `p8-tests-cpp.log`): 502 tests, 88 skipped (Linux-only signals, /proc, tmpfs,
PAL2NAL and playwright gates), 2 failures on both engines, neither about COAL's output: `test_cpp_vs_reference`
needs the symlink privilege (V-14) and `test_resume_with_workers` counts the venv launcher (V-15); the exit-code
tests (`test_decisions`, `test_robustness`) pass, so SKILL.md's "expected, not yet verified there" is settled (V-6).
`tests/test_docs_explainer.py` is skipped on Windows by its Linux-only browser lookup (V-13); run by hand against the
Windows Chromium its 2 tests pass (`verify/logs/p8-explainer.log`).

**Phase 9 (distribution and versions).** `verify/phase9.py`: the v0.5.0 assets download and match `SHA256SUMS.txt`; in
clean venvs on Python 3.10.22 and 3.13.13 the pure wheel plus the app wheel install, `coal --version` is `coal 0.5.0`,
the clean example runs, the `pal2nal` entry point and `python -m coal_app --check` work, and `coal --cpp-path` exits 1
with a reason (no bundled binary); the manylinux platform wheel installs, `--cpp-path` finds `coal/bin/coal-cpp`, whose
`--version` is `coal 0.5.0 (C++)` and whose output equals the Python engine's; the release sdists install;
`python -m build --sdist` from the checkout builds and installs; `pip install git+https://github.com/Plumat-COAL/COAL`
works. Version strings (`verify/p9/versions.json`): `src/coal/__init__.py`, both `pyproject.toml`, the app's
`coal-align==` pin, `cpp/src/main.cpp`, `cpp/src/report.cpp`, `CMakeLists.txt`, `CITATION.cff`, `site/index.html`,
`site/cite.html`, `tests/report_example.json`, the CHANGELOG's latest release heading, `coal --version`,
`coal-cpp --version` and the release tag all say 0.5.0; `apps/coal-app/desktop/README.md` says 0.4.2 (V-5).

**Phase 10 (docs and website).** `verify/phase10.py`: every `coal` / `pal2nal` command in README, SKILL.md, the
desktop README, docs/*.md and site/*.html that names a file in `examples/` runs with the exit code shown, and the
outputs the site quotes match (the mito example: 20 warnings then "... and 159 more warnings", the code-2 summary;
one_broken_gene: the OG0001_sp14 error, "Frameshift at column 102", "6 sequences were not checked", exit 1, and 19
genes with `--keep-going`, exit 2; the clean example's 8 genes). The live site (https://plumat-coal.github.io/COAL/)
serves exactly the repository's `site/` files, shows 0.5.0, every link answers (the download links are release
assets, the sdist checksums match `SHA256SUMS.txt`), and the install page has the Windows / macOS / Linux steps and
checksum commands. Numbers: the test counts (285 / 238 / 24 / 23 / 39 / 4) are the printed line of phase 1; the
compatibility paragraph is identical in README, SKILL.md, tests/README.md and the site (`test_docs_numbers`); the
speed-ups (74× / 223× / 297× / 1,115×, 12.5 s, 172 MB, 62 min, 18.6 GB) trace to `docs/benchmarks.md` and
`benchmarks/results/windows_47f57f1/SUMMARY.md` (phase 11 compares this machine's numbers); the "0 unexplained
differences" claim holds for the documented seeds (V-1 for new ones). None of the claims the team decided not to make
appears (NVIDIA Brev, a hosted web version, PyPI availability, code signing, DIAMOND on real genomes); "C++ via pip" is
now true in the form the docs state it (platform wheels on the release page; phase 9 installed one). Stale text: V-4,
V-5, V-6.

**Phase 11 (performance claims).** `benchmarks/run_bench.py` at 10,000 × 500, 100 × 10,000 and 100,000 × 500, both
engines (`--impl both`, 1 process/thread and all cores, 3 repeats after a warm-up), then the same script against
PAL2NAL v14 at the two smaller sizes, first with the data on the Windows drive (DrvFs, as the rest of the run) and
again on the Linux file system, where file I/O no longer dominates (the 100,000-sequence Python run drops from 23.8 s
to 4.1 s). Every output pair was identical. The numbers and the comparison with the committed tables are in section 3:
the shapes hold (the Python engine ~20× faster than PAL2NAL on many sequences and ~60× on long genes here; coal-cpp
5–7× faster than the Python engine; all-core Python helps on 100,000 sequences, 4.1 s → 2.2 s, and not on long genes).
The 1,000,000-sequence row (62 minutes of PAL2NAL) was not re-run. Peak memory at 500,000 × 150 codons in every output
format with `tests/tools/stress_run.py` (summed PSS, `-j 1`, output to a file): Python 190 / 191 / 190 / 197 MB and
coal-cpp 112 / 137 / 113 / 149 MB for fasta / clustal / paml / codon, identical outputs, 49–62 s and 38–45 s; the
figures `docs/usage.md` gives for this shape (203 / 203 / 209 and 112 / 136 / 148) are met or bettered, so the
streaming claim and the per-sequence figures (about 350 and 190 bytes) stand (B02, K02 PASS). The usage.md
"100,000 unrelated genes: 118 / 69 / 885 MB" figures were measured by the Task A work on this machine two hours
earlier (`docs/stress-test-2026-10-23.md`) and were not repeated.

**Phase 12 (demos).** `demo/race.py` no longer exists (removed with the hackathon material): SKIPPED.
`demo/diamond_pipeline.py`: `diamond` is absent, so only its stand-in test ran (`test_diamond_pipeline`, 9 tests, OK).
`demo/mcp_demo.md`: regenerating it with `make_mcp_demo.py` changes 27 lines of the recorded tool output (V-4).

**Phase 13 (can the tests catch real bugs?).** `verify/mutations.py`: 15 single-line bugs, each in a fresh rsync copy
of the checkout outside the repository (`~/coal-kit/coal-mut`, deleted afterwards), then the relevant suites with
`PYTHONPATH` pointing at the copy (11 Python-engine bugs against the whole engine suite; 2 C++ bugs rebuilt and run
through the differential, parity, pool, CLI, decisions and reference suites with `COAL_IMPL=cpp`; 2 app bugs against
the app suite from both app environments, `verify/mutations_app.py`, after the first pass had left the app package
resolving to the unmutated checkout). Results (`verify/p13/results.json`, `results_app.json`): **12 of 15 caught, 3 not
caught, of which 2 are equivalent mutants and 1 is a test gap.** Caught: mismatch column off by one (171 failures),
`-nogap` keeping TAG stop columns (32), table 2 TGA as a stop (42), exit 0 after an error (140), report `ok` count off
by one (29), `-nomismatch` dropping the wrong column (24), the genetic-code hint picking the worst tied code (17), PAML
header off by one (49), `--keep-going` exiting 0 (25), C++ rows emitted in reverse batch order with `-j 4` (the
differential suite hangs rather than fails: the pool's ordered writer waits for a batch that never comes, so the run
is caught by the 25-minute timeout, not by an assertion), C++ accepting only ATG as the first Met (9 failures in the
reference suite), MCP `all_ok` true with a failed gene (4 failures). Not caught: (a) removing CTG from the standard
code's start codons in the Python engine: equivalent, the start set is read only by `_quick_offset0`, a fast path
whose miss falls back to the general search (checked by hand: identical output, report and exit with the mutant), so
no test can see it, while the C++ engine's start set does decide the alignment and its mutant was caught; (b)
`-blockonly`'s `keep_col` lambda shifted by one column: equivalent as far as probes reached (identical output, stderr
and report with a mismatch in a dropped column next to kept ones; the lambda filters warnings only in a path these
inputs did not take); (c) **the GUI's "Download codon alignment" button handing out the report file instead of the
alignment: not caught by any test (V-16)**, the one real gap. Every mutation applied cleanly; the suites took 40–56 s
per Python mutant, 1,510 s for the hanging C++ one.

## 8. Coverage matrix

`verify/MATRIX.md` (174 rows; the result column names the finding or the reason a row was skipped; the rows marked "Windows pass" are filled by the second, native-Windows run of phase 8).

## 9. Files

- `verify/LOG.md`: every decision, timestamped.
- `verify/matrix.py`, `MATRIX.md`: the coverage matrix.
- `verify/adversarial.py`, `phase4.py` … `phase10.py`, `gui_followup.py`, `gui_probe*.py`, `mutations.py`: the drivers (each prints its checks; results in `verify/p4` … `p13`, `adv`).
- `verify/logs/`: the drivers' and jobs' logs, the memory measurements (`p11-memory.jsonl`).
- `verify/p11/`: the benchmark CSVs of this machine.
- `verify/screenshots/`: the GUI at 1366×768, 1920×1080, 200 % zoom, dark mode, the 100,000-sequence run, the bad upload; `verify/p8/`: the Windows pass (installed and portable app screenshots, install logs, the machine state before, between and after the install cycles).
- `verify/repro/V-n/`: one folder per finding with the inputs and steps.
