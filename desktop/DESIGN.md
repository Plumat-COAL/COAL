# COAL for Windows: design

Status: design for the hackathon prototype and what comes after it. Written before any desktop
code. Everything in this folder (`desktop/`) is a thin wrapper: the science lives in the engine
(`src/coal`, run as `python -m coal`), the screens live in the existing local app
(`src/coal_app`, Gradio). If the desktop wrapper broke, results would be byte-for-byte the same,
because the wrapper never touches sequences.

Contents

1. What the desktop app is (one paragraph)
2. Three people and the moment each feature removes
3. Main flows, step by step
4. Screens: layout in words and wireframes
5. Windows decisions, with reasons
6. Large data
7. In scope for the hackathon demo, and later
8. Changes this design needs outside `desktop/` (for a human to decide)

---

## 1. What the desktop app is

A Windows program called **COAL** that you install without admin rights, start from the Start
menu, and use in its own window. Under the hood it starts the same local COAL app the
`python -m coal_app` command starts, on `127.0.0.1` (this computer only) on a random free port,
and shows it in a native window using the Edge WebView2 component that ships with Windows. No
browser tab, no terminal, no internet connection, no account, no AI calls at run time. Files
are read from your disk and results are written to your disk; nothing is uploaded anywhere.

Three promises shape every decision below:

- **No Linux, no terminal, no installer drama.** Double-click, click through SmartScreen once,
  done.
- **Your data stays on your machine.** The server only listens on 127.0.0.1, nothing is
  fetched from the internet, there is no update check, no telemetry.
- **Errors explained in plain English, and every run is reproducible.** Each run shows the
  equivalent command, a run receipt (report JSON with input checksums) and a ready-to-paste
  methods sentence.

---

## 2. Three people and the moment each feature removes

### Priya, bench scientist (molecular ecologist, Windows laptop from the university)

Has 20,000 orthogroups from OrthoFinder and wants dN/dS with HyPhy. A collaborator ran PAL2NAL
on a Linux server last year; that person has left. She pasted alignments into the PAL2NAL web
server one at a time for an afternoon, then 3 genes errored with "ERROR: inconsistency between
the following pep and nuc seqs" and she could not tell why.

Moments of pain and the feature that removes each:

| Moment | Feature |
|---|---|
| "I don't have Linux and IT won't install anything." | Per-user installer, no admin rights, no Python, no WSL. |
| "The download is blocked / Windows says it's dangerous." | Exact SmartScreen wording and the two clicks, in the README and the download page (section 5.7). |
| "3 genes out of 20,000 failed and the error is gibberish." | Health table, problems first, one plain sentence per gene ("one DNA letter is missing near codon 102, so every later codon is read out of step"); `--keep-going` so the other 19,997 are still written. |
| "Did I pick the right genetic code?" | Hint banner: "Most mismatches are TGA read as W. These look like vertebrate mitochondrial genes: re-run with code 2", with a one-click re-run button. |
| "Am I allowed to put this on OneDrive / is it sent somewhere?" | Statement on the home screen and About: "Nothing leaves this computer"; results written next to the input. |

### Sam, bioinformatician (core facility, builds pipelines for others)

Lives in a terminal on the cluster, but the people who bring data use Windows. Sam is asked
"which button did you click?" weekly and needs runs that can be reproduced on the cluster.

| Moment | Feature |
|---|---|
| "I need to rerun exactly what the student did, on the cluster." | Equivalent command for every run, for both `coal` and `pal2nal.pl`, with a copy button; the report JSON records version, options and SHA-256 of both inputs. |
| "Is the GUI doing something different from the CLI?" | The desktop app calls the same engine via the same command; the About screen says which engine (Python or C++) and version ran. |
| "Our pipeline breaks when the IDs have `.1` suffixes." | ID-pairing hint and "paired with DNA AT1G01010" shown per row, so mispairing is visible, not silent. |
| "100,000 sequences in the browser version just hangs." | Files read by path from disk (no upload), progress with time remaining, a Cancel button, health table that pages and leads with problems (section 6). |

### Lin, MSc student (first evolution project, never used a terminal)

Was told "make codon alignments for PAML". Doesn't know what a codon alignment is beyond the
words, and is anxious about doing something wrong in front of the supervisor.

| Moment | Feature |
|---|---|
| "What files do I even need?" | Home screen names the two inputs in words, with a "Try an example" button for each of the four bundled datasets. |
| "What does 'nogap' mean? Will ticking it break my thesis?" | Every option has a one-line explanation; a glossary ("What do these words mean?"); safe defaults. |
| "My supervisor asks how I did it." | Methods paragraph with a copy button, citing the version and options. |
| "I use a screen reader / I need big text." | Keyboard-only use, Narrator labels, follows Windows text size, dark mode and high-contrast mode (section 5.13 to 5.16). |

The PI writing the methods section is served by the same methods paragraph and the run receipt.

---

## 3. Main flows

Each flow lists what the person does (**bold**) and what the app does. "Today" marks what the
current local app (`coal_app/ui.py`) already does; the rest is target behaviour.

### 3.1 First launch

1. **Downloads `COAL-Setup-0.4.0.exe` and double-clicks it.** Edge may first say the file "isn't
   commonly downloaded" (section 5.7 has the exact clicks).
2. Windows shows "Windows protected your PC". **Clicks "More info", then "Run anyway".**
3. Installer: one page "Install COAL for you only (no administrator needed)", with two
   unticked boxes: "Create a desktop shortcut" and "Open .aln and .fasta files with COAL
   (adds COAL to 'Open with'; doesn't change your defaults)". **Clicks Install**, then Finish
   with "Launch COAL" ticked.
4. The COAL window opens at a sensible size (1200 x 800, or smaller on small screens) with a
   neutral splash: "Starting COAL… (this computer only, no internet needed)". The wrapper
   picks a free port, starts the local app, and swaps the splash for the home screen when the
   local app says it is ready (usually 3 to 8 seconds; Gradio import is the slow part).
5. Home screen: the two drop zones, the "Try an example" row, and one sentence: "Nothing
   leaves this computer." No sign-up, no tour, no "check for updates" prompt.
6. If WebView2 is missing (rare: very old Windows 10), a plain dialog says so and offers to
   open COAL in the default browser instead, still on 127.0.0.1 (section 5.1).
7. If COAL is already open, starting it again brings the existing window to the front instead
   of starting a second server (single instance).

### 3.2 Drag and drop a protein alignment and a DNA file, and run

1. **Drags `OG0001.aln.fasta` from Explorer onto drop zone 1** ("Protein alignment") and
   **`OG0001.cds.fasta` onto drop zone 2** ("Coding DNA"). Today: the Gradio file boxes accept
   drops inside the window. Target: dropping both files anywhere on the window works; COAL
   looks inside each file (gaps and amino-acid letters vs only ACGTN) and puts each in the
   right box, saying which it chose: "Protein alignment: OG0001.aln.fasta (20 sequences,
   aligned). DNA: OG0001.cds.fasta (20 sequences)." A wrong guess is fixed with a "Swap"
   button.
2. Options keep their defaults: output FASTA, genetic code 1, no column removal.
3. **Presses "Start alignment"** (or Enter when the button has focus, or Ctrl+Enter
   anywhere).
4. While running: the button becomes "Cancel", a progress line shows "Checking 12,400 of
   20,000 sequences, about 40 s left" (target; needs engine progress, section 8).
5. Done: banner "✔ Done: every sequence matched. 20 sequences, genetic code 1, 0.3 s."; health
   table; alignment viewer; downloads; methods paragraph; equivalent command. Today all of
   these exist except progress and Cancel.
6. Target: results are also written next to the input, in a folder
   `OG0001_coal\` (section 6.4), and the banner has "Open folder" which opens Explorer there.

### 3.3 A run with problems, and how the app guides the fix

**Wrong genetic code (example `mito`).**

1. **Runs with the default code 1.** Banner: "⚠ Done, with some problems (see the table)."
2. Hint banner (orange border, icon and words, not colour alone): "Tip: 24 of 24 sequences
   have codons that don't match, almost all TGA where the protein has W. That fits the
   vertebrate mitochondrial code (2)." Button: **"Re-run with genetic code 2"**.
3. **Clicks it.** The dropdown changes to 2, the run repeats, banner turns to "✔ Done". The
   methods paragraph and the equivalent command now say `-codontable 2`. Today: works.

**Mismatched IDs (example `messy_ids`).**

1. **Runs.** Today, with the current engine: "✖ COAL stopped: no alignment was written",
   because protein `AT1G01010.1` has no DNA called `AT1G01010.1`.
2. Hint: "The protein IDs end in .1, .2 … but the DNA IDs don't. They match once the version
   suffix is removed." Once the engine has ID normalisation, the run succeeds and each row
   says "paired with DNA AT1G01010", so the pairing is visible and checkable.
3. Target: a "Re-run matching IDs without version suffixes" button, mirroring the
   genetic-code one (needs the engine option, section 8).

**One broken gene (example `one_broken_gene`).**

1. **Runs.** With `--keep-going` (engine branch): "⚠ Done, with some problems. 19 of 20
   sequences OK, 1 failed." The codon alignment of the 19 good genes is written.
2. The failed row is at the top: "✖ Failed. OG0001_sp14: from codon 102 on, the DNA no longer
   spells the protein. One DNA letter seems to be missing there (a frameshift), often a
   sequencing or assembly error." A link "Show in alignment" jumps to column 102 in the
   viewer, where the cells are double-outlined and hatched.
3. What the person can do, listed under the row: fix or drop that gene and re-run; or keep
   the 19 and note it in the methods ("1 sequence was excluded because of a frameshift").
   The methods paragraph already counts excluded sequences.

### 3.4 Opening a large dataset (100,000 sequences)

1. **Clicks "Choose files…"** (native Windows file dialog) or drops two files of 120 MB and
   400 MB. Target: the app receives the *paths*, not the contents; nothing is copied
   (section 6.1). Today the Gradio file box copies the file into its temp folder over
   127.0.0.1: correct, but slow and doubles disk use for big files.
2. The home screen shows file facts before running: "100,000 sequences, 1.2 billion DNA
   letters, alignment 2,400 columns." It warns once if the destination disk has less free
   space than about 3 times the DNA file.
3. **Runs.** Progress line with time remaining, updated about twice a second; **Cancel** stops
   the engine within a second and deletes the partial output (a cancelled run leaves no
   half-written file that looks finished).
4. Results: the summary counts first, then only the problem rows (paged, 100 per page), then
   "Download the full table (TSV)". The alignment viewer opens on the first problem and
   shows at most 200 sequences x 2,000 columns around it, saying so ("showing sequences 1 to
   200 of 100,000"). Today the viewer already caps at these sizes.
5. The window stays responsive throughout: the engine runs in a separate process.

### 3.5 Exporting results and copying the methods sentence

1. After a run, **clicks "Download codon alignment"**. Desktop target: a native "Save as"
   dialog, starting in the input's folder, with the suggested name `OG0001.codon.fasta`.
   (Today: Gradio download buttons; inside WebView2 these show WebView2's own download
   behaviour. pywebview has a documented `create_file_dialog`, which the wrapper can use
   once the GUI hands it a path; section 8.)
2. **Clicks "Download report (JSON)"**: the run receipt (inputs with SHA-256, options, engine
   version, per-sequence status).
3. **Clicks the copy icon on "Methods paragraph"**, pastes into Word. The paragraph names
   COAL's version, the genetic code by number and name, counts, removed columns, and the
   command with file names only (no `C:\Users\priya\…` paths leak into a paper).
4. **Clicks the copy icon on "Equivalent commands"** and sends it to Sam, who reruns it on the
   cluster and gets the same file.

### 3.6 Re-running after changing an option

1. **Ticks "Remove gap columns (-nogap)"** or picks another output format or genetic code.
2. Target: the results area dims and says "Options changed since this run. Press Make codon
   alignment to update." Results never silently mix old output with new options.
3. **Presses "Start alignment".** The new run replaces the results on screen. Desktop
   target: the previous run's files are kept (each run gets its own folder,
   `OG0001_coal\run-2026-10-17_0912\`), so nothing is overwritten without asking.

---

## 4. Screens

The hackathon build shows the existing single-page Gradio app in the window (sections
4.1 to 4.3 are one scrolling page today). The wireframes show the target layout; differences
from today are listed under each screen. Window minimum size 900 x 600; below 1000 px wide the
two-column areas stack.

### 4.1 Home / drop zone

```
+-------------------------------------------------------------------------+
| COAL  codon alignment                              [Settings] [About]   |
+-------------------------------------------------------------------------+
|  Give COAL a protein alignment and the DNA of the same genes. It lines  |
|  the DNA up codon by codon. Nothing leaves this computer.               |
|                                                                         |
|  +-------------------------------+   +-------------------------------+  |
|  | 1. Protein alignment          |   | 2. Coding DNA (not aligned)   |  |
|  |                               |   |                               |  |
|  |   Drop a file here, or        |   |   Drop a file here, or        |  |
|  |   [ Choose file... ]          |   |   [ Choose file... ]          |  |
|  |   FASTA or CLUSTAL            |   |   FASTA                       |  |
|  +-------------------------------+   +-------------------------------+  |
|   OG0001.aln.fasta  20 seqs  [x]      OG0001.cds.fasta  20 seqs  [x]    |
|                                       [ Swap ]                          |
|                                                                         |
|  Try an example:  [Clean nuclear genes] [Mitochondrial (wrong code)]    |
|                   [Messy IDs] [One broken gene]                         |
|                                                                         |
|  > Options  (FASTA, genetic code 1, no columns removed)                 |
|                                                                         |
|                       [  Start alignment  ]                        |
+-------------------------------------------------------------------------+
```

- Drop zones are large targets (at least 200 x 120 px) and real buttons for keyboard users
  (Enter or Space opens the file dialog). Each has a visible label, not only placeholder text.
- After a file is chosen the zone shows name, sequence count and a remove button.
- Today: same two zones (Gradio `gr.File`), example buttons, intro; no Settings/About bar, no
  Swap, no sequence count, options always visible.

### 4.2 Options

Collapsed by default under one summary line, because defaults are right for most people.

```
|  v Options                                                              |
|  +-------------------------------------------------------------------+  |
|  | Genetic code     [ 1: Standard (nuclear genes of most organisms) v]|  |
|  |                  Code 1 for nuclear genes. Mitochondria need their |  |
|  |                  own code (2 for vertebrates).                     |  |
|  | Output format    [ FASTA v ]   Which file format you want back.    |  |
|  | [ ] Remove gap columns (-nogap)                                    |  |
|  |     Drop columns with a gap or stop codon in any sequence.         |  |
|  | [ ] Remove mismatched columns (-nomismatch)                        |  |
|  |     Drop columns where a codon doesn't match its amino acid.       |  |
|  | Save results in  (•) a folder next to the protein file             |  |
|  |                  ( ) another folder... [ Browse ]                  |  |
|  +-------------------------------------------------------------------+  |
```

- Each control's explanation is its accessible description (read by Narrator after the label).
- Today: genetic code, format, the two checkboxes, each with an explanation. No "Save results
  in" (results go to a temp folder and are downloaded).

### 4.3 Results

```
+-------------------------------------------------------------------------+
| ⚠ Done, with some problems                                [Open folder] |
| 19 of 20 sequences OK, 1 failed. Genetic code 1. 0.4 s.                 |
+-------------------------------------------------------------------------+
| ! Tip: ...these look like vertebrate mitochondrial genes...             |
|   [ Re-run with genetic code 2 ]                                        |
+-------------------------------------------------------------------------+
| Sequence health, problems first           [ Show: problems | all ]      |
|  Status      Sequence        What this means                            |
|  ✖ Failed    OG0001_sp14     From codon 102 the DNA no longer spells    |
|                              the protein: one letter missing...         |
|                              [Show in alignment]                        |
|  ✔ OK        OG0001_sp01     Every codon matches.                       |
|  ...                                         page 1 of 1  [TSV]         |
+-------------------------------------------------------------------------+
| v View alignment                                                        |
|   name        |  98   99  100  101  102  103 ...                        |
|   OG0001_sp01 |  GCT  AAA  TGG  CTG  GAT  ...                           |
|   OG0001_sp14 |  GCT  AAA  TGG  CTG ##G#A# ...  (frameshift)            |
|   [ Download alignment view (HTML) ]                                    |
+-------------------------------------------------------------------------+
| Downloads:  [Codon alignment]  [Report (JSON)]                          |
+-------------------------------------------------------------------------+
| Methods paragraph                                             [Copy]    |
|  Codon alignments were built with COAL 0.4.0 (genetic code 1, ...)      |
+-------------------------------------------------------------------------+
| Equivalent commands                                           [Copy]    |
|  coal OG0001.aln.fasta OG0001.cds.fasta --keep-going -o ...             |
|  pal2nal.pl OG0001.aln.fasta OG0001.cds.fasta -output fasta > ...       |
+-------------------------------------------------------------------------+
| > Raw COAL messages     > What do these words mean?                     |
+-------------------------------------------------------------------------+
```

- Status is always icon + word + sentence (Okabe-Ito colours are a fourth cue, never the only
  one). The banner is an ARIA live region so Narrator announces "Done, with some problems"
  without the person hunting for it.
- After a run, keyboard focus moves to the banner heading.
- Today: all blocks exist except Open folder, the problems/all toggle, paging and TSV.

### 4.4 Settings

Small, because every setting is a support question. Stored in
`%LOCALAPPDATA%\COAL\settings.json`.

```
+--------------------------- Settings ------------------------------------+
| Results folder     (•) next to the protein file                          |
|                    ( ) always in: [ C:\Users\priya\Documents\COAL ] [..] |
| Default genetic code   [ 1: Standard v ]                                 |
| Text size          follows Windows (Settings > Accessibility > Text size)|
| Appearance         follows Windows light / dark / contrast themes        |
| For AI agents      [x] Let AI tools on this computer use COAL (MCP)      |
|                        On by default. Only programs on this computer.    |
| Logs               [ Open log folder ]  (for bug reports; no sequences)  |
|                                                       [ Close ]          |
+--------------------------------------------------------------------------+
```

- No "check for updates", no telemetry switch (there is no telemetry to switch).
- Hackathon build: no Settings screen; the defaults above are fixed.

### 4.5 About and licences

```
+---------------------------- About COAL ---------------------------------+
| COAL 0.4.0 (engine: Python 0.4.0)          MIT licence                   |
| A fast, open reimplementation of PAL2NAL. Works offline. Nothing you     |
| open in COAL leaves this computer.                                       |
| How to cite: ... [Copy]                                                  |
| Updates: there is no automatic update. New versions are on               |
|          github.com/plumat-coal/COAL/releases (copy link).                |
| Third-party software included:  [ View licences ]                        |
|   Python (PSF), Gradio (Apache-2.0), pywebview (BSD-3-Clause),           |
|   Microsoft Edge WebView2 (part of Windows), PyInstaller bootloader      |
|   (GPL-2.0 with bootloader exception), and Gradio's dependencies         |
|   (THIRD_PARTY_LICENSES.txt in the install folder).                      |
| Built: 2026-10-17, commit a2351d0.  Unsigned test build.                 |
+--------------------------------------------------------------------------+
```

- The link is shown as text with a copy button rather than opened, so nothing reaches the
  internet unless the person chooses to paste it.
- Hackathon build: licences are files in the install folder (`licenses\`), listed in the
  README; no About screen.

---

## 5. Windows decisions

### 5.1 Wrapper technology: pywebview + Edge WebView2

**Decision:** pywebview (Python, BSD licence) with its `edgechromium` backend, which uses the
Microsoft Edge WebView2 runtime built into Windows 11 and installed on Windows 10 by Edge
updates. The window loads `http://127.0.0.1:<random port>/` served by the existing local app.
Pinned version: see `requirements.txt`.

| | pywebview (chosen) | Tauri | Electron | Just open the browser |
|---|---|---|---|---|
| Language of the wrapper | Python, same as COAL | Rust + JS | JS (Node) | none |
| Reuses the Gradio app as is | yes | yes, but Python still needs bundling as a "sidecar" | yes, Python as a child process | yes |
| Extra runtime shipped | none (WebView2 is in Windows) | none (WebView2) | a whole Chromium, ~150 MB | none |
| Toolchain the team already has | yes (pip) | no (Rust, cargo, MSVC) | no (Node, electron-builder) | yes |
| Native window, no browser tab | yes | yes | yes | no: a tab among 40, easy to close by accident, looks like a website (which undermines "nothing is uploaded") |
| Native file dialogs, drag-drop paths | yes (`create_file_dialog`, DOM `drop` event with full paths) | yes | yes | no paths, only uploads |
| Risk | Python GUI threading quirks, pythonnet dependency on Windows | two languages, sidecar packaging | size, Chromium security updates become our job | none technically, worst UX |

Reasons: one language (Python) for the team; nothing new to install for builders; the
smallest download; WebView2 is serviced by Microsoft Update, so we don't ship a browser we
must patch. Tauri would be the choice if we rewrote the GUI in JS later. Electron's only
advantage (identical Chromium everywhere) doesn't matter for one Windows target.

Fallback if WebView2 is missing or fails to start: a native message box explains it and offers
"Open COAL in your web browser instead" (still 127.0.0.1 only), plus the link to Microsoft's
WebView2 runtime installer as text. We do not bundle the WebView2 runtime for the hackathon
(pywebview 6.1+ can point to a bundled "fixed version" runtime via
`webview.settings['WEBVIEW2_RUNTIME_PATH']`, adding ~150 MB; a later option for locked-down
machines).

### 5.2 Process model and shutdown

```
COAL.exe (wrapper: window + single-instance lock)        main thread: pywebview window
  └─ COAL.exe --serve --port 53817   (local app, Gradio)  stdin = pipe from the wrapper
       └─ COAL.exe -m coal ...        (engine, one per run)
```

- One executable plays three roles chosen by its arguments, because in a PyInstaller build
  `sys.executable` is `COAL.exe`, not `python.exe`, and the local app starts the engine as
  `sys.executable -m coal`. The wrapper recognises `-m coal` and runs the engine.
- The local app runs in a child process, not a thread, so a crash or a hang in Gradio cannot
  freeze the window, and shutting it down is certain (terminate, then kill).
- **Lifeline:** the child watches its stdin; when the wrapper exits for any reason, including
  a crash or Task Manager "End task", the pipe closes and the child exits by itself.
- On close: the window's `closed` event stops the child in a background thread (close stdin,
  wait up to 5 s, terminate, wait 3 s, kill), so the window never waits for it. Closing during
  "Starting COAL…" is not an error. On Windows the child is also put in a Job Object with
  "kill on job close", so engine runs (grandchildren) die with it; if that fails, stopping
  uses `taskkill /T /F` on the child's process tree instead.
  <!-- CHANGED(WIN-W2, WIN-R9, WIN-R10). See overnight/windows-app/CHANGES.md -->
- **Readiness without guessing:** the child prints one line `COAL-READY http://127.0.0.1:N/`
  after the server is listening; the wrapper waits for that line with a timeout (180 s; first
  launch after install can be slow while Defender scans the new files, and after 10 s the
  splash says so) and shows a clear error with the log path if it doesn't come.
  <!-- CHANGED(WIN-R8): was 60 s. See overnight/windows-app/CHANGES.md -->
- **Port:** the wrapper asks Windows for a free port on 127.0.0.1 (bind to port 0), then
  starts the child on it; if the child fails to bind (another program took it in the gap),
  the wrapper tries a new port, up to 3 times.
- **Single instance:** a lock file in `%LOCALAPPDATA%\COAL\` held for the lifetime of the
  wrapper. A second launch finds the lock taken, connects to the first instance's control
  port (127.0.0.1, random, protected by a random token stored next to the lock) and asks it
  to come to the front (and, later, to open the files it was given by a file association).
  Then the second launch exits. Locks are released by Windows when a process dies, so a crash
  never leaves COAL "already running" forever (the release can lag, so a launch retries for
  3 s before saying COAL is running but not responding). Only the instance holding the lock
  writes or rotates logs. <!-- CHANGED(WIN-R5, WIN-W3). See overnight/windows-app/CHANGES.md -->

### 5.3 Security of the local server

- Binds to `127.0.0.1` only, never `0.0.0.0`: other computers on the network cannot reach it,
  and Windows Defender Firewall does not show its "allow access" prompt (that prompt appears
  only for programs listening on network interfaces), so no admin-looking dialog.
- Random port each launch, so nothing else on the machine can rely on finding it.
- MCP (AI-agent tools) **off** in the desktop window by default (IMP-M-3: on a random port nothing can register
  it); `COAL.exe --mcp` turns it on at `http://127.0.0.1:<port>/gradio_api/mcp/sse` (the URL is in
  `local_app.log`). For Claude, use the stdio server (`desktop\setup_mcp.bat`).
  <!-- CHANGED(WIN-M): was off by default with opt-in --mcp (this line said --no-mcp, the code
  had --mcp). User decision. See overnight/windows-app/CHANGES.md -->
- Known gap: any program running as any user on the same machine can talk to 127.0.0.1 ports;
  Gradio has no per-launch secret. On a shared lab PC another logged-in user could in
  principle call the app while it is open. Proposed fix (section 8): a per-launch token the
  wrapper passes to the child and puts in the window's URL. Acceptable for the hackathon.
- Proxies: the wrapper sets `NO_PROXY=127.0.0.1,localhost` for the child, because university
  machines often have a system proxy and the local app must never send its own traffic
  through it.
- Offline: `GRADIO_ANALYTICS_ENABLED=False`; the GUI uses system fonts only (no CDN). The
  smoke test checklist includes running with Wi-Fi off.

### 5.4 Installer: per-user, no admin rights (Inno Setup)

- Inno Setup 6 (free, mature, one script file, `PrivilegesRequired=lowest`), installing to
  `%LOCALAPPDATA%\Programs\COAL` (the per-user Program Files used by VS Code, Teams, etc.).
  No UAC prompt, works on managed university laptops where users are not admins.
- Alternative considered: MSIX. Cleaner install/uninstall, but needs signing even to test-
  install, and sideloading policies vary on managed machines. Later, once signed.
- Alternative considered: a zip of the PyInstaller folder. Kept as a second download for
  people whose IT blocks installers ("portable": unzip, double-click `COAL.exe`).
- The installer is one page with two optional ticks (desktop shortcut, file associations),
  both off by default.

### 5.5 Start menu, file associations, uninstall, updates

- **Start menu:** one shortcut "COAL" in the user's Start menu (`{userprograms}`), no folder
  full of "Uninstall COAL" and "Readme" shortcuts (uninstall lives in Settings > Apps, where
  Windows users look for it).
- **File associations (optional, off by default):** registered per user under
  `HKCU\Software\Classes`, and only as an *Open with* option (`OpenWithProgids`) for `.aln`,
  `.fasta`, `.fa`, `.fas`, `.faa`, `.fna`. We never take over the default program for
  `.fasta` (people open those in Notepad++, Jalview, MEGA, and would be rightly annoyed).
  Opening a file with COAL passes its path to COAL; if COAL is open, the running window
  receives it via the single-instance channel. COAL comes to the front and loads the file by
  path into the right box (sniffed: protein alignment or DNA) through the GUI's
  `?protein=&dna=` hook, keeping the other box's file. <!-- CHANGED(WIN-W4b/R2-A9): was "shows a
  message ... can't fill the box yet". See overnight/integration/CHANGES.md -->
- **Uninstall:** Settings > Apps > COAL > Uninstall. Removes the install folder, the shortcut(s),
  the per-user registry keys, and COAL's own app data (`%LOCALAPPDATA%\COAL`: logs, lock,
  WebView2 cache) and `%TEMP%\coal_app_runs` (the GUI's run folders). Never touches results
  the person saved. Setup and uninstall ask to close a running COAL using Inno Setup's
  `AppMutex` (the wrapper holds a named mutex `COAL-desktop-running`); they don't check the
  lock file. Upgrades delete the old `_internal` folder first.
  <!-- CHANGED(WIN-W6, WIN-W8): was "plus a check of the lock". See overnight/windows-app/CHANGES.md -->
- **Updates:** manual. No update check, no auto-updater, nothing phoning home. About shows the
  version and the releases address as copyable text. Installing a new version over an old one
  works (same AppId) and keeps settings.

### 5.6 Code signing and SmartScreen

**Hackathon: unsigned.** Signing needs an identity-validated certificate or a service account,
which takes days to weeks. What people will see, in order:

1. **Edge, after download** (sometimes): "COAL-Setup-0.4.0.exe isn't commonly downloaded.
   Make sure you trust COAL-Setup-0.4.0.exe before you open it." Click the **"…"** (More
   actions) next to the warning, then **"Keep"**; if asked again, **"Show more"**, then **"Keep
   anyway"**.
2. **Running the installer:** a blue box titled **"Windows protected your PC"**: "Microsoft
   Defender SmartScreen prevented an unrecognized app from starting. Running this app might
   put your PC at risk." Click **"More info"**: it shows "App: COAL-Setup-0.4.0.exe,
   Publisher: Unknown publisher". Click **"Run anyway"**.
3. No UAC prompt follows, because the installer is per-user.
4. **Windows 11 with Smart App Control on** (on by default only on fresh installs of 22H2 and
   later, and it switches itself off on many machines after an evaluation period): unsigned
   apps are blocked with "Smart App Control blocked an app that may be unsafe" and **there is
   no Run anyway**. The only fix is signing (or the person turning Smart App Control off,
   which can't be turned back on without reinstalling Windows; we do not recommend that).
   This is the strongest reason to sign before any public release.
5. Managed laptops with AppLocker/WDAC may block unsigned programs from `%LOCALAPPDATA%`;
   then IT must allow it. The README says so plainly.

**Later:** sign `COAL.exe`, the uninstaller and the installer with Microsoft's Trusted
Signing service (low monthly cost, identity validation required) or SignPath's free plan for
open-source projects. Note that since 2024 even EV certificates no longer give instant
SmartScreen reputation; reputation builds with downloads of signed files. Sign with a
timestamp so signatures outlive the certificate.

### 5.7 Antivirus false positives

PyInstaller programs are sometimes flagged by antivirus heuristics, because malware also uses
PyInstaller. Mitigations, in order of effect:

1. **One-folder build, not one-file.** One-file builds unpack themselves to `%TEMP%` at every
   start, which looks like a dropper and is slower. One-folder is also faster to start.
2. **No UPX compression** (`upx=False`), a common trigger.
3. **Build with an official python.org Python**, a recent PyInstaller, and keep build
   reproducible (pinned requirements).
4. **Sign the files** (section 5.6): the biggest single improvement.
5. If Defender still flags a release: submit the file at Microsoft's Security Intelligence
   file-submission page as a false positive (usually cleared in days). Other vendors: via
   VirusTotal's vendor contacts.
6. Later: build the PyInstaller bootloader from source (unique binary, fewer heuristic hits).

### 5.8 Long paths, OneDrive and spaces

- Install path is short (`%LOCALAPPDATA%\Programs\COAL`), well under Windows' classic 260
  character limit. Python 3.6+ supports long paths only if the system setting
  `LongPathsEnabled` is on, which we can't assume; so the app writes results with short names
  (`OG0001_coal\OG0001.codon.fasta`) and, if a result path would exceed 250 characters, falls
  back to `Documents\COAL results\` and says so.
- Data often lives in `C:\Users\Priya Sharma\OneDrive - University of X\Thesis data\`: spaces,
  a hyphen, and OneDrive "Files On-Demand" placeholders. Checked here: the engine handles a
  path with spaces and CRLF files correctly (test in `test_smoke.py`). The equivalent command
  quotes paths with spaces using double quotes on Windows (the GUI's existing `_q`), which
  works in both cmd and PowerShell for paths without `$` or backticks.
- OneDrive placeholders: reading a file that isn't downloaded yet makes OneDrive download it
  first, which can take minutes for big files. The app shows "Reading file…" rather than
  appearing frozen, and Cancel works during it (target).
- Writing results into a OneDrive folder is allowed (it syncs, which is often what people
  want). The output folder setting lets people pick a local folder instead.

### 5.9 CRLF line endings and other Windows file quirks

- Windows tools (Notepad, MEGA, Excel exports) write CRLF line endings. The engine reads text
  in universal-newline mode: tested here, CRLF input gives byte-identical output to LF input.
- Output files are written with LF endings (the engine's choice; every bioinformatics tool and
  modern Notepad read them). Not changed per platform, so results are identical on all
  systems and checksums in reports match across machines.
- **Found while writing this design:** a FASTA file saved with a UTF-8 byte-order mark (as
  old Notepad and some Windows tools do) makes the engine stop with
  `'ascii' codec can't decode byte 0xef in position 0`. That message is not plain English, and
  the BOM is harmless. Fixed in the engine (ENG-h skips a BOM); `test_smoke.py`
  `test_utf8_bom_is_accepted` is now a normal test (INT-2). <!-- CHANGED(R2-D8) -->

### 5.10 High-DPI scaling

- Laptops commonly run at 125 % to 200 % display scaling. WebView2 renders the page at the
  monitor's scale. Whether pywebview's Windows backend declares the process per-monitor
  DPI-aware (sharp text rather than a blurry stretched bitmap) is not confirmed from its docs;
  the PyInstaller spec therefore does not override the manifest, and this is a smoke-test
  check. To verify on hardware: dragging the window between a
  100 % and a 150 % monitor should re-render at the new scale (smoke test step).
- Window size is given in logical pixels, capped to the screen's working area so the window
  never opens bigger than a small laptop screen.
- The page uses relative units (the GUI already uses an 18 px base and `rem`), so Windows'
  "Text size" accessibility setting and Ctrl+plus/minus zoom both work. Zoom is left enabled
  in the window.

### 5.11 Light and dark mode

- WebView2 passes the Windows "app mode" (Settings > Personalization > Colors) to the page as
  `prefers-color-scheme`. Gradio 5 follows it automatically. Whether the native title bar follows
  the theme depends on pywebview's backend (not confirmed; smoke-test check).
- The status colours (Okabe-Ito) were chosen to stay distinguishable on dark backgrounds; the
  warning text colour needs checking in dark mode (smoke test step).

### 5.12 Keyboard use

- Everything reachable with Tab / Shift+Tab in reading order; Enter/Space activate; visible
  3 px focus outline (already in the GUI CSS).
- Shortcuts (target, implemented in the page so they work identically in the browser):
  Ctrl+O choose protein file, Ctrl+Shift+O choose DNA file, Ctrl+Enter run, Esc cancel a
  run, Ctrl+Shift+C copy methods paragraph. Shown in tooltips and the About screen.
- Window-level: Alt+F4 closes and shuts down cleanly; F11 is not used (no full screen
  needed).

### 5.13 Narrator (screen reader)

- WebView2 exposes the page through Windows UI Automation, so Narrator reads the Gradio page
  like a web page in Edge. The window title is "COAL codon alignment" so Alt+Tab and Narrator
  announce it.
- Requirements on the page (most already true): every input has a label; status is words, not
  colour; the health table has a caption and column headers; the result banner is a live
  region (target); after a run focus moves to the banner (target); the alignment viewer is a
  table with row and column headers and a text summary before it ("1 problem: frameshift in
  OG0001_sp14 at column 102"), because a 2,000-column grid is not usable by ear.
- Smoke test includes a Narrator pass of one run (Ctrl+Win+Enter starts Narrator).

### 5.14 Windows high-contrast mode (Contrast themes)

- WebView2 supports `forced-colors: active`: in a Contrast theme, the page's colours are
  replaced by the theme's few system colours. Consequences: the coloured left borders in the
  health table and the tinted cells in the alignment viewer may disappear. Because every
  status also has an icon and a word, nothing is lost; the viewer's mismatch outline and
  frameshift hatching need `forced-color-adjust` care so they remain visible as borders
  (proposed CSS in section 8).
- Smoke test: switch to "Night sky" contrast theme (Settings > Accessibility > Contrast
  themes) and check the failed row and the frameshift cells are still identifiable.

### 5.15 Temporary files and logs

- Gradio's own temp files (`%TEMP%\gradio`, shared with other Gradio apps) and the GUI's run
  folders (`%TEMP%\coal_app_runs`, also used by runs from AI agents over MCP) go to `%TEMP%`.
  <!-- CHANGED(R2-D8): was "MCP runs: %TEMP%\coal_*". See overnight/integration/CHANGES.md -->
- The wrapper's logs go to `%LOCALAPPDATA%\COAL\logs\`: `coal_desktop.log` (newest 5 kept)
  and `local_app.log` (newest 3), rotated only by the instance holding the lock. They never
  contain sequence data, but `local_app.log` holds everything the local app prints (including
  file paths in error messages); the wrapper itself logs dropped/opened files by name only.
- WebView2 profile: the wrapper calls `webview.start(private_mode=True,
  storage_path=%LOCALAPPDATA%\COAL\webview)`. Private mode means no cookies or local storage
  persist between launches (COAL needs none). Where pywebview then puts WebView2's data folder
  (that `storage_path`, or a temporary folder) is **not verified**: check on Windows. Either
  way, uninstall removes `%LOCALAPPDATA%\COAL`.
  <!-- CHANGED(WIN-W6): was "the profile goes to %LOCALAPPDATA%\COAL\webview" as a fact.
  See overnight/windows-app/CHANGES.md -->

---

## 6. Large data

### 6.1 Read from disk, never upload

- In a desktop app there is no reason to "upload": the engine can read the file where it is.
  Target: the wrapper's native file dialog (`create_file_dialog`) or the window-wide drop
  handler (pywebview's DOM `drop` event gives full paths) produces a path, which the GUI
  passes straight to the engine. Zero copies, works for files larger than RAM.
- Done (R2-A9): files dropped anywhere on the window, opened with COAL or given to a second
  launch are passed by path (`?protein=&dna=`) and read in place. Only the box's own file
  picker still uploads a copy into Gradio's temp folder. <!-- CHANGED(R2-A9) -->

### 6.2 Progress with time remaining

- The engine processes sequences one by one; it can report "done k of n" cheaply. Proposed:
  `--progress` makes it print `coal: progress k/n` to stderr at most every 0.5 s (section 8).
  The GUI shows "Checking 12,400 of 20,000 sequences, about 40 s left", estimating time left
  from the rate so far, and "Reading file…" while it counts sequences.
- Today: Gradio's built-in spinner/timer only.

### 6.3 Cancel

- Runs happen in an engine subprocess, so cancel is reliable: kill the process (and its tree,
  via the Job Object), delete the partial output, say "Cancelled. Nothing was written." Today
  the GUI calls the engine with a blocking `subprocess.run` and no cancel; proposed change in
  section 8 (`Popen` + Gradio's `cancels=`).

### 6.4 Where results go

- Default: a new folder next to the protein file, `<protein stem>_coal\`, containing the
  codon alignment, the report JSON, and the alignment view HTML; a second run makes a
  timestamped sub-folder instead of overwriting. If the input folder isn't writable (a USB
  stick, a read-only share), fall back to `Documents\COAL results\` and say so in the banner.
- Today: a temp folder plus Download buttons (fine for the web version, but on Windows the
  downloads end up in Downloads with no link to the input). Needs a GUI option (section 8).

### 6.5 Showing big results

- Health table: summary counts, then problems first, paged at 100 rows; full table as TSV.
  A 100,000-row HTML table would freeze the page.
- Viewer: already capped at 200 sequences x 2,000 columns and 200,000 cells with a note;
  target is to centre that window on the selected problem.
- Raw messages: last 4,000 characters shown, full log downloadable.

---

## 7. Scope

### Hackathon demo (this branch)

- `desktop/coal_desktop.py`: native window around the existing GUI; random port on 127.0.0.1;
  local app in a child process with lifeline and clean shutdown; single instance with
  bring-to-front; engine dispatch for frozen builds; offline; drag-and-drop paths captured
  where pywebview supports it (a window-wide drop guard, registered each time the app's page
  loads, prevents the window from navigating away to a dropped file; dropped file names are
  logged); the Gradio drop zones accept drops as they do in a browser.
  <!-- CHANGED(WIN-W1): the guard used to sit on the splash page only. See overnight/windows-app/CHANGES.md -->
- `desktop/build.ps1` + `coal_desktop.spec`: PyInstaller one-folder build bundling Python,
  COAL engine, GUI, examples, and `coal-cpp.exe` if it exists.
- `desktop/installer.iss`: per-user Inno Setup installer, Start menu shortcut, optional
  desktop shortcut and Open-with entries, clean uninstall.
- Unsigned. SmartScreen documented.
- `desktop/SMOKE_TEST.md` and `desktop/test_smoke.py`.
- The four examples run from the "Try an example" buttons.

### Later (in rough priority order)

1. Code signing (removes Smart App Control blocks and most SmartScreen friction).
2. Paths instead of uploads; native Save-as dialogs; results next to the input; Open folder.
3. Engine progress + Cancel; paged health table + TSV.
4. Engine fixes: UTF-8 BOM, ID normalisation re-run button.
5. Window-wide drop that sorts protein vs DNA automatically; Swap button.
6. Settings and About screens; THIRD_PARTY_LICENSES generated at build time.
7. Per-launch token on the local server.
8. Forced-colors CSS for the viewer; full Narrator pass with a screen-reader user.
9. Bundled WebView2 fixed runtime as an option for locked-down machines; MSIX.
10. Batch mode in the GUI (a folder of 20,000 orthogroups in one go, which is Priya's real
    task; today the CLI and `demo/race.py` cover it).

---

## 8. Changes this design needs outside `desktop/` (for a human to decide)

The desktop branch may only edit `desktop/` and `overnight/windows-app/`. These are proposals;
exact wording is also in `overnight/windows-app/LOG.md` under "Needs a human decision".

1. **Engine: accept a UTF-8 BOM.** In `src/coal/seqio.py`, open text files with
   `encoding="utf-8-sig"` and then check `seq.isascii()` (or decode the first line with
   `utf-8-sig` and continue in ASCII), so a BOM is skipped and any other non-ASCII byte still
   gives a clear error naming the file and line.
2. **GUI: examples folder in frozen builds.** `core.EXAMPLES_DIR` is computed as
   `<coal_app>/../../examples`, which doesn't exist inside a PyInstaller bundle. The wrapper
   currently sets `core.EXAMPLES_DIR` itself before building the UI (works, but reaches into
   a module). Proposed: `EXAMPLES_DIR = os.environ.get("COAL_EXAMPLES_DIR") or <current>`.
3. **GUI: path inputs.** A hidden `gr.Textbox` pair (protein path, DNA path) or a
   `gr.api` endpoint `load_paths`, so the wrapper's file dialog / window-wide drop can hand
   paths to the GUI without uploading.
4. **GUI: results next to the input** (output-folder option, Open folder button calling the
   wrapper through pywebview's `js_api`).
5. **Engine: `--progress`** lines on stderr; **GUI:** `Popen` streaming those lines to a
   progress text, and a Cancel button (`cancels=[run_event]`) that kills the engine.
6. **GUI: per-launch token.** Accept `--token` and reject requests without it (Gradio `auth`
   or a middleware); the wrapper passes it in the window URL.
7. **GUI: forced-colors CSS**, e.g.
   `@media (forced-colors: active) { .coal-health tr td:first-child { border-left-color: CanvasText; } }`
   and the equivalent for the viewer's mismatch/frameshift outlines.
8. **GUI: result banner as a live region** (`elem_id` + `aria-live="polite"`), and focus moved
   to it after a run.
