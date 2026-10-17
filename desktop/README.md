# COAL for Windows

COAL turns a **protein alignment** and the **DNA of the same genes** into a **codon
alignment**: the DNA lined up three letters at a time, so that each codon sits under its amino
acid. You need one before measuring selection (dN/dS with PAML or HyPhy) or building trees from
DNA. COAL does what PAL2NAL does, but runs on Windows with no terminal, and tells you in plain
words when a gene doesn't fit.

- **Nothing leaves your computer.** COAL works offline. It doesn't upload your files, check for
  updates or collect statistics.
- **No administrator rights needed.** It installs just for you.
- **Free and open source** (MIT licence).

*(Screenshots to be added from the first Windows build.)*

## Install

1. Download `COAL-Setup-0.4.0.exe`.
2. Edge may say the file **"isn't commonly downloaded"**. Click the **…** next to the message,
   then **Keep**, then **Show more** and **Keep anyway**.
3. Double-click the file. Windows shows a blue box, **"Windows protected your PC"**. Click
   **More info**, then **Run anyway**.
   *Why:* this test version of COAL isn't digitally signed yet, so Windows doesn't recognise
   the publisher. It will say "Publisher: Unknown publisher". That is expected for this build.
4. Choose whether you want a desktop shortcut and an "Open with" entry for .fasta and .aln
   files (both optional; COAL never becomes your default program for those files). Click
   **Install**.
   <!-- CHANGED(WIN-W4b/R2-A9): was "it can't load the file into the box by itself yet"; it now does. See overnight/integration/CHANGES.md -->
   *About "Open with":* it opens COAL (or brings it to the front) and puts the file into the
   right box (COAL looks inside to tell a protein alignment from DNA). The file is read where it
   is, never copied. Open the other file the same way, or type its path.

COAL now appears in your Start menu.

**If Windows says "Smart App Control blocked an app"** with no "Run anyway" button, your PC
only allows signed apps, and this test build can't run on it yet. Please don't turn Smart App
Control off for COAL (it can't be turned back on without reinstalling Windows); wait for a
signed version.

**On a university or company laptop** the install may be blocked by IT policy. Send IT this
page; COAL installs in your own user folder (`%LOCALAPPDATA%\Programs\COAL`) and only listens
on this computer (127.0.0.1).

## Use

1. Start **COAL** from the Start menu. A window opens (the first start can take ten seconds or
   so).
2. Drop your **protein alignment** (FASTA or CLUSTAL, from MAFFT, MUSCLE, Clustal…) and the
   **coding DNA** of the same genes (FASTA, not aligned) anywhere on the window: COAL puts each
   one in the right box and reads it where it is (no copy, so big files are fine). You can also
   click a box to choose a file, or type or paste a file's path under the box.
   <!-- CHANGED(R2-A9): drops anywhere are now loaded by path. See overnight/integration/CHANGES.md -->
3. Check the **genetic code**: 1 for nuclear genes of most organisms, 2 for vertebrate
   mitochondria, 11 for bacteria and plastids.
4. Click **Make codon alignment**.

Not sure what your files should look like? Click one of the **Try an example** buttons.

### Reading the results

- The **banner** says whether everything matched.
- The **health table** lists problem genes first, each with a sentence explaining what's wrong
  (✔ OK, ⚠ Warning, ✖ Failed: always a symbol and a word, never just a colour).
- A **Tip** box appears when COAL can guess the fix, for example "these look like
  mitochondrial genes", with a button to re-run with the suggested setting.
- **View alignment** shows exactly where a gene breaks.
- **Download** buttons save the codon alignment, a report (JSON, a record of exactly what was
  run, including checksums of your files) and the alignment view.
- **Methods paragraph**: a ready-to-paste sentence for your paper, with a copy button.
- **Equivalent commands**: the same run as a command line, for a bioinformatician who wants to
  repeat it on a server (with COAL or with PAL2NAL).

### Common problems

| COAL says | What it means | What to do |
|---|---|---|
| many codons "don't match", tip suggests another genetic code | the genes use a different genetic code (common for mitochondrial genes) | click the re-run button in the tip |
| protein and DNA IDs don't match | e.g. proteins called `AT1G01010.1`, DNA called `AT1G01010` | rename so the IDs match (the tip says how they differ) |
| a gene failed with a "frameshift" | a DNA letter is missing or extra in that gene (often a sequencing or assembly error) | fix or remove that gene and run again: COAL stops at the first failed gene, so the genes after it are shown as "not checked", not as fine |
<!-- CHANGED(R2-D8): BOM row removed (the engine now skips a byte order mark, ENG-h); frameshift advice: later genes are "not checked". See overnight/integration/CHANGES.md -->

## Uninstall

Settings > Apps > Installed apps > **COAL** > **Uninstall**. This removes COAL, its shortcuts,
its "Open with" entries, its own settings and logs (`%LOCALAPPDATA%\COAL`) and its run folders
(`%TEMP%\coal_app_runs`). Your sequence files and anything you downloaded from COAL are not
touched.

Two kinds of temporary files are left for Windows' own cleanup (Settings > System > Storage >
Temporary files), or you can delete them yourself: `%TEMP%\gradio` (uploaded copies of files you
chose with a box's file picker; shared with any other Gradio-based program, so the uninstaller
leaves it). Runs started by AI tools over MCP use the same `%TEMP%\coal_app_runs` folder as the
window. <!-- CHANGED(R2-D8): was "%TEMP%\coal_* folders (MCP runs)"; MCP runs go to coal_app_runs (core.RUNS_DIR) -->
<!-- CHANGED(WIN-W6): temp folders listed. See overnight/windows-app/CHANGES.md -->

## Updates

COAL never checks for updates by itself. New versions are published at
`github.com/plumat-coal/COAL/releases`. Installing a new version over the old one is fine.

## Privacy and what's inside

- COAL runs a small web page **on your own computer only** (address 127.0.0.1, a different
  random port each time) and shows it in a window, using the Microsoft Edge WebView2 component
  that comes with Windows. Other computers can't reach it, and nothing is sent to the internet.
- **But any program, and any user logged in to this same computer, can use it while COAL is
  open**: the page has no password, and it can read any file your account can read if given its
  path. On a shared lab PC, close COAL when you're done.
  <!-- CHANGED(WIN-W6, decision 3): security line made plain. See overnight/windows-app/CHANGES.md -->
- **AI tools on this computer** (MCP clients) can use COAL's tools while it is open, at
  `http://127.0.0.1:<port>/gradio_api/mcp/sse`; the port changes every launch and is written in
  `%LOCALAPPDATA%\COAL\logs\local_app.log`. To turn this off, start COAL with `COAL.exe --no-mcp`.
  <!-- CHANGED(WIN-M): MCP on by default. See overnight/windows-app/CHANGES.md -->
- Logs are in `%LOCALAPPDATA%\COAL\logs`. They never contain your sequences, but they do contain
  file names, and `local_app.log` holds everything COAL's local app prints, which can include
  full file paths (for example in error messages). Look through them before sending; if
  something goes wrong, sending that folder to the COAL team helps.
  <!-- CHANGED(WIN-W6): was "start-up messages only". See overnight/windows-app/CHANGES.md -->
- COAL includes Python, Gradio (the page), pywebview (the window) and their dependencies; their
  licences are in the install folder.

## Help

Open an issue at `github.com/plumat-coal/COAL/issues`, with the log folder zipped and, if
possible, a small file that shows the problem.
