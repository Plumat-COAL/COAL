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

1. Download `COAL-Setup-0.4.1.exe`.
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
only allows signed apps, and this unsigned build can't run on it: Smart App Control has no way to
allow a single app. You can turn Smart App Control off (Windows Security > App & browser control),
but on some Windows versions it can then only be turned back on by resetting or reinstalling
Windows, so we suggest leaving it on and installing COAL with pip instead (see "Install with pip"
in the release notes or on the website).
<!-- CHANGED(FIN-C18): was "it can't be turned back on without reinstalling Windows; wait for a signed version"
(inconsistent with the site and the release notes, which recommend pip) -->

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
4. Click **Start alignment**.  <!-- CHANGED(IMP-M-4): button renamed in POL-7 -->

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
| a gene failed with a "frameshift" | a DNA letter is missing or extra in that gene (often a sequencing or assembly error) | fix or remove that gene and run again: COAL stops at the first failed gene, so no alignment is written: the genes before it are shown as "Not written" (with their own problems, if any) and the genes after it as "Not checked", not as fine | <!-- CHANGED(FIN-C7) -->
<!-- CHANGED(R2-D8): BOM row removed (the engine now skips a byte order mark, ENG-h); frameshift advice: later genes are "not checked". See overnight/integration/CHANGES.md -->

## Run from source (no installer)

<!-- CHANGED(IMP-M-2): new; the launchers are in this folder -->
With Python 3.10 or newer (3.12 recommended, from https://www.python.org/downloads/windows/) and a copy of the COAL
repository (`git clone` or the source zip):

- **The window:** double-click `desktop\run_from_source.bat`. The first start creates `.venv-gui` in the repository
  folder and installs COAL with its GUI (a few minutes); later starts open the window straight away. If an install
  was interrupted, the next start notices (it tries to import COAL, Gradio and pywebview) and repairs it.
- **COAL in Claude (MCP):** double-click `desktop\setup_mcp.bat`. It creates `.venv-mcp` (a separate environment:
  the GUI's and the MCP server's packages conflict), checks that the server imports, and prints a
  `claude mcp add coal -- ...` line for Claude Code and a ready-to-paste `claude_desktop_config.json` entry for
  Claude Desktop. Run it again to repair the install or to see the lines again. Give Claude full paths to your files.

Both use the Python launcher (`py -3.12`, then `py -3`), or `python` if there is no launcher, and stop with a message
if Python is missing or older than 3.10.

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
- **AI tools (MCP)**: the window does **not** serve COAL's MCP tools by default, because its port changes every
  launch and an MCP client could not be pointed at it. To use COAL from Claude, set up the stand-alone MCP server
  once with `desktop\setup_mcp.bat` (see "Run from source" above). `COAL.exe --mcp` still turns the endpoint on
  for a session, at `http://127.0.0.1:<port>/gradio_api/mcp/sse` (the port is written in
  `%LOCALAPPDATA%\COAL\logs\local_app.log`).
  <!-- CHANGED(IMP-M-3): was on by default (WIN-M) on a random port. See overnight/integration/CHANGES.md -->
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
