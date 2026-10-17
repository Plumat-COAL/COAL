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

1. Download `COAL-Setup-0.3.0.exe`.
2. Edge may say the file **"isn't commonly downloaded"**. Click the **…** next to the message,
   then **Keep**, then **Show more** and **Keep anyway**.
3. Double-click the file. Windows shows a blue box, **"Windows protected your PC"**. Click
   **More info**, then **Run anyway**.
   *Why:* this test version of COAL isn't digitally signed yet, so Windows doesn't recognise
   the publisher. It will say "Publisher: Unknown publisher". That is expected for this build.
4. Choose whether you want a desktop shortcut and an "Open with" entry for .fasta and .aln
   files (both optional; COAL never becomes your default program for those files). Click
   **Install**.

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
2. Drop your **protein alignment** (FASTA or CLUSTAL, from MAFFT, MUSCLE, Clustal…) on box 1,
   and the **coding DNA** of the same genes (FASTA, not aligned) on box 2. Or click a box to
   choose a file.
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
| a gene failed with a "frameshift" | a DNA letter is missing or extra in that gene (often a sequencing or assembly error) | fix or remove that gene; the others are fine |
| "ascii codec can't decode byte 0xef" | the file was saved with a hidden "byte order mark" (old Notepad does this) | open it in Notepad, File > Save as, Encoding: **UTF-8** (not "UTF-8 with BOM"), save, try again |

## Uninstall

Settings > Apps > Installed apps > **COAL** > **Uninstall**. This removes COAL, its shortcuts,
its "Open with" entries and its own settings and logs. Your sequence files and anything you
downloaded from COAL are not touched.

## Updates

COAL never checks for updates by itself. New versions are published at
`github.com/chloe-newt/COAL/releases`. Installing a new version over the old one is fine.

## Privacy and what's inside

- COAL runs a small web page **on your own computer only** (address 127.0.0.1, a different
  random port each time) and shows it in a window, using the Microsoft Edge WebView2 component
  that comes with Windows. Other computers can't reach it, and nothing is sent to the internet.
- Logs (start-up messages only, never your sequences) are in `%LOCALAPPDATA%\COAL\logs`. If
  something goes wrong, sending that folder to the COAL team helps.
- COAL includes Python, Gradio (the page), pywebview (the window) and their dependencies; their
  licences are in the install folder.

## Help

Open an issue at `github.com/chloe-newt/COAL/issues`, with the log folder zipped and, if
possible, a small file that shows the problem.
