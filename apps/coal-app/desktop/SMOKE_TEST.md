# COAL for Windows: manual smoke test

For a teammate with a **clean Windows 10 or 11 laptop** (or a fresh VM, or Windows Sandbox).
About 30 minutes, plus 10 for the optional accessibility pass. Tick each box; for anything that
fails, write what you saw in the "Notes" column and keep going.

You need: `COAL-Setup-<version>.exe` (from `apps\coal-app\desktop\build.ps1`, see `BUILD.md`), sent to you
**through a browser download** (email, Teams, a GitHub release, or a shared drive opened in
Edge), not copied from a USB stick: the browser download is what triggers SmartScreen, which is
part of what we're testing. No admin account needed; please test as a normal user.

The automated part (`apps\coal-app\desktop\test_smoke.py`) checks that every example gives the right result
through the local app **only if it ran with Gradio installed** (`build.ps1` does this, and runs it
again against `COAL.exe`): look for `RealLocalApp ... ok` in the build output, not `skipped`.
Without Gradio it checks the same route with a stand-in for Gradio's HTTP layer. Either way, this
checklist is for what only a person can see.

| # | Step | Expected | ✓ | Notes |
|---|---|---|---|---|
| 1 | **Prepare.** Note Windows version (Win+R, `winver`), display scaling (Settings > System > Display), whether Smart App Control is on (Windows Security > App & browser control). Make a folder `Documents\COAL test data` and copy the four `examples\*` folders into it (from the repo zip). | You know your setup. | ☐ | Win __ build __, scaling __ %, SAC on/off |
| 2 | **Download and SmartScreen.** Download the installer in Edge. If Edge says "isn't commonly downloaded": "…" > **Keep** > **Show more** > **Keep anyway**. Double-click it. | Blue box "Windows protected your PC". **More info** shows "Publisher: Unknown publisher". **Run anyway** starts setup. *If Smart App Control blocks it with no Run anyway: note it, that's a known limit of unsigned builds; stop here.* | ☐ | Which warnings appeared? |
| 3 | **Install.** | **No UAC (admin) prompt.** One page with two unticked options (desktop shortcut, "Open with"). Tick both. Install takes under a minute. "Open COAL now" ticked on the last page. Start menu has one "COAL" entry with the COAL mark icon (squares, charcoal to orange, on a light tile; a 3 x 3 block at small sizes). | ☐ | |
| 4 | **First launch.** Turn **Wi-Fi off** first. Click Finish (or Start > COAL). | A window titled "COAL codon alignment" shows "Starting COAL…", then the app within ~10 s (first launch may be slower). **No browser tab**, no console window, **no Windows Firewall prompt**. Start COAL again from the Start menu: the same window comes to the front, no second window. | ☐ | Seconds to ready: __ |
| 5 | **Run each example** (buttons under "Try an example"). | **Clean nuclear genes:** "✔ All 8 genes match codon for codon". **Mitochondrial (wrong code):** "ℹ Wrong genetic code?" with code 1 vs code 2 in one sentence and the button **Re-run with genetic code 2** in the banner; the 24 warnings are one closed row; the button gives "✔ All 24 genes match codon for codon". **Messy IDs:** explains that protein IDs end in `.1` but DNA IDs don't. **One broken gene:** `OG0001_sp14` shown as ✖ Failed with a sentence about a missing letter (frameshift) near codon 102. Each message is understandable without a manual. | ☐ | Anything confusing? |
| 6 | **Your own files.** In Explorer, open `Documents\COAL test data\clean`. Drag `protein.aln.fasta` onto box 1 and `dna.fasta` onto box 2, press **Start alignment**. Then drag a file onto the page *outside* the boxes. Also try right-click `protein.aln.fasta` > Open with > COAL. | Drops into the boxes work; result as in step 5. A drop outside the boxes puts the file into the right box (its path appears under the box; the file is read where it is, not copied) and does **not** replace the app with the file's text (also after the app has finished loading). "Open with" brings COAL to the front and puts the file into the right box the same way; open `dna.fasta` with COAL too and both boxes are filled. Do it once with COAL closed and once with it open. Selecting two DNA files and choosing Open with > COAL says only one DNA file can be loaded at a time and names the one used. If you have OneDrive, repeat the first drag from a OneDrive folder with spaces in its name. | ☐ | |
| 7 | **View the alignment.** After "One broken gene", open "View alignment" and click the jump link in the failed row. | The viewer jumps to column ~102; the broken cells are outlined/hatched and readable; names stay visible while scrolling sideways; text is sharp (not blurry) at your scaling. | ☐ | |
| 8 | **Export.** Click **Download codon alignment**, **Download report (JSON)**, **Download alignment view (HTML)**. | Each asks where to save or lands in Downloads (note which). The FASTA opens in Notepad and looks like FASTA; the HTML view opens in Edge, also with Wi-Fi off. | ☐ | Where did files go? |
| 9 | **Copy the methods sentence and the command.** Click the copy icon on "Methods paragraph", paste into Word or Notepad; same for "Equivalent commands". | The paragraph names COAL's version, the genetic code and counts, and contains **no** `C:\Users\...` path. The command uses file names only. | ☐ | |
| 9b | **AI tools (MCP).** Open `%LOCALAPPDATA%\COAL\logs\local_app.log` and check there is **no** `MCP tools for AI agents` line (off by default, IMP-M-3). Quit COAL, start `COAL.exe --mcp` from the install folder, find the line `MCP tools for AI agents at http://127.0.0.1:…/gradio_api/mcp/sse` and open that address in Edge. | No MCP line by default. With `--mcp` the address answers (an event stream that stays open; close the tab). | ☐ | |
| 9c | **Close while starting.** Quit COAL, start it, and close the window while it still says "Starting COAL…". | The window closes; **no** "COAL could not start" message; no `COAL.exe` left in Task Manager after ~10 s. | ☐ | |
| 10 | **Close and uninstall.** Close the window (X, or Alt+F4). Open Task Manager (Ctrl+Shift+Esc) > Details: no `COAL.exe` left after ~10 s. Then Settings > Apps > Installed apps > COAL > **Uninstall**. | Window closes at once, no leftover processes. Uninstall asks no admin rights and finishes. Afterwards: no Start menu or desktop entry; `%LOCALAPPDATA%\Programs\COAL`, `%LOCALAPPDATA%\COAL` and `%TEMP%\coal_app_runs` are gone; right-click a .fasta > Open with no longer lists COAL; your files in `Documents\COAL test data` and in Downloads are untouched. | ☐ | |

Before uninstalling, if anything failed: zip `%LOCALAPPDATA%\COAL\logs` (no sequence data in
there, but file names and possibly full paths) and send it with this sheet.

## Frozen-build checks (for whoever made the build; see BUILD.md)

- [ ] Hidden `COAL.exe --serve` child gets working stdout/stdin pipes (if not: "stopped while starting (exit code 4)").
- [ ] `COAL.exe -m coal --version` prints the version (build.ps1 checks it).
- [ ] A large dataset starts worker processes (spawn via `freeze_support`) with no extra windows, and they exit afterwards.
- [ ] A forced error lands in `local_app.log`, not in an "Unhandled exception" box.
- [ ] Where WebView2 put its profile (`%LOCALAPPDATA%\COAL\webview` or a temp folder): note it in LOG.md.

## Optional: accessibility pass (10 minutes, before step 10)

| # | Check | Expected | ✓ |
|---|---|---|---|
| A1 | **Keyboard only.** Unplug or ignore the mouse. Tab through the page, run "Clean nuclear genes" with Enter, reach the copy buttons. | Every control reachable, visible focus outline, nothing traps focus. | ☐ |
| A2 | **Narrator** (Ctrl+Win+Enter). Run "One broken gene" and listen. | Labels of the two file boxes and options are read; the result banner and the failed row can be found and are read as words ("Failed"), not only colours. | ☐ |
| A3 | **Scaling.** Set display scaling to 150 % (or move the window to a monitor with different scaling). | Text re-renders sharp; nothing cut off at 900 x 600 window size. Ctrl+plus / Ctrl+minus zoom works. | ☐ |
| A4 | **Dark mode.** Settings > Personalization > Colors > Choose your mode: Dark. Restart COAL. | Page is dark; warning and failed colours still readable. Title bar: note whether it turns dark. | ☐ |
| A5 | **Contrast theme.** Settings > Accessibility > Contrast themes > Night sky. Run "One broken gene". | Status is still identifiable by icon and word; note whether the viewer's broken cells are still visible. | ☐ |

## Result

Tester: ________  Date: ________  Build/version: ________

Overall: ☐ ready for the demo  ☐ ready with notes  ☐ blocked (by step __)
