# Overnight C (Windows desktop app): morning summary

All 5 tasks done on branch `keanu_overnight/windows-app`, **local commits only, not pushed**
(per the launching session). Written on a **Linux cloud machine, not Windows**: no PyPI, so
pywebview, Gradio, PyInstaller and Inno Setup were never installed or run. Everything that needs
them is written from their documentation and marked untested below.

## Done
| Task | Commit |
|---|---|
| 1. `desktop/DESIGN.md`: personas, 6 flows, 5 wireframed screens, Windows decisions, large data, scope, proposed changes | 6cefb22 |
| 2. `desktop/coal_desktop.py` wrapper + `requirements.txt` + `test_wrapper.py` (23 tests) | 8bc2308 |
| 3. `coal_desktop.spec`, `build.ps1`, `installer.iss`, `BUILD.md`, `coal.ico` (+ `make_icon.py`), `.gitignore` (not built) | 8e1dcde |
| 4. `SMOKE_TEST.md` (10 steps + accessibility pass) + `test_smoke.py` | da4748d |
| 5. `desktop/README.md` for users (no screenshots: no GUI here) + this summary | 5ab0470 (+ this hash fix) |

## Tested on this machine (Linux, Python 3.13)
- `python -m unittest desktop/test_wrapper.py`: 23 OK. Port picking; local app child start ->
  `COAL-READY` line -> HTTP 200 -> clean stop (exit 0) and port freed; retry when the port is
  taken, give up after 3; start timeout kills the child; crash and missing-Gradio messages;
  **lifeline**: window process SIGKILLed -> local app exits by itself; single instance (lock,
  hand-over of files, wrong token refused, lock freed when the holder is killed); dropped-file
  sniffing (protein vs DNA, CLUSTAL, BOM, CRLF); `coal_desktop.py -m coal` byte-identical to
  `python -m coal`.
- `python -m unittest desktop/test_smoke.py`: 6 OK, 1 expected failure (UTF-8 BOM, engine),
  4 skipped (need Gradio / pywebview). All four examples (5 runs) go through the GUI's logic with
  the engine started the frozen-build way, and through a real child process started by the
  wrapper with a stdlib stand-in for Gradio's HTTP API.
- Engine tests 11 OK, app tests 52 OK (2 skipped, Gradio) unchanged.
- Spec file's collection logic dry-run with a fake PyInstaller (examples, sample report,
  LICENSE, version). Icon checked visually.

## Needs a Windows machine (untested)
- Whole window part of `coal_desktop.py`: `run_window` (pywebview window, splash -> `load_url`,
  `closed` event, document drop handler with `pywebviewFullPath`, bring-to-front,
  `ALLOW_DOWNLOADS`, `gui="edgechromium"`), `_run_in_browser`, `_error_dialog` (MessageBoxW).
- `launch_gradio` against real Gradio (`ui.launch(prevent_thread_lock=True)`, `demo.close`).
- Windows-only code: msvcrt lock, Job Object (ctypes), CREATE_NO_WINDOW, parent-handle
  lifeline, named mutex for Inno's AppMutex.
- `build.ps1` (never parsed: no PowerShell here), the spec (never run), `installer.iss` (never
  compiled), and the frozen build itself (`sys._MEIPASS` examples, `freeze_support`, a windowed
  exe's stdout/stdin pipes in the `--serve` child).
- `test_smoke.py` RealLocalApp (real Gradio HTTP API) and WindowApi; `SMOKE_TEST.md` by a person.

## Parked
- **No PyPI / no Windows tools here** (pip download of pywebview: DNS failure, tried once as
  instructed). Version pins `pywebview==6.2.1`, `pyinstaller==6.22.3` come from the projects'
  GitHub release pages, not PyPI.

## Needs a human decision
1. **Engine: UTF-8 BOM.** A FASTA saved with a BOM stops the engine with
   "'ascii' codec can't decode byte 0xef in position 0". Proposed: in `src/coal/seqio.py` read
   text with `encoding="utf-8-sig"` and keep rejecting other non-ASCII bytes with a message
   naming file and line. `test_smoke.py::test_utf8_bom_is_accepted` is an expected failure
   that will flip when fixed (then remove the decorator).
2. **GUI: examples in frozen builds.** `core.EXAMPLES_DIR` = `<coal_app>/../../examples` doesn't
   exist in a bundle; the wrapper overwrites `core.EXAMPLES_DIR` itself. Proposed in
   `src/coal_app/core.py`: `EXAMPLES_DIR = os.environ.get("COAL_EXAMPLES_DIR") or _find_examples()`.
3. **GUI: path inputs** (read from disk, no upload; needed for 100,000-sequence files and for
   "Open with" / window-wide drop to fill the boxes): a `gr.api` endpoint or hidden textboxes
   taking local paths. Until then, "Open with" only brings COAL to the front.
4. **GUI: results next to the input** + "Open folder"; **engine `--progress`** + GUI Cancel
   (`Popen` instead of `subprocess.run`); paged health table (DESIGN.md section 8).
5. **Security: per-launch token** for the local server (any local user can reach a 127.0.0.1
   port while COAL is open). Acceptable for the demo?
6. **MCP is off in the desktop app** (`--no-mcp` equivalent; `COAL.exe --mcp` turns it on).
   OK, or should the desktop app serve MCP for local AI agents by default?
7. **Code signing** before any public release (Smart App Control blocks unsigned apps with no
   override). Options: Microsoft Trusted Signing or SignPath (free for open source). Needs
   someone to own the identity validation.
8. **C++ engine**: the spec bundles `coal-cpp.exe` if built, but nothing calls it; the GUI
   always runs `python -m coal`. Decide how the GUI should choose it.
9. **Gradio pin**: still unverified (from overnight/app). Check `gradio[mcp]==5.49.1` installs
   with Python 3.12 on Windows before the first build.

## How to build and try it (Windows, PowerShell, repo root)
```powershell
powershell -ExecutionPolicy Bypass -File desktop\build.ps1     # venv, tests, PyInstaller, installer if Inno Setup 6 is installed
desktop\dist\COAL\COAL.exe                                     # portable app
& "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe" /DAppVersion=0.3.0 desktop\installer.iss   # installer by hand
```
From source without building: `py -3.12 -m venv desktop\.venv`,
`desktop\.venv\Scripts\python -m pip install -e ".[gui]" -r desktop\requirements.txt`,
`desktop\.venv\Scripts\python desktop\coal_desktop.py`. Then follow `desktop\SMOKE_TEST.md`.
On any OS: `python -m unittest desktop/test_wrapper.py desktop/test_smoke.py`.

---

## Progress log

- Sat 02:38 start. Machine: Linux cloud sandbox (not Windows), Python 3.13, on branch keanu_overnight/windows-app at a2351d0. Brief adapted by the launching session: skip the Windows git step; PyPI unreachable so pywebview/PyInstaller/Inno Setup can't be installed or run; commit locally only, never push.
- Sat 02:38 Note: the system Python has an editable install of coal/coal_app pointing at /home/claude/coal-app (another session's folder). Every command here sets PYTHONPATH=/home/claude/coal-win/src so this worktree's code is used; coal_desktop.py also puts ../src first on sys.path and PYTHONPATH in a source checkout.
- Sat 02:38 Problem "install pywebview": attempt 1 `pip download pywebview` (timeout 40) -> DNS failure for pypi.org. Not retried (instruction: don't try more than once). PyInstaller and Inno Setup not attempted for the same reason. Gradio is not installed either (known from overnight/app/LOG.md).
- Sat 02:38 pywebview docs were reachable instead (pywebview.flowrl.com guide/api and guide/dom, github examples/drag_drop.py, github releases page): latest release 6.2.1; DOM `drop` event gives full paths in `e['dataTransfer']['files'][i]['pywebviewFullPath']` via `window.dom.document.events.drop += DOMEventHandler(fn, prevent_default, stop_propagation)`; window events closing/closed/loaded/shown; methods restore/show/evaluate_js/destroy; `on_top` is a property; webview.start(func, args, gui, debug, private_mode, storage_path, icon, ...). No platform limitation for drop is documented. The installed-version check the brief asks for is impossible here; the pin is 6.2.1 from the releases page.
- Sat 02:40 Finding (tested here): engine reads CRLF files and paths with spaces correctly (output byte-identical to LF input), but a UTF-8 BOM makes it stop with "'ascii' codec can't decode byte 0xef in position 0". See Needs a human decision.
- Sat 02:41 Task 1 done: desktop/DESIGN.md (3 personas with moment-of-pain tables, 6 flows, 5 wireframed screens, 15 Windows decisions incl. wrapper comparison, process model, SmartScreen wording, Smart App Control, AV, OneDrive/long paths, CRLF/BOM, DPI, dark/contrast, keyboard/Narrator; large data; scope; 8 proposed changes outside desktop/).
- Sat 02:47 Task 2: first test run 22/23 OK; the failure was the test itself (binding the freed port hit TIME_WAIT from the HTTP request); fixed with SO_REUSEADDR (still fails if a listener remained). Also closed pipes of failed port attempts (ResourceWarning). Rerun with ResourceWarning as error: 23/23 OK.
- Sat 02:47 Task 2 done: desktop/coal_desktop.py + desktop/requirements.txt (pywebview==6.2.1, pyinstaller==6.22.3, unverified against PyPI) + desktop/test_wrapper.py (23 tests, all with timeouts, 5 s total here). TESTED HERE (Linux): free-port picking; LocalServer start/READY line/HTTP 200/clean stop and port release; retry when the port is taken (and giving up after 3); start timeout kills the child; crash and missing-Gradio messages; lifeline (window process SIGKILLed -> local app exits by itself); serve() protocol; NO_PROXY/analytics env; single instance (lock, hand-over message with files, wrong token refused, lock freed when the holder is killed); dropped-file sniffing (protein vs DNA, CLUSTAL, BOM, CRLF); engine dispatch `coal_desktop.py -m coal` byte-identical to `python -m coal` on 3 examples. NOT TESTED HERE: everything in run_window()/_run_in_browser()/_error_dialog() (pywebview window, splash -> load_url, closed event, drop handler, bring-to-front, ALLOW_DOWNLOADS, edgechromium-only, WebView2-missing dialog); launch_gradio() against real Gradio (Gradio not installable); Windows-only code: msvcrt lock, Job Object (ctypes), CREATE_NO_WINDOW, parent-handle lifeline, MessageBoxW; frozen-mode paths (sys._MEIPASS examples, freeze_support).
- Sat 02:47 Drag and drop decision: Gradio's drop zones keep accepting drops (in WebView2 they work like in Edge, to be confirmed on Windows). A document-level pywebview drop handler (documented API, prevent_default) stops the window navigating away to a file dropped elsewhere and logs the full paths, sorted protein/DNA. Handing those paths to the GUI needs a GUI hook (Needs a human decision 3).
- Sat 02:51 Task 3: PyInstaller, Inno Setup and PowerShell are not available here (no pwsh, no wine, no PyPI); per instructions, not attempted. The spec's own logic was dry-run with a fake PyInstaller module: version 0.3.0 parsed; datas = sample_report.json, LICENSE, the 4 examples (README + 2 FASTA each) under examples/<name>; no cpp/build so no coal-cpp.exe (message printed). build.ps1 not parsed (no PowerShell here).
- Sat 02:51 Task 3 done (written, NOT built): desktop/coal_desktop.spec (one-folder, windowed, upx off, version info, icon, collect_all for gradio/gradio_client/safehttpx/groovy/uvicorn/fastapi/starlette/webview, gradio kept as .py), desktop/build.ps1 (venv in desktop\.venv, every step with a timeout via Start-Process + WaitForExit, tests, PyInstaller, frozen checks incl. smoke test against COAL.exe, Inno Setup if found else prints the ISCC command), desktop/installer.iss (per-user, PrivilegesRequired=lowest, {autopf}\COAL, one Start menu shortcut, optional desktop shortcut and Open-with entries under HKCU only, AppMutex so setup/uninstall ask to close COAL, WebView2 runtime check with Microsoft's documented registry key, uninstall removes %LOCALAPPDATA%\COAL), desktop/BUILD.md (exact Windows commands), desktop/coal.ico + make_icon.py (Pillow; checked visually), desktop/.gitignore. Wrapper gained the named mutex (Windows only, untested).
- Sat 02:54 Task 4 done: desktop/SMOKE_TEST.md (10-step manual checklist: prepare, download + SmartScreen, install, offline first launch + single instance, each example, own files incl. drag-drop/OneDrive/Open with, viewer, export, methods sentence, close + uninstall with leftover checks; optional 5-step accessibility pass) and desktop/test_smoke.py. Ran here: 11 tests, 6 OK, 1 expected failure (UTF-8 BOM, engine), 4 skipped (RealLocalApp x3: no Gradio; WindowApi: no pywebview). What ran: all 4 examples (5 runs incl. mito code 1 and 2) through coal_app.tools with the engine started as `coal_desktop.py -m coal` (frozen-build route); the same 5 runs through a real child process started by LocalServer + serve() with a stdlib stand-in for Gradio's HTTP API (/gradio_api/call/codon_align + SSE), then clean shutdown (exit 0) and port freed; CRLF files in "OneDrive - Uni of X\Thesis data" byte-identical output. Expectations accept both today's engine and the engine branch (--keep-going, ID normalisation). With COAL_SMOKE_EXE set, the same tests target the PyInstaller COAL.exe (build.ps1 does this).
- Sat 02:55 Task 5 done: desktop/README.md (install incl. exact SmartScreen/Edge clicks and Smart App Control, use, reading results, common problems incl. the BOM workaround, uninstall, updates, privacy). No screenshots (no GUI on this machine). Final run: engine 11 OK; app 52 OK/2 skipped; desktop 34 tests, 29 OK/4 skipped/1 expected failure. Summary written. Committed locally; not pushed (per launching session). Stopping.
