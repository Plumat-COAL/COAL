# Overnight C (Windows desktop app): log

(Summary goes here at the end.)

---

## Progress log

- Sat 02:38 start. Machine: Linux cloud sandbox (not Windows), Python 3.13, on branch keanu_overnight/windows-app at a2351d0. Brief adapted by the launching session: skip the Windows git step; PyPI unreachable so pywebview/PyInstaller/Inno Setup can't be installed or run; commit locally only, never push.
- Sat 02:38 Note: the system Python has an editable install of coal/coal_app pointing at /home/claude/coal-app (another session's folder). Every command here sets PYTHONPATH=/home/claude/coal-win/src so this worktree's code is used; coal_desktop.py also puts ../src first on sys.path and PYTHONPATH in a source checkout.
- Sat 02:38 Problem "install pywebview": attempt 1 `pip download pywebview` (timeout 40) -> DNS failure for pypi.org. Not retried (instruction: don't try more than once). PyInstaller and Inno Setup not attempted for the same reason. Gradio is not installed either (known from overnight/app/LOG.md).
- Sat 02:38 pywebview docs were reachable instead (pywebview.flowrl.com guide/api and guide/dom, github examples/drag_drop.py, github releases page): latest release 6.2.1; DOM `drop` event gives full paths in `e['dataTransfer']['files'][i]['pywebviewFullPath']` via `window.dom.document.events.drop += DOMEventHandler(fn, prevent_default, stop_propagation)`; window events closing/closed/loaded/shown; methods restore/show/evaluate_js/destroy; `on_top` is a property; webview.start(func, args, gui, debug, private_mode, storage_path, icon, ...). No platform limitation for drop is documented. The installed-version check the brief asks for is impossible here; the pin is 6.2.1 from the releases page.
- Sat 02:40 Finding (tested here): engine reads CRLF files and paths with spaces correctly (output byte-identical to LF input), but a UTF-8 BOM makes it stop with "'ascii' codec can't decode byte 0xef in position 0". See Needs a human decision.
- Sat 02:41 Task 1 done: desktop/DESIGN.md (3 personas with moment-of-pain tables, 6 flows, 5 wireframed screens, 15 Windows decisions incl. wrapper comparison, process model, SmartScreen wording, Smart App Control, AV, OneDrive/long paths, CRLF/BOM, DPI, dark/contrast, keyboard/Narrator; large data; scope; 8 proposed changes outside desktop/).
