# Building COAL for Windows

For teammates making the Windows app. (People who just want to use COAL: see `README.md`.)

Status: these files were written on Linux without PyPI, PyInstaller or Inno Setup, so **the
build has never been run**. Expect to fix small things on the first run; note them in
`overnight/windows-app/LOG.md`.

## Once per machine

1. Windows 10 or 11, 64-bit.
2. Python 3.12 from python.org (tick "Add python.exe to PATH" is not needed; the `py` launcher
   is enough). Check: `py -3.12 --version`.
3. Git, and this repo checked out on branch `keanu_overnight/windows-app`.
4. Optional, for the installer: Inno Setup 6.3 or newer from https://jrsoftware.org/isdl.php
   (choose "Install for me only" if you have no admin rights).

## Build (one command)

In PowerShell, from the repo root:

```powershell
powershell -ExecutionPolicy Bypass -File desktop\build.ps1
```

It will, each step with a timeout:

1. create `desktop\.venv` and install `.[gui]` and `desktop\requirements.txt` (pinned);
2. run the engine, app, wrapper and smoke tests (`-SkipTests` to skip);
3. run PyInstaller with `desktop\coal_desktop.spec` into `desktop\dist\COAL\`;
4. check the frozen app: `COAL.exe -m coal --version` and the smoke test against `COAL.exe`;
5. if Inno Setup is installed, build `desktop\dist\COAL-Setup-<version>.exe`.

Options: `-Clean` (start from scratch), `-NoInstaller`, `-Python -3.13`.

## The same, step by step

```powershell
py -3.12 -m venv desktop\.venv
desktop\.venv\Scripts\python -m pip install -e ".[gui]" -r desktop\requirements.txt
desktop\.venv\Scripts\python -m unittest desktop\test_wrapper.py
desktop\.venv\Scripts\python -m unittest desktop\test_smoke.py -v
desktop\.venv\Scripts\python desktop\coal_desktop.py          # try the window from source
desktop\.venv\Scripts\python -m PyInstaller --noconfirm --clean desktop\coal_desktop.spec --distpath desktop\dist --workpath desktop\build
desktop\dist\COAL\COAL.exe                                       # try the built app
$env:COAL_SMOKE_EXE = "$PWD\desktop\dist\COAL\COAL.exe"; desktop\.venv\Scripts\python -m unittest desktop\test_smoke.py -v; Remove-Item Env:\COAL_SMOKE_EXE
& "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe" /DAppVersion=0.3.0 desktop\installer.iss
```

(For a machine-wide Inno Setup the last line is
`& "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe" /DAppVersion=0.3.0 desktop\installer.iss`.)

## C++ engine

If `cpp\build\...\coal-cpp.exe` exists when PyInstaller runs, the spec copies it into the
bundle (`_internal\bin\coal-cpp.exe`) and says so in the build output. Nothing calls it yet:
the GUI always runs `python -m coal`. Using it needs a GUI/engine decision (LOG.md).

## Sharing a test build

- Installer: send `desktop\dist\COAL-Setup-<version>.exe`.
- Portable: zip `desktop\dist\COAL\` and send the zip; unzip anywhere, run `COAL.exe`.
- Both are unsigned: testers will see SmartScreen (see `README.md` for the clicks).

## Likely first-build problems (and where to look)

| Symptom | Probable cause | Fix to try |
|---|---|---|
| `COAL.exe` opens nothing; `%LOCALAPPDATA%\COAL\logs\local_app.log` shows a missing module | a module Gradio/uvicorn imports by name | add it to `hiddenimports` in the spec |
| log shows a missing file inside `gradio\` | a Gradio data file not collected | `collect_all` should cover it; check `module_collection_mode` |
| window says "Starting COAL…" forever then an error after 60 s | local app crashed: see `local_app.log` | as above |
| "Try an example" buttons missing | examples not found in the bundle | check `_internal\examples\clean\dna.fasta` exists |
| a console window flashes on each run | engine started without `CREATE_NO_WINDOW` | should not happen in the windowed build; report it |
| Defender quarantines `COAL.exe` | heuristic false positive | DESIGN.md 5.7; submit the file to Microsoft |
