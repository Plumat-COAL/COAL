# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for the COAL desktop app: ONE-FOLDER build (not one-file, see DESIGN.md 5.7).
#
# Build on Windows with desktop\build.ps1, or by hand from the repo root:
#   desktop\.venv\Scripts\pyinstaller --noconfirm --clean desktop\coal_desktop.spec ^
#       --distpath desktop\dist --workpath desktop\build
# Result: desktop\dist\COAL\COAL.exe plus its _internal\ folder.
#
# Written for PyInstaller 6.x. NOT BUILT YET: written on a Linux machine without PyPI access,
# so the Gradio/pywebview collection below is untested (see overnight/windows-app/LOG.md).
#
# What goes in:
#   - Python + COAL engine (src/coal) + local app (src/coal_app) + the wrapper
#   - Gradio and its runtime pieces (collected whole: it loads templates, JS and .pyi at run time)
#   - pywebview (+ pythonnet / WebView2 interop DLLs, via pywebview's hooks)
#   - examples/ (the four "Try an example" datasets) -> <bundle>\_internal\examples
#   - coal-cpp.exe if it has been built (cpp\build\**\coal-cpp.exe) -> <bundle>\_internal\bin
#   - LICENSE
import glob
import os
import re

from PyInstaller.utils.hooks import collect_all, collect_submodules

REPO = os.path.abspath(os.path.join(SPECPATH, ".."))  # noqa: F821 - SPECPATH is set by PyInstaller
SRC = os.path.join(REPO, "src")

with open(os.path.join(SRC, "coal", "__init__.py"), encoding="utf-8") as fh:
    VERSION = re.search(r'__version__\s*=\s*"([^"]+)"', fh.read()).group(1)

datas, binaries, hiddenimports = [], [], []

# Gradio and friends: data files (frontend build, templates), metadata (importlib.metadata
# version lookups), and every submodule (uvicorn and fastapi import some by name).
for pkg in ("gradio", "gradio_client", "safehttpx", "groovy", "uvicorn", "fastapi", "starlette",
            "webview"):
    try:
        d, b, h = collect_all(pkg)
    except Exception as e:  # a package missing from the venv must fail loudly, not silently
        raise SystemExit(f"coal_desktop.spec: could not collect {pkg}: {e}")
    datas += d
    binaries += b
    hiddenimports += h

hiddenimports += collect_submodules("coal") + collect_submodules("coal_app")

datas.append((os.path.join(SRC, "coal_app", "sample_report.json"), "coal_app"))
datas.append((os.path.join(REPO, "LICENSE"), "."))

# Example datasets, keeping their folder names (core.example_paths expects <dir>/<name>/...).
for f in glob.glob(os.path.join(REPO, "examples", "*", "*")):
    if os.path.isfile(f) and f.endswith((".fasta", ".md")):
        datas.append((f, os.path.join("examples", os.path.basename(os.path.dirname(f)))))

# Optional C++ engine.
cpp = sorted(glob.glob(os.path.join(REPO, "cpp", "build", "**", "coal-cpp.exe"), recursive=True))
if cpp:
    binaries.append((cpp[0], "bin"))
    print(f"coal_desktop.spec: bundling C++ engine {cpp[0]}")
else:
    print("coal_desktop.spec: no cpp\\build\\**\\coal-cpp.exe found; building without the C++ engine")

# Windows file properties (Explorer > Properties > Details). Version info also makes the .exe
# look less like anonymous malware to antivirus heuristics.
version_info = None
try:
    from PyInstaller.utils.win32.versioninfo import (FixedFileInfo, StringFileInfo, StringStruct,
                                                     StringTable, VarFileInfo, VarStruct,
                                                     VSVersionInfo)
    nums = tuple(int(x) for x in (VERSION.split(".") + ["0", "0", "0"])[:4] if x.isdigit())
    nums = (nums + (0, 0, 0, 0))[:4]
    version_info = VSVersionInfo(
        ffi=FixedFileInfo(filevers=nums, prodvers=nums),
        kids=[
            StringFileInfo([StringTable("040904B0", [
                StringStruct("CompanyName", "COAL contributors"),
                StringStruct("FileDescription", "COAL codon alignment"),
                StringStruct("FileVersion", VERSION),
                StringStruct("InternalName", "COAL"),
                StringStruct("LegalCopyright", "MIT licence"),
                StringStruct("OriginalFilename", "COAL.exe"),
                StringStruct("ProductName", "COAL"),
                StringStruct("ProductVersion", VERSION),
            ])]),
            VarFileInfo([VarStruct("Translation", [1033, 1200])]),
        ],
    )
except ImportError:
    pass

a = Analysis(  # noqa: F821
    [os.path.join(SPECPATH, "coal_desktop.py")],  # noqa: F821
    pathex=[SRC],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=["tkinter", "matplotlib", "IPython", "pytest"],
    # Gradio reads its own source files at run time (e.g. inspect, .pyi event stubs); keep it as
    # plain .py files on disk rather than only bytecode inside the archive.
    module_collection_mode={"gradio": "py", "gradio_client": "py"},
    noarchive=False,
)
pyz = PYZ(a.pure)  # noqa: F821

exe = EXE(  # noqa: F821
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="COAL",
    icon=os.path.join(SPECPATH, "coal.ico"),  # noqa: F821
    version=version_info,
    console=False,          # windowed: no console window
    # CHANGED(WIN-R1): was False (modal traceback box in hidden processes); now errors go to the logs.
    # See overnight/windows-app/CHANGES.md
    disable_windowed_traceback=True,
    upx=False,              # UPX is a common antivirus false-positive trigger
    debug=False,
    strip=False,
)

coll = COLLECT(  # noqa: F821
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="COAL",
)
