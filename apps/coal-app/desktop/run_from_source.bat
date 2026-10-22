@echo off
rem COAL desktop window, run from a source checkout (git clone). Double-click it, or run it in cmd.
rem First start: creates .venv-gui in the repo folder and installs the COAL engine and app, with the app's GUI
rem extra (a few minutes). This file is in apps\coal-app\desktop; it works from the repository root.
rem Later starts: checks the environment by importing COAL, Gradio and pywebview, repairs it if that fails,
rem then opens the window. Arguments are passed on (e.g. --browser, --mcp).
rem The GUI and the MCP server need separate environments (gradio[mcp] and the mcp extra conflict):
rem for Claude, use setup_mcp.bat, which uses .venv-mcp.
setlocal EnableExtensions
cd /d "%~dp0..\..\.." || goto :fail
set "VENV=%CD%\.venv-gui"
set "PY=%VENV%\Scripts\python.exe"
set "CHECK=import coal_app.ui, gradio, webview"

if exist "%PY%" (
  "%PY%" -c "%CHECK%" >nul 2>&1 && goto :run
  echo COAL: the GUI environment in "%VENV%" is incomplete; repairing it.
)
call :find_python || goto :fail
if not exist "%PY%" (
  echo COAL: creating "%VENV%"
  %BASEPY% -m venv "%VENV%" || goto :fail
)
rem pywebview: same pin as apps\coal-app\desktop\requirements.txt (apps/coal-app/tests/test_launchers.py keeps them equal)
"%PY%" -m pip install -e . -e "apps\coal-app[gui]" "pywebview==6.2.1" || goto :fail
"%PY%" -c "%CHECK%" || goto :fail

:run
"%PY%" apps\coal-app\desktop\coal_desktop.py %* || goto :fail
exit /b 0

:find_python
rem Prefer the Python launcher with 3.12 (the version COAL's Windows build uses), then any Python 3,
rem then "python" on PATH. "python" can be the Microsoft Store stub, which fails the test below.
set "BASEPY="
py -3.12 -c "import sys" >nul 2>&1 && set "BASEPY=py -3.12"
if not defined BASEPY py -3 -c "import sys" >nul 2>&1 && set "BASEPY=py -3"
if not defined BASEPY python -c "import sys" >nul 2>&1 && set "BASEPY=python"
if not defined BASEPY (
  echo COAL: Python was not found. Install Python 3.12 from https://www.python.org/downloads/windows/
  echo       ^(tick "Add python.exe to PATH" or keep the "py launcher" option^), then run this file again.
  exit /b 1
)
%BASEPY% -c "import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)" >nul 2>&1
if errorlevel 1 (
  echo COAL: needs Python 3.10 or newer; "%BASEPY%" is older.
  echo       Install Python 3.12 from https://www.python.org/downloads/windows/ and run this file again.
  exit /b 1
)
echo COAL: using %BASEPY%
exit /b 0

:fail
echo.
echo COAL: setup did not finish (see the messages above). Run this file again to repair it.
pause
exit /b 1
