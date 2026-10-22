@echo off
rem Set up COAL's MCP server for Claude (Claude Code or Claude Desktop), from a source checkout.
rem Creates .venv-mcp in the repo folder (separate from the GUI's .venv-gui: the two extras conflict),
rem installs the COAL engine and app with the app's mcp extra, checks it imports, and prints what to register with
rem Claude. This file is in apps\coal-app\desktop; it works from the repository root.
rem Run it again any time: it repairs a half-finished install and prints the registration again.
setlocal EnableExtensions
cd /d "%~dp0..\..\.." || goto :fail
set "VENV=%CD%\.venv-mcp"
set "PY=%VENV%\Scripts\python.exe"
set "CHECK=import coal_app.mcp_stdio, mcp.server.mcpserver"

if exist "%PY%" (
  "%PY%" -c "%CHECK%" >nul 2>&1 && goto :register
  echo COAL: the MCP environment in "%VENV%" is incomplete; repairing it.
)
call :find_python || goto :fail
if not exist "%PY%" (
  echo COAL: creating "%VENV%"
  %BASEPY% -m venv "%VENV%" || goto :fail
)
"%PY%" -m pip install -e . -e "apps\coal-app[mcp]" || goto :fail
"%PY%" -c "%CHECK%" || goto :fail

:register
"%PY%" -m coal_app --check || goto :fail
echo.
"%PY%" apps\coal-app\desktop\mcp_config.py || goto :fail
echo.
pause
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
