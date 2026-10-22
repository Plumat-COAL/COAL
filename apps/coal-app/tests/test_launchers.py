"""IMP-M-2: the Windows launchers in desktop/ (static checks; cmd.exe itself only runs on Windows)."""
import importlib.util
import json
import os
import re
import unittest

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # apps/coal-app
DESKTOP = os.path.join(APP, "desktop")
# both install the engine (the repo root) and the app (apps\coal-app) with one extra, from the repo root
LAUNCHERS = {"run_from_source.bat": (".venv-gui", '. -e "apps\\coal-app[gui]"', "webview"),
             "setup_mcp.bat": (".venv-mcp", '. -e "apps\\coal-app[mcp]"', "mcp.server.mcpserver")}


def _read(name):
    with open(os.path.join(DESKTOP, name), "rb") as fh:
        return fh.read().decode("ascii")


def _load_mcp_config():
    spec = importlib.util.spec_from_file_location("coal_mcp_config", os.path.join(DESKTOP, "mcp_config.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestLaunchers(unittest.TestCase):
    def test_crlf_line_endings(self):
        for name in LAUNCHERS:
            text = _read(name)
            self.assertEqual(text.count("\n"), text.count("\r\n"), f"{name} must use CRLF (cmd mis-reads labels)")

    def test_separate_environments_and_extras(self):
        for name, (venv, extra, module) in LAUNCHERS.items():
            text = _read(name)
            self.assertIn(f'set "VENV=%CD%\\{venv}"', text, name)
            self.assertIn(f"pip install -e {extra}", text, name)
            self.assertIn(module, re.search(r'set "CHECK=([^"]*)"', text).group(1), name)
            self.assertNotRegex(text, r"\\\.venv\\", name)  # never the shared .venv

    def test_python_choice_and_version_check(self):
        for name in LAUNCHERS:
            text = _read(name)
            pos = [text.index(s) for s in ('py -3.12 -c', 'py -3 -c', 'python -c "import sys" ')]
            self.assertEqual(pos, sorted(pos), name)
            self.assertIn("sys.version_info >= (3, 10)", text, name)
            self.assertIn("https://www.python.org/downloads/windows/", text, name)

    def test_repair_tests_imports_not_folders(self):
        for name in LAUNCHERS:
            text = _read(name)
            self.assertIn('"%PY%" -c "%CHECK%" >nul 2>&1 && goto :', text, name)
            self.assertIn('"%PY%" -c "%CHECK%" || goto :fail', text, name)
            self.assertNotIn('if exist "%VENV%"', text, name)

    def test_every_step_fails_loudly(self):
        for name in LAUNCHERS:
            for line in _read(name).splitlines():
                s = line.strip()
                if s.lower().startswith(("rem", "echo", "set ", "if ", "exit", ":", "pause", "@echo", "setlocal",
                                         "py -3", "python -c")):
                    continue  # comments, messages, probes (their failure is handled on the next line)
                if s.endswith(">nul 2>&1") or ">nul 2>&1 && goto :" in s:
                    continue  # silent probes: a failure means "repair", handled by the following lines
                if s and not s.startswith(")"):
                    self.assertTrue(s.endswith("|| goto :fail"), f"{name}: {s}")

    def test_pywebview_pin_matches_requirements(self):
        with open(os.path.join(DESKTOP, "requirements.txt"), encoding="utf-8") as fh:
            pin = re.search(r"^pywebview==\S+", fh.read(), re.M).group(0)
        self.assertIn(f'"{pin}"', _read("run_from_source.bat"))

    def test_setup_mcp_prints_the_registration(self):
        text = _read("setup_mcp.bat")
        self.assertIn('"%PY%" apps\\coal-app\\desktop\\mcp_config.py || goto :fail', text)
        self.assertIn('cd /d "%~dp0..\\..\\.." || goto :fail', text)  # the repository root


class TestMcpConfig(unittest.TestCase):
    def test_windows_path_is_json_escaped(self):
        mc = _load_mcp_config()
        py = r"C:\Users\A B\COAL\.venv-mcp\Scripts\python.exe"
        text = mc.claude_desktop_config(py)
        self.assertIn(r'"C:\\Users\\A B\\COAL\\.venv-mcp\\Scripts\\python.exe"', text)
        self.assertEqual(json.loads(text), {"mcpServers": {"coal": {"command": py,
                                                                    "args": ["-m", "coal_app", "--mcp"]}}})

    def test_claude_code_command(self):
        mc = _load_mcp_config()
        self.assertEqual(mc.claude_code_command(r"C:\A B\python.exe"),
                         'claude mcp add coal -- "C:\\A B\\python.exe" -m coal_app --mcp')
        self.assertEqual(mc.claude_code_command("/x/python"), "claude mcp add coal -- /x/python -m coal_app --mcp")


if __name__ == "__main__":
    unittest.main()
