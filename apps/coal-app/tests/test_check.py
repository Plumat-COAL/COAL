"""`coal-gui --check`: what the Ka/Ks and pipeline tabs can do on this machine (new-user audit 2026-11-14, U19)."""
import io
import unittest
from contextlib import redirect_stdout
from unittest import mock

from coal_app import dnds
from coal_app.__main__ import main


def check(status: dnds.ToolStatus) -> str:
    out = io.StringIO()
    with mock.patch.object(dnds, "tool_status", return_value=status), redirect_stdout(out):
        rc = main(["--check"])
    assert rc == 0, rc
    return out.getvalue()


class TestCheck(unittest.TestCase):
    def test_tools_found(self):
        text = check(dnds.ToolStatus("/x/coal-cpp", "/y/diamond", "2.2.8"))
        self.assertIn("coal-cpp: /x/coal-cpp (Ka/Ks and dN/dS pipeline tabs enabled)", text)
        self.assertIn("DIAMOND: /y/diamond (version 2.2.8)", text)

    def test_tools_missing(self):
        text = check(dnds.ToolStatus(None, None))
        self.assertIn("coal-cpp: not found (Ka/Ks and dN/dS pipeline tabs disabled;", text)
        self.assertIn("docs/install.md", text)
        self.assertIn(f"DIAMOND: not found (the dN/dS pipeline needs it: {dnds.install_command()})", text)
        self.assertIn("examples: ", text)   # the lines before stay


if __name__ == "__main__":
    unittest.main()
