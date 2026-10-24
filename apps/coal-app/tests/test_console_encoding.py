"""The GUI's console output must not die on a non-UTF-8 console (stress-test finding: `python -m coal_app` on
Windows with cp1252 stdout crashed on Gradio's MCP banner emoji)."""
import io
import sys
import unittest


class TolerantConsole(unittest.TestCase):
    def test_cp1252_stdout_takes_an_emoji_after_tolerant_console(self):
        from coal_app.__main__ import tolerant_console
        old_out, old_err = sys.stdout, sys.stderr
        try:
            sys.stdout = io.TextIOWrapper(io.BytesIO(), encoding="cp1252", errors="strict")
            sys.stderr = io.TextIOWrapper(io.BytesIO(), encoding="cp1252", errors="strict")
            with self.assertRaises(UnicodeEncodeError):
                print("\U0001f528 Launching MCP server", file=sys.stdout)
            tolerant_console()
            print("\U0001f528 Launching MCP server", file=sys.stdout)   # no exception
            sys.stdout.flush()
            self.assertEqual(sys.stdout.errors, "replace")
        finally:
            sys.stdout, sys.stderr = old_out, old_err

    def test_utf8_console_is_left_alone(self):
        from coal_app.__main__ import tolerant_console
        old = sys.stdout
        try:
            sys.stdout = io.TextIOWrapper(io.BytesIO(), encoding="utf-8", errors="strict")
            tolerant_console()
            self.assertEqual(sys.stdout.errors, "strict")
        finally:
            sys.stdout = old


if __name__ == "__main__":
    unittest.main()
