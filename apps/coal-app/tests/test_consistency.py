"""Nothing drifts: the app's genetic codes equal the engine's tables, UTF-8 is named everywhere, the demo wording matches what the app does."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # helpers.py

from helpers import APP, ROOT  # noqa: E402


class TestNoDrift(unittest.TestCase):
    def test_app_genetic_codes_equal_the_engine_tables(self):
        from coal import tables
        from coal_app import gcode
        self.assertEqual(gcode._RAW, tables._RAW)

    def test_viewer_defines_everything_before_main_block(self):
        import ast
        with open(os.path.join(APP, "src", "coal_app", "viewer.py"), encoding="utf-8") as fh:
            body = ast.parse(fh.read()).body
        main_at = next(i for i, n in enumerate(body) if isinstance(n, ast.If))
        self.assertEqual([n.name for n in body[main_at:] if isinstance(n, ast.FunctionDef)], [])
        self.assertEqual(main_at, len(body) - 1)

    def test_test_classes_defined_before_unittest_main(self):
        import ast
        import glob
        for path in glob.glob(os.path.join(APP, "tests", "test_*.py")):
            with open(path, encoding="utf-8") as fh:
                body = ast.parse(fh.read()).body
            ifs = [i for i, n in enumerate(body) if isinstance(n, ast.If)]
            if ifs:
                self.assertEqual(ifs[-1], len(body) - 1, path)  # finding 16: nothing after the main block


class TestUtf8Everywhere(unittest.TestCase):
    FILES = ["demo/diamond_pipeline.py", "demo/make_mcp_demo.py", "apps/coal-app/tests/diamond_stub.py",
             "apps/coal-app/src/coal_app/core.py", "apps/coal-app/src/coal_app/viewer.py", "apps/coal-app/src/coal_app/tools.py", "apps/coal-app/src/coal_app/ui.py",
             "src/coal_compat/cli.py", "apps/coal-app/space/build_space.py", "examples/make_examples.py"]

    def test_text_files_and_pipes_name_utf8(self):
        import ast
        bad = []
        for rel in self.FILES:
            with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
                tree = ast.parse(fh.read())
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                fn = node.func
                name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
                kw = {k.arg for k in node.keywords}
                if name == "open":
                    mode = node.args[1].value if len(node.args) > 1 and isinstance(node.args[1], ast.Constant) \
                        else next((k.value.value for k in node.keywords if k.arg == "mode"), "r")
                    if "b" not in str(mode) and "encoding" not in kw:
                        bad.append(f"{rel}:{node.lineno} open")
                elif name == "run" and "text" in kw and "encoding" not in kw:
                    bad.append(f"{rel}:{node.lineno} subprocess.run(text=True)")
        self.assertEqual(bad, [])


class TestDemoWording(unittest.TestCase):
    def test_location_is_attributed_to_the_app_estimate(self):
        for rel in ("demo/mcp_demo.md", "demo/make_mcp_demo.py"):
            with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
                text = fh.read()
            self.assertNotIn("says *where*: a missing base", text, rel)
            self.assertNotIn("what is wrong and where: a missing base", text, rel)
        with open(os.path.join(ROOT, "demo", "mcp_demo.md"), encoding="utf-8") as fh:
            self.assertIn("COAL's app estimates", fh.read())


if __name__ == "__main__":
    unittest.main()
