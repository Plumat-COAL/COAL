import os
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestPackaging(unittest.TestCase):
    def test_ci_workflow_exists(self):
        self.assertTrue(os.path.isfile(os.path.join(ROOT, ".github", "workflows", "ci.yml")))

    def test_app_imports_without_gradio(self):
        import coal_app  # noqa: F401  (must not need gradio)


if __name__ == "__main__":
    unittest.main()
