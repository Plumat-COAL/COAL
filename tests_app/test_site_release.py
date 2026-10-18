"""WEB-18: the deployed site matches the latest release (scripts/site_for_release.py, run by pages.yml)."""
import importlib.util
import os
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _script():
    spec = importlib.util.spec_from_file_location("site_for_release", os.path.join(ROOT, "scripts", "site_for_release.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _read(name):
    with open(os.path.join(ROOT, "site", name), encoding="utf-8") as fh:
        return fh.read()


class TestSiteForRelease(unittest.TestCase):
    def setUp(self):
        self.m = _script()

    def test_older_release_hides_newer_features_and_shows_its_version(self):
        out = self.m.rewrite(_read("install.html"), "0.4.1")
        self.assertNotIn("Connect COAL to Claude", out)
        self.assertNotIn("data-since", out)
        self.assertIn("<li>Start <strong>COAL</strong>", out)  # the next step is untouched
        idx = self.m.rewrite(_read("index.html"), "0.4.1")
        self.assertIn('<span id="version" data-version="0.4.1">0.4.1</span>', idx)
        self.assertIn("<span data-version-ref>0.4.1</span>", self.m.rewrite(_read("cite.html"), "0.4.1"))

    def test_same_or_newer_release_keeps_features(self):
        for v in ("0.4.2", "0.5.0", "1.0"):
            out = self.m.rewrite(_read("install.html"), v)
            self.assertIn("Connect COAL to Claude", out)
            self.assertNotIn("data-since", out)

    def test_idempotent(self):
        once = self.m.rewrite(_read("install.html"), "0.4.1")
        self.assertEqual(self.m.rewrite(once, "0.4.1"), once)

    def test_bad_arguments(self):
        self.assertEqual(self.m.main(["latest", "x.html"]), 2)
        self.assertEqual(self.m.main(["0.4.1"]), 2)

    def test_pages_workflow_runs_it_on_the_deployed_copy(self):
        with open(os.path.join(ROOT, ".github", "workflows", "pages.yml"), encoding="utf-8") as fh:
            wf = fh.read()
        self.assertIn('python scripts/site_for_release.py "${tag#v}" _site/*.html', wf)
        self.assertLess(wf.index("site_for_release"), wf.index("upload-pages-artifact"))


if __name__ == "__main__":
    unittest.main()
