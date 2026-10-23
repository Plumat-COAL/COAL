"""site/explainer.html in a real browser: columns line up at every step, labels match captions, no horizontal
overflow at phone width. Needs Playwright + Chromium; skipped otherwise."""
import glob
import importlib.util
import os
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAGE = "file://" + os.path.join(ROOT, "site", "explainer.html")


def _chromium():
    roots = [os.environ.get("PLAYWRIGHT_BROWSERS_PATH", ""), "/opt/pw-browsers"]
    for r in roots:
        for p in glob.glob(os.path.join(r, "chromium-*", "chrome-linux", "chrome")) if r else []:
            return p
    return None


HAVE = importlib.util.find_spec("playwright") is not None and _chromium() is not None

GEOMETRY = """() => {
  const rows = [...document.querySelectorAll('.row')].map(r => [...r.querySelectorAll('.cell')]
      .map(c => Math.round(c.getBoundingClientRect().left)));
  const bad = document.querySelector('.bad');
  return {rows, scroll: document.documentElement.scrollWidth, width: window.innerWidth,
          after: getComputedStyle(bad, '::after').content, step: document.getElementById('stage').className};
}"""


@unittest.skipUnless(HAVE, "playwright with a Chromium build not available")
class TestExplainer(unittest.TestCase):
    def _steps(self, width, height):
        from playwright.sync_api import sync_playwright
        out = []
        with sync_playwright() as pw:
            b = pw.chromium.launch(executable_path=_chromium())
            try:
                pg = b.new_page(viewport={"width": width, "height": height}, reduced_motion="reduce")
                pg.goto(PAGE, timeout=30000)
                pg.evaluate("document.getElementById('play').click()")  # pause autoplay (the button is hidden), step by hand
                for i in range(5):
                    if i:
                        pg.click("#next")
                    pg.wait_for_timeout(100)  # reduced motion: no transitions
                    out.append(pg.evaluate(GEOMETRY))
            finally:
                b.close()
        return out

    def test_columns_line_up_at_every_step_and_labels_match_captions(self):
        steps = self._steps(1280, 720)
        for i, g in enumerate(steps):
            self.assertEqual(g["step"], f"stage s{i}")
            self.assertEqual(g["rows"][0], g["rows"][1], f"step {i}: Human and Mouse cells not aligned")
        self.assertIn("STOP", steps[3]["after"])
        self.assertIn("W (mito)", steps[4]["after"])

    def test_phone_width_has_no_horizontal_overflow(self):
        for i, g in enumerate(self._steps(390, 844)):
            self.assertLessEqual(g["scroll"], g["width"], f"step {i}")
            self.assertEqual(g["rows"][0], g["rows"][1])


if __name__ == "__main__":
    unittest.main()
