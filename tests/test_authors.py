"""COAL's two authors, Chloe Newton and Keanu Murati, are named in every place that declares authorship."""
import os
import re
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUTHORS = ["Chloe Newton", "Keanu Murati"]
HANDLES = ["chloe-newt", "Keanu-MUR"]


def _read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as fh:
        return fh.read()


class TestAuthors(unittest.TestCase):
    def test_python_packages(self):
        for parts in (("pyproject.toml",), ("apps", "coal-app", "pyproject.toml")):
            line = re.search(r"^authors = \[(.*)\]$", _read(*parts), re.M)  # no tomllib on Python 3.10
            self.assertIsNotNone(line, os.path.join(*parts))
            names = re.findall(r'\{ name = "([^"]+)" \}', line.group(1))
            self.assertEqual(names, AUTHORS, os.path.join(*parts))

    def test_citation_cff(self):
        head = _read("CITATION.cff").split("references:")[0]
        names = [f"{g} {f}" for g, f in re.findall(r"- given-names: (\S+)\n\s+family-names: (\S+)", head)]
        self.assertEqual(names, AUTHORS)
        for name, handle in zip(AUTHORS, HANDLES):  # each author's GitHub profile, so the credit links to them
            given, family = name.split()
            self.assertIn(f'  - given-names: {given}\n    family-names: {family}\n    alias: {handle}\n'
                          f'    website: "https://github.com/{handle}"\n', head)

    def test_readme_authors_section(self):
        readme = _read("README.md")
        section = readme[readme.index("\n## Authors\n"):]
        for name, handle in zip(AUTHORS, HANDLES):
            self.assertIn(f"**{name}** ([@{handle}](https://github.com/{handle}))", section)

    def test_coalr_description(self):
        desc = _read("packaging", "r", "coalr", "DESCRIPTION")
        for given, family in (a.split() for a in AUTHORS):
            self.assertRegex(desc, rf'person\("{given}", "{family}"[^)]*role = c\("aut"')

    def test_copyright_lines(self):
        for parts in (("NOTICE",), ("packaging", "r", "coalr", "inst", "COPYRIGHTS")):
            self.assertIn("Copyright (C) 2026 Chloe Newton and Keanu Murati", _read(*parts), os.path.join(*parts))
        self.assertIn("Developed by Chloe Newton and Keanu Murati.", _read("site", "index.html"))
        self.assertIn('"LegalCopyright", "Copyright (C) 2026 Chloe Newton and Keanu Murati.',
                      _read("apps", "coal-app", "desktop", "coal_desktop.spec"))


if __name__ == "__main__":
    unittest.main()
