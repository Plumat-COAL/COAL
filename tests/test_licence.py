"""The licence is declared the same everywhere: GPL-3.0-only from the release after 0.5.0 (0.5.0 and earlier: MIT)."""
import hashlib
import os
import re
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as fh:
        return fh.read()


class TestLicence(unittest.TestCase):
    def test_licence_file_is_the_gpl3_text(self):
        text = _read("LICENSE")
        self.assertIn("GNU GENERAL PUBLIC LICENSE", text.splitlines()[0])
        self.assertIn("Version 3, 29 June 2007", text.splitlines()[1])
        # the FSF's text, unchanged (the same file KaKs_Calculator 2.0 and PAML ship); line endings aside
        digest = hashlib.md5(text.replace("\r\n", "\n").encode("utf-8")).hexdigest()
        self.assertEqual(digest, "1ebbd3e34237af26da5dc08a4e440464")
        for copy in ("LICENSE", "NOTICE"):  # the app is built from its own folder and carries copies
            self.assertEqual(_read(copy), _read("apps", "coal-app", copy), f"apps/coal-app/{copy} is not a copy")

    def test_notice_has_the_copyright_and_the_earlier_licence(self):
        notice = _read("NOTICE")
        for must in ("Copyright (C) 2026 Chloe Newton and Keanu Murati", "version 3 of the License only",
                     "GPL-3.0-only", "0.5.0 and earlier releases were published under the MIT License",
                     "OFL.txt", "LICENSE-hyphy.txt", "LICENSE-biopython.rst"):
            self.assertIn(must, notice)

    def test_every_declaration_says_gpl3_only(self):
        for rel in (("pyproject.toml",), ("apps", "coal-app", "pyproject.toml")):
            text = _read(*rel)
            self.assertRegex(text, r'(?m)^license = \{text = "GPL-3\.0-only"\}', rel)
            self.assertIn("License :: OSI Approved :: GNU General Public License v3 (GPLv3)", text, rel)
        self.assertRegex(_read("CITATION.cff"), r"(?m)^license: GPL-3\.0-only$")
        meta = _read("packaging", "bioconda", "meta.yaml")
        self.assertRegex(meta, r"(?m)^  license: GPL-3\.0-only$")
        self.assertRegex(meta, r"(?m)^  license_family: GPL3$")
        self.assertRegex(_read("apps", "coal-app", "space", "README.md"), r"(?m)^license: gpl-3\.0$")
        self.assertIn("GPL-3.0-only", _read("README.md"))
        self.assertIn("0.5.0 and earlier were\nreleased under the MIT License and stay MIT", _read("README.md"))
        self.assertIn("install(FILES LICENSE NOTICE", _read("CMakeLists.txt"))

    def test_kaks_engine_attribution(self):
        """The Ka/Ks engine is derived from KaKs_Calculator 2.0 and PAML yn00 (GPL v3), and says so everywhere."""
        notices = ("Copyright (C) 2009, Beijing Institute of Genomics of CAS", "Da-Peng Wang", "Yu-Bin Zhang",
                   "Zhang Zhang", "Kristian K Ullrich", "Copyright, 1998, Ziheng Yang", "yn00.c",
                   "GNU General Public License, version 3")
        for rel in (("cpp", "src", "kaks.cpp"), ("NOTICE",), ("packaging", "r", "coalr", "inst", "COPYRIGHTS"),
                    ("packaging", "r", "coalr", "src", "coal", "kaks.cpp")):
            text = _read(*rel)
            for must in notices:
                self.assertIn(must, text, rel)
        head = _read("cpp", "src", "kaks.cpp").split("#include")[0]
        self.assertIn("SPDX-License-Identifier: GPL-3.0-only", head)
        self.assertIn("is derived from", head)
        self.assertIn("Modified for COAL by Chloe Newton and Keanu Murati, 2026. What was changed:", head)
        for rel in (("cpp", "src", "kaks.hpp"), ("cpp", "src", "kaks_cli.cpp")):
            flat = " ".join(_read(*rel).replace("//", " ").split())
            self.assertTrue("Beijing Institute of Genomics of CAS" in flat, rel)
        desc = _read("packaging", "r", "coalr", "DESCRIPTION")
        for must in ('person("Ziheng", "Yang", role = c("ctb", "cph")', "Beijing Institute of Genomics of CAS",
                     "Copyright: See the file inst/COPYRIGHTS"):
            self.assertIn(must, desc)
        # the claim the independent review of 2026-11-12 found to be false is gone, everywhere
        for rel in (("docs", "kaks.md"), ("README.md",), ("CHANGELOG.md",), ("docs", "r.md"), ("docs", "pipeline.md"),
                    ("packaging", "r", "coalr", "README.md"), ("cpp", "src", "kaks.cpp"), ("cpp", "src", "kaks.hpp"),
                    ("cpp", "src", "kaks_cli.cpp")):
            text = " ".join(_read(*rel).split())
            for false_claim in ("independent implementation of the published methods",
                                "none of their code is copied", "Independent implementation",
                                "a drop-in for `KaKs_Calculator", "is a drop-in for KaKs_Calculator",
                                "KaKs_Calculator drop-in"):
                self.assertNotIn(false_claim, text, rel)
        kaks_doc = " ".join(_read("docs", "kaks.md").split())
        self.assertIn("The engine is not an independent implementation", kaks_doc)
        self.assertIn("derived from the source code of KaKs_Calculator 2.0", kaks_doc)

    def test_nothing_still_calls_coal_mit(self):
        """'MIT' is left only where it is true: the earlier releases, and third-party test data."""
        allowed = re.compile(r"0\.5\.0|earlier|MIT at the time|HyPhy|hyphy|license_family: MIT|license: MIT")
        skip_dirs = {".git", "build", "node_modules", "__pycache__", "realdata", "PAL2NAL"}
        left = []
        for d, dirs, files in os.walk(ROOT):
            dirs[:] = [x for x in dirs if x not in skip_dirs and not x.endswith(".egg-info")
                       and not x.startswith("build")]
            for f in files:
                if not f.endswith((".md", ".toml", ".cff", ".yaml", ".yml", ".html", ".py", ".spec", ".iss", ".txt",
                                   ".cpp", ".hpp", ".in", ".R", ".Rd")) and f not in ("LICENSE", "NOTICE",
                                                                                       "DESCRIPTION"):
                    continue
                path = os.path.join(d, f)
                if os.path.abspath(path) == os.path.abspath(__file__):
                    continue
                try:
                    with open(path, encoding="utf-8") as fh:
                        lines = fh.read().splitlines()
                except (UnicodeDecodeError, OSError):
                    continue
                for i, line in enumerate(lines):
                    if re.search(r"\bMIT\b", line):
                        context = " ".join(lines[max(0, i - 2):i + 3])
                        if not allowed.search(context):
                            left.append(f"{os.path.relpath(path, ROOT)}:{i + 1}: {line.strip()[:100]}")
        self.assertEqual(left, [])


if __name__ == "__main__":
    unittest.main()
