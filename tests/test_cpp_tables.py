"""C++ genetic code tables == the reviewed Python tables.

Reference: src/coal/tables.py (checked against NCBI: 64/64 codons and every start set; table 3 includes GTG),
from COAL_REF_TABLES, this repo's src/coal, or REF_BRANCH via `git show`. The C++
tables come from `coal-cpp-dump tables`. Every supported code must have the same 64 amino acids (NCBI
codon order TCAG) and exactly the same start codons, and the set of codes must be the same.
"""
from __future__ import annotations

import sys
import tempfile
import unittest
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_cpp_readers import REF_BRANCH, cpp_dump, load_reference, run_dump

CODONS = ["".join(c) for c in product("TCAG", repeat=3)]


class TestTables(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dump = cpp_dump()
        if cls.dump is None:
            raise unittest.SkipTest("coal-cpp-dump not built (cmake --build build)")
        cls._tmp = tempfile.TemporaryDirectory(prefix="coal_tables_")
        cls.ref = load_reference("tables", Path(cls._tmp.name))
        if cls.ref is None:
            raise unittest.SkipTest(f"reference tables not available ({REF_BRANCH} not fetched; or set COAL_REF_TABLES)")
        cls.cpp = run_dump(cls.dump, "tables")

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_same_codes(self):
        self.assertEqual(sorted(int(k) for k in self.cpp), list(self.ref.SUPPORTED))
        self.assertEqual(list(self.ref.SUPPORTED), [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23])

    def test_amino_acids_and_starts(self):
        for n in self.ref.SUPPORTED:
            with self.subTest(code=n):
                code = self.ref.get_code(n)
                cpp = self.cpp[str(n)]
                self.assertEqual(cpp["aa"], "".join(code.code[c] for c in CODONS))
                cpp_starts = {CODONS[i] for i, f in enumerate(cpp["starts"]) if f == "M"}
                self.assertEqual(cpp_starts, set(code.starts))

    def test_table3_starts_include_gtg(self):
        self.assertIn("GTG", {CODONS[i] for i, f in enumerate(self.cpp["3"]["starts"]) if f == "M"})


if __name__ == "__main__":
    unittest.main()
