"""C++ input readers == the reviewed Python readers (CPP-2).

Edge-case files are built in a temp dir and read by both:
  * Python: src/coal/seqio.py from the reviewed engine on origin/keanu_overnight/integration (read with
    `git show`, not merged), or the file named by COAL_REF_SEQIO;
  * C++: build/coal-cpp-dump (COAL_CPP_DUMP), the same reader code coal-cpp uses.
Both must give the same IDs, the same sequences/rows, or the same error message (with file:line).
Protein files are read as the pipeline reads them (CLUSTAL if the first line starts with CLUSTAL,
else FASTA one record at a time); DNA files through DnaIndex (IDs, then each cleaned sequence).
Skipped if the C++ helper is not built or the reference reader is not available.
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF_BRANCH = "origin/keanu_overnight/integration"
TIMEOUT = 60


def cpp_dump() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_DUMP", ROOT / "build" / "coal-cpp-dump"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file():
            return cand
    return None


def load_reference(module: str, tmp: Path):
    """The reviewed Python module src/coal/<module>.py, as a standalone module (or None)."""
    env = os.environ.get(f"COAL_REF_{module.upper()}")
    if env:
        src = Path(env).read_text(encoding="utf-8")
    else:
        try:
            r = subprocess.run(["git", "show", f"{REF_BRANCH}:src/coal/{module}.py"], cwd=ROOT,
                               capture_output=True, text=True, timeout=TIMEOUT)
        except (OSError, subprocess.SubprocessError):
            return None
        if r.returncode != 0:
            return None
        src = r.stdout
    path = tmp / f"ref_{module}.py"
    path.write_text(src, encoding="utf-8")
    spec = importlib.util.spec_from_file_location(f"ref_{module}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_dump(dump: Path, *args, cwd=None) -> dict:
    r = subprocess.run([str(dump), *map(str, args)], cwd=cwd, capture_output=True, timeout=TIMEOUT)
    out = r.stdout.decode("utf-8", "replace").strip()
    if not out:
        raise AssertionError(f"coal-cpp-dump printed nothing (exit {r.returncode}): {r.stderr!r}")
    return json.loads(out)


# name -> bytes. Protein alignments.
PROTEIN_CASES = {
    "plain": b">s1 desc\nMKV-L\n>s2\nMKVAL\n",
    "multiline_spaces": b">s1\nMK V\n-L\n\n>s2\tdesc\nMK\nVAL\n",
    "bom": b"\xef\xbb\xbf>s1\nMKV-L\n>s2\nMKVAL\n",
    "crlf": b">s1\r\nMKV-L\r\n>s2\r\nMKVAL\r\n",
    "cr_only": b">s1\rMKV-L\r>s2\rMKVAL\r",
    "no_final_newline": b">s1\nMKV-L\n>s2\nMKVAL",
    "lowercase_dots": b">s1\nmkv.l\n>s2\nmkval\n",
    "latin1_line3": b">s1\nMKVAL\nMK\xe9AL\n>s2\nMKVAL\n",
    "utf8_line2": ">s1\nMKéAL\n".encode("utf-8"),
    "nonascii_header": ">sé1\nMKVAL\n".encode("utf-8"),
    "duplicate_id": b">s1\nMKVAL\n>s2\nMKVAL\n>s1\nMKVAL\n",
    "empty_header": b">s1\nMKVAL\n>\nMKVAL\n",
    "empty_header_spaces": b">s1\nMKVAL\n>   \nMKVAL\n",
    "data_before_header": b"MKVAL\n>s1\nMKVAL\n",
    "blank_lines_first": b"\n\n  \n>s1\nMKVAL\n",
    "empty_file": b"",
    "only_header": b">s1\n",
    "clustal": b"CLUSTAL W (1.83) multiple sequence alignment\n\ns1   MKV-L 5\ns2   MKVAL 5\n     **.*:\n\n"
               b"s1   AAAA\ns2   AAAA\n",
    "clustal_bom_crlf": b"\xef\xbb\xbfCLUSTAL W\r\n\r\ns1  MKV\r\ns2  MKV\r\n\r\ns1  -L\r\ns2  AL\r\n",
    "clustal_malformed": b"CLUSTAL\n\ns1   MKVAL\nlonely\n",
    "clustal_tab_conservation": b"CLUSTAL\n\ns1 MKV\ns2 MKV\n\t***\n",
    "clustal_nonascii": "CLUSTAL\n\ns1 MKV\nsé2 MKV\n".encode("utf-8"),
}

# name -> list of file contents (several files = multi-file DNA input)
DNA_CASES = {
    "plain": [b">s1 desc\nATGAAA\n>s2\nATGCCC\n"],
    "bom": [b"\xef\xbb\xbf>s1\nATGAAA\n"],
    "crlf": [b">s1\r\nATG\r\nAAA\r\n>s2\r\nCCC\r\n"],
    "cr_only": [b">s1\rATG\rAAA\r>s2\rCCC\r"],
    "cr_only_cr_at_end": [b">s1\rATGAAA\r"],
    "mixed_crlf_and_lone_cr": [b">s1\r\nATG\rAAA\r\n>s2\rCCC\n"],
    "hash_before_header": [b"# a comment\n#another\n>s1\nATGAAA\n"],
    "hash_inside_record": [b">s1\nATG\n# comment GGG\nAAA\n>s2\nCCC\n"],
    "non_letters_removed": [b">s1\nAT-G 12 aa*a.\n\tccc\n"],
    "nonascii_in_sequence_dropped": [b">s1\nATG\xe9AAA\n"],
    "nonascii_id": [">sé1\nATG\n".encode("utf-8")],
    "empty_header": [b">s1\nATG\n>\nAAA\n"],
    "data_before_header": [b"\n  \nATG\n>s1\nAAA\n"],
    "duplicate_in_file": [b">s1\nATG\n>s1\nAAA\n"],
    "duplicate_across_files": [b">s1\nATG\n", b">s2\nAAA\n>s1\nCCC\n"],
    "two_files": [b">s1\nATG\n", b">s2\nAAA\n>s3\nCCC"],
    "record_without_sequence": [b">s1\n>s2\nATG\n"],
    "no_final_newline": [b">s1\nATG\n>s2\nAAA"],
    "empty_file": [b""],
    "header_tab_desc": [b">s1\tdesc text\nATG\n"],
    # a lone \r only after the first 64 kB: both readers decide line endings from the start of the file
    "late_lone_cr": [b">s0\n" + b"A" * 70000 + b"\n>s1\rATG\n>s2\nCCC\n"],
}


class ReaderCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dump = cpp_dump()
        if cls.dump is None:
            raise unittest.SkipTest("coal-cpp-dump not built (cmake --build build)")
        cls._tmp = tempfile.TemporaryDirectory(prefix="coal_readers_")
        cls.d = Path(cls._tmp.name)
        cls.ref = load_reference("seqio", cls.d)
        if cls.ref is None:
            raise unittest.SkipTest(f"reference reader not available ({REF_BRANCH} not fetched; "
                                    f"or set COAL_REF_SEQIO)")

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def py_protein(self, path: str) -> dict:
        ref = self.ref
        try:
            recs = ref.read_clustal(path) if ref.is_clustal(path) else list(ref.iter_fasta(path))
            return {"records": [list(r) for r in recs]}
        except (ref.InputFormatError, OSError) as e:
            return {"error": str(e)}

    def py_dna(self, paths: list[str]) -> dict:
        ref = self.ref
        try:
            idx = ref.DnaIndex(paths)
            try:
                return {"records": [[idx.ids[i], idx.get(i)] for i in range(len(idx))]}
            finally:
                idx.close()
        except (ref.InputFormatError, OSError) as e:
            return {"error": str(e)}


class TestProteinReaders(ReaderCase):
    def test_every_case(self):
        for name, data in PROTEIN_CASES.items():
            with self.subTest(case=name):
                (self.d / f"{name}.fa").write_bytes(data)
                # relative names, so the error messages (which name the file) are comparable
                cwd = os.getcwd()
                os.chdir(self.d)
                try:
                    want = self.py_protein(f"{name}.fa")
                finally:
                    os.chdir(cwd)
                got = run_dump(self.dump, "protein", f"{name}.fa", cwd=self.d)
                self.assertEqual(got, want)

    def test_missing_file(self):
        want = self.py_protein(str(self.d / "nope.fa"))
        self.assertEqual(run_dump(self.dump, "protein", self.d / "nope.fa"), want)


class TestDnaReaders(ReaderCase):
    def test_every_case(self):
        for name, files in DNA_CASES.items():
            with self.subTest(case=name):
                names = []
                for k, data in enumerate(files):
                    fn = f"{name}_{k}.fna"
                    (self.d / fn).write_bytes(data)
                    names.append(fn)
                cwd = os.getcwd()
                os.chdir(self.d)
                try:
                    want = self.py_dna(names)
                finally:
                    os.chdir(cwd)
                got = run_dump(self.dump, "dna", *names, cwd=self.d)
                self.assertEqual(got, want)

    def test_missing_file(self):
        want = self.py_dna([str(self.d / "nope.fna")])
        self.assertEqual(run_dump(self.dump, "dna", self.d / "nope.fna"), want)


if __name__ == "__main__":
    unittest.main()
