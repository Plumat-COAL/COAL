"""COAL on real biological data (tests/realdata/, see tests/README.md "Real data").

Run:  python -m unittest tests.test_realdata -v        (or: pytest tests/test_realdata.py)

Engine under test, as in test_differential.py: COAL_IMPL=python (default, in-process) or COAL_IMPL=cpp (the coal-cpp
executable; path from COAL_CPP_BIN, default ./build/coal-cpp; skipped if not built).

* roundtrip / partial datasets (published codon alignments from HyPhy's test data): the protein alignment and the
  unaligned DNA are derived from the codon alignment (tests/realdata/realdata_lib.py translates with its own copy of
  the NCBI tables, not COAL's); COAL, given them and the dataset's -codontable, must rebuild the codon alignment
  exactly in all four output formats, and -nogap must follow COAL's documented rule (gap or stop of that code).
* pair datasets (Biopython's codonalign test data: a real protein alignment + real mRNA/genes with UTRs and other
  IDs): COAL's output must translate back to the protein alignment, with the expected number of warnings.
* messy variants of one mitochondrial dataset: what users really feed COAL (UTRs, versioned IDs in another order,
  lower case, CRLF, another aligner's gap placement, the default code on code-2 genes).

Each dataset is skipped if its files are absent. No network, no PAL2NAL needed.
"""
from __future__ import annotations

import io
import json
import os
import random
import re
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RD = ROOT / "tests" / "realdata"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(RD))
from coal.cli import main as coal_main  # noqa: E402
import realdata_lib as L  # noqa: E402

IMPL = os.environ.get("COAL_IMPL", "python").lower()
FORMATS = ("fasta", "clustal", "paml", "codon")
MANIFEST = json.loads((RD / "datasets.json").read_text())["datasets"]


def cpp_binary() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


def run_coal(args: list[str]) -> tuple[int, str, str]:
    if IMPL == "cpp":
        r = subprocess.run([str(cpp_binary()), *args], capture_output=True, timeout=300)
        return r.returncode, r.stdout.decode("ascii", "replace"), \
            r.stderr.decode("ascii", "replace").replace("\r\n", "\n")
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        try:
            rc = coal_main(list(args))
        except SystemExit as e:
            rc = e.code
    return rc, out.getvalue(), err.getvalue()


def data_file(name: str, key: str) -> Path:
    return RD / name / MANIFEST[name][key]


def text_of(path: Path) -> str:
    return L._open_text(path)


class RealDataCase(unittest.TestCase):
    def setUp(self):
        if IMPL == "cpp" and cpp_binary() is None:
            self.skipTest("coal-cpp not built (set COAL_CPP_BIN or build into ./build)")
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name: str, text: str, newline: str = "\n") -> str:
        p = self.dir / name
        p.write_bytes(text.replace("\n", newline).encode())
        return str(p)

    def rows_of(self, name: str) -> list[tuple[str, str]]:
        spec = MANIFEST[name]
        path = data_file(name, "file")
        if not path.is_file():
            self.skipTest(f"real data not present: {path}")
        return L.clean_rows(L.load(path, spec["format"]))


class TestRoundTrip(RealDataCase):
    """Codon alignment -> (protein alignment, unaligned DNA) -> COAL -> the same codon alignment."""

    def check_dataset(self, name: str):
        spec = MANIFEST[name]
        rows = self.rows_of(name)
        code = spec["code"]
        pep, nuc, expected = L.derive(rows, code, partial="gap" if spec["kind"] == "partial" else "error")
        pep_p, nuc_p = self.write("pep.fa", pep), self.write("nuc.fa", nuc)
        for fmt in FORMATS:
            with self.subTest(dataset=name, output=fmt):
                rc, out, err = run_coal([pep_p, nuc_p, "-codontable", str(code), "-output", fmt])
                self.assertEqual(rc, 0, err[-2000:])
                got = L.parse_output(fmt, out)
                if fmt == "fasta":
                    self.assertEqual(got, expected)  # IDs and rows
                else:
                    self.assertEqual([s for _, s in got], [s for _, s in expected])
                    self.assertEqual(len(got), len(expected))
                self.assertNotIn("pepAlnPos", err)  # no per-codon warnings on correct data
        with self.subTest(dataset=name, option="-nogap"):
            rc, out, err = run_coal([pep_p, nuc_p, "-codontable", str(code), "-nogap"])
            self.assertEqual(rc, 0, err[-2000:])
            self.assertEqual(L.fasta_rows(out), L.nogap_rule(expected, code))


for _name, _spec in MANIFEST.items():
    if _spec["kind"] in ("roundtrip", "partial"):
        setattr(TestRoundTrip, f"test_{_name}", lambda self, n=_name: self.check_dataset(n))


class TestRealPairs(RealDataCase):
    """Real protein alignments + real unaligned DNA (Biopython's codonalign test data)."""

    def inputs(self, name: str, mapped: bool = True) -> tuple[str, str, list]:
        spec = MANIFEST[name]
        pep_path, nuc_path = data_file(name, "pep"), data_file(name, "nuc")
        if not (pep_path.is_file() and nuc_path.is_file()):
            self.skipTest(f"real data not present: {pep_path.parent}")
        pep_text, nuc_text = text_of(pep_path), text_of(nuc_path)
        if mapped and "id_map" in spec:  # protein ID -> DNA ID file: rename the DNA records to the protein IDs
            pairs = [line.split() for line in text_of(data_file(name, "id_map")).splitlines() if line.strip()]
            to_pep = {dna_id: pep_id for pep_id, dna_id in pairs}
            nuc_text = re.sub(r"^>(\S+)", lambda m: ">" + to_pep.get(m.group(1), m.group(1)), nuc_text, flags=re.M)
        pep_rows = L.read_clustal(pep_text) if pep_text.startswith("CLUSTAL") else L.parse_fasta(pep_text)
        return self.write("pep.aln", pep_text), self.write("nuc.fa", nuc_text), pep_rows

    def check_pair(self, name: str):
        spec = MANIFEST[name]
        pep_p, nuc_p, pep_rows = self.inputs(name)
        rc, out, err = run_coal([pep_p, nuc_p, "-codontable", str(spec["code"])])
        self.assertEqual(rc, 0, err[-2000:])
        got = L.fasta_rows(out)
        self.assertEqual([n for n, _ in got], [n for n, _ in pep_rows])
        self.assertEqual(err.count("pepAlnPos"), spec["warnings"], err)
        # every codon is its residue, except the ones COAL warned about
        self.assertEqual(L.back_translation_mismatches(pep_rows, got, spec["code"]), spec["warnings"])

    def test_egfr_ids_differ_without_map(self):
        """egfr: protein (NP_/gb) and mRNA (NM_/gb) IDs all differ and the files are in another order, so pairing
        is by file order (as in PAL2NAL): the first protein gets another gene's mRNA and COAL refuses the run."""
        pep_p, nuc_p, _ = self.inputs("egfr", mapped=False)
        rc, out, err = run_coal([pep_p, nuc_p])
        self.assertNotEqual(rc, 0)
        self.assertEqual(out, "")
        self.assertIn("gi|17136534|ref|NP_476758.1|: protein and DNA are inconsistent", err)


for _name, _spec in MANIFEST.items():
    if _spec["kind"] == "pair":
        setattr(TestRealPairs, f"test_{_name}", lambda self, n=_name: self.check_pair(n))


class MessyBase(RealDataCase):
    """Inputs built from one real codon alignment (NAME), changed the way users' files differ."""
    NAME = ""

    def setUp(self):
        super().setUp()
        rows = self.rows_of(self.NAME)
        self.code = MANIFEST[self.NAME]["code"]
        self.rows = [(n.split()[0], s) for n, s in rows]
        self.pep = [(n, "".join(L.translate_codon(s[k:k + 3], self.code) for k in range(0, len(s), 3)))
                    for n, s in self.rows]
        self.dna = [(n, s.replace("-", "")) for n, s in self.rows]

    @staticmethod
    def fa(rows) -> str:
        return "".join(f">{n}\n{s}\n" for n, s in rows)

    def run_rows(self, pep, dna, *opts, newline="\n"):
        return run_coal([self.write("pep.fa", self.fa(pep), newline), self.write("nuc.fa", self.fa(dna), newline),
                         *opts])

    def assert_rebuilt(self, res, expected=None):
        rc, out, err = res
        self.assertEqual(rc, 0, err[-2000:])
        self.assertEqual(L.fasta_rows(out), expected or self.rows)
        self.assertNotIn("pepAlnPos", err)


class TestMessyRealInput(MessyBase):
    """What users really do, on the vertebrate mitochondrial COX1 alignment (code 2)."""
    NAME = "coxI_mammals"

    def test_utr_stop_and_polya_around_each_gene(self):
        rng = random.Random(1)
        flank = lambda: "".join(rng.choice("ACGT") for _ in range(rng.randint(5, 60)))  # noqa: E731
        dna = [(n, flank() + s + "AGA" + flank() + "A" * 20) for n, s in self.dna]
        self.assert_rebuilt(self.run_rows(self.pep, dna, "-codontable", "2"))

    def test_stop_codon_in_dna_and_protein(self):
        res = self.run_rows([(n, p + "*") for n, p in self.pep], [(n, s + "TAA") for n, s in self.dna],
                            "-codontable", "2")
        self.assert_rebuilt(res, [(n, s + "TAA") for n, s in self.rows])

    def test_versioned_protein_ids_dna_in_other_order(self):
        dna = self.dna[::-1]
        rc, out, err = self.run_rows([(n + ".1", p) for n, p in self.pep], dna, "-codontable", "2")
        self.assert_rebuilt((rc, out, err), [(n + ".1", s) for n, s in self.rows])
        self.assertIn("match a DNA ID only after normalising", err)  # a notice, not a per-codon warning

    def test_dna_headers_with_descriptions_in_other_order(self):
        dna = [(f"{n} cytochrome c oxidase subunit I, mitochondrial", s) for n, s in self.dna[::-1]]
        self.assert_rebuilt(self.run_rows(self.pep, dna, "-codontable", "2"))

    def test_lower_case_and_crlf(self):
        rc, out, err = self.run_rows([(n, p.lower()) for n, p in self.pep], [(n, s.lower()) for n, s in self.dna],
                                     "-codontable", "2", newline="\r\n")
        self.assert_rebuilt((rc, out, err))  # fasta_rows upper-cases the rows
        self.assertEqual(L.parse_fasta(out)[0][1], self.rows[0][1].lower())  # DNA case is kept, as in PAL2NAL

    def test_nuclear_code_on_mitochondrial_genes_names_the_right_code(self):
        rc, out, err = self.run_rows(self.pep, self.dna)  # no -codontable: code 1
        self.assertEqual(rc, 0)  # like PAL2NAL: output plus warnings (DEC-1)
        self.assertIn("does not correspond to TGA", err)
        self.assertRegex(err, r"Under genetic code 2 \(Vertebrate mitochondrial\), all 21 sequences match")
        # --strict does not refuse it: far fewer than half of each gene's codons differ (TGA/ATA/AGR only)

    def test_one_gene_truncated(self):
        dna = [(n, s[:-30] if i == 4 else s) for i, (n, s) in enumerate(self.dna)]
        rc, out, err = self.run_rows(self.pep, dna, "-codontable", "2")
        self.assertNotEqual(rc, 0)
        self.assertEqual(out, "")
        self.assertIn(f"{self.dna[4][0]}: protein and DNA are inconsistent", err)
        self.assertIn("The DNA ends before the protein does", err)


class TestOtherAlignerGaps(MessyBase):
    """The same protein alignment with gaps placed differently (as another aligner would), on REM2 (many gaps)."""
    NAME = "rem2_mammals"

    def test_other_aligner_gap_placement(self):
        rng = random.Random(3)

        def shift(p):  # move internal gap runs by one residue: an equally valid alignment from another aligner
            s = list(p)
            for m in list(re.finditer(r"-+", p)):
                a, b = m.start(), m.end()
                if 0 < a and b < len(p) and s[a - 1] != "-" and rng.random() < 0.7:
                    s[a - 1:b] = s[a:b] + [s[a - 1]]
            return "".join(s)
        pep = [(n, shift(p)) for n, p in self.pep]
        self.assertTrue(pep != self.pep, "no gap was moved")
        expected = []
        for (n, p), (_, d) in zip(pep, self.dna):
            it = iter(d[k:k + 3] for k in range(0, len(d), 3))
            expected.append((n, "".join("---" if a == "-" else next(it) for a in p)))
        self.assert_rebuilt(self.run_rows(pep, self.dna, "-codontable", str(self.code)), expected)


if __name__ == "__main__":
    unittest.main()
