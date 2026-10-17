"""Run report (--report / --json, schema coal.report/1) and the plain-language
end-of-run summary with genetic-code hints. Inputs are built in a temp dir."""
from __future__ import annotations

import hashlib
import io
import json
import os
import random
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from coal import __version__  # noqa: E402
from coal.cli import main as coal_main  # noqa: E402
from coal.tables import get_code  # noqa: E402

TOP_KEYS = ["schema", "coal_version", "command", "pal2nal_equivalent", "options", "inputs", "output",
            "summary", "hints", "sequences", "runtime_seconds", "exit_code"]


def run(d: Path, args: list[str]):
    out, err = io.StringIO(), io.StringIO()
    cwd = os.getcwd()
    os.chdir(d)
    try:
        with redirect_stdout(out), redirect_stderr(err):
            try:
                rc = coal_main(args)
            except SystemExit as e:
                rc = e.code
    finally:
        os.chdir(cwd)
    return rc, out.getvalue(), err.getvalue()


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write_genes(d: Path, n: int, codons: int = 60, table: int = 1, seed: int = 3, mito: bool = False,
                broken: set[int] = frozenset(), wrong: set[int] = frozenset()):
    """n genes; mito=True uses vertebrate-mitochondrial codons (TGA = W, ATA = M)."""
    rng = random.Random(seed)
    code = get_code(table)
    aas = sorted(a for a in code.codons if a not in "*")
    pep, nuc = [], []
    for i in range(n):
        prot = ["M"] + [rng.choice(aas) for _ in range(codons - 1)]
        dna = []
        for k, aa in enumerate(prot):
            if k == 0:
                c = "ATG"
            elif mito and aa == "W":
                c = "TGA"
            elif mito and aa == "M":
                c = "ATA"
            else:
                c = rng.choice(sorted(code.codons[aa]))
            dna.append(c)
        if mito and "W" not in prot[1:]:
            prot[5], dna[5] = "W", "TGA"   # every mito gene needs at least one TGA
        if i in wrong:
            dna[20] = "TGG" if prot[20] != "W" else "AAA"
        seq = "".join(dna)
        if i in broken:
            seq = seq[:40] + seq[41:]          # one-base deletion: no fit in any frame
        pep.append(f">g{i}\n{''.join(prot)}\n")
        nuc.append(f">g{i}\n{seq}\n")
    (d / "pep.fasta").write_text("".join(pep), newline="\n")
    (d / "nuc.fasta").write_text("".join(nuc), newline="\n")


class TestReport(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def report(self):
        return json.loads((self.d / "r.json").read_text())

    def test_clean_run_report(self):
        write_genes(self.d, 5)
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "-o", "out.paml", "-output", "paml",
                                    "--report", "r.json"])
        self.assertEqual(rc, 0, err)
        r = self.report()
        for k in TOP_KEYS:
            self.assertIn(k, r)
        self.assertEqual(r["schema"], "coal.report/1")
        self.assertEqual(r["coal_version"], __version__)
        # MRG-2: "command" reproduces the output (as the app writes it); the exact argv is "command_line"
        self.assertEqual(r["command"], "coal pep.fasta nuc.fasta -output paml -o out.paml")
        self.assertEqual(r["command_line"], "coal pep.fasta nuc.fasta -o out.paml -output paml --report r.json")
        self.assertEqual(r["pal2nal_equivalent"], "pal2nal.pl pep.fasta nuc.fasta -output paml")
        self.assertEqual(r["options"], {"output_format": "paml", "codontable": 1, "nogap": False,
                                        "nomismatch": False, "keep_going": False, "strict": False,
                                        "blockonly": False})
        self.assertEqual([(i["role"], i["path"], i["sequences"]) for i in r["inputs"]],
                         [("protein", "pep.fasta", 5), ("dna", "nuc.fasta", 5)])
        self.assertEqual(r["inputs"][0]["sha256"], sha(self.d / "pep.fasta"))
        self.assertEqual(r["output"], {"path": "out.paml", "format": "paml",
                                       "sha256": sha(self.d / "out.paml"), "sequences": 5})
        self.assertEqual({k: r["summary"][k] for k in ("total", "ok", "warning", "failed", "pairing")},
                         {"total": 5, "ok": 5, "warning": 0, "failed": 0, "pairing": "id"})
        self.assertEqual(r["hints"], [])
        self.assertEqual([s["status"] for s in r["sequences"]], ["ok"] * 5)
        self.assertEqual(r["sequences"][0]["dna_id"], "g0")
        self.assertEqual(r["sequences"][0]["pairing"], "id")
        self.assertEqual(r["exit_code"], 0)
        self.assertIsInstance(r["runtime_seconds"], float)
        self.assertNotIn("coal: summary:", err)       # no problems -> no summary

    def test_codontable_in_pal2nal_equivalent(self):
        write_genes(self.d, 2, table=2)
        rc, _, err = run(self.d, ["pep.fasta", "nuc.fasta", "-codontable", "2", "-nogap", "-o", "o.fa",
                                  "--report", "r.json"])
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.report()["pal2nal_equivalent"],
                         "pal2nal.pl pep.fasta nuc.fasta -output fasta -codontable 2 -nogap")

    def test_json_to_stdout_requires_o(self):
        write_genes(self.d, 2)
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "--json"])
        self.assertEqual(rc, 2)
        self.assertIn("--json requires -o", err)
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "--json", "-o", "o.fa"])
        self.assertEqual(rc, 0)
        self.assertEqual(json.loads(out)["summary"]["total"], 2)

    def test_mismatch_is_a_warning(self):
        write_genes(self.d, 3, wrong={1})
        rc, _, err = run(self.d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
        self.assertEqual(rc, 0)
        r = self.report()
        s = r["sequences"][1]
        self.assertEqual(s["status"], "warning")
        self.assertEqual(s["messages"][0]["code"], "mismatch")
        self.assertEqual(s["messages"][0]["column"], 21)
        self.assertIn("Column 21:", s["messages"][0]["text"])
        self.assertEqual(r["summary"]["warning"], 1)
        self.assertIn("coal: summary: 1 of 3 sequences has codons that do not match their amino acid (kept in the "
                      "alignment).", err)

    def test_failed_run_still_writes_report(self):
        write_genes(self.d, 4, broken={2})
        rc, out, err = run(self.d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
        self.assertEqual(rc, 1)
        self.assertFalse((self.d / "o.fa").exists())
        r = self.report()
        self.assertEqual(r["exit_code"], 1)
        self.assertEqual(r["error"]["code"], "inconsistent")  # MRG-2: the app's message codes
        failed = [s for s in r["sequences"] if s["status"] == "failed"]
        self.assertEqual([s["id"] for s in failed], ["g2"])
        # MRG-2: where the gene breaks comes first (as the app's locate_breaks), then the engine's reason
        self.assertEqual([m["code"] for m in failed[0]["messages"]], ["frameshift", "inconsistent"])
        self.assertIn("Frameshift at column 14", failed[0]["messages"][0]["text"])
        # g0 and g1 were checked (not written), g3 was never reached
        self.assertEqual([s["status"] for s in r["sequences"]], ["ok", "ok", "failed", "not_checked"])
        self.assertIsNone(r["output"])
        self.assertEqual(r["summary"]["failed"], 1)
        self.assertIn("coal: summary: 1 sequence of 4 could not be converted (first: g2:", err)

    def test_input_error_report(self):
        write_genes(self.d, 3)
        (self.d / "nuc.fasta").write_text(">g0\nATG\n>g1\nATG\n")
        rc, _, err = run(self.d, ["pep.fasta", "nuc.fasta", "--report", "r.json"])
        self.assertEqual(rc, 1)
        r = self.report()
        self.assertEqual(r["error"]["code"], "count_mismatch")
        self.assertIn("number of sequences differs (protein: 3, DNA: 2)", r["error"]["text"])
        # MRG-2: as the app's report: every protein listed, none checked, no pairing
        self.assertEqual([(s["id"], s["status"], s["dna_id"]) for s in r["sequences"]],
                         [(f"g{i}", "not_checked", None) for i in range(3)])
        self.assertIsNone(r["summary"]["pairing"])
        self.assertEqual(r["errors"][0]["code"], "count_mismatch")
        self.assertEqual(r["inputs"][1]["sequences"], 2)

    def test_example_file_is_valid(self):
        ex = json.loads((ROOT / "tests" / "report_example.json").read_text())
        for k in TOP_KEYS:
            self.assertIn(k, ex)
        self.assertEqual(ex["schema"], "coal.report/1")
        for s in ex["sequences"]:
            self.assertIn(s["status"], ("ok", "warning", "failed"))
            self.assertIn(s["pairing"], ("id", "normalised", "order"))


class TestSummaryAndHints(unittest.TestCase):
    def test_mito_genes_with_standard_code_suggest_code_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_genes(d, 30, table=2, mito=True)
            rc, _, err = run(d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
            r = json.loads((d / "r.json").read_text())
            hints = [h for h in r["hints"] if h["code"] == "genetic_code"]
            self.assertEqual(len(hints), 1, err)
            self.assertEqual(hints[0]["suggested_codontable"], 2)
            self.assertIn("vertebrate mitochondrial", hints[0]["text"])
            self.assertIn("Re-run with code 2 (-codontable 2).", err)
            self.assertIn("coal: summary: 30 of 30 sequences have codons that do not match their amino acid (kept "
                          "in the alignment). With genetic code 2 (vertebrate mitochondrial), all 30 sequences match "
                          "their protein codon for codon; with code 1 (standard), the one used, none do. Re-run with "
                          "code 2 (-codontable 2).", err)
            self.assertEqual((hints[0]["checked"], hints[0]["fit"], hints[0]["fit_current"]), (30, 30, 0))
            # with the right code: clean run, no summary, no hint
            rc, _, err = run(d, ["pep.fasta", "nuc.fasta", "-codontable", "2", "-o", "o.fa", "--report", "r.json"])
            self.assertEqual(rc, 0)
            self.assertEqual(json.loads((d / "r.json").read_text())["hints"], [])
            self.assertNotIn("coal: summary:", err)

    def test_no_hint_below_threshold(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_genes(d, 10, wrong={4})          # MRG-2: no genetic code fixes this codon: no hint
            rc, _, err = run(d, ["pep.fasta", "nuc.fasta", "-o", "o.fa", "--report", "r.json"])
            self.assertEqual(json.loads((d / "r.json").read_text())["hints"], [])
            self.assertIn("coal: summary: 1 of 10 sequences has codons that do not match their amino acid (kept in "
                          "the alignment).", err)

    def test_quiet_suppresses_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_genes(d, 3, wrong={0})
            rc, _, err = run(d, ["pep.fasta", "nuc.fasta", "-q"])
            self.assertEqual(err, "")


if __name__ == "__main__":
    unittest.main()
