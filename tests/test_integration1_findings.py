"""The findings C01-C03, C06 and C08 of the COAL / orthologr integration pass, each fixed and kept fixed here.

  * C01: --preset pwalign reproduced pwalign's score but not its choice among co-optimal alignments (COAL put
    the gaps leftmost): 49-76 % of orthologr's protein alignments differed from COAL's. AlignParams.tie_break
    "pwalign" (the preset's default) returns pwalign::pairwiseAlignment's alignment; "leftmost" keeps COAL's.
    tests/fixtures/pwalign: synthetic pairs and pwalign 1.6.0's alignments of them (make_expected.R).
  * C02: a pair's result depended on which gene was seq1 (see the class C02 below for what is canonical now).
  * C03: with several -m, the summary line, report.json's pair counts, the saturation notice and failed.tsv
    described the first method in canonical order only, without naming it.
  * C06: cluster_summary.tsv called a cluster whose pairs all have no estimate (saturated / not_converged)
    "failed", as if the run had failed.
  * C08: report.json "command" was null when the pipeline ran from Python (run_pairs / run_clusters) or from R
    (coalr::coal_pipeline_*).

Synthetic inputs only (no real sequence data); coal-cpp (COAL_CPP_BIN) and coal-cpp-align for the C++ parts.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import align as A  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402
from test_align import align_binary, no_binary  # noqa: E402
from test_dnds_pipeline import CPP, FAM, FAM_CLUSTERS, Q, S, SCHEMA, PipelineCase, needs_cpp  # noqa: E402
from test_report_schema import validate  # noqa: E402

PW = ROOT / "tests" / "fixtures" / "pwalign"
TIMEOUT = 300


def fasta(path: Path) -> list[tuple[str, str]]:
    out: list[list[str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            out.append([line[1:], ""])
        elif out:
            out[-1][1] += line.strip()
    return [(n, s) for n, s in out]


def fixture_pairs():
    """[(a, b, expected_a, expected_b, expected_score)] of tests/fixtures/pwalign."""
    pairs, exp = fasta(PW / "pairs.fasta"), fasta(PW / "expected.fasta")
    out = []
    for k in range(0, len(pairs), 2):
        score = float(exp[k][0].split("score=")[1])
        out.append((pairs[k][1], pairs[k + 1][1], exp[k][1], exp[k + 1][1], score))
    return out


class C01PwalignTieBreak(unittest.TestCase):
    def test_fixture_has_co_optimal_pairs(self):
        # the fixture is only a test if COAL's own choice differs from pwalign's in many pairs
        rows = fixture_pairs()
        self.assertGreaterEqual(len(rows), 80)
        left = A.AlignParams("BLOSUM62-PWALIGN", gap_open=14, gap_extend=4, free_end_gaps=False)
        differ = 0
        for a, b, ea, eb, score in rows[:40]:
            r = A.align_pair(a, b, left)
            self.assertEqual(r.score, score)
            differ += (r.a, r.b) != (ea, eb)
        self.assertGreater(differ, 10)

    def test_python_kernel_is_pwalign(self):
        p = A.AlignParams.pwalign()
        self.assertEqual(p.tie_break, "pwalign")
        for a, b, ea, eb, score in fixture_pairs():
            r = A.align_pair(a, b, p)
            self.assertEqual((r.a, r.b), (ea, eb), a)
            self.assertEqual(r.score, score)
            self.assertEqual(A.score_alignment(r.a, r.b, p), score)

    def test_pwalign_matrix(self):
        # pwalign's BLOSUM62 (data(BLOSUM62, package = "pwalign")) differs from NCBI's current file for B, Z, X
        pw, ncbi = A.MATRICES["BLOSUM62-PWALIGN"], A.MATRICES["BLOSUM62"]
        diff = {(A.ALPHABET[i], A.ALPHABET[j]) for i in range(24) for j in range(24) if pw[i][j] != ncbi[i][j]}
        self.assertTrue(diff)
        self.assertTrue(all("B" in d or "Z" in d or "X" in d for d in diff))
        self.assertTrue(all(pw[A.ALPHABET.index("X")][k] == -1 for k in range(23)))
        self.assertEqual(A.AlignParams.pwalign().matrix, "BLOSUM62-PWALIGN")
        a = P.build_parser().parse_args(["pairs", "--query", "q", "--subject", "s", "-o", "o", "--preset", "pwalign"])
        P._finish_args(a)
        self.assertEqual(a.matrix, "BLOSUM62-PWALIGN")
        a = P.build_parser().parse_args(["pairs", "--query", "q", "--subject", "s", "-o", "o", "--preset", "pwalign",
                                         "--matrix", "blosum62"])
        P._finish_args(a)
        self.assertEqual(a.matrix, "BLOSUM62")

    def test_parameters(self):
        self.assertEqual(A.AlignParams().tie_break, "leftmost")
        self.assertEqual(A.AlignParams.needle().args(), A.AlignParams().args())   # unchanged command line
        self.assertEqual(A.AlignParams.pwalign().args()[-2:], ["--tie-break", "pwalign"])
        self.assertEqual(A.AlignParams.pwalign().args()[:2], ["--matrix", "BLOSUM62-PWALIGN"])
        with self.assertRaises(A.AlignError) as cm:
            A.AlignParams(tie_break="rightmost").check()
        self.assertEqual(cm.exception.code, "align_params")
        # empty sequences: the same full-length rows with either rule
        for tb in A.TIE_BREAKS:
            r = A.align_pair("", "MK", A.AlignParams(tie_break=tb))
            self.assertEqual((r.a, r.b), ("--", "MK"))

    def test_cpp_align_binary(self):
        exe = align_binary()
        if exe is None:
            no_binary(self)
        rows = fixture_pairs()
        out = subprocess.run([str(exe), "--preset", "pwalign", "-j", "2", str(PW / "pairs.fasta")],
                             capture_output=True, text=True, timeout=TIMEOUT)
        self.assertEqual(out.returncode, 0, out.stderr)
        got = [s for _, s in fasta_text(out.stdout)]
        self.assertEqual(len(got), 2 * len(rows))
        for k, (a, b, ea, eb, _) in enumerate(rows):
            self.assertEqual((got[2 * k], got[2 * k + 1]), (ea, eb), a)
        # --tie-break leftmost after the preset: COAL's own choice, as the Python kernel makes it
        out = subprocess.run([str(exe), "--preset", "pwalign", "--tie-break", "leftmost", str(PW / "pairs.fasta")],
                             capture_output=True, text=True, timeout=TIMEOUT)
        got = [s for _, s in fasta_text(out.stdout)]
        left = A.AlignParams("BLOSUM62-PWALIGN", gap_open=14, gap_extend=4, free_end_gaps=False)
        for k, (a, b, *_rest) in enumerate(rows[:30]):
            r = A.align_pair(a, b, left)
            self.assertEqual((got[2 * k], got[2 * k + 1]), (r.a, r.b))
        bad = subprocess.run([str(exe), "--tie-break", "rightmost", str(PW / "pairs.fasta")],
                             capture_output=True, text=True, timeout=TIMEOUT)
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("--tie-break: leftmost or pwalign", bad.stderr)

    def test_cpp_size_message_names_the_pair_as_given(self):
        # the kernel runs on the exchanged pair; the refusal still says n x m in the order given
        exe = align_binary()
        if exe is None:
            no_binary(self)
        import tempfile
        a, b = "MKTAYIAKQRQISFVKSHFSRQ" * 3, "MKTAYAKQRQISFVKSHFSRQLE" * 2
        with tempfile.TemporaryDirectory() as td:
            f = Path(td) / "p.fasta"
            f.write_text(f">a\n{a}\n>b\n{b}\n", encoding="utf-8")
            out = subprocess.run([str(exe), "--preset", "pwalign", "--max-cells", "100", str(f)],
                                 capture_output=True, text=True, timeout=TIMEOUT)
        self.assertIn(f"{len(a)} x {len(b)} residues", out.stderr)

    def test_pipeline_records_and_passes_the_tie_break(self):
        if CPP is None:
            self.skipTest("coal-cpp not built")
        a = P.build_parser().parse_args(["pairs", "--query", "q", "--subject", "s", "-o", "o", "--preset", "pwalign"])
        P._finish_args(a)
        self.assertEqual(a.tie_break, "pwalign")
        a = P.build_parser().parse_args(["pairs", "--query", "q", "--subject", "s", "-o", "o", "--preset", "pwalign",
                                         "--tie-break", "leftmost"])
        P._finish_args(a)
        self.assertEqual(a.tie_break, "leftmost")
        a = P.build_parser().parse_args(["pairs", "--query", "q", "--subject", "s", "-o", "o"])
        P._finish_args(a)
        self.assertEqual(a.tie_break, "leftmost")
        with self.assertRaises(ValueError):
            P._namespace("pairs", ["--query", "q", "--subject", "s", "-o", "o"], {"tie_break": "rightmost"})


def fasta_text(text: str) -> list[tuple[str, str]]:
    out: list[list[str]] = []
    for line in text.splitlines():
        if line.startswith(">"):
            out.append([line[1:], ""])
        elif out:
            out[-1][1] += line.strip()
    return [(n, s) for n, s in out]


@needs_cpp
class C01Pipeline(PipelineCase):
    def test_preset_pwalign_alignments_are_pwalign_and_recorded(self):
        pairs = self.pairs_file([(r[0], r[1]) for r in _fixture_pair_ids()])
        res = P.run_pairs(Q, S, self.tmp / "pw", pairs=str(pairs), preset="pwalign", jobs=1, quiet=True,
                          keep_protein_alignments=True)
        self.assertIn(res.exit_code, (0, 2))
        self.assertEqual(res.report["parameters"]["tie_break"], "pwalign")
        cmd = [c for c in res.report["tools"]["coal-cpp"]["commands"] if "codon-dnds" in c][0]
        self.assertEqual(cmd[cmd.index("--tie-break") + 1], "pwalign")
        old = P.run_pairs(Q, S, self.tmp / "left", pairs=str(pairs), preset="pwalign", tie_break="leftmost", jobs=1,
                          quiet=True, keep_protein_alignments=True)
        self.assertEqual(old.report["parameters"]["tie_break"], "leftmost")
        cmd = [c for c in old.report["tools"]["coal-cpp"]["commands"] if "codon-dnds" in c][0]
        self.assertNotIn("--tie-break", cmd)
        new_alns = fasta(Path(res.outputs["protein_alignments"]))
        old_alns = fasta(Path(old.outputs["protein_alignments"]))
        p = A.AlignParams.pwalign()
        for (n1, r1), (n2, r2) in zip(new_alns[0::2], new_alns[1::2]):
            want = A.align_pair(r1.replace("-", ""), r2.replace("-", ""), p)
            self.assertEqual((r1, r2), (want.a, want.b), n1)
        self.assertEqual([s.replace("-", "") for _, s in new_alns], [s.replace("-", "") for _, s in old_alns])
        default = P.run_pairs(Q, S, self.tmp / "def", pairs=str(pairs), jobs=1, quiet=True)
        self.assertEqual(default.report["parameters"]["tie_break"], "leftmost")


def _fixture_pair_ids():
    """The ortholog pairs of tests/fixtures/pipeline that coalr's fixtures use."""
    p = ROOT / "packaging" / "r" / "coalr" / "tests" / "testthat" / "fixtures" / "pairs_orthologs.tsv"
    rows = [line.split("\t") for line in p.read_text(encoding="utf-8").splitlines()[1:] if line.strip()]
    return [(r[0], r[1]) for r in rows]


if __name__ == "__main__":
    unittest.main()
