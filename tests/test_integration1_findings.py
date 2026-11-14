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


SENSE = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG" if a + b + c not in ("TAA", "TAG", "TGA")]


def related_cds(rng, n_pairs=12):
    """Pairs of related CDS with substitutions and indels in repeats: co-optimal protein alignments, and Comeron
    values that depend on which gene is first."""
    out = []
    for k in range(n_pairs):
        unit = [rng.choice(SENSE) for _ in range(2)]
        a = ["ATG"] + [rng.choice(SENSE) if rng.random() < 0.7 else unit[rng.random() < 0.5] for _ in range(150)]
        b = list(a)
        for i in range(1, len(b)):
            if rng.random() < 0.2:
                b[i] = rng.choice(SENSE)
        del b[40:40 + 1 + k % 3]
        b[90:90] = [unit[0]] * (1 + k % 2)
        out.append(("".join(a) + "TAA", "".join(b) + "TGA"))
    return out


@needs_cpp
class C02InputOrder(PipelineCase):
    """A pair is computed in input order (the gene of the earlier --cds file, then the earlier record, first) in
    every route: `pipeline clusters` lists the representative first but aligns and estimates the pair as the
    pairs route and orthologr do."""

    def files(self):
        import random
        rng = random.Random(7)
        pairs = related_cds(rng)
        f1 = self.tmp / "sp1.fa"
        f2 = self.tmp / "sp2.fa"
        f1.write_text("".join(f">y{k}\n{b}\n" for k, (_, b) in enumerate(pairs)), encoding="utf-8")
        f2.write_text("".join(f">x{k}\n{a}\n" for k, (a, _) in enumerate(pairs)), encoding="utf-8")
        cl = self.tmp / "clusters.tsv"   # the representative (x, from the second file) first
        cl.write_text("".join(f"x{k}\tx{k}\nx{k}\ty{k}\n" for k in range(len(pairs))), encoding="utf-8")
        return f1, f2, cl, len(pairs)

    @staticmethod
    def values(path, key_cols=("seq1", "seq2")):
        import csv
        with open(path, encoding="utf-8") as fh:
            return {(frozenset((r["seq1"], r["seq2"])), r["method"]): (r["dN"], r["dS"], r["status"])
                    for r in csv.DictReader(fh, delimiter="\t")}

    def test_clusters_equal_pairs_in_input_order(self):
        f1, f2, cl, n = self.files()
        methods = ["YN", "LWL", "Comeron"]
        for preset in ("needle", "pwalign"):
            clu = P.run_clusters([f1, f2], self.tmp / f"clu_{preset}", clusters=cl, method=methods, jobs=1,
                                 quiet=True, preset=preset)
            self.assertIn(clu.exit_code, (0, 2))
            rows = (self.tmp / f"clu_{preset}" / "dnds.tsv").read_text(encoding="utf-8").splitlines()[1:]
            self.assertTrue(all(r.split("\t")[1].startswith("x") for r in rows))   # still the representative first
            pp = self.pairs_file([(f"y{k}", f"x{k}") for k in range(n)])
            pr = P.run_pairs(f1, f2, self.tmp / f"pairs_{preset}", pairs=str(pp), method=methods, jobs=1, quiet=True,
                             preset=preset)
            self.assertIn(pr.exit_code, (0, 2))
            c, p = self.values(clu.outputs["dnds"]), self.values(pr.outputs["dnds"])
            self.assertEqual(len(p), 3 * n)
            self.assertEqual(c, p, preset)
            # the other input order: the representative's file first, no exchange, the pairs route with x first
            clu2 = P.run_clusters([f2, f1], self.tmp / f"clu2_{preset}", clusters=cl, method=methods, jobs=1,
                                  quiet=True, preset=preset)
            pp2 = self.pairs_file([(f"x{k}", f"y{k}") for k in range(n)])
            pr2 = P.run_pairs(f2, f1, self.tmp / f"pairs2_{preset}", pairs=str(pp2), method=methods, jobs=1,
                              quiet=True, preset=preset)
            c2, p2 = self.values(clu2.outputs["dnds"]), self.values(pr2.outputs["dnds"])
            self.assertEqual(c2, p2, preset)
            # the test can tell: Comeron depends on which gene is first
            self.assertTrue(any(c[k] != c2[k] for k in c if k[1] == "Comeron"))

    def test_engine_swap_field(self):
        f1, f2, cl, n = self.files()
        x = f2.read_text(encoding="utf-8").split()[1]
        y = f1.read_text(encoding="utf-8").split()[1]

        def run(text):
            out = self.tmp / "o.tsv"
            r = subprocess.run([str(CPP), "codon-dnds", "-o", str(out), "-m", "YN", "-m", "Comeron", "-q"],
                               input=text, capture_output=True, text=True, timeout=TIMEOUT)
            return r, (out.read_text(encoding="utf-8").splitlines()[1:] if r.returncode in (0, 2) else [])
        _, swapped = run(f"G\tg\tall\tswap\nM\tx\t{x}\nM\ty\t{y}\nE\n")
        _, yx = run(f"G\tg\tall\nM\ty\t{y}\nM\tx\t{x}\nE\n")
        _, plain = run(f"G\tg\tall\tnocodon\nM\tx\t{x}\nM\ty\t{y}\nE\n")
        self.assertEqual([r.split("\t")[:2] for r in swapped], [["x", "y"]] * 2)
        self.assertEqual([r.split("\t")[2:] for r in swapped], [r.split("\t")[2:] for r in yx])
        self.assertNotEqual(swapped, plain)
        bad, _ = run(f"G\tg\tall\tswapped\nM\tx\t{x}\nM\ty\t{y}\nE\n")
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("nocodon or swap", bad.stderr)


class _StubRun:
    def __init__(self):
        self.counts, self.failed, self.warnings = {}, [], []
        self.a = None

    def fail(self, level, code, message, cluster="", seq1="", seq2=""):
        self.failed.append((level, code, message, cluster, seq1, seq2))

    def warn(self, text):
        self.warnings.append(text)


class C03SeveralMethods(unittest.TestCase):
    HEAD = "seq1\tseq2\tmethod\tdN\tdS\tdNdS\tstatus\tnote\n"

    def table(self, rows, tmp):
        p = Path(tmp) / "dnds.tsv"
        p.write_text(self.HEAD + "".join("\t".join(r) + "\n" for r in rows), encoding="utf-8")
        return str(p)

    def test_counts_failures_and_notice_per_method(self):
        import tempfile
        rows = []
        for i in range(30):
            sat_yn = i < 25           # YN: 25 of 30 saturated, LWL: 5 of 30
            rows.append((f"a{i}", f"b{i}", "LWL", "0.1", "NA" if i >= 25 else "0.2", "NA",
                         "saturated" if i >= 25 else "ok", ""))
            rows.append((f"a{i}", f"b{i}", "YN", "0.1", "NA" if sat_yn else "0.2", "NA",
                         "saturated" if sat_yn else "ok", ""))
        rows.append(("x", "y", "LWL", "NA", "NA", "NA", "error", "protein alignment failed: too large"))
        rows.append(("x", "y", "YN", "NA", "NA", "NA", "error", "protein alignment failed: too large"))
        rows.append(("u", "v", "LWL", "0.1", "0.2", "0.5", "ok", ""))
        rows.append(("u", "v", "YN", "NA", "NA", "NA", "error", "made-up YN-only error"))
        with tempfile.TemporaryDirectory() as td:
            run = _StubRun()
            P._count_pairs(run, self.table(rows, td), False)
        c = run.counts
        self.assertEqual(c["pairs"]["method"], "LWL")
        self.assertEqual(set(c["pairs_by_method"]), {"LWL", "YN"})
        self.assertEqual(c["pairs_by_method"]["YN"]["saturated"], 25)
        self.assertEqual(c["pairs_by_method"]["YN"]["error"], 2)
        self.assertEqual(c["pairs"]["saturated"], 5)
        # one entry for the pair that failed in every method, the YN-only error named
        self.assertEqual([(f[4], f[5], f[2]) for f in run.failed],
                         [("x", "y", "protein alignment failed: too large"), ("u", "v", "YN: made-up YN-only error")])
        # the notice names each method (only YN is above 30%, all are listed)
        self.assertEqual(len(run.warnings), 1)
        self.assertIn("LWL 5 of 32 (16%), YN 25 of 32 (78%)", run.warnings[0])

    def test_one_method_as_before(self):
        import tempfile
        rows = [(f"a{i}", f"b{i}", "YN", "0.1", "NA", "NA", "saturated", "") for i in range(25)]
        with tempfile.TemporaryDirectory() as td:
            run = _StubRun()
            P._count_pairs(run, self.table(rows, td), False)
        self.assertNotIn("method", run.counts["pairs"])
        self.assertNotIn("pairs_by_method", run.counts)
        self.assertTrue(run.warnings[0].startswith("25 of the 25 pairs (100%) are saturated or not_converged: dS"))


@needs_cpp
class C03Pipeline(PipelineCase):
    def test_report_and_summary_name_the_methods(self):
        import contextlib
        import io
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            res = P.run_pairs(Q, S, self.tmp / "o", method=["YN", "NG", "LWL"], jobs=1)
        self.assertIn(res.exit_code, (0, 2))
        self.check_report(res)
        by = res.report["counts"]["pairs_by_method"]
        self.assertEqual(list(by), ["NG", "LWL", "YN"])   # the order of the rows
        self.assertEqual(res.report["counts"]["pairs"]["method"], "NG")
        line = [x for x in err.getvalue().splitlines() if "pairs: NG " in x]
        self.assertEqual(len(line), 1, err.getvalue())
        self.assertIn("; LWL ", line[0])
        self.assertIn("; YN ", line[0])


class C06ClusterStatus(unittest.TestCase):
    def test_status_and_note(self):
        import csv
        import tempfile
        head = "cluster_id\tseq1\tseq2\tmethod\tdN\tdS\tdNdS\tstatus\tnote\n"
        rows = [("c1", "a", "b", "YN", "0.5", "NA", "NA", "saturated"), ("c1", "a", "c", "YN", "NA", "NA", "NA",
                                                                          "not_converged"),
                ("c2", "d", "e", "YN", "NA", "NA", "NA", "not_converged"),
                ("c3", "f", "g", "YN", "NA", "NA", "NA", "no_data"),
                ("c4", "h", "i", "YN", "NA", "NA", "NA", "error"), ("c4", "h", "j", "YN", "0.5", "NA", "NA",
                                                                     "saturated"),
                ("c5", "k", "l", "YN", "0.1", "0.2", "0.5", "ok"), ("c5", "k", "m", "YN", "0.5", "NA", "NA",
                                                                    "saturated"),
                ("c6", "n", "o", "YN", "0.5", "NA", "NA", "saturated")]
        plans = [{"cluster_id": c, "representative": c, "n_members": 3, "n_used": 3, "pair_selection": "all",
                  "aligner": "nw", "note": "an MSA note" if c == "c6" else ""} for c in ("c1", "c2", "c3", "c4", "c5",
                                                                                         "c6")]
        with tempfile.TemporaryDirectory() as td:
            d, o = Path(td) / "dnds.tsv", Path(td) / "summary.tsv"
            d.write_text(head + "".join("\t".join(r) + "\t\n" for r in rows), encoding="utf-8")
            P.summarise(str(d), plans, str(o))
            with open(o, encoding="utf-8") as fh:
                got = {r["cluster_id"]: (r["status"], r["note"]) for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(got["c1"], ("saturated", "no pair has an estimate (2 pairs: 1 saturated, 1 not converged)"))
        self.assertEqual(got["c2"], ("not_converged", "no pair has an estimate (1 pair: 1 not converged)"))
        self.assertEqual(got["c3"], ("no_data", "no pair has an estimate (1 pair: 1 without codons to compare)"))
        self.assertEqual(got["c4"], ("failed", "no pair has an estimate (2 pairs: 1 saturated, 1 failed)"))
        self.assertEqual(got["c5"], ("ok", ""))
        self.assertEqual(got["c6"], ("saturated", "an MSA note"))
