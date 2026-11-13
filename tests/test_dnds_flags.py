"""--min-codons N / --max-ds X: flags for estimates that need care (docs/kaks.md "Filtering before interpretation").

In `coal dnds`, `coal codon-dnds` and `coal pipeline`: a result with an estimate (status ok or warning) and fewer than
N compared codons, or a dS above X, gets a note and status `warning`. No row is dropped and no number changes: the
table with the flags is the table without them but for `status` and `note` of the flagged rows. Off by default (the
output is then byte-identical to a run without the options). Needs build/coal-cpp (COAL_CPP_BIN).
"""
from __future__ import annotations

import io
import os
import subprocess
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import dnds_pipeline as P  # noqa: E402
from test_dnds_pipeline import (CPP, FAM, FAM_CLUSTERS, PLAIN, Q, S, PipelineCase, needs_cpp,  # noqa: E402
                                read_tsv, run_cpp)

AXT = ROOT / "tests" / "kaks" / "comeron_li_cases.axt"   # 499 pairs of many lengths and divergences


def differing(plain: list[dict], flagged: list[dict]) -> list[int]:
    """Rows that differ; asserts that nothing but status and note ever differs."""
    assert len(plain) == len(flagged)
    out = []
    for k, (a, b) in enumerate(zip(plain, flagged)):
        for col in a:
            if col not in ("status", "note"):
                assert a[col] == b[col], (k, col, a[col], b[col])
        if (a["status"], a["note"]) != (b["status"], b["note"]):
            out.append(k)
    return out


@needs_cpp
class DndsFlags(PipelineCase):
    def table(self, *args):
        out = self.tmp / "o.tsv"
        r = run_cpp("dnds", AXT, "-o", out, *args)
        self.assertEqual(r.returncode, 0, r.stderr)
        return read_tsv(out), r.stderr.decode()

    def test_off_by_default_and_zero_is_off(self):
        a = run_cpp("dnds", AXT, "-q", "-m", "YN", "-m", "NG")
        b = run_cpp("dnds", AXT, "-q", "-m", "YN", "-m", "NG", "--min-codons", "0")
        self.assertEqual(a.returncode, 0, a.stderr)
        self.assertEqual(a.stdout, b.stdout)
        self.assertNotIn(b"codons compared", a.stdout)

    def test_flags_set_status_and_note_only(self):
        plain, _ = self.table("-m", "YN", "-m", "NG")
        flagged, err = self.table("-m", "YN", "-m", "NG", "--min-codons", "60", "--max-ds", "1.5")
        changed = set(differing(plain, flagged))
        n_short = n_high = 0
        for k, (a, b) in enumerate(zip(plain, flagged)):
            has = a["status"] in ("ok", "warning")
            short = has and int(a["codons_compared"]) < 60
            high = has and a["dS"] != "NA" and float(a["dS"]) > 1.5
            n_short += short
            n_high += high
            want = [x for x in (a["note"],) if x] + (["fewer than 60 codons compared"] if short else []) + \
                (["dS above 1.5"] if high else [])
            self.assertEqual(b["note"], "; ".join(want), k)
            self.assertEqual(b["status"], "warning" if (short or high) else a["status"], k)
            self.assertEqual(k in changed, bool(short or high), k)
        self.assertGreater(n_short, 5)
        self.assertGreater(n_high, 5)
        self.assertLess(len(changed), len(plain))
        # rows without an estimate keep their status and note
        self.assertTrue(any(a["status"] in ("saturated", "no_data", "not_converged", "error") for a in plain))
        self.assertIn(f"coal dnds: {n_short + n_high} flag(s) set (status warning, numbers unchanged): {n_short} "
                      f"result(s) with fewer than 60 codons compared, {n_high} with dS above 1.5", err)

    def test_one_flag_at_a_time_and_number_format(self):
        flagged, err = self.table("--max-ds", "2")
        self.assertTrue(any(r["note"].endswith("dS above 2") for r in flagged))
        self.assertNotIn("codons compared", err)
        self.assertRegex(err, r"flag\(s\) set \(status warning, numbers unchanged\): \d+ with dS above 2\n")
        flagged, err = self.table("--max-ds=0.25", "--min-codons=3")
        self.assertTrue(any("dS above 0.25" in r["note"] for r in flagged))

    def test_bad_values_and_kaks_format(self):
        for args in (["--min-codons", "-1"], ["--min-codons", "x"], ["--min-codons", "1.5"], ["--max-ds", "0"],
                     ["--max-ds", "-2"], ["--max-ds", "nan"], ["--max-ds", "inf"], ["--max-ds"]):
            r = run_cpp("dnds", AXT, *args)
            self.assertEqual(r.returncode, 2, (args, r.stderr))
            self.assertRegex(r.stderr.decode(), "--min-codons|--max-ds")
        # KaKs_Calculator's format has no status column: the options belong to `coal dnds`
        r = run_cpp("kaks", "-i", AXT, "-o", self.tmp / "o.kaks", "--min-codons", "50")
        self.assertEqual(r.returncode, 2)

    def test_help_names_the_options(self):
        for sub in ("dnds", "codon-dnds"):
            text = run_cpp(sub, "--help").stdout.decode()
            self.assertIn("--min-codons N", text)
            self.assertIn("--max-ds X", text)


@needs_cpp
class PipelineFlags(PipelineCase):
    def run_pairs(self, name, **kw):
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_pairs(Q, S, self.tmp / name, pairs=self.pairs_file(PLAIN), jobs=1, **kw)
        return res, err.getvalue()

    def test_pairs_mode(self):
        plain, perr = self.run_pairs("plain")
        res, err = self.run_pairs("flag", min_codons=100, max_ds=0.5, method=["YN", "NG"])
        base, _ = self.run_pairs("base", method=["YN", "NG"])
        a, b = read_tsv(base.outputs["dnds"]), read_tsv(res.outputs["dnds"])
        changed = differing(a, b)
        self.assertGreater(len(changed), 0)
        self.assertLess(len(changed), len(a))
        for k in changed:
            self.assertEqual(b[k]["status"], "warning")
            self.assertRegex(b[k]["note"], "fewer than 100 codons compared|dS above 0.5")
        first = [r for r in b if r["method"] == b[0]["method"]]
        want = {"min_codons": sum("fewer than 100 codons compared" in r["note"] for r in first),
                "max_ds": sum("dS above 0.5" in r["note"] for r in first)}
        self.assertEqual(res.report["counts"]["flagged"], want)
        self.assertEqual(res.report["parameters"]["min_codons"], 100)
        self.assertEqual(res.report["parameters"]["max_ds"], 0.5)
        self.assertIn(f"summary: flagged (status warning, numbers unchanged): {want['min_codons']} pair(s) with "
                      f"fewer than 100 codons compared, {want['max_ds']} with dS above 0.5", err)
        self.check_report(res)
        self.assertEqual(res.exit_code, base.exit_code)   # a flag is not a failure
        # off: nothing in the report, the summary or the table
        self.assertNotIn("flagged", plain.report["counts"])
        self.assertNotIn("flagged", perr)
        self.assertEqual((plain.report["parameters"]["min_codons"], plain.report["parameters"]["max_ds"]), (0, None))

    def test_clusters_mode_summary_is_unchanged(self):
        with redirect_stderr(io.StringIO()):
            plain = P.run_clusters(FAM, self.tmp / "plain", clusters=FAM_CLUSTERS, max_cluster_size=2, jobs=1)
            res = P.run_clusters(FAM, self.tmp / "flag", clusters=FAM_CLUSTERS, max_cluster_size=2, jobs=1,
                                 min_codons=10000)
        self.assertEqual(Path(res.outputs["cluster_summary"]).read_bytes(),
                         Path(plain.outputs["cluster_summary"]).read_bytes())
        a, b = read_tsv(plain.outputs["dnds"]), read_tsv(res.outputs["dnds"])
        changed = differing(a, b)
        self.assertEqual(len(changed), sum(r["status"] in ("ok", "warning") for r in a))
        self.assertEqual(res.report["counts"]["flagged"], {"min_codons": len(changed)})
        self.assertEqual(res.report["counts"]["pairs"]["ok"], 0)

    def test_bad_values(self):
        for kw in ({"min_codons": -1}, {"min_codons": 2.5}, {"max_ds": 0}, {"max_ds": float("nan")},
                   {"max_ds": "2"}):
            with self.assertRaises(ValueError, msg=kw):
                P.run_pairs(Q, S, self.tmp / "x", **kw)
        for args in (["--min-codons", "-3"], ["--max-ds", "0"]):
            with redirect_stderr(io.StringIO()) as err, self.assertRaises(SystemExit) as cm:
                P.main(["pairs", "--query", str(Q), "--subject", str(S), "-o", str(self.tmp / "x"), *args])
            self.assertEqual(cm.exception.code, 2)
            self.assertIn(args[0], err.getvalue())


if __name__ == "__main__":
    unittest.main()
