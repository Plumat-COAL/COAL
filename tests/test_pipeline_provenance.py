"""report.json of `coal pipeline`: what a methods section needs (docs/pipeline.md "report.json").

Start and end time (UTC), platform, Python and COAL versions, the thread counts as resolved, every input file with its
size and SHA-256 (skipped above a size limit, with the reason; --hash-inputs always / never), the exact command lines
of DIAMOND and coal-cpp and the template of the aligner's, the citation keys. The tables are not touched: they are
byte-identical whatever --hash-inputs is.
"""
from __future__ import annotations

import hashlib
import io
import os
import re
import sys
import time
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
import coal  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402
from test_dnds_pipeline import (CLUSTALO, CPP, FAM, FAM_CLUSTERS, PLAIN, Q, S, PipelineCase,  # noqa: E402
                                needs_clustalo, needs_cpp, needs_diamond)

ISO = re.compile(r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$")


def quiet(fn, *a, **k):
    with redirect_stderr(io.StringIO()):
        return fn(*a, **k)


@needs_cpp
class Provenance(PipelineCase):
    def test_times_platform_threads_inputs(self):
        pairs = self.pairs_file(PLAIN)
        t0 = time.time()
        res = quiet(P.run_pairs, Q, S, self.tmp / "o", pairs=pairs, jobs=1)
        rep = res.report
        self.check_report(res)
        for k in ("started", "finished"):
            self.assertRegex(rep[k], ISO)
            got = time.mktime(time.strptime(rep[k], "%Y-%m-%dT%H:%M:%SZ")) - time.timezone
            self.assertLess(abs(got - t0), 120, k)   # UTC, now
        self.assertLessEqual(rep["started"], rep["finished"])
        import platform
        self.assertEqual(rep["platform"]["system"], platform.system())
        self.assertEqual(rep["platform"]["python"], platform.python_version())
        self.assertEqual(rep["platform"]["machine"], platform.machine())
        self.assertEqual(rep["coal_version"], coal.__version__)
        self.assertEqual(rep["threads"], {"requested": 1, "engine": 1, "tools": 1})
        self.assertEqual([i["option"] for i in rep["inputs"]], ["--query", "--subject", "--pairs"])
        for rec, path in zip(rep["inputs"], (Q, S, pairs)):
            self.assertEqual(rec["path"], os.path.abspath(path))
            self.assertEqual(rec["bytes"], path.stat().st_size)
            self.assertEqual(rec["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertRegex(rec["modified"], ISO)
            self.assertNotIn("sha256_skipped", rec)
        # the engine's command line, as it was run
        cmds = rep["tools"]["coal-cpp"]["commands"]
        self.assertEqual([c[1] for c in cmds], ["translate", "translate", "codon-dnds"])
        self.assertTrue(all(c[0] == str(CPP) for c in cmds))
        self.assertIn("-m", cmds[2])
        stage = {s["name"]: s for s in rep["stages"]}
        self.assertEqual(stage["translate_query"]["commands"], [cmds[0]])

    def test_default_jobs_are_resolved(self):
        res = quiet(P.run_pairs, Q, S, self.tmp / "o", pairs=self.pairs_file(PLAIN))
        cpus = os.cpu_count() or 1
        self.assertEqual(res.report["threads"], {"requested": 0, "engine": cpus, "tools": cpus})
        self.assertEqual(res.report["platform"]["cpus"], os.cpu_count())

    def test_hash_limit_and_option(self):
        pairs = self.pairs_file(PLAIN)
        base = quiet(P.run_pairs, Q, S, self.tmp / "base", pairs=pairs)
        with mock.patch.object(P, "HASH_LIMIT", Q.stat().st_size - 1):   # the query is "too large", the rest not
            auto = quiet(P.run_pairs, Q, S, self.tmp / "auto", pairs=pairs)
            forced = quiet(P.run_pairs, Q, S, self.tmp / "forced", pairs=pairs, hash_inputs="always")
        q = auto.report["inputs"][0]
        self.assertIsNone(q["sha256"])
        self.assertIn("--hash-inputs always", q["sha256_skipped"])
        self.assertIn("larger than", q["sha256_skipped"])
        self.assertEqual(q["bytes"], Q.stat().st_size)
        self.assertEqual(auto.report["inputs"][1]["sha256"], base.report["inputs"][1]["sha256"])
        self.assertEqual([i["sha256"] for i in forced.report["inputs"]], [i["sha256"] for i in base.report["inputs"]])
        never = quiet(P.run_pairs, Q, S, self.tmp / "never", pairs=pairs, hash_inputs="never")
        self.assertEqual([(i["sha256"], i["sha256_skipped"]) for i in never.report["inputs"]],
                         [(None, "--hash-inputs never")] * 3)
        for r in (auto, forced, never):
            self.check_report(r)
            self.assertEqual(Path(r.outputs["dnds"]).read_bytes(), Path(base.outputs["dnds"]).read_bytes())
            self.assertEqual(Path(r.outputs["genes"]).read_bytes(), Path(base.outputs["genes"]).read_bytes())
        with self.assertRaises(ValueError):
            P.run_pairs(Q, S, self.tmp / "x", hash_inputs="sometimes")

    def test_a_failed_run_still_says_when_and_where(self):
        res = quiet(P.run_pairs, Q, S, self.tmp / "o", pairs=self.tmp / "no_such_pairs.tsv")
        self.assertEqual(res.exit_code, 1)
        rep = res.report
        self.assertRegex(rep["started"], ISO)
        self.assertEqual(rep["inputs"][2]["option"], "--pairs")
        self.assertIn("error", rep["inputs"][2])
        self.assertEqual(rep["cite"][:2], ["coal", "suyama2006"])

    def test_gz_input_is_described_as_given(self):
        import gzip
        qz = self.tmp / "q.fa.gz"
        qz.write_bytes(gzip.compress(Q.read_bytes()))
        res = quiet(P.run_pairs, qz, S, self.tmp / "o", pairs=self.pairs_file(PLAIN))
        rec = res.report["inputs"][0]
        self.assertEqual((rec["compression"], rec["bytes"]), ("gzip", qz.stat().st_size))
        self.assertEqual(rec["sha256"], hashlib.sha256(qz.read_bytes()).hexdigest())
        self.assertNotIn("compression", res.report["inputs"][1])


@needs_diamond
class DiamondCommands(PipelineCase):
    def test_commands_are_recorded_and_survive_resume(self):
        res = quiet(P.run_pairs, Q, S, self.tmp / "o", jobs=2)
        cmds = res.report["tools"]["diamond"]["commands"]
        self.assertEqual([c[1] for c in cmds], ["makedb", "blastp", "makedb", "blastp"])
        self.assertTrue(all(c[0] == res.report["tools"]["diamond"]["path"] for c in cmds))
        self.assertIn("--max-target-seqs", cmds[1])
        log = (self.tmp / "o" / "work" / "logs" / "diamond.log").read_text(encoding="utf-8")
        for c in cmds:   # the report has what the log shows was run
            self.assertIn("$ " + " ".join(c), log)
        again = quiet(P.run_pairs, Q, S, self.tmp / "o", jobs=2, resume=True)
        self.assertTrue(next(s for s in again.report["stages"] if s["name"] == "diamond_search")["cached"])
        self.assertEqual(again.report["tools"]["diamond"]["commands"], cmds)
        self.check_report(again)


@needs_clustalo
class AlignerTemplate(PipelineCase):
    def test_cluster_run_has_the_aligner_command_template(self):
        res = quiet(P.run_clusters, FAM, self.tmp / "o", clusters=FAM_CLUSTERS, jobs=2, msa_threads=1)
        tool = res.report["tools"]["clustalo"]
        self.assertEqual(tool["command_template"][0], str(CLUSTALO))
        self.assertEqual(tool["command_template"][1:5], ["-i", "IN.fasta", "-o", "OUT.fasta"])
        self.assertIn("--threads=1", tool["command_template"])
        self.assertIn("sievers2011", res.report["cite"])
        self.check_report(res)


if __name__ == "__main__":
    unittest.main()
