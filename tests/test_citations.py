"""What to cite for what was run: one list (src/coal/citations.py), the same in every place that shows it.

`coal cite`, the `cite` list of the pipeline's report.json, README.md, docs/kaks.md, CITATION.cff and site/cite.html.
"""
from __future__ import annotations

import io
import json
import os
import re
import subprocess
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import citations as C  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402
from test_dnds_pipeline import FAM, FAM_CLUSTERS, PLAIN, Q, S, PipelineCase, needs_cpp, needs_diamond  # noqa: E402


def coal(*args):
    return subprocess.run([sys.executable, "-m", "coal", *map(str, args)], capture_output=True, text=True,
                          timeout=60, stdin=subprocess.DEVNULL, env={**os.environ, "PYTHONPATH": str(ROOT / "src")})


def title(key: str) -> str:
    """The title of a reference: between "(year). " and the next ". "."""
    return re.search(r"\(\d{4}\)\. (.+?)[.?] ", C.REFERENCES[key] + " ").group(1)


class References(unittest.TestCase):
    def test_every_key_is_complete(self):
        self.assertEqual(set(C.LABELS), set(C.REFERENCES))
        for key, text in C.REFERENCES.items():
            self.assertRegex(text, r"\(\d{4}\)\. ", key)
            self.assertNotIn("doi", text.lower(), key)   # none was verified: none is given
        for m, keys in C.METHOD_KEYS.items():
            self.assertIn(m, P.METHODS)
            for k in keys:
                self.assertIn(k, C.REFERENCES)
        self.assertEqual(set(C.METHOD_KEYS), set(P.METHODS))
        for _what, keys in C.TABLE:
            for k in keys:
                self.assertIn(k, C.REFERENCES)
        used = {k for _w, keys in C.TABLE for k in keys}
        self.assertEqual(used, set(C.REFERENCES), "a reference that no row of the table names")

    def test_keys_for(self):
        self.assertEqual(C.keys_for(), ["coal", "suyama2006"])
        self.assertEqual(C.keys_for(methods=["YN"], diamond=True),
                         ["coal", "suyama2006", "yang2000", "wang2010", "yang2007", "buchfink2021"])
        got = C.keys_for(methods=["MYN", "YN", "Li", "Comeron", "MLPB"], deepclust=True, clustalo=True, mafft=True)
        self.assertEqual(len(got), len(set(got)))
        for k in ("zhang2006", "comeron1995", "thornton2003", "charif2007", "tzeng2004", "pamilo1993",
                  "buchfink2021", "buchfink2023", "sievers2011", "katoh2013"):
            self.assertIn(k, got)
        self.assertNotIn("nei1986", got)

    def test_documents_hold_the_table_and_every_reference(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        kaks = (ROOT / "docs" / "kaks.md").read_text(encoding="utf-8")
        cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        site = (ROOT / "site" / "cite.html").read_text(encoding="utf-8")
        self.assertIn(C.markdown_table(), readme)
        self.assertIn(C.markdown_table(), kaks)
        self.assertIn(C.markdown_references(), kaks)
        plain_site = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", site)).replace("&amp;", "&")
        for key in C.REFERENCES:
            if key == "coal":
                continue
            t = title(key)
            self.assertIn(t, kaks, key)
            self.assertIn(t, cff, f"CITATION.cff lacks {key}")
            self.assertIn(f"# {key}:", cff)
            self.assertIn(t, plain_site, f"site/cite.html lacks {key}")
        for what, keys in C.TABLE:
            for k in keys:
                self.assertIn(C.LABELS[k], plain_site, f"site/cite.html: the row {what!r} lacks {k}")

    def test_citation_cff_parses(self):
        try:
            import yaml
        except ImportError:
            self.skipTest("PyYAML is not installed")
        d = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
        self.assertEqual(len(d["references"]), len(C.REFERENCES) - 1)
        for ref in d["references"]:
            self.assertIn(ref["type"], ("article", "generic"))
            self.assertTrue(ref["authors"] and ref["title"] and ref["year"])
            self.assertNotIn("doi", ref)


class CiteCommand(unittest.TestCase):
    def test_table_and_all_references(self):
        r = coal("cite")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("What you ran -> what to cite", r.stdout)
        for key, text in C.REFERENCES.items():
            self.assertIn(f"[{key}] " + text.split("{version}")[0], r.stdout)

    def test_methods_keys_and_errors(self):
        r = coal("cite", "yn")
        self.assertEqual([ln.split("]")[0][1:] for ln in r.stdout.splitlines()],
                         ["coal", "suyama2006", "yang2000", "wang2010", "yang2007"])
        self.assertIn(f"version {__import__('coal').__version__}", r.stdout)
        r = coal("cite", "katoh2013", "C95", "katoh2013")
        self.assertEqual([ln.split("]")[0][1:] for ln in r.stdout.splitlines()],
                         ["katoh2013", "coal", "suyama2006", "comeron1995", "thornton2003"])
        r = coal("cite", "nobody2020")
        self.assertEqual(r.returncode, 2)
        self.assertIn("is not a citation key", r.stderr)
        self.assertEqual(coal("cite", "--help").returncode, 0)
        r = coal("cite", ROOT / "pyproject.toml")
        self.assertEqual(r.returncode, 1)
        self.assertIn("is not a `coal pipeline` report", r.stderr)


@needs_cpp
class ReportCite(PipelineCase):
    def test_pairs_given_and_cite_from_the_report(self):
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_pairs(Q, S, self.tmp / "o", pairs=self.pairs_file(PLAIN), method=["NG", "YN"])
        self.assertEqual(res.report["cite"], ["coal", "suyama2006", "nei1986", "wang2010", "yang2000", "yang2007"])
        self.assertIn("summary: what to cite for this run: coal cite ", err.getvalue())
        r = coal("cite", res.outputs["report"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual([ln.split("]")[0][1:] for ln in r.stdout.splitlines()], res.report["cite"])
        on_disk = json.loads(Path(res.outputs["report"]).read_text(encoding="utf-8"))
        self.assertEqual(on_disk["cite"], res.report["cite"])

    def test_clusters_given_with_fallback_aligner(self):
        with redirect_stderr(io.StringIO()):
            res = P.run_clusters(FAM, self.tmp / "o", clusters=FAM_CLUSTERS, max_cluster_size=2, method=["Li"])
        self.assertEqual(res.report["cite"], ["coal", "suyama2006", "li1993", "charif2007"])


@needs_diamond
class ReportCiteDiamond(PipelineCase):
    def test_search_and_deepclust(self):
        with redirect_stderr(io.StringIO()):
            res = P.run_pairs(Q, S, self.tmp / "o", jobs=2)
            cl = P.run_clusters(FAM, self.tmp / "c", max_cluster_size=2, jobs=2)
        self.assertEqual(res.report["cite"][-1], "buchfink2021")
        self.assertEqual(cl.report["cite"][-2:], ["buchfink2021", "buchfink2023"])


if __name__ == "__main__":
    unittest.main()
