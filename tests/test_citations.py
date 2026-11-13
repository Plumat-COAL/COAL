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
            self.assertNotIn("doi", text.lower(), key)   # the DOI is in DOIS, printed as a link after the text
        # every cited work has a DOI; COAL itself has none (the owner's rule for the website and CITATION.cff)
        self.assertEqual(set(C.DOIS), set(C.REFERENCES) - {"coal"})
        for key, doi in C.DOIS.items():
            self.assertRegex(doi, r"^10\.\d{4,5}/\S+$", key)
            self.assertEqual(C.doi_url(key), "https://doi.org/" + doi)
        self.assertEqual(len(set(C.DOIS.values())), len(C.DOIS))
        self.assertIsNone(C.doi_url("coal"))
        for m, keys in C.METHOD_KEYS.items():
            self.assertIn(m, P.METHODS)
            for k in keys:
                self.assertIn(k, C.REFERENCES)
        self.assertEqual(set(C.METHOD_KEYS), set(P.METHODS))
        for _what, keys in C.TABLE:
            for k in keys:
                self.assertIn(k, C.REFERENCES)
        used = {k for _w, keys in C.TABLE for k in keys}
        self.assertEqual(used, set(C.REFERENCES) - set(C.SUPERSEDED), "a reference that no row of the table names")

    def test_deepclust_paper_replaces_its_preprint(self):
        """buchfink2026 (Nature Methods) is the DeepClust citation; buchfink2023 (bioRxiv) only resolves."""
        self.assertEqual(C.SUPERSEDED, {"buchfink2023": "buchfink2026"})
        self.assertIn("Nature Methods 23:724-727", C.REFERENCES["buchfink2026"])
        self.assertIn("Barb\u00e9 \u00c9", C.REFERENCES["buchfink2026"])
        self.assertEqual(C.DOIS["buchfink2026"], "10.1038/s41592-026-03030-z")
        self.assertEqual(C.DOIS["buchfink2023"], "10.1101/2023.01.24.525373")
        self.assertNotIn("buchfink2023", C.listed_keys())
        self.assertNotIn("buchfink2023", C.keys_for(deepclust=True))
        self.assertEqual(C.current(["coal", "buchfink2023", "buchfink2026"]), ["coal", "buchfink2023", "buchfink2026"])
        for doc in ("README.md", "docs/kaks.md", "docs/pipeline.md", "CITATION.cff", "site/cite.html"):
            text = (ROOT / doc).read_text(encoding="utf-8")
            self.assertNotIn("DeepClust, preprint", text, doc)
            self.assertNotIn("Sensitive clustering of protein sequences", text, doc)

    def test_keys_for(self):
        self.assertEqual(C.keys_for(), ["coal", "suyama2006"])
        self.assertEqual(C.keys_for(methods=["YN"], diamond=True),
                         ["coal", "suyama2006", "yang2000", "wang2010", "yang2007", "buchfink2021"])
        got = C.keys_for(methods=["MYN", "YN", "Li", "Comeron", "MLPB"], deepclust=True, clustalo=True, mafft=True)
        self.assertEqual(len(got), len(set(got)))
        for k in ("zhang2006", "comeron1995", "thornton2003", "charif2007", "tzeng2004", "pamilo1993",
                  "buchfink2021", "buchfink2026", "sievers2011", "katoh2013"):
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
        for key in C.listed_keys():
            t = title(key)
            self.assertIn(t, kaks, key)
            self.assertIn(t, cff, f"CITATION.cff lacks {key}")
            self.assertIn(f"# {key}:", cff)
            self.assertIn(t, plain_site, f"site/cite.html lacks {key}")
            # the DOI of the work in every copy: a link on the website and in the documents, `doi:` in CITATION.cff
            self.assertIn(f'<a href="{C.doi_url(key)}">', site, key)
            self.assertIn(f"<{C.doi_url(key)}>", kaks, key)
            self.assertIn(f'doi: "{C.DOIS[key]}"', cff, key)
        self.assertEqual(len(re.findall(r"https://doi\.org/", site)), len(C.listed_keys()))
        self.assertIn("<https://doi.org/10.1093/nar/gkl315>", readme)
        for what, keys in C.TABLE:
            for k in keys:
                self.assertIn(C.LABELS[k], plain_site, f"site/cite.html: the row {what!r} lacks {k}")

    def test_citation_cff_parses(self):
        try:
            import yaml
        except ImportError:
            self.skipTest("PyYAML is not installed")
        d = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
        self.assertEqual(len(d["references"]), len(C.listed_keys()))
        self.assertNotIn("doi", d)             # COAL itself has no DOI; these are the DOIs of the cited works
        self.assertNotIn("identifiers", d)
        dois = set()
        for ref in d["references"]:
            self.assertIn(ref["type"], ("article", "generic"))
            self.assertTrue(ref["authors"] and ref["title"] and ref["year"])
            dois.add(ref["doi"])
            by_title = [k for k in C.listed_keys() if title(k) == ref["title"]]
            self.assertEqual(len(by_title), 1, ref["title"])
            self.assertEqual(ref["doi"], C.DOIS[by_title[0]])
            self.assertIn(f"({ref['year']})", C.REFERENCES[by_title[0]])
            self.assertIn(ref["authors"][0]["family-names"], C.REFERENCES[by_title[0]])
        self.assertEqual(dois, {C.DOIS[k] for k in C.listed_keys()})


class CiteCommand(unittest.TestCase):
    def test_table_and_all_references(self):
        r = coal("cite")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("What you ran -> what to cite", r.stdout)
        for key in ["coal", *C.listed_keys()]:
            self.assertIn(f"[{key}] " + C.REFERENCES[key].split("{version}")[0], r.stdout)
            if key != "coal":
                self.assertIn(f" https://doi.org/{C.DOIS[key]}\n", r.stdout)
        self.assertNotIn("[buchfink2023]", r.stdout)
        self.assertNotIn("doi.org", r.stdout.split("[suyama2006]")[0].split("[coal]")[1])   # none for COAL

    def test_methods_keys_and_errors(self):
        r = coal("cite", "yn")
        self.assertEqual([ln.split("]")[0][1:] for ln in r.stdout.splitlines()],
                         ["coal", "suyama2006", "yang2000", "wang2010", "yang2007"])
        self.assertIn(f"version {__import__('coal').__version__}", r.stdout)
        r = coal("cite", "katoh2013", "C95", "katoh2013")
        self.assertEqual([ln.split("]")[0][1:] for ln in r.stdout.splitlines()],
                         ["katoh2013", "coal", "suyama2006", "comeron1995", "thornton2003"])
        r = coal("cite", "buchfink2023")   # the key of the preprint still resolves, and names the paper
        self.assertEqual([ln.split("]")[0][1:] for ln in r.stdout.splitlines()], ["buchfink2023", "buchfink2026"])
        self.assertIn("https://doi.org/10.1101/2023.01.24.525373", r.stdout)
        r = coal("cite", "nobody2020")
        self.assertEqual(r.returncode, 2)
        self.assertIn("is not a citation key", r.stderr)
        self.assertEqual(coal("cite", "--help").returncode, 0)
        r = coal("cite", ROOT / "pyproject.toml")
        self.assertEqual(r.returncode, 1)
        self.assertIn("is not a `coal pipeline` report", r.stderr)


class OldReport(unittest.TestCase):
    def test_a_report_that_names_the_preprint_still_resolves(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            rep = Path(d) / "report.json"
            rep.write_text(json.dumps({"cite": ["coal", "suyama2006", "buchfink2021", "buchfink2023"]}),
                           encoding="utf-8", newline="\n")
            r = coal("cite", rep)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual([ln.split("]")[0][1:] for ln in r.stdout.splitlines()],
                         ["coal", "suyama2006", "buchfink2021", "buchfink2023", "buchfink2026"])

    def test_ascii_only_terminal(self):
        """An author's name with an accent does not end `coal cite` with a traceback."""
        r = subprocess.run([sys.executable, "-m", "coal", "cite", "buchfink2026"], capture_output=True, timeout=60,
                           stdin=subprocess.DEVNULL,
                           env={**os.environ, "PYTHONPATH": str(ROOT / "src"), "PYTHONIOENCODING": "ascii",
                                "PYTHONUTF8": "0"})
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn(b"Barb? ?, Ashkenazy", r.stdout)
        self.assertIn(b"https://doi.org/10.1038/s41592-026-03030-z", r.stdout)


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
        self.assertEqual(cl.report["cite"][-2:], ["buchfink2021", "buchfink2026"])


if __name__ == "__main__":
    unittest.main()
