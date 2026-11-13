"""`coal pipeline clusters`: the clusters do not depend on the order of the input (docs/pipeline.md
"Reproducibility").

DIAMOND DeepClust's result depends on the order of the sequences it is given, so the pipeline gives it the proteins
sorted by ID. Units: the sorter. With coal-cpp and a stand-in DIAMOND (a shell script that records what it was
given): the bytes DeepClust receives are the same for shuffled records and for the same records split over several
--cds files. With the real DIAMOND: clusters.tsv, dnds.tsv and cluster_summary.tsv are byte-identical.
"""
from __future__ import annotations

import os
import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import dnds_pipeline as P  # noqa: E402
from test_dnds_pipeline import FAM, PipelineCase, needs_cpp, needs_diamond  # noqa: E402


def records(path) -> list[str]:
    return [">" + r for r in Path(path).read_text(encoding="utf-8").split(">")[1:]]


def shuffled_inputs(tmp: Path, seed: int) -> list[Path]:
    """The family fixture's records (first record of an ID) in another order, split over two files."""
    seen, recs = set(), []
    for r in records(FAM):
        rid = r[1:].split()[0] if r[1:].split() else ""
        if rid not in seen:
            seen.add(rid)
            recs.append(r if r.endswith("\n") else r + "\n")
    random.Random(seed).shuffle(recs)
    a, b = tmp / f"part_a_{seed}.fasta", tmp / f"part_b_{seed}.fasta"
    a.write_text("".join(recs[: len(recs) // 3]), encoding="utf-8")
    b.write_text("".join(recs[len(recs) // 3:]), encoding="utf-8")
    return [b, a]


class SortFasta(unittest.TestCase):
    def test_sorted_by_id_bytes_and_complete(self):
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            src, dst = Path(t) / "in.faa", Path(t) / "out.faa"
            src.write_bytes(b">b2 desc\nMKV\n>B1\nMK\nLL\n>a10\nM\n>a9\nMM\n>\xe9x\nMW\n>a9\nSECOND\n>z")
            self.assertEqual(P.sort_fasta_by_id(str(src), str(dst)), 7)
            self.assertEqual(dst.read_bytes(), b">B1\nMK\nLL\n>a10\nM\n>a9\nMM\n>a9\nSECOND\n>b2 desc\nMKV\n>z\n"
                                               b">\xe9x\nMW\n")
            self.assertFalse((Path(t) / "out.faa.tmp").exists())
            src.write_bytes(b"")
            self.assertEqual(P.sort_fasta_by_id(str(src), str(dst)), 0)
            self.assertEqual(dst.read_bytes(), b"")


@needs_cpp
@unittest.skipIf(os.name == "nt", "the stand-in DIAMOND is a POSIX shell script")
class DeepClustInput(PipelineCase):
    def fake_diamond(self) -> Path:
        """Records the --db file of `deepclust` and writes one cluster per sequence."""
        fake = self.tmp / "diamond"
        fake.write_text(
            "#!/bin/sh\n"
            'if [ "$1" = version ]; then echo "diamond version 2.2.8"; exit 0; fi\n'
            'db=; out=\n'
            'while [ $# -gt 0 ]; do case "$1" in --db) db=$2; shift;; --out) out=$2; shift;; esac; shift; done\n'
            f'cp "$db" "{self.tmp}/given_$(basename "$(dirname "$out")").faa"\n'
            "grep '^>' \"$db\" | sed 's/^>//; s/[[:space:]].*//' | awk '{print $1 \"\\t\" $1}' > \"$out\"\n",
            encoding="utf-8")
        fake.chmod(0o755)
        return fake

    def test_deepclust_gets_the_same_bytes_for_any_input_order(self):
        fake = self.fake_diamond()
        given = []
        for name, cds in (("o1", [FAM]), ("o2", shuffled_inputs(self.tmp, 1)), ("o3", shuffled_inputs(self.tmp, 2))):
            res = P.run_clusters(cds, self.tmp / name, diamond=str(fake), quiet=True, jobs=1)
            self.assertIn(res.exit_code, (0, 2), res.report.get("error"))
            given.append((self.tmp / f"given_{name}.faa").read_bytes())
            stage = next(s for s in res.report["stages"] if s["name"] == "deepclust")
            self.assertEqual(stage["input_order"], "sorted by ID")
        self.assertEqual(given[0], given[1])
        self.assertEqual(given[0], given[2])
        ids = [ln[1:].split()[0] for ln in given[0].decode().splitlines() if ln.startswith(">")]
        self.assertEqual(ids, sorted(ids, key=lambda x: x.encode()))
        self.assertGreater(len(ids), 15)
        # the genes table still follows the input (only DeepClust's input is sorted)
        g1 = [ln.split("\t")[1] for ln in (self.tmp / "o1" / "genes.tsv").read_text().splitlines()[1:]]
        g2 = [ln.split("\t")[1] for ln in (self.tmp / "o2" / "genes.tsv").read_text().splitlines()[1:]]
        self.assertNotEqual(g1, g2)
        self.assertEqual(sorted(set(g1)), sorted(set(g2)))


@needs_diamond
class DeepClustInvariance(PipelineCase):
    def test_clusters_and_tables_do_not_depend_on_file_or_record_order(self):
        runs = []
        for name, cds in (("o1", [FAM]), ("o2", shuffled_inputs(self.tmp, 3)), ("o3", shuffled_inputs(self.tmp, 4))):
            res = P.run_clusters(cds, self.tmp / name, max_cluster_size=2, quiet=True, jobs=2)
            self.assertIn(res.exit_code, (0, 2), res.report.get("error"))
            runs.append({k: Path(res.outputs[k]).read_bytes() for k in ("clusters", "dnds", "cluster_summary")})
        for k in ("clusters", "dnds", "cluster_summary"):
            self.assertEqual(runs[0][k], runs[1][k], k)
            self.assertEqual(runs[0][k], runs[2][k], k)
        self.assertGreater(runs[0]["dnds"].count(b"\n"), 10)


if __name__ == "__main__":
    unittest.main()
