"""Several transcripts per gene (docs/pipeline.md "One record is one gene"): the warning, --gene-map and
--longest-isoform.

Fixtures: tests/fixtures/isoforms (tests/tools/make_isoform_fixture.py; synthetic sequences): the same six ortholog
pairs as 11 transcripts per file, with Ensembl-style headers (gene:), NCBI-style headers ([gene=]) and TAIR-style
IDs (AT1G01010.1 / .2, no gene field); expected.tsv names the record that is the longest of each gene.

  * Units: the header fields, the suffix rules and when they do not apply (transcript versions, numbered genes), the
    gene map reader, the choice of the longest record and its tie rule.
  * With coal-cpp: the default run is unchanged but for one warning; --longest-isoform keeps the expected records in
    both modes, lists the others in genes.tsv with status `dropped`, and the report records the choice.
  * With DIAMOND: with isoforms the reciprocal best hits name transcripts; with --longest-isoform there is one pair
    per gene, between the longest records.
"""
from __future__ import annotations

import io
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import dnds_pipeline as P  # noqa: E402
from coal import isoforms as I  # noqa: E402
from test_dnds_pipeline import PLAIN, Q, S, PipelineCase, needs_cpp, needs_diamond, read_tsv  # noqa: E402

FIX = ROOT / "tests" / "fixtures" / "isoforms"
ENS, NCBI, TAIR, GMAP = (FIX / n for n in ("ensembl_query_cds.fasta", "ncbi_subject_cds.fasta",
                                           "tair_query_cds.fasta", "gene_map.tsv"))
EXPECTED = {(r["file"], r["gene"]): r["kept"] for r in read_tsv(FIX / "expected.tsv")}
KEPT = {f: sorted(k for (ff, _), k in EXPECTED.items() if ff == f) for f in {f for f, _ in EXPECTED}}


def run(fn, *args, **kw):
    with redirect_stderr(io.StringIO()) as err:
        res = fn(*args, **kw)
    return res, err.getvalue()


class Detect(unittest.TestCase):
    def det(self, lines):
        return I.detect([(ln.split(None, 1) + [""])[:2] for ln in lines])

    def test_fixture_files(self):
        for path, source in ((ENS, "header field gene:"), (NCBI, "header field [gene="), (TAIR, "ID suffix .N")):
            d = I.detect(I.read_headers([path]))
            self.assertEqual((d.source, d.records, d.genes, d.multi, d.extra), (source, 11, 6, 4, 5), path.name)
            self.assertEqual(len(d.example[1]), 3)
        self.assertEqual(I.detect(I.read_headers([Q])).source, None)   # plain IDs: nothing to find

    def test_header_fields(self):
        d = self.det(["t1 gene:G1 x", "t2 gene:G1", "t3 gene:G2"])
        self.assertEqual((d.source, d.gene_of), ("header field gene:", {"t1": "G1", "t2": "G1", "t3": "G2"}))
        d = self.det(["a [locus_tag=B0001] [protein=x y]", "b [locus_tag=B0001]", "c [locus_tag=B0002]"])
        self.assertEqual((d.source, d.multi), ("header field [locus_tag=", 1))
        d = self.det(["a ID=x;gene=AT1G01010;len=3", "b gene=AT1G01010", "c gene=AT1G01020"])
        self.assertEqual((d.source, d.gene_of["a"], d.multi), ("header field gene=", "AT1G01010", 1))
        d = self.det(["a gene_symbol:TP53", "b gene_symbol:TP53"])
        self.assertEqual(d.source, "header field gene_symbol:")
        # the gene field wins over the symbol, and a field most records lack is not used
        d = self.det(["a gene:G1 gene_symbol:S", "b gene:G2 gene_symbol:S"])
        self.assertEqual((d.source, d.multi), ("header field gene:", 0))
        d = self.det(["a gene:G1"] + [f"r{k}" for k in range(20)])
        self.assertIsNone(d.source)
        # words that only contain the field name are not the field
        d = self.det(["a pseudogene:yes transgene=1", "b pseudogene:yes transgene=1"])
        self.assertIsNone(d.source)

    def test_suffix_rules_are_conservative(self):
        d = self.det(["AT1G01010.1", "AT1G01010.2", "AT1G01020.1", "AT1G01030.1"])
        self.assertEqual((d.source, d.genes, d.multi, d.extra), ("ID suffix .N", 3, 1, 1))
        self.assertEqual(d.gene_of, {"AT1G01010.1": "AT1G01010", "AT1G01010.2": "AT1G01010"})
        # transcript versions (one record per stem): not isoforms
        self.assertIsNone(self.det(["ENST00000380152.8", "ENST00000544455.6", "NM_000059.4"]).source)
        # one record per gene, all .1 (a representative gene model set): nothing
        self.assertIsNone(self.det([f"AT1G{k:05d}.1" for k in range(50)]).source)
        # genes numbered within a scaffold: many records per stem, not isoforms
        self.assertIsNone(self.det([f"scaffold{s}.{k}" for s in (1, 2) for k in range(1, 40)]).source)
        d = self.det(["g1.t1", "g1.t2", "g2.t1"])
        self.assertEqual((d.source, d.multi), ("ID suffix .tN", 1))
        d = self.det(["TRINITY_DN1_c0_g1_i1", "TRINITY_DN1_c0_g1_i2", "TRINITY_DN1_c0_g2_i1"])
        self.assertEqual((d.source, d.gene_of["TRINITY_DN1_c0_g1_i2"]), ("ID suffix _iN (Trinity)", "TRINITY_DN1_c0_g1"))
        self.assertIsNone(self.det([]).source)

    def test_gene_map(self):
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "m.tsv"
            p.write_text("# comment\ntranscript_id\tgene_id\nt1\tG1\nt2\tG1\tignored\n\nt3 G2\nt1\tOTHER\n", encoding="utf-8")
            m = I.read_gene_map(str(p))
            self.assertEqual(m, {"t1": "G1", "t2": "G1", "t3": "G2"})
            d = I.from_gene_map(m, ["t1", "t2", "t3", "t4"])
            self.assertEqual((d.source, d.genes, d.multi, d.extra), ("gene map", 3, 1, 1))
            p.write_text("t1\tG1\nonly_one_column\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "line 2: expected two columns"):
                I.read_gene_map(str(p))

    def test_choose_longest(self):
        d = I.Detection(gene_of={"a": "G", "b": "G", "c": "G", "x": "H", "y": "H", "z": "K"})
        got = I.choose_longest(d, {"a": 300, "b": 900, "c": 900, "x": 30, "y": 30, "z": 3, "free": 5})
        self.assertEqual(got, {"a": ("G", "b"), "c": ("G", "b"), "y": ("H", "x")})
        # a record that could not be used is no candidate: the longest usable one stays
        self.assertEqual(I.choose_longest(d, {"a": 300, "c": 900}), {"a": ("G", "c")})


@needs_cpp
class Pipeline(PipelineCase):
    def pairs(self, qf, sf):
        """One pair per gene between the kept records, and one that names a dropped record."""
        q = [EXPECTED[(qf.name, g)] for (f, g) in EXPECTED if f == qf.name]
        s = [EXPECTED[(sf.name, g)] for (f, g) in EXPECTED if f == sf.name]
        return list(zip(q, s))

    def test_default_is_unchanged_but_for_one_warning(self):
        pairs = self.pairs(ENS, NCBI)
        res, err = run(P.run_pairs, ENS, NCBI, self.tmp / "o", pairs=self.pairs_file(pairs))
        self.assertEqual(res.exit_code, 0, err)
        lines = [ln for ln in err.splitlines() if "warning:" in ln]
        self.assertEqual(len(lines), 1, err)
        self.assertIn("the input looks like several transcripts per gene", lines[0])
        self.assertIn("query: 11 records of 6 genes, 4 genes with several records", lines[0])
        self.assertIn("read from the header field gene:", lines[0])
        self.assertIn("subject: 11 records of 6 genes", lines[0])
        self.assertIn("header field [gene=", lines[0])
        self.assertIn("--longest-isoform", lines[0])
        self.assertIn("--gene-map", lines[0])
        genes = read_tsv(res.outputs["genes"])
        self.assertEqual(len(genes), 22)
        self.assertEqual({g["status"] for g in genes}, {"ok"})        # nothing is dropped by default
        self.assertEqual(res.report["isoforms"]["longest_isoform"], False)
        self.assertEqual(res.report["isoforms"]["sets"]["query"]["genes_with_several_records"], 4)
        self.assertEqual(res.report["counts"]["genes"], {"total": 22, "ok": 22, "warning": 0, "failed": 0})
        self.check_report(res)
        # and input without transcripts gets no warning and the same tables as ever
        res, err = run(P.run_pairs, Q, S, self.tmp / "plain", pairs=self.pairs_file(PLAIN))
        self.assertNotIn("transcripts per gene", err)
        self.assertNotIn("warnings", res.report)
        self.assertEqual(res.report["isoforms"]["sets"]["query"]["gene_source"], None)

    def test_longest_isoform_from_headers(self):
        pairs = self.pairs(ENS, NCBI)
        dropped_q = sorted(set(r for r, _ in I.read_headers([ENS])) - set(KEPT[ENS.name]))[0]
        res, err = run(P.run_pairs, ENS, NCBI, self.tmp / "o", longest_isoform=True,
                       pairs=self.pairs_file(pairs + [(dropped_q, pairs[0][1])]))
        self.assertIn("--longest-isoform: query: 5 of 11 records dropped, 6 genes (from header field gene:); "
                      "subject: 5 of 11 records dropped, 6 genes (from header field [gene=); 10 in all, listed in "
                      "genes.tsv with status dropped", err)
        self.assertNotIn("looks like several transcripts", err)
        genes = read_tsv(res.outputs["genes"])
        for name, set_name in ((ENS.name, "query"), (NCBI.name, "subject")):
            mine = [g for g in genes if g["set"] == set_name]
            self.assertEqual(sorted(g["id"] for g in mine if g["status"] == "ok"), KEPT[name])
            drops = [g for g in mine if g["status"] == "dropped"]
            self.assertEqual(len(drops), 5)
            for g in drops:
                self.assertEqual(g["code"], "not_longest_isoform")
                self.assertRegex(g["note"], r"gene \S+: \S+ is kept \(--longest-isoform\)")
                self.assertIn(g["note"].split(": ")[1].split(" is kept")[0], KEPT[name])
        self.assertEqual(res.report["counts"]["genes"], {"total": 22, "ok": 12, "warning": 0, "failed": 0,
                                                         "dropped": 10})
        rows = read_tsv(res.outputs["dnds"])
        self.assertEqual([(r["seq1"], r["seq2"]) for r in rows[:6]], pairs)
        self.assertTrue(all(r["status"] == "ok" for r in rows[:6]))
        # the pair that names a dropped record says so, and is the run's one failure
        self.assertEqual(rows[6]["status"], "error")
        self.assertIn("dropped by --longest-isoform", rows[6]["note"])
        self.assertEqual(res.exit_code, 2)
        iso = res.report["isoforms"]
        self.assertEqual((iso["longest_isoform"], iso["gene_map"]), (True, None))
        self.assertEqual(iso["sets"]["subject"]["dropped"], 5)
        self.check_report(res)

    def test_tair_style_ids_and_the_gene_map_agree(self):
        pairs = self.pairs(TAIR, NCBI)
        by_suffix, err1 = run(P.run_pairs, TAIR, NCBI, self.tmp / "a", longest_isoform=True,
                              pairs=self.pairs_file(pairs))
        by_map, err2 = run(P.run_pairs, TAIR, NCBI, self.tmp / "b", longest_isoform=True, gene_map=GMAP,
                           pairs=self.pairs_file(pairs))
        self.assertIn("(from ID suffix .N)", err1)
        self.assertIn("(from --gene-map)", err2)
        self.assertEqual(by_suffix.exit_code, 0, err1)
        for name in ("genes", "dnds"):
            self.assertEqual(Path(by_suffix.outputs[name]).read_bytes(), Path(by_map.outputs[name]).read_bytes())
        ok = sorted(g["id"] for g in read_tsv(by_map.outputs["genes"]) if g["status"] == "ok" and g["set"] == "query")
        self.assertEqual(ok, KEPT[TAIR.name])
        self.assertEqual(by_map.report["isoforms"]["gene_map"], str(GMAP))
        self.assertIn("--gene-map", [i["option"] for i in by_map.report["inputs"]])

    def test_gene_map_alone_warns_and_partial_maps_are_reported(self):
        res, err = run(P.run_pairs, TAIR, Q, self.tmp / "o", gene_map=GMAP,
                       pairs=self.pairs_file([(KEPT[TAIR.name][0], "g01")]))
        self.assertIn("several transcripts per gene", err)
        self.assertIn("read from the --gene-map", err)
        res, err = run(P.run_pairs, TAIR, Q, self.tmp / "o2", gene_map=GMAP, longest_isoform=True,
                       pairs=self.pairs_file([(KEPT[TAIR.name][0], "g01")]))
        self.assertRegex(err, r"\d+ record\(s\) are not in --gene-map")   # the second file's records
        self.assertEqual(res.report["isoforms"]["sets"]["subject"]["dropped"], 0)

    def test_no_gene_information_keeps_everything_and_says_so(self):
        res, err = run(P.run_pairs, Q, S, self.tmp / "o", longest_isoform=True, pairs=self.pairs_file(PLAIN))
        self.assertIn("the headers of the query and subject file(s) name no gene", err)
        base, _ = run(P.run_pairs, Q, S, self.tmp / "base", pairs=self.pairs_file(PLAIN))
        self.assertEqual(Path(res.outputs["dnds"]).read_bytes(), Path(base.outputs["dnds"]).read_bytes())
        self.assertEqual(Path(res.outputs["genes"]).read_bytes(), Path(base.outputs["genes"]).read_bytes())

    def test_bad_gene_map(self):
        bad = self.tmp / "bad.tsv"
        bad.write_text("justone\n", encoding="utf-8")
        res, err = run(P.run_pairs, TAIR, Q, self.tmp / "o", gene_map=bad, pairs=self.pairs_file([("a", "b")]))
        self.assertEqual(res.exit_code, 1)
        self.assertIn("line 1: expected two columns (record ID, gene ID)", err)
        res, err = run(P.run_pairs, TAIR, Q, self.tmp / "o", gene_map=self.tmp / "none.tsv",
                       pairs=self.pairs_file([("a", "b")]))
        self.assertEqual(res.exit_code, 1)
        self.assertIn("cannot open", err)

    def test_cluster_mode(self):
        cl = self.tmp / "cl.tsv"
        kept = KEPT[TAIR.name]
        every = [r for r, _ in I.read_headers([TAIR])]
        # the clusters of a run made before the isoforms were removed: every transcript is a member
        cl.write_text("".join(f"{every[0]}\t{m}\n" for m in every), encoding="utf-8")
        res, err = run(P.run_clusters, TAIR, self.tmp / "o", clusters=cl, longest_isoform=True, max_cluster_size=2)
        self.assertIn("--longest-isoform: cds: 5 of 11 records dropped, 6 genes (from ID suffix .N)", err)
        summ = read_tsv(res.outputs["cluster_summary"])
        self.assertEqual((summ[0]["n_members"], summ[0]["n_used"]), ("11", "6"))
        used = {r["seq1"] for r in read_tsv(res.outputs["dnds"])} | {r["seq2"] for r in read_tsv(res.outputs["dnds"])}
        self.assertEqual(sorted(used), kept)
        failed = read_tsv(res.outputs["failed"])
        self.assertEqual(sum("dropped by --longest-isoform" in f["message"] for f in failed), 5)
        res, err = run(P.run_clusters, TAIR, self.tmp / "o2", clusters=cl, max_cluster_size=2)
        self.assertIn("cds: 11 records of 6 genes, 4 genes with several records (e.g. QT1G00010: QT1G00010.1, "
                      "QT1G00010.2, QT1G00010.3)", err)

    def test_command_line(self):
        out = self.tmp / "o"
        with redirect_stderr(io.StringIO()) as err:
            rc = P.main(["pairs", "--query", str(TAIR), "--subject", str(NCBI), "-o", str(out), "--longest-isoform",
                         "--gene-map", str(GMAP), "--pairs", str(self.pairs_file(self.pairs(TAIR, NCBI)))])
        self.assertEqual(rc, 0, err.getvalue())
        self.assertEqual(sum(g["status"] == "dropped" for g in read_tsv(out / "genes.tsv")), 10)


@needs_diamond
class WithDiamond(PipelineCase):
    def test_search_with_and_without_isoforms(self):
        plain, err = run(P.run_pairs, ENS, NCBI, self.tmp / "all", jobs=2)
        self.assertIn("several transcripts per gene", err)
        res, err = run(P.run_pairs, ENS, NCBI, self.tmp / "one", jobs=2, longest_isoform=True)
        self.assertEqual(res.exit_code, 0, err)
        pairs = [(r["seq1"], r["seq2"]) for r in read_tsv(res.outputs["dnds"])]
        want = [(EXPECTED[(ENS.name, g)], EXPECTED[(NCBI.name, h)])
                for (f, g), (f2, h) in zip([k for k in EXPECTED if k[0] == ENS.name],
                                           [k for k in EXPECTED if k[0] == NCBI.name])]
        self.assertEqual(sorted(pairs), sorted(want))        # one pair per gene, between the longest records
        self.assertEqual(res.report["counts"]["unpaired"], {"no_hit": 0, "not_reciprocal": 0})
        # without the option the transcripts compete: every record of a gene with several is paired or left over
        n_all = len(read_tsv(plain.outputs["dnds"]))
        self.assertGreaterEqual(n_all, 6)
        self.assertGreater(sum(plain.report["counts"]["unpaired"].values()), 0)
        # --resume reuses the search of the kept proteins
        again, _ = run(P.run_pairs, ENS, NCBI, self.tmp / "one", jobs=2, longest_isoform=True, resume=True)
        self.assertTrue(next(s for s in again.report["stages"] if s["name"] == "diamond_search")["cached"])
        self.assertEqual(Path(again.outputs["dnds"]).read_bytes(), Path(res.outputs["dnds"]).read_bytes())

    def test_deepclust_clusters_the_kept_records_only(self):
        res, err = run(P.run_clusters, TAIR, self.tmp / "o", jobs=2, longest_isoform=True, max_cluster_size=2)
        self.assertIn(res.exit_code, (0, 2), err)
        clusters, _ = P.read_clusters(res.outputs["clusters"])
        members = sorted(m for mem in clusters.values() for m in mem)
        self.assertEqual(members, KEPT[TAIR.name])


if __name__ == "__main__":
    unittest.main()
