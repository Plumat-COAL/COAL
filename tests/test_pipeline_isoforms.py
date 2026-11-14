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
        # the NCBI-style file has [gene=] and [db_xref=GeneID:]: the stable identifier wins over the symbol
        for path, source in ((ENS, "header field gene:"), (NCBI, "header field GeneID:"), (TAIR, "ID suffix .N")):
            d = I.detect(I.read_headers([path]))
            self.assertEqual((d.source, d.records, d.genes, d.multi, d.extra), (source, 11, 6, 4, 5), path.name)
            self.assertEqual(len(d.example[1]), 3)
        self.assertEqual(I.detect(I.read_headers([Q])).source, None)   # plain IDs: nothing to find

    def test_header_fields(self):
        d = self.det(["t1 gene:G1 x", "t2 gene:G1", "t3 gene:G2"])
        self.assertEqual((d.source, d.gene_of), ("header field gene:", {"t1": "G1", "t2": "G1", "t3": "G2"}))
        d = self.det(["a [locus_tag=B0001] [protein=x y]", "b [locus_tag=B0001]", "c [locus_tag=B0002]"])
        self.assertEqual((d.source, d.multi), ("header field [locus_tag=", 1))
        # stable identifiers before symbols: a locus tag, a GeneID, Ensembl's gene:, FlyBase's parent=FBgn
        d = self.det(["a [gene=tuf] [locus_tag=X_001]", "b [gene=tuf] [locus_tag=X_900]", "c [gene=rpoB] [locus_tag=X_2]"])
        self.assertEqual((d.source, d.kind, d.genes, d.multi), ("header field [locus_tag=", "header", 3, 0))
        d = self.det(["a [gene=A] [db_xref=CCDS:CCDS1.1,GeneID:11]", "b [gene=A] [db_xref=GeneID:11]", "c GeneID:12"])
        self.assertEqual((d.source, d.gene_of["a"], d.gene_of["c"], d.multi), ("header field GeneID:", "GeneID:11",
                                                                              "GeneID:12", 1))
        d = self.det(["p1 type=CDS; parent=FBgn0000001,FBtr0000001;", "p2 type=CDS; parent=FBgn0000001,FBtr0000002;"])
        self.assertEqual((d.source, d.gene_of["p2"]), ("header field parent=FBgn", "FBgn0000001"))
        d = self.det(["a gene=SYM locus_tag=L1", "b gene=SYM locus_tag=L2"])
        self.assertEqual((d.source, d.multi), ("header field locus_tag=", 0))
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

    def test_a_symbol_shared_by_different_loci_is_not_merged(self):
        # two strains in one file (another contig in the record ID)
        d = I.detect(I.read_headers([FIX / "multi_strain_cds.fasta"]))
        self.assertEqual((d.source, d.kind, d.records, d.genes, d.multi), ("header field [gene=", "symbol", 6, 6, 0))
        self.assertEqual(d.conflicts, 3)
        self.assertEqual(d.conflict_example[0], "dnaA")
        self.assertIn("contigs STRAIN1, STRAIN2", d.conflict_example[2])
        self.assertEqual(d.gene_of, {})
        # one contig, the symbol on most records, locus tags on too few to be the field: they still tell loci apart
        d = self.det(["lcl|c_cds_1 [gene=tuf] [locus_tag=X_001]", "lcl|c_cds_2 [gene=tuf] [locus_tag=X_900]"]
                     + [f"lcl|c_cds_{k} [gene=g{k}]" for k in range(3, 30)]
                     + ["lcl|c_cds_30 [gene=iso]", "lcl|c_cds_31 [gene=iso]"])
        self.assertEqual((d.source, d.conflicts, d.multi), ("header field [gene=", 1, 1))
        self.assertIn("locus tags X_001, X_900", d.conflict_example[2])
        self.assertEqual(d.example[0], "iso")          # transcripts of one locus are still merged
        # the same symbol on one contig with nothing that tells the records apart: transcripts
        d = self.det(["lcl|c_cds_1 [gene=A]", "lcl|c_cds_2 [gene=A]"])
        self.assertEqual((d.conflicts, d.multi), (0, 1))

    def test_every_documented_style(self):
        """One fixture per row of docs/pipeline.md "Header styles" (coal.isoforms.STYLES): the source, the gene,
        and the record kept (the second: the first is the shorter transcript)."""
        rows = read_tsv(FIX / "styles" / "expected.tsv")
        self.assertEqual([(r["style"], r["source"]) for r in rows], [(s, src) for s, _, src in I.STYLES])
        doc = (ROOT / "docs" / "pipeline.md").read_text(encoding="utf-8")
        for r, (style, (h1, h2), source) in zip(rows, I.STYLES):
            path = FIX / "styles" / r["file"]
            self.assertNotIn(b"\r", path.read_bytes())
            heads = I.read_headers([path])
            self.assertEqual([h for h, _ in heads[:2]], [h1.split()[0], h2.split()[0]])
            d = I.detect(heads)
            self.assertEqual((d.source, d.records, d.genes, d.multi, d.extra), (source, 3, 2, 1, 1), style)
            self.assertEqual(d.gene_of[heads[0][0]], r["gene"])
            lengths = {}
            for line in path.read_text(encoding="utf-8").splitlines():
                if line.startswith(">"):
                    cur = line[1:].split()[0]
                    lengths[cur] = 0
                else:
                    lengths[cur] += len(line)
            self.assertEqual(I.choose_longest(d, lengths), {heads[0][0]: (r["gene"], r["kept"])}, style)
            self.assertIn(f"| {style} |", doc, "docs/pipeline.md lacks the style")
            self.assertIn("`" + h1.split()[0].replace("|", "\\|") + "`", doc, style)
        rules = [name for name, _ in I.SUFFIX_RULES]
        self.assertLess(rules.index("ID suffix _iN.pN (TransDecoder on Trinity)"), rules.index("ID suffix _iN (Trinity)"))
        self.assertEqual({src for _, _, src in I.STYLES if src.startswith("ID suffix")}, set(rules))
        self.assertLessEqual({src[len("header field "):] for _, _, src in I.STYLES if src.startswith("header")},
                             {name for name, _, _ in I.HEADER_FIELDS})

    def test_refseq_accessions_hold_no_gene(self):
        d = I.detect(I.read_headers([FIX / "refseq_rna.fasta"]))
        self.assertEqual((d.source, d.records, d.refseq_ids), (None, 3, 3))
        self.assertEqual(self.det(["g1", "g2"]).refseq_ids, 0)

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
            # a record named twice with the same gene is fine
            p.write_text("# comment\ntranscript_id\tgene_id\nt1\tG1\nt2\tG1\tignored\n\nt3 G2\nt1\tG1\n", encoding="utf-8")
            m = I.read_gene_map(str(p))
            self.assertEqual(m, {"t1": "G1", "t2": "G1", "t3": "G2"})
            d = I.from_gene_map(m, ["t1", "t2", "t3", "t4"])
            self.assertEqual((d.source, d.genes, d.multi, d.extra, d.in_map), ("gene map", 3, 1, 1, 3))
            # with two different genes it is an error that lists the records
            p.write_text("t1\tG1\nt2\tG1\nt1\tOTHER\nt2\tG1\nt2\tG9\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, r"2 record\(s\) are given two different genes: t1 \(line 3: OTHER; "
                                                    r"before: G1\); t2 \(line 5: G9; before: G1\)"):
                I.read_gene_map(str(p))
            p.write_text("t1\tG1\nonly_one_column\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "line 2: expected two columns"):
                I.read_gene_map(str(p))

    def test_gene_map_wins_per_record_and_the_rest_falls_back(self):
        heads = [("AT1G01010.1", ""), ("AT1G01010.2", ""), ("AT1G01020.1", ""), ("AT1G01020.2", ""), ("AT1G01030.1", "")]
        # the map names one record: its gene is the map's; the others are found by their suffix as without a map
        d = I.from_gene_map({"AT1G01010.1": "MINE"}, heads)
        self.assertEqual((d.source, d.kind, d.fallback, d.in_map), ("gene map", "map", "ID suffix .N", 1))
        self.assertEqual(d.gene_of["AT1G01010.1"], "MINE")
        self.assertEqual(d.gene_of["AT1G01020.2"], "AT1G01020")
        self.assertEqual((d.genes, d.multi), (4, 1))
        self.assertEqual(d.by_suffix, {"AT1G01010.2", "AT1G01020.1", "AT1G01020.2"})
        # a map gene with the name of a detected one is that gene
        d = I.from_gene_map({"AT1G01030.1": "AT1G01020"}, heads)
        self.assertEqual((d.genes, d.multi, d.extra), (2, 2, 3))
        # a map that names no record of the set: the set is read as without a map
        d = I.from_gene_map({"other": "G"}, heads)
        self.assertEqual((d.source, d.in_map, d.multi), ("ID suffix .N", 0, 2))
        # --isoform-rule map: nothing is read from headers or IDs
        d = I.from_gene_map({"AT1G01010.1": "MINE"}, heads, detection=False)
        self.assertEqual((d.source, d.fallback, d.genes, d.multi), ("gene map", None, 5, 0))

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
        self.assertIn("header field GeneID:", lines[0])
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
                      "subject: 5 of 11 records dropped, 6 genes (from header field GeneID:); 10 in all, listed in "
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
        self.assertEqual(Path(by_suffix.outputs["dnds"]).read_bytes(), Path(by_map.outputs["dnds"]).read_bytes())

        def without_gene_names(path):   # the subject's genes are GeneID:9001 from the header, HOM1 from the map
            return [{k: (v.split(": ", 1)[-1] if k == "note" else v) for k, v in g.items()} for g in read_tsv(path)]
        self.assertEqual(without_gene_names(by_suffix.outputs["genes"]), without_gene_names(by_map.outputs["genes"]))
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
        # the map names no record of the second file, and nothing else does: said once, everything kept there
        self.assertIn("the headers of the subject file(s) name no gene", err)
        self.assertIn(f"and --gene-map {GMAP} names none of their records", err)
        sets = res.report["isoforms"]["sets"]
        self.assertEqual((sets["subject"]["dropped"], sets["subject"]["records_in_gene_map"]), (0, 0))
        self.assertEqual((sets["query"]["dropped"], sets["query"]["records_in_gene_map"]), (5, 11))

    def test_gene_map_for_one_set_leaves_detection_on_for_the_other(self):
        """A map of the query's transcripts only: the subject's isoforms are still found from its headers (the
        map used to switch detection off for every set)."""
        qmap = self.tmp / "q.tsv"
        qmap.write_text("".join(ln + "\n" for ln in GMAP.read_text(encoding="utf-8").splitlines() if ln[:2] == "QT"),
                        encoding="utf-8", newline="\n")
        pairs = self.pairs(TAIR, NCBI)
        res, err = run(P.run_pairs, TAIR, NCBI, self.tmp / "o", longest_isoform=True, gene_map=qmap,
                       pairs=self.pairs_file(pairs))
        self.assertEqual(res.exit_code, 0, err)
        self.assertIn("query: 5 of 11 records dropped, 6 genes (from --gene-map); subject: 5 of 11 records dropped, "
                      "6 genes (from header field GeneID:)", err)
        self.assertIn(f"--gene-map {qmap} names no record of the subject file(s): their genes are read from the "
                      "header field GeneID:", err)
        base, _ = run(P.run_pairs, TAIR, NCBI, self.tmp / "base", longest_isoform=True, pairs=self.pairs_file(pairs))
        self.assertEqual(Path(res.outputs["dnds"]).read_bytes(), Path(base.outputs["dnds"]).read_bytes())
        # a map of some records of a set: the map for those, the suffix rule for the rest
        part = self.tmp / "part.tsv"
        part.write_text("QT1G00010.1\tQT1G00010\nQT1G00010.2\tQT1G00010\n", encoding="utf-8", newline="\n")
        res, err = run(P.run_pairs, TAIR, NCBI, self.tmp / "p", longest_isoform=True, gene_map=part,
                       pairs=self.pairs_file(pairs))
        self.assertIn("query: 5 of 11 records dropped, 6 genes (from --gene-map, and ID suffix .N for the records it "
                      "does not name)", err)
        self.assertIn("9 record(s) are not in --gene-map", err)
        q = res.report["isoforms"]["sets"]["query"]
        self.assertEqual((q["records_in_gene_map"], q["records_not_in_gene_map"],
                          q["gene_source_of_records_not_in_gene_map"]), (2, 9, "ID suffix .N"))
        self.assertEqual(Path(res.outputs["dnds"]).read_bytes(), Path(base.outputs["dnds"]).read_bytes())
        # --isoform-rule map: the old behaviour, nothing but the map
        res, err = run(P.run_pairs, TAIR, NCBI, self.tmp / "m", longest_isoform=True, gene_map=part,
                       isoform_rule="map", pairs=self.pairs_file(pairs[:1]))
        self.assertIn("query: 1 of 11 records dropped, 10 genes (from --gene-map)", err)   # .3 is not in the map
        self.assertEqual(res.report["isoforms"]["isoform_rule"], "map")

    def test_gene_map_that_matches_nothing_is_an_error(self):
        m = self.tmp / "m.tsv"
        m.write_text("ENST1\tG1\nENST2\tG1\n", encoding="utf-8", newline="\n")
        for extra in ({}, {"longest_isoform": True}):
            res, err = run(P.run_pairs, TAIR, NCBI, self.tmp / "o", gene_map=m, pairs=self.pairs_file([("a", "b")]),
                           **extra)
            self.assertEqual(res.exit_code, 1, err)
            self.assertIn(f"--gene-map {m} names none of the records of the input", err)
            self.assertIn("e.g. ENST1, ENST2", err)                       # the map's IDs
            self.assertIn("QT1G00010.1, QT1G00010.2, QT1G00010.3 (query)", err)      # and the input's, of both sets
            self.assertIn("lcl|NC_000001.1_cds_XP_000101.1_1", err)
            self.assertNotIn("Traceback", err)
            self.assertFalse((self.tmp / "o" / "dnds.tsv").exists())
        res, err = run(P.run_clusters, TAIR, self.tmp / "c", gene_map=m, clusters=self.tmp / "m.tsv")
        self.assertEqual(res.exit_code, 1, err)
        self.assertIn("names none of the records of the input", err)

    def test_gene_map_duplicates(self):
        m = self.tmp / "m.tsv"
        lines = GMAP.read_text(encoding="utf-8").splitlines()
        m.write_text("\n".join(lines + lines[1:4]) + "\n", encoding="utf-8", newline="\n")   # identical: fine
        pairs = self.pairs(TAIR, NCBI)
        res, err = run(P.run_pairs, TAIR, NCBI, self.tmp / "o", gene_map=m, longest_isoform=True,
                       pairs=self.pairs_file(pairs))
        self.assertEqual(res.exit_code, 0, err)
        m.write_text("\n".join(lines + ["QT1G00010.1\tELSE", "QT1G00010.2\tQT1G00010"]) + "\n", encoding="utf-8",
                     newline="\n")
        res, err = run(P.run_pairs, TAIR, NCBI, self.tmp / "o2", gene_map=m, longest_isoform=True,
                       pairs=self.pairs_file(pairs))
        self.assertEqual(res.exit_code, 1)
        self.assertIn("1 record(s) are given two different genes: QT1G00010.1 (line 24: ELSE; before: QT1G00010)", err)

    def test_suffix_alone_never_drops_most_of_a_file(self):
        """Genes numbered per scaffold look like transcripts to the .N rule (three per stem passes the guard):
        --longest-isoform refuses instead of merging paralogs, unless --isoform-rule suffix says so."""
        scaf = FIX / "scaffold_numbered_cds.fasta"
        Q = ENS            # a subject without unusable records, so that the exit status is the run's own
        pairs = self.pairs_file([("s1.1", KEPT[ENS.name][0])])
        res, err = run(P.run_pairs, scaf, Q, self.tmp / "o", longest_isoform=True, pairs=pairs)
        self.assertEqual(res.exit_code, 1, err)
        self.assertIn("--longest-isoform would drop 4 of the 6 query records on the evidence of their IDs alone "
                      "(ID suffix .N; e.g. s1.1, s1.2, s1.3 would be one gene, s1): more than half", err)
        self.assertIn("--isoform-rule suffix", err)
        self.assertIn("--gene-map TSV", err)
        self.assertIn("Nothing was dropped", err)
        self.assertFalse((self.tmp / "o" / "dnds.tsv").exists())
        self.assertEqual(res.report["status"], "failed")
        # without --longest-isoform nothing is dropped, so the run goes on with the usual warning
        res, err = run(P.run_pairs, scaf, Q, self.tmp / "w", pairs=pairs)
        self.assertEqual(res.exit_code, 0, err)
        self.assertIn("looks like several transcripts per gene", err)
        # said explicitly, the suffixes are transcripts
        res, err = run(P.run_pairs, scaf, Q, self.tmp / "s", longest_isoform=True, isoform_rule="suffix", pairs=pairs)
        self.assertNotEqual(res.exit_code, 1, err)     # (2 when the pair names a record that was dropped)
        self.assertIn("query: 4 of 6 records dropped, 2 genes (from ID suffix .N)", err)
        q = res.report["isoforms"]["sets"]["query"]
        self.assertEqual((q["dropped"], q["dropped_by_id_suffix"]), (4, 4))
        self.assertEqual(res.report["isoforms"]["isoform_rule"], "suffix")
        # a gene map takes the guess out: no refusal
        m = self.tmp / "m.tsv"
        m.write_text("".join(f"s{a}.{b}\tscaffold{a}_gene{b}\n" for a in (1, 2) for b in (1, 2, 3)), encoding="utf-8",
                     newline="\n")
        res, err = run(P.run_pairs, scaf, Q, self.tmp / "m", longest_isoform=True, gene_map=m, pairs=pairs)
        self.assertEqual(res.exit_code, 0, err)
        self.assertEqual(res.report["isoforms"]["sets"]["query"]["dropped"], 0)
        # the same guard in cluster mode
        cl = self.tmp / "cl.tsv"
        cl.write_text("s1.1\ts1.1\ns1.1\ts2.1\n", encoding="utf-8", newline="\n")
        res, err = run(P.run_clusters, scaf, self.tmp / "c", clusters=cl, longest_isoform=True)
        self.assertEqual(res.exit_code, 1)
        self.assertIn("would drop 4 of the 6 cds records", err)

    def test_documented_failure_mode_of_the_suffix_rule(self):
        """HSP70.1 / HSP70.2 among other genes: one of six records, so no refusal; they are merged (docs/pipeline.md
        says so, and the summary line names the rule so that the user can see it)."""
        res, err = run(P.run_pairs, FIX / "paralog_numbered_cds.fasta", ENS, self.tmp / "o", longest_isoform=True,
                       pairs=self.pairs_file([("ACT1", KEPT[ENS.name][0])]))
        self.assertEqual(res.exit_code, 0, err)
        self.assertIn("query: 1 of 6 records dropped, 5 genes (from ID suffix .N)", err)
        self.assertIn("HSP70.1 / HSP70.2", (ROOT / "docs" / "pipeline.md").read_text(encoding="utf-8"))

    def test_symbols_at_different_loci_and_refseq_ids_are_said(self):
        res, err = run(P.run_pairs, FIX / "multi_strain_cds.fasta", FIX / "refseq_rna.fasta", self.tmp / "o",
                       longest_isoform=True, pairs=self.pairs_file([("lcl|STRAIN1_cds_WP_1001.1_1", "NM_000001.1")]))
        self.assertEqual(res.exit_code, 0, err)
        self.assertIn("query: 3 gene symbol(s) are shared by records at different loci (e.g. dnaA: "
                      "lcl|STRAIN1_cds_WP_1001.1_1, lcl|STRAIN2_cds_WP_2001.1_1; contigs STRAIN1, STRAIN2)", err)
        self.assertIn("they are not merged", err)
        self.assertIn("RefSeq transcript accessions (NM_, XM_) hold no gene in the ID, and the headers of the subject "
                      "file(s) have no gene field", err)
        sets = res.report["isoforms"]["sets"]
        self.assertEqual((sets["query"]["dropped"], sets["query"]["gene_symbols_at_different_loci"]), (0, 3))
        self.assertEqual((sets["subject"]["dropped"], sets["subject"]["refseq_transcript_ids"]), (0, 3))
        self.assertEqual({g["status"] for g in read_tsv(res.outputs["genes"])}, {"ok"})
        self.check_report(res)
        # without --longest-isoform the RefSeq note is a warning too (nothing else can say it)
        res, err = run(P.run_pairs, FIX / "multi_strain_cds.fasta", FIX / "refseq_rna.fasta", self.tmp / "w",
                       pairs=self.pairs_file([("lcl|STRAIN1_cds_WP_1001.1_1", "NM_000001.1")]))
        self.assertIn("RefSeq transcript accessions", err)

    def test_records_not_genes_and_counts_that_agree(self):
        """Before anything is collapsed the lines count records; a repeated ID is one record that is not usable in
        the first line and in the summary alike."""
        dup = self.tmp / "dup.fasta"
        text = TAIR.read_text(encoding="utf-8")
        first = text.split(">")[1]
        dup.write_text(text + ">" + first, encoding="utf-8", newline="\n")
        res, err = run(P.run_pairs, dup, NCBI, self.tmp / "o", longest_isoform=True,
                       pairs=self.pairs_file(self.pairs(TAIR, NCBI)))
        self.assertIn("translated 12 query and 11 subject records (1 not usable)", err)
        self.assertIn("summary: 23 records (1 not usable, 10 dropped by --longest-isoform: 12 genes); ", err)
        res, err = run(P.run_pairs, dup, NCBI, self.tmp / "p", pairs=self.pairs_file(self.pairs(TAIR, NCBI)))
        self.assertIn("summary: 23 records (1 not usable); ", err)
        self.assertNotRegex(err, r"summary: \d+ genes")
        res, err = run(P.run_clusters, dup, self.tmp / "c", clusters=self.stop_after_translation())
        self.assertIn("translated 12 records (1 not usable)", err)

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
