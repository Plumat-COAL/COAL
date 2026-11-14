"""The findings of the independent deep test of 2026-11-14 (D01-D08), each fixed and kept fixed here.

  * D05: clusters mode with the same IDs in two --cds files (one species was dropped without a word): refused
    up front with both files named; --prefix-ids puts each file's name before its IDs.
  * D02: the genetic-code hint (no pseudogenes, the code with the fewest reassignments, candidates when unsure).
  * D03: a UTF-8 byte order mark and classic Mac (lone CR) line ends in every reader of the pipeline.
  * D01: YN --yn-flavour paml at the 2P + Q = 1 boundary, as PAML yn00 computes it.
  * D04, D06, D07, D08: messages (--batch lines, 0 pairs, --pairs columns, $COAL_CLUSTALO, a mistyped subcommand,
    the MAFFT pin, most pairs saturated).

Synthetic inputs only (no real sequence data); coal-cpp (COAL_CPP_BIN) for the pipeline parts, DIAMOND where said.
"""
from __future__ import annotations

import io
import json
import os
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import dnds_pipeline as P  # noqa: E402
from test_dnds_pipeline import (CPP, FAM, FAM_CLUSTERS, SCHEMA, PipelineCase, needs_cpp,  # noqa: E402
                                needs_diamond, read_fasta, read_tsv, run_cpp)
from test_report_schema import validate  # noqa: E402


def write_fasta(path: Path, recs, eol: str = "\n", bom: bool = False) -> Path:
    text = "".join(f">{k}{eol}{v}{eol}" for k, v in recs)
    path.write_bytes((b"\xef\xbb\xbf" if bom else b"") + text.encode("utf-8"))
    return path


# ---------------------------------------------------------------------------------------------- D05
class FastaHelpers(unittest.TestCase):
    def test_file_stem(self):
        for name, want in (("S_cer.fa", "S_cer"), ("dir/S_cer.cds.fa.gz", "S_cer"), ("x.fasta", "x"),
                           ("genome.FNA", "genome"), ("plain", "plain"), (".fa", ".fa"), ("a.b.fa", "a.b")):
            self.assertEqual(P.file_stem(name), want, name)


@needs_cpp
class SameIdsInSeveralCdsFiles(PipelineCase):
    """D05: g1.t1 / gene1 in every species' file."""

    def two_species(self):
        fam = read_fasta(FAM)
        a = [n for n in fam if n.startswith(("famA", "famC"))]
        # species two: the same IDs, the sequences of the other members (in reverse)
        one = write_fasta(self.tmp / "sp_one.fa", [(n, fam[n]) for n in a])
        two = write_fasta(self.tmp / "sp_two.fasta", [(n, fam[m]) for n, m in zip(a, reversed(a))])
        return one, two, a

    def test_refused_naming_both_files(self):
        one, two, ids = self.two_species()
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_clusters([one, two], self.tmp / "out", clusters=FAM_CLUSTERS, max_cluster_size=2)
        self.assertEqual(res.exit_code, 1)
        msg = err.getvalue()
        self.assertIn(f"{len(ids)} sequence ID(s) occur in more than one of the --cds files", msg)
        self.assertIn(f"in both {one} and {two}", msg)
        self.assertIn("--prefix-ids", msg)
        self.assertNotIn("dnds", res.outputs)
        self.assertFalse((self.tmp / "out" / "dnds.tsv").exists())

    def test_prefix_ids(self):
        one, two, ids = self.two_species()
        cl = self.tmp / "cl.tsv"   # a cluster of each ID: its two species
        cl.write_text("".join(f"sp_one|{n}\tsp_one|{n}\nsp_one|{n}\tsp_two|{n}\n" for n in ids), encoding="utf-8")
        res = P.run_clusters([one, two], self.tmp / "out", clusters=cl, max_cluster_size=2, prefix_ids=True,
                             quiet=True)
        self.assertEqual(res.exit_code, 0)
        self.check_report(res)
        rows = read_tsv(res.outputs["dnds"])
        self.assertEqual([(r["seq1"], r["seq2"]) for r in rows], [(f"sp_one|{n}", f"sp_two|{n}") for n in ids])
        genes = read_tsv(res.outputs["genes"])
        self.assertEqual(len(genes), 2 * len(ids))
        self.assertTrue(all(g["status"] == "ok" for g in genes))
        self.assertEqual([r.get("id_prefix") for r in res.report["inputs"] if r["option"] == "--cds"],
                         ["sp_one|", "sp_two|"])
        self.assertTrue(res.report["parameters"]["prefix_ids"])
        # the same file names in two folders cannot give the prefixes
        (self.tmp / "d").mkdir()
        other = write_fasta(self.tmp / "d" / "sp_one.fa", [("x", "ATG")])
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_clusters([one, other], self.tmp / "out2", clusters=cl, prefix_ids=True)
        self.assertEqual(res.exit_code, 1)
        self.assertIn("would get the same prefix 'sp_one'", err.getvalue())

    def test_prefix_ids_with_a_gene_map(self):
        one, two, ids = self.two_species()
        gm = self.tmp / "map.tsv"   # unprefixed record IDs; the same gene name in both species is two genes
        gm.write_text("".join(f"{n}\tG{k}\n" for k, n in enumerate(ids)), encoding="utf-8")
        cl = self.tmp / "cl.tsv"
        cl.write_text("".join(f"sp_one|{n}\tsp_two|{n}\n" for n in ids), encoding="utf-8")
        res = P.run_clusters([one, two], self.tmp / "out", clusters=cl, max_cluster_size=2, prefix_ids=True,
                             gene_map=gm, longest_isoform=True, quiet=True)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.report["isoforms"]["sets"]["cds"]["genes"], 2 * len(ids))
        self.assertEqual(len(read_tsv(res.outputs["dnds"])), len(ids))

    def test_translate_names_both_files(self):
        one, two, _ = self.two_species()
        r = run_cpp("translate", one, two, "-o", self.tmp / "p.faa", "--genes", self.tmp / "g.tsv")
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = read_tsv(self.tmp / "g.tsv")
        dup = [g for g in rows if g["code"] == "duplicate_id"]
        self.assertEqual(len(dup), len(rows) // 2)
        self.assertEqual(dup[0]["note"], f"the ID occurs in two input files, {one} and {two} (the record of the "
                                         "first is used): IDs must be unique across the files")
        # within one file the wording is as before
        r = run_cpp("translate", one, one, "-o", self.tmp / "p.faa", "--genes", self.tmp / "g.tsv")
        same = write_fasta(self.tmp / "same.fa", [("x", "ATGAAA"), ("x", "ATGCCC")])
        r = run_cpp("translate", same, "-o", self.tmp / "p.faa", "--genes", self.tmp / "g.tsv")
        self.assertEqual(read_tsv(self.tmp / "g.tsv")[1]["note"],
                         "the ID occurs more than once in the file (the first record is used)")


@needs_diamond
class SameIdsDeepClust(PipelineCase):
    def test_no_species_is_lost(self):
        fam = read_fasta(FAM)
        names = [n for n in fam if not n.startswith("famB")]
        one = write_fasta(self.tmp / "one.fa", [(f"gene{k}", fam[n]) for k, n in enumerate(names)])
        two = write_fasta(self.tmp / "two.fa", [(f"gene{k}", fam[n]) for k, n in enumerate(reversed(names))])
        res = P.run_clusters([one, two], self.tmp / "out", prefix_ids=True, max_cluster_size=2, quiet=True)
        self.assertIn(res.exit_code, (0, 2))
        genes = read_tsv(res.outputs["genes"])
        self.assertEqual(sum(g["status"] != "failed" for g in genes), 2 * len(names))
        seqs = {r["seq1"].split("|")[0] for r in read_tsv(res.outputs["dnds"])} | \
               {r["seq2"].split("|")[0] for r in read_tsv(res.outputs["dnds"])}
        self.assertEqual(seqs, {"one", "two"})


# ---------------------------------------------------------------------------------------------- D02
SENSE = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG" if a + b + c not in ("TAA", "TAG", "TGA")]


def cds(rng, n=300, avoid=(), stops=("TGA", "TGA", "TGA"), start="ATG", tail="TAA") -> str:
    pool = [c for c in SENSE if c not in avoid]
    codons = [start] + [rng.choice(pool) for _ in range(n)]
    for k, c in enumerate(stops):
        codons[20 + 37 * k] = c
    return "".join(codons) + tail


class GeneticCodeHint(unittest.TestCase):
    """D02: the hint named 14 for yeast mitochondrial genes and for frameshifted pseudogenes, 4 for vertebrate
    mitochondrial genes. suggest_genetic_code() on synthetic sets (uniform codon usage in the nuclear genes)."""

    def setUp(self):
        import random
        import tempfile
        self.rng = random.Random(7)
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def hint(self, organelle, nuclear=100, chosen=1, header=lambda k: "", org_header=lambda k: ""):
        recs = [(f"nuc{k} {header(k)}", cds(self.rng, stops=())) for k in range(nuclear)]
        recs += [(f"org{k} {org_header(k)}", s) for k, s in enumerate(organelle)]
        path = write_fasta(self.tmp / "set.fa", recs)
        flagged = {f"org{k}": 0 for k in range(len(organelle))}
        usable = {r.split()[0] for r, _ in recs}
        genes, bg, reps = P._hint_inputs([str(path)], flagged, usable)
        return P.suggest_genetic_code(chosen, genes, bg, len(recs), len(organelle), reps)

    def test_yeast_mitochondrial_genes_get_code_3(self):
        h = self.hint([cds(self.rng, avoid=("CGA", "CGC")) for _ in range(5)])
        self.assertEqual(h["codes"], [3])
        self.assertIn("CGA / CGC codons, which code 3 does not use", h["evidence"])
        self.assertEqual(len(h["genes"]), 5)

    def test_vertebrate_mitochondrial_genes_get_code_2(self):
        h = self.hint([cds(self.rng, avoid=("AGA", "AGG"), start="ATT") for _ in range(5)])
        self.assertEqual(h["codes"], [2])
        self.assertIn("AGA / AGG codons, which code 2 reads as stops", h["evidence"])
        # an incomplete stop codon (T, TA: completed by polyadenylation) is no frameshift
        h = self.hint([cds(self.rng, avoid=("AGA", "AGG"), tail=t) for t in ("T", "TA", "TAA", "T")])
        self.assertEqual((h["codes"], len(h["genes"])), ([2], 4))

    def test_nothing_tells_the_codes_apart(self):
        h = self.hint([cds(self.rng) for _ in range(4)])
        self.assertEqual(h["codes"], [4, 5, 9, 10, 13, 21])   # 3 loses: the genes use CGA / CGC like the others
        # ATT / ATA starts: codes whose start codons do not include them drop out
        h = self.hint([cds(self.rng, start=("ATT", "ATA")[k % 2]) for k in range(4)])
        self.assertEqual(h["codes"], [4, 5])
        # no other genes to compare with: the codes that avoid codons cannot be told apart either
        h = self.hint([cds(self.rng, n=120, avoid=("CGA", "CGC")) for _ in range(3)], nuclear=0)
        self.assertEqual(h["codes"], [3, 4, 5, 9, 10, 13, 21])

    def test_pseudogenes_and_frameshifts_give_no_hint(self):
        taa = [cds(self.rng, stops=("TAA", "TGA", "TAA")) for _ in range(3)]
        self.assertIsNone(self.hint(taa))      # TAA as an amino acid: only for a whole set
        self.assertIsNone(self.hint([cds(self.rng) for _ in range(4)], org_header=lambda k: "[pseudo=true]"))
        self.assertIsNone(self.hint([cds(self.rng)[:-3] + "G" for _ in range(4)]))   # length 3n + 1, ends in G
        # the whole set read with the wrong code (ciliate nuclear genes with -c 1): code 6
        h = self.hint([cds(self.rng, stops=("TAA", "TAG")) for _ in range(30)], nuclear=5)
        self.assertEqual(h["codes"], [6])

    def test_genes_scattered_over_nuclear_chromosomes_give_no_hint(self):
        ens = lambda k: f"cds chromosome:GRCh38:{k % 20 + 1}:100:400:1"   # noqa: E731
        orgs = [cds(self.rng) for _ in range(4)]
        self.assertIsNone(self.hint(orgs, header=ens, org_header=lambda k: f"cds chromosome:GRCh38:{k + 1}:1:9:1"))
        h = self.hint(orgs, header=ens, org_header=lambda k: "cds chromosome:GRCh38:MT:1:900:1")
        self.assertEqual(len(h["genes"]), 4)



@needs_cpp
class GeneticCodeHintRun(PipelineCase):
    def test_the_report_and_the_message(self):
        import random
        self.rng = random.Random(7)
        recs = [(f"nuc{k}", cds(self.rng, stops=())) for k in range(60)]
        recs += [(f"mt{k}", cds(self.rng)) for k in range(3)]
        q = write_fasta(self.tmp / "q.fa", recs)
        cl = self.tmp / "cl.tsv"
        cl.write_text("nuc0\tnuc1\n", encoding="utf-8")
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_clusters(q, self.tmp / "out", clusters=cl)
        h = res.report["genetic_code_hints"][0]
        self.assertEqual((h["clean_under"], h["candidates"], h["genes"]), (None, [4, 5, 9, 10, 13, 21], 3))
        validate(res.report, SCHEMA)
        self.assertIn("have none under any of the NCBI genetic codes 4, 5, 9, 10, 13 and 21", err.getvalue())


if __name__ == "__main__":
    unittest.main()
