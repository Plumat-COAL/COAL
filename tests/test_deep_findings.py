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
from unittest import mock

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
        # the map in prefixed form, with the same gene IDs in both species (retest R02: one record of each gene
        # was dropped as "not the longest isoform" of the other species' record, exit 0); a gene ID that already
        # carries the prefix stays as it is
        for name, lines in (("pref", [f"{x}|{n}\tG{k}\n" for x in ("sp_one", "sp_two") for k, n in enumerate(ids)]),
                            ("both", [f"{x}|{n}\t{x}|G{k}\n" for x in ("sp_one", "sp_two")
                                      for k, n in enumerate(ids)])):
            gm.write_text("".join(lines), encoding="utf-8")
            res = P.run_clusters([one, two], self.tmp / name, clusters=cl, max_cluster_size=2, prefix_ids=True,
                                 gene_map=gm, longest_isoform=True, quiet=True)
            self.assertEqual(res.exit_code, 0, name)
            iso = res.report["isoforms"]["sets"]["cds"]
            self.assertEqual((iso["genes"], iso["dropped"]), (2 * len(ids), 0), name)
            self.assertEqual(len(read_tsv(res.outputs["dnds"])), len(ids), name)
        # two transcripts of one gene of one species, named in prefixed form: still one gene
        gm.write_text(f"sp_one|{ids[0]}\tG0\nsp_one|{ids[1]}\tG0\n", encoding="utf-8")
        res = P.run_clusters([one, two], self.tmp / "iso", clusters=cl, max_cluster_size=2, prefix_ids=True,
                             gene_map=gm, longest_isoform=True, quiet=True)
        self.assertEqual(res.report["isoforms"]["sets"]["cds"]["dropped"], 1)

    def test_the_same_file_twice_is_named(self):   # retest R04 (cosmetic)
        one, _, _ = self.two_species()
        link = self.tmp / "also_one.fa"
        try:
            os.symlink(one, link)
        except OSError:
            link = None
        for prefix in (False, True):
            with redirect_stderr(io.StringIO()) as err:
                res = P.run_clusters([one, one], self.tmp / f"o{prefix}", clusters=FAM_CLUSTERS, prefix_ids=prefix)
            self.assertEqual(res.exit_code, 1)
            self.assertIn(f"--cds names {one} twice. Each file", err.getvalue())
            self.assertNotIn("two input files", err.getvalue())
        if link is not None:
            with redirect_stderr(io.StringIO()) as err:
                res = P.run_clusters([one, link], self.tmp / "o3", clusters=FAM_CLUSTERS)
            self.assertIn(f"--cds names the same file twice: {one} and {link}", err.getvalue())

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


# ---------------------------------------------------------------------------------------------- D03
BOM = b"\xef\xbb\xbf"


class ByteOrderMarkAndCarriageReturns(unittest.TestCase):
    """D03: a UTF-8 byte order mark (Windows Notepad) and classic Mac line ends (a lone CR) are read by every
    Python-side reader of the pipeline as coal-cpp reads them: the first record / line is not lost."""

    def setUp(self):
        import tempfile
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def variants(self, text: str):
        """name -> bytes of text in each style."""
        b = text.encode("utf-8")
        return {"plain": b, "bom": BOM + b, "crlf": b.replace(b"\n", b"\r\n"), "cr": b.replace(b"\n", b"\r"),
                "bom_cr": BOM + b.replace(b"\n", b"\r"), "mixed": b.replace(b"\n", b"\r", 2)}

    def test_fasta_lines(self):
        data = BOM + b">a x\rAC\r\nGT\n>b\rTT"
        got = list(P.fasta_lines(io.BytesIO(data)))
        self.assertEqual([ln for _, ln in got], [b">a x\r", b"AC\r\n", b"GT\n", b">b\r", b"TT"])
        self.assertEqual([data[o:o + len(ln)] for o, ln in got], [ln for _, ln in got])   # offsets into the file

    def test_fasta_readers(self):
        for name, data in self.variants(">s1 desc\nATGAAA\nCCC\n>s2\nATGGGG\n").items():
            f = self.tmp / f"{name}.fa"
            f.write_bytes(data)
            idx = P.FastaIndex(str(f))
            self.assertEqual(sorted(idx.pos), ["s1", "s2"], name)
            self.assertEqual((idx.get("s1"), idx.get("s2")), (b"ATGAAACCC", b"ATGGGG"), name)
            idx.close()
            from coal import isoforms
            self.assertEqual(isoforms.read_headers([str(f)]), [("s1", "desc"), ("s2", "")], name)
            self.assertEqual(P.fasta_ids(str(f)), [b"s1", b"s2"], name)

    def test_table_readers(self):
        from coal import batch, isoforms
        for name, data in self.variants("t1\tG1\nt2\tG1\n").items():
            f = self.tmp / f"{name}.tsv"
            f.write_bytes(data)
            self.assertEqual(isoforms.read_gene_map(str(f)), {"t1": "G1", "t2": "G1"}, name)
            self.assertEqual(P._read_pairs_file(str(f)), [("t1", "G1"), ("t2", "G1")], name)
            self.assertEqual(P.read_clusters(str(f)), ({"t1": ["t1", "G1"], "t2": ["t2"]}, 1), name)
            self.assertEqual(batch.read_jobs(str(f)), [(1, "t1\tG1"), (2, "t2\tG1")], name)


@needs_cpp
class ByteOrderMarkRuns(PipelineCase):
    def test_the_same_results_with_a_bom_or_lone_cr(self):
        from test_dnds_pipeline import PLAIN, Q, S
        pairs = "".join(f"{a}\t{b}\n" for a, b in PLAIN).encode()
        want = None
        for name, (fa, pr) in {"plain": (b"", b"\n"), "bom": (BOM, b"\n"), "cr": (b"", b"\r"),
                               "bom_cr": (BOM, b"\r")}.items():
            q = self.tmp / f"{name}.fa"
            q.write_bytes(fa + Q.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", pr))
            p = self.tmp / f"{name}.tsv"
            p.write_bytes(fa + pairs.replace(b"\n", pr))
            res = P.run_pairs(q, S, self.tmp / name, pairs=p, quiet=True)
            got = (res.exit_code, Path(res.outputs["dnds"]).read_bytes(),
                   Path(res.outputs["failed"]).read_bytes().replace(str(q).encode(), b"Q"))
            if want is None:
                want = got
            self.assertEqual(got, want, name)
            got = got[1]
            self.assertEqual(len(got.splitlines()), 1 + len(PLAIN))

# ---------------------------------------------------------------------------------------------- D01
@needs_cpp
class YnPamlBoundary(unittest.TestCase):
    """D01: tests/kaks/yn_paml_boundary_pair.fasta is an S. cerevisiae / C. glabrata ortholog pair (NCBI RefSeq,
    the codon alignment of `coal pipeline`, stop codons and ambiguous columns removed) whose 4-fold sites have
    P = 11/34 and Q = 12/34, so 1 - 2P - Q is exactly 0. PAML yn00 4.10.10 (weighting 1): S 128.0, N 529.0,
    dN 0.4123 (K80 fails, kappa 0.7704 from the non-degenerate sites); KaKs_Calculator 2.0 YN: S 194.22, N 462.78,
    Ka 0.502369, Ks 2.15976 (a rounding residue of 1e-17 gives kappa 60 at the 4-fold sites)."""

    def rows(self, *opts):
        import subprocess
        r = subprocess.run([str(CPP), "dnds", str(ROOT / "tests" / "kaks" / "yn_paml_boundary_pair.fasta"),
                            "--pairs", "consecutive", "-m", "YN", "-o", "-", "-q", *opts], capture_output=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        lines = r.stdout.decode().splitlines()
        return dict(zip(lines[0].split("\t"), lines[1].split("\t")))

    def test_paml_flavour_is_yn00(self):
        for w in ("1", "0"):
            r = self.rows("--yn-flavour", "paml", "--yn-weighting", w)
            self.assertEqual((round(float(r["S_sites"]), 1), round(float(r["N_sites"]), 1)), (128.0, 529.0), w)
            self.assertAlmostEqual(float(r["dN"]), {"1": 0.4123, "0": 0.4278}[w], places=4)
            # yn00 prints dS 3.6389; the Ks correction is undefined, which COAL reports as saturated (docs/kaks.md)
            self.assertEqual(r["status"], "saturated")

    def test_default_flavour_is_still_kaks_calculator(self):
        r = self.rows()
        self.assertEqual((r["S_sites"], r["N_sites"], r["dN"], r["dS"]),
                         ("194.22024", "462.77976", "0.50236926", "2.1597608"))

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

    # retest of 2026-11-14, R01: the code of the WHOLE set
    def mycoplasma(self, n=200, nuclear=60):
        """A Mycoplasma-like genome on one molecule: most genes read TGA as Trp, 15% start with GTG / TTG."""
        starts = ["ATG"] * 17 + ["GTG", "TTG", "TTG"]
        return [cds(self.rng, start=starts[k % 20]) for k in range(n - nuclear)], nuclear

    def test_a_genome_read_with_the_wrong_code_gets_its_code(self):
        orgs, nuc = self.mycoplasma()
        lcl = lambda k: ""   # noqa: E731
        h = self.hint(orgs, nuclear=nuc, chosen=11, header=lcl)
        # no molecule named: the codes whose start codons do not fit 15% of the genes drop out, 4, 5, 13 stay
        self.assertEqual((h["codes"], h["whole_set"]), ([4, 5, 13], True))
        # one molecule with 200 records is a genome, not an organelle: code 4 (the old replicon rule, "at least
        # half of the molecule's records are evidence", dropped every gene: 140 of 200 here, 253 of 524 in
        # M. genitalium, and the hint was silent)
        recs = [(f"lcl|NC_000908.2_cds_N{k}_{k}", cds(self.rng, stops=())) for k in range(nuc)]
        recs += [(f"lcl|NC_000908.2_cds_O{k}_{nuc + k}", s) for k, s in enumerate(orgs)]
        path = write_fasta(self.tmp / "myco.fa", recs)
        flagged = {r: 0 for r, _ in recs[nuc:][::2]}   # (half of them, as a sample would be)
        genes, bg, reps = P._hint_inputs([str(path)], flagged, {r for r, _ in recs[:nuc]})
        h = P.suggest_genetic_code(11, genes, bg, len(recs), len(orgs), reps, fraction=0.5)
        self.assertEqual((h["codes"], h["whole_set"]), ([4], True))
        # a mitochondrial set alone (few records on its molecule): organelle codes stay candidates
        h = self.hint([cds(self.rng, avoid=("CGA", "CGC")) for _ in range(5)], nuclear=2,
                      org_header=lambda k: "", header=lambda k: "")
        self.assertEqual((h["codes"], h["whole_set"]), ([3], True))

    def test_organelle_genes_in_a_nuclear_set_are_not_the_whole_set(self):
        h = self.hint([cds(self.rng, avoid=("CGA", "CGC")) for _ in range(5)])
        self.assertFalse(h["whole_set"])



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
        self.assertFalse(h["whole_set"])

    def run_set(self, recs, **options):
        q = write_fasta(self.tmp / "q.fa", recs)
        cl = self.tmp / "cl.tsv"
        cl.write_text(f"{recs[0][0]}\t{recs[1][0]}\n", encoding="utf-8")
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_clusters(q, self.tmp / "out", clusters=cl, **options)
        validate(res.report, SCHEMA)
        return res.report["genetic_code_hints"], " ".join(err.getvalue().split())

    def test_a_genome_with_the_wrong_code_is_told_to_run_again(self):
        """R01: M. genitalium with -c 11 got no hint at all; now 'the whole set looks like code 4'."""
        import random
        self.rng = random.Random(7)
        starts = ["ATG"] * 17 + ["GTG", "TTG", "TTG"]
        recs = [(f"lcl|NC_000908.2_cds_N{k}_{k}", cds(self.rng, stops=())) for k in range(60)]
        recs += [(f"lcl|NC_000908.2_cds_O{k}_{60 + k}", cds(self.rng, start=starts[k % 20])) for k in range(140)]
        hints, err = self.run_set(recs, codontable=11)
        h = hints[0]
        self.assertEqual((h["whole_set"], h["candidates"], h["clean_under"], h["records"],
                          h["genes_with_internal_stops"], h["genes"]), (True, [4], 4, 200, 140, 140))
        self.assertIn("cds: most genes (140 of the 200) have internal stop codons under genetic code 11, and 140 of the "
                      "140 with two or more have none under NCBI genetic code 4", err)
        self.assertIn("The whole set looks like code 4: its rows are estimates from mistranslated codons. Run "
                      "again with -c 4", err)
        self.assertNotIn("take them out", err)
        self.assertNotIn("nuclear set", err)

    def test_a_large_set_is_sampled_evenly(self):
        """R01: the first CODE_HINT_MAX_GENES genes of a record-shuffled ciliate set lay on many scaffolds (81 of
        26,927 genes, organelle wording); now a sample spread over the set, and the whole set's code."""
        import random
        self.rng = random.Random(7)
        recs = [(f"lcl|NW_{k:05d}.1_cds_X{k}_1", cds(self.rng, stops=("TAA", "TAG"))) for k in range(40)]
        recs += [(f"lcl|NW_{k:05d}.1_cds_N{k}_2", cds(self.rng, stops=())) for k in range(40, 50)]
        with mock.patch.object(P._Run, "CODE_HINT_MAX_GENES", 8):
            hints, err = self.run_set(recs)
        h = hints[0]
        self.assertEqual((h["whole_set"], h["candidates"], h["genes_looked_at"], h["genes"]), (True, [6], 8, 8))
        self.assertIn("most genes (40 of the 50) have internal stop codons under genetic code 1, and 8 of the 8 "
                      "looked at (of 40 with two or more) have none under NCBI genetic code 6 (ciliate nuclear)", err)
        self.assertIn("Run again with -c 6", err)
        # spread over the set: the 8 genes are every fifth, not the first 8
        self.assertIn("lcl|NW_00004.1_cds_X4_1", h["examples"])
        self.assertNotIn("lcl|NW_00001.1_cds_X1_1", h["examples"])


# ---------------------------------------------------------------------------------------------- D04, D06, D07, D08
class Messages(unittest.TestCase):
    def test_batch_line_with_the_program_name(self):   # D04
        import subprocess
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            jobs = Path(td) / "jobs.txt"
            jobs.write_text("coal pep.aln nuc.fa -o out.fa\n", encoding="utf-8")
            r = subprocess.run([sys.executable, "-m", "coal", "--batch", str(jobs)], capture_output=True, text=True,
                               cwd=td, env=dict(os.environ, PYTHONPATH=str(ROOT / "src")))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("the line starts with 'coal': a job line is a coal command line without the program name",
                      r.stderr)
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn("without the program name", readme)

    def test_mistyped_command(self):   # D06
        import subprocess
        from coal import cli
        self.assertEqual([cli.mistyped_command(w) for w in ("dnd", "piepline", "kak", "Dnds", "codon-dns", "data",
                                                             "x", "pep.aln", "cite")],
                         ["dnds", "pipeline", "kaks", "dnds", "codon-dnds", "", "", "", ""])
        want = ("coal: error: cannot open 'dnd': no such file, and no command of that name (did you mean "
                "`coal dnds`?)")
        r = subprocess.run([sys.executable, "-m", "coal", "dnd", "x.fa"], capture_output=True, text=True,
                           env=dict(os.environ, PYTHONPATH=str(ROOT / "src")))
        self.assertEqual((r.returncode, r.stderr.strip()), (1, want))
        if CPP:   # the same in coal-cpp
            r = run_cpp("dnd", "x.fa")
            self.assertEqual((r.returncode, r.stderr.decode().strip()), (1, want))

    def test_mafft_is_pinned(self):   # D07
        from unittest import mock
        from coal import msa
        self.assertEqual(msa.MAFFT_VERSION, "7.526")
        calls = []

        def fake_install(tool, spec, conda, quiet, timeout):
            calls.append(spec)
            return "/x/mafft"
        with mock.patch.object(msa, "_conda_install", fake_install), \
                mock.patch.object(msa, "_conda_exe", return_value="/x/conda"), \
                mock.patch.object(msa, "_trusted_cached", return_value=False), \
                mock.patch.object(msa.platform, "system", return_value="Linux"), \
                mock.patch.object(msa, "_install_lock"), \
                mock.patch.object(msa, "tool_version", return_value="v7.526 (2024/Apr/26)"), \
                redirect_stderr(io.StringIO()) as err:
            self.assertEqual(msa.install_mafft(), "/x/mafft")
        self.assertEqual(calls, ["mafft=7.526"])
        self.assertIn("MAFFT 7.526 installed: /x/mafft", err.getvalue())

    @unittest.skipIf(os.name == "nt", "a shell script stands in for MAFFT")
    def test_an_unpinned_mafft_in_the_cache_is_named(self):   # retest R03
        """An earlier COAL installed MAFFT unpinned: install kept it without a word, and runs used it."""
        import tempfile
        from coal import msa
        with tempfile.TemporaryDirectory() as td, mock.patch.dict(os.environ, {"COAL_CACHE_DIR": td}):
            exe = fake_mafft(msa, "v7.525 (2024/Mar/13)")
            calls = []
            with mock.patch.object(msa, "_conda_install", lambda *a: calls.append(a)), \
                    redirect_stderr(io.StringIO()) as err:
                self.assertEqual(msa.install_mafft(), exe)
            self.assertEqual(calls, [])     # kept, not replaced (--force replaces it)
            self.assertIn(f"coal: note: MAFFT 7.525 in COAL's cache ({exe}) is not the pinned 7.526 (an earlier "
                          "COAL installed it unpinned); it is used as it is. To install 7.526: coal pipeline "
                          "install mafft --force", err.getvalue())
            with redirect_stderr(io.StringIO()) as err:
                msa.install_mafft(quiet=True)
            self.assertEqual(err.getvalue(), "")
            # the pinned version, and a MAFFT of your own (not in the cache): nothing to say
            self.assertIsNone(msa.pin_note("mafft", exe, "v7.526 (2024/Apr/26)"))
            self.assertIsNone(msa.pin_note("mafft", "/usr/bin/mafft", "v7.505"))
            self.assertIn("Clustal Omega 1.2.3", msa.pin_note("clustalo", msa._cached("clustalo")[1], "1.2.3"))
            fake_mafft(msa, "v7.526 (2024/Apr/26)")
            with redirect_stderr(io.StringIO()) as err:
                msa.install_mafft()
            self.assertEqual(err.getvalue(), "")


def fake_mafft(msa, version: str) -> str:
    """A 'MAFFT' in COAL's cache, with the checksum record an install writes: --version prints `version`, and an
    alignment is its input (the test's sequences have one length)."""
    exe = msa._cached("mafft")[1]
    os.makedirs(os.path.dirname(exe), exist_ok=True)
    with open(exe, "w", encoding="utf-8") as fh:
        fh.write(f'#!/bin/sh\nif [ "$1" = --version ]; then echo "{version}" >&2; exit 0; fi\n'
                 'for a; do last=$a; done\nexec /bin/cat "$last"\n')
    os.chmod(exe, 0o755)
    msa._write_manifest(exe, "mafft=7.525 (test)", "test")
    return exe


@needs_cpp
class PipelineMessages(PipelineCase):
    def test_pairs_file_swapped_repeated_and_empty(self):   # D06
        from test_dnds_pipeline import PLAIN, Q, S
        out = io.StringIO()
        with redirect_stderr(out):
            res = P.run_pairs(Q, S, self.tmp / "sw", pairs=self.pairs_file([(b, a) for a, b in PLAIN[:4]]))
        self.assertIn("look swapped: in 4 of the 4 pairs the first ID is a record of --subject", out.getvalue())
        self.assertTrue(any("look swapped" in w for w in res.report["warnings"]))
        out = io.StringIO()
        with redirect_stderr(out):
            P.run_pairs(Q, S, self.tmp / "dup", pairs=self.pairs_file([*PLAIN[:3], PLAIN[0], PLAIN[0]]))
        self.assertIn("2 line(s) of --pairs", out.getvalue())
        self.assertNotIn("look swapped", out.getvalue())
        empty = self.tmp / "empty.tsv"
        empty.write_text("query\tsubject\n", encoding="utf-8")
        out = io.StringIO()
        with redirect_stderr(out):
            res = P.run_pairs(Q, S, self.tmp / "none", pairs=empty)
        self.assertIn("no pairs: the --pairs file", out.getvalue())

    def test_a_bad_tool_variable_is_named(self):   # D06
        from unittest import mock
        with mock.patch.dict(os.environ, {"COAL_CLUSTALO": str(self.tmp / "nope")}), \
                redirect_stderr(io.StringIO()) as err:
            res = P.run_clusters(FAM, self.tmp / "o", clusters=FAM_CLUSTERS, max_cluster_size=8)
        self.assertEqual(res.exit_code, 1)
        self.assertIn(f"$COAL_CLUSTALO is set to '{self.tmp / 'nope'}', which is not a file", err.getvalue())
        with mock.patch.dict(os.environ, {"COAL_DIAMOND": str(self.tmp)}), redirect_stderr(io.StringIO()) as err:
            res = P.run_clusters(FAM, self.tmp / "o2", max_cluster_size=8)
        self.assertIn("$COAL_DIAMOND is set to", err.getvalue())
        self.assertIn("(a folder)", err.getvalue())

    @unittest.skipIf(os.name == "nt", "a shell script stands in for MAFFT")
    def test_a_run_names_an_unpinned_mafft(self):   # retest R03
        import random
        from coal import msa
        rng = random.Random(5)
        a = cds(rng, stops=())

        def mutant():   # (the same length: the stand-in's "alignment" is its input)
            b = "".join(c if rng.random() > 0.1 else rng.choice("ACGT") for c in a[:-3])
            return "".join(b[k:k + 3] if b[k:k + 3] not in ("TAA", "TAG", "TGA") else a[k:k + 3]
                           for k in range(0, len(b), 3)) + "TAA"
        q = write_fasta(self.tmp / "q.fa", [("g1", a), ("g2", mutant()), ("g3", mutant())])
        cl = self.tmp / "cl.tsv"   # (three members: two are aligned without an MSA program)
        cl.write_text("g1\tg2\ng1\tg3\n", encoding="utf-8")
        env = {"COAL_CACHE_DIR": str(self.tmp / "cache"), "PATH": str(self.tmp / "nothing")}
        with mock.patch.dict(os.environ, env):
            os.environ.pop("COAL_MAFFT", None)
            exe = fake_mafft(msa, "v7.525 (2024/Mar/13)")
            with redirect_stderr(io.StringIO()) as err:
                res = P.run_clusters(q, self.tmp / "o", clusters=cl, aligner="mafft")
        self.assertEqual(res.exit_code, 0, err.getvalue())
        self.assertEqual(res.report["tools"]["mafft"]["path"], exe)
        note = f"MAFFT 7.525 in COAL's cache ({exe}) is not the pinned 7.526"
        self.assertIn(note, err.getvalue())
        self.assertTrue(any(note in w for w in res.report["warnings"]))

    def test_most_pairs_saturated(self):   # D08
        import random
        rng = random.Random(3)
        recs_q, recs_s, pairs = [], [], []
        for k in range(25):   # unrelated sequences of one length: synonymous sites saturated
            recs_q.append((f"q{k}", cds(rng, n=150, stops=())))
            recs_s.append((f"s{k}", cds(rng, n=150, stops=())))
            pairs.append((f"q{k}", f"s{k}"))
        q, s = write_fasta(self.tmp / "q.fa", recs_q), write_fasta(self.tmp / "s.fa", recs_s)
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_pairs(q, s, self.tmp / "o", pairs=self.pairs_file(pairs))
        c = res.report["counts"]["pairs"]
        lost = c["saturated"] + c["not_converged"]
        self.assertGreater(lost, 0.3 * c["total"])
        self.assertIn(f"{lost} of the 25 pairs", err.getvalue())
        self.assertIn("Filtering before interpretation", err.getvalue())
        # the closely related pairs of the fixture: no such line
        from test_dnds_pipeline import PLAIN, Q, S
        with redirect_stderr(io.StringIO()) as err:
            P.run_pairs(Q, S, self.tmp / "ok", pairs=self.pairs_file(PLAIN))
        self.assertNotIn("too divergent", err.getvalue())

# ---------------------------------------------------------------------------------------------- memory
@needs_cpp
class CheckpointedTraceback(PipelineCase):
    """A pair whose traceback is above AlignParams::traceback_bytes is aligned with checkpoints (titin: 630 MB of
    traceback before). The test hook COAL_TEST_TRACEBACK_BYTES makes every pair of the fixture take that route,
    with blocks of one row and of a few rows: every output byte must equal the plain run's."""

    def test_the_same_bytes(self):
        from unittest import mock
        from cpp_hooks import need_test_hooks
        from test_dnds_pipeline import PLAIN, Q, S
        need_test_hooks(self, CPP)
        outs = {}
        for budget in (None, "1", "100", "4096"):
            env = {"COAL_TEST_HOOKS": "1", **({"COAL_TEST_TRACEBACK_BYTES": budget} if budget else {})}
            with mock.patch.dict(os.environ, env):
                if budget is None:
                    os.environ.pop("COAL_TEST_TRACEBACK_BYTES", None)
                res = P.run_pairs(Q, S, self.tmp / f"o{budget}", pairs=self.pairs_file(PLAIN), quiet=True,
                                  keep_alignments=True, keep_protein_alignments=True, preset="pwalign")
            outs[budget] = [Path(res.outputs[k]).read_bytes() for k in ("dnds", "codon_alignments",
                                                                        "protein_alignments")]
        for budget in ("1", "100", "4096"):
            self.assertEqual(outs[budget], outs[None], budget)

if __name__ == "__main__":
    unittest.main()
