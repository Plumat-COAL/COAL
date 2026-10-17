"""FIN-*: last review round (methods text, genetic-code hint, viewer, commands, MCP options)."""
import os
import random
import re
import shutil
import tempfile
import unittest

from coal_app import core

CLEAN_P, CLEAN_D = core.example_paths("clean")


def _read_fasta(path):
    recs, name = [], None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
                recs.append([name, ""])
            elif name is not None:
                recs[-1][1] += line.strip()
    return recs


def _write_fasta(path, recs):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.writelines(f">{n}\n{s}\n" for n, s in recs)


class Case(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="coal_fin_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def variant(self, prot=None, dna=None, name="v"):
        """examples/clean with the protein / DNA records changed by the given functions (list -> list)."""
        p = _read_fasta(CLEAN_P)
        d = _read_fasta(CLEAN_D)
        if prot:
            p = prot(p)
        if dna:
            d = dna(d)
        pp, dd = os.path.join(self.d, f"{name}.pep.fa"), os.path.join(self.d, f"{name}.dna.fa")
        _write_fasta(pp, p)
        _write_fasta(dd, d)
        return pp, dd

    def run_app(self, p, d, probe=False, **kw):
        return core.run_coal(p, d, out_dir=os.path.join(self.d, "run_" + str(random.random())[2:]),
                             probe_codes=probe, **kw)


def _randomise(i):
    def f(recs):
        rng = random.Random(7)
        recs[i][1] = "ATG" + "".join(rng.choice("ACGT") for _ in recs[i][1][3:])
        return recs
    return f


class TestFinC1MethodsText(Case):  # FIN-C1 / FIN-A2
    def test_suspect_gene_is_written_not_excluded(self):
        r = self.run_app(*self.variant(dna=_randomise(2)))
        self.assertEqual(r.report["summary"]["failed"], 1)
        with open(r.output_path) as fh:
            self.assertEqual(fh.read().count(">"), 8)  # the suspect gene IS in the output
        m = core.methods_text(r.report)
        self.assertNotIn("excluded", m)
        self.assertIn("1 was written to the output but flagged as probably wrong", m)
        self.assertIn("remove it before analysis", m)

    def test_non_standard_residue_is_not_a_mismatch(self):
        def ubz(recs):
            s = recs[0][1]
            i = s.index("W")  # gene01 has a W: replace by U (selenocysteine), B and Z further on
            s = s[:i] + "U" + s[i + 1:]
            recs[0][1] = s
            return recs
        r = self.run_app(*self.variant(prot=ubz))
        self.assertEqual(r.report["summary"]["warning"], 1, r.stderr)
        m = core.methods_text(r.report)
        self.assertNotIn("did not match", m)
        self.assertIn("non-standard amino acids", m)

    def test_pairing_only_warning_is_not_a_mismatch(self):
        def ren(recs):
            recs[2][0], recs[4][0] = "geneX3", "geneX5"
            return recs
        r = self.run_app(*self.variant(dna=ren))
        m = core.methods_text(r.report)
        self.assertNotIn("did not match", m)
        self.assertIn("2 were paired with their DNA only by their order in the files", m)

    def test_real_mismatch_still_counted(self):
        def mm(recs):
            s = recs[0][1]
            recs[0][1] = s[:9] + ("TTT" if s[9:12] != "TTT" else "GGG") + s[12:]  # 4th codon changed
            return recs
        r = self.run_app(*self.variant(dna=mm))
        m = core.methods_text(r.report)
        self.assertIn("1 contained codons that did not match the aligned amino acid", m)

    def test_proper_names_keep_their_capitals(self):  # FIN-C11
        rep = core.load_sample_report()
        for n, frag in ((22, "Scenedesmus obliquus mitochondrial genetic code"),
                        (23, "Thraustochytrium mitochondrial genetic code"),
                        (15, "Blepharisma nuclear genetic code"),
                        (6, "ciliate, dasycladacean and Hexamita nuclear genetic code"),
                        (4, "Mycoplasma/Spiroplasma"),
                        (1, "standard genetic code (NCBI translation table 1)")):
            r = dict(rep, options=dict(rep["options"], codontable=n))
            self.assertIn(frag, core.methods_text(r))


if __name__ == "__main__":
    unittest.main()
