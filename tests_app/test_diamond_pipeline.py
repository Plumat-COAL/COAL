import importlib.util
import os
import random
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "demo", "diamond_pipeline.py")
STUB = os.path.join(ROOT, "tests_app", "diamond_stub.py")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mx = load("make_examples", os.path.join(ROOT, "examples", "make_examples.py"))
dp = load("diamond_pipeline", SCRIPT)


def make_input(d):
    """Species A (with isoforms .1/.2) and B: 5 orthologous genes + 1 A-only gene."""
    rng = random.Random(3)
    cmap = mx.codon_map(mx.CODE1)
    ap, ac, bp, bc = [], [], [], []
    for g in range(5):
        a_row, b_row = mx.make_family(rng, 2, 120)
        a, b = a_row.replace("-", ""), b_row.replace("-", "")
        a_short = a[:60]
        for tid, prot in ((f"Hs{g}.1", a_short), (f"Hs{g}.2", a)):
            ap.append((tid, prot))
            ac.append((tid, mx.back_translate(rng, prot, cmap)))
        bp.append((f"Mm{g}-PA", b))
        bc.append((f"Mm{g}-PA", mx.back_translate(rng, b, cmap)))
    lone = "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(80))
    ap.append(("HsLone.1", lone))
    ac.append(("HsLone.1", mx.back_translate(rng, lone, cmap)))
    paths = {}
    for k, recs in (("ap", ap), ("ac", ac), ("bp", bp), ("bc", bc)):
        paths[k] = os.path.join(d, k + ".fa")
        mx.write_fasta(paths[k], recs)
    return paths


class TestDiamondPipeline(unittest.TestCase):
    def test_helpers(self):
        self.assertEqual(dp.gene_of("AT1G01010.1", "AT1G01010.1 x"), "AT1G01010")
        self.assertEqual(dp.gene_of("ENSP1", "ENSP1 pep gene:ENSG9 transcript:ENST1"), "ENSG9")
        self.assertEqual(dp.gene_of("CG1-PA", "CG1-PA"), "CG1")
        self.assertEqual(dp.reciprocal_best_hits({"a": "x", "b": "y"}, {"x": "a", "y": "c"}), [("a", "x")])
        ra, rb = dp.needleman_wunsch("MKTAYIAKQRQISFVK", "MKTAYIAKQRISFVK")
        self.assertEqual((ra.replace("-", ""), rb.replace("-", "")), ("MKTAYIAKQRQISFVK", "MKTAYIAKQRISFVK"))
        self.assertEqual(len(ra), len(rb))
        self.assertEqual(rb.count("-"), 1)

    def test_no_diamond_clear_error(self):
        with tempfile.TemporaryDirectory() as t:
            p = make_input(t)
            r = subprocess.run([sys.executable, SCRIPT, "--a-protein", p["ap"], "--a-cds", p["ac"],
                                "--b-protein", p["bp"], "--b-cds", p["bc"], "--out", os.path.join(t, "o"),
                                "--diamond", os.path.join(t, "no_such_diamond")],
                               capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 1)
        self.assertIn("DIAMOND was not found", r.stderr)

    def test_pipeline_with_stand_in(self):
        with tempfile.TemporaryDirectory() as t:
            p = make_input(t)
            out = os.path.join(t, "o")
            r = subprocess.run([sys.executable, SCRIPT, "--a-protein", p["ap"], "--a-cds", p["ac"], "--a-name", "human",
                                "--b-protein", p["bp"], "--b-cds", p["bc"], "--b-name", "mouse",
                                "--out", out, "--diamond", STUB, "--threads", "1"],
                               capture_output=True, text=True, timeout=300)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertIn("6 genes", r.stdout)            # 11 proteins -> 6 genes
            self.assertIn("reciprocal best hits: 5", r.stdout)
            pairs = sorted(os.listdir(os.path.join(out, "pairs")))
            self.assertEqual(len(pairs), 5)
            self.assertIn("human__Hs0.2__Mm0-PA.codon.fasta", pairs)  # the longest isoform was used
            # each pair: two codon rows of equal length that translate back to the proteins
            prots = {**{i: s for i, _, s in dp.read_fasta(p["ap"])}, **{i: s for i, _, s in dp.read_fasta(p["bp"])}}
            code = dict(zip(mx.CODONS, mx.CODE1))
            for f in pairs:
                recs = dp.read_fasta(os.path.join(out, "pairs", f))
                self.assertEqual(len(recs), 2)
                self.assertEqual(len(recs[0][2]), len(recs[1][2]))
                for rid, _, seq in recs:
                    tid = rid.split("__", 1)[1]
                    aa = "".join(code[seq[i:i + 3]] for i in range(0, len(seq), 3) if seq[i:i + 3] != "---")
                    self.assertEqual(aa, prots[tid])


if __name__ == "__main__":
    unittest.main()
