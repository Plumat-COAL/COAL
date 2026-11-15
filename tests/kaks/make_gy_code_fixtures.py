"""Make tests/kaks/gy_codes.axt and tests/kaks/golden/gy_codes.codeml.tsv: pairs for method GY under the genetic
codes other than the standard one that PAML's codeml has, with codeml's results (tests/test_kaks_gy.py,
TestGeneticCodes).

    python tests/kaks/make_gy_code_fixtures.py --codeml /path/to/codeml

codeml's icode and the NCBI table it is: 1 = 2 (vertebrate mitochondrial), 2 = 3 (yeast mt), 3 = 4 (mold mt), 4 = 5
(invertebrate mt), 5 = 6 (ciliate nuclear), 6 = 9 (echinoderm mt), 7 = 10 (euplotid nuclear), 8 = 12 (alternative
yeast nuclear), 9 = 13 (ascidian mt), 10 = 15 (blepharisma nuclear). Two pairs per code: random sense codons of
that code (so that the codons whose meaning differs from the standard code occur: TGA, AGA / AGG, ATA, CTN, TAA /
TAG, AAA, CTG where they are sense codons), changed at random positions, never into a stop codon of the code.
codeml: runmode = -2, CodonFreq = 2, getSE = 1 (PAML 4.10.10 made the stored file). The NCBI tables codeml does not
have (11, 14, 16, 21, 22, 23) are checked in the test against the likelihood recomputed from the model.
"""
import argparse
import os
import random
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "benchmarks", "kaks"))
sys.path.insert(0, os.path.join(ROOT, "tests"))
import test_kaks_properties as TP  # noqa: E402  (the NCBI tables)
import validate_gy as V  # noqa: E402

ICODE = {1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 9, 7: 10, 8: 12, 9: 13, 10: 15}
BASES = "TCAG"


def make_pair(rng, code, n, div, keep_nonsyn):
    aa = TP.AA[code]
    a = [rng.choice(TP.SENSE[code]) for _ in range(n)]
    b = []
    for c in a:
        d = c
        for _ in range(3):
            if rng.random() < div:
                p = rng.randrange(3)
                e = d[:p] + rng.choice(BASES) + d[p + 1:]
                if aa[e] == "*" or (aa[e] != aa[d] and rng.random() > keep_nonsyn):
                    continue
                d = e
        b.append(d)
    return "".join(a), "".join(b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--codeml", required=True)
    a = ap.parse_args()
    base = V.CTL
    rows, pairs_all = [], []
    with tempfile.TemporaryDirectory() as td:
        for icode, ncbi in ICODE.items():
            rng = random.Random(7000 + ncbi)
            pairs = [(f"code{ncbi}_{k}",) + make_pair(rng, ncbi, 150, div, kn)
                     for k, (div, kn) in enumerate(((0.08, 0.3), (0.3, 0.6)))]
            V.CTL = base.replace("icode = 0", f"icode = {icode}")
            res, _ = V.run_codeml(a.codeml, pairs, td)
            for (nm, x, y), r in zip(pairs, res):
                if r is None:
                    raise SystemExit(f"codeml printed nothing for {nm}")
                pairs_all.append((nm, x, y))
                rows.append(f"{nm}\t{ncbi}\t{icode}\t{r['lnL']:.6f}\t{r['t']:.5f}\t{r['kappa']:.5f}\t{r['omega']:.5f}\t"
                            f"{r['dN']:.4f}\t{r['dS']:.4f}\t{r['S']:.1f}\t{r['N']:.1f}\t{r['se_dN']:.5f}\t"
                            f"{r['se_dS']:.5f}")
    with open(os.path.join(HERE, "gy_codes.axt"), "w", newline="\n") as fh:
        for nm, x, y in pairs_all:
            fh.write(f"{nm}&ref\n{x}\n{y}\n\n")
    with open(os.path.join(HERE, "golden", "gy_codes.codeml.tsv"), "w", newline="\n") as fh:
        fh.write("name\tncbi\ticode\tlnL\tt\tkappa\tomega\tdN\tdS\tS\tN\tSE_dN\tSE_dS\n" + "\n".join(rows) + "\n")
    print(len(rows), "pairs")


if __name__ == "__main__":
    main()
