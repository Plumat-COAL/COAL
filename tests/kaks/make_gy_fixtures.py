"""Make the input and the reference file of tests/test_kaks_gy.py: pairs for method GY with PAML codeml's results.

    python tests/kaks/make_gy_fixtures.py --codeml /path/to/codeml --sim SIM.axt --real REAL.axt [--real ...]

SIM.axt: pairs simulated under GY94 / F3x4 (benchmarks/kaks/validate_gy.py sim writes <out>/coal_in.axt: 8
combinations of codon usage, kappa and omega x dS 0.1 ... 3); REAL.axt: real ortholog pairs. A few pairs of each
are taken (fixed positions), cleaned to the codon columns COAL compares, and given to codeml (runmode = -2,
CodonFreq = 2 (F3x4), one data set per pair; PAML 4.10.10 made the stored file). Writes

    tests/kaks/gy_cases.axt                    the pairs
    tests/kaks/golden/gy_cases.codeml.tsv      name, icode, lnL, t, kappa, omega, dN, dS, S, N as codeml prints them

The last pairs are run a second time under the vertebrate mitochondrial code (codeml icode = 1, COAL -c 2), with
the columns that hold one of its stop codons (AGA, AGG) removed first: rows with icode 1.
"""
import argparse
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "benchmarks", "kaks"))
import validate_gy as V  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--codeml", required=True)
    ap.add_argument("--sim", required=True)
    ap.add_argument("--real", action="append", default=[])
    a = ap.parse_args()
    pairs = []
    sim = V.read_axt(a.sim)
    for i in range(0, len(sim), max(1, len(sim) // 10)):
        x, y = V.clean_pair(sim[i][1], sim[i][2])
        pairs.append((f"sim{i}", x[:600], y[:600]))        # 200 codons
    for k, path in enumerate(a.real):
        real = V.read_axt(path, 40)
        took = 0
        for j, (_, x, y) in enumerate(real):
            x, y = V.clean_pair(x, y)
            if 120 * 3 <= len(x) <= 260 * 3 and x != y and took < 3:
                pairs.append((f"real{k}_{j}", x, y))
                took += 1
    mito = []
    for nm, x, y in pairs[:3]:
        ox, oy = [], []
        for i in range(0, len(x), 3):
            if x[i:i + 3] not in ("AGA", "AGG") and y[i:i + 3] not in ("AGA", "AGG"):
                ox.append(x[i:i + 3])
                oy.append(y[i:i + 3])
        mito.append((nm + "_mt", "".join(ox), "".join(oy)))
    with tempfile.TemporaryDirectory() as td:
        std, _ = V.run_codeml(a.codeml, pairs, td)
        V.CTL = V.CTL.replace("icode = 0", "icode = 1")
        mt, _ = V.run_codeml(a.codeml, mito, td)
    with open(os.path.join(HERE, "gy_cases.axt"), "w", newline="\n") as fh:
        for nm, x, y in pairs + mito:
            fh.write(f"{nm}&ref\n{x}\n{y}\n\n")
    with open(os.path.join(HERE, "golden", "gy_cases.codeml.tsv"), "w", newline="\n") as fh:
        fh.write("name\ticode\tlnL\tt\tkappa\tomega\tdN\tdS\tS\tN\n")
        for icode, ps, rs in ((0, pairs, std), (1, mito, mt)):
            for (nm, _, _), r in zip(ps, rs):
                if r is None:
                    raise SystemExit(f"codeml printed nothing for {nm}")
                fh.write(f"{nm}\t{icode}\t{r['lnL']:.6f}\t{r['t']:.5f}\t{r['kappa']:.5f}\t{r['omega']:.5f}\t"
                         f"{r['dN']:.4f}\t{r['dS']:.4f}\t{r['S']:.1f}\t{r['N']:.1f}\n")
    print(len(pairs), "pairs,", len(mito), "under the mitochondrial code")


if __name__ == "__main__":
    main()
