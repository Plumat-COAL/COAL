"""Generate the simulation grid as AXT files (resumable: existing complete files are skipped).

One file per (freqs, kappa, omega, length) holding every dS x replicate:
  data/f{freqs}_k{kappa}_w{omega}_L{len}.axt, name line  d{dS}_r{rep}&b
Truth table: data/truth.tsv (one row per file x dS: t, dS, dN, omega, S, N, rhoS).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from codonsim import FREQS, CodonModel, sample_pair_P, to_str  # noqa: E402

DS = [0.01, 0.03, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
OMEGA = [0.05, 0.2, 0.5, 1.0, 2.0]
KAPPA = [1.0, 2.0, 4.0]
LENS = [30, 100, 300, 1000]
REPS = 200
OUT = os.path.join(os.environ.get("COAL_SIM_DIR", "work/sim_truth"), "data")


def fname(f, k, w, L):
    return f"{OUT}/f{f}_k{k:g}_w{w:g}_L{L}.axt"


def main():
    os.makedirs(OUT, exist_ok=True)
    truth_rows = []
    for fi, f in enumerate(FREQS):
        for ki, k in enumerate(KAPPA):
            for wi, w in enumerate(OMEGA):
                m = CodonModel(k, w, FREQS[f])
                for L in LENS:
                    for d in DS:
                        t = m.t_for_dS(d)
                        truth_rows.append((os.path.basename(fname(f, k, w, L)), f, k, w, L, d, t, d,
                                           m.dN_of_t(t), m.S, m.N, m.rhoS))
                    path = fname(f, k, w, L)
                    if os.path.exists(path):
                        continue
                    seed = [20261114, fi, ki, wi, L]
                    rng = np.random.default_rng(seed)
                    tmp = path + ".tmp"
                    with open(tmp, "w") as fh:
                        for d in DS:
                            t = m.t_for_dS(d)
                            for r in range(REPS):
                                a, b = sample_pair_P(m, t, L, rng)
                                fh.write(f"d{d:g}_r{r}&b\n{to_str(a)}\n{to_str(b)}\n\n")
                    os.replace(tmp, path)
                print(f, k, w, flush=True)
    with open(f"{OUT}/truth.tsv", "w") as fh:
        fh.write("file\tfreqs\tkappa\tomega\tlen\tdS_label\tt\tdS\tdN\tS\tN\trhoS\n")
        for r in truth_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


if __name__ == "__main__":
    main()
