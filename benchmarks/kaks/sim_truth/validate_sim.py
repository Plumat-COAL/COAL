"""Validate codonsim: model algebra, P(t) sampler vs Gillespie path sampler, realised substitution counts vs theory."""
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from codonsim import FREQS, CodonModel, sample_pair_P, sample_pair_path, SENSE  # noqa: E402


def ndiff(a, b):
    out = np.zeros(len(a), dtype=int)
    for k, div in enumerate((16, 4, 1)):
        ca = np.array([SENSE[x] for x in a]) // div % 4
        cb = np.array([SENSE[x] for x in b]) // div % 4
        out += ca != cb
    return out


def main():
    rng = np.random.default_rng(20261114)
    n = 200_000
    print("model checks | kappa omega freqs | S N rhoS | max|piQ| max|detailed balance| rowsum(Q)")
    rows = []
    for kappa in (1, 4):
        for omega in (0.2, 2):
            for fname in ("equal", "skewed"):
                m = CodonModel(kappa, omega, FREQS[fname])
                piQ = np.abs(m.pi @ m.Q).max()
                db = np.abs(m.pi[:, None] * m.Q - (m.pi[:, None] * m.Q).T).max()
                rs = np.abs(m.Q.sum(1)).max()
                scale = -(m.pi * np.diag(m.Q)).sum()
                print(f"  {kappa} {omega} {fname:6s} | S={m.S:.4f} N={m.N:.4f} rhoS={m.rhoS:.4f} | {piQ:.1e} {db:.1e} {rs:.1e} rate={scale:.6f}")
                for dS in (0.1, 1.0, 3.0):
                    t = m.t_for_dS(dS)
                    r0, d0, ns, nn = sample_pair_path(m, t, n, rng)
                    r1, d1 = sample_pair_P(m, t, n, rng)
                    # realised substitutions per site vs truth
                    dS_real = ns.sum() / (n * m.S)
                    dN_real = nn.sum() / (n * m.N)
                    se_s = np.sqrt(ns.var() * n) / (n * m.S)
                    se_n = np.sqrt(nn.var() * n) / (n * m.N)
                    zS = (dS_real - dS) / se_s
                    zN = (dN_real - m.dN_of_t(t)) / se_n
                    # endpoint: number of differing positions, path vs P(t)
                    k0 = np.bincount(ndiff(r0, d0), minlength=4)
                    k1 = np.bincount(ndiff(r1, d1), minlength=4)
                    chi = stats.chi2_contingency(np.vstack([k0, k1]) + 0.5)
                    # identical-codon fraction vs theory sum_i pi_i P_ii(t)
                    pid_th = float((m.pi * np.diag(m.P(t))).sum())
                    pid_P = float((r1 == d1).mean())
                    se_pid = np.sqrt(pid_th * (1 - pid_th) / n)
                    # root distribution vs pi
                    cnt = np.bincount(r1, minlength=m.n)
                    chi_pi = stats.chisquare(cnt, m.pi * n)
                    rows.append((kappa, omega, fname, dS, t, dS_real, zS, m.dN_of_t(t), dN_real, zN, chi.pvalue,
                                 pid_th, pid_P, (pid_P - pid_th) / se_pid, chi_pi.pvalue))
    print()
    print("kappa omega freqs dS_true t | dS_realised z | dN_true dN_realised z | p(chi2 path vs P: #diff positions) |"
          " P(identical) theory sampled z | p(root ~ pi)")
    for r in rows:
        print("  %g %g %-6s %4.2f %6.3f | %.4f %+5.2f | %.4f %.4f %+5.2f | %.3f | %.4f %.4f %+5.2f | %.3f" % r)
    zs = np.array([[r[6], r[9], r[13]] for r in rows])
    ps = np.array([[r[10], r[14]] for r in rows])
    print(f"\nmax |z| = {np.abs(zs).max():.2f} over {zs.size} z-scores; min p = {ps.min():.3f} over {ps.size} tests"
          f" (Bonferroni 0.05 threshold {0.05/ps.size:.4f})")


if __name__ == "__main__":
    main()
