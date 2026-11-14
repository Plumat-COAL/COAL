"""Figures for REPORT.md from res/cells.tsv (static PNGs; every figure has a table counterpart in the report)."""
import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

W = os.environ.get("COAL_SIM_DIR", "work/sim_truth")  # the work folder: data/, res/, plots/
d = pd.read_csv(f"{W}/res/cells.tsv", sep="\t")
d = d[d.run == "main"]
METHODS = ["YN", "MYN", "NG", "LWL", "MLWL", "LPB", "Comeron", "Li"]
COL = dict(zip(METHODS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]))
MRK = dict(zip(METHODS, ["o", "s", "^", "v", "D", "P", "X", "*"]))
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#888", "axes.grid": True, "grid.color": "#e6e6e6",
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb", "lines.linewidth": 2})


def panel_grid(metric, omega, L, fname, ylabel, ylim, ref=0.0, transform=None):
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), sharex=True, sharey=True)
    for i, f in enumerate(["equal", "skewed"]):
        for j, k in enumerate([1.0, 2.0, 4.0]):
            ax = axes[i, j]
            x = d[(d.freqs == f) & (d.kappa == k) & (d.omega == omega) & (d.len == L)]
            for m in METHODS:
                g = x[x.method == m].sort_values("dS_true")
                y = g[metric] if transform is None else transform(g)
                ax.plot(g.dS_true, y, color=COL[m], marker=MRK[m], markersize=5, label=m)
            ax.axhline(ref, color="#52514e", lw=1)
            ax.axhspan(ref - 0.1, ref + 0.1, color="#e6e6e6", alpha=0.5, lw=0)
            ax.set_xscale("log")
            ax.set_ylim(*ylim)
            ax.set_title(f"{f} codon usage, kappa = {k:g}", fontsize=9)
            if j == 0:
                ax.set_ylabel(ylabel)
            if i == 1:
                ax.set_xlabel("true dS")
    axes[0, 0].legend(ncol=4, fontsize=7.5, loc="lower left", frameon=False)
    fig.suptitle(f"{ylabel}: omega = {omega:g}, {L} codons, 200 pairs per point (grey band: +-10%)", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{W}/plots/{fname}", dpi=130)
    plt.close(fig)


def main():
    panel_grid("dS_medrel", 0.2, 1000, "fig1_dS_bias_L1000_w0.2.png", "relative error of median dS", (-0.6, 0.6))
    panel_grid("w_medrel", 0.2, 1000, "fig2_omega_bias_L1000_w0.2.png", "relative error of median dN/dS", (-0.6, 0.6))
    panel_grid("w_medrel", 1.0, 1000, "fig2b_omega_bias_L1000_w1.png", "relative error of median dN/dS", (-0.6, 0.6))
    # fraction with an estimate (dN and dS), pooled over freqs and kappa, omega 0.2
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for ax, L in zip(axes, [100, 300, 1000]):
        x = d[(d.omega == 0.2) & (d.len == L)].groupby(["method", "dS_true"]).frac_both.mean().reset_index()
        for m in METHODS:
            g = x[x.method == m]
            ax.plot(g.dS_true, g.frac_both, color=COL[m], marker=MRK[m], markersize=5, label=m)
        ax.set_xscale("log")
        ax.set_title(f"{L} codons", fontsize=9)
        ax.set_xlabel("true dS")
    axes[0].set_ylabel("fraction of pairs with dN and dS")
    axes[0].legend(ncol=2, fontsize=7.5, loc="lower left", frameon=False)
    fig.suptitle("Pairs with an estimate (omega 0.2; mean over both codon usages and kappa 1, 2, 4)", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{W}/plots/fig3_estimable.png", dpi=130)
    plt.close(fig)
    # length: relative RMSE of omega at moderate divergence
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for ax, dS in zip(axes, [0.1, 0.5, 1.0]):
        x = d[(d.omega == 0.2) & (d.dS_true == dS)].groupby(["method", "len"]).w_relrmse.median().reset_index()
        for m in METHODS:
            g = x[x.method == m]
            ax.plot(g.len, g.w_relrmse, color=COL[m], marker=MRK[m], markersize=5, label=m)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axhline(0.5, color="#52514e", lw=1)
        ax.set_title(f"true dS = {dS:g}", fontsize=9)
        ax.set_xlabel("codons")
        ax.set_xticks([30, 100, 300, 1000], ["30", "100", "300", "1000"])
    axes[0].set_ylabel("relative RMSE of dN/dS")
    axes[0].legend(ncol=2, fontsize=7.5, frameon=False)
    fig.suptitle("Gene length: relative RMSE of dN/dS (omega 0.2; median over codon usage x kappa; line: 50%)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{W}/plots/fig4_length_rmse.png", dpi=130)
    plt.close(fig)
    # SE coverage
    fig, axes = plt.subplots(2, 4, figsize=(12, 5.4), sharey=True, sharex=True)
    for i, c in enumerate(["dS", "dN"]):
        for j, L in enumerate([30, 100, 300, 1000]):
            ax = axes[i, j]
            x = d[(d.len == L) & d.method.isin(["YN", "MYN", "NG", "Li"])].groupby(["method", "dS_true"])[
                f"{c}_cov95"].mean().reset_index()
            for m in ["YN", "MYN", "NG", "Li"]:
                g = x[x.method == m]
                ax.plot(g.dS_true, g[f"{c}_cov95"], color=COL[m], marker=MRK[m], markersize=5, label=m)
            ax.axhline(0.95, color="#52514e", lw=1)
            ax.set_xscale("log")
            ax.set_ylim(0, 1.02)
            ax.set_title(f"{c}, {L} codons", fontsize=9)
            if i == 1:
                ax.set_xlabel("true dS")
        axes[i, 0].set_ylabel("coverage of est +- 1.96 SE")
    axes[0, 0].legend(fontsize=7.5, frameon=False, loc="lower left")
    fig.suptitle("Coverage of the nominal 95% interval from SE_dN / SE_dS (mean over omega, kappa, codon usage;"
                 " estimated pairs only)", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{W}/plots/fig5_se_coverage.png", dpi=130)
    plt.close(fig)
    print("figures written")


if __name__ == "__main__":
    main()
