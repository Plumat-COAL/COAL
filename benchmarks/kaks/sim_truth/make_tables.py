"""Markdown tables for REPORT.md from res/cells.tsv, res/ref_*.tsv, res/ml_accuracy.tsv -> res/tables.md"""
import os
import numpy as np
import pandas as pd

W = os.environ.get("COAL_SIM_DIR", "work/sim_truth")  # the work folder: data/, res/, plots/
c = pd.read_csv(f"{W}/res/cells.tsv", sep="\t")
M9 = ["YN", "MYN", "NG", "LWL", "MLWL", "LPB", "MLPB", "Comeron", "Li"]
out = []


def md(df, floatfmt="{:.2f}"):
    cols = list(df.columns)
    s = "| " + " | ".join(str(x) for x in cols) + " |\n|" + "---|" * len(cols) + "\n"
    for _, r in df.iterrows():
        cells = []
        for x in r:
            if isinstance(x, (float, np.floating)):
                cells.append("" if np.isnan(x) else floatfmt.format(x))
            else:
                cells.append(str(x))
        s += "| " + " | ".join(cells) + " |\n"
    return s


def pct(x):
    return "" if pd.isna(x) else f"{100 * x:+.0f}%"


def h(t):
    out.append(f"\n### {t}\n")


main = c[c.run == "main"]

# A: range of median relative error of dS and omega across the 30 (freqs x kappa x omega) combos, L = 1000
for L in (1000, 300):
    for col, nm in (("dS_medrel", "dS"), ("w_medrel", "dN/dS")):
        h(f"A. Median relative error of {nm}, {L} codons: median [min, max] over 30 combinations (2 codon usages x 3 kappa x 5 omega)")
        rows = []
        for m in M9:
            r = {"method": m}
            for ds in (0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0):
                x = main[(main.method == m) & (main.len == L) & (main.dS_true == ds)][col]
                r[f"dS={ds:g}"] = f"{pct(x.median())} [{pct(x.min())}, {pct(x.max())}]" if x.notna().any() else "NA"
            rows.append(r)
        out.append(md(pd.DataFrame(rows)))

# B: showcase cells
for f, k, w, L in (("skewed", 4, 0.2, 300), ("equal", 1, 0.2, 1000), ("skewed", 2, 1.0, 300), ("equal", 2, 0.05, 1000)):
    h(f"B. {f} codon usage, kappa {k}, omega {w}, {L} codons: median relative error of dS / of dN/dS (fraction with both dN and dS)")
    x = main[(main.freqs == f) & (main.kappa == k) & (main.omega == w) & (main.len == L)]
    rows = []
    for m in M9:
        r = {"method": m}
        for ds in (0.03, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0):
            g = x[(x.method == m) & (x.dS_true == ds)]
            if len(g):
                g = g.iloc[0]
                r[f"dS={ds:g}"] = f"{pct(g.dS_medrel)} / {pct(g.w_medrel)} ({g.frac_both:.2f})"
        rows.append(r)
    out.append(md(pd.DataFrame(rows)))

# C: reliability onset
h("C. Saturation onset: largest true dS up to which every grid dS from 0.1 up has |median rel. error of dS| <= 15% "
  "and >= 90% of pairs with dN and dS [second number: also relative RMSE of dS <= 50%]: median (min..max) over the "
  "15 combinations of kappa x omega per codon usage; none = fails already at dS 0.1")
DSG = sorted(main.dS_true.unique())


def onset(g, rmse=None):
    g = g.set_index("dS_true").reindex(DSG)
    good = (g.dS_medrel.abs() <= 0.15) & (g.frac_both >= 0.9)
    if rmse is not None:
        good &= g.dS_relrmse <= rmse
    last = 0.0
    for ds in [x for x in DSG if x >= 0.1]:
        if bool(good.get(ds, False)):
            last = ds
        else:
            break
    return last


fmt = lambda v: "none" if v == 0 else f"{v:g}"
med = lambda s: fmt(float(np.percentile(s, 50, method="lower")))
rows = []
for m in M9:
    r = {"method": m}
    for f in ("equal", "skewed"):
        for L in (100, 300, 1000):
            x = main[(main.method == m) & (main.len == L) & (main.freqs == f)]
            on = x.groupby(["kappa", "omega"]).apply(onset, include_groups=False)
            on2 = x.groupby(["kappa", "omega"]).apply(lambda g: onset(g, 0.5), include_groups=False)
            r[f"{f}, {L}"] = f"{med(on)} ({fmt(on.min())}..{fmt(on.max())}) [{med(on2)}]"
    rows.append(r)
out.append(md(pd.DataFrame(rows)))

h("C2. Same criterion for dN/dS: |median rel. error of dN/dS| <= 15% and >= 90% of pairs estimable")


def onset_w(g):
    g = g.set_index("dS_true").reindex(DSG)
    good = (g.w_medrel.abs() <= 0.15) & (g.frac_both >= 0.9)
    first = [ds for ds in DSG if bool(good.get(ds, False))]
    if not first:
        return (np.nan, np.nan)
    lo = first[0]
    hi = lo
    for ds in DSG[DSG.index(lo):]:
        if bool(good.get(ds, False)):
            hi = ds
        else:
            break
    return (lo, hi)


rows = []
for m in M9:
    r = {"method": m}
    for L in (100, 300, 1000):
        x = main[(main.method == m) & (main.len == L)]
        on = x.groupby(["freqs", "kappa", "omega"]).apply(onset_w, include_groups=False)
        lo = pd.Series([v[0] for v in on]); hi = pd.Series([v[1] for v in on])
        r[f"{L}: combos with a reliable range"] = f"{lo.notna().sum()}/30"
        r[f"{L}: median range"] = f"{lo.median():g}-{hi.median():g}" if lo.notna().any() else "none"
    rows.append(r)
out.append(md(pd.DataFrame(rows)))

# E: gene length
for ds in (0.1, 0.5):
    h(f"E. Gene length at true dS {ds}, omega 0.2: relative RMSE of dN/dS / share of dN/dS within +-50% of truth / pairs estimable (median over 2 codon usages x 3 kappa)")
    rows = []
    for m in M9:
        r = {"method": m}
        for L in (30, 100, 300, 1000):
            x = main[(main.method == m) & (main.len == L) & (main.dS_true == ds) & (main.omega == 0.2)]
            r[f"{L} codons"] = f"{x.w_relrmse.median():.2f} / {x.w_within50.median():.2f} / {x.frac_both.median():.2f}"
        rows.append(r)
    out.append(md(pd.DataFrame(rows)))
h("E2. YN only, all omega: share of pairs whose dN/dS is within +-50% of the truth (median over codon usage x kappa), by length and true dS")
rows = []
for w in sorted(main.omega.unique()):
    for L in (30, 100, 300, 1000):
        r = {"omega": w, "codons": L}
        for ds in (0.01, 0.03, 0.1, 0.5, 1.0, 2.0):
            x = main[(main.method == "YN") & (main.len == L) & (main.dS_true == ds) & (main.omega == w)]
            r[f"dS={ds:g}"] = f"{(x.w_within50 * x.frac_both).median():.2f}"
        rows.append(r)
out.append(md(pd.DataFrame(rows)))

# F: SE coverage
h("F. Coverage of est +- 1.96 SE (dS / dN), estimated pairs, mean over the 30 combinations")
rows = []
for m in ("YN", "MYN", "NG", "Li"):
    for L in (30, 100, 300, 1000):
        r = {"method": m, "codons": L}
        for ds in (0.01, 0.1, 0.5, 1.0, 2.0):
            x = main[(main.method == m) & (main.len == L) & (main.dS_true == ds)]
            r[f"dS={ds:g}"] = f"{x.dS_cov95.mean():.2f} / {x.dN_cov95.mean():.2f}"
        rows.append(r)
out.append(md(pd.DataFrame(rows)))
h("F2. SE calibration: mean reported SE_dS / empirical SD of dS (estimated pairs), mean over combinations")
rows = []
for m in ("YN", "MYN", "NG", "Li"):
    for L in (100, 1000):
        r = {"method": m, "codons": L}
        for ds in (0.1, 0.5, 1.0, 1.5, 2.0):
            x = main[(main.method == m) & (main.len == L) & (main.dS_true == ds)]
            r[f"dS={ds:g}"] = f"{(x.dS_se_median / x.dS_sd_est).median():.2f}"
        rows.append(r)
out.append(md(pd.DataFrame(rows)))

# G: Comeron min-sites
h("G. Comeron: --comeron-min-sites 0 / 3 (default) / 10 and the gestimator flavour; fraction of pairs with dN, "
  "and median relative error of dN (omega 0.2, pooled over codon usage and kappa)")
rows = []
for L in (30, 100, 300, 1000):
    for ds in (0.1, 0.5, 1.0, 1.5, 2.0):
        r = {"codons": L, "true dS": ds}
        for m, lab in (("Comeron_cms0", "min 0"), ("Comeron", "min 3"), ("Comeron_cms10", "min 10"), ("Comeron_cgest", "gestimator")):
            x = c[(c.method == m) & (c.len == L) & (c.dS_true == ds) & (c.omega == 0.2)]
            r[lab] = f"{x.dN_frac_est.mean():.3f} / {pct(x.dN_medrel.median())} / RMSE {x.dN_relrmse.median():.2f}"
        rows.append(r)
out.append(md(pd.DataFrame(rows)))
h("G2. Comeron Ka (dN) NA although dN is small: fraction of pairs without dN, omega 0.05, by true dS (pooled)")
rows = []
for L in (100, 300, 1000):
    r = {"codons": L}
    for ds in (0.5, 1.0, 1.5, 2.0, 3.0):
        x = main[(main.method == "Comeron") & (main.len == L) & (main.dS_true == ds) & (main.omega == 0.05)]
        y = main[(main.method == "YN") & (main.len == L) & (main.dS_true == ds) & (main.omega == 0.05)]
        r[f"dS={ds:g} (dN={0.05*ds:g})"] = f"{1 - x.dN_frac_est.mean():.2f} (YN {1 - y.dN_frac_est.mean():.2f})"
    rows.append(r)
out.append(md(pd.DataFrame(rows)))

# H: flavours
h("H. YN flavours (kappa 2 only): median relative error of dS / dN/dS, pooled over codon usage and omega (median of cells)")
rows = []
for m in ("YN", "YN_paml", "YN_w0"):
    for L in (100, 1000):
        r = {"variant": m, "codons": L}
        for ds in (0.1, 0.5, 1.0, 1.5, 2.0):
            x = c[(c.method == m) & (c.len == L) & (c.dS_true == ds) & (c.kappa == 2)]
            r[f"dS={ds:g}"] = f"{pct(x.dS_medrel.median())} / {pct(x.w_medrel.median())} ({x.frac_both.mean():.2f})"
        rows.append(r)
out.append(md(pd.DataFrame(rows)))

# I: kappa
h("I. kappa estimated by YN (coal kaks, Substitution-Rate-Ratio): median over pairs, by true kappa, 1000 codons, omega 0.2")
rows = []
for f in ("equal", "skewed"):
    for k in (1, 2, 4):
        r = {"codon usage": f, "true kappa": k}
        for ds in (0.1, 0.5, 1.0, 2.0):
            x = c[(c.method == "YN_kaks") & (c.len == 1000) & (c.dS_true == ds) & (c.kappa == k) & (c.freqs == f) & (c.omega == 0.2)]
            r[f"dS={ds:g}"] = f"{x.kappa_med.median():.2f}"
        rows.append(r)
out.append(md(pd.DataFrame(rows)))

# J: references
try:
    a = pd.read_csv(f"{W}/res/ref_agreement.tsv", sep="\t")
    h("J. Agreement on the reference subset (8 files x 500 pairs): pairs where both give a number / within tolerance; max |difference|")
    g = a.groupby(["a", "b", "col"]).agg(both=("both", "sum"), within=("within_tol", "sum"), only_a=("only_a", "sum"),
                                         only_b=("only_b", "sum"), max_abs=("max_abs", "max"), tol=("tol", "first")).reset_index()
    out.append(md(g, "{:.2g}"))
    r = pd.read_csv(f"{W}/res/ref_accuracy.tsv", sep="\t")
    h("J2. Accuracy of every tool on the reference subset (kappa 2): median relative error of dS / dN/dS (fraction with dN/dS), pooled over codon usage, omega 0.2 and 1, 100 and 1000 codons (median of cells)")
    rows = []
    for m in sorted(r.method.unique()):
        row = {"tool_method": m}
        for ds in (0.1, 0.5, 1.0, 1.5, 2.0, 3.0):
            x = r[(r.method == m) & (r.dS_label == ds)]
            row[f"dS={ds:g}"] = f"{pct(x.dS_medrel.median())} / {pct(x.dNdS_medrel.median())} ({x.dNdS_frac_est.mean():.2f})"
        rows.append(row)
    out.append(md(pd.DataFrame(rows)))
    ml = pd.read_csv(f"{W}/res/ml_accuracy.tsv", sep="\t")
    h("J3. codeml ML (GY94 F3x4, pairwise) vs COAL methods, 1000 codons, omega 0.2, 20 pairs per cell: median relative error of dS / dN/dS")
    rows = []
    for m in ["codeml_ML", "yn00_YN", "yn00_NG86", "COAL_YN", "COAL_MYN", "COAL_NG", "COAL_LPB", "COAL_Comeron", "COAL_Li"]:
        for f in ("equal", "skewed"):
            for k in (1, 2, 4):
                row = {"method": m, "usage": f, "kappa": k}
                for ds in (0.1, 0.5, 1.0, 1.5, 2.0, 3.0):
                    x = ml[(ml.method == m) & (ml.dS_label == ds) & (ml.freqs == f) & (ml.kappa == k)]
                    row[f"dS={ds:g}"] = f"{pct(x.dS_medrel.median())} / {pct(x.dNdS_medrel.median())}" if len(x) else ""
                rows.append(row)
    out.append(md(pd.DataFrame(rows)))
except FileNotFoundError as e:
    out.append(f"(missing: {e})\n")

open(f"{W}/res/tables.md", "w").write("".join(out))
print("tables written")

# D (per pair)
pp = pd.read_csv(f"{W}/res/perpair.tsv", sep="\t")
for L in (100, 1000):
    out.append(f"\n### D. {L} codons: share of pairs flagged (status not ok/warning) / share with dN and dS whose dS is >30% off the truth (pooled over the 30 combinations)\n")
    rows = []
    for m in M9:
        r = {"method": m}
        for ds in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
            x = pp[(pp.method == m) & (pp.len == L) & (pp.dS_true == ds)]
            r[f"dS={ds:g}"] = f"{(x.flagged * x.n).sum() / x.n.sum():.2f} / {x.ok_off30.sum() / x.n.sum():.2f}"
        rows.append(r)
    out.append(md(pd.DataFrame(rows)))
pf = pd.read_csv(f"{W}/res/perpair_filter.tsv", sep="\t")
out.append("\n### D2. The doc's filter 'dS at most X': of the pairs that pass it (estimate exists, estimated dS <= X), share whose TRUE dS is above X / above 1.5 X (pooled over everything; the grid has equal numbers of pairs at true dS 0.01-3)\n")
rows = []
for thr in (1.0, 2.0):
    for m in M9:
        r = {"X": thr, "method": m}
        for L in (100, 300, 1000):
            x = pf[(pf.method == m) & (pf.len == L) & (pf.thr == thr)]
            r[f"{L} codons"] = f"{x.true_above.sum() / x.kept.sum():.3f} / {x.true_above_15thr.sum() / x.kept.sum():.3f}"
        rows.append(r)
out.append(md(pd.DataFrame(rows)))
x = pp[pp.method == "MYN"]
out.append(f"\nMYN rows with status no_data and the note 'no synonymous sites; no nonsynonymous sites' (all simulated pairs have sites; cause: kappa estimated as infinite): {int(x.kinf.sum())} of {int(x.n.sum())}\n")
x = pp[pp.method == "YN"]
out.append(f"\nYN rows with that note: {int(x.kinf.sum())} of {int(x.n.sum())}\n")
x = pp[pp.method == "Comeron"]
out.append(f"\nComeron rows with the min-sites note: {int(x.minsites.sum())} of {int(x.n.sum())}; by length: " +
           ", ".join(f"{L}: {int(x[x.len == L].minsites.sum())}" for L in (30, 100, 300, 1000)) + "\n")
open(f"{W}/res/tables.md", "w").write("".join(out))
print("tables + D written")
