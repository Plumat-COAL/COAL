"""Aggregate COAL estimates against the truth. Writes res/cells.tsv (one row per file x run x method x dS)."""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

W = os.environ.get("COAL_SIM_DIR", "work/sim_truth")  # the work folder: data/, res/, plots/
truth = pd.read_csv(f"{W}/data/truth.tsv", sep="\t")
truth["key"] = truth.file.str.replace(".axt", "", regex=False)
TR = {(r.key, round(r.dS_label, 6)): r for r in truth.itertuples()}

DN_COLS = ["seq1", "method", "dN", "dS", "dNdS", "SE_dN", "SE_dS", "codons_compared", "status", "note"]


def load_dnds(path):
    d = pd.read_csv(path, sep="\t", usecols=DN_COLS, na_values=["NA"], keep_default_na=False)
    m = d.seq1.str.extract(r"d([0-9.]+)_r(\d+)")
    d["dS_label"] = m[0].astype(float).round(6)
    d["rep"] = m[1].astype(int)
    return d


def load_kaks(path):
    d = pd.read_csv(path, sep="\t", na_values=["NA", "nan", "-nan", "inf"], keep_default_na=False)
    d = d.rename(columns={"Sequence": "seq1", "Method": "method", "Ka": "dN", "Ks": "dS", "Ka/Ks": "dNdS"})
    rr = d["Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)"].astype(str).str.split(":", expand=True)
    d["kappa_hat"] = pd.to_numeric(rr[0], errors="coerce")
    d["kappa_hat2"] = pd.to_numeric(rr[1], errors="coerce")
    m = d.seq1.str.extract(r"d([0-9.]+)_r(\d+)")
    d["dS_label"] = m[0].astype(float).round(6)
    return d


def summarise(est, true):
    est = np.asarray(est, float)
    ok = np.isfinite(est)
    n = len(est)
    if ok.sum() == 0:
        return dict(n_est=0, frac_est=0.0, mean=np.nan, median=np.nan, bias=np.nan, relbias=np.nan,
                    medrel=np.nan, rmse=np.nan, relrmse=np.nan, mae_rel=np.nan)
    e = est[ok]
    return dict(n_est=int(ok.sum()), frac_est=ok.mean(), mean=e.mean(), median=np.median(e),
                bias=e.mean() - true, relbias=(e.mean() - true) / true, medrel=(np.median(e) - true) / true,
                rmse=np.sqrt(((e - true) ** 2).mean()), relrmse=np.sqrt(((e - true) ** 2).mean()) / true,
                mae_rel=np.median(np.abs(e - true)) / true)


def cell_rows(key, run, d):
    out = []
    for (meth, dl), g in d.groupby(["method", "dS_label"]):
        tr = TR[(key, dl)]
        row = dict(key=key, run=run, freqs=tr.freqs, kappa=tr.kappa, omega=tr.omega, len=tr.len, dS_true=tr.dS,
                   dN_true=tr.dN, method=meth, n=len(g))
        if "status" in g:
            st = g.status.value_counts()
            for s in ("ok", "warning", "saturated", "not_converged", "no_data", "error"):
                row["st_" + s] = st.get(s, 0) / len(g)
            row["warn_minsites"] = g.note.fillna("").str.contains("fewer than").mean()
        both = np.isfinite(g.dN) & np.isfinite(g.dS)
        row["frac_both"] = both.mean()
        for nm, col, tv in (("dS", "dS", tr.dS), ("dN", "dN", tr.dN)):
            for k, v in summarise(g[col], tv).items():
                row[f"{nm}_{k}"] = v
        w = g.dNdS.where(both)
        wt = tr.omega
        for k, v in summarise(w, wt).items():
            row[f"w_{k}"] = v
        # omega robust summaries: fraction within +-20% / 50% of truth among estimated
        we = w[np.isfinite(w)]
        row["w_within20"] = (np.abs(we - wt) <= 0.2 * wt).mean() if len(we) else np.nan
        row["w_within50"] = (np.abs(we - wt) <= 0.5 * wt).mean() if len(we) else np.nan
        row["w_q25"] = we.quantile(0.25) if len(we) else np.nan
        row["w_q75"] = we.quantile(0.75) if len(we) else np.nan
        # ratio of means estimator (pooled), used by some studies
        row["w_ratio_of_means"] = g.dN[both].mean() / g.dS[both].mean() if both.any() else np.nan
        if "SE_dN" in g:
            for nm, col, se, tv in (("dS", "dS", "SE_dS", tr.dS), ("dN", "dN", "SE_dN", tr.dN)):
                okk = np.isfinite(g[col]) & np.isfinite(g[se])
                row[f"{nm}_se_n"] = int(okk.sum())
                if okk.any():
                    e, s = g[col][okk], g[se][okk]
                    row[f"{nm}_cov95"] = (np.abs(e - tv) <= 1.96 * s).mean()
                    row[f"{nm}_se_mean"] = s.mean()
                    row[f"{nm}_sd_est"] = e.std()
                    row[f"{nm}_se_median"] = s.median()
        if "kappa_hat" in g:
            kh = g.kappa_hat[np.isfinite(g.kappa_hat) & both]
            row["kappa_med"] = kh.median() if len(kh) else np.nan
            row["kappa_mean"] = kh.mean() if len(kh) else np.nan
            kh2 = g.kappa_hat2[np.isfinite(g.kappa_hat2) & both]
            row["kappa2_med"] = kh2.median() if len(kh2) else np.nan
            row["kappa_inf"] = (~np.isfinite(g.kappa_hat2)).mean()
        out.append(row)
    return out


def main():
    rows = []
    files = sorted(glob.glob(f"{W}/res/coal/*.tsv.gz"))
    for p in files:
        b = os.path.basename(p)
        key, run = re.match(r"(.+)\.(\w+)\.tsv\.gz", b).groups()
        d = load_kaks(p) if run == "kaks" else load_dnds(p)
        if run in ("cgest", "cms0", "cms10"):
            d["method"] = d.method + "_" + run
        elif run == "liseq":
            d["method"] = "Li_seqinr"
        elif run == "ynpaml":
            d["method"] = "YN_paml"
        elif run == "ynw0":
            d["method"] = "YN_w0"
        elif run == "kaks":
            d["method"] = d.method + "_kaks"
        rows += cell_rows(key, run, d)
    df = pd.DataFrame(rows)
    df.to_csv(f"{W}/res/cells.tsv", sep="\t", index=False)
    print(len(files), "files", len(df), "rows")


if __name__ == "__main__":
    main()
