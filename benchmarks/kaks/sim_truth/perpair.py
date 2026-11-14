"""Per-pair statistics over the main runs: flags vs truth, filter leakage. -> res/perpair.tsv, res/perpair_filter.tsv"""
import glob, os, re
import numpy as np, pandas as pd
W = os.environ.get("COAL_SIM_DIR", "work/sim_truth")  # the work folder: data/, res/, plots/
truth = pd.read_csv(f"{W}/data/truth.tsv", sep="\t")
truth["key"] = truth.file.str.replace(".axt", "", regex=False)
truth["dS_label"] = truth.dS_label.round(6)
rows, filt, notes = [], [], []
for p in sorted(glob.glob(f"{W}/res/coal/*.main.tsv.gz")):
    key = os.path.basename(p).split(".main")[0]
    d = pd.read_csv(p, sep="\t", na_values=["NA"], keep_default_na=False, usecols=["seq1", "method", "dN", "dS", "dNdS", "codons_compared", "status", "note"])
    d["dS_label"] = d.seq1.str.extract(r"d([0-9.]+)_")[0].astype(float).round(6)
    t = truth[truth.key == key][["dS_label", "dS", "dN", "omega", "len", "freqs", "kappa"]].rename(columns={"dS": "dS_true", "dN": "dN_true"})
    d = d.merge(t, on="dS_label")
    d["has"] = np.isfinite(d.dS) & np.isfinite(d.dN)
    d["okst"] = d.status.isin(["ok", "warning"])
    d["off30"] = d.has & (np.abs(d.dS - d.dS_true) > 0.3 * d.dS_true)
    d["woff50"] = d.has & np.isfinite(d.dNdS) & (np.abs(d.dNdS - d.omega) > 0.5 * d.omega)
    g = d.groupby(["method", "len", "dS_true"]).agg(n=("has", "size"), flagged=("okst", lambda x: 1 - x.mean()),
        has=("has", "mean"), ok_off30=("off30", "sum"), ok_woff50=("woff50", "sum"),
        kinf=("note", lambda x: x.fillna("").str.contains("no synonymous sites; no nonsynonymous").sum()),
        minsites=("note", lambda x: x.fillna("").str.contains("fewer than 3 sites").sum())).reset_index()
    g["key"] = key
    rows.append(g)
    # filter leakage: pairs with estimate and estimated dS <= 2: distribution of true dS
    e = d[d.has]
    for thr in (1.0, 2.0):
        k = e[e.dS <= thr]
        f2 = k.groupby(["method", "len"]).agg(kept=("dS", "size"), true_above=("dS_true", lambda x: (x > thr).sum()),
            true_above_15thr=("dS_true", lambda x: (x > 1.5 * thr).sum())).reset_index()
        f2["thr"] = thr; f2["key"] = key
        filt.append(f2)
pd.concat(rows).to_csv(f"{W}/res/perpair.tsv", sep="\t", index=False)
pd.concat(filt).to_csv(f"{W}/res/perpair_filter.tsv", sep="\t", index=False)
print("ok")
