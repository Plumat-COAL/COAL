"""Check a `coal pipeline clusters` run independently (docs/pipeline.md "Validation").

  1. Per-pair values: for every cluster with a kept codon alignment (RUN/codon_alignments/*.fasta), `coal-cpp dnds`
     on that file must give the rows of RUN/dnds.tsv (text of the estimates, SEs, sites, differences, codons).
  2. Clusters above --max-cluster-size (aligned pairwise, no file kept): up to --pairwise N of their pairs are
     redone by hand (translate in Python -> coal-cpp-align -> coal-cpp -> coal-cpp dnds).
  3. RUN/cluster_summary.tsv is recomputed from RUN/dnds.tsv with pandas / numpy (counts by status; median, mean,
     quartiles, min, max over the pairs with status ok or warning) and from the clusters file (members, singletons).

  python benchmarks/pipeline/validate_clusters.py --cds all_cds.fa [more.fa] --run RUN [--clusters clusters.tsv] \
      --coal-cpp build/coal-cpp --coal-cpp-align build/coal-cpp-align

Needs pandas. The run must have been made with --keep-alignments.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_pairs import read_fasta, read_fasta_list, translate  # noqa: E402

VALUE_COLS = ["dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "Nd", "Sd", "codons_compared"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cds", nargs="+", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--clusters", help="default: RUN/clusters.tsv (DeepClust's)")
    ap.add_argument("--coal-cpp", required=True)
    ap.add_argument("--coal-cpp-align", required=True)
    ap.add_argument("--pairwise", type=int, default=300)
    a = ap.parse_args()
    bad = 0
    cds = {}
    for p in a.cds:
        for k, v in read_fasta(p).items():
            cds.setdefault(k, v)
    cfile = a.clusters or os.path.join(a.run, "clusters.tsv")
    member_of, members = {}, {}
    for line in open(cfile).read().splitlines():
        rep, mem = line.split("\t")[:2]
        if mem not in member_of:
            member_of[mem] = rep
            members.setdefault(rep, []).append(mem)
    d = pd.read_csv(os.path.join(a.run, "dnds.tsv"), sep="\t", dtype=str, keep_default_na=False)
    methods = list(dict.fromkeys(d["method"]))
    ours = {(r.cluster_id, r.seq1, r.seq2, r.method): r for r in d.itertuples(index=False)}
    summ = pd.read_csv(os.path.join(a.run, "cluster_summary.tsv"), sep="\t", dtype=str, keep_default_na=False)
    margs = [x for m in methods for x in ("-m", m)]

    # 1. coal dnds on every kept codon MSA
    adir = os.path.join(a.run, "codon_alignments")
    n_files = n_rows = n_same = 0
    seen = set()
    with tempfile.TemporaryDirectory() as td:
        for name in sorted(os.listdir(adir)):
            path = os.path.join(adir, name)
            first = open(path).readline()[1:].split()[0]
            cid = member_of[first]
            out = os.path.join(td, "o.tsv")
            subprocess.run([a.coal_cpp, "dnds", path, "-o", out, "-q", *margs], check=True)
            n_files += 1
            for line in open(out).read().splitlines()[1:]:
                f = line.split("\t")
                o = ours.get((cid, f[0], f[1], f[2]))
                n_rows += 1
                seen.add((cid, f[0], f[1], f[2]))
                if o is not None and [getattr(o, c) for c in VALUE_COLS] == f[3:13] and \
                        (o.status == f[13] or (o.status == "warning" and "seq" in o.note)):
                    n_same += 1
                else:
                    bad += 1
                    print("differs:", cid, f, o)
        print(f"1. {n_files} codon MSAs, {n_rows} rows of `coal dnds` on them: {n_same} equal to dnds.tsv")
        rest = [k for k in ours if k not in seen]
        pw_clusters = set(summ.loc[summ["aligner"] == "nw_pairwise", "cluster_id"])
        stray = [k for k in rest if k[0] not in pw_clusters]
        if stray:
            bad += len(stray)
            print(f"   {len(stray)} rows of dnds.tsv belong to no kept alignment and no pairwise cluster:", stray[:3])

        # 2. pairwise clusters by hand
        todo = sorted({k[:3] for k in rest if k[0] in pw_clusters})[:a.pairwise]
        if todo:
            A, B = os.path.join(td, "a.faa"), os.path.join(td, "b.faa")
            with open(A, "w") as fa, open(B, "w") as fb:
                for _, x, y in todo:
                    fa.write(f">{x}\n{translate(cds[x])}\n")
                    fb.write(f">{y}\n{translate(cds[y])}\n")
            aln = os.path.join(td, "aln.fasta")
            subprocess.run([a.coal_cpp_align, A, B, "-o", aln], check=True)
            rows = read_fasta_list(aln)
            stacked = os.path.join(td, "stacked.fasta")
            with open(stacked, "w") as fs:
                for k, (_, x, y) in enumerate(todo):
                    pep, nuc = os.path.join(td, "pep.fa"), os.path.join(td, "nuc.fa")
                    open(pep, "w").write(f">a\n{rows[2 * k][1]}\n>b\n{rows[2 * k + 1][1]}\n")
                    open(nuc, "w").write(f">a\n{cds[x]}\n>b\n{cds[y]}\n")
                    r = subprocess.run([a.coal_cpp, pep, nuc, "-q"], capture_output=True, text=True, check=True)
                    r1, r2 = ("".join(p.splitlines()[1:]) for p in r.stdout.split(">")[1:])
                    fs.write(f">{x}\n{r1}\n>{y}\n{r2}\n")
            out = os.path.join(td, "o.tsv")
            subprocess.run([a.coal_cpp, "dnds", stacked, "--pairs", "consecutive", "-o", out, "-q", *margs],
                           check=True)
            same = n = 0
            for k, line in enumerate(open(out).read().splitlines()[1:]):
                f = line.split("\t")
                cid = todo[k // len(methods)][0]
                o = ours.get((cid, f[0], f[1], f[2]))
                n += 1
                if o is not None and [getattr(o, c) for c in VALUE_COLS] == f[3:13]:
                    same += 1
                else:
                    bad += 1
                    print("differs (pairwise):", cid, f, o)
            print(f"2. {len(todo)} pairs of clusters above the size cap redone by hand: {same} of {n} rows equal")

    # 3. the summary, with pandas
    num = d.copy()
    for c in ("dN", "dS", "dNdS"):
        num[c] = pd.to_numeric(num[c].replace("NA", np.nan))
    n_checked = n_cells = 0
    by = {k: g for k, g in num.groupby(["cluster_id", "method"], sort=False)}
    for r in summ.itertuples(index=False):
        g = by.get((r.cluster_id, r.method))
        mem = members.get(r.cluster_id, [])
        n_checked += 1
        problems = []
        if int(r.n_members) != len(mem):
            problems.append(f"n_members {r.n_members} != {len(mem)}")
        if g is None:
            if int(r.n_pairs) != 0 or r.dN_median != "NA":
                problems.append("no rows in dnds.tsv but n_pairs > 0")
            if len(mem) == 1 and r.pair_selection != "singleton":
                problems.append("singleton not marked")
        else:
            ok = g[g["status"].isin(["ok", "warning"])]
            want = {"n_pairs": len(g), "n_ok": len(ok), "n_saturated": int((g["status"] == "saturated").sum()),
                    "n_no_data": int((g["status"] == "no_data").sum()), "n_failed": int((g["status"] == "error").sum()),
                    "n_dNdS": int(ok["dNdS"].notna().sum())}
            for k, v in want.items():
                if int(getattr(r, k)) != v:
                    problems.append(f"{k} {getattr(r, k)} != {v}")
            for c in ("dN", "dS", "dNdS"):
                v = ok[c].dropna().to_numpy()
                stats = {"median": np.median, "mean": np.mean, "q25": lambda x: np.quantile(x, 0.25),
                         "q75": lambda x: np.quantile(x, 0.75), "min": np.min, "max": np.max}
                for s, fn in stats.items():
                    got = getattr(r, f"{c}_{s}")
                    n_cells += 1
                    if len(v) == 0:
                        if got != "NA":
                            problems.append(f"{c}_{s} {got} != NA")
                    elif got == "NA" or abs(float(got) - fn(v)) > 1e-7 * max(1.0, abs(fn(v))):
                        problems.append(f"{c}_{s} {got} != {fn(v)}")
            if r.pair_selection == "all" and len(g) != int(r.n_used) * (int(r.n_used) - 1) // 2:
                problems.append("pair count is not n(n-1)/2")
            if r.pair_selection == "representative" and len(g) != int(r.n_used) - 1:
                problems.append("pair count is not n-1")
        if problems:
            bad += 1
            print("summary differs:", r.cluster_id, r.method, problems)
    missing = set(members) - set(summ["cluster_id"])
    if missing:
        bad += len(missing)
        print("clusters missing from the summary:", sorted(missing)[:5])
    singles = summ[summ["pair_selection"] == "singleton"]
    print(f"3. {n_checked} summary rows ({len(set(summ['cluster_id']))} clusters, {len(singles)} singletons with "
          f"n_pairs 0), {n_cells} statistics recomputed with pandas/numpy: {'all equal' if not bad else 'see above'}")
    print("OK" if not bad else f"{bad} difference(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
