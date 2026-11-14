#!/usr/bin/env python3
"""Filter and summarise a COAL dnds.tsv (coal pipeline pairs / clusters, coal dnds, coalr tables written out).

Used by the COAL tutorials (docs/tutorials/). Standard library only.

  python3 dnds_summary.py out/dnds.tsv [--cds a.fna.gz b.fna.gz] [--method YN]
          [--min-codons 50] [--min-coverage 0.5] [--min-identity 40] [--min-ds 0.01] [--max-ds 2]
          [--top 10] [--write-kept kept.tsv] [--only-ids ids.txt]

It prints: the rows per status; how many rows each filter removes, in order; the quartiles of dN, dS and dN/dS
of the rows that pass; and the rows with the highest and lowest dN/dS among those, with gene names read from
the [gene=...] / [protein=...] fields of NCBI cds_from_genomic headers when --cds is given. No number in the
table is changed; the thresholds are the starting points of docs/kaks.md, "Filtering before interpretation".
"""
import argparse
import csv
import gzip
import re
import statistics
import sys


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_names(paths):
    names = {}
    for path in paths or []:
        with open_text(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    rid = line[1:].split()[0]
                    gene = re.search(r"\[gene=([^\]]+)\]", line)
                    prot = re.search(r"\[protein=([^\]]+)\]", line)
                    tag = re.search(r"\[locus_tag=([^\]]+)\]", line)
                    label = gene.group(1) if gene else (tag.group(1) if tag else "")
                    names[rid] = (label, prot.group(1) if prot else "")
    return names


def num(x):
    return None if x in ("NA", "") else float(x)


def quartiles(v):
    v = sorted(v)
    if len(v) < 2:
        return (v[0],) * 3 if v else (None,) * 3
    q = statistics.quantiles(v, n=4, method="inclusive")  # R type 7
    return q[0], q[1], q[2]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dnds")
    ap.add_argument("--cds", nargs="*", help="the CDS files, for gene names (NCBI headers)")
    ap.add_argument("--method", default=None, help="only rows of this method (default: the first method)")
    ap.add_argument("--min-codons", type=int, default=50)
    ap.add_argument("--min-coverage", type=float, default=0.5, help="codons_compared / prot_aln_length")
    ap.add_argument("--min-identity", type=float, default=40)
    ap.add_argument("--min-ds", type=float, default=0.01)
    ap.add_argument("--max-ds", type=float, default=2)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--write-kept", metavar="TSV", help="write the rows that pass every filter")
    ap.add_argument("--only-ids", metavar="FILE", help="only rows whose seq1 or seq2 is one of the IDs in FILE (one per line)")
    a = ap.parse_args()

    with open(a.dnds, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        rows = list(rd)
        fields = rd.fieldnames
    method = a.method or (rows[0]["method"] if rows else "YN")
    rows = [r for r in rows if r["method"] == method]
    if a.only_ids:
        with open(a.only_ids) as fh:
            ids = {line.strip() for line in fh if line.strip()}
        rows = [r for r in rows if r["seq1"] in ids or r["seq2"] in ids]
    print(f"{a.dnds}: {len(rows)} rows of method {method}" + (f" (only the IDs in {a.only_ids})" if a.only_ids else ""))

    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    print("status: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1])))

    has_ident = "prot_identity" in fields
    steps = [
        ("has an estimate (status ok or warning)", lambda r: r["status"] in ("ok", "warning")),
        (f"codons_compared >= {a.min_codons}", lambda r: int(r["codons_compared"]) >= a.min_codons),
    ]
    if has_ident:
        steps += [
            (f"codons_compared / prot_aln_length >= {a.min_coverage}",
             lambda r: int(r["codons_compared"]) >= a.min_coverage * int(r["prot_aln_length"])),
            (f"prot_identity >= {a.min_identity}", lambda r: float(r["prot_identity"]) >= a.min_identity),
        ]
    steps += [
        (f"dS <= {a.max_ds}", lambda r: num(r["dS"]) is not None and num(r["dS"]) <= a.max_ds),
        (f"dS >= {a.min_ds}", lambda r: num(r["dS"]) >= a.min_ds),
    ]
    kept = rows
    print("filters (each applied to what the previous ones kept):")
    for label, f in steps:
        before = len(kept)
        kept = [r for r in kept if f(r)]
        print(f"  {label:<45} {len(kept):>7}  (-{before - len(kept)})")

    if not kept:
        return
    print(f"quartiles over the {len(kept)} rows kept (25% / median / 75%):")
    for col in ("dN", "dS", "dNdS"):
        v = [num(r[col]) for r in kept if num(r[col]) is not None]
        q = quartiles(v)
        print(f"  {col:<5} {q[0]:.4f} / {q[1]:.4f} / {q[2]:.4f}")
    above1 = sum(1 for r in kept if num(r["dNdS"]) is not None and num(r["dNdS"]) > 1)
    print(f"  dN/dS > 1: {above1} of {len(kept)}")

    names = read_names(a.cds)
    def show(r):
        g1 = names.get(r["seq1"], ("", ""))
        g2 = names.get(r["seq2"], ("", ""))
        label = g1[0] if g1[0] == g2[0] or not g2[0] else f"{g1[0]}/{g2[0]}"
        cl = (r["cluster_id"] + "  ") if "cluster_id" in r else ""
        print(f"  {cl}{label:<10} dN {float(r['dN']):.4f}  dS {float(r['dS']):.4f}  dN/dS {float(r['dNdS']):.4f}"
              f"  codons {r['codons_compared']:>4}  {g1[1][:45]}")
    ranked = sorted((r for r in kept if num(r["dNdS"]) is not None), key=lambda r: float(r["dNdS"]))
    if a.top and ranked:
        print(f"highest dN/dS among the rows kept:")
        for r in ranked[::-1][: a.top]:
            show(r)
        print(f"lowest dN/dS among the rows kept:")
        for r in ranked[: a.top]:
            show(r)
    if a.write_kept:
        with open(a.write_kept, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(kept)
        print(f"wrote {len(kept)} rows to {a.write_kept}")


if __name__ == "__main__":
    sys.exit(main())
