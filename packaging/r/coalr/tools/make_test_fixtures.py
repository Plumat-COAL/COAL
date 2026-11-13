#!/usr/bin/env python3
"""Expected results for coalr's tests, written by COAL's command line (tests/testthat/fixtures).

    python packaging/r/coalr/tools/make_test_fixtures.py        # needs `coal` with coal-cpp, DIAMOND, Clustal Omega

Inputs: tests/fixtures/pipeline of the repository (synthetic sequences). The R tests compare the in-process
results of coalr with these tables, so they need neither the command line nor DIAMOND / Clustal Omega.
"""
from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SRC = ROOT / "tests" / "fixtures" / "pipeline"
DST = HERE.parent / "tests" / "testthat" / "fixtures"
sys.path.insert(0, str(ROOT / "src"))
from coal import dnds_pipeline as P   # noqa: E402
from coal.codon import fast_code      # noqa: E402

ALL = ["NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN"]


def fasta(path):
    out, rid = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if line.startswith(">"):
            rid = line[1:].split()[0]
            out[rid] = ""
        elif rid is not None:
            out[rid] += line
    return out


def main() -> int:
    DST.mkdir(parents=True, exist_ok=True)
    for name in ("query_cds.fasta", "subject_cds.fasta", "family_cds.fasta", "family_clusters.tsv"):
        shutil.copyfile(SRC / name, DST / name)
    code = fast_code(1)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        # pair mode: reciprocal best hits by DIAMOND, every method
        r = P.run_pairs(SRC / "query_cds.fasta", SRC / "subject_cds.fasta", tmp / "pairs", method=ALL, jobs=1,
                        quiet=True)
        assert r.exit_code in (0, 2), r.report.get("error")
        shutil.copyfile(r.outputs["orthologs"], DST / "pairs_orthologs.tsv")
        shutil.copyfile(r.outputs["dnds"], DST / "pairs_dnds.tsv")
        shutil.copyfile(r.outputs["genes"], DST / "pairs_genes.tsv")
        # the same pairs with pwalign's scoring
        r = P.run_pairs(SRC / "query_cds.fasta", SRC / "subject_cds.fasta", tmp / "pw", method=["YN"], jobs=1,
                        quiet=True, pairs=str(DST / "pairs_orthologs.tsv"), preset="pwalign")
        assert r.exit_code in (0, 2), r.report.get("error")
        shutil.copyfile(r.outputs["dnds"], DST / "pairs_pwalign_dnds.tsv")
        # cluster mode: Clustal Omega per cluster, all pairs
        r = P.run_clusters([SRC / "family_cds.fasta"], tmp / "cl", clusters=SRC / "family_clusters.tsv",
                           method=["YN", "NG"], jobs=1, quiet=True, keep_alignments=True)
        assert r.exit_code in (0, 2), r.report.get("error")
        shutil.copyfile(r.outputs["dnds"], DST / "clusters_dnds.tsv")
        shutil.copyfile(r.outputs["cluster_summary"], DST / "clusters_summary.tsv")
        with open(DST / "clusters_msa.tsv", "w", encoding="utf-8", newline="\n") as out:
            out.write("cluster_id\tid\tprotein_row\tcodon_row\n")
            adir = tmp / "cl" / "codon_alignments"
            for f in sorted(os.listdir(adir)):
                for rid, row in fasta(adir / f).items():
                    prot = "".join("-" if row[i:i + 3] == "---" else code.get(row[i:i + 3].upper().replace("U", "T"), "X")
                                   for i in range(0, len(row), 3))
                    out.write(f"{f[:-6]}\t{rid}\t{prot}\t{row}\n")
        # cluster mode: clusters above the size cap against their representative, pairwise
        r = P.run_clusters([SRC / "family_cds.fasta"], tmp / "rep", clusters=SRC / "family_clusters.tsv",
                           method=["YN"], jobs=1, quiet=True, max_cluster_size=5)
        assert r.exit_code in (0, 2), r.report.get("error")
        shutil.copyfile(r.outputs["dnds"], DST / "clusters_rep_dnds.tsv")
        shutil.copyfile(r.outputs["cluster_summary"], DST / "clusters_rep_summary.tsv")
        r = P.run_clusters([SRC / "family_cds.fasta"], tmp / "skip", clusters=SRC / "family_clusters.tsv",
                           method=["YN"], jobs=1, quiet=True, max_cluster_size=5, large_clusters="skip")
        shutil.copyfile(r.outputs["cluster_summary"], DST / "clusters_skip_summary.tsv")
    ver = subprocess.run([sys.executable, "-m", "coal", "--version"], capture_output=True, text=True,
                         env={**os.environ, "PYTHONPATH": str(ROOT / "src")}).stdout.strip()
    (DST / "README.txt").write_text(
        "Written by packaging/r/coalr/tools/make_test_fixtures.py with " + ver + " (coal pipeline).\n"
        "Inputs are copies of tests/fixtures/pipeline (synthetic sequences, no real data).\n", encoding="utf-8")
    print("wrote", DST)
    return 0


if __name__ == "__main__":
    sys.exit(main())
