# coalr

R bindings to [COAL](https://github.com/Plumat-COAL/COAL): codon alignment and dN/dS (Ka/Ks) in process, with
COAL's C++ core compiled into the package. No external programs, no temporary files.

```r
remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")

library(coalr)
res <- coal_dnds_pairs("query_cds.fasta", "subject_cds.fasta", pairs = ortholog_table, threads = 4)
cl  <- coal_dnds_clusters(cds, clusters)             # every pair of each cluster + a per-cluster summary
coal_kaks(codon_row_1, codon_row_2, method = "YN")   # KaKs_Calculator's numbers, without running it
coal_codon_align(protein_rows, cds)                  # PAL2NAL's job
```

Documentation: [docs/r.md](https://github.com/Plumat-COAL/COAL/blob/main/docs/r.md) and the help pages
(`?coal_dnds_pairs`). `src/coal/` is a copy of COAL's `cpp/src`, made by `scripts/sync_r_sources.py`; do not edit
it here. Licence: GPL-3 (version 3 only), like COAL; copyright and third-party notices are in `inst/COPYRIGHTS`.
