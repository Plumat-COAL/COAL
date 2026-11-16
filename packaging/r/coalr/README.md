# coalr

R bindings to [COAL](https://github.com/Plumat-COAL/COAL): codon alignment and dN/dS (Ka/Ks) in process, with
COAL's C++ core compiled into the package. No external programs, no temporary files.

```r
# a release (COAL 0.6.1 and later have coalr_<version>.tar.gz), or the same version from its tag
install.packages("Rcpp"); install.packages("https://github.com/Plumat-COAL/COAL/releases/download/v0.6.1/coalr_0.6.1.tar.gz", repos = NULL, type = "source")
remotes::install_github("Plumat-COAL/COAL@v0.6.1", subdir = "packaging/r/coalr")
# development version (main)
remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")

library(coalr)
res <- coal_dnds_pairs("query_cds.fasta", "subject_cds.fasta", pairs = ortholog_table, threads = 4)
cl  <- coal_dnds_clusters(cds, clusters)             # every pair of each cluster + a per-cluster summary
coal_kaks(codon_row_1, codon_row_2, method = "YN")   # KaKs_Calculator's numbers, without running it
coal_codon_align(protein_rows, cds)                  # PAL2NAL's job
```

Protein alignments are scored by default as `coal pipeline` and EMBOSS needle score them (BLOSUM62, open 10,
extend 0.5, free end gaps); `preset = "pwalign"` gives the alignments of `pwalign::pairwiseAlignment()`, which
orthologr uses. For some pairs the two differ and dS with them (by more than 10% for 118 of 2,958 pairs in COAL's
tutorial 4): use `preset = "pwalign"` to reproduce orthologr's numbers. `attr(res, "alignment")` says which scoring
a result used.

Documentation: [docs/r.md](https://github.com/Plumat-COAL/COAL/blob/main/docs/r.md) and the help pages
(`?coal_dnds_pairs`). `src/coal/` is a copy of COAL's `cpp/src`, made by `scripts/sync_r_sources.py`; do not edit
it here. Licence: GPL-3 (version 3 only), like COAL; copyright and third-party notices are in `inst/COPYRIGHTS`.
