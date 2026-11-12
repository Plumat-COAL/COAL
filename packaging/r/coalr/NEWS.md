# coalr 0.5.0

First version. In-process bindings (Rcpp) to COAL 0.5's C++ core:

- `coal_dnds_pairs()`: CDS of ortholog pairs -> protein alignment -> codon alignment -> dN, dS, dN/dS; the columns
  of `coal pipeline pairs`' `dnds.tsv`.
- `coal_dnds_clusters()`: every pair of each cluster, with the per-cluster summary of `coal pipeline clusters`;
  protein alignments from given rows, an aligner function (`coal_msa_aligner()` for Clustal Omega / MAFFT) or the
  built-in pairwise aligner.
- `coal_kaks()`, `coal_codon_align()`, `coal_align_pairs()`, `coal_translate()`, `coal_read_fasta()`,
  `coal_info()`.
- `coal_pipeline_pairs()`, `coal_pipeline_clusters()`, `coal_python_available()`: the complete pipeline through
  Python (reticulate); optional.
- Licence: GPL-3, as COAL from this release (the Ka/Ks engine in `src/coal/kaks.cpp` is derived from
  KaKs_Calculator 2.0 and PAML's yn00; `inst/COPYRIGHTS`, `Authors@R`).
- One string that looks like a path but names no file is the error "file not found", not a sequence
  (`?coalr`, "Sequence inputs").
- Ka/Ks: YN and MYN estimate dS (dN) when there are no nonsynonymous (synonymous) sites; MYN no longer returns
  `saturated` for pairs whose TN93 term is 0/0 (COAL's docs/kaks.md, differences 8 to 10).
