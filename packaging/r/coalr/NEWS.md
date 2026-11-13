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
- One string that is written like a file name but names no file is the error "file not found", not a sequence.
  An aligned row is never taken for a file name because of what it ends in or contains (`"MKV..LFA.FA"`,
  `"MK/VL"`, `"~~MKV"`); `?coalr`, "Sequence inputs", has the rule and what it cannot tell apart.
- Ka/Ks: a YN or MYN iteration that does not converge gives status `not_converged` and no estimate (`dN`, `dS`,
  `dNdS`, the SEs, `Nd`, `Sd` are `NA`; `note` keeps the last round). The values of such a last round differ
  between platforms in the third digit (compiler, C library), so they are not a result. `coal_dnds_clusters()`'s
  summary counts these pairs in the new column `n_not_converged`, not in `n_ok`.
- The engine is compiled without floating-point contraction (no fused multiply-add) on every compiler.
- Ka/Ks: YN and MYN estimate dS (dN) when there are no nonsynonymous (synonymous) sites; MYN no longer returns
  `saturated` for pairs whose TN93 term is 0/0 (COAL's docs/kaks.md, differences 8 to 10).
