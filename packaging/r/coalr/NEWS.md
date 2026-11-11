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
