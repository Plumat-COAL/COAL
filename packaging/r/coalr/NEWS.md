# coalr 0.6.0 (in development: 0.6.0.9000)

First version. In-process bindings (Rcpp) to the C++ core of COAL 0.6:

- `preset = "pwalign"` now returns the protein alignment `pwalign::pairwiseAlignment()` returns, not only one
  with the same score: pwalign's choice among co-optimal alignments (`tie_break = "pwalign"`, the preset's
  default). Before, 49% (A. thaliana) to 76% (yeast) of orthologr's protein alignments differed from COAL's
  under the same scoring; now 100% are identical on five genome pairs, and orthologr's `engine = "coal"`
  agrees with its native engine to print precision. New argument `tie_break` (`"leftmost"`, COAL's own
  choice and the default otherwise, or `"pwalign"`) in `coal_align_pairs()`, `coal_dnds_pairs()` and
  `coal_dnds_clusters()`.

- A vignette, "From two CDS files to dN/dS" (`vignette("from-two-cds-files-to-dnds")`), on the example that
  ships in `inst/extdata`; it runs in process, in seconds.
- `min_identity` in `coal_dnds_pairs()` and `coal_dnds_clusters()`: the third flag next to `min_codons` and
  `max_ds` (`coal pipeline --min-identity`), for pairs whose protein identity is below a percentage.
- `method` and `methods` name the same argument and both are accepted by `coal_dnds_pairs()`,
  `coal_dnds_clusters()`, `coal_kaks()`, `coal_pipeline_pairs()` and `coal_pipeline_clusters()` (a `method =`
  given to the pipeline functions used to be ignored). `?coalr` has a table of the command line's options and
  the R arguments; it also states the thread defaults (one thread in process, all CPUs through Python).

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
- Ka/Ks methods `"Comeron"` (alias `"C95"`; Comeron 1995 as libsequence's gestimator, the default method of the
  R package orthologr) and `"Li"` (Li 1993 as `seqinr::kaks()`), in `coal_kaks()`, `coal_dnds_pairs()` and
  `coal_dnds_clusters()`, for the standard genetic code. Arguments `comeron_flavour` and `li_flavour`: COAL's
  policy (default), or the reference program's reading of the pair (`"gestimator"`: the numbers of orthologr
  0.4.3; `"seqinr"`). `src/coal/kaks_comeron.cpp` is derived from libsequence (GPL >= 3; orthologr's copy) and
  `src/coal/kaks_li.cpp` from seqinr's `kaks.c` (GPL >= 2): `inst/COPYRIGHTS`, `Authors@R`. `"ALL"` still
  means KaKs_Calculator's seven methods.
- One string that is written like a file name but names no file is the error "file not found", not a sequence.
  An aligned row is never taken for a file name because of what it ends in or contains (`"MKV..LFA.FA"`,
  `"MK/VL"`, `"~~MKV"`); `?coalr`, "Sequence inputs", has the rule and what it cannot tell apart.
- Ka/Ks: a YN or MYN iteration that does not converge gives status `not_converged` and no estimate (`dN`, `dS`,
  `dNdS`, the SEs, `Nd`, `Sd` are `NA`; `note` keeps the last round). The values of such a last round differ
  between platforms in the third digit (compiler, C library), so they are not a result. `coal_dnds_clusters()`'s
  summary counts these pairs in the new column `n_not_converged`, not in `n_ok`.
- The engine is compiled without floating-point contraction (no fused multiply-add) on every compiler: `Comeron`
  and `Li` included, and the quartiles of `coal_dnds_clusters()`'s summary. A build for a CPU with FMA (arm64,
  `-march=x86-64-v3`) gives the same bytes as one without (COAL's docs/kaks.md, "Platform reproducibility").
- Ka/Ks: YN and MYN estimate dS (dN) when there are no nonsynonymous (synonymous) sites; MYN no longer returns
  `saturated` for pairs whose TN93 term is 0/0 (COAL's docs/kaks.md, differences 8 to 10).
