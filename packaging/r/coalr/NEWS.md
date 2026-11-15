# coalr 0.6.0 (in development: 0.6.0.9000)

First version. In-process bindings (Rcpp) to the C++ core of COAL 0.6:

- `coal_dnds_pairs()` without `pairs` still pairs the i-th query with the i-th subject, but now warns when both
  inputs have IDs and the i-th IDs differ: how many pairs differ (with an example) and how many of them have a
  protein identity below 20%, as unrelated genes do. Two CDS files are rarely in ortholog order, and such pairs
  used to get status `ok` without a word (the command line never pairs two CDS files by position). New value
  `pairs = "position"`: the same pairing, asked for, without the warning. No number changes.
- `coal_dnds_pairs()`: when every pair fails, a warning names the common cause: the columns of `pairs` are not
  the query and subject IDs (or are the other way round), or the sequences are protein, not coding DNA.
- A single string with a sequence-file extension in upper case (`"CDS.FASTA"`) that names no file is still read
  as a sequence, but the error it then leads to says that no such file exists in the working folder.
- `R CMD check --as-cran`: no error or warning (notes: the development version number; checks that need the
  network or HTML Tidy). DESCRIPTION cites the methods with their DOIs; the examples of `coal_msa_aligner()` run
  on the shipped data where Clustal Omega is installed, and those of `coal_pipeline_pairs()` (which need Python,
  `coal-cpp` and DIAMOND) use the shipped data; the tests use at most two threads unless `NOT_CRAN=true`.
- Comeron: `comeron_min_sites` applies only where it gives a positive distance; a 0 caused by a small class
  stays `NA` (`saturated`). Numbers: COAL's CHANGELOG (D9, Refined).
- `coal_dnds_clusters()`: the pairs of a cluster with `msa` rows or an `aligner` are computed in input order
  (the earlier gene of `cds` first), as `coal pipeline clusters` does since COAL's PR #71. Changes only Comeron
  with `comeron_twosubs = "libsequence"` (the default of `comeron_flavour = "gestimator"`).
- Ka/Ks engine fixes of COAL (vendored; numbers in COAL's CHANGELOG): MYN with an infinite kappa is `saturated`
  (it was `no_data`); a multiple-hit correction at an exact boundary is `saturated`, not a rounding residue
  printed as an estimate; the special cases of LWL, LPB, MLWL and MLPB follow the genetic code; YN/MYN rows
  without an estimate no longer have `Sd` / `Nd` that depend on the order of the sequences. Without
  COAL's Python pipeline, `coal_pipeline_pairs()` / `coal_pipeline_clusters()` of a development coalr name the
  build from a COAL checkout instead of the latest release.

- Comeron: new arguments `comeron_twosubs` and `comeron_min_sites` in `coal_kaks()`, `coal_dnds_pairs()` and
  `coal_dnds_clusters()` (the command line's `--comeron-twosubs`, `--comeron-min-sites`). By default Comeron
  now counts Comeron's two pathways for codons that differ at two positions (libsequence's branches made the
  estimate depend on which sequence is first; `comeron_twosubs = "libsequence"`, the default with
  `comeron_flavour = "gestimator"`, keeps them), and a site class with fewer than 3 sites can no longer make dN
  or dS `NA` (`comeron_min_sites = 0` restores that). Numbers: COAL's CHANGELOG.

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
