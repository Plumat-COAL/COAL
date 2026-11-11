# COAL from R: the `coalr` package

`coalr` (in `packaging/r/coalr`) gives R two ways to COAL's dN/dS pipeline ([pipeline.md](pipeline.md)):

| Route | Functions | What runs | Needs |
|---|---|---|---|
| **In process** (Rcpp) | `coal_dnds_pairs()`, `coal_dnds_clusters()`, `coal_kaks()`, `coal_codon_align()`, `coal_align_pairs()`, `coal_translate()` | COAL's C++ core compiled into the package: translation, protein alignment (NW, BLOSUM62), codon threading, Ka/Ks. No process is started and no file is written | a C++17 compiler at install time; `Rcpp` |
| **Through Python** (reticulate) | `coal_pipeline_pairs()`, `coal_pipeline_clusters()` | `coal.dnds_pipeline.run_pairs()` / `run_clusters()`: the whole `coal pipeline`, with the DIAMOND search or DeepClust and Clustal Omega; the tables it writes are read back as data frames | `reticulate`, a Python with `coal-align` (and its `coal-cpp`), DIAMOND, Clustal Omega for clusters |

**Which one.** When you already have the ortholog pairs or the clusters (from orthologr's `diamond_rec()`,
`deepclust()`, OrthoFinder, your own table), use the in-process functions: they take R objects, return data
frames, and need nothing outside R. When COAL should also find the orthologs or cluster the proteins, or you want
its report, `--resume`, or the files on disk, use the Python route (or the command line). Both give the same
numbers: the in-process functions compile the same C++ sources as `coal-cpp`.

## Install

```r
# install.packages("remotes")
remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")
```

or from a clone: `R CMD INSTALL packaging/r/coalr`. The package is self-contained: `src/coal/` holds a copy of
`cpp/src` (tables, core, kaks, align, dnds, pool), made by `python scripts/sync_r_sources.py` and committed
(`--check` verifies it; `tests/test_r_bindings.py` runs that check). It needs a C++17 compiler (Rtools on Windows,
Xcode command-line tools on macOS) and `Rcpp`. `reticulate`, `Biostrings` and `testthat` are optional.

For the Python route also:

```
pip install coal-align          # the platform wheels include coal-cpp; then DIAMOND: coal pipeline install diamond
```

```r
reticulate::use_virtualenv("~/venvs/coal")     # or use_python(), RETICULATE_PYTHON, an active virtualenv
coalr::coal_python_available()                 # TRUE, or FALSE with attr(, "reason")
```

`coal_python_available()` does not let reticulate download a Python of its own to answer; call
`reticulate::py_require("coal-align")` first if you want reticulate to manage one.

## Pairs

```r
library(coalr)
res <- coal_dnds_pairs("thaliana_cds.fasta", "lyrata_cds.fasta",
                       pairs = data.frame(query = c("AT1G01010.1", "AT1G01020.1"),
                                          subject = c("333554|PACid:16033839", "470181|PACid:16064328")),
                       methods = c("YN", "NG"), threads = 4)
```

- `query_cds`, `subject_cds`: a named character vector, a `DNAStringSet`, a `seqinr` list, a two-column data frame
  (ID, sequence: `orthologr::read.cds()`'s table works) or a FASTA path.
- `pairs`: two columns, query ID and subject ID (for example `orthologr::diamond_rec()`'s first two columns).
  Without it, the i-th query goes with the i-th subject.
- Result: a data frame with the columns of [`dnds.tsv`](pipeline.md#dndstsv): `seq1 seq2 method dN dS dNdS SE_dN
  SE_dS N_sites S_sites Nd Sd codons_compared prot_aln_length prot_identity status note`, rows in the order of
  `pairs`, `NA` for a value that is not available. A gene that cannot be used, an unknown ID or a pair that cannot
  be aligned is a row with status `error` and the reason, never an R error.
- Alignment options are those of the command line: `matrix`, `gap_open`, `gap_extend`, `end_gaps`, `max_cells`,
  `preset = "pwalign"`; `code` is the NCBI genetic code; `yn_flavour`, `yn_weighting` as in [kaks.md](kaks.md).
- `keep_alignments = TRUE` adds the protein and codon alignment of each pair as `attr(res, "alignments")`.

## Clusters

```r
cl <- orthologr::deepclust(...)             # or any table: representative, member
res <- coal_dnds_clusters(cds, cl[, 1:2], aligner = coal_msa_aligner("clustalo"), threads = 4)
res$pairs      # cluster_id + the pair columns, every pair of each cluster
res$summary    # one row per cluster and method: the columns of cluster_summary.tsv
res$failed     # members and clusters that could not be used
```

- `clusters`: a two-column table (representative, member) as `diamond deepclust` writes it, its file, or a named
  list of member IDs.
- Protein alignment of clusters of three or more: `msa =` rows you already have (a named list, cluster ID -> rows
  named by member), or `aligner =` a function that aligns one cluster's proteins. `coal_msa_aligner("clustalo")`
  and `coal_msa_aligner("mafft")` make such a function; they start the program once per cluster, the only place
  where the package starts a process. With neither, every pair is aligned on its own by the built-in aligner
  (`aligner` column `nw_pairwise`): no external program, but the values come from pairwise alignments.
- `max_cluster_size` (200), `large_clusters` (`"representative"`, `"sample"`, `"skip"`), `pair_sample`, `seed` as on
  the command line. The sample is reproducible and leaves R's random number stream alone, but it is not the draw
  the command line makes for the same seed (another generator).

## The building blocks

```r
coal_translate(cds)                                  # proteins (terminal stop dropped); details = TRUE: genes.tsv
coal_align_pairs(prot_a, prot_b)                     # global NW, BLOSUM62, full length; data frame of rows + score
coal_codon_align(protein_rows, cds)                  # PAL2NAL's job: protein alignment + DNA -> codon alignment
coal_kaks(codon_row_1, codon_row_2, method = "YN")   # Ka/Ks of aligned codon sequences (KaKs_Calculator's numbers)
coal_kaks(codon_msa, method = c("NG", "YN"))         # every pair of one alignment
coal_info()                                          # versions, methods, genetic codes, matrices
```

`coal_codon_align()` takes `nogap`, `nomismatch` and `strict` (PAL2NAL's `-nogap`, `-nomismatch`; COAL's
`--strict`) and matches DNA to protein rows by name when both are named.

## The Python route

```r
res <- coal_pipeline_pairs("q_cds.fa", "s_cds.fa", methods = "YN", threads = 8)       # DIAMOND RBH + everything
res$pairs; res$orthologs; res$unpaired; res$genes; res$failed; res$report$counts; res$exit_code
cl <- coal_pipeline_clusters("all_cds.fa", threads = 8, max_cluster_size = 100L, out_dir = "run1")
cl$summary
```

Further arguments are the options of `coal pipeline` with underscores (`ortho = "bh"`, `aligner = "diamond"`,
`pairs = "pairs.tsv"`, `codontable = 2L`, `install_tools = TRUE`, `resume = TRUE`, ...). Without `out_dir` the
run's folder is temporary and removed. Exit code 2 (finished, something is listed in `failed`) is returned, not
raised; a run that could not be done (exit code 1: an unreadable input, a missing tool) is an R error with the
pipeline's message. Without reticulate or a Python that has `coal-align`, both functions stop with a message that
says what to install; nothing else in the package needs Python.

## Threads, interrupts, errors

- `threads` (default `getOption("coalr.threads", 1L)`; `0` = one per hardware thread). Results are identical for
  every value.
- Worker threads only see C++ copies of the inputs. R's API is used on the calling thread only; R objects are
  built after the workers have been joined.
- A user interrupt (Ctrl-C, Esc) is noticed within about 0.1 s plus the pair each thread is on; the workers are
  stopped and joined before the interrupt is signalled.
- No C++ exception reaches R as a crash: bad data is a row's `status`/`note`; a bad argument is an R error.
- Memory: the inputs stay in R; the engine copies one batch (at most 1,024 groups, 16 Mb of CDS or 32,768 pairs)
  at a time, plus one alignment table per thread (`max_cells`).

## Validation

2026-11-11, R 4.5.3, Linux x86-64 (2-core cloud VM), gcc 16.2:

- `testthat`: 30 tests, 759 expectations. The package's fixtures are tables written by the command line
  (`packaging/r/coalr/tools/make_test_fixtures.py`, from `tests/fixtures/pipeline`): pair mode with all seven
  methods, the `pwalign` preset, cluster mode with Clustal Omega rows, clusters above the size cap (representative,
  skip), the summary, translation, the kept codon alignments. Every numeric column agrees to the 8 digits the
  command line prints; every text column is identical.
- Run against a live command line (`COALR_TEST_COAL=path/to/coal`): orthologr's 1,000-gene example, 999 reciprocal
  best hits, YN and NG: all 1,998 rows equal `coal pipeline pairs`. The reticulate route gives the same tables.
- `R CMD check --as-cran --no-manual`: no ERROR, no WARNING. NOTEs: the compile flag `-march=nocona` (conda's R
  adds it; not from the package) and "unable to verify current time" (no network).
- Time for the 999 pairs, YN, sequences read from FASTA: 1.5 s on one thread (0.3 s alignment + codon threading,
  1.1 s Ka/Ks); all seven methods 3.0 s.

Not tested: Windows and macOS (the C++ core builds there in COAL's CI, the R package was built on Linux only);
R older than 4.5; clusters or pair sets large enough to matter for memory.

## Maintenance

- After changing `cpp/src`: `python scripts/sync_r_sources.py`, commit the copies.
- The package version follows COAL's (`DESCRIPTION` = `pyproject.toml`; `tests/test_r_bindings.py` checks it).
- After changing the pipeline's output: `python packaging/r/coalr/tools/make_test_fixtures.py` (needs `coal` with
  `coal-cpp`, DIAMOND and Clustal Omega), then the R tests.
- `Rscript -e 'Rcpp::compileAttributes("packaging/r/coalr"); roxygen2::roxygenise("packaging/r/coalr")'` after
  changing exported C++ functions or the roxygen comments. Delete `src/*.o` afterwards: roxygen builds them
  without optimisation, and `R CMD INSTALL` would reuse them (ten times slower).
