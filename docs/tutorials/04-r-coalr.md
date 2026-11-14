# Tutorial 4: the same from R (coalr, and orthologr's `engine = "coal"`)

**Goal:** repeat tutorials 2 and 3 from R, with the `coalr` package (in process, or the whole pipeline through
Python), and with the fork of orthologr that has COAL built in; and check that every route gives the same numbers.

| | |
|---|---|
| Data | the files of tutorials 2 and 3 (`get_enterobacteria.sh`, `get_chlamydia.sh`) and their `orthologs.tsv` / `clusters.tsv` |
| Needs | R 4.1 or later with a C++17 compiler and `Rcpp`; `reticulate` and the COAL platform wheel for the Python route; Clustal Omega for clusters; DIAMOND and `remotes` for orthologr |
| Runtime | in-process pairs 7 s; the pipeline through Python 14 s; in-process clusters with Clustal Omega 71 s; orthologr 20 s (one thread each) |
| Tested with | R 4.5.3, coalr 0.6.0.9000 (COAL commit 009bdd5), orthologr 0.5.0 of Plumat-COAL/orthologr `integration/all` (df3ea87), DIAMOND 2.2.8, Clustal Omega 1.2.4, Linux x86_64 |

Run the R code from the scratch folder of the earlier tutorials (`coal-tutorials/`), after tutorials 2 and 3 have
written `ecoli_vs_salmonella/` and `chlamydia/` (skip their clean-up steps, or run them again).

## 0. Install

```r
install.packages(c("remotes", "Rcpp"))
remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")
```

That is all the in-process functions need: COAL's C++ core is compiled into the package, so no Python, no
`coal-cpp` and no external program (apart from Clustal Omega or MAFFT, if you want them for clusters).

**orthologr with COAL built in** is a fork, [Plumat-COAL/orthologr](https://github.com/Plumat-COAL/orthologr),
not the upstream [drostlab/orthologr](https://github.com/drostlab/orthologr); the COAL changes have not been
merged there. Its branch `integration/all` has `dNdS(engine = "coal")`:

```r
remotes::install_github("Plumat-COAL/orthologr", ref = "integration/all")
```

(orthologr itself needs Biostrings and its other dependencies; see its README.)

## 1. Pairs, in process: `coal_dnds_pairs()`

Given the sequences and the pairs, `coal_dnds_pairs()` aligns, threads and estimates inside R. Here the pairs are
the DIAMOND reciprocal best hits that `coal pipeline pairs` found in tutorial 2:

```r
library(coalr)
q <- coal_read_fasta("coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz")
s <- coal_read_fasta("coal-tutorial-data/styphimurium_LT2.cds.fna.gz")
pairs <- read.delim("ecoli_vs_salmonella/orthologs.tsv")
res <- coal_dnds_pairs(q, s, pairs = pairs[, 1:2], methods = "YN", threads = 1)
dim(res)
table(res$status)
head(res[, c("dN", "dS", "dNdS", "codons_compared", "prot_identity", "status")], 3)
```

```
[1] 3142   17

not_converged            ok     saturated       warning 
           11          2973           155             3 
          dN       dS       dNdS codons_compared prot_identity status
1 0.03146997 1.909781 0.01647831             820         94.51     ok
2 0.03568057 2.016680 0.01769273             309         92.90     ok
3 0.03896516 1.238365 0.03146501             428         93.46     ok
```

The columns are those of `dnds.tsv` (`seq1 seq2 method dN dS dNdS SE_dN SE_dS N_sites S_sites Nd Sd
codons_compared prot_aln_length prot_identity status note`), and so are the numbers:

```r
cli <- read.delim("ecoli_vs_salmonella/dnds.tsv")
identical(res$status, cli$status)
identical(signif(res$dS, 8), signif(cli$dS, 8))    # the TSV has 8 significant digits
```

```
[1] TRUE
[1] TRUE
```

The filters of tutorial 2, in R (the line in [kaks.md](../kaks.md#filtering-before-interpretation)):

```r
keep <- res$status %in% c("ok", "warning") & res$codons_compared >= 50 &
        res$codons_compared / res$prot_aln_length >= 0.5 & res$prot_identity >= 40 &
        !is.na(res$dS) & res$dS >= 0.01 & res$dS <= 2
sum(keep)
quantile(res$dS[keep], c(0.25, 0.5, 0.75))
summary(res$dNdS[keep])
```

```
[1] 2468
      25%       50%       75% 
0.8539104 1.2227868 1.4945717 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00000 0.02578 0.04083 0.05720 0.06682 3.96510 
```

The same 2,468 pairs and quartiles as the summary script of tutorial 2. The pairs can come from anywhere: a
two-column data frame of query and subject IDs (OrthoFinder, `orthologr::diamond_rec()`, your own table). The
inputs can be FASTA paths, named character vectors, `DNAStringSet`s or seqinr lists. More options
(`min_codons`, `max_ds`, `min_identity`, `keep_alignments`, `code`, the methods): `?coal_dnds_pairs` and
[r.md](../r.md).

## 2. The whole pipeline from R: `coal_pipeline_pairs()`

When COAL should also find the orthologs, call the pipeline through reticulate. This needs a Python with the
COAL platform wheel (which has `coal-cpp`) and DIAMOND, as for tutorial 2:

```r
library(coalr)
reticulate::use_virtualenv("~/venvs/coal", required = TRUE)    # the Python that has coal-align
coal_python_available()
run <- coal_pipeline_pairs("coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz",
                           "coal-tutorial-data/styphimurium_LT2.cds.fna.gz",
                           threads = 1, out_dir = "ecsa_from_r")
names(run)
table(run$pairs$status)
run$report$entry
```

```
[1] TRUE
[1] "exit_code" "pairs"     "orthologs" "unpaired"  "genes"     "failed"   
[7] "report"    "out_dir"  

not_converged            ok     saturated       warning 
           11          2973           155             3 
[1] "coalr 0.6.0.9000 coal_pipeline_pairs()"
```

It is `coal pipeline pairs`, run from R: `ecsa_from_r/` has the same files, and its `dnds.tsv` is byte-identical
to tutorial 2's (`sha256sum`: `7ed75371...e1b8e8` for both). The report records the R entry point and the
equivalent command line (`run$report$command`). Without `out_dir` the folder is temporary.

## 3. Clusters, in process: `coal_dnds_clusters()`

Tutorial 3's clusters, with one Clustal Omega alignment per cluster:

```r
library(coalr)
files <- file.path("coal-tutorial-data",
                   c("ctrachomatis_D", "ctrachomatis_L2", "cmuridarum_Nigg", "cpneumoniae_CWL029"))
cds <- unlist(lapply(paste0(files, ".cds.fna.gz"), coal_read_fasta))
cl  <- read.delim("chlamydia/clusters.tsv", header = FALSE)    # representative, member
res <- coal_dnds_clusters(cds, cl, aligner = coal_msa_aligner("clustalo"), threads = 1)
names(res)
table(res$pairs$status)
head(res$summary[order(-res$summary$n_members),
                 c("n_members", "n_pairs", "n_ok", "n_saturated", "dS_median", "dNdS_median", "status")], 3)
```

```
[1] "pairs"   "summary" "failed" 

      no_data not_converged            ok     saturated       warning 
           58           105          3386          3290            34 
    n_members n_pairs n_ok n_saturated dS_median dNdS_median status
623        40     780  103         672 2.7301650   0.2895782     ok
382        16     120   43          77 2.0384366   0.3843151     ok
292        15     105   71           7 0.7758856   0.5640490     ok
```

`res$pairs` and `res$summary` are `dnds.tsv` and `cluster_summary.tsv` of tutorial 3 (same rows, statuses and
values: `identical(res$summary$status, read.delim("chlamydia/cluster_summary.tsv")$status)` is `TRUE`). The first
row is the 40-member Pmp family of tutorial 3. Without `aligner =`, every pair is aligned on its own by the
built-in aligner (no external program, other numbers); `coal_msa_aligner("mafft")` uses MAFFT. A cluster table
from `orthologr::deepclust()` works the same way.

## 4. orthologr: `dNdS(engine = "coal")`

orthologr's `dNdS()` does the search itself (DIAMOND reciprocal best hits) and, with `engine = "coal"`, hands
all pairs to coalr at once. orthologr's default search is `sensitivity_mode = "fast"`; `"default"` matches
`coal pipeline`:

```r
library(orthologr)
d <- dNdS("ecoli.fna", "salmonella.fna", dnds_est.method = "YN", engine = "coal",
          sensitivity_mode = "default", comp_cores = 1)
```

(`ecoli.fna` and `salmonella.fna` are the two files of tutorial 2, uncompressed:
`zcat coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz > ecoli.fna`, the same for *Salmonella*.)

```
dNdS_coal(): reciprocal best hits of 'ecoli.fna' against 'salmonella.fna' (diamond, sensitivity_mode 'default', e-value 1E-5, best_hit_rule 'bitscore'), then dN/dS by method 'YN' with the COAL engine (aa_aln_scoring 'pwalign', BLOSUM62, comp_cores = 1).
dNdS_coal(): 3143 gene pairs - 2969 with dN/dS, 4 with dS = 0 (no ratio), 170 saturated (dN or dS not estimable). Median dS: 1.33. 13 not converged (see the notes of the engine).
1 coding sequence has a length that is not a multiple of 3 and was trimmed at the 3' end (delete_corrupt_cds = FALSE; ...). They are listed in attr(<result>, "trimmed_sequences"); ...
```

The result is orthologr's table (`query_id subject_id dN dS dNdS method perc_identity ... status`, 3,143 rows).
**It is not tutorial 2's table, for two documented reasons**, and both can be made to match:

| Difference | Here | Why | To match |
|---|---|---|---|
| the pairs | 3,143 against 3,142: 2 pairs only in orthologr, 1 only in COAL | orthologr takes the best hit by bit score, COAL takes DIAMOND's first hit (lowest e-value) | give orthologr's pairs to COAL (`--pairs`), or COAL's to coalr (section 1) |
| the protein alignment | for 118 of the 2,958 pairs with dS in both, dS differs by more than 10% | the fork aligns as `pwalign::pairwiseAlignment()` does (`aa_aln_scoring = "pwalign"`: end gaps penalised, open 14, extend 4); `coal pipeline` uses free end gaps, open 10, extend 0.5 | `coal pipeline pairs ... --preset pwalign` |

With `--preset pwalign`, the dS of all 2,972 pairs that both routes estimate agree to the 8 significant digits of
the TSV (largest difference 5e-8):

```
coal pipeline pairs --query coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz \
    --subject coal-tutorial-data/styphimurium_LT2.cds.fna.gz -o ecsa_pwalign -j 1 --preset pwalign -q
```

The alignment choice matters more than one would guess: *thrB* (the second row above) has the same 309 compared
codons and the same 92.9% identity under both alignments, but dS 2.02 with free end gaps and 3.83 with pwalign's,
because one gap sits elsewhere and this pair is close to saturation. Neither is "the" alignment; choose one,
keep it for everything you compare, and say which ([pipeline.md, "How much the aligner
matters"](../pipeline.md#how-much-the-aligner-matters)). orthologr's `"Comeron"` (its default method) and `"Li"`
are in COAL's engine too (`methods = "Comeron"` in coalr, with `comeron_flavour = "gestimator"` for orthologr
0.4.3's numbers; [kaks.md](../kaks.md#comeron-and-li)).

## Which route

| You have | Use |
|---|---|
| pairs or clusters already (orthologr, OrthoFinder, your own) | `coal_dnds_pairs()` / `coal_dnds_clusters()`: in process, nothing outside R |
| two CDS files, and want COAL's search, report and files | `coal_pipeline_pairs()` (or the command line) |
| an orthologr workflow | `dNdS(engine = "coal")` of the fork, knowing the two differences above |

## Clean up

```
rm -rf ecsa_from_r ecsa_pwalign ecoli.fna salmonella.fna
```

## What to cite

COAL and what its run used (`coal cite`; for the in-process functions, the method and PAL2NAL); for orthologr,
also Drost et al. 2015, which `dNdS()` prints once per session.
