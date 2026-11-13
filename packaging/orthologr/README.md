# COAL in orthologr

[orthologr](https://github.com/drostlab/orthologr) (`dNdS()`, `codon_aln()`) builds one codon alignment per
ortholog pair with PAL2NAL v14. It runs the `pal2nal.pl` bundled inside the installed R package, once per pair:

```
perl <orthologr>/pal2nal/pal2nal.v14/pal2nal.pl <pair>_AA.aln <pair>.fasta -output fasta > <pair>_pal2nal.aln
```

`codon_aln()` accepts only `tool = "pal2nal"` and has no argument, option or environment variable for another
program (orthologr `codon_aln.R`, checked October 2026), so COAL's `pal2nal` command on `PATH` is not picked up.
Three ways in. (Licences: COAL is GPL-3.0-only from the release after 0.5.0; orthologr's DESCRIPTION says GPL-2.
Way 1 and the `coal-cpp` program start COAL as a separate program, so the two are not combined. Way 2 loads `coalr`
(GPL-3) into the same R process as orthologr: see the licence note there.)

## 1. Replace orthologr's bundled script with the stand-in (works today)

`pal2nal.pl` in this folder is a small Perl script that hands the unchanged PAL2NAL command line to COAL:
to `coal-cpp` when it is installed (it starts in a few milliseconds, which is what matters when a program runs
once per gene pair), otherwise to COAL's `pal2nal` command (Python); `COAL_PAL2NAL=<path>` forces one. A command
line with a PAL2NAL option that `coal-cpp` does not take (`-nostderr`, `-h`) always goes to `pal2nal`, which
accepts all of PAL2NAL's. Output is
PAL2NAL v14's, byte for byte, apart from the [documented differences](../../docs/differences.md).

In R, with COAL installed (`pip install coal-align`, or `conda install -c bioconda coal-align` once the recipe is
accepted; both put `coal`, `pal2nal` and `coal-cpp` on `PATH`):

```r
standin <- "path/to/COAL/packaging/orthologr/pal2nal.pl"   # this file
target  <- system.file("pal2nal", "pal2nal.v14", "pal2nal.pl", package = "orthologr", mustWork = TRUE)
file.copy(target, paste0(target, ".orig"))                 # keep the original
file.copy(standin, target, overwrite = TRUE)
```

Reinstalling or updating orthologr restores the original file. To check which program runs:

```r
system2("perl", c(target, "--version"))   # "coal 0.4.2 (C++)" or "coal 0.4.2"; PAL2NAL v14 prints its usage
```

## 2. The fork with COAL built in (Plumat-COAL/orthologr, branch `feat/coal`)

The fork adds COAL to orthologr as an optional dependency (the R package [`coalr`](../r/coalr), see
[docs/r.md](../../docs/r.md); or the `coal-cpp` program when `coalr` is not installed). Not merged upstream.

```r
remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")
remotes::install_github("Plumat-COAL/orthologr", ref = "feat/coal")

dNdS(q, s, codon_aln_tool = "coal")                                         # PAL2NAL's step without Perl
dNdS(q, s, dnds_est.method = "YN", codon_aln_tool = "coal", kaks_tool = "coal")   # ... and without KaKs_Calculator
dNdS(q, s, dnds_est.method = "YN", engine = "coal", comp_cores = 8)         # all pairs in one call to COAL
dNdS_clusters(c(q, s), aa_aln_tool = "clustalo")                            # every pair of each DeepClust cluster
```

- `codon_aln_tool = "coal"` gives the same codon alignments as `"pal2nal"` (so the same dN/dS), one pair at a
  time; the gain is robustness (no Perl, no length limit, a checked result).
- `kaks_tool = "coal"` computes NG, LWL, MLWL, LPB, MLPB, YN and MYN with COAL's engine ([docs/kaks.md](../../docs/kaks.md)).
- `engine = "coal"` (`dNdS_coal()`) keeps orthologr's DIAMOND/BLAST reciprocal-best-hit search and hands every
  pair to `coalr::coal_dnds_pairs()`: no file and no process per pair. It returns `dNdS()`'s table.
- `dNdS_clusters()` uses `deepclust()` / `deepclust_realign()` tables with `coalr::coal_dnds_clusters()`.
- orthologr's own two methods, `"Comeron"` (its default) and `"Li"`, are in COAL's engine too (`methods =
  "Comeron"`, `"Li"` in coalr, with `comeron_flavour = "gestimator"` / `li_flavour = "seqinr"` for orthologr's
  numbers). The fork does not use them yet: its `engine = "coal"` and `kaks_tool = "coal"` accept
  KaKs_Calculator's seven methods. [docs/kaks.md](../../docs/kaks.md), "Comeron and Li", has the validation
  against orthologr and the five things an orthologr build has to do to offer them.

The fork's `NEWS.md` also lists its fixes to orthologr's own path (BLOSUM62 scoring and full-length alignments in
`pairwise_aln()`, exit-status checks, one failing pair no longer ends a run, hyphens in gene IDs, input files of
the same name).

**Licence note (open question for the maintainers).** `coalr` is GPL-3 (it contains COAL's Ka/Ks engine, derived
from KaKs_Calculator 2.0 and PAML) and orthologr's DESCRIPTION says `GPL-2`. The fork contains no COAL code and
loads `coalr` only when it is installed, but a GPL-2-only package calling a GPL-3-only one in the same R process is
a combination the two licences do not clearly allow. Before the fork is merged upstream or distributed with `coalr`
as a dependency, orthologr's authors should decide (for instance by licensing orthologr as `GPL (>= 2)` or
`GPL-3`). Using the `coal-cpp` program instead (a separate process) does not raise the question.

## 3. Add `tool = "coal"` to upstream orthologr yourself


In `R/codon_aln.R`, accept `"coal"` in the `tool` check and, in the `tool == "coal"` branch, run COAL directly
(no Perl, no shell redirection, and a checked exit status):

```r
# untested sketch; COAL takes PAL2NAL's -output and -codontable options and writes the file itself with -o
exe <- Sys.which("coal-cpp"); if (exe == "") exe <- Sys.which("coal")
if (exe == "") stop("COAL is not installed: pip install coal-align (or conda install -c bioconda coal-align)")
status <- system2(exe, c(shQuote(file_aln), shQuote(file_nuc), "-output", format, "-o", shQuote(file.out), params))
if (status != 0) stop("COAL failed on ", file_aln, " (exit status ", status, "); see its messages above")
```

and let `dNdS(codon_aln_tool = "coal")` pass it through (`compute_dnds()` already forwards `codon_aln_tool`
unchecked). COAL exits non-zero when a gene cannot be converted (PAL2NAL v14 always exits 0, and orthologr does
not check its status, so a failed pair surfaces only later, if at all).

## What to expect

Measured here (Linux, one pair per process, median of 5–20 runs; outputs identical in every case):

| Pair | PAL2NAL v14 (Perl) | stand-in → coal-cpp |
|---|---:|---:|
| 2 × 150 codons | 10.7 ms | 13.3 ms |
| 2 × 500 codons | 16.2 ms | 16.2 ms |
| 2 × 2,000 codons | 50.4 ms | 12.5 ms |
| 2 × 5,000 codons | 208.7 ms | 13.4 ms |

- **Typical gene lengths:** the time per pair is process start-up (Perl's or coal-cpp's), so COAL is no faster.
  The stand-in via COAL's Python `pal2nal` is slower (about 90 ms per pair) because of Python start-up: that is
  why the stand-in prefers `coal-cpp`.
- **Long genes:** PAL2NAL's time grows faster than the gene length; COAL's does not. From about 2,000 codons COAL
  is faster per pair, 15× at 5,000 codons.
- **Correctness:** on 2,595 real Ensembl Compara trees COAL matched the truth on all; PAL2NAL v14 failed 22
  (19 with selenocysteine, 1 over Perl's regex capture-group limit, which stops it on proteins longer than about
  22,000 residues; `docs/benchmarks.md`). Those pairs come back as missing dN/dS values in orthologr today.
- **A real speed-up for whole-genome runs** needs one COAL call for all pairs instead of one process per pair:
  that is `dNdS(engine = "coal")` in the fork (section 2), built on `coalr::coal_dnds_pairs()`.
