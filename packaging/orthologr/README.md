# COAL in orthologr

[orthologr](https://github.com/drostlab/orthologr) (`dNdS()`, `codon_aln()`) builds one codon alignment per
ortholog pair with PAL2NAL v14. It runs the `pal2nal.pl` bundled inside the installed R package, once per pair:

```
perl <orthologr>/pal2nal/pal2nal.v14/pal2nal.pl <pair>_AA.aln <pair>.fasta -output fasta > <pair>_pal2nal.aln
```

`codon_aln()` accepts only `tool = "pal2nal"` and has no argument, option or environment variable for another
program (orthologr `codon_aln.R`, checked October 2026), so COAL's `pal2nal` command on `PATH` is not picked up.
Two ways in:

## 1. Replace orthologr's bundled script with the stand-in (works today)

`pal2nal.pl` in this folder is a 20-line Perl script that hands the unchanged PAL2NAL command line to COAL:
to `coal-cpp` when it is installed (it starts in a few milliseconds, which is what matters when a program runs
once per gene pair), otherwise to COAL's `pal2nal` command (Python); `COAL_PAL2NAL=<path>` forces one. Output is
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

## 2. Add `tool = "coal"` to orthologr (the proper fix; needs a pull request there)

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
- **A real speed-up for whole-genome runs** needs one COAL process for all pairs instead of one per pair, which is
  a change to orthologr's loop (a batch mode in COAL would make that call a one-liner; not built yet).
