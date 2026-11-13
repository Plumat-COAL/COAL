# SPDX-License-Identifier: GPL-3.0-only
# Reference values for Method::Li: seqinr::kaks() on every pair of an AXT file, as orthologr's
# substitutionrate(est.method = "Li") calls it (an alignment read by seqinr::read.alignment holds lower-case
# rows; kaks() turns them to upper case), written with 17 significant digits.
#
#   Rscript benchmarks/kaks/ref/seqinr_ref.R pairs.axt ref.tsv
#
# Output, one line per pair: name ka ks vka vks l0 l2 l4   (9.999999 = seqinr's "saturated"; "ERR" where kaks()
# refuses the pair: rows of unequal length or not a multiple of 3, or an empty row).
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
suppressPackageStartupMessages(library(seqinr))
lines <- sub("\r$", "", readLines(args[1], warn = FALSE))
out <- file(args[2], "w")
fmt <- function(x) formatC(x, digits = 17, format = "g")
rec <- character(0)
flush_rec <- function() {
    if (length(rec) == 0) return(invisible())
    if (length(rec) != 3 || nchar(rec[2]) == 0 || nchar(rec[2]) != nchar(rec[3]) || nchar(rec[2]) %% 3 != 0) {
        writeLines(paste(rec[1], "ERR", sep = "\t"), out)
    } else {
        aln <- seqinr::as.alignment(nb = 2, nam = c("a", "b"), seq = tolower(rec[2:3]))
        r <- suppressWarnings(seqinr::kaks(aln, verbose = TRUE))
        v <- c(r$ka[1], r$ks[1], r$vka[1], r$vks[1], r$l0[1], r$l2[1], r$l4[1])
        writeLines(paste(c(rec[1], fmt(v)), collapse = "\t"), out)
    }
    rec <<- character(0)
}
for (l in lines) {
    if (nchar(l) == 0) flush_rec() else rec <- c(rec, l)
}
flush_rec()
close(out)
cat("seqinr", as.character(packageVersion("seqinr")), "\n")
