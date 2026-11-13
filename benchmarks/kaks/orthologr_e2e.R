# SPDX-License-Identifier: GPL-3.0-only
# End to end: orthologr's dNdS(dnds_est.method = "Comeron") against coalr on the same gene pairs
# (docs/kaks.md, "Comeron and Li": "With orthologr").
#
#   Rscript benchmarks/kaks/orthologr_e2e.R <orthologr-lib> <coalr-lib> <query.fasta> <subject.fasta> <outdir>
#
# 1. orthologr::dNdS() with its default path (reciprocal best hits, pairwise protein alignment, PAL2NAL,
#    gestimator); timed. Its per-pair protein and codon alignment files are copied from the session's temporary
#    folder before R ends. (A fork of orthologr that marks saturated pairs as NA is asked for the 0.4.3 values
#    with options(orthologr.legacy_comeron_saturation = TRUE); orthologr 0.4.3 ignores the option.)
# 2. coalr::coal_kaks(method = "Comeron", comeron_flavour = "gestimator") on orthologr's own codon alignments:
#    the same input, so every value must be orthologr's (which gestimator prints with 4 significant digits).
# 3. coalr::coal_dnds_pairs(methods = "Comeron", comeron_flavour = "gestimator") on the same pairs, from the
#    coding sequences (COAL translates, aligns the proteins with preset = "pwalign", the scoring of orthologr's
#    default pairwise alignment, and threads the codons itself); timed. Pairs are split by whether COAL's
#    protein alignment is orthologr's (orthologr's rows without their last column when that is '*' in both).
# 4. coalr::coal_codon_align() on orthologr's protein alignments and the coding sequences, then coal_kaks():
#    the same pairs and the same protein alignments, so every value must be orthologr's.
# Writes <outdir>/orthologr.tsv, pairs.tsv (one row per pair with both results) and summary.txt.
a <- commandArgs(TRUE)
stopifnot(length(a) == 5)
olib <- a[1]; clib <- a[2]; q <- a[3]; s <- a[4]; out <- a[5]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(olib, clib, .libPaths()))
suppressMessages(library(orthologr))
options(orthologr.legacy_comeron_saturation = TRUE)
say <- function(...) cat(sprintf(...), "\n", sep = "", file = file.path(out, "summary.txt"), append = TRUE)
unlink(file.path(out, "summary.txt"))
say("orthologr %s (%s); coalr %s (%s)", as.character(packageVersion("orthologr")), olib,
    as.character(packageVersion("coalr", lib.loc = clib)), clib)

t0 <- proc.time()[["elapsed"]]
res <- suppressWarnings(suppressMessages(dNdS(query_file = q, subject_file = s, dnds_est.method = "Comeron", comp_cores = 1,
                             quiet = TRUE, print_citation = FALSE)))
t_orth <- proc.time()[["elapsed"]] - t0
res <- as.data.frame(res)
res <- res[, c("query_id", "subject_id", "dN", "dS", "dNdS")]
write.table(res, file.path(out, "orthologr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("orthologr dNdS(dnds_est.method = \"Comeron\"), default path, 1 core: %.1f s, %d pairs", t_orth, nrow(res))

# orthologr's alignment files of this session
files <- list.files(tempdir(), recursive = TRUE, full.names = TRUE)
aa_files <- grep("pairwise_aln/.*_AA\\.aln$", files, value = TRUE)
co_files <- grep("codon_aln/.*\\.aln$", files, value = TRUE)
read_fa <- function(f) {
  l <- readLines(f, warn = FALSE)
  h <- grep("^>", l)
  if (length(h) != 2) return(NULL)
  nm <- sub("^>\\s*", "", l[h]); nm <- sub("\\s.*$", "", nm)
  sq <- c(paste(l[(h[1] + 1):(h[2] - 1)], collapse = ""), paste(l[(h[2] + 1):length(l)], collapse = ""))
  list(names = nm, seq = gsub("\\s", "", sq))
}
aa <- Filter(Negate(is.null), lapply(aa_files, read_fa))
co <- Filter(Negate(is.null), lapply(co_files, read_fa))
key <- function(x) vapply(x, function(z) paste(z$names, collapse = "\t"), "")
names(aa) <- key(aa); names(co) <- key(co)
say("alignment files kept by orthologr: %d protein, %d codon", length(aa), length(co))
saveRDS(list(aa = aa, co = co), file.path(out, "orthologr_alignments.rds"))

pk <- paste(res$query_id, res$subject_id, sep = "\t")
fmt4 <- function(x) suppressWarnings(as.numeric(sprintf("%.4g", x)))   # gestimator prints 4 significant digits
same <- function(x, y) (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)

# 2. the engine on orthologr's own codon alignments
library(coalr, lib.loc = clib)
have <- pk %in% names(co)
c1 <- vapply(co[pk[have]], function(z) z$seq[1], ""); c2 <- vapply(co[pk[have]], function(z) z$seq[2], "")
k <- coal_kaks(unname(c1), unname(c2), method = "Comeron", comeron_flavour = "gestimator")
eng <- data.frame(dN = rep(NA_real_, nrow(res)), dS = NA_real_, dNdS = NA_real_, status = NA_character_)
eng[have, ] <- data.frame(k$dN, k$dS, k$dNdS, k$status)
ok_dn <- same(fmt4(eng$dN), res$dN); ok_ds <- same(fmt4(eng$dS), res$dS); ok_w <- same(fmt4(eng$dNdS), res$dNdS)
say("")
say("coal_kaks(Comeron, comeron_flavour = \"gestimator\") on orthologr's codon alignments (%d of %d pairs have a file):",
    sum(have), nrow(res))
say("  equal to orthologr at its 4 significant digits: dN %d, dS %d, dN/dS %d of %d", sum(ok_dn[have]), sum(ok_ds[have]),
    sum(ok_w[have]), sum(have))
dflt <- coal_kaks(unname(c1), unname(c2), method = "Comeron")
rel <- function(x, y) abs(x - y) / abs(y)
ch_dn <- !same(dflt$dN, k$dN); ch_ds <- !same(dflt$dS, k$dS)
say("  default flavour on the same alignments: dN changes in %d pairs (to NA: %d; by > 1%%: %d; by > 20%%: %d), dS in %d (to NA: %d; by > 1%%: %d; by > 20%%: %d)",
    sum(ch_dn), sum(is.na(dflt$dN) & !is.na(k$dN)), sum(rel(dflt$dN, k$dN) > 0.01, na.rm = TRUE),
    sum(rel(dflt$dN, k$dN) > 0.20, na.rm = TRUE), sum(ch_ds), sum(is.na(dflt$dS) & !is.na(k$dS)),
    sum(rel(dflt$dS, k$dS) > 0.01, na.rm = TRUE), sum(rel(dflt$dS, k$dS) > 0.20, na.rm = TRUE))
say("  pairs whose codon alignment holds a stop codon column that both sequences fill: %d",
    sum(dflt$codons_compared != k$codons_compared))

# 3. coalr from the coding sequences
qs <- coal_read_fasta(q); ss <- coal_read_fasta(s)
t0 <- proc.time()[["elapsed"]]
cp <- coal_dnds_pairs(qs, ss, pairs = res[, c("query_id", "subject_id")], methods = "Comeron",
                      comeron_flavour = "gestimator", preset = "pwalign", threads = 1L, keep_alignments = TRUE)
t_coal <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
invisible(coal_dnds_pairs(qs, ss, pairs = res[, c("query_id", "subject_id")], methods = "Comeron",
                          comeron_flavour = "gestimator", preset = "pwalign", threads = 1L))
t_coal2 <- proc.time()[["elapsed"]] - t0
al <- attr(cp, "alignments")
stopifnot(identical(paste(cp$seq1, cp$seq2, sep = "\t"), pk))
same_aln <- rep(NA, nrow(res))
has_aa <- pk %in% names(aa)
same_aln[has_aa] <- vapply(which(has_aa), function(i) {
  o <- aa[[pk[i]]]$seq
  # orthologr translates the terminal stop codon to '*' and aligns it; COAL's protein rows end before it
  if (all(grepl("\\*$", o))) o <- sub("\\*$", "", o)
  identical(toupper(o[1]), toupper(al$protein1[i])) && identical(toupper(o[2]), toupper(al$protein2[i]))
}, NA)
e_dn <- same(fmt4(cp$dN), res$dN); e_ds <- same(fmt4(cp$dS), res$dS); e_w <- same(fmt4(cp$dNdS), res$dNdS)
say("")
say("coal_dnds_pairs(methods = \"Comeron\", comeron_flavour = \"gestimator\", preset = \"pwalign\") from the coding sequences, 1 thread: %.2f s with the alignments returned, %.2f s without (orthologr: %.1f s; %.0f times)",
    t_coal, t_coal2, t_orth, t_orth / t_coal2)
for (grp in list(list("COAL's protein alignment is orthologr's", which(same_aln %in% TRUE)),
                 list("the protein alignments differ", which(same_aln %in% FALSE)),
                 list("no orthologr protein alignment file", which(is.na(same_aln))))) {
  i <- grp[[2]]
  if (!length(i)) next
  say("  %s: %d pairs; equal to orthologr at 4 digits: dN %d, dS %d, dN/dS %d", grp[[1]], length(i), sum(e_dn[i]),
      sum(e_ds[i]), sum(e_w[i]))
  d <- i[!e_dn[i] | !e_ds[i]]
  if (length(d)) {
    rdn <- rel(cp$dN[d], res$dN[d]); rds <- rel(cp$dS[d], res$dS[d])
    say("    of the %d that differ: median |rel. diff| dN %.3g, dS %.3g; max dN %.3g, dS %.3g; > 1%%: dN %d, dS %d; NA in one only: %d",
        length(d), median(rdn, na.rm = TRUE), median(rds, na.rm = TRUE), max(rdn, na.rm = TRUE), max(rds, na.rm = TRUE),
        sum(rdn > 0.01, na.rm = TRUE), sum(rds > 0.01, na.rm = TRUE),
        sum(xor(is.na(cp$dN[d]), is.na(res$dN[d])) | xor(is.na(cp$dS[d]), is.na(res$dS[d]))))
  }
}
# 4. coalr from the coding sequences and orthologr's protein alignments: coal_codon_align() threads the codons
#    onto orthologr's rows (their '*' column takes the stop codons, as PAL2NAL does), coal_kaks() estimates
t0 <- proc.time()[["elapsed"]]
idx <- which(has_aa)
cds_all <- c(qs, ss)
crow <- t(vapply(idx, function(i) {
  o <- aa[[pk[i]]]
  r <- tryCatch(coal_codon_align(stats::setNames(toupper(o$seq), o$names), cds_all[o$names]), error = function(e) c(NA, NA))
  unname(as.character(r))
}, c("", "")))
kk <- coal_kaks(crow[, 1], crow[, 2], method = "Comeron", comeron_flavour = "gestimator")
t_cl <- proc.time()[["elapsed"]] - t0
g_dn <- g_ds <- g_w <- rep(NA, nrow(res))
g_dn[idx] <- same(fmt4(kk$dN), res$dN[idx]); g_ds[idx] <- same(fmt4(kk$dS), res$dS[idx])
g_w[idx] <- same(fmt4(kk$dNdS), res$dNdS[idx])
same_codon <- vapply(seq_along(idx), function(n) {
  o <- co[[pk[idx[n]]]]
  !is.null(o) && identical(toupper(o$seq[1]), toupper(crow[n, 1])) && identical(toupper(o$seq[2]), toupper(crow[n, 2]))
}, NA)
say("")
say("coal_codon_align(orthologr's protein alignment, CDS) + coal_kaks(Comeron, comeron_flavour = \"gestimator\"): %.2f s",
    t_cl)
say("  codon alignment identical to orthologr's PAL2NAL output: %d of %d pairs", sum(same_codon), length(idx))
say("  equal to orthologr at 4 digits: dN %d, dS %d, dN/dS %d of %d pairs", sum(g_dn[idx]), sum(g_ds[idx]), sum(g_w[idx]),
    length(idx))
bad <- which(!(g_dn[idx] & g_ds[idx] & g_w[idx]))
for (n in head(bad, 10)) {
  i <- idx[n]
  say("    differs: %s / %s: orthologr %s %s %s, coalr %.6g %.6g %.6g (%s; %s)", res$query_id[i], res$subject_id[i],
      format(res$dN[i]), format(res$dS[i]), format(res$dNdS[i]), kk$dN[n], kk$dS[n], kk$dNdS[n], kk$status[n], kk$note[n])
}

both <- !is.na(cp$dN) & !is.na(res$dN) & !is.na(cp$dS) & !is.na(res$dS)
say("  all pairs: Pearson r dN %.6f, dS %.6f (n = %d)", cor(cp$dN[both], res$dN[both]), cor(cp$dS[both], res$dS[both]), sum(both))
write.table(data.frame(query_id = res$query_id, subject_id = res$subject_id, orthologr_dN = res$dN, orthologr_dS = res$dS,
                       orthologr_dNdS = res$dNdS, engine_dN = eng$dN, engine_dS = eng$dS, engine_dNdS = eng$dNdS,
                       coalr_dN = cp$dN, coalr_dS = cp$dS, coalr_dNdS = cp$dNdS, coalr_status = cp$status,
                       same_protein_alignment = same_aln),
            file.path(out, "pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(readLines(file.path(out, "summary.txt")), sep = "\n")
