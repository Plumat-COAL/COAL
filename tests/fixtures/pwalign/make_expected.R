# expected.fasta: pwalign::pairwiseAlignment(type = "global", substitutionMatrix = "BLOSUM62") (pwalign's default
# gapOpening = 10, gapExtension = 4) of the pairs in pairs.fasta, the full-length rows; the score in the name line.
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
x <- Biostrings::readAAStringSet(file.path(here, "pairs.fasta"))
out <- character(length(x))
nm <- names(x)
for (k in seq(1, length(x), by = 2)) {
  al <- pwalign::pairwiseAlignment(x[[k]], x[[k + 1]], type = "global", substitutionMatrix = "BLOSUM62")
  out[k] <- as.character(pwalign::alignedPattern(al))
  out[k + 1] <- as.character(pwalign::alignedSubject(al))
  nm[k] <- paste0(nm[k], " score=", pwalign::score(al))
}
writeLines(paste0(">", nm, "\n", out), file.path(here, "expected.fasta"))
cat("pwalign", as.character(utils::packageVersion("pwalign")), "\n")
