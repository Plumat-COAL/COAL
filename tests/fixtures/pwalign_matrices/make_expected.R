# expected_<MATRIX>.fasta: pwalign::pairwiseAlignment(type = "global", substitutionMatrix = MATRIX) (pwalign's
# default gapOpening = 10, gapExtension = 4) of the pairs in pairs.fasta, for each matrix COAL has and pwalign ships:
# the full-length rows, the score in the name line. pwalign refuses a letter its matrix has no row for (U and O; J
# in PAM250): for such a pair the alignment is pwalign's of the pair with those letters as X (which is how COAL
# reads them), and the name line says so ("as_x=UO").
#     Rscript make_expected.R
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
x <- Biostrings::readAAStringSet(file.path(here, "pairs.fasta"))
nm <- sub(" .*", "", names(x))
for (m in c("BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250")) {
  e <- new.env()
  utils::data(list = m, package = "pwalign", envir = e)
  letters <- rownames(e[[m]])
  out <- character(length(x))
  name <- nm
  for (k in seq(1, length(x), by = 2)) {
    a <- as.character(x[[k]])
    b <- as.character(x[[k + 1]])
    foreign <- setdiff(unique(strsplit(paste0(a, b), "")[[1]]), letters)
    note <- ""
    if (length(foreign)) {
      re <- paste0("[", paste(foreign, collapse = ""), "]")
      a <- gsub(re, "X", a)
      b <- gsub(re, "X", b)
      note <- paste0(" as_x=", paste(sort(foreign), collapse = ""))
    }
    al <- pwalign::pairwiseAlignment(Biostrings::AAString(a), Biostrings::AAString(b), type = "global",
                                     substitutionMatrix = m)
    out[k] <- as.character(pwalign::alignedPattern(al))
    out[k + 1] <- as.character(pwalign::alignedSubject(al))
    name[k] <- paste0(nm[k], " score=", pwalign::score(al), note)
  }
  writeLines(paste0(">", name, "\n", out), file.path(here, paste0("expected_", m, ".fasta")))
}
cat("pwalign", as.character(utils::packageVersion("pwalign")), "\n")
