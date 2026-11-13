q <- coal_read_fasta(fx("query_cds.fasta"))
s <- coal_read_fasta(fx("subject_cds.fasta"))

test_that("pairs equal `coal pipeline pairs` for every method", {
  want <- read_expected("pairs_dnds.tsv")
  got <- coal_dnds_pairs(q, s, pairs = orthologs(), methods = "ALL")
  expect_same_table(got, want)
  expect_identical(unique(got$method), c("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN"))
  expect_true(all(c("ok", "warning") %in% got$status))
})

test_that("file paths, lists and data frames are accepted and give the same table", {
  want <- coal_dnds_pairs(q, s, pairs = orthologs())
  attr(want, "seconds") <- NULL
  strip <- function(x) {
    attr(x, "seconds") <- NULL
    x
  }
  expect_identical(strip(coal_dnds_pairs(fx("query_cds.fasta"), fx("subject_cds.fasta"), orthologs())), want)
  ql <- lapply(q, function(x) strsplit(x, "")[[1]])
  expect_identical(strip(coal_dnds_pairs(ql, as.list(s), as.matrix(orthologs()))), want)
  qd <- data.frame(geneids = names(q), seqs = unname(q))
  expect_identical(strip(coal_dnds_pairs(qd, s, orthologs())), want)
  skip_if_not_installed("Biostrings")
  keep <- !duplicated(names(q)) & grepl("^[ACGTNRYacgtnry]+$", q)
  dss <- Biostrings::DNAStringSet(q[keep])
  o <- orthologs()
  o <- o[o$query %in% names(dss) & !o$query %in% c("rna_q"), ]
  a <- strip(coal_dnds_pairs(dss, s, o))
  b <- strip(coal_dnds_pairs(q, s, o))
  expect_equal(a[, 4:13], b[, 4:13])
})

test_that("the result does not depend on the number of threads", {
  one <- coal_dnds_pairs(q, s, pairs = orthologs(), methods = c("YN", "NG"), threads = 1)
  for (t in c(2, 4, 0)) {
    many <- coal_dnds_pairs(q, s, pairs = orthologs(), methods = c("YN", "NG"), threads = t)
    attr(one, "seconds") <- attr(many, "seconds") <- NULL
    expect_identical(many, one)
  }
})

test_that("the pwalign preset equals the command line's", {
  want <- read_expected("pairs_pwalign_dnds.tsv")
  got <- coal_dnds_pairs(q, s, pairs = orthologs(), preset = "pwalign")
  expect_same_table(got, want)
  same <- coal_dnds_pairs(q, s, pairs = orthologs(), gap_open = 14, gap_extend = 4, end_gaps = "penalized")
  expect_equal(got$dN, same$dN)
})

test_that("bad genes and unknown IDs are rows, not errors", {
  p <- data.frame(q = c("g01", "bad_q", "empty_q", "nobody", "g02", NA), s = c("h01", "h02", "h03", "h04", "ghost", "h05"))
  got <- coal_dnds_pairs(q, s, pairs = p)
  expect_identical(nrow(got), 6L)
  expect_identical(got$status, c("ok", rep("error", 5)))
  expect_true(all(is.na(got$dN[-1])))
  expect_identical(got$codons_compared[-1], rep(0L, 5))
  expect_match(got$note[2], "seq1: the CDS holds a character that is no nucleotide code")
  expect_match(got$note[3], "seq1: the CDS is empty")
  expect_match(got$note[4], "seq1: no CDS record with this ID")
  expect_match(got$note[5], "seq2: no CDS record with this ID")
  expect_identical(got$seq1[4], "nobody")
  # NA sequences
  q2 <- c(a = NA_character_, b = "ATGGCTAAAGGTCTGGAACGTTAA")
  s2 <- c(x = "ATGGCTAAAGGTCTGGAACGTTAA", y = "ATGGCTAAAGGTCTGGAACGTTAA")
  got <- coal_dnds_pairs(q2, s2)
  expect_identical(got$status, c("error", "ok"))
  expect_match(got$note[1], "the CDS is NA")
  expect_identical(got$dN[2], 0)
  expect_true(is.na(got$dNdS[2]))
})

test_that("a pair above max_cells is that pair's error", {
  got <- coal_dnds_pairs(q, s, pairs = orthologs()[1:3, ], max_cells = 2000)
  expect_true(all(got$status == "error"))
  expect_match(got$note[1], "protein alignment failed")
})

test_that("alignments are returned when asked for, and are consistent", {
  o <- orthologs()[1:5, ]
  got <- coal_dnds_pairs(q, s, pairs = o, keep_alignments = TRUE)
  al <- attr(got, "alignments")
  expect_identical(al$seq1, o$query)
  expect_identical(nchar(al$codon1), 3L * nchar(al$protein1))
  expect_identical(nchar(al$codon1), nchar(al$codon2))
  again <- coal_kaks(al$codon1, al$codon2)
  expect_equal(again$dN, got$dN)
  expect_equal(again$dS, got$dS)
  expect_identical(gsub("-", "", al$protein1), unname(coal_translate(q[o$query])))
  thread <- coal_codon_align(c(al$protein1[1], al$protein2[1]), c(q[[o$query[1]]], s[[o$subject[1]]]))
  expect_identical(unname(thread), c(al$codon1[1], al$codon2[1]))
})

test_that("arguments are checked", {
  expect_error(coal_dnds_pairs(q, s[1:3]), "same length")
  expect_error(coal_dnds_pairs(unname(q), s, orthologs()), "needs names")
  expect_error(coal_dnds_pairs(q, s, orthologs(), methods = "GY"), "not one COAL computes")
  expect_error(coal_dnds_pairs(q, s, orthologs(), code = 7), "genetic code 7 is not supported")
  expect_error(coal_dnds_pairs(q, s, orthologs(), matrix = "BLOSUM99"), "BLOSUM99")
  expect_error(coal_dnds_pairs(q, s, orthologs(), gap_open = 1, gap_extend = 5))
  expect_error(coal_dnds_pairs(q, s, orthologs(), threads = -1), "threads")
  expect_error(coal_dnds_pairs(q, s, pairs = 1:3), "two columns")
  expect_identical(nrow(coal_dnds_pairs(q, s, pairs = orthologs()[0, ])), 0L)
  expect_identical(nrow(coal_dnds_pairs(character(0), character(0))), 0L)
})

test_that("other genetic codes are used", {
  a <- c(x = "ATGTGAAAAGGTCTGGAACGTCTTGCA")   # TGA: stop in code 1, Trp in code 2 (vertebrate mitochondrial)
  b <- c(y = "ATGTGGAAAGGCCTGGAACGACTTGCA")
  one <- coal_dnds_pairs(a, b, code = 1)
  two <- coal_dnds_pairs(a, b, code = 2)
  expect_identical(one$status, "warning")
  expect_identical(two$status, "ok")
  expect_identical(two$codons_compared, 9L)
})

test_that("min_codons and max_ds flag rows without dropping them or changing a number", {
  q <- c(g1 = "ATGGCTAAAGGTCTGGAACGTTAA", g2 = "ATGAAACCCGGGTTTTAG", g3 = "ATGAAAGTTCTGGCAGGTCCGTTA")
  s <- c(h1 = "ATGGCGAAAGGCCTGGATCGTTAA", h2 = "ATGAAGCCAGGGTTCTAG", h3 = "ATGAAGGTTTTGGCAGGACCGCTA")
  plain <- coal_dnds_pairs(q, s, methods = c("YN", "NG"))
  expect_null(attr(plain, "flagged"))
  f <- coal_dnds_pairs(q, s, methods = c("YN", "NG"), min_codons = 8, max_ds = 0.5)
  num <- setdiff(names(plain), c("status", "note"))
  expect_identical(f[num], plain[num])
  expect_identical(nrow(f), nrow(plain))
  has <- plain$status %in% c("ok", "warning")
  short <- has & plain$codons_compared < 8
  high <- has & !is.na(plain$dS) & plain$dS > 0.5
  expect_true(any(short) && any(high) && any(!short))
  expect_identical(grepl("fewer than 8 codons compared", f$note, fixed = TRUE), short)
  expect_identical(grepl("dS above 0.5", f$note, fixed = TRUE), high)
  expect_true(all(f$status[short | high] == "warning"))
  expect_identical(f$status[!(short | high)], plain$status[!(short | high)])
  expect_identical(f$note[!(short | high)], plain$note[!(short | high)])
  both <- short & high
  if (any(both)) expect_true(all(endsWith(f$note[both], "fewer than 8 codons compared; dS above 0.5")))
  expect_identical(attr(f, "flagged"), c(min_codons = sum(short), max_ds = sum(high)))
  # rows without an estimate are left alone; the notes are appended to an existing one
  sat <- coal_kaks(strrep("AAACCCGGGTTT", 5), strrep("CCCGGGTTTAAA", 5), method = "NG", min_codons = 1000, max_ds = 0.01)
  if (!sat$status %in% c("ok", "warning")) expect_false(grepl("fewer than", sat$note))
  k <- coal_kaks(c(a = "ATGAAAGTTCTGGCAGGTCCGTTA"), c(b = "ATGAAGGTTTTGGCAGGACCGCTA"), min_codons = 50)
  expect_identical(k$status, "warning")
  expect_identical(k$note, "fewer than 50 codons compared")
  expect_error(coal_dnds_pairs(q, s, min_codons = -1), "min_codons")
  expect_error(coal_dnds_pairs(q, s, max_ds = 0), "max_ds")
  expect_error(coal_kaks("ATGAAA", "ATGAAG", max_ds = "2"), "max_ds")
})

test_that("the cluster summary does not change with the flags", {
  cds <- c(a = "ATGAAAGTTCTGGCAGGTCCGTTA", b = "ATGAAGGTTTTGGCAGGACCGCTA", c = "ATGAAAGTTCTTGCAGGTCCGCTA")
  cl <- data.frame(rep = c("a", "a", "a"), member = c("a", "b", "c"))
  rows <- c(a = "MKVLAGPL", b = "MKVLAGPL", c = "MKVLAGPL")
  plain <- coal_dnds_clusters(cds, cl, msa = list(a = rows))
  f <- coal_dnds_clusters(cds, cl, msa = list(a = rows), min_codons = 100)
  expect_identical(f$summary, plain$summary)
  expect_true(all(grepl("fewer than 100 codons compared", f$pairs$note[f$pairs$status == "warning"])))
  expect_identical(f$pairs$dN, plain$pairs$dN)
  expect_identical(unname(attr(f$pairs, "flagged")), sum(plain$pairs$status %in% c("ok", "warning")))
})
