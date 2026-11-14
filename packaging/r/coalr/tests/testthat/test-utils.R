test_that("preset pwalign takes pwalign's copy of the matrix named, as the command line does", {
  ap <- function(matrix, ...) {
    coalr:::align_params(matrix = matrix, gap_open = 10, gap_extend = 0.5, end_gaps = "free", max_cells = 4e9, ...)
  }
  expect_identical(ap("BLOSUM45", preset = "pwalign")$matrix, "BLOSUM45-PWALIGN")
  expect_identical(ap("blosum80", preset = "pwalign")$matrix, "BLOSUM80-PWALIGN")
  expect_identical(ap("BLOSUM62", preset = "pwalign")$matrix, "BLOSUM62-PWALIGN")
  expect_identical(ap(NULL, preset = "pwalign")$matrix, "BLOSUM62-PWALIGN")
  # pwalign's PAM250 is NCBI's: no copy
  expect_identical(ap("pam250", preset = "pwalign")$matrix, "PAM250")
  # a *-PWALIGN name is used as given
  expect_identical(ap("BLOSUM45-PWALIGN", preset = "pwalign")$matrix, "BLOSUM45-PWALIGN")
  # without the preset (or with needle) a plain name is NCBI's matrix
  expect_identical(ap("BLOSUM45")$matrix, "BLOSUM45")
  expect_identical(ap("BLOSUM45", preset = "needle")$matrix, "BLOSUM45")
  expect_identical(ap(NULL)$matrix, "BLOSUM62")
})

test_that("preset pwalign with a plain matrix name gives pwalign's alignments", {
  q <- coal_read_fasta(fx("query_cds.fasta"))
  s <- coal_read_fasta(fx("subject_cds.fasta"))
  pr <- utils::read.delim(fx("pairs_orthologs.tsv"), stringsAsFactors = FALSE)[, 1:2]
  plain <- coal_dnds_pairs(q, s, pairs = pr, preset = "pwalign", matrix = "BLOSUM45", keep_alignments = TRUE)
  named <- coal_dnds_pairs(q, s, pairs = pr, preset = "pwalign", matrix = "BLOSUM45-PWALIGN", keep_alignments = TRUE)
  ncbi <- coal_dnds_pairs(q, s, pairs = pr, matrix = "BLOSUM45", gap_open = 14, gap_extend = 4,
                          end_gaps = "penalized", tie_break = "pwalign", keep_alignments = TRUE)
  attr(plain, "seconds") <- attr(named, "seconds") <- attr(ncbi, "seconds") <- NULL
  expect_identical(plain, named)
  expect_identical(nrow(plain), nrow(ncbi))
})
