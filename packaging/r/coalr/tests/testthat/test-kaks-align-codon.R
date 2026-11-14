msa_tbl <- utils::read.delim(fx("clusters_msa.tsv"), stringsAsFactors = FALSE)
famA <- msa_tbl[msa_tbl$cluster_id == "famA_001", ]
cds <- coal_read_fasta(fx("family_cds.fasta"))

test_that("coal_kaks on a codon alignment equals the pipeline's rows", {
  want <- read_expected("clusters_dnds.tsv")
  want <- want[want$cluster_id == "famA_001", ]
  rows <- structure(famA$codon_row, names = famA$id)
  got <- coal_kaks(rows, method = c("YN", "NG"))
  cols <- c("seq1", "seq2", "method", "dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "Nd", "Sd",
            "codons_compared", "status", "note")
  expect_same_table(got, want[cols])
  ij <- utils::combn(nrow(famA), 2)
  two <- coal_kaks(rows[ij[1, ]], rows[ij[2, ]], method = c("NG", "YN"), threads = 2)
  expect_equal(two[, -(1:2)], got[, -(1:2)])
  d <- coal_kaks(rows[1], rows[2], details = TRUE)
  expect_true(all(c("p_fisher", "kappa", "length", "gc") %in% names(d)))
  expect_identical(d$length, 3L * d$codons_compared)
})

test_that("coal_kaks reports bad pairs in status", {
  r <- coal_kaks(c("ATGAAA", "ATGAAAC", NA, "ATG---", "ATGTAAAAA"), c("ATGAAG", "ATGAAA", "ATG", "---ATG", "ATGAAAAAA"))
  expect_identical(r$status, c("ok", "error", "error", "no_data", "warning"))
  expect_match(r$note[3], "NA")
  expect_error(coal_kaks("ATG"), "two sequences")
  expect_error(coal_kaks("ATG", c("ATG", "ATG")), "same length")
  expect_error(coal_kaks("ATG", "ATG", method = "MA"), "not one COAL computes")
  paml <- coal_kaks(famA$codon_row[1], famA$codon_row[2], yn_flavour = "paml", yn_weighting = FALSE)
  expect_identical(paml$status, "ok")
})

test_that("coal_codon_align equals the pipeline's codon alignment", {
  rows <- structure(famA$protein_row, names = famA$id)
  got <- coal_codon_align(rows, cds)
  expect_identical(unname(got), famA$codon_row)
  expect_identical(names(got), famA$id)
  # by position, with shuffled names ignored
  expect_identical(unname(coal_codon_align(unname(rows), unname(cds[famA$id]))), famA$codon_row)
  # -nogap removes the gapped columns
  ng <- coal_codon_align(rows, cds, nogap = TRUE)
  expect_identical(length(unique(nchar(ng))), 1L)
  expect_false(any(grepl("-", ng)))
  gap_cols <- Reduce(`|`, lapply(strsplit(famA$protein_row, ""), function(r) r == "-"))
  expect_identical(unname(nchar(ng[1])), 3L * sum(!gap_cols))
})

test_that("coal_codon_align explains what it cannot thread", {
  rows <- c(a = "MAKG-LER", b = "MAKGALDR")
  dna <- c(a = "ATGGCTAAAGGTCTGGAACGTTAA", b = "ATGGCGAAAGGCGCACTGGATCGTTAA")
  expect_identical(unname(coal_codon_align(rows, dna)), c("ATGGCTAAAGGT---CTGGAACGT", "ATGGCGAAAGGCGCACTGGATCGT"))
  expect_error(coal_codon_align(rows, dna[1]), "no coding sequence named 'b'")
  expect_error(coal_codon_align(c(rows, c = "MAK"), c(dna, c = "ATGGCTAAA")), "differ in length")
  wrong <- c(a = dna[["a"]], b = "CCCCCCCCCCCCCCCCCCCCCCCCCCC")
  expect_error(coal_codon_align(rows, wrong), "b")
  expect_warning(res <- coal_codon_align(rows, wrong, keep_going = TRUE), "could not be threaded")
  expect_identical(unname(is.na(res)), c(FALSE, TRUE))
  expect_identical(attr(res, "problems")$id, "b")
})

test_that("coal_align_pairs is a full-length global alignment", {
  r <- coal_align_pairs("MAKGLERWW", "MAKGLDR")
  expect_identical(r$aligned_a, "MAKGLERWW")
  expect_identical(r$aligned_b, "MAKGLDR--")
  expect_identical(r$identities, 6L)
  expect_true(is.na(r$error))
  pw <- coal_align_pairs("MAKGLERWW", "MAKGLDR", preset = "pwalign")
  expect_lt(pw$score, r$score)       # end gaps cost there
  many <- coal_align_pairs(rep(c("MKTAYIAKQR", "MSSHEGGKKK"), 50), "MKTAIAKQR", threads = 2)
  expect_identical(nrow(many), 100L)
  expect_identical(many[1, ], many[99, ], ignore_attr = TRUE)
  expect_identical(gsub("-", "", many$aligned_a[2]), "MSSHEGGKKK")
  bad <- coal_align_pairs(c("MK", NA, "M-K"), c("MK", "MK", "MK"))
  expect_identical(is.na(bad$aligned_a), c(FALSE, TRUE, TRUE))
  expect_match(bad$error[3], "not an amino acid")
  big <- coal_align_pairs(strrep("MKTAYIAKQR", 20), strrep("MKTAYIAKQR", 20), max_cells = 100)
  expect_match(big$error, "cells")
  expect_error(coal_align_pairs(c("MK", "MK"), c("MK", "MK", "MK")), "same length")
})

test_that("coal_translate follows the pipeline's rules", {
  genes <- utils::read.delim(fx("pairs_genes.tsv"), stringsAsFactors = FALSE, quote = "")
  q <- coal_read_fasta(fx("query_cds.fasta"))
  q <- q[!duplicated(names(q))]
  d <- coal_translate(q, details = TRUE)
  g <- genes[genes$set == "query" & genes$code != "duplicate_id", ]
  g <- g[match(d$id, g$id), ]
  expect_identical(d$nt_length, g$nt_length)
  expect_identical(d$aa_length, g$aa_length)
  expect_identical(d$status, g$status)
  expect_identical(d$code, ifelse(is.na(g$code), "", g$code))
  expect_identical(d$internal_stops, g$internal_stops)
  p <- coal_translate(c(a = "ATGGCTAAATAA", b = "atgnnngcn", c = "ATG-CC", d = NA, e = "AUGUGA", f = "ATGTAGAAA"))
  expect_identical(p, c(a = "MAK", b = "MXA", c = NA, d = NA, e = "M", f = "M*K"))
  expect_identical(unname(coal_translate("ATGTGA", code = 2)), "MW")
})

test_that("coal_read_fasta reads what it should", {
  f <- tempfile(fileext = ".fa")
  on.exit(unlink(f))
  writeLines(c("# comment", ">a desc", "ATG", "AAA", "", ">b", ">c\tx", "AT G"), f)
  expect_identical(coal_read_fasta(f), c(a = "ATGAAA", b = "", c = "ATG"))
  expect_error(coal_read_fasta(tempfile()), "no such file")
  writeLines("ATG", f)
  expect_error(coal_read_fasta(f), "not a FASTA file")
  # the package version is COAL's; a development tree (COAL 0.6.0.dev0) is 0.6.0.9000 in R's form
  v <- coal_info()$coal_version
  dev <- regmatches(v, regexec("^([0-9.]+)\\.dev([0-9]+)$", v))[[1]]
  if (length(dev) == 3) v <- paste0(dev[2], ".", 9000 + as.integer(dev[3]))
  expect_identical(v, as.character(utils::packageVersion("coalr")))
})

test_that("one string that looks like a path but names no file is 'file not found', never a sequence", {
  missing <- file.path(tempdir(), "no_such_folder", "cds.fasta")
  expect_error(coal_translate(missing), "`cds`: file not found: '.*cds\\.fasta'")
  expect_error(coal_translate("cds.fasta"), "file not found")             # an extension is enough
  expect_error(coal_translate("genes.fa.gz"), "file not found")
  expect_error(coal_translate("~/no/such/file"), "file not found")
  expect_error(coal_translate("C:\\data\\cds"), "file not found")           # a Windows path
  expect_error(coal_dnds_pairs("query_cds.fasta", c(a = "ATGAAA")), "`query_cds`: file not found")
  expect_error(coal_codon_align("aln.fasta", c(a = "ATGAAA")), "`protein_rows`: file not found")
  expect_error(coal_align_pairs("a.faa", "MKV"), "`a`: file not found")
  expect_error(coal_translate(tempdir()), "is a folder, not a FASTA file")
  # the message says how the string was read and what to do for a single sequence
  expect_error(coal_translate("x/y"), "read as a FASTA file name; for a single sequence give a named vector")
  # a plain single string is still one sequence, and a named one always is
  expect_identical(unname(coal_translate("ATGAAATAA")), "MK")
  expect_identical(unname(coal_translate(c(x = "ATGAAATAA"))), "MK")
  expect_no_error(coal_align_pairs("MKVL", "MKL"))
  # an existing file is read, whatever its name
  f <- tempfile("noext")
  writeLines(c(">g", "ATGAAATAA"), f)
  expect_identical(coal_translate(f), c(g = "MK"))
  expect_true(coalr:::looks_like_path("data/cds"))
  expect_false(coalr:::looks_like_path("ATGAAA"))
  expect_false(coalr:::looks_like_path("MKV-L.A*"))
})

test_that("an aligned row is a sequence whatever it ends in or contains (never 'file not found')", {
  # independent review 2026-11-13, finding 7: dot-gapped rows that end in the letters of a file extension, rows
  # with a frameshift mark (/ or \) and rows with GCG's ~ gaps were called missing files
  rows <- c("MKV..LFA.FA", "MKV.FA", "ACGT.AA", "ACG.NT", "MKL.SEQ", "MK.PEP", "AC.TXT", "MKVLA.FASTA", "MK.CDS",
            "MKV/LA", "MK\\VL", "~~MKV", "~MKV..FA", "MKV~~", "mk..v.fa", "acgt-ac.nt", "ac.gt.aa", "mkv*.fa",
            "MKV.FA.GZ", "M/K.FA")
  for (r in rows) expect_false(coalr:::looks_like_path(r), info = r)
  # ... and is used as one
  expect_no_error(coal_align_pairs("MKV.FA", "MKVFA"))
  expect_no_error(coal_align_pairs("MKV..LFA.FA", "MKVLFAFA"))
  # what is written like a file name still is one
  paths <- c("cds.fasta", "genes.fa.gz", "aln.fasta", "a.faa", "x/y", "~/no/such/file", "data/cds", "Data/Cds.fasta",
             "C:\\data\\cds", "query_cds.fasta", "run2/cds", "genes2.fa", "my genes.txt", "CDS_1.FASTA", "~/CDS1",
             "..\\cds.fa", "./cds.fa")
  for (p in paths) expect_true(coalr:::looks_like_path(p), info = p)
  # the two readings the rule cannot tell apart, as documented
  expect_false(coalr:::looks_like_path("CDS.FASTA"))     # a mistyped upper-case file name: taken as a sequence
  expect_true(coalr:::looks_like_path("acgtn.nt"))       # a lower-case row with one dot before "nt": a file name
})

test_that("an iteration that does not converge gives no estimate", {
  # AAA against ATG: the YN (10 rounds) and MYN (100 rounds) iterations wander. The last round's values differ
  # between platforms (Ka 3.169, 3.174 or saturated), so nothing of them is tested: only the status and the NAs.
  r <- coal_kaks("AAA", "ATG", method = c("NG", "YN", "MYN"), details = TRUE)
  expect_identical(r$method, c("NG", "YN", "MYN"))
  nc <- r[r$method != "NG", ]
  expect_identical(nc$status, c("not_converged", "not_converged"))
  for (v in c("dN", "dS", "dNdS", "SE_dN", "SE_dS", "Nd", "Sd", "p_fisher", "divergence")) {
    expect_true(all(is.na(nc[[v]])), info = v)
  }
  expect_identical(nc$N_sites, c(2, 2))
  expect_identical(nc$S_sites, c(1, 1))
  expect_match(nc$note[1], "^YN00 iteration did not converge in 10 rounds: no estimate \\(last round: ")
  expect_match(nc$note[2], "^MYN iteration did not converge in 100 rounds: no estimate \\(last round: ")
  expect_false(r$status[1] == "not_converged")
  # without pathway weighting there is no iteration
  expect_false(coal_kaks("AAA", "ATG", yn_flavour = "paml", yn_weighting = FALSE)$status == "not_converged")
})

test_that("preset pwalign returns pwalign::pairwiseAlignment's alignment among co-optimal ones", {
  # each pair has several alignments with the best score; COAL's own tie-break puts the gaps leftmost
  a <- c("MKTAYIAKQRQISFVKSHFSRQ", "MAAAKAAA", "MKVLAAGG", "GGGGAGGG", "MSTNPKPQRKTKRNTNRRPQDVKFPGG")
  b <- c("MKTAYAKQRQISFVKSHFSRQLE", "MAAKAAAA", "MKVLAGG", "GGGAGGGG", "MSTNPKPQKKNKRNTNRRPQDVKFPG")
  pw <- coal_align_pairs(a, b, preset = "pwalign")
  left <- coal_align_pairs(a, b, preset = "pwalign", tie_break = "leftmost")
  expect_identical(pw$score, left$score)
  expect_false(identical(pw$aligned_a, left$aligned_a) && identical(pw$aligned_b, left$aligned_b))
  expect_identical(coal_align_pairs(a, b, tie_break = "pwalign", gap_open = 14, gap_extend = 4,
                                    end_gaps = "penalized")$aligned_a, pw$aligned_a)
  expect_error(coal_align_pairs("MK", "MK", tie_break = "rightmost"), "tie_break")
  skip_if_not_installed("pwalign")
  for (i in seq_along(a)) {
    ref <- pwalign::pairwiseAlignment(Biostrings::AAString(a[i]), Biostrings::AAString(b[i]), type = "global",
                                      substitutionMatrix = "BLOSUM62")
    expect_identical(pw$aligned_a[i], as.character(pwalign::alignedPattern(ref)), info = a[i])
    expect_identical(pw$aligned_b[i], as.character(pwalign::alignedSubject(ref)), info = a[i])
    expect_equal(pw$score[i], pwalign::score(ref), info = a[i])
  }
})
