cds <- coal_read_fasta(fx("family_cds.fasta"))
msa_tbl <- utils::read.delim(fx("clusters_msa.tsv"), stringsAsFactors = FALSE)
msa <- lapply(split(msa_tbl, msa_tbl$cluster_id), function(d) structure(d$protein_row, names = d$id))

test_that("clusters with given MSA rows equal `coal pipeline clusters` (Clustal Omega)", {
  got <- coal_dnds_clusters(cds, fx("family_clusters.tsv"), msa = msa, methods = c("YN", "NG"))
  expect_same_table(got$pairs, read_expected("clusters_dnds.tsv"))
  want <- read_expected("clusters_summary.tsv")
  got$summary$aligner[got$summary$aligner == "msa"] <- "clustalo"   # the rows came from Clustal Omega
  want$note <- sub("failed.tsv", "failed", want$note, fixed = TRUE)
  expect_same_table(got$summary, want)
  expect_setequal(got$failed$seq1, c("famB_007", "famB_ghost"))
  expect_identical(got$failed$code[got$failed$seq1 == "famB_ghost"], "no_cds")
  expect_identical(got$failed$code[got$failed$seq1 == "famB_007"], "gene_failed")
  expect_true(all(got$failed$level == "gene"))
})

test_that("an aligner callback gives the same as the rows it returns", {
  calls <- 0L
  cb <- structure(function(prot) {
    calls <<- calls + 1L
    cid <- msa_tbl$cluster_id[match(names(prot)[1], msa_tbl$id)]
    expect_identical(gsub("-", "", msa[[cid]][names(prot)]), prot)
    rev(msa[[cid]][names(prot)])   # named, in another order
  }, name = "clustalo")
  got <- coal_dnds_clusters(cds, fx("family_clusters.tsv"), aligner = cb, methods = c("YN", "NG"))
  expect_identical(calls, 2L)   # famA and famB; famC has two members, famD one
  expect_same_table(got$pairs, read_expected("clusters_dnds.tsv"))
  expect_identical(unique(got$summary$aligner), c("clustalo", "nw", "none"))
})

test_that("a failing aligner falls back to pairwise alignment, and says so", {
  cb <- function(prot) if (length(prot) == 6L) stop("boom") else rep("MK", length(prot) - 1L)
  got <- coal_dnds_clusters(cds, fx("family_clusters.tsv"), aligner = cb)
  sm <- got$summary
  expect_identical(sm$aligner, c("nw_pairwise", "nw_pairwise", "nw", "none"))
  expect_match(sm$note[1], "custom failed \\(boom\\); pairs aligned pairwise")
  expect_match(sm$note[2], "returned 10 rows for 11 sequences")
  expect_identical(sm$n_pairs, c(15L, 55L, 1L, 0L))
  expect_identical(sum(got$failed$code == "msa_failed"), 2L)
  # and that is what no aligner at all gives
  plain <- coal_dnds_clusters(cds, fx("family_clusters.tsv"))
  expect_equal(plain$pairs, got$pairs)
})

test_that("clusters above the size cap: representative, skip, sample", {
  cl <- utils::read.delim(fx("family_clusters.tsv"), header = FALSE, stringsAsFactors = FALSE)
  got <- coal_dnds_clusters(cds, cl, max_cluster_size = 5)
  expect_same_table(got$pairs, read_expected("clusters_rep_dnds.tsv"))
  want <- read_expected("clusters_rep_summary.tsv")
  want$note <- sub("failed.tsv", "failed", sub("--max-cluster-size", "max_cluster_size", want$note, fixed = TRUE),
                   fixed = TRUE)
  expect_same_table(got$summary, want)

  skip <- coal_dnds_clusters(cds, cl, max_cluster_size = 5, large_clusters = "skip")
  want <- read_expected("clusters_skip_summary.tsv")
  want$note <- sub("failed.tsv", "failed", sub("--max-cluster-size", "max_cluster_size", want$note, fixed = TRUE),
                   fixed = TRUE)
  expect_same_table(skip$summary, want)
  expect_identical(skip$pairs$cluster_id, "famC_001")

  set.seed(42)
  before <- .Random.seed
  a <- coal_dnds_clusters(cds, cl, max_cluster_size = 5, large_clusters = "sample", pair_sample = 7, seed = 3)
  expect_identical(.Random.seed, before)   # R's random numbers are left alone
  b <- coal_dnds_clusters(cds, cl, max_cluster_size = 5, large_clusters = "sample", pair_sample = 7, seed = 3,
                          threads = 2)
  d <- coal_dnds_clusters(cds, cl, max_cluster_size = 5, large_clusters = "sample", pair_sample = 7, seed = 4)
  expect_equal(a$pairs, b$pairs)
  expect_identical(a$summary$pair_selection, c("sample:7", "sample:7", "all", "singleton"))
  expect_identical(a$summary$n_pairs, c(7L, 7L, 1L, 0L))
  expect_false(identical(a$pairs[, 1:3], d$pairs[, 1:3]))
  key <- paste(a$pairs$cluster_id, a$pairs$seq1, a$pairs$seq2)
  expect_identical(anyDuplicated(key), 0L)
  # every sampled pair has the value of the all-pairs pairwise run
  all_pw <- coal_dnds_clusters(cds, cl)$pairs
  m <- match(key, paste(all_pw$cluster_id, all_pw$seq1, all_pw$seq2))
  expect_false(anyNA(m))
  expect_equal(a$pairs$dN, all_pw$dN[m])
  # more pairs asked for than there are: all of them
  full <- coal_dnds_clusters(cds, cl, max_cluster_size = 5, large_clusters = "sample", pair_sample = 1000)
  expect_identical(full$summary$n_pairs, c(15L, 55L, 1L, 0L))
})

test_that("sample_pairs draws distinct ordered pairs", {
  p <- coalr:::sample_pairs(50, 200, 1, "c1")
  expect_identical(dim(p), c(2L, 200L))
  expect_true(all(p[1, ] < p[2, ] & p[1, ] >= 1 & p[2, ] <= 50))
  expect_identical(anyDuplicated(t(p)), 0L)
  expect_identical(p, coalr:::sample_pairs(50, 200, 1, "c1"))
  expect_false(identical(p, coalr:::sample_pairs(50, 200, 1, "c2")))
  expect_identical(ncol(coalr:::sample_pairs(5, 100, 1, "x")), 10L)
})

test_that("cluster tables are read as the command line reads them", {
  cl <- coalr:::as_clusters(data.frame(representative = c("representative", "a", "a", "b", "b", "c", "a"),
                                       member = c("member", "a", "x", "b", "x", "a", "y")))
  expect_identical(cl, list(a = c("a", "x", "y"), b = "b", c = "c"))
  expect_identical(coalr:::as_clusters(list(k = c("a", "b"))), list(k = c("a", "b")))
  expect_identical(unname(lengths(coalr:::as_clusters(fx("family_clusters.tsv")))), c(6L, 13L, 2L, 1L))
  expect_error(coalr:::as_clusters(1:3), "two-column")
})

test_that("rows that do not fit the sequences are that cluster's error", {
  bad <- msa
  bad$famA_001[2] <- paste0(bad$famA_001[2], "-")
  got <- coal_dnds_clusters(cds, fx("family_clusters.tsv"), msa = bad)
  a <- got$pairs[got$pairs$cluster_id == "famA_001", ]
  expect_true(all(a$status == "error"))
  expect_match(a$note[1], "differ in length")
  expect_identical(got$summary$status[1], "failed")
  expect_true(all(got$pairs$status[got$pairs$cluster_id == "famB_001"] == "ok"))
  expect_true("alignment_ragged" %in% got$failed$code)
})

test_that("kept alignments are the command line's codon alignments", {
  got <- coal_dnds_clusters(cds, fx("family_clusters.tsv"), msa = msa, keep_alignments = TRUE)
  al <- got$alignments
  m <- match(paste(msa_tbl$cluster_id, msa_tbl$id), paste(al$cluster_id, al$id))
  expect_false(anyNA(m))
  expect_identical(al$codon[m], msa_tbl$codon_row)
  expect_identical(al$protein[m], msa_tbl$protein_row)
})

test_that("Clustal Omega through coal_msa_aligner() reproduces the command line", {
  skip_on_cran()
  skip_if(!nzchar(Sys.which("clustalo")), "clustalo is not on the PATH")
  got <- coal_dnds_clusters(cds, fx("family_clusters.tsv"), aligner = coal_msa_aligner("clustalo"),
                            methods = c("YN", "NG"))
  expect_same_table(got$pairs, read_expected("clusters_dnds.tsv"))
  expect_identical(got$summary$aligner[1], "clustalo")
})

test_that("a missing MSA program is explained", {
  expect_error(coal_msa_aligner("clustalo", path = "/no/such/clustalo"), "was not found")
})

test_that("pairs without an estimate are counted apart in the summary", {
  # AAA against ATG: the YN iteration does not converge (rows are given: three members, so they are used). Such a
  # pair is not in n_ok nor in the statistics. Nothing of the last round is tested: it differs between platforms.
  seqs <- c(k1 = "AAA", m1 = "ATG", k2 = "AAA", cds[c("famA_001", "famA_002")])
  got <- coal_dnds_clusters(seqs, list(nc = c("k1", "m1", "k2"), fam = c("famA_001", "famA_002")),
                            msa = list(nc = c(k1 = "K", m1 = "M", k2 = "K")))
  expect_identical(names(got$summary), coalr:::SUMMARY_COLUMNS)
  expect_true("n_not_converged" %in% names(got$summary))
  expect_identical(got$pairs$status, c("not_converged", "no_data", "not_converged", "ok"))
  nc <- got$pairs[got$pairs$status == "not_converged", ]
  expect_true(all(is.na(nc[, c("dN", "dS", "dNdS", "SE_dN", "SE_dS", "Nd", "Sd")])))
  expect_true(all(grepl("did not converge", nc$note)))
  sm <- got$summary
  expect_identical(sm$n_pairs, c(3L, 1L))
  expect_identical(sm$n_ok, c(0L, 1L))
  expect_identical(sm$n_not_converged, c(2L, 0L))
  expect_identical(sm$n_no_data, c(1L, 0L))
  expect_identical(sm$n_saturated + sm$n_failed, c(0L, 0L))
  expect_true(is.na(sm$dN_median[1]))
  # no pair with an estimate and none failed: the pairs' status (saturated, else no_data, else not_converged)
  expect_identical(sm$status, c("no_data", "ok"))
  expect_identical(sm$note[1], "no pair has an estimate (3 pairs: 1 without codons to compare, 2 not converged)")
})

test_that("a two-member cluster is computed in input order, as coal_dnds_pairs computes the pair", {
  q <- coal_read_fasta(fx("query_cds.fasta"))
  s <- coal_read_fasta(fx("subject_cds.fasta"))
  pr <- orthologs()
  m <- c("YN", "LWL", "Comeron")
  pairs <- coal_dnds_pairs(q, s, pairs = pr, methods = m)
  # the representative is the subject's gene, which comes after the query's in `cds`
  cl <- data.frame(rep = rep(pr[[2]], each = 2), mem = as.vector(rbind(pr[[2]], pr[[1]])))
  clu <- coal_dnds_clusters(c(q, s), cl, methods = m)$pairs
  expect_identical(clu$seq1, pairs$seq2)
  for (v in c("dN", "dS", "status")) expect_identical(clu[[v]], pairs[[v]], info = v)
  # the other way round in `cds`: the subject's gene first, as coal_dnds_pairs with the sets exchanged
  clu2 <- coal_dnds_clusters(c(s, q), cl, methods = m)$pairs
  back <- coal_dnds_pairs(s, q, pairs = pr[, 2:1], methods = m)
  for (v in c("dN", "dS", "status")) expect_identical(clu2[[v]], back[[v]], info = v)
})
