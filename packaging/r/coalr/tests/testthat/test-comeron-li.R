# The methods Comeron (libsequence's gestimator, as in orthologr) and Li (seqinr::kaks): the reference flavours
# against stored values of the two references, the default flavours, names, arguments, the three entry points.

refs <- function() {
  utils::read.delim(fx("comeron_li_refs.tsv"), sep = "\t", quote = "", stringsAsFactors = FALSE,
                    colClasses = "character")
}
num <- function(x) suppressWarnings(as.numeric(x))
expect_reference <- function(got, want, info) {
  expect_identical(is.na(got), is.na(want), info = info)
  ok <- !is.na(want)
  expect_true(all(abs(got[ok] - want[ok]) <= 1e-12 * pmax(abs(want[ok]), 1e-300)), info = info)
}

test_that("Comeron, flavour gestimator, gives the numbers of orthologr's gestimator", {
  r <- refs()
  got <- coal_kaks(r$seq1, r$seq2, method = "Comeron", comeron_flavour = "gestimator")
  expect_identical(nrow(got), nrow(r))
  expect_true(all(got$method == "Comeron"))
  refused <- r$g_dN == "ERR"                  # gestimator cannot take the pair
  expect_true(all(got$status[refused] == "error"))
  for (v in c("dN", "dS", "dNdS")) {
    want <- num(r[[paste0("g_", v)]])
    want[want == 999 | refused] <- NA         # 999: gestimator's "cannot be calculated"
    expect_reference(got[[v]], want, v)
  }
  expect_gt(sum(!is.na(got$dN)), 100)
  # a value that gestimator gets by setting an undefined logarithm to 0 is flagged
  flagged <- grepl("gestimator's value", got$note)
  expect_true(any(flagged))
  expect_true(all(got$status[flagged] %in% c("warning", "saturated", "no_data")))
})

test_that("Li, flavour seqinr, gives the numbers of seqinr::kaks", {
  r <- refs()
  got <- coal_kaks(r$seq1, r$seq2, method = "Li", li_flavour = "seqinr")
  refused <- r$s_ka == "ERR"
  expect_true(all(got$status[refused] == "error"))
  for (v in list(c("dN", "s_ka", "s_vka", "SE_dN"), c("dS", "s_ks", "s_vks", "SE_dS"))) {
    want <- num(r[[v[2]]])
    sentinel <- !is.na(want) & want == 9.999999   # seqinr's "saturated": never a number in coalr
    want[sentinel] <- NA
    expect_reference(got[[v[1]]], want, v[1])
    var <- num(r[[v[3]]])
    se <- suppressWarnings(sqrt(var))
    se[sentinel | is.na(want) | !is.finite(se)] <- NA
    expect_reference(got[[v[4]]], se, v[4])
  }
  expect_gt(sum(!is.na(got$dS)), 60)
})

test_that("seqinr::kaks itself agrees (when seqinr is installed)", {
  skip_if_not_installed("seqinr")
  r <- refs()
  r <- r[r$s_ka != "ERR", ][1:60, ]
  got <- coal_kaks(r$seq1, r$seq2, method = "Li", li_flavour = "seqinr")
  for (i in seq_len(nrow(r))) {
    k <- suppressWarnings(seqinr::kaks(seqinr::as.alignment(nb = 2, nam = c("a", "b"),
                                                           seq = tolower(c(r$seq1[i], r$seq2[i])))))
    ka <- k$ka[1]; ks <- k$ks[1]
    if (ka == 9.999999) expect_true(is.na(got$dN[i])) else expect_equal(got$dN[i], ka, tolerance = 1e-12)
    if (ks == 9.999999) expect_true(is.na(got$dS[i])) else expect_equal(got$dS[i], ks, tolerance = 1e-12)
  }
})

test_that("the default flavours apply COAL's policy", {
  r <- refs()
  row <- function(n) r[r$name == n, ]
  a0 <- row("typical")$seq1
  b0 <- row("typical")$seq2
  a <- paste0(a0, "TAA")
  b <- paste0(b0, "TGA")
  d <- coal_kaks(a, b, method = c("Comeron", "Li"))
  g <- coal_kaks(a, b, method = c("Comeron", "Li"), comeron_flavour = "gestimator", li_flavour = "seqinr")
  expect_identical(d$method, c("Comeron", "Li"))
  expect_identical(d$status, c("ok", "ok"))
  expect_identical(d$codons_compared, c(40L, 40L))     # the terminal stop codons are dropped
  expect_identical(g$codons_compared, c(41L, 41L))     # the references compare TAA with TGA
  expect_true(all(d$dN != g$dN))
  # without the stop codons the two flavours are one
  d2 <- coal_kaks(a0, b0, method = c("Comeron", "Li"))
  g2 <- coal_kaks(a0, b0, method = c("Comeron", "Li"), comeron_flavour = "gestimator", li_flavour = "seqinr")
  expect_identical(d2$dN, g2$dN)
  expect_identical(d2$dS, g2$dS)
  expect_identical(d$dN, d2$dN)
  expect_identical(d$dS, d2$dS)
  # Li has standard errors, Comeron has none
  expect_true(is.na(d$SE_dN[1]) && !is.na(d$SE_dN[2]) && !is.na(d$SE_dS[2]))
  # an internal stop: skipped with a warning by default, compared by the references
  i <- row("internal_stop_one_sequence")
  x <- coal_kaks(i$seq1, i$seq2, method = c("Comeron", "Li"))
  expect_identical(x$status, c("warning", "warning"))
  expect_match(x$note[1], "internal stop")
  expect_identical(x$codons_compared, c(39L, 39L))
  y <- coal_kaks(i$seq1, i$seq2, method = c("Comeron", "Li"), comeron_flavour = "gestimator", li_flavour = "seqinr")
  expect_identical(y$codons_compared, c(40L, 40L))
  expect_false(any(grepl("internal stop", y$note)))
  # saturation: NA and a status by default; gestimator's finite number is flagged in its flavour
  u <- row("unrelated_saturated")
  z <- coal_kaks(u$seq1, u$seq2, method = "Comeron")
  expect_identical(z$status, "saturated")
  expect_true(is.na(z$dS) || is.na(z$dN))
})

test_that("not_converged is a status of the iterating methods only", {
  # AAA against ATG: YN's iteration does not converge (no estimate); Comeron and Li do not iterate. The status
  # names come from the engine's enum, so a status added there cannot shift the names of the others here.
  r <- coal_kaks("AAA", "ATG", method = c("YN", "Comeron", "Li"))
  expect_identical(r$method, c("YN", "Comeron", "Li"))
  expect_identical(r$status[1], "not_converged")
  expect_true(all(r$status[2:3] %in% c("ok", "warning", "saturated", "no_data")))
  rf <- refs()
  row <- function(n) rf[rf$name == n, ]
  u <- row("unrelated_saturated")
  expect_identical(coal_kaks(u$seq1, u$seq2, method = "Li")$status, "saturated")
  g <- row("all_gap")
  expect_identical(coal_kaks(g$seq1, g$seq2, method = c("Comeron", "Li"))$status, c("no_data", "no_data"))
})

test_that("names, aliases and arguments", {
  a <- "ATGGCTAAAGGTCTGTTTGATCCGGAATGG"
  b <- "ATGGCAAAGGGTTTGTTCGACCCGGAGTGG"
  want <- coal_kaks(a, b, method = "Comeron")
  for (m in c("comeron", "COMERON", "C95", "c95")) expect_identical(coal_kaks(a, b, method = m), want)
  expect_identical(coal_kaks(a, b, method = "li"), coal_kaks(a, b, method = "Li"))
  expect_true(all(c("Comeron", "Li") %in% coal_info()$methods))
  # "ALL" stays KaKs_Calculator's seven methods; rows of several methods come in the engine's order
  expect_identical(coal_kaks(a, b, method = "ALL")$method, c("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN"))
  expect_identical(coal_kaks(a, b, method = c("Li", "Comeron", "YN"))$method, c("YN", "Comeron", "Li"))
  expect_error(coal_kaks(a, b, method = "Comeron", code = 2L), "standard genetic code only")
  expect_error(coal_kaks(a, b, method = "Li", code = 4L), "standard genetic code only")
  expect_error(coal_kaks(a, b, method = "Comeron", comeron_flavour = "x"), "should be one of")
  expect_error(coal_kaks(a, b, method = "Li", li_flavour = "x"), "should be one of")
  expect_error(coal_kaks(a, b, method = "C96"), "not one COAL computes")
})

test_that("coal_dnds_pairs and coal_dnds_clusters compute Comeron and Li", {
  q <- coal_read_fasta(fx("query_cds.fasta"))
  s <- coal_read_fasta(fx("subject_cds.fasta"))
  got <- coal_dnds_pairs(q, s, pairs = orthologs(), methods = c("Comeron", "Li", "YN"), keep_alignments = TRUE)
  expect_identical(unique(got$method), c("YN", "Comeron", "Li"))
  aln <- attr(got, "alignments")
  # the same numbers as coal_kaks on the codon alignments the function made
  ok <- !is.na(aln$codon1)
  direct <- coal_kaks(aln$codon1[ok], aln$codon2[ok], method = c("Comeron", "Li"))
  sub <- got[got$method %in% c("Comeron", "Li") & paste(got$seq1, got$seq2) %in% paste(aln$seq1, aln$seq2)[ok], ]
  expect_identical(sub$dN, direct$dN)
  expect_identical(sub$dS, direct$dS)
  # terminal stop codons: dropped also by the reference flavour, unless terminal_stops = "compare" (orthologr
  # 0.4.3's pipeline): the pair is then compared with its two stop codons as the last column
  ends <- function(x) toupper(substr(x, nchar(x) - 2L, nchar(x)))
  o <- orthologs()
  stop1 <- ends(q[o[[1]]]); stop2 <- ends(s[o[[2]]])
  both <- stop1 %in% c("TAA", "TAG", "TGA") & stop2 %in% c("TAA", "TAG", "TGA") & ok
  expect_gt(sum(both), 5)
  leg <- coal_dnds_pairs(q, s, pairs = o, methods = "Comeron", comeron_flavour = "gestimator")
  cmp <- coal_dnds_pairs(q, s, pairs = o, methods = "Comeron", comeron_flavour = "gestimator",
                         terminal_stops = "compare")
  k0 <- coal_kaks(aln$codon1[both], aln$codon2[both], method = "Comeron", comeron_flavour = "gestimator")
  k1 <- coal_kaks(paste0(aln$codon1[both], stop1[both]), paste0(aln$codon2[both], stop2[both]), method = "Comeron",
                  comeron_flavour = "gestimator")
  expect_identical(leg$dN[both], k0$dN)
  expect_identical(cmp$dN[both], k1$dN)
  expect_identical(cmp$dS[both], k1$dS)
  expect_identical(cmp$codons_compared[both], leg$codons_compared[both] + 1L)
  # the option has no effect on the default flavours
  expect_identical(coal_dnds_pairs(q, s, pairs = o, methods = c("Comeron", "YN"), terminal_stops = "compare")$dN,
                   coal_dnds_pairs(q, s, pairs = o, methods = c("Comeron", "YN"))$dN)
  expect_error(coal_dnds_pairs(q, s, pairs = o, methods = "Comeron", terminal_stops = "x"), "should be one of")
  expect_error(coal_dnds_pairs(q, s, pairs = orthologs(), methods = "Comeron", code = 11L),
               "standard genetic code only")
  # threads do not change the result
  expect_identical(coal_dnds_pairs(q, s, pairs = orthologs(), methods = c("Comeron", "Li"), threads = 1L)$dN,
                   coal_dnds_pairs(q, s, pairs = orthologs(), methods = c("Comeron", "Li"), threads = 4L)$dN)

  cds <- coal_read_fasta(fx("family_cds.fasta"))
  msa_tbl <- utils::read.delim(fx("clusters_msa.tsv"), stringsAsFactors = FALSE)
  rows <- lapply(split(msa_tbl, msa_tbl$cluster_id), function(d) structure(d$protein_row, names = d$id))
  res <- coal_dnds_clusters(cds, fx("family_clusters.tsv"), msa = rows, methods = c("Comeron", "Li"))
  res2 <- coal_dnds_clusters(cds, fx("family_clusters.tsv"), msa = rows, methods = "Comeron",
                             comeron_flavour = "gestimator", terminal_stops = "compare")
  expect_identical(nrow(res2$pairs), sum(res$pairs$method == "Comeron"))
  expect_true(all(res$pairs$method %in% c("Comeron", "Li")))
  expect_true(all(c("Comeron", "Li") %in% res$summary$method))
  expect_gt(sum(!is.na(res$pairs$dN)), 10)
})

test_that("Comeron's two-position pathways (comeron_twosubs) and the minimum class size (comeron_min_sites)", {
  set.seed(20261114)
  bases <- c("T", "C", "A", "G")
  sense <- setdiff(apply(expand.grid(bases, bases, bases)[, 3:1], 1, paste0, collapse = ""), c("TAA", "TAG", "TGA"))
  mutate2 <- function(c) {
    repeat {
      d <- strsplit(c, "")[[1]]
      for (p in sample(3, 2)) d[p] <- sample(setdiff(bases, d[p]), 1)
      d <- paste0(d, collapse = "")
      if (!d %in% c("TAA", "TAG", "TGA")) return(d)
    }
  }
  pair <- function(n) {
    x <- sample(sense, n, replace = TRUE)
    y <- vapply(x, function(c) if (runif(1) < 0.3) mutate2(c) else c, "")
    c(paste0(x, collapse = ""), paste0(y, collapse = ""))
  }
  p <- replicate(40, pair(150))
  ab <- coal_kaks(p[1, ], p[2, ], method = "Comeron")
  ba <- coal_kaks(p[2, ], p[1, ], method = "Comeron")
  # the default (comeron's pathways) does not depend on which sequence is first; libsequence's branches do
  expect_identical(ab$dN, ba$dN)
  expect_identical(ab$dS, ba$dS)
  lab <- coal_kaks(p[1, ], p[2, ], method = "Comeron", comeron_twosubs = "libsequence")
  lba <- coal_kaks(p[2, ], p[1, ], method = "Comeron", comeron_twosubs = "libsequence")
  expect_gt(sum(lab$dN != lba$dN), 20)
  # NULL = by the flavour: comeron with "coal", libsequence with "gestimator"
  expect_identical(coal_kaks(p[1, ], p[2, ], method = "Comeron", comeron_twosubs = "comeron"), ab)
  expect_identical(coal_kaks(p[1, ], p[2, ], method = "Comeron", comeron_flavour = "gestimator"),
                   coal_kaks(p[1, ], p[2, ], method = "Comeron", comeron_flavour = "gestimator",
                             comeron_twosubs = "libsequence"))
  # the pair functions take the options too (codon alignments of coal_dnds_pairs = coal_kaks on them)
  q <- coal_read_fasta(fx("query_cds.fasta"))
  s <- coal_read_fasta(fx("subject_cds.fasta"))
  for (ts in c("comeron", "libsequence")) {
    got <- coal_dnds_pairs(q, s, pairs = orthologs(), methods = "Comeron", comeron_twosubs = ts,
                           keep_alignments = TRUE)
    aln <- attr(got, "alignments")
    ok <- !is.na(aln$codon1)
    k <- coal_kaks(aln$codon1[ok], aln$codon2[ok], method = "Comeron", comeron_twosubs = ts)
    expect_identical(got$dS[ok], k$dS)
  }
  # one codon, a nonsynonymous change: the 2-fold transversional class has less than one site and an undefined
  # correction. Default (3): gestimator's value, status warning; 0: NA, saturated (the rule before 2026-11-14)
  d3 <- coal_kaks("GCT", "GAT", method = "Comeron")
  d0 <- coal_kaks("GCT", "GAT", method = "Comeron", comeron_min_sites = 0)
  expect_false(is.na(d3$dN))
  expect_identical(d3$status, "warning")
  expect_match(d3$note, "fewer than 3 sites has an undefined correction", fixed = TRUE)
  expect_true(is.na(d0$dN))
  expect_identical(d0$status, "saturated")
  # never changes a value that is not NA otherwise
  r <- refs()
  a <- coal_kaks(r$seq1, r$seq2, method = "Comeron", comeron_min_sites = 0)
  b <- coal_kaks(r$seq1, r$seq2, method = "Comeron")
  expect_identical(a$dN[!is.na(a$dN)], b$dN[!is.na(a$dN)])
  expect_identical(a$dS[!is.na(a$dS)], b$dS[!is.na(a$dS)])
  expect_gt(sum(is.na(a$dN) & !is.na(b$dN)), 0)
  expect_error(coal_kaks("GCT", "GAT", method = "Comeron", comeron_twosubs = "x"), "comeron_twosubs")
  expect_error(coal_kaks("GCT", "GAT", method = "Comeron", comeron_min_sites = -1), "comeron_min_sites")
  expect_error(coal_dnds_pairs(q, s, pairs = orthologs(), methods = "Comeron", comeron_min_sites = NA),
               "comeron_min_sites")
})
