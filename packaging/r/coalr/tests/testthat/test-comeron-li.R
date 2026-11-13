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
  # the flavour reaches the engine: the CDS of the fixtures end in stop codons
  leg <- coal_dnds_pairs(q, s, pairs = orthologs(), methods = "Comeron", comeron_flavour = "gestimator")
  dfl <- got[got$method == "Comeron", ]
  expect_identical(nrow(leg), nrow(dfl))
  expect_error(coal_dnds_pairs(q, s, pairs = orthologs(), methods = "Comeron", code = 11L),
               "standard genetic code only")
  # threads do not change the result
  expect_identical(coal_dnds_pairs(q, s, pairs = orthologs(), methods = c("Comeron", "Li"), threads = 1L)$dN,
                   coal_dnds_pairs(q, s, pairs = orthologs(), methods = c("Comeron", "Li"), threads = 4L)$dN)

  cds <- coal_read_fasta(fx("family_cds.fasta"))
  msa_tbl <- utils::read.delim(fx("clusters_msa.tsv"), stringsAsFactors = FALSE)
  rows <- lapply(split(msa_tbl, msa_tbl$cluster_id), function(d) structure(d$protein_row, names = d$id))
  res <- coal_dnds_clusters(cds, fx("family_clusters.tsv"), msa = rows, methods = c("Comeron", "Li"))
  expect_true(all(res$pairs$method %in% c("Comeron", "Li")))
  expect_true(all(c("Comeron", "Li") %in% res$summary$method))
  expect_gt(sum(!is.na(res$pairs$dN)), 10)
})
