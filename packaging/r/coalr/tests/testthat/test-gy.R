# The method "GY": maximum likelihood under Goldman & Yang (1994), as PAML codeml's pairwise estimate
# (COAL cpp/src/kaks_gy.cpp, docs/kaks.md "GY: maximum likelihood").

gy_a <- "ATGGCTAAAGGTCCCGAACTGTTTACCGGCGTTGTGCCGATTCTGGTTGAACTGGATGGCGATGTTAACGGCCACAAATTC"
gy_b <- "ATGGCAAAAGGCCACGAACTATTTACCGGCGTTGTGCCGATTCTAGCTAAACTAGATGGCGATGTTAACGGCCACAAATTC"

test_that("GY is a method of coalr and coal_info() says so", {
  info <- coal_info()
  expect_true("GY" %in% info$methods)
  expect_true("gy94_ml" %in% info$features)
  expect_identical(info$methods, c("NG", "LWL", "MLWL", "LPB", "MLPB", "GY", "YN", "MYN", "Comeron", "Li"))
  # "ALL" stays KaKs_Calculator's seven counting and approximate methods: GY is asked for by name
  expect_false("GY" %in% coal_kaks(gy_a, gy_b, method = "ALL")$method)
})

test_that("GY gives PAML codeml's numbers", {
  # codeml 4.10.10, runmode = -2, CodonFreq = 2 on this pair:
  #   lnL = -135.208131; t 0.37388, kappa 4.98710, omega 0.23571; dN = 0.0584, dS = 0.2476, S = 28.4, N = 52.6
  r <- coal_kaks(gy_a, gy_b, method = c("GY", "YN"), details = TRUE)
  expect_identical(r$method, c("GY", "YN"))
  g <- r[1, ]
  expect_identical(g$status, "ok")
  expect_equal(g$dN, 0.0584, tolerance = 1e-3)
  expect_equal(g$dS, 0.2476, tolerance = 1e-3)
  expect_equal(g$dNdS, 0.23571, tolerance = 1e-4)
  expect_equal(g$kappa, 4.98710, tolerance = 1e-4)
  expect_equal(g$divergence * 3, 0.37388, tolerance = 1e-4)
  expect_equal(g$S_sites, 28.4, tolerance = 2e-3)
  expect_equal(g$N_sites, 52.6, tolerance = 2e-3)
  expect_true(is.na(g$SE_dN) && is.na(g$SE_dS))
  # names in any case, and KaKs_Calculator's name of its GY rows
  expect_identical(coal_kaks(gy_a, gy_b, method = "gy")$dS, g$dS)
})

test_that("GY: identical sequences, saturation, another genetic code, threads", {
  same <- coal_kaks(gy_a, gy_a, method = "GY")
  expect_identical(c(same$dN, same$dS), c(0, 0))
  expect_identical(same$status, "ok")
  expect_true(is.na(same$dNdS))
  # a real, saturated ortholog pair (S. cerevisiae / C. glabrata; codeml: t = 50, its bound, dS 59.9, dN 0.3613)
  x <- paste0(
      "ATGTTACGCACTACTAGACTATGGACCACCCGCATGAAATTGTTTTTGAGAAACAGCTCCGGCAATGCCCTAAACAAGAATAAACTACCATTTTTGTACT",
      "CATCCCAAGGACCTCAAGCCGTGAGGTACACTTCCCAACATGAGTGGATAGCTGTGCATCAGGACAAGACTGCCTTTGTCGGAATTACAAAATACGCCAC",
      "TGATGCCTTAGGGGACGCTACCTATGTTGAGTTGCCAGAAGTGGGCACTGAGATTGCCCAAGGTGAGTCGCTAGGGTCCATTGAGTCCGTCAAGTCAGCC",
      "TCCGAGATCTACCAGCCTGCCGATGGTACCGTAGAGGAAATTAACACTAATCTTGAGGAAAATCCAGGTGTGGTGAACGAAGATCCTATGGGTGACGGCT",
      "GGCTAGTCAAAATGAAGCTTGGTGAGGGCGTTAATGTGGAACAGGTCGAGGGTCTAATGTCCTTAGAACAGTACGAAAAGACACTGGTTCATGATGAC")
  y <- paste0(
      "GTTGTTAGATCAAGCAGGTTTGTAGCTTCTCCACTACGTTTGGCTCTCAGATGTCAATCCACAAACGCTTTAAACAAGAGCCAGTTGCCATTTACTCACT",
      "CTGGTGAAGGACCATTGGCCGTCAAGTACACTTCTCAACACGAATGGATCGCTGCACACAGTGACGGTACCGCATTCTTAGGTATAACGAAATATGCTGC",
      "CGATGCCCTAGGTGATGCCACTTACATTGAGATCCCAGAAGCTGACACTGAGGTCTCCGCCGGTGACTCTGTCTCATCTGTTGAGTCAGTCAAGTCTGCC",
      "TCCGAAGTATACACCCCATTGTCAGGTACGGTGGTGGAAGGTAACGCAGCTCTAGAAGAGTCTCCACAATTGATCAACACTGACCCACTGGGTGAAGGTT",
      "GGATCGCCAAGATAAAATTAGAGAACCCAGAAGAGATCAACACCACAGAAGGTCTACTGACTTTAGAACAATACGAGCAATCACTAGAGCACGATGAG")
  sat <- coal_kaks(x, y, method = "GY")
  expect_identical(sat$status, "saturated")
  expect_true(is.na(sat$dS) && is.na(sat$dNdS))
  expect_equal(sat$dN, 0.3613, tolerance = 1e-3)
  expect_match(sat$note, "Ks saturated")
  mt <- coal_kaks(gy_a, gy_b, method = "GY", code = 2L)
  expect_identical(mt$status, "ok")
  expect_false(isTRUE(all.equal(mt$dS, coal_kaks(gy_a, gy_b, method = "GY")$dS)))
  a <- rep(c(gy_a, x), 6)
  b <- rep(c(gy_b, y), 6)
  expect_identical(coal_kaks(a, b, method = "GY", threads = 1L), coal_kaks(a, b, method = "GY", threads = 2L))
})

test_that("coal_dnds_pairs() takes GY", {
  q <- c(g1 = gy_a)
  s <- c(h1 = gy_b)
  r <- coal_dnds_pairs(q, s, pairs = data.frame(query = "g1", subject = "h1"), methods = c("YN", "GY"))
  expect_identical(r$method, c("GY", "YN"))
  expect_identical(r$status, c("ok", "ok"))
  expect_equal(r$dS[1], 0.2476, tolerance = 1e-3)
})
