fx <- function(name) testthat::test_path("fixtures", name)

read_expected <- function(name) {
  x <- utils::read.delim(fx(name), sep = "\t", quote = "", na.strings = "NA", comment.char = "",
                         stringsAsFactors = FALSE, check.names = FALSE)
  for (v in intersect(c("note"), names(x))) x[[v]][is.na(x[[v]])] <- ""
  x
}

# The command line prints 8 significant digits; coalr returns doubles.
expect_same_table <- function(got, want, tolerance = 2e-7) {
  testthat::expect_identical(names(got), names(want))
  testthat::expect_identical(nrow(got), nrow(want))
  for (v in names(want)) {
    g <- got[[v]]
    w <- want[[v]]
    if (is.numeric(w) || (is.logical(w) && all(is.na(w)))) {
      testthat::expect_identical(is.na(g), is.na(w), info = v)
      ok <- !is.na(w)
      scale <- pmax(abs(as.numeric(w[ok])), 1e-12)
      testthat::expect_true(all(abs(as.numeric(g[ok]) - as.numeric(w[ok])) / scale <= tolerance |
                                  abs(as.numeric(g[ok]) - as.numeric(w[ok])) < 1e-12), info = v)
    } else {
      testthat::expect_identical(as.character(g), as.character(w), info = v)
    }
  }
}

orthologs <- function() utils::read.delim(fx("pairs_orthologs.tsv"), stringsAsFactors = FALSE)[, 1:2]

# The `coal` command for the live comparisons: COALR_TEST_COAL, or `coal` on the PATH. "" when it cannot run the
# pipeline's engine (no coal-cpp).
coal_cli <- function() {
  exe <- Sys.getenv("COALR_TEST_COAL", unname(Sys.which("coal")))
  if (!nzchar(exe) || !file.exists(exe)) return("")
  ok <- tryCatch(suppressWarnings(system2(exe, "--cpp-path", stdout = TRUE, stderr = FALSE)), error = function(e) "")
  if (length(ok) && nzchar(ok[1L]) && is.null(attr(ok, "status"))) exe else ""
}

# orthologr's 1,000-gene example: COALR_TEST_ORTHOLOGR_SEQS, or an installed orthologr.
orthologr_1000 <- function() {
  dir <- Sys.getenv("COALR_TEST_ORTHOLOGR_SEQS", system.file("seqs", package = "orthologr"))
  f <- file.path(dir, c("ortho_thal_cds_1000.fasta", "ortho_lyra_cds_1000.fasta"))
  if (nzchar(dir) && all(file.exists(f))) f else NULL
}
