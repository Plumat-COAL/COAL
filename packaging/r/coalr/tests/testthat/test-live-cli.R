# Comparisons with a COAL command line run during the test (skipped without one): set COALR_TEST_COAL to the
# `coal` program (it needs coal-cpp), and COALR_TEST_ORTHOLOGR_SEQS to orthologr's inst/seqs for the 1,000-gene set.

run_cli <- function(exe, args) {
  out <- suppressWarnings(system2(exe, args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status") %||% 0L
  if (!status %in% c(0L, 2L)) stop("coal failed (", status, "): ", paste(utils::tail(out, 5), collapse = "\n"))
  invisible(status)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

test_that("pairs equal a fresh `coal pipeline pairs --pairs` run", {
  skip_on_cran()
  exe <- coal_cli()
  skip_if(!nzchar(exe), "no coal command line with coal-cpp (COALR_TEST_COAL)")
  out <- tempfile("coalr_cli_")
  on.exit(unlink(out, recursive = TRUE))
  run_cli(exe, c("pipeline", "pairs", "--query", fx("query_cds.fasta"), "--subject", fx("subject_cds.fasta"),
                 "--pairs", fx("pairs_orthologs.tsv"), "-o", out, "-m", "YN", "-m", "MYN", "-c", "1", "-q"))
  want <- utils::read.delim(file.path(out, "dnds.tsv"), quote = "", stringsAsFactors = FALSE)
  want$note[is.na(want$note)] <- ""
  got <- coal_dnds_pairs(fx("query_cds.fasta"), fx("subject_cds.fasta"), orthologs(), methods = c("YN", "MYN"))
  expect_same_table(got, want)
})

test_that("orthologr's 1,000-gene example equals the command line (999 pairs)", {
  skip_on_cran()
  exe <- coal_cli()
  skip_if(!nzchar(exe), "no coal command line with coal-cpp (COALR_TEST_COAL)")
  files <- orthologr_1000()
  skip_if(is.null(files), "orthologr's 1,000-gene example was not found (COALR_TEST_ORTHOLOGR_SEQS)")
  skip_if(!nzchar(Sys.which("diamond")) && !nzchar(Sys.getenv("COAL_DIAMOND")), "DIAMOND is not installed")
  out <- tempfile("coalr_cli_")
  on.exit(unlink(out, recursive = TRUE))
  run_cli(exe, c("pipeline", "pairs", "--query", files[1], "--subject", files[2], "-o", out, "-m", "YN", "-m", "NG",
                 "-q"))
  want <- utils::read.delim(file.path(out, "dnds.tsv"), quote = "", stringsAsFactors = FALSE)
  want$note[is.na(want$note)] <- ""
  pairs <- utils::read.delim(file.path(out, "orthologs.tsv"), stringsAsFactors = FALSE)[, 1:2]
  expect_gt(nrow(pairs), 900)
  got <- coal_dnds_pairs(files[1], files[2], pairs, methods = c("YN", "NG"), threads = 2)
  expect_same_table(got, want)
  expect_identical(sum(got$status == "saturated" & got$method == "YN"), sum(want$status == "saturated" & want$method == "YN"))
})

test_that("the reticulate route gives the command line's tables", {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  py <- Sys.getenv("COALR_TEST_PYTHON", "")
  ok <- coal_python_available(if (nzchar(py)) py)
  skip_if(!ok, paste("Python route not available:", attr(ok, "reason")))
  res <- coal_pipeline_pairs(fx("query_cds.fasta"), fx("subject_cds.fasta"), methods = "ALL", threads = 1,
                             pairs = fx("pairs_orthologs.tsv"))
  if (res$exit_code == 1L) skip("the pipeline could not run here")
  expect_same_table(res$pairs, read_expected("pairs_dnds.tsv"), tolerance = 1e-12)
  expect_identical(res$exit_code, 2L)          # the fixture has genes that cannot be used
  expect_true(all(c("genes", "failed", "report") %in% names(res)))
  expect_identical(res$report$schema, "coal.pipeline/1")
  inproc <- coal_dnds_pairs(fx("query_cds.fasta"), fx("subject_cds.fasta"), orthologs(), methods = "ALL")
  expect_same_table(inproc, res$pairs)

  cl <- coal_pipeline_clusters(fx("family_cds.fasta"), clusters = fx("family_clusters.tsv"), threads = 1,
                               max_cluster_size = 5L)
  expect_same_table(cl$pairs, read_expected("clusters_rep_dnds.tsv"), tolerance = 1e-12)
  expect_same_table(cl$summary, read_expected("clusters_rep_summary.tsv"), tolerance = 1e-12)
  dir <- tempfile("coalr_keep_")
  on.exit(unlink(dir, recursive = TRUE))
  kept <- coal_pipeline_clusters(fx("family_cds.fasta"), clusters = fx("family_clusters.tsv"), out_dir = dir,
                                 threads = 1, max_cluster_size = 5L, large_clusters = "skip")
  expect_true(file.exists(file.path(dir, "cluster_summary.tsv")))
  expect_identical(kept$out_dir, dir)
})

test_that("the Python route is refused cleanly when it is not there", {
  skip_if_not_installed("reticulate")
  skip_if(isTRUE(coal_python_available()), "Python with coal-align is available")
  expect_false(coal_python_available())
  expect_error(coal_pipeline_pairs("a.fa", "b.fa"), "not available")
})
