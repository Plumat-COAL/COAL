# The reticulate route: COAL's complete pipeline (DIAMOND search or DeepClust, Clustal Omega, the engine) through
# the Python package `coal-align`. Optional: needs reticulate, a Python with coal-align and its coal-cpp binary.

#' Is the Python route available?
#'
#' `TRUE` when \pkg{reticulate} is installed and its Python can import `coal.dnds_pipeline` (the Python package
#' `coal-align`, 0.6 or later). Never errors.
#'
#' Which Python: the one \pkg{reticulate} already runs; else `python`; else the one reticulate finds
#' (`RETICULATE_PYTHON`, an active virtual environment, `reticulate::use_virtualenv()`, the `PATH`). coalr does not
#' let reticulate download a Python for this check unless you asked for one with
#' `reticulate::py_require("coal-align")`.
#'
#' @param python Optional path of the Python (or of a virtual environment's Python) to use; passed to
#'   `reticulate::use_python()` before Python is first started.
#' @return `TRUE` or `FALSE`, with attribute `"reason"` when `FALSE`.
#' @examples
#' coal_python_available()
#' @export
coal_python_available <- function(python = NULL) {
  no <- function(reason) structure(FALSE, reason = reason)
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(no("the R package 'reticulate' is not installed: install.packages(\"reticulate\")"))
  }
  ok <- tryCatch({
    if (!reticulate::py_available(initialize = FALSE)) {
      if (!is.null(python)) {
        reticulate::use_python(python, required = TRUE)
      } else if (is.na(Sys.getenv("RETICULATE_USE_MANAGED_VENV", NA)) && !py_required_coal()) {
        # Only a Python that is already there: do not let reticulate download one just to answer this question.
        Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")
        on.exit(Sys.unsetenv("RETICULATE_USE_MANAGED_VENV"), add = TRUE)
      }
    }
    found <- withCallingHandlers(
      reticulate::py_module_available("coal.dnds_pipeline"),
      warning = function(w) {
        # reticulate's note that its default requirement (numpy) is not in this environment: COAL does not need it
        if (grepl("py_require", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
      })
    if (!found) {
      stop("this Python has no coal.dnds_pipeline: install COAL 0.6 or later (a platform wheel from https://github.com/plumat-coal/COAL/releases; docs/install.md)")
    }
    TRUE
  }, error = function(e) no(conditionMessage(e)))
  ok
}

# Did the user ask reticulate for coal-align (reticulate::py_require("coal-align"))?
py_required_coal <- function() {
  req <- tryCatch(reticulate::py_require()$packages, error = function(e) character(0))
  any(grepl("^coal[-_]align", req))
}

py_pipeline <- function(python = NULL) {
  ok <- coal_python_available(python)
  if (!ok) {
    stop("COAL's Python pipeline is not available: ", attr(ok, "reason"),
         "\nUse coal_dnds_pairs() / coal_dnds_clusters() (no Python needed), or install it: ",
         "python -m pip install <the platform wheel coal_align-<version>-py3-none-<platform>.whl from ",
         "https://github.com/plumat-coal/COAL/releases/latest> (the source package has no coal-cpp), then ",
         "reticulate::use_python() / use_virtualenv().", call. = FALSE)
  }
  reticulate::import("coal.dnds_pipeline", convert = TRUE, delay_load = FALSE)
}

read_tsv <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  utils::read.delim(path, sep = "\t", quote = "", na.strings = "NA", comment.char = "", check.names = FALSE,
                    stringsAsFactors = FALSE, colClasses = c(seq1 = "character", seq2 = "character",
                    cluster_id = "character", representative = "character", id = "character",
                    query = "character", subject = "character", note = "character", message = "character")[
                      intersect(c("seq1", "seq2", "cluster_id", "representative", "id", "query", "subject", "note",
                                  "message"),
                                strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]])])
}

pipeline_result <- function(res, out_dir, keep) {
  outputs <- res$outputs
  tables <- lapply(outputs[setdiff(names(outputs), c("report", "dnds_kaks", "codon_alignments", "protein_alignments"))], read_tsv)
  for (v in intersect(c("dnds", "cluster_summary", "failed", "unpaired"), names(tables))) {
    if (!is.null(tables[[v]]) && "note" %in% names(tables[[v]])) tables[[v]]$note[is.na(tables[[v]]$note)] <- ""
  }
  out <- list(exit_code = as.integer(res$exit_code), pairs = tables$dnds, summary = tables$cluster_summary,
              orthologs = tables$orthologs, unpaired = tables$unpaired, genes = tables$genes, failed = tables$failed,
              report = res$report, out_dir = if (keep) out_dir else NULL)
  out[!vapply(out, is.null, NA)]
}

#' The complete COAL pipeline from R, through Python
#'
#' `coal_pipeline_pairs()` runs `coal pipeline pairs` (translate, DIAMOND reciprocal best hits, protein
#' alignment, codon alignment, Ka/Ks) and `coal_pipeline_clusters()` runs `coal pipeline clusters` (DIAMOND
#' DeepClust, Clustal Omega per cluster, Ka/Ks for every pair, per-cluster summary), by calling the Python
#' functions `coal.dnds_pipeline.run_pairs()` / `run_clusters()` through \pkg{reticulate} and reading the tables
#' they write. Use them when COAL should also do the ortholog search or the clustering; when you already have
#' pairs or clusters, [coal_dnds_pairs()] and [coal_dnds_clusters()] need neither Python nor files.
#'
#' Needs: the R package \pkg{reticulate}; a Python with `coal-align` 0.6 or later installed from a platform wheel
#' (`coal_align-<version>-py3-none-<platform>.whl` of a COAL release: it includes the `coal-cpp` program; the
#' package is not on PyPI yet); DIAMOND, and for clusters Clustal Omega (on the `PATH`, or `install_tools = TRUE`).
#' [coal_python_available()] tells whether Python is ready.
#'
#' @param query,subject FASTA files of coding sequences.
#' @param cds One or more FASTA files of coding sequences.
#' @param clusters Optional clusters file (representative, member); `NULL` runs DIAMOND DeepClust.
#' @param out_dir Output folder; `NULL` uses a temporary folder that is removed afterwards.
#' @param methods Ka/Ks methods (see [coal_dnds_pairs()]).
#' @param method The same as `methods`; give one of the two.
#' @param threads Threads; `0` (the default here, as `coal pipeline -j 0`) = one per CPU the process may use. The
#'   in-process functions default to `getOption("coalr.threads", 1L)` instead.
#' @param ... Further options of `coal pipeline pairs` / `clusters`, with underscores: for example
#'   `ortho = "bh"`, `aligner = "diamond"`, `pairs = "pairs.tsv"`, `codontable = 2L`, `max_cluster_size = 100L`,
#'   `large_clusters = "sample"`, `install_tools = TRUE`, `resume = TRUE`, `longest_isoform = TRUE`,
#'   `min_codons = 50L`, `max_ds = 2`, `min_identity = 30`. The section "Command line and R arguments" of
#'   [coalr][coalr-package] has the names side by side.
#' @param python Optional Python to use (see [coal_python_available()]).
#' @return A list: `exit_code` (0 finished; 2 finished with entries in `failed`; 1 the run could not be done),
#'   `pairs` (`dnds.tsv` as a data frame), `summary` (clusters: `cluster_summary.tsv`), `orthologs` and `unpaired`
#'   (pairs), `genes`, `failed`, `report` (the JSON report as a list) and `out_dir` (when it was given). A run that
#'   could not be done is an R error with the pipeline's message.
#' @examples
#' \dontrun{
#' res <- coal_pipeline_pairs("thaliana_cds.fa", "lyrata_cds.fa", threads = 4)
#' head(res$pairs)
#' cl <- coal_pipeline_clusters("all_cds.fa", methods = c("YN", "NG"))
#' head(cl$summary)
#' }
#' @export
coal_pipeline_pairs <- function(query, subject, out_dir = NULL, methods = "YN", threads = 0L, ..., python = NULL,
                                method = NULL) {
  methods <- pick_methods(methods, method, !missing(methods), "methods", "method")
  P <- py_pipeline(python)
  run_pipeline(function(dir, opts) do.call(P$run_pairs, c(list(path.expand(query), path.expand(subject), dir),
                                                         opts)),
               out_dir, methods, threads, list(...))
}

#' @rdname coal_pipeline_pairs
#' @export
coal_pipeline_clusters <- function(cds, clusters = NULL, out_dir = NULL, methods = "YN", threads = 0L, ...,
                                   python = NULL, method = NULL) {
  methods <- pick_methods(methods, method, !missing(methods), "methods", "method")
  P <- py_pipeline(python)
  cds <- as.list(path.expand(cds))
  run_pipeline(function(dir, opts) {
    if (!is.null(clusters)) opts$clusters <- path.expand(clusters)
    do.call(P$run_clusters, c(list(cds, dir), opts))
  }, out_dir, methods, threads, list(...))
}

run_pipeline <- function(call, out_dir, methods, threads, opts) {
  keep <- !is.null(out_dir)
  if (!keep) {
    out_dir <- tempfile("coalr_pipeline_")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  }
  out_dir <- path.expand(out_dir)
  opts$method <- as.list(check_methods(methods))
  opts$jobs <- check_count(threads, "threads", 0L)
  if (is.null(opts$quiet)) opts$quiet <- TRUE
  # whole numbers must reach Python as int
  opts <- lapply(opts, function(v) if (is.numeric(v) && length(v) == 1L && v == floor(v)) as.integer(v) else v)
  res <- call(out_dir, opts)
  code <- as.integer(res$exit_code)
  if (!code %in% c(0L, 2L) || is.null(res$outputs$dnds) || !file.exists(res$outputs$dnds)) {
    stop("the COAL pipeline could not be run (exit code ", code, "): ",
         res$report$error %||% "see the report", call. = FALSE)
  }
  pipeline_result(res, out_dir, keep)
}
