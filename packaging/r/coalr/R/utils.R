# Internal helpers: input coercion and option checking.

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Read a FASTA file into a named character vector
#'
#' A small reader for the functions of this package: one string per record, named by the header up to the first
#' white space. Plain or gzip-compressed files. Lines starting with `#` before a record are skipped.
#'
#' @param path Path of a FASTA file.
#' @return A named character vector of sequences (white space removed).
#' @examples
#' f <- tempfile(fileext = ".fasta")
#' writeLines(c(">a some text", "ATG", "AAA", ">b", "ATGCCC"), f)
#' coal_read_fasta(f)
#' @export
coal_read_fasta <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) stop("`path` must be one file name", call. = FALSE)
  if (!file.exists(path)) stop("cannot open '", path, "': no such file", call. = FALSE)
  con <- file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  hdr <- startsWith(lines, ">")
  if (!length(lines)) return(structure(character(0), names = character(0)))
  if (!hdr[1L]) stop(path, ": text before the first '>' header line: not a FASTA file", call. = FALSE)
  ids <- sub("^>\\s*", "", lines[hdr])
  ids <- sub("\\s.*$", "", ids)
  grp <- cumsum(hdr)[!hdr]
  seqs <- character(length(ids))
  if (length(grp)) {
    body <- gsub("[[:space:]]", "", lines[!hdr])
    joined <- vapply(split(body, grp), paste0, "", collapse = "", USE.NAMES = TRUE)
    seqs[as.integer(names(joined))] <- joined
  }
  names(seqs) <- ids
  seqs
}

# Sequences from what the user gave: a named character vector, a FASTA path, an XStringSet (Biostrings), a list
# (seqinr: one vector of characters per record) or a two-column data.frame (id, sequence).
as_seqs <- function(x, arg, need_names = TRUE) {
  if (inherits(x, "XStringSet") || inherits(x, "XString") || inherits(x, "BStringSet")) {
    x <- as.character(x)
  } else if (is.data.frame(x)) {
    if (ncol(x) < 2L) stop("`", arg, "`: a data frame needs two columns (ID, sequence)", call. = FALSE)
    x <- structure(as.character(x[[2L]]), names = as.character(x[[1L]]))
  } else if (is.list(x)) {
    x <- vapply(x, function(s) if (length(s) == 1L) as.character(s) else paste0(as.character(s), collapse = ""), "")
  } else if (is.character(x) && length(x) == 1L && is.null(names(x)) && !is.na(x) && nchar(x) < 4096L &&
             file.exists(x) && !dir.exists(x)) {
    x <- coal_read_fasta(x)
  } else if (is.factor(x)) {
    x <- structure(as.character(x), names = names(x))
  }
  if (!is.character(x)) {
    stop("`", arg, "` must be a character vector, a FASTA file, a list of sequences or a DNAStringSet/AAStringSet",
         call. = FALSE)
  }
  if (need_names) {
    nm <- names(x)
    if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) {
      stop("`", arg, "` needs names: the sequence IDs", call. = FALSE)
    }
  }
  x
}

check_count <- function(x, arg, min = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min || x != floor(x)) {
    stop("`", arg, "` must be one whole number >= ", min, call. = FALSE)
  }
  as.integer(x)
}

check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop("`", arg, "` must be TRUE or FALSE", call. = FALSE)
  x
}

check_methods <- function(methods) {
  if (!is.character(methods) || !length(methods) || anyNA(methods)) {
    stop("`methods` must name at least one Ka/Ks method: ", paste(coalr_info()$methods, collapse = ", "),
         call. = FALSE)
  }
  m <- toupper(methods)
  if (identical(m, "ALL")) return(coalr_info()$methods)
  bad <- setdiff(m, coalr_info()$methods)
  if (length(bad)) {
    stop("method '", bad[1L], "' is not one COAL computes; choose from ",
         paste(coalr_info()$methods, collapse = ", "), call. = FALSE)
  }
  unique(m)
}

# The alignment parameters as the C++ side reads them. preset: NULL, "needle" (the defaults) or "pwalign"
# (pwalign::pairwiseAlignment(type = "global", substitutionMatrix = "BLOSUM62")'s scoring: open 14, extend 4,
# end gaps penalised).
align_params <- function(matrix, gap_open, gap_extend, end_gaps, max_cells, preset = NULL) {
  end_gaps <- match.arg(end_gaps, c("free", "penalized"))
  if (!is.null(preset)) {
    preset <- match.arg(preset, c("needle", "pwalign"))
    if (preset == "pwalign") {
      gap_open <- 14
      gap_extend <- 4
      end_gaps <- "penalized"
    } else {
      gap_open <- 10
      gap_extend <- 0.5
      end_gaps <- "free"
    }
  }
  if (!is.character(matrix) || length(matrix) != 1L || is.na(matrix)) {
    stop("`matrix` must be one of ", paste(coalr_info()$matrices, collapse = ", "), call. = FALSE)
  }
  for (v in list(gap_open, gap_extend, max_cells)) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop("`gap_open`, `gap_extend` and `max_cells` must be single numbers", call. = FALSE)
    }
  }
  list(matrix = matrix, gap_open = as.numeric(gap_open), gap_extend = as.numeric(gap_extend),
       free_end_gaps = end_gaps == "free", max_cells = as.numeric(max_cells))
}

default_threads <- function() getOption("coalr.threads", 1L)

PAIR_COLUMNS <- c("seq1", "seq2", "method", "dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "Nd", "Sd",
                  "codons_compared", "prot_aln_length", "prot_identity", "status", "note")
STATS <- c("median", "mean", "q25", "q75", "min", "max")
SUMMARY_COLUMNS <- c("cluster_id", "representative", "n_members", "n_used", "method", "pair_selection", "aligner",
                     "n_pairs", "n_ok", "n_saturated", "n_no_data", "n_failed",
                     paste0(rep(c("dN", "dS", "dNdS"), each = 6L), "_", STATS), "n_dNdS", "status", "note")
