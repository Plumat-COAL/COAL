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

# Whether one unnamed string is meant as a file name although no such file exists. The first version said yes to
# every string with "/", "\\", a leading "~" or a sequence-file extension, and so called aligned rows files:
# "MKV..LFA.FA" (a dot-gapped row ending in the residues F, A), "MK/VL" (a frameshift mark), "~~MKV" (GCG's end
# gaps). A string is a file name here only if
#   1. it has a path separator, a leading "~" or a sequence-file extension, and a character that no sequence row
#      has (a digit, "_", ":", a space; anything but letters and - . * ~ ? / \): "C:\\data\\cds", "run2/cds",
#      "query_cds.fasta"; or
#   2. it is written like a file name and not like a row: it has a lower-case letter and none of - * ?, and either
#      a path separator ("data/cds", "~/x/cds") or a lower-case extension after a name without a dot ("cds.fasta",
#      "genes.fa.gz").
# Every other string is a sequence: rows in upper case whatever their gaps ("MKV.FA", "ACGT.NT", "MK/VL"), and
# lower-case rows with dashes or with more dots than an extension has ("ac.gt.aa").
looks_like_path <- function(x) {
  ext <- "\\.(fa|fasta|fas|fna|ffn|faa|frn|cds|pep|aa|nt|aln|afa|mfa|seq|txt)(\\.(gz|bz2|xz))?$"
  sep <- grepl("[/\\\\]", x)
  pathish <- sep || startsWith(x, "~") || grepl(ext, x, ignore.case = TRUE)
  if (!pathish) return(FALSE)
  if (grepl("[^A-Za-z.*~?/\\\\-]", x)) return(TRUE)                       # 1: no sequence row has this character
  if (!grepl("[a-z]", x) || grepl("[*?-]", x)) return(FALSE)              # an upper-case row, or one with gaps
  if (sep) return(!grepl("~", sub("^~", "", x)))                           # 2: lower case with a separator
  grepl(ext, x) && !grepl("[.~]", sub(ext, "", x))                         # 2: name.ext, the name without a dot
}

# Sequences from what the user gave: a named character vector, a FASTA path, an XStringSet (Biostrings), a list
# (seqinr: one vector of characters per record) or a two-column data.frame (id, sequence).
# One unnamed string is a FASTA file name when that file exists or when it is written like a file name and not
# like a sequence row (looks_like_path): then a missing file is an error, and the string is never taken for a
# sequence. Any other single string is one sequence.
as_seqs <- function(x, arg, need_names = TRUE) {
  if (inherits(x, "XStringSet") || inherits(x, "XString") || inherits(x, "BStringSet")) {
    x <- as.character(x)
  } else if (is.data.frame(x)) {
    if (ncol(x) < 2L) stop("`", arg, "`: a data frame needs two columns (ID, sequence)", call. = FALSE)
    x <- structure(as.character(x[[2L]]), names = as.character(x[[1L]]))
  } else if (is.list(x)) {
    x <- vapply(x, function(s) if (length(s) == 1L) as.character(s) else paste0(as.character(s), collapse = ""), "")
  } else if (is.character(x) && length(x) == 1L && is.null(names(x)) && !is.na(x) && nchar(x) < 4096L &&
             (file.exists(x) || looks_like_path(x))) {
    if (dir.exists(x)) stop("`", arg, "`: '", x, "' is a folder, not a FASTA file", call. = FALSE)
    if (!file.exists(x)) {
      stop("`", arg, "`: file not found: '", x, "' (one string written like a file name is read as a FASTA file ",
           "name; for a single sequence give a named vector, c(id = \"ATG...\"))",
           call. = FALSE)
    }
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

# The seven methods of KaKs_Calculator that COAL computes: what "ALL" means. Comeron and Li (which KaKs_Calculator
# does not have, and which exist for the standard genetic code only) are asked for by name.
kaks_calculator_methods <- c("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN")

# Method names as the engine writes them, whatever the case given ("comeron", "C95" -> "Comeron"; "li" -> "Li").
check_methods <- function(methods) {
  known <- coalr_info()$methods
  if (!is.character(methods) || !length(methods) || anyNA(methods)) {
    stop("`methods` must name at least one Ka/Ks method: ", paste(known, collapse = ", "), call. = FALSE)
  }
  m <- toupper(methods)
  if (identical(m, "ALL")) return(kaks_calculator_methods)
  m[m == "C95"] <- "COMERON"
  i <- match(m, toupper(known))
  if (anyNA(i)) {
    stop("method '", methods[is.na(i)][1L], "' is not one COAL computes; choose from ",
         paste(known, collapse = ", "), call. = FALSE)
  }
  unique(known[i])
}

# Comeron and Li are implemented for the standard genetic code only, as their references (gestimator, seqinr).
check_method_code <- function(methods, code) {
  m <- intersect(c("Comeron", "Li"), methods)
  if (length(m) && code != 1L) {
    stop("method '", m[1L], "' is implemented for the standard genetic code only (code = 1)", call. = FALSE)
  }
  invisible(NULL)
}

# The alignment parameters as the C++ side reads them. preset: NULL, "needle" (the defaults) or "pwalign"
# (pwalign::pairwiseAlignment(type = "global", substitutionMatrix = "BLOSUM62")'s scoring: open 14, extend 4,
# end gaps penalised, pwalign's tie-break). tie_break: NULL (the preset's: pwalign for preset "pwalign", else
# leftmost), "leftmost" or "pwalign".
align_params <- function(matrix, gap_open, gap_extend, end_gaps, max_cells, preset = NULL, tie_break = NULL) {
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
  if (is.null(matrix)) {
    matrix <- if (identical(preset, "pwalign")) "BLOSUM62-PWALIGN" else "BLOSUM62"
  }
  if (!is.character(matrix) || length(matrix) != 1L || is.na(matrix)) {
    stop("`matrix` must be one of ", paste(coalr_info()$matrices, collapse = ", "), call. = FALSE)
  }
  for (v in list(gap_open, gap_extend, max_cells)) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop("`gap_open`, `gap_extend` and `max_cells` must be single numbers", call. = FALSE)
    }
  }
  if (is.null(tie_break)) {
    tie_break <- if (identical(preset, "pwalign")) "pwalign" else "leftmost"
  }
  if (!is.character(tie_break) || length(tie_break) != 1L || !tie_break %in% c("leftmost", "pwalign")) {
    stop("`tie_break` must be \"leftmost\" or \"pwalign\"", call. = FALSE)
  }
  list(matrix = matrix, gap_open = as.numeric(gap_open), gap_extend = as.numeric(gap_extend),
       free_end_gaps = end_gaps == "free", max_cells = as.numeric(max_cells), tie_break = tie_break)
}

default_threads <- function() getOption("coalr.threads", 1L)

# min_codons / max_ds of coal_dnds_pairs(), coal_dnds_clusters() and coal_kaks(): the flags of
# `coal dnds --min-codons N --max-ds X` (cpp/src/kaks_cli.cpp PairFlags), with the same notes. A row that has an
# estimate (status ok or warning) and fewer than min_codons compared codons, or a dS above max_ds, gets the note
# appended and status "warning" when it was "ok". No row is dropped and no number changes. Attribute "flagged":
# how many rows each flag marked.
check_flags <- function(min_codons, max_ds, min_identity = NULL) {
  min_codons <- check_count(min_codons, "min_codons", 0L)
  if (!is.null(min_identity) && (!is.numeric(min_identity) || length(min_identity) != 1L || is.na(min_identity) ||
                                 !is.finite(min_identity) || min_identity <= 0 || min_identity > 100)) {
    stop("`min_identity` must be NULL (off) or one percentage, above 0 and at most 100", call. = FALSE)
  }
  if (!is.null(max_ds) && (!is.numeric(max_ds) || length(max_ds) != 1L || is.na(max_ds) || !is.finite(max_ds) ||
                           max_ds <= 0)) {
    stop("`max_ds` must be NULL (off) or one number greater than 0", call. = FALSE)
  }
  list(min_codons = min_codons, max_ds = max_ds, min_identity = min_identity)
}

# Comeron's options (docs/kaks.md "Comeron"): the pathway order for codons that differ at two positions
# (`comeron_twosubs`: NULL = by the flavour, "comeron" with "coal", "libsequence" with "gestimator") and the
# smallest site class whose undefined correction makes dN or dS NA (`comeron_min_sites`). Returns the values the
# C++ glue takes ("" for NULL).
check_comeron <- function(comeron_twosubs, comeron_min_sites) {
  if (!is.null(comeron_twosubs) && (!is.character(comeron_twosubs) || length(comeron_twosubs) != 1L ||
                                    !comeron_twosubs %in% c("comeron", "libsequence"))) {
    stop("`comeron_twosubs` must be NULL (by `comeron_flavour`), \"comeron\" or \"libsequence\"", call. = FALSE)
  }
  if (!is.numeric(comeron_min_sites) || length(comeron_min_sites) != 1L || is.na(comeron_min_sites) ||
      !is.finite(comeron_min_sites) || comeron_min_sites < 0) {
    stop("`comeron_min_sites` must be one number, 0 or more (0: any site class can make dN or dS NA)",
         call. = FALSE)
  }
  list(twosubs = comeron_twosubs %||% "", min_sites = as.numeric(comeron_min_sites))
}

# `methods` and `method` name the same argument (the command line's -m / --method; coal_kaks() has `method`, the
# other functions `methods`): either is accepted everywhere. given: whether the function's own argument was given.
pick_methods <- function(own, other, given, own_name, other_name) {
  if (is.null(other)) return(own)
  if (given && !identical(toupper(as.character(own)), toupper(as.character(other)))) {
    stop("give the Ka/Ks methods once, as `", own_name, "` or as `", other_name, "` (both were given and differ)",
         call. = FALSE)
  }
  other
}

flag_pairs <- function(out, flags) {
  if (flags$min_codons == 0L && is.null(flags$max_ds) && is.null(flags$min_identity)) return(out)
  has <- out$status %in% c("ok", "warning")
  add <- function(rows, text) {
    if (!any(rows)) return(invisible())
    out$note[rows] <<- ifelse(nzchar(out$note[rows]), paste0(out$note[rows], "; ", text), text)
    out$status[rows & out$status == "ok"] <<- "warning"
  }
  flagged <- integer(0)
  if (flags$min_codons > 0L) {
    rows <- has & !is.na(out$codons_compared) & out$codons_compared < flags$min_codons
    add(rows, sprintf("fewer than %d codons compared", flags$min_codons))
    flagged <- c(flagged, min_codons = sum(rows))
  }
  if (!is.null(flags$max_ds)) {
    rows <- has & !is.na(out$dS) & out$dS > flags$max_ds
    add(rows, paste0("dS above ", sprintf("%.6g", flags$max_ds)))
    flagged <- c(flagged, max_ds = sum(rows))
  }
  if (!is.null(flags$min_identity) && "prot_identity" %in% names(out)) {
    rows <- has & !is.na(out$prot_identity) & out$prot_identity < flags$min_identity
    add(rows, paste0("protein identity below ", sprintf("%.6g", flags$min_identity), "%"))
    flagged <- c(flagged, min_identity = sum(rows))
  }
  attr(out, "flagged") <- flagged
  out
}

PAIR_COLUMNS <- c("seq1", "seq2", "method", "dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "Nd", "Sd",
                  "codons_compared", "prot_aln_length", "prot_identity", "status", "note")
STATS <- c("median", "mean", "q25", "q75", "min", "max")
SUMMARY_COLUMNS <- c("cluster_id", "representative", "n_members", "n_used", "method", "pair_selection", "aligner",
                     "n_pairs", "n_ok", "n_not_converged", "n_saturated", "n_no_data", "n_failed",
                     paste0(rep(c("dN", "dS", "dNdS"), each = 6L), "_", STATS), "n_dNdS", "status", "note")
