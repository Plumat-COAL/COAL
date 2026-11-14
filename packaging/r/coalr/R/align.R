#' Global protein alignment of pairs (Needleman-Wunsch, BLOSUM62)
#'
#' Aligns `a[i]` against `b[i]` with COAL's exact global aligner (affine gaps). The result is always full length:
#' every residue of both sequences is in the alignment, overhanging ends are gapped, never cut off. The defaults
#' are those of EMBOSS needle.
#'
#' @param a,b Protein sequences: character vectors (or `AAStringSet`s) of equal length; a length-one argument is
#'   used for every pair. Letters and `*`; anything else is that pair's error.
#' @inheritParams coal_dnds_pairs
#' @return A data frame, one row per pair: `aligned_a`, `aligned_b` (the gapped rows), `score`, `length`,
#'   `identities`, `similarities`, `gaps`, and `error` (`NA` unless the pair could not be aligned; its rows are
#'   then `NA`).
#' @examples
#' coal_align_pairs("MAKGLERWW", "MAKGLDR")
#' coal_align_pairs(c("MKTAYIAK", "MSSHEGGK"), c("MKTAIAK", "MSHEGGKK"), preset = "pwalign")
#' @export
coal_align_pairs <- function(a, b, matrix = "BLOSUM62", gap_open = 10, gap_extend = 0.5,
                             end_gaps = c("free", "penalized"), preset = NULL, tie_break = NULL, max_cells = 4e9,
                             threads = getOption("coalr.threads", 1L)) {
  a <- as_seqs(a, "a", need_names = FALSE)
  b <- as_seqs(b, "b", need_names = FALSE)
  if (length(a) == 1L && length(b) > 1L) a <- rep(a, length(b))
  if (length(b) == 1L && length(a) > 1L) b <- rep(b, length(a))
  if (length(a) != length(b)) {
    stop("`a` and `b` must have the same length (", length(a), " and ", length(b), ")", call. = FALSE)
  }
  threads <- check_count(threads, "threads", 0L)
  ap <- align_params(matrix, gap_open, gap_extend, end_gaps, max_cells, preset, tie_break)
  as.data.frame(coalr_align(unname(a), unname(b), ap, threads), stringsAsFactors = FALSE)
}
