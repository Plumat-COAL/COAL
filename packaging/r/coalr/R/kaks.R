#' Ka/Ks (dN/dS) for codon-aligned sequences
#'
#' Estimates dN (Ka), dS (Ks) and dN/dS for pairs of sequences that are already aligned codon by codon (for
#' example by [coal_codon_align()] or PAL2NAL), with COAL's native engine: the methods, conventions and numbers of
#' KaKs_Calculator 2.0, without starting it.
#'
#' Codon columns with a gap, an ambiguous base or an internal stop codon in either sequence are skipped; a stop
#' codon at the end is dropped. A pair that cannot be used gets a `status` instead of an error.
#'
#' @param aligned1,aligned2 Character vectors (or `DNAStringSet`s) of equal length: the i-th elements are one
#'   pair. Each pair's two rows must have the same length, a multiple of 3; `-` or `.` for gaps.
#'   Alternatively give one alignment in `aligned1` (two or more rows) and leave `aligned2` out: every pair of its
#'   rows is computed (1-2, 1-3, ..., 2-3, ...).
#' @param method One or more of `"YN"`, `"NG"`, `"LWL"`, `"MLWL"`, `"LPB"`, `"MLPB"`, `"MYN"` (or `"ALL"` for
#'   these seven), `"Comeron"` (alias `"C95"`) and `"Li"`; see [coal_dnds_pairs()].
#' @param comeron_flavour,li_flavour Which reading of the pair `"Comeron"` and `"Li"` follow: COAL's policy
#'   (`"coal"`, the default) or the reference program's (`"gestimator"`, `"seqinr"`); see [coal_dnds_pairs()].
#' @param details `TRUE`: add KaKs_Calculator's further columns (`p_fisher`, `length`, `substitutions`, `kappa`,
#'   `divergence`, `gc`) and the counts of skipped columns.
#' @inheritParams coal_dnds_pairs
#' @return A data frame, one row per pair and method: `seq1`, `seq2` (the names, or the positions), `method`,
#'   `dN`, `dS`, `dNdS`, `SE_dN`, `SE_dS`, `N_sites`, `S_sites`, `Nd`, `Sd`, `codons_compared`, `status`, `note`.
#' @examples
#' coal_kaks("ATGGCTAAAGGTCTGGAACGT", "ATGGCGAAAGGCCTGGATCGT")
#' coal_kaks(c(a = "ATGGCTAAAGGTCTGGAACGT", b = "ATGGCGAAAGGCCTGGATCGT", c = "ATGGCTAAGGGT---GAACGA"),
#'           method = c("NG", "YN"))
#' @export
coal_kaks <- function(aligned1, aligned2 = NULL, method = "YN", code = 1L,
                      threads = getOption("coalr.threads", 1L), yn_flavour = c("kaks", "paml"),
                      yn_weighting = TRUE, details = FALSE, comeron_flavour = c("coal", "gestimator"),
                      li_flavour = c("coal", "seqinr")) {
  a <- as_seqs(aligned1, "aligned1", need_names = FALSE)
  methods <- check_methods(method)
  code <- check_count(code, "code", 1L)
  threads <- check_count(threads, "threads", 0L)
  check_method_code(methods, code)
  yn_flavour <- match.arg(yn_flavour)
  comeron_flavour <- match.arg(comeron_flavour)
  li_flavour <- match.arg(li_flavour)
  check_flag(yn_weighting, "yn_weighting")
  if (is.null(aligned2)) {
    if (length(a) < 2L) stop("give two sequences, or one alignment of two or more rows", call. = FALSE)
    ij <- utils::combn(length(a), 2L)
    id <- names(a) %||% as.character(seq_along(a))
    id1 <- id[ij[1L, ]]
    id2 <- id[ij[2L, ]]
    b <- a[ij[2L, ]]
    a <- a[ij[1L, ]]
  } else {
    b <- as_seqs(aligned2, "aligned2", need_names = FALSE)
    if (length(a) != length(b)) {
      stop("`aligned1` and `aligned2` must have the same length (", length(a), " and ", length(b), ")",
           call. = FALSE)
    }
    id1 <- names(a) %||% as.character(seq_along(a))
    id2 <- names(b) %||% as.character(seq_along(b))
  }
  r <- coalr_kaks(unname(a), unname(b), methods, code, threads, yn_flavour, yn_weighting, comeron_flavour,
                  li_flavour)
  out <- data.frame(seq1 = id1[r$pair + 1L], seq2 = id2[r$pair + 1L], method = r$method_names[r$method + 1L],
                    dN = r$dN, dS = r$dS, dNdS = r$dNdS, SE_dN = r$SE_dN, SE_dS = r$SE_dS, N_sites = r$N_sites,
                    S_sites = r$S_sites, Nd = r$Nd, Sd = r$Sd, codons_compared = r$codons_compared,
                    status = r$status_names[r$status + 1L], note = r$note, stringsAsFactors = FALSE)
  if (check_flag(details, "details")) {
    for (v in c("p_fisher", "length", "substitutions", "kappa", "divergence", "gc", "internal_stops", "gap_codons",
                "ambiguous_codons")) out[[v]] <- r[[v]]
  }
  out
}
