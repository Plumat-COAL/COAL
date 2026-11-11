#' coalr: codon alignment and dN/dS (Ka/Ks) with the COAL C++ engine
#'
#' In-process bindings to the C++ core of COAL (<https://github.com/Plumat-COAL/COAL>): no external programs, no
#' temporary files. The main functions:
#'
#' * [coal_dnds_pairs()]: coding sequences of ortholog pairs -> protein alignment -> codon alignment -> dN, dS,
#'   dN/dS, one row per pair and method (the columns of `coal pipeline pairs`' `dnds.tsv`).
#' * [coal_dnds_clusters()]: the same for every pair of the genes in each cluster, with a per-cluster summary
#'   (`dnds.tsv` and `cluster_summary.tsv` of `coal pipeline clusters`).
#' * [coal_kaks()]: Ka/Ks for sequences that are already codon-aligned.
#' * [coal_codon_align()]: protein alignment + coding DNA -> codon alignment (what PAL2NAL does).
#' * [coal_align_pairs()]: global protein alignment (Needleman-Wunsch, BLOSUM62).
#' * [coal_translate()]: coding DNA -> protein.
#' * [coal_pipeline_pairs()], [coal_pipeline_clusters()]: the complete COAL pipeline (DIAMOND, Clustal Omega)
#'   through Python and \pkg{reticulate}; optional.
#'
#' @section Threads and interrupts:
#' `threads` sets the number of worker threads (default: `getOption("coalr.threads", 1L)`; `0` = one per hardware
#' thread). Results do not depend on it. Worker threads never call into R; a user interrupt (Ctrl-C, Esc) is
#' noticed within about 0.1 s plus the time of the pair each thread is working on.
#'
#' @keywords internal
#' @useDynLib coalr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

#' Versions and capabilities of the engine
#'
#' @return A list: `coalr` (package version), `coal_version` (the COAL release the C++ sources come from),
#'   `methods` (Ka/Ks methods), `genetic_codes` (NCBI table numbers), `matrices` (substitution matrices) and
#'   `hardware_threads`.
#' @examples
#' coal_info()$methods
#' @export
coal_info <- function() {
  c(list(coalr = as.character(utils::packageVersion("coalr"))), coalr_info())
}
