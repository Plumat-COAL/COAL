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
#' @section Sequence inputs:
#' Where a function takes sequences (`cds`, `query_cds`, `protein_rows`, ...) it accepts a named character vector,
#' a `DNAStringSet` / `AAStringSet`, a list as `seqinr::read.fasta()` returns it, a two-column data frame (ID,
#' sequence), or **one string naming a FASTA file**. One unnamed string is taken as a file name when a file of
#' that name exists. When no such file exists it is still taken as a file name, and is then the error "file not
#' found" instead of being used as a sequence, if it is written like a file name and not like a sequence row:
#' \itemize{
#' \item it has a path separator (`/` or `\`), a leading `~` or a sequence-file extension (`.fa`, `.fasta`,
#'   `.fas`, `.fna`, `.ffn`, `.faa`, `.frn`, `.cds`, `.pep`, `.aa`, `.nt`, `.aln`, `.afa`, `.mfa`, `.seq`, `.txt`,
#'   each optionally followed by `.gz`, `.bz2` or `.xz`) **and** a character that no sequence row has (a digit,
#'   `_`, `:`, a space, ...): `"run2/cds"`, `"query_cds.fasta"`, `"C:\\data\\cds"`; or
#' \item it has a lower-case letter, none of `-`, `*`, `?`, and either a path separator (`"data/cds"`) or a
#'   lower-case extension after a name without a dot (`"cds.fasta"`, `"genes.fa.gz"`).
#' }
#' Any other single string is one sequence, where a function accepts sequences without names (for example
#' `coal_translate("ATGAAA")`). In particular an aligned row is a sequence whatever it ends in or contains:
#' `"MKV..LFA.FA"`, `"ACGT.NT"`, `"MK/VL"` (a frameshift mark), `"~~MKV"`; so is a lower-case row with dashes or
#' with several dots. What this cannot tell apart: a mistyped file name in upper case (`"CDS.FASTA"`) is taken
#' as a sequence, and a lower-case row whose only gap is one dot just before `fa`, `aa`, `nt`, ... (or that
#' uses `/`) is taken as a file name. To pass one sequence that could be mistaken for a file name, name it:
#' `c(id = "...")`.
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
