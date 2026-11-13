#' dN/dS for pairs of coding sequences
#'
#' For every pair: translate both coding sequences (CDS), align the proteins globally (Needleman-Wunsch,
#' BLOSUM62, full length), thread the codons onto the alignment and estimate dN, dS and dN/dS. Everything runs
#' inside the R process on `threads` threads; no program is started and no file is written. The result has the
#' columns of `dnds.tsv` of `coal pipeline pairs`.
#'
#' One bad gene or pair never stops the call: the pair's row gets `status = "error"` and the reason in `note`.
#'
#' @param query_cds,subject_cds Coding sequences: a named character vector (names = IDs), a `DNAStringSet`, a list
#'   as `seqinr::read.fasta()` returns it, a two-column data frame (ID, sequence) or the path of a FASTA file.
#'   When an ID occurs twice the first record is used.
#' @param pairs Which query goes with which subject: a data frame or matrix whose first two columns are query ID
#'   and subject ID. `NULL`: the i-th query with the i-th subject.
#' @param methods Ka/Ks methods: any of `"YN"` (Yang & Nielsen 2000, the default), `"NG"`, `"LWL"`, `"MLWL"`,
#'   `"LPB"`, `"MLPB"`, `"MYN"` (KaKs_Calculator's methods; `"ALL"` for these seven), `"Comeron"` (Comeron 1995
#'   as libsequence's gestimator computes it, the default method of the R package orthologr; alias `"C95"`) and
#'   `"Li"` (Li 1993 as `seqinr::kaks()` computes it). Names are read in any case. With several methods a pair's
#'   rows come in the order NG, LWL, MLWL, LPB, MLPB, YN, MYN, Comeron, Li. Comeron and Li are implemented for
#'   the standard genetic code only (`code = 1`).
#' @param code NCBI genetic code number (1 = standard), used for translation, codon threading and Ka/Ks.
#' @param threads Worker threads; `0` = one per hardware thread. Default `getOption("coalr.threads", 1L)`.
#' @param matrix Substitution matrix: `"BLOSUM62"`, `"BLOSUM45"`, `"BLOSUM80"` or `"PAM250"`.
#' @param gap_open,gap_extend Cost of the first position of a gap and of each further position (EMBOSS needle's
#'   convention; defaults 10 and 0.5).
#' @param end_gaps `"free"`: gaps at the ends cost nothing (needle's default); `"penalized"`: as internal gaps.
#' @param preset `NULL`, `"needle"` (the defaults above) or `"pwalign"` (the scoring of
#'   `pwalign::pairwiseAlignment(type = "global", substitutionMatrix = "BLOSUM62")`: open 14, extend 4, end gaps
#'   penalised). A preset overrides `gap_open`, `gap_extend` and `end_gaps`.
#' @param max_cells Refuse (status `error`) a pair whose alignment table would have more cells than this; 0 = no
#'   limit. The table takes half a byte per cell and thread.
#' @param yn_flavour,yn_weighting Where KaKs_Calculator and PAML's yn00 differ for YN/MYN: `"kaks"` (default) or
#'   `"paml"`; `yn_weighting = FALSE` gives equal pathway weights (yn00's default control file).
#' @param comeron_flavour For `"Comeron"`: `"coal"` (default) applies COAL's policy (a terminal stop codon is
#'   dropped, columns with an internal stop are skipped) and gives `NA` with status `saturated` where the
#'   correction is undefined; `"gestimator"` gives the numbers of orthologr 0.4.3's `gestimator()` (stop codons
#'   are compared like other codons, an undefined logarithm is replaced by 0; status `warning` then), for
#'   reproducing published orthologr results.
#' @param li_flavour For `"Li"`: `"coal"` (default) applies COAL's policy; `"seqinr"` compares stop codons as
#'   `seqinr::kaks()` does. In both, seqinr's sentinel 9.999999 is `NA`.
#' @param terminal_stops Only with `comeron_flavour = "gestimator"` or `li_flavour = "seqinr"`: `"drop"` (default)
#'   ends every pair before the terminal stop codons, as everywhere else in COAL; `"compare"` takes the two
#'   terminal stop codons as the last compared codon column when both genes end in a stop codon. That is what
#'   orthologr 0.4.3 does (it translates the stop codon to `*`, aligns it and lets PAL2NAL thread it), so
#'   `"compare"` is needed to reproduce its numbers; pairs whose two stop codons differ (`TAA` against `TGA`)
#'   get another dN by it.
#' @param keep_alignments `TRUE`: also return the protein and codon alignment of every pair, as attribute
#'   `"alignments"` (a data frame: `seq1`, `seq2`, `protein1`, `protein2`, `codon1`, `codon2`).
#' @param min_codons,max_ds Flags for estimates that need care, as `coal dnds --min-codons N --max-ds X`; both
#'   off by default (`0`, `NULL`). A row that has an estimate (status `ok` or `warning`) and fewer than
#'   `min_codons` compared codons gets the note `"fewer than N codons compared"`; one with a dS above `max_ds`
#'   (a saturated synonymous distance) gets `"dS above X"`; its status becomes `warning` if it was `ok`. No row is
#'   dropped and no number changes, so the cluster summary is unaffected. Attribute `"flagged"` counts the rows
#'   each flag marked. See "Filtering before interpretation" in COAL's `docs/kaks.md`.
#'
#' @return A data frame with one row per pair and method: `seq1`, `seq2`, `method`, `dN`, `dS`, `dNdS`, `SE_dN`,
#'   `SE_dS`, `N_sites`, `S_sites`, `Nd`, `Sd`, `codons_compared`, `prot_aln_length`, `prot_identity`, `status`
#'   (`ok`, `warning`, `not_converged`, `saturated`, `no_data`, `error`) and `note`. Values that are not available
#'   are `NA` (`dNdS` is `NA` when dS is below 1e-6). `not_converged`: the YN / MYN iteration ended at its round
#'   limit without settling, which happens with very short or nearly saturated pairs; there is no estimate (`dN`,
#'   `dS`, `dNdS`, the SEs, `Nd` and `Sd` are `NA`), and `note` has the values of the last round, which are not
#'   estimates and can differ between platforms. Rows are in the order of `pairs`. Attribute `"seconds"`: time spent in
#'   alignment + codon threading and in Ka/Ks.
#'
#' @examples
#' q <- c(g1 = "ATGGCTAAAGGTCTGGAACGTTAA", g2 = "ATGAAACCCGGGTTTTAG")
#' s <- c(h1 = "ATGGCGAAAGGCCTGGATCGTTAA", h2 = "ATGAAGCCAGGGTTCTAG")
#' coal_dnds_pairs(q, s, pairs = data.frame(query = c("g1", "g2"), subject = c("h1", "h2")))
#' coal_dnds_pairs(q, s, methods = c("YN", "NG"))[, 1:7]
#' @seealso [coal_dnds_clusters()], [coal_kaks()], [coal_pipeline_pairs()]
#' @export
coal_dnds_pairs <- function(query_cds, subject_cds, pairs = NULL, methods = "YN", code = 1L,
                            threads = getOption("coalr.threads", 1L),
                            matrix = "BLOSUM62", gap_open = 10, gap_extend = 0.5, end_gaps = c("free", "penalized"),
                            preset = NULL, max_cells = 4e9, yn_flavour = c("kaks", "paml"), yn_weighting = TRUE,
                            comeron_flavour = c("coal", "gestimator"), li_flavour = c("coal", "seqinr"),
                            terminal_stops = c("drop", "compare"),
                            keep_alignments = FALSE, min_codons = 0L, max_ds = NULL) {
  flags <- check_flags(min_codons, max_ds)
  q <- as_seqs(query_cds, "query_cds", need_names = !is.null(pairs))
  s <- as_seqs(subject_cds, "subject_cds", need_names = !is.null(pairs))
  methods <- check_methods(methods)
  code <- check_count(code, "code", 1L)
  check_method_code(methods, code)
  threads <- check_count(threads, "threads", 0L)
  ap <- align_params(matrix, gap_open, gap_extend, end_gaps, max_cells, preset)
  yn_flavour <- match.arg(yn_flavour)
  comeron_flavour <- match.arg(comeron_flavour)
  li_flavour <- match.arg(li_flavour)
  terminal_stops <- match.arg(terminal_stops)
  check_flag(yn_weighting, "yn_weighting")
  check_flag(keep_alignments, "keep_alignments")

  if (is.null(pairs)) {
    if (length(q) != length(s)) {
      stop("without `pairs`, `query_cds` and `subject_cds` must have the same length (", length(q), " and ",
           length(s), ")", call. = FALSE)
    }
    n <- length(q)
    qi <- si <- seq_len(n)
    id1 <- names(q) %||% as.character(qi)
    id2 <- names(s) %||% as.character(si)
  } else {
    if (is.matrix(pairs)) pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
    if (!is.list(pairs) || length(pairs) < 2L) {
      stop("`pairs` must be a data frame or matrix with two columns: query ID, subject ID", call. = FALSE)
    }
    id1 <- as.character(pairs[[1L]])
    id2 <- as.character(pairs[[2L]])
    if (length(id1) != length(id2)) stop("the two columns of `pairs` differ in length", call. = FALSE)
    n <- length(id1)
    qi <- match(id1, names(q))
    si <- match(id2, names(s))
  }
  nq <- length(q)
  members <- as.vector(rbind(qi - 1L, si + nq - 1L))  # interleaved; NA where the ID is unknown
  member_ids <- as.vector(rbind(id1, id2))
  raw <- coalr_dnds(c(unname(q), unname(s)), member_ids, as.integer(members), NULL,
                    as.integer(2L * (0:n)), rep.int(0L, n), integer(n + 1L), integer(0), integer(0),
                    methods, code, threads, ap, yn_flavour, yn_weighting, comeron_flavour, li_flavour,
                    terminal_stops == "compare", keep_alignments)
  out <- flag_pairs(pair_table(raw, member_ids), flags)
  attr(out, "seconds") <- raw$seconds
  if (keep_alignments) {
    odd <- seq.int(1L, by = 2L, length.out = n)
    na_empty <- function(x) ifelse(nzchar(x), x, NA_character_)
    attr(out, "alignments") <- data.frame(
      seq1 = id1, seq2 = id2,
      protein1 = na_empty(raw$protein_rows[odd]), protein2 = na_empty(raw$protein_rows[odd + 1L]),
      codon1 = na_empty(raw$codon_rows[odd]), codon2 = na_empty(raw$codon_rows[odd + 1L]),
      stringsAsFactors = FALSE)
  }
  out
}

# The pair table of coalr_dnds() as a data frame with the columns of dnds.tsv.
pair_table <- function(raw, member_ids, cluster_id = NULL) {
  p <- raw$pairs
  out <- data.frame(
    seq1 = member_ids[p$m1 + 1L], seq2 = member_ids[p$m2 + 1L], method = raw$method_names[p$method + 1L],
    dN = p$dN, dS = p$dS, dNdS = p$dNdS, SE_dN = p$SE_dN, SE_dS = p$SE_dS, N_sites = p$N_sites,
    S_sites = p$S_sites, Nd = p$Nd, Sd = p$Sd, codons_compared = p$codons_compared,
    prot_aln_length = p$prot_aln_length, prot_identity = p$prot_identity,
    status = raw$status_names[p$status + 1L], note = p$note, stringsAsFactors = FALSE)
  if (!is.null(cluster_id)) out <- cbind(data.frame(cluster_id = cluster_id, stringsAsFactors = FALSE), out)
  out
}
