#' Codon alignment from a protein alignment and coding DNA (PAL2NAL replacement)
#'
#' Threads each sequence's coding DNA onto its row of a protein alignment, giving the codon alignment that
#' `pal2nal.pl pep.aln nuc.fasta -output fasta` gives, in process. The DNA may have untranslated ends, a final stop
#' codon or lower case; the protein is located in it.
#'
#' @param protein_rows The rows of the protein alignment (equal length, `-` for gaps): a character vector, an
#'   `AAStringSet`, a list or a FASTA file.
#' @param cds The coding DNA of the same sequences. When both arguments are named, DNA is matched to rows by name
#'   (every row needs one); otherwise by position.
#' @param code NCBI genetic code number.
#' @param nogap `TRUE`: remove columns with a gap or an in-frame stop codon (PAL2NAL's `-nogap`).
#' @param nomismatch `TRUE`: remove columns where protein and DNA disagree (PAL2NAL's `-nomismatch`).
#' @param strict `TRUE`: probably-wrong data (many mismatches) is an error instead of a warning.
#' @param keep_going `FALSE` (default): stop when a sequence cannot be threaded. `TRUE`: return `NA` for it, with
#'   a warning.
#' @return A character vector of codon-aligned rows, named like `protein_rows`. Attribute `"problems"`: a data
#'   frame (`id`, `type`, `message`) of mismatches and failures, when there are any.
#' @examples
#' rows <- c(a = "MAKG-LER", b = "MAKGALDR")
#' cds <- c(a = "ATGGCTAAAGGTCTGGAACGTTAA", b = "ATGGCGAAAGGCGCACTGGATCGTTAA")
#' coal_codon_align(rows, cds)
#' @export
coal_codon_align <- function(protein_rows, cds, code = 1L, nogap = FALSE, nomismatch = FALSE, strict = FALSE,
                             keep_going = FALSE) {
  rows <- as_seqs(protein_rows, "protein_rows", need_names = FALSE)
  dna <- as_seqs(cds, "cds", need_names = FALSE)
  code <- check_count(code, "code", 1L)
  if (!is.null(names(rows)) && !is.null(names(dna))) {
    miss <- setdiff(names(rows), names(dna))
    if (length(miss)) {
      stop("no coding sequence named '", miss[1L], "'", if (length(miss) > 1L) sprintf(" (and %d more)",
           length(miss) - 1L), call. = FALSE)
    }
    dna <- dna[names(rows)]
  } else if (length(rows) != length(dna)) {
    stop("`protein_rows` and `cds` must have the same length (", length(rows), " and ", length(dna),
         ") or both be named", call. = FALSE)
  }
  ids <- names(rows) %||% names(dna) %||% as.character(seq_along(rows))
  r <- coalr_codon_align(unname(rows), unname(dna), as.character(ids), code, check_flag(strict, "strict"),
                         check_flag(nogap, "nogap"), check_flag(nomismatch, "nomismatch"))
  bad <- !is.na(r$error)
  if (any(bad) && !check_flag(keep_going, "keep_going")) {
    stop(r$error[bad][1L], if (sum(bad) > 1L) sprintf(" (and %d more sequences failed)", sum(bad) - 1L),
         call. = FALSE)
  }
  out <- structure(r$aligned, names = names(rows) %||% names(dna))
  warned <- !is.na(r$warning)
  if (any(bad) || any(warned)) {
    attr(out, "problems") <- data.frame(
      id = c(ids[bad], ids[warned]), type = c(rep("error", sum(bad)), rep("warning", sum(warned))),
      message = c(r$error[bad], r$warning[warned]), stringsAsFactors = FALSE)
    if (any(bad)) {
      warning(sum(bad), " sequence(s) could not be threaded and are NA; see attr(, \"problems\")", call. = FALSE)
    }
  }
  out
}

#' Translate coding sequences
#'
#' Translates from the first base, as COAL's pipeline does: a stop codon at the end is left out, an internal stop
#' is `*`, a codon with an ambiguous base is `X` unless every reading gives the same residue, and bases after the
#' last complete codon are ignored.
#'
#' @param cds Coding sequences (see `query_cds` of [coal_dnds_pairs()]); names are optional.
#' @param code NCBI genetic code number.
#' @param details `FALSE`: return the proteins. `TRUE`: return a data frame with the columns of the pipeline's
#'   `genes.tsv` (`id`, `protein`, `nt_length`, `aa_length`, `terminal_stop`, `internal_stops`,
#'   `ambiguous_codons`, `trailing_nt`, `status`, `code`, `note`).
#' @return A character vector of proteins, named like `cds`; `NA` for a sequence that cannot be translated (empty,
#'   shorter than a codon, only a stop codon, or with a character that is no nucleotide code). Or the data frame.
#' @examples
#' coal_translate(c(a = "ATGGCTAAATAA", b = "ATGNNNGCN", c = "ATG-CC"))
#' coal_translate("ATGGCTAAATAAGG", details = TRUE)
#' @export
coal_translate <- function(cds, code = 1L, details = FALSE) {
  x <- as_seqs(cds, "cds", need_names = FALSE)
  r <- coalr_translate(unname(x), check_count(code, "code", 1L))
  if (check_flag(details, "details")) {
    return(data.frame(id = names(x) %||% as.character(seq_along(x)), protein = r$protein, nt_length = r$nt_length,
                      aa_length = r$aa_length, terminal_stop = r$terminal_stop,
                      internal_stops = r$internal_stops, ambiguous_codons = r$ambiguous_codons,
                      trailing_nt = r$trailing_nt, status = r$status, code = r$code, note = r$note,
                      stringsAsFactors = FALSE))
  }
  structure(r$protein, names = names(x))
}
