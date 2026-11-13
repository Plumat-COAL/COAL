#' dN/dS for every pair of the genes in each cluster
#'
#' For each cluster (for example of `diamond deepclust`, or of orthologr's `deepclust()`): translate the members,
#' take or build their protein alignment, thread the codons and estimate dN, dS and dN/dS for the pairs of the
#' cluster; then summarise each cluster. The result matches `dnds.tsv` and `cluster_summary.tsv` of
#' `coal pipeline clusters`.
#'
#' @section Protein alignments:
#' Clusters of two are aligned by the built-in global aligner. For clusters of three or more:
#' * `msa`: the rows of a protein multiple alignment per cluster, when you already have them; or
#' * `aligner`: a function that gets a cluster's proteins (a named character vector) and returns the aligned rows
#'   (same length; named, or in the same order). [coal_msa_aligner()] makes one that runs Clustal Omega or MAFFT.
#'   A cluster whose aligner fails is aligned pairwise instead, and the summary says so; or
#' * neither: every pair is aligned on its own by the built-in aligner (`aligner` column `nw_pairwise`). This needs
#'   no external program, but the values then come from pairwise alignments, not from one MSA.
#'
#' A stop codon at the end of a CDS is not part of its protein, so alignment rows must not hold it.
#'
#' @section Large clusters:
#' The number of pairs grows with the square of the cluster size. A cluster with more than `max_cluster_size`
#' usable members is computed as `large_clusters` says: `"representative"` (the representative against every other
#' member), `"sample"` (`pair_sample` pairs drawn without replacement; the draw depends only on `seed` and the
#' cluster ID, and does not touch R's random number stream) or `"skip"`. Such clusters are aligned pairwise.
#' The sample is not the one `coal pipeline clusters` draws for the same seed (another generator).
#'
#' @param cds Coding sequences of all members (see `query_cds` of [coal_dnds_pairs()]).
#' @param clusters The clusters: a data frame or matrix whose first two columns are representative and member, one
#'   row per member, as `diamond deepclust` writes them (a member named in two clusters stays in its first one);
#'   the path of such a tab-separated file; or a named list (cluster ID -> member IDs, the representative first).
#' @param msa Optional protein alignment rows: a named list, cluster ID -> character vector of rows named by
#'   member ID (extra clusters and clusters without an entry are allowed).
#' @param aligner Optional function(proteins) -> aligned rows, used for clusters of three or more without `msa`
#'   rows. Called in the main R session, one cluster at a time.
#' @param max_cluster_size,large_clusters,pair_sample,seed See "Large clusters".
#' @inheritParams coal_dnds_pairs
#'
#' @return A list:
#' * `pairs`: data frame, one row per pair and method: `cluster_id` and the columns of [coal_dnds_pairs()].
#'   Clusters in the order of `clusters`, pairs 1-2, 1-3, ..., 2-3, ...
#' * `summary`: data frame, one row per cluster and method: `cluster_id`, `representative`, `n_members`, `n_used`,
#'   `method`, `pair_selection` (`all`, `representative`, `sample:K`, `skipped`, `singleton`, `none`), `aligner`,
#'   `n_pairs`, `n_ok`, `n_saturated`, `n_no_data`, `n_failed`, then median, mean, q25, q75, min and max of dN, dS
#'   and dN/dS over the pairs with an estimate (quartiles of R's type 7), `n_dNdS`, `status`, `note`. A cluster
#'   without pairs has one row with `method` `NA`.
#' * `failed`: data frame of what could not be used (`level`, `cluster_id`, `seq1`, `seq2`, `code`, `message`).
#' * `alignments` (with `keep_alignments = TRUE`): data frame `cluster_id`, `id`, `protein`, `codon` for clusters
#'   aligned as a whole.
#'
#' @examples
#' cds <- c(a = "ATGGCTAAAGGTCTGGAACGTTAA", b = "ATGGCGAAAGGCCTGGATCGTTAA",
#'          c = "ATGGCTAAGGGTCTGGAACGATAA", d = "ATGTTTCCCTAA")
#' cl <- data.frame(representative = c("a", "a", "a", "d"), member = c("a", "b", "c", "d"))
#' res <- coal_dnds_clusters(cds, cl)
#' res$pairs[, 1:7]
#' res$summary[, c("cluster_id", "n_members", "n_pairs", "dNdS_median")]
#' @seealso [coal_dnds_pairs()], [coal_msa_aligner()], [coal_pipeline_clusters()]
#' @export
coal_dnds_clusters <- function(cds, clusters, msa = NULL, aligner = NULL, methods = "YN", code = 1L,
                               threads = getOption("coalr.threads", 1L),
                               max_cluster_size = 200L, large_clusters = c("representative", "sample", "skip"),
                               pair_sample = 1000L, seed = 1L,
                               matrix = "BLOSUM62", gap_open = 10, gap_extend = 0.5,
                               end_gaps = c("free", "penalized"), preset = NULL, max_cells = 4e9,
                               yn_flavour = c("kaks", "paml"), yn_weighting = TRUE, keep_alignments = FALSE,
                               comeron_flavour = c("coal", "gestimator"), li_flavour = c("coal", "seqinr")) {
  cds <- as_seqs(cds, "cds")
  cl <- as_clusters(clusters)
  methods <- check_methods(methods)
  code <- check_count(code, "code", 1L)
  check_method_code(methods, code)
  threads <- check_count(threads, "threads", 0L)
  cap <- check_count(max_cluster_size, "max_cluster_size", 2L)
  large_clusters <- match.arg(large_clusters)
  pair_sample <- check_count(pair_sample, "pair_sample", 1L)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) stop("`seed` must be one number", call. = FALSE)
  ap <- align_params(matrix, gap_open, gap_extend, end_gaps, max_cells, preset)
  yn_flavour <- match.arg(yn_flavour)
  comeron_flavour <- match.arg(comeron_flavour)
  li_flavour <- match.arg(li_flavour)
  check_flag(yn_weighting, "yn_weighting")
  check_flag(keep_alignments, "keep_alignments")
  if (!is.null(msa) && (!is.list(msa) || is.null(names(msa)))) {
    stop("`msa` must be a named list: cluster ID -> named character vector of alignment rows", call. = FALSE)
  }
  if (!is.null(aligner) && !is.function(aligner)) stop("`aligner` must be a function or NULL", call. = FALSE)
  aligner_name <- if (is.null(aligner)) NA_character_ else (attr(aligner, "name") %||% "custom")

  failed <- list()
  fail <- function(level, cid, seq1, code, message) {
    failed[[length(failed) + 1L]] <<- list(level, cid, seq1, "", code, message)
  }
  # ---- genes: the first record of an ID is used; translate once to know which members are usable
  dup <- duplicated(names(cds))
  if (any(dup)) {
    for (id in unique(names(cds)[dup])) {
      fail("gene", "", id, "duplicate_id", "the ID occurs more than once (the first record is used)")
    }
  }
  all_members <- unlist(cl, use.names = FALSE)
  idx <- match(all_members, names(cds))                # per member: position of its CDS (NA: none)
  used <- sort(unique(idx[!is.na(idx)]))
  tr <- coalr_translate(unname(cds[used]), code)
  tr_pos <- match(idx, used)                           # per member: row of tr
  gene_ok <- !is.na(idx) & tr$status[tr_pos] != "failed"
  sizes <- lengths(cl)
  ends <- cumsum(sizes)
  starts <- ends - sizes

  nc <- length(cl)
  cids <- names(cl)
  n_used <- integer(nc)
  selection <- rep("all", nc)
  aligner_col <- rep("none", nc)
  notes <- character(nc)
  # per unit (grown in lists, one entry per cluster)
  u_members <- vector("list", nc)   # pool indices (0-based), unit after unit
  u_ids <- vector("list", nc)
  u_rows <- vector("list", nc)
  u_sizes <- vector("list", nc)     # members per unit
  u_mode <- vector("list", nc)
  u_whole <- logical(nc)            # the cluster is one unit (an alignment of all its members)

  for (k in seq_len(nc)) {
    cid <- cids[k]
    span <- starts[k] + seq_len(sizes[k])
    ids_all <- all_members[span]
    ok <- gene_ok[span]
    if (!all(ok)) {
      for (m in which(!ok)) {
        if (is.na(idx[span[m]])) {
          fail("gene", cid, ids_all[m], "no_cds", "left out of the cluster: no CDS record with this ID")
        } else {
          fail("gene", cid, ids_all[m], "gene_failed",
               paste0("left out of the cluster: ", tr$note[tr_pos[span[m]]]))
        }
      }
    }
    ids <- ids_all[ok]
    pool <- idx[span][ok] - 1L
    n <- length(ids)
    n_used[k] <- n
    note <- character(0)
    if (length(ids_all) != n) note <- c(note, sprintf("%d member(s) left out (see failed)", length(ids_all) - n))
    pairwise <- NULL   # a two-row matrix of member numbers (from 1) when the cluster is aligned pair by pair
    if (n < 2L) {
      selection[k] <- if (length(ids_all) < 2L) "singleton" else "none"
    } else if (n > cap) {
      if (large_clusters == "skip") {
        selection[k] <- "skipped"
        note <- c(note, sprintf("%d members: above max_cluster_size %d, not computed", n, cap))
      } else {
        aligner_col[k] <- "nw_pairwise"
        if (large_clusters == "representative") {
          selection[k] <- "representative"
          pairwise <- rbind(1L, 2:n)
          if (ids[1L] != cid) {
            note <- c(note, sprintf("the representative is not usable; %s stands in for it", ids[1L]))
          }
        } else {
          pairwise <- sample_pairs(n, pair_sample, seed, cid)
          selection[k] <- sprintf("sample:%d", ncol(pairwise))
        }
        note <- c(note, sprintf("%d members: above max_cluster_size %d", n, cap))
      }
    } else if (n == 2L) {
      aligner_col[k] <- "nw"
      u_members[[k]] <- pool
      u_ids[[k]] <- ids
      u_rows[[k]] <- c("", "")
      u_sizes[[k]] <- 2L
      u_mode[[k]] <- 0L
      u_whole[k] <- TRUE
    } else {
      rows <- NULL
      label <- NA_character_
      given <- if (!is.null(msa)) msa[[cid]] else NULL
      if (!is.null(given)) {
        given <- as_seqs(given, "msa", need_names = FALSE)
        rows <- if (!is.null(names(given))) unname(given[ids]) else if (length(given) == length(ids_all)) given[ok]
        if (is.null(rows) || anyNA(rows)) {
          rows <- NULL
          msg <- "the rows given in `msa` do not cover the cluster's members; pairs aligned pairwise"
          note <- c(note, msg)
          fail("cluster", cid, "", "msa_failed", msg)
        } else {
          label <- "msa"
        }
      } else if (!is.null(aligner)) {
        prot <- structure(tr$protein[tr_pos[span][ok]], names = ids)
        res <- tryCatch(check_rows(aligner(prot), prot), error = function(e) e)
        if (inherits(res, "error")) {
          msg <- sprintf("%s failed (%s); pairs aligned pairwise", aligner_name,
                         gsub("[\t\r\n]+", " ", conditionMessage(res)))
          note <- c(note, msg)
          fail("cluster", cid, "", "msa_failed", msg)
        } else {
          rows <- res
          label <- aligner_name
        }
      }
      if (is.null(rows)) {
        aligner_col[k] <- "nw_pairwise"
        pairwise <- utils::combn(n, 2L)
      } else {
        aligner_col[k] <- label
        u_members[[k]] <- pool
        u_ids[[k]] <- ids
        u_rows[[k]] <- rows
        u_sizes[[k]] <- n
        u_mode[[k]] <- 0L
        u_whole[k] <- TRUE
      }
    }
    if (!is.null(pairwise) && ncol(pairwise)) {
      np <- ncol(pairwise)
      u_members[[k]] <- pool[as.vector(pairwise)]
      u_ids[[k]] <- ids[as.vector(pairwise)]
      u_rows[[k]] <- rep("", 2L * np)
      u_sizes[[k]] <- rep(2L, np)
      u_mode[[k]] <- rep(0L, np)
    }
    notes[k] <- paste(note, collapse = "; ")
  }

  unit_sizes <- unlist(u_sizes, use.names = FALSE) %||% integer(0)
  unit_cluster <- rep.int(seq_len(nc), lengths(u_sizes))
  members <- as.integer(unlist(u_members, use.names = FALSE) %||% integer(0))
  member_ids <- as.character(unlist(u_ids, use.names = FALSE) %||% character(0))
  rows <- as.character(unlist(u_rows, use.names = FALSE) %||% character(0))
  nu <- length(unit_sizes)
  raw <- coalr_dnds(unname(cds), member_ids, members, rows, as.integer(c(0L, cumsum(unit_sizes))),
                    as.integer(unlist(u_mode, use.names = FALSE) %||% integer(0)), integer(nu + 1L), integer(0),
                    integer(0), methods, code, threads, ap, yn_flavour, yn_weighting, comeron_flavour, li_flavour,
                    keep_alignments)

  row_cluster <- unit_cluster[raw$pairs$unit + 1L]
  pairs <- pair_table(raw, member_ids, cluster_id = cids[row_cluster])
  # ---- summary
  sm <- coalr_summary(as.integer(row_cluster - 1L), raw$pairs$method, raw$pairs$status, raw$pairs$dN, raw$pairs$dS,
                      raw$pairs$dNdS, length(raw$method_names))
  have <- unique(sm$cluster + 1L)
  none <- setdiff(seq_len(nc), have)
  k_all <- c(sm$cluster + 1L, none)
  ord <- order(k_all, c(sm$method, rep(0L, length(none))))
  pad <- function(x, fill) c(x, rep(fill, length(none)))[ord]
  k_ord <- k_all[ord]
  n_ok <- pad(sm$n_ok, 0L)
  n_pairs <- pad(sm$n_pairs, 0L)
  summary <- data.frame(
    cluster_id = cids[k_ord], representative = cids[k_ord], n_members = unname(sizes)[k_ord], n_used = n_used[k_ord],
    method = pad(raw$method_names[sm$method + 1L], NA_character_), pair_selection = selection[k_ord],
    aligner = aligner_col[k_ord], n_pairs = n_pairs, n_ok = n_ok, n_saturated = pad(sm$n_saturated, 0L),
    n_no_data = pad(sm$n_no_data, 0L), n_failed = pad(sm$n_failed, 0L), stringsAsFactors = FALSE)
  for (v in paste0(rep(c("dN", "dS", "dNdS"), each = 6L), "_", STATS)) summary[[v]] <- pad(sm[[v]], NA_real_)
  summary$n_dNdS <- pad(sm$n_dNdS, 0L)
  summary$status <- ifelse(n_pairs == 0L, ifelse(selection[k_ord] == "skipped", "skipped", "ok"),
                           ifelse(n_ok > 0L, "ok", "failed"))
  summary$note <- notes[k_ord]
  rownames(summary) <- NULL

  # ---- failures of the engine (a member that could not be threaded, an alignment that does not fit)
  if (length(raw$fail_member)) {
    seen <- character(0)
    for (i in seq_along(raw$fail_member)) {
      pos <- raw$fail_member[i] + 1L
      unit <- findInterval(pos - 1L, c(0L, cumsum(unit_sizes)))
      key <- paste(unit_cluster[unit], member_ids[pos], raw$fail_code[i])
      if (key %in% seen) next
      seen <- c(seen, key)
      fail("gene", cids[unit_cluster[unit]], member_ids[pos], raw$fail_code[i], raw$fail_message[i])
    }
  }
  bad_units <- which(nzchar(raw$unit_error) & u_whole[unit_cluster])
  for (u in bad_units) fail("cluster", cids[unit_cluster[u]], "", raw$unit_error_code[u], raw$unit_error[u])
  failed_df <- if (length(failed)) {
    as.data.frame(do.call(rbind, lapply(failed, function(f) unlist(f))), stringsAsFactors = FALSE)
  } else {
    as.data.frame(matrix(character(0), 0L, 6L), stringsAsFactors = FALSE)
  }
  names(failed_df) <- c("level", "cluster_id", "seq1", "seq2", "code", "message")

  out <- list(pairs = pairs, summary = summary[SUMMARY_COLUMNS], failed = failed_df)
  if (keep_alignments) {
    whole_units <- which(u_whole[unit_cluster])
    offs <- c(0L, cumsum(unit_sizes))
    pos <- unlist(lapply(whole_units, function(u) offs[u] + seq_len(unit_sizes[u])), use.names = FALSE)
    na_empty <- function(x) ifelse(nzchar(x), x, NA_character_)
    out$alignments <- data.frame(
      cluster_id = cids[unit_cluster[rep.int(whole_units, unit_sizes[whole_units])]], id = member_ids[pos],
      protein = na_empty(raw$protein_rows[pos]), codon = na_empty(raw$codon_rows[pos]), stringsAsFactors = FALSE)
  }
  attr(out, "seconds") <- raw$seconds
  out
}

# Rows returned by an aligner callback, in the order of `prot`. Errors when they do not fit.
check_rows <- function(rows, prot) {
  rows <- as_seqs(rows, "the aligner's result", need_names = FALSE)
  if (length(rows) != length(prot)) {
    stop(sprintf("the aligner returned %d rows for %d sequences", length(rows), length(prot)))
  }
  if (!is.null(names(rows))) {
    if (!all(names(prot) %in% names(rows))) stop("the aligner's rows are not named by the sequence IDs")
    rows <- rows[names(prot)]
  }
  if (anyNA(rows)) stop("the aligner returned NA")
  if (length(unique(nchar(rows))) != 1L) stop("the aligner's rows differ in length")
  unname(rows)
}

# Clusters as a named list (cluster ID -> member IDs, representative first), with DIAMOND's two-column table
# read as `coal pipeline clusters` reads it.
as_clusters <- function(x) {
  if (is.character(x) && length(x) == 1L && is.null(names(x))) {
    if (!file.exists(x)) stop("cannot open '", x, "': no such file", call. = FALSE)
    lines <- readLines(x, warn = FALSE)
    lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#")]
    parts <- strsplit(lines, "\t", fixed = TRUE)
    if (any(lengths(parts) < 2L)) {
      stop(x, ": expected two tab-separated columns (representative, member), as `diamond deepclust` writes them",
           call. = FALSE)
    }
    x <- data.frame(rep = trimws(vapply(parts, `[`, "", 1L)), mem = trimws(vapply(parts, `[`, "", 2L)),
                    stringsAsFactors = FALSE)
  }
  if (is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (is.data.frame(x)) {
    if (ncol(x) < 2L) stop("`clusters` needs two columns: representative, member", call. = FALSE)
    rep <- as.character(x[[1L]])
    mem <- as.character(x[[2L]])
    headers <- c("centroid member", "representative member", "representative_id member_id", "cluster member",
                 "cluster_id member_id", "rep member")
    if (length(rep) && paste(tolower(rep[1L]), tolower(mem[1L])) %in% headers) {
      rep <- rep[-1L]
      mem <- mem[-1L]
    }
    if (anyNA(rep) || anyNA(mem)) stop("`clusters` holds NA", call. = FALSE)
    # the representative is its cluster's first member (unless it already belongs to an earlier cluster); a member
    # named again stays where it was first
    first <- !duplicated(rep)
    key <- c(seq_along(rep), which(first) - 0.5)
    o <- order(key)
    cl <- c(rep, rep[first])[o]
    id <- c(mem, rep[first])[o]
    keep <- !duplicated(id)
    out <- split(id[keep], factor(cl[keep], levels = unique(rep)))
    return(out)
  }
  if (is.list(x) && !is.null(names(x))) {
    out <- lapply(x, as.character)
    if (anyDuplicated(names(out))) stop("`clusters`: a cluster ID occurs twice", call. = FALSE)
    return(out)
  }
  stop("`clusters` must be a two-column data frame (representative, member), a file of that table, or a named ",
       "list of member IDs", call. = FALSE)
}

# k distinct pairs (i < j) of n members as a 2 x k matrix, sorted; the same for the same (seed, cluster_id),
# without touching R's random number stream. All pairs when k >= n(n-1)/2.
sample_pairs <- function(n, k, seed, cluster_id) {
  total <- n * (n - 1) / 2
  if (k >= total) return(utils::combn(n, 2L))
  h <- as.numeric(seed) %% 2147483647
  for (ch in utf8ToInt(enc2utf8(cluster_id))) h <- (h * 31 + ch) %% 2147483647
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv(), inherits = FALSE)
      }
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  })
  suppressWarnings(set.seed(as.integer(h), kind = "Mersenne-Twister", normal.kind = "Inversion",
                            sample.kind = "Rejection"))
  t <- sort(sample.int(total, k)) - 1                    # 0-based rank in the order (1,2), (1,3), ..., (2,3), ...
  i <- floor((2 * n - 1 - sqrt((2 * n - 1)^2 - 8 * t)) / 2)
  row_start <- function(i) i * n - i * (i + 1) / 2
  i <- i - (row_start(i) > t)
  i <- i + (row_start(i + 1) <= t)
  j <- i + 1 + t - row_start(i)
  rbind(as.integer(i + 1), as.integer(j + 1))
}

#' An aligner for [coal_dnds_clusters()] that runs Clustal Omega or MAFFT
#'
#' Returns a function for the `aligner` argument of [coal_dnds_clusters()]. Unlike the rest of the package it
#' starts an external program (one process and two temporary files per cluster), so the program must be
#' installed. Internal stop codons (`*`) are aligned as `X` and put back.
#'
#' @param tool `"clustalo"` (Clustal Omega) or `"mafft"`.
#' @param path The executable; `NULL` looks on the `PATH`.
#' @param threads Threads of the program per cluster.
#' @param args Further command-line arguments (a character vector).
#' @return A function(proteins) -> aligned rows, with attribute `"name"` (the `aligner` column of the summary).
#' @examples
#' \dontrun{
#' res <- coal_dnds_clusters(cds, clusters, aligner = coal_msa_aligner("clustalo"))
#' }
#' @export
coal_msa_aligner <- function(tool = c("clustalo", "mafft"), path = NULL, threads = 1L, args = NULL) {
  tool <- match.arg(tool)
  exe <- path %||% unname(Sys.which(tool))
  if (!nzchar(exe) || !file.exists(exe)) {
    stop(tool, " was not found", if (is.null(path)) " on the PATH" else paste0(" at '", path, "'"),
         "; install it (for example `conda install -c bioconda ", if (tool == "clustalo") "clustalo" else "mafft",
         "`) or give `path`", call. = FALSE)
  }
  threads <- check_count(threads, "threads", 1L)
  fn <- function(proteins) {
    inp <- tempfile("coalr_msa_", fileext = ".fasta")
    outp <- tempfile("coalr_msa_", fileext = ".aln")
    errp <- tempfile("coalr_msa_", fileext = ".err")
    on.exit(unlink(c(inp, outp, errp)))
    n <- length(proteins)
    tmp_ids <- sprintf("s%d", seq_len(n))
    writeLines(as.vector(rbind(paste0(">", tmp_ids), gsub("*", "X", proteins, fixed = TRUE))), inp)
    if (tool == "clustalo") {
      rc <- system2(exe, c("-i", shQuote(inp), "-o", shQuote(outp), "--outfmt=fa", "--seqtype=Protein",
                           "--output-order=input-order", paste0("--threads=", threads), "--force", args),
                    stdout = FALSE, stderr = errp)
    } else {
      rc <- system2(exe, c("--auto", "--amino", "--anysymbol", "--quiet", "--thread", threads, args, shQuote(inp)),
                    stdout = outp, stderr = errp)
    }
    if (!identical(rc, 0L) || !file.exists(outp) || file.size(outp) == 0) {
      msg <- if (file.exists(errp)) paste(utils::tail(readLines(errp, warn = FALSE), 3L), collapse = " ") else ""
      stop(sprintf("exit status %s%s", rc, if (nzchar(msg)) paste0(": ", msg) else ""))
    }
    rows <- coal_read_fasta(outp)
    if (!all(tmp_ids %in% names(rows))) stop("the alignment lacks some of the sequences")
    rows <- unname(rows[tmp_ids])
    # give every row its own residues back (the program saw X for '*', and may change case)
    for (i in seq_len(n)) {
      r <- strsplit(rows[i], "", fixed = TRUE)[[1L]]
      res <- strsplit(proteins[[i]], "", fixed = TRUE)[[1L]]
      is_res <- !(r %in% c("-", "."))
      if (sum(is_res) != length(res)) stop("the alignment does not hold the residues of ", names(proteins)[i])
      r[is_res] <- res
      r[!is_res] <- "-"
      rows[i] <- paste0(r, collapse = "")
    }
    structure(rows, names = names(proteins))
  }
  structure(fn, name = tool)
}
