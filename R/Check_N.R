#' Check proportion of N bases in each sequence.
#'
#' Calculates the proportion of "N" bases (ambiguous bases) in each sequence and flags if above the given threshold.
#'
#' @param rep_seqs Character vector, list (e.g., from seqinr::read.fasta(as.string=TRUE)),
#'   or (optionally) a DNAStringSet.
#' @param cutoff Numeric, percent threshold (default 1).
#' @return Data frame with columns: qseqid, N_percent, N_flag.
#' @examples
#' seqs <- c(seq1 = "ATGCNNNN", seq2 = "NNNNATGC")
#' check_N(seqs)
#' check_N(seqs, cutoff = 10)
#' @export
check_N <- function(rep_seqs, cutoff = 1) {
  # Normalize input to a named character vector of sequences
  if (is.list(rep_seqs)) {
    seqs <- unlist(rep_seqs, use.names = TRUE)
  } else {
    seqs <- as.character(rep_seqs)
  }
  
  if (is.null(names(seqs)) || anyNA(names(seqs)) || any(names(seqs) == "")) {
    names(seqs) <- paste0("seq", seq_along(seqs))
  }
  
  n_perc <- vapply(seqs, function(s) {
    s <- as.character(s)
    len <- nchar(s)
    if (is.na(len) || len == 0) return(NA_real_)
    chars <- strsplit(s, "", fixed = TRUE)[[1]]
    sum(chars %in% c("N", "n")) / len * 100
  }, numeric(1))
  
  data.frame(
    qseqid = names(seqs),
    N_percent = unname(n_perc),
    N_flag = !is.na(n_perc) & (unname(n_perc) > cutoff),
    stringsAsFactors = FALSE
  )
}
