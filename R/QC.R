#' Load and check BLAST results and rep-seq FASTA
#'
#' @param blast_file Path to BLAST results TSV file.
#' @param rep_fasta Path to representative sequences FASTA file.
#' @param taxonomy_col The column in BLAST file containing taxonomy strings (default "stitle").
#' @param verbose Logical; if TRUE, emit progress messages. Default FALSE.
#' @return List with BLAST dataframe (kingdom-filtered) and rep_seqs as a named list of DNA strings.
#' @importFrom utils read.table
#' @importFrom seqinr read.fasta
#' @export
load_and_check <- function(blast_file, rep_fasta, taxonomy_col = "stitle", verbose = FALSE) {
  expected_names <- c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
  )
  
  blast <- utils::read.table(
    blast_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
  )
  
  if (!all(expected_names %in% colnames(blast))) {
    blast <- utils::read.table(
      blast_file, sep = "\t", header = FALSE, col.names = expected_names,
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  
  if (taxonomy_col %in% colnames(blast)) {
    tax_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    already_parsed <- all(tax_ranks %in% colnames(blast))
    
    if (!already_parsed) {
      if (isTRUE(verbose)) message("Parsing taxonomy from '", taxonomy_col, "' column...")
      
      pattern <- "([a-z]__[^;]+)"
      tax_split <- regmatches(blast[[taxonomy_col]], gregexpr(pattern, blast[[taxonomy_col]]))
      
      tax_df <- do.call(rbind, lapply(tax_split, function(x) {
        v <- rep(NA_character_, 7)
        names(v) <- tax_ranks
        
        for (seg in x) {
          rank_code <- substr(seg, 1, 3)
          value <- sub("^[a-z]__", "", seg)
          
          if (grepl("_sp$", value) ||
              grepl("Incertae_sedis", value, ignore.case = TRUE) ||
              is.na(value) || value == "") {
            value <- "Unclassified"
          }
          
          rank <- switch(rank_code,
                         "k__" = "kingdom",
                         "p__" = "phylum",
                         "c__" = "class",
                         "o__" = "order",
                         "f__" = "family",
                         "g__" = "genus",
                         "s__" = "species",
                         NA
          )
          
          if (!is.na(rank)) {
            if (rank == "species" && !is.na(value) && value != "Unclassified") {
              value <- gsub("_", " ", value)
            }
            v[[rank]] <- value
          }
        }
        
        v
      }))
      
      tax_df <- as.data.frame(tax_df, stringsAsFactors = FALSE)
      blast <- cbind(blast, tax_df)
    } else {
      if (isTRUE(verbose)) message("Taxonomy columns already present. Skipping parsing.")
    }
  }
  
  # Drop undefined kingdom and phylum
  if (all(c("kingdom", "phylum") %in% colnames(blast))) {
    blast <- blast[
      !(is.na(blast$kingdom) | blast$kingdom == "" | blast$kingdom == "Unclassified" |
          is.na(blast$phylum)  | blast$phylum  == "" | blast$phylum  == "Unclassified"),
      ,
      drop = FALSE
    ]
  }
  
  if (!file.exists(rep_fasta)) stop("FASTA file not found.")
  if (!grepl("\\.(fa|fasta)$", rep_fasta, ignore.case = TRUE)) stop(".FASTA format required.")
  
  rep_seqs <- seqinr::read.fasta(
    file = rep_fasta, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE
  )
  
  list(blast = blast, rep_seqs = rep_seqs)
}

#' Trim BLAST alignments by minimum length
#'
#' @param blast BLAST data frame.
#' @param rep_seqs Named list/character vector of DNA strings (from seqinr::read.fasta(as.string = TRUE)).
#' @param fraction Numeric; fraction of the median rep-seq length used as the cutoff. Default 0.6.
#' @return Filtered BLAST data frame.
#' @importFrom stats median
#' @export
trim_alignments <- function(blast, rep_seqs, fraction = 0.6) {
  seq_lengths <- nchar(unlist(rep_seqs, use.names = FALSE))
  cutoff <- stats::median(seq_lengths) * fraction
  blast[blast$length > cutoff, , drop = FALSE]
}

#' Create and return alignment length histogram (ggplot object)
#'
#' @param blast BLAST data frame.
#' @param rep_seqs Named list/character vector of DNA strings (from seqinr::read.fasta(as.string = TRUE)).
#' @param cutoff_fraction Numeric; fraction of median alignment length for cutoff line. Default 0.6.
#' @return A ggplot object.
#' @importFrom stats median setNames
#' @export
plot_alignment_hist <- function(blast, rep_seqs, cutoff_fraction = 0.6) {
  median_length <- stats::median(blast$length)
  cutoff_length <- median_length * cutoff_fraction
  mean_fasta_length <- mean(nchar(unlist(rep_seqs, use.names = FALSE)))
  
  vline_data <- data.frame(
    x = c(median_length, cutoff_length, mean_fasta_length),
    label = c("Median BLAST alignment length",
              "Alignment cutoff applied",
              "Mean FASTA length"),
    color = c("red", "green", "purple"),
    stringsAsFactors = FALSE
  )
  
  ggplot2::ggplot(blast, ggplot2::aes(x = length)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    ggplot2::geom_vline(
      data = vline_data,
      ggplot2::aes(xintercept = x, color = label),
      linetype = "dashed", linewidth = 1.2, show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      name = "Reference line",
      values = stats::setNames(vline_data$color, vline_data$label)
    ) +
    ggplot2::labs(
      title = "Distribution of alignment lengths",
      x = "Alignment length (bp)",
      y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
}