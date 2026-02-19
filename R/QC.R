#' Load and check BLAST results and rep-seq FASTA
#'
#' @param blast_file Path to BLAST results TSV file
#' @param rep_fasta Path to representative sequences FASTA file
#' @return List with BLAST dataframe (kingdom-filtered) and rep_seqs DNAStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange
#' @export
load_and_check <- function(blast_file, rep_fasta, taxonomy_col = "stitle") {
  expected_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle")
  blast <- read.table(blast_file, sep = "\t", header = TRUE)
  if (!all(expected_names %in% colnames(blast))) {
    blast <- read.table(blast_file, sep = "\t", header = FALSE, col.names = expected_names)
  }
  
  if (taxonomy_col %in% colnames(blast)) {
    tax_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    already_parsed <- all(tax_ranks %in% colnames(blast))
    if (!already_parsed) {
      cat("Parsing taxonomy from", taxonomy_col, "column...\n")
      pattern <- "([a-z]__[^;]+)"
      tax_split <- regmatches(blast[[taxonomy_col]], gregexpr(pattern, blast[[taxonomy_col]]))
      tax_df <- do.call(rbind, lapply(tax_split, function(x) {
        v <- rep(NA, 7)
        names(v) <- tax_ranks
        for (seg in x) {
          rank_code <- substr(seg, 1, 3)
          value <- sub("^[a-z]__", "", seg)
          if (grepl("_sp$", value) || grepl("Incertae_sedis", value, ignore.case=TRUE) || is.na(value) || value == "") value <- "Unclassified"
          rank <- switch(rank_code,
                         "k__"="kingdom",
                         "p__"="phylum",
                         "c__"="class",
                         "o__"="order",
                         "f__"="family",
                         "g__"="genus",
                         "s__"="species",
                         NA)
          if (rank == "species" && !is.na(value) && value != "Unclassified") {
            value <- gsub("_", " ", value)
          }
          if (!is.na(rank)) v[rank] <- value
        }
        return(v)
      }))
      tax_df <- as.data.frame(tax_df, stringsAsFactors=FALSE)
      blast <- cbind(blast, tax_df)
    } else {
      cat("Taxonomy columns already present. Skipping parsing.\n")
    }
  }
  
  # Drop undefined kingdom and phylum
  if (all(c("kingdom", "phylum") %in% colnames(blast))) {
    blast <- blast[
      !(is.na(blast$kingdom) | blast$kingdom == "" | blast$kingdom == "Unclassified" |
          is.na(blast$phylum)  | blast$phylum  == "" | blast$phylum  == "Unclassified"),
    ]
  }
  
  if (!file.exists(rep_fasta)) stop("FASTA file not found.")
  if (!grepl("\\.(fa|fasta)$", rep_fasta, ignore.case = TRUE)) stop(".FASTA format required.")
  rep_seqs <- Biostrings::readDNAStringSet(rep_fasta)
  list(blast = blast, rep_seqs = rep_seqs)
}

#' Trim BLAST alignments by minimum length (default fraction: 0.6 of median rep-seq length)
#'
#' @param blast BLAST dataframe
#' @param rep_seqs DNAStringSet
#' @param fraction Numeric fraction of median rep-seq length (default: 0.6)
#' @return Filtered BLAST dataframe
#' @export
trim_alignments <- function(blast, rep_seqs, fraction = 0.6) {
  med_len <- median(nchar(as.character(rep_seqs)))
  cutoff <- med_len * fraction
  dplyr::filter(blast, length > cutoff)
}

#' Create and return alignment length histogram (ggplot object)
#'
#' @param blast BLAST dataframe
#' @param rep_seqs DNAStringSet
#' @param cutoff_fraction Numeric, fraction of median for cutoff (default: 0.6)
#' @return ggplot object (for later use in PDF or display)
#' @export
plot_alignment_hist <- function(blast, rep_seqs, cutoff_fraction = 0.6) {
  median_length <- median(blast$length)
  cutoff_length <- median_length * cutoff_fraction
  mean_fasta_length <- mean(nchar(as.character(rep_seqs)))
  
  vline_data <- data.frame(
    x = c(median_length, cutoff_length, mean_fasta_length),
    label = c("Median BLAST Alignment length",
              "Alignment Cutoff Applied",
              "Mean FASTA length"),
    color = c("red", "green", "purple")
  )
  
  p <- ggplot2::ggplot(blast, ggplot2::aes(x = length)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    ggplot2::geom_vline(
      data = vline_data,
      ggplot2::aes(xintercept = x, color = label),
      linetype = "dashed", linewidth = 1.2, show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      name = "Reference Line",
      values = setNames(vline_data$color, vline_data$label)
    ) +
    ggplot2::labs(
      title = "Distribution of Alignment Lengths in User Provided BLAST Results",
      x = "Alignment Length (bp)",
      y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )
  
  return(p)
}
