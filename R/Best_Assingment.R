#' Parse taxonomy cutoffs file
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange
#' @export
parse_taxonomy_cutoffs <- function(cutoffs_file = NULL) {
  # Robust path resolution for cutoffs file
  if (is.null(cutoffs_file) || !file.exists(cutoffs_file)) {
    cutoffs_file <- system.file("extdata", "taxonomy_cutoffs.csv", package = "ClassifyITS")
    if (!file.exists(cutoffs_file) || nchar(cutoffs_file) == 0) {
      cutoffs_file <- file.path("inst", "extdata", "taxonomy_cutoffs.csv")
    }
  }
  if (!file.exists(cutoffs_file) || nchar(cutoffs_file) == 0) {
    stop("Could not locate taxonomy_cutoffs.csv. Please reinstall the package or supply a custom file.")
  }
  cutoffs_raw <- read.csv(cutoffs_file, stringsAsFactors = FALSE)
  tax_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  long_cutoffs <- list()
  for (i in seq_len(nrow(cutoffs_raw))) {
    r <- cutoffs_raw[i,]
    for (rank in tax_ranks) {
      type <- ifelse(rank %in% c("genus", "species"), "percent_identity", "evalue")
      col <- switch(rank,
                    kingdom = "e.value.kingdom",
                    phylum  = "e.value.phylum",
                    class   = "e.value.class",
                    order   = "e.value.order",
                    family  = "e.value.family",
                    genus   = if(type == "evalue") "e.value.genus" else "per.ident.genus",
                    species = "per.ident.species"
      )
      val <- r[[col]]
      if (!is.null(val) && !is.na(val) && val != "") {
        long_cutoffs[[length(long_cutoffs)+1]] <- data.frame(
          rank = rank,
          kingdom = r$Kingdom,
          phylum = r$Phylum,
          class = r$Class,
          order = r$Order,
          family = r$Family,
          genus = r$Genus,
          species = r$Species,
          cutoff_type = type,
          cutoff_value = val,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  cutoffs_long <- do.call(rbind, long_cutoffs)
  list(long = cutoffs_long, ranks = tax_ranks)
}

#' Hierarchical best-hit taxonomy assignment with per-rank fallback rule
#'
#' Pass ONLY those OTUs that haven't been assigned already!
#' For each rank, if the best e-value hit is undefined and the second-best hit is defined
#' and at least 60% as good, use the second-best hit's value for that rank.
#'
#' @export
best_hit_taxonomy_assignment <- function(blast_qc, cutoffs_long, defaults) {
  tax_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  find_cutoff <- function(rank, hit, type, defaults, cutoffs_long) {
    match_order <- tax_ranks
    sel <- cutoffs_long[cutoffs_long$rank == rank & cutoffs_long$cutoff_type == type, ]
    for (lev in seq(length(match_order), 1)) {
      tmp_sel <- sel
      for (m in seq_len(lev)) {
        col <- match_order[m]
        v <- hit[[col]]
        tmp_sel <- tmp_sel[is.na(tmp_sel[[col]]) | tmp_sel[[col]] == "" | tmp_sel[[col]] == v, ]
      }
      if (nrow(tmp_sel) > 0) {
        cv <- tmp_sel$cutoff_value[1]
        n <- suppressWarnings(as.numeric(cv))
        if (is.na(n) && grepl("^e-", cv)) n <- as.numeric(sub("^e-", "1e-", cv))
        return(n)
      }
    }
    return(defaults[[rank]])
  }
  
  final_assignments <- list()
  for (otu in unique(blast_qc$qseqid)) {
    hits <- blast_qc[blast_qc$qseqid == otu, ]
    if (nrow(hits) == 0) next
    
    # Sort hits by evalue (ascending)
    hits <- hits[order(as.numeric(hits$evalue)), , drop=FALSE]
    best_hit <- as.list(hits[1, ])
    best_evalue <- as.numeric(best_hit$evalue)
    second_hit <- NULL
    second_evalue <- Inf
    if (nrow(hits) > 1) {
      second_hit <- as.list(hits[2, ])
      second_evalue <- as.numeric(second_hit$evalue)
    }
    
    taxonomy <- list()
    e_val <- as.numeric(best_hit$evalue)
    pident <- as.numeric(best_hit$pident)
    
    for (rank in tax_ranks) {
      # Try to assign from best hit
      val <- best_hit[[rank]]
      
      # Fallback: use second-best hit for this rank if best is undefined and qualifies
      if (is.na(val) || val == "" || val == "Unclassified") {
        use_second_this_rank <- FALSE
        if (!is.null(second_hit) && !is.na(second_evalue) && best_evalue > 0 && second_evalue <= best_evalue / 0.6) {
          second_val <- second_hit[[rank]]
          if (!is.na(second_val) && second_val != "" && second_val != "Unclassified") {
            use_second_this_rank <- TRUE
            val <- second_val
            # For percent_identity, update if genus/species
            if (rank %in% c("genus", "species")) pident <- as.numeric(second_hit$pident)
            if (rank == "kingdom") e_val <- as.numeric(second_hit$evalue)
          }
        }
      }
      
      if (rank %in% c("genus", "species")) {
        cutoff <- find_cutoff(rank, best_hit, "percent_identity", defaults, cutoffs_long)
        value <- pident
        if (is.na(val) || val == "" || val == "Unclassified" || is.na(value) || value < cutoff*100) {
          taxonomy[[rank]] <- "Unclassified"
        } else {
          taxonomy[[rank]] <- val
        }
      } else {
        cutoff <- find_cutoff(rank, best_hit, "evalue", defaults, cutoffs_long)
        value <- e_val
        if (is.na(val) || val == "" || val == "Unclassified" || is.na(value) || value > cutoff) {
          taxonomy[[rank]] <- "Unclassified"
        } else {
          taxonomy[[rank]] <- val
        }
      }
    }
    
    final_assignments[[length(final_assignments) + 1]] <- c(
      qseqid = best_hit$qseqid,
      percent_identity = best_hit$pident,
      alignment_length = best_hit$length,
      e_value = best_hit$evalue,
      taxonomy
    )
  }
  as.data.frame(do.call(rbind, lapply(final_assignments, unlist)), stringsAsFactors = FALSE)
}