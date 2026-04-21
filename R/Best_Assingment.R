#' Parse taxonomy cutoffs file
#'
#' Reads and processes a taxonomy cutoffs CSV for assignment thresholds at various ranks.
#'
#' @param cutoffs_file Path to a taxonomy cutoffs CSV file. If not supplied or invalid,
#'   attempts to locate the default file in the package.
#' @return A list with two elements: \code{long}, a data frame of parsed cutoffs, and
#'   \code{ranks}, the vector of taxonomic ranks.
#' @importFrom utils read.csv
#' @export
parse_taxonomy_cutoffs <- function(cutoffs_file = NULL) {
  if (is.null(cutoffs_file) || !file.exists(cutoffs_file)) {
    cutoffs_file <- system.file("extdata", "taxonomy_cutoffs.csv", package = "ClassifyITS")
  }
  if (!nzchar(cutoffs_file) || !file.exists(cutoffs_file)) {
    stop("Could not locate taxonomy_cutoffs.csv. Please reinstall the package or supply a custom file.")
  }
  
  cutoffs_raw <- utils::read.csv(cutoffs_file, stringsAsFactors = FALSE)
  
  tax_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # rank+type -> column name in CSV
  col_map <- list(
    kingdom = list(evalue = "e.value.kingdom"),
    phylum  = list(evalue = "e.value.phylum"),
    class   = list(evalue = "e.value.class"),
    order   = list(evalue = "e.value.order"),
    family  = list(evalue = "e.value.family"),
    genus   = list(evalue = "e.value.genus", percent_identity = "per.ident.genus"),
    species = list(percent_identity = "per.ident.species")
  )
  
  long_cutoffs <- list()
  
  for (i in seq_len(nrow(cutoffs_raw))) {
    r <- cutoffs_raw[i, , drop = FALSE]
    
    for (rank in tax_ranks) {
      types_for_rank <- names(col_map[[rank]])
      
      for (type in types_for_rank) {
        col <- col_map[[rank]][[type]]
        if (!col %in% names(cutoffs_raw)) next
        
        val <- r[[col]]
        if (!is.null(val) && !is.na(val) && val != "") {
          long_cutoffs[[length(long_cutoffs) + 1]] <- data.frame(
            rank = rank,
            kingdom = r$Kingdom,
            phylum  = r$Phylum,
            class   = r$Class,
            order   = r$Order,
            family  = r$Family,
            genus   = r$Genus,
            species = r$Species,
            cutoff_type  = type,
            cutoff_value = val,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  cutoffs_long <- if (length(long_cutoffs) == 0) {
    data.frame(
      rank = character(),
      kingdom = character(),
      phylum  = character(),
      class   = character(),
      order   = character(),
      family  = character(),
      genus   = character(),
      species = character(),
      cutoff_type  = character(),
      cutoff_value = character(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, long_cutoffs)
  }
  
  list(long = cutoffs_long, ranks = tax_ranks)
}


#' Hierarchical best-hit taxonomy assignment with per-rank fallback rule
#'
#' Pass ONLY those OTUs that haven't been assigned already!
#' For each rank, if the best e-value hit is undefined and the second-best hit is defined
#' and at least 60% as good, use the second-best hit's value for that rank.
#'
#' @param blast_qc A data.frame of BLAST results for query sequences. Must include qseqid,
#'   evalue, pident, length, and taxonomy columns: kingdom/phylum/class/order/family/genus/species.
#' @param cutoffs_long Long-form cutoffs (parse_taxonomy_cutoffs()$long).
#' @param defaults Named list of default cutoff values. Required names:
#'   - kingdom/phylum/class/order/family: evalue defaults
#'   - species: percent_identity default
#'   - genus_evalue and/or genus_percent_identity defaults (see genus_cutoff_mode)
#' @param genus_cutoff_mode One of: "prefer_evalue", "prefer_pident", "both".
#' @return A data.frame containing hierarchical taxonomy assignment for each query sequence.
#' @export
best_hit_taxonomy_assignment <- function(blast_qc, cutoffs_long, defaults,
                                         genus_cutoff_mode = c("prefer_evalue", "prefer_pident", "both")) {
  genus_cutoff_mode <- match.arg(genus_cutoff_mode)
  
  tax_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # FIXED fallback threshold: "second is at least 60% as good" => allow up to 1/0.6 ≈ 1.6667x worse e-value
  F <- 1 / 0.6  # 1.6666667
  
  # helper to parse e-notation like "e-50" -> 1e-50
  parse_cutoff_numeric <- function(cv) {
    n <- suppressWarnings(as.numeric(cv))
    if (is.na(n) && is.character(cv) && grepl("^e-", cv)) n <- as.numeric(sub("^e-", "1e-", cv))
    n
  }
  
  # find most specific matching cutoff for rank+type and hit lineage
  find_cutoff <- function(rank, hit, type, defaults, cutoffs_long) {
    match_order <- tax_ranks
    sel <- cutoffs_long[cutoffs_long$rank == rank & cutoffs_long$cutoff_type == type, , drop = FALSE]
    
    for (lev in seq(length(match_order), 1)) {
      tmp_sel <- sel
      for (m in seq_len(lev)) {
        col <- match_order[m]
        v <- hit[[col]]
        tmp_sel <- tmp_sel[is.na(tmp_sel[[col]]) | tmp_sel[[col]] == "" | tmp_sel[[col]] == v, , drop = FALSE]
      }
      if (nrow(tmp_sel) > 0) return(parse_cutoff_numeric(tmp_sel$cutoff_value[1]))
    }
    
    key <- if (rank == "genus") paste0("genus_", type) else rank
    defaults[[key]]
  }
  
  # detect whether ANY matching cutoff exists (vs. falling back to default)
  has_matching_cutoff <- function(rank, hit, type, cutoffs_long) {
    match_order <- tax_ranks
    sel <- cutoffs_long[cutoffs_long$rank == rank & cutoffs_long$cutoff_type == type, , drop = FALSE]
    if (nrow(sel) == 0) return(FALSE)
    
    for (lev in seq(length(match_order), 1)) {
      tmp_sel <- sel
      for (m in seq_len(lev)) {
        col <- match_order[m]
        v <- hit[[col]]
        tmp_sel <- tmp_sel[is.na(tmp_sel[[col]]) | tmp_sel[[col]] == "" | tmp_sel[[col]] == v, , drop = FALSE]
      }
      if (nrow(tmp_sel) > 0) return(TRUE)
    }
    FALSE
  }
  
  final_assignments <- vector("list", length = 0)
  
  for (otu in unique(blast_qc$qseqid)) {
    hits <- blast_qc[blast_qc$qseqid == otu, , drop = FALSE]
    if (nrow(hits) == 0) next
    
    hits$evalue_num <- suppressWarnings(as.numeric(hits$evalue))
    hits$pident_num <- suppressWarnings(as.numeric(hits$pident))
    
    hits <- hits[order(hits$evalue_num, na.last = TRUE), , drop = FALSE]
    
    best_hit   <- as.list(hits[1, ])
    second_hit <- if (nrow(hits) > 1) as.list(hits[2, ]) else NULL
    
    best_e <- best_hit$evalue_num
    sec_e  <- if (!is.null(second_hit)) suppressWarnings(as.numeric(second_hit$evalue)) else NA_real_
    
    taxonomy <- list()
    
    for (rank in tax_ranks) {
      # start with best hit
      hit_for_rank <- best_hit
      val <- hit_for_rank[[rank]]
      
      # fallback to second hit if best has undefined label and second is within F-fold of best (fixed)
      if (is.na(val) || val == "" || val == "Unclassified") {
        if (!is.null(second_hit) && is.finite(sec_e)) {
          qualifies <- if (is.finite(best_e) && best_e > 0) {
            sec_e <= best_e * F
          } else {
            # if best is 0, only allow second if also 0
            is.finite(best_e) && best_e == 0 && sec_e == 0
          }
          
          if (qualifies) {
            second_val <- second_hit[[rank]]
            if (!is.na(second_val) && second_val != "" && second_val != "Unclassified") {
              hit_for_rank <- second_hit
              val <- second_val
            }
          }
        }
      }
      
      # apply cutoffs using hit_for_rank
      if (rank == "species") {
        cutoff <- find_cutoff("species", hit_for_rank, "percent_identity", defaults, cutoffs_long)
        pid <- suppressWarnings(as.numeric(hit_for_rank$pident))
        taxonomy[[rank]] <- if (is.na(val) || val == "" || val == "Unclassified" || is.na(pid) || pid < cutoff * 100) "Unclassified" else val
        
      } else if (rank == "genus") {
        pid <- suppressWarnings(as.numeric(hit_for_rank$pident))
        ev  <- suppressWarnings(as.numeric(hit_for_rank$evalue))
        
        has_e  <- has_matching_cutoff("genus", hit_for_rank, "evalue", cutoffs_long)
        has_pi <- has_matching_cutoff("genus", hit_for_rank, "percent_identity", cutoffs_long)
        
        use_e  <- (genus_cutoff_mode == "prefer_evalue" && has_e) ||
          (genus_cutoff_mode == "prefer_pident" && !has_pi && has_e) ||
          (genus_cutoff_mode == "both" && has_e)
        
        use_pi <- (genus_cutoff_mode == "prefer_pident" && has_pi) ||
          (genus_cutoff_mode == "prefer_evalue" && !has_e && has_pi) ||
          (genus_cutoff_mode == "both" && has_pi)
        
        pass_e  <- TRUE
        pass_pi <- TRUE
        
        if (use_e) {
          cutoff_e <- find_cutoff("genus", hit_for_rank, "evalue", defaults, cutoffs_long)
          pass_e <- !(is.na(ev) || is.na(cutoff_e) || ev > cutoff_e)
        }
        if (use_pi) {
          cutoff_pi <- find_cutoff("genus", hit_for_rank, "percent_identity", defaults, cutoffs_long)
          pass_pi <- !(is.na(pid) || is.na(cutoff_pi) || pid < cutoff_pi * 100)
        }
        
        pass <- pass_e && pass_pi
        taxonomy[[rank]] <- if (is.na(val) || val == "" || val == "Unclassified" || !pass) "Unclassified" else val
        
      } else {
        cutoff <- find_cutoff(rank, hit_for_rank, "evalue", defaults, cutoffs_long)
        ev <- suppressWarnings(as.numeric(hit_for_rank$evalue))
        taxonomy[[rank]] <- if (is.na(val) || val == "" || val == "Unclassified" || is.na(ev) || ev > cutoff) "Unclassified" else val
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