#' Easy taxonomy assignment for OTUs using BLAST QC output & phylum-specific thresholds.
#'
#' @param blast_filtered QC-filtered BLAST dataframe (with parsed taxonomy columns!)
#' @param cutoffs_file Path to taxonomy cutoffs CSV file. If not supplied or invalid,
#'   attempts to locate the default file in the package.
#' @param default_cutoff Default percent identity cutoff (kept for API compatibility)
#' @return List with assigned_otus_df and remaining_otus_df
#' @importFrom stats setNames
#' @export
easy_assignments <- function(
    blast_filtered,
    cutoffs_file = NULL,
    default_cutoff = 98
) {
  taxonomic_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  expected_columns <- c("qseqid", taxonomic_levels)
  assigned_otus <- list()
  easy_otus <- character()
  all_otus <- unique(blast_filtered$qseqid)
  
  # Locate cutoffs file (installed package path) [kept for compatibility]
  if (is.null(cutoffs_file) || !file.exists(cutoffs_file)) {
    cutoffs_file <- system.file("extdata", "taxonomy_cutoffs.csv", package = "ClassifyITS")
  }
  if (!nzchar(cutoffs_file) || !file.exists(cutoffs_file)) {
    stop("Could not locate taxonomy_cutoffs.csv. Please reinstall the package or supply a custom file.")
  }
  
  # Kept: read cutoffs (even if not used by the new decision rules)
  cutoffs_data <- parse_taxonomy_cutoffs(cutoffs_file)$long
  phylum_species_cutoff <- cutoffs_data[
    cutoffs_data$rank == "species" & cutoffs_data$cutoff_type == "percent_identity",
    , drop = FALSE
  ]
  phylum_species_cutoff_vec <- stats::setNames(
    as.numeric(phylum_species_cutoff$cutoff_value),
    phylum_species_cutoff$phylum
  )
  
  norm_val <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | x == "" | tolower(x) == "unclassified"] <- "Unclassified"
    x
  }
  
  # NEW: safe extractor so missing columns never create length-0 values
  best_vals_from_row <- function(row, ranks) {
    out <- vector("list", length(ranks))
    names(out) <- ranks
    for (r in ranks) {
      if (r %in% names(row)) {
        v <- row[[r]]
        out[[r]] <- if (length(v) == 0) "Unclassified" else norm_val(v)[1]
      } else {
        out[[r]] <- "Unclassified"
      }
    }
    out
  }
  
  # UPDATED: treat missing/length-0 as not assigned (prevents vapply error)
  fully_assigned <- function(vals) {
    all(vapply(taxonomic_levels, function(r) {
      v <- vals[[r]]
      if (length(v) == 0) return(FALSE)
      norm_val(v)[1] != "Unclassified"
    }, logical(1)))
  }
  
  family_supported <- function(hits_sorted) {
    fam <- norm_val(hits_sorted$family)
    fam_def <- fam[fam != "Unclassified"]
    if (length(fam_def) == 0) return(FALSE)
    
    if (nrow(hits_sorted) >= 2) {
      top2 <- fam[1:2]
      if (all(top2 != "Unclassified") && length(unique(top2)) == 1) return(TRUE)
    }
    
    tab <- table(fam_def)
    (max(tab) / sum(tab)) > 0.5
  }
  
  consensus_to_best <- function(hits_sorted, best_vals) {
    out <- best_vals
    for (r in taxonomic_levels) {
      orig <- best_vals[[r]]
      if (orig == "Unclassified") {
        out[[r]] <- "Unclassified"
        next
      }
      
      v <- norm_val(hits_sorted[[r]])
      v <- v[v != "Unclassified"]
      if (length(v) == 0) {
        out[[r]] <- "Unclassified"
        next
      }
      
      if (nrow(hits_sorted) >= 2) {
        top2 <- norm_val(hits_sorted[1:2, r, drop = TRUE])
        if (all(top2 != "Unclassified") && length(unique(top2)) == 1 && top2[1] == orig) {
          out[[r]] <- orig
          next
        }
      }
      
      tab <- table(v)
      maj <- names(tab)[which.max(tab)]
      if ((max(tab) / sum(tab)) > 0.5 && maj == orig) {
        out[[r]] <- orig
      } else {
        out[[r]] <- "Unclassified"
      }
    }
    out
  }
  
  # --- Assignment loop over OTUs (Category A/B rules) ---
  for (otu in all_otus) {
    otu_hits_all <- blast_filtered[blast_filtered$qseqid == otu, , drop = FALSE]
    if (nrow(otu_hits_all) == 0) next
    
    # fungi-only for easy assignment logic
    otu_hits <- otu_hits_all[tolower(otu_hits_all$kingdom) == "fungi", , drop = FALSE]
    if (nrow(otu_hits) == 0) next
    
    otu_hits$evalue_num <- suppressWarnings(as.numeric(otu_hits$evalue))
    otu_hits$pident_num <- suppressWarnings(as.numeric(otu_hits$pident))
    otu_hits <- otu_hits[order(otu_hits$evalue_num, na.last = TRUE), , drop = FALSE]
    
    best_hit <- otu_hits[1, , drop = FALSE]
    
    # ---------- Category A: any 100% hits ----------
    hits100 <- otu_hits[!is.na(otu_hits$pident_num) & otu_hits$pident_num == 100, , drop = FALSE]
    if (nrow(hits100) > 0) {
      hits995 <- otu_hits[!is.na(otu_hits$pident_num) & otu_hits$pident_num >= 99.5, , drop = FALSE]
      hits995 <- hits995[order(hits995$evalue_num, na.last = TRUE), , drop = FALSE]
      
      if (!family_supported(hits995)) next
      
      best_vals <- best_vals_from_row(best_hit, taxonomic_levels)
      if (!fully_assigned(best_vals)) next
      
      # No other >=99.5 hit may disagree at ANY rank
      any_disagree <- FALSE
      if (nrow(hits995) > 1) {
        others <- hits995[-1, , drop = FALSE]
        for (r in taxonomic_levels) {
          other_vals <- norm_val(others[[r]])
          if (any(other_vals != best_vals[[r]])) { any_disagree <- TRUE; break }
        }
      }
      if (any_disagree) next
      
      best_vals$qseqid <- otu
      assigned_otus[[otu]] <- best_vals
      easy_otus <- c(easy_otus, otu)
      next
    }
    
    # ---------- Category B: no 100% hits ----------
    if (is.na(best_hit$pident_num) || best_hit$pident_num < 98.5) next
    if (nrow(otu_hits) < 2) next
    
    second_hit <- otu_hits[2, , drop = FALSE]
    if (is.na(second_hit$pident_num)) next
    
    # best must be >= 0.5% better than second
    if (best_hit$pident_num < (second_hit$pident_num + 0.5)) next
    
    best_vals <- best_vals_from_row(best_hit, taxonomic_levels)
    
    ranks_to_check <- c("kingdom","phylum","class","order","family","genus")
    best_second_agree_k2g <- TRUE
    for (r in ranks_to_check) {
      if (best_vals[[r]] == "Unclassified") { best_second_agree_k2g <- FALSE; break }
      if (norm_val(best_hit[[r]]) != norm_val(second_hit[[r]])) { best_second_agree_k2g <- FALSE; break }
    }
    
    final_vals <- if (best_second_agree_k2g) best_vals else consensus_to_best(otu_hits, best_vals)
    
    if (!fully_assigned(final_vals)) next
    
    final_vals$qseqid <- otu
    assigned_otus[[otu]] <- final_vals
    easy_otus <- c(easy_otus, otu)
    next
  }
  
  # -- Build assigned_otus_df, guaranteeing all columns exist and are ordered --
  if (length(assigned_otus) == 0) {
    assigned_otus_df <- data.frame(
      qseqid = character(),
      kingdom = character(),
      phylum = character(),
      class = character(),
      order = character(),
      family = character(),
      genus = character(),
      species = character(),
      stringsAsFactors = FALSE
    )
  } else {
    assigned_otus_df <- do.call(rbind, lapply(assigned_otus, function(x) {
      vals <- stats::setNames(rep(NA, length(expected_columns)), expected_columns)
      for (nm in names(x)) {
        if (nm %in% expected_columns) vals[[nm]] <- x[[nm]]
      }
      as.data.frame(as.list(vals), stringsAsFactors = FALSE)
    }))
    rownames(assigned_otus_df) <- NULL
  }
  
  still_remaining <- setdiff(all_otus, assigned_otus_df$qseqid)
  remaining_otus_df <- blast_filtered[blast_filtered$qseqid %in% still_remaining, , drop = FALSE]
  
  # === FINAL KINGDOM CHECK (kept) ===
  bad_otus <- character()
  for (otu in assigned_otus_df$qseqid) {
    otu_hits <- blast_filtered[blast_filtered$qseqid == otu, , drop = FALSE]
    n_hits <- nrow(otu_hits)
    n_fungal <- sum(tolower(otu_hits$kingdom) == "fungi")
    if (n_hits == 0) next
    if (n_fungal / n_hits <= 0.5) bad_otus <- c(bad_otus, otu)
  }
  if (length(bad_otus) > 0) {
    assigned_otus_df <- assigned_otus_df[!assigned_otus_df$qseqid %in% bad_otus, , drop = FALSE]
    still_remaining <- union(still_remaining, bad_otus)
    remaining_otus_df <- blast_filtered[blast_filtered$qseqid %in% still_remaining, , drop = FALSE]
  }
  
  # === FINAL PHYLUM, CLASS, ORDER DOUBLE CHECK (kept; scalar fix) ===
  bad_ranks <- character()
  ranks_to_check <- c("phylum", "class", "order")
  for (otu in assigned_otus_df$qseqid) {
    otu_hits <- blast_filtered[blast_filtered$qseqid == otu, , drop = FALSE]
    for (rank in ranks_to_check) {
      assigned_val <- assigned_otus_df[[rank]][assigned_otus_df$qseqid == otu]
      if (is.na(assigned_val) || assigned_val == "" || assigned_val == "Unclassified") {
        blast_vals <- otu_hits[[rank]]
        blast_defined <- blast_vals[!is.na(blast_vals) & blast_vals != "" & blast_vals != "Unclassified"]
        if (length(blast_defined) > 0) {
          bad_ranks <- c(bad_ranks, otu)
          break
        }
      }
    }
  }
  if (length(bad_ranks) > 0) {
    assigned_otus_df <- assigned_otus_df[!assigned_otus_df$qseqid %in% bad_ranks, , drop = FALSE]
    still_remaining <- union(still_remaining, bad_ranks)
    remaining_otus_df <- blast_filtered[blast_filtered$qseqid %in% still_remaining, , drop = FALSE]
  }
  
  # === FINAL: enforce fully classified at ALL ranks (kingdom->species) ===
  is_unclassified <- function(x) is.na(x) | x == "" | x == "Unclassified"
  incomplete <- Reduce(`|`, lapply(taxonomic_levels, function(r) is_unclassified(assigned_otus_df[[r]])))
  
  if (any(incomplete)) {
    bad2 <- assigned_otus_df$qseqid[incomplete]
    assigned_otus_df <- assigned_otus_df[!assigned_otus_df$qseqid %in% bad2, , drop = FALSE]
    still_remaining <- union(still_remaining, bad2)
    remaining_otus_df <- blast_filtered[blast_filtered$qseqid %in% still_remaining, , drop = FALSE]
  }
  
  # === NEW FINAL: enforce best fungal pident >= 98.5 for ANY easy assignment ===
  if (nrow(assigned_otus_df) > 0) {
    bad_pident <- character()
    for (otu in assigned_otus_df$qseqid) {
      hits <- blast_filtered[blast_filtered$qseqid == otu & tolower(blast_filtered$kingdom) == "fungi", , drop = FALSE]
      if (nrow(hits) == 0) { bad_pident <- c(bad_pident, otu); next }
      pid <- suppressWarnings(as.numeric(hits$pident))
      best_pid <- max(pid, na.rm = TRUE)
      if (!is.finite(best_pid) || best_pid < 98.5) bad_pident <- c(bad_pident, otu)
    }
    if (length(bad_pident) > 0) {
      assigned_otus_df <- assigned_otus_df[!assigned_otus_df$qseqid %in% bad_pident, , drop = FALSE]
      still_remaining <- union(still_remaining, bad_pident)
      remaining_otus_df <- blast_filtered[blast_filtered$qseqid %in% still_remaining, , drop = FALSE]
    }
  }
  
  assigned_otus_df <- assigned_otus_df[, expected_columns, drop = FALSE]
  
  list(
    assigned_otus_df = assigned_otus_df,
    remaining_otus_df = remaining_otus_df
  )
}