easy_assignments <- function(
    blast_filtered,
    cutoffs_file = NULL,
    default_cutoff = 98,
    package_name = NULL
) {
  taxonomic_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  expected_columns <- c("qseqid", taxonomic_levels)
  assigned_otus <- list()
  easy_otus <- c()
  all_otus <- unique(blast_filtered$qseqid)
  
  # --- Path resolution for cutoffs_file ---
  if (is.null(cutoffs_file)) {
    # Try to find default file via system.file (works if package is installed)
    cutoffs_file <- system.file("extdata", "taxonomy_cutoffs.csv", package = "ClassifyITS")
    # Fallback for development mode (not installed, just source code)
    if (!file.exists(cutoffs_file) || nchar(cutoffs_file) == 0) {
      cutoffs_file <- "inst/extdata/taxonomy_cutoffs.csv"
    }
  }
  if (!file.exists(cutoffs_file) || nchar(cutoffs_file) == 0) {
    stop(sprintf(
      "Cutoffs file '%s' not found. Please provide a valid path or ensure the package includes inst/extdata/taxonomy_cutoffs.csv.",
      cutoffs_file
    ))
  }
  
  cutoffs_data <- parse_taxonomy_cutoffs(cutoffs_file)$long
  phylum_species_cutoff <- cutoffs_data[cutoffs_data$rank == "species" &
                                          cutoffs_data$cutoff_type == "percent_identity", ]
  phylum_species_cutoff_vec <- setNames(as.numeric(phylum_species_cutoff$cutoff_value), phylum_species_cutoff$phylum)
  
  # --- STEP 1: Assign easy OTUs (100% identity, kingdom Fungi, matching rules)
  for (otu in all_otus) {
    hits <- blast_filtered[blast_filtered$qseqid == otu & as.numeric(blast_filtered$pident) == 100, ]
    if (nrow(hits) == 0) next
    fungi_hits <- hits[tolower(hits$kingdom) == "fungi", ]
    if (nrow(fungi_hits) == 0) next
    
    if (nrow(fungi_hits) == 1) {
      taxonomy <- as.list(fungi_hits[1, taxonomic_levels])
      taxonomy$qseqid <- otu
      assigned_otus[[otu]] <- taxonomy
      easy_otus <- c(easy_otus, otu)
      next
    }
    other_levels <- taxonomic_levels[1:6]
    ref_tax <- as.character(fungi_hits[1, other_levels])
    all_other_same <- all(apply(fungi_hits[, other_levels, drop=FALSE], 1, function(vals) all(as.character(vals) == ref_tax)))
    unique_species <- unique(fungi_hits$species)
    unique_species <- unique_species[!is.na(unique_species) & unique_species != "Unclassified" & unique_species != ""]
    if (all_other_same) {
      if (length(unique_species) == 1) {
        taxonomy <- as.list(fungi_hits[1, taxonomic_levels])
        taxonomy$qseqid <- otu
        assigned_otus[[otu]] <- taxonomy
        easy_otus <- c(easy_otus, otu)
      } else if (length(unique_species) > 1) {
        taxonomy <- as.list(fungi_hits[1, other_levels])
        taxonomy$species <- "Unclassified"
        taxonomy$qseqid <- otu
        assigned_otus[[otu]] <- taxonomy
        easy_otus <- c(easy_otus, otu)
      }
    }
    # If other ranks differ, skip assignment for now
  }
  
  remaining_otus <- setdiff(all_otus, easy_otus)
  
  # --- STEP 2: Assign using best evalue & phylum-specific percent cutoff
  for (otu in remaining_otus) {
    otu_hits <- blast_filtered[blast_filtered$qseqid == otu, ]
    if (nrow(otu_hits) == 0) next
    
    best_hit_idx <- which.min(as.numeric(otu_hits$evalue))
    best_hit <- otu_hits[best_hit_idx, ]
    if (tolower(best_hit$kingdom) != "fungi") next
    
    phylum_val <- best_hit$phylum
    cutoff <- if (!is.na(phylum_species_cutoff_vec[phylum_val])) phylum_species_cutoff_vec[phylum_val] else default_cutoff
    
    high_hits <- otu_hits[as.numeric(otu_hits$pident) >= cutoff, ]
    if (nrow(high_hits) == 0) next
    
    genus_check <- function(hit, others) {
      if (nrow(others) == 0) return(TRUE)
      all_genus_disagree <- all(others$genus != hit$genus | is.na(others$genus) | others$genus == "Unclassified" | others$genus == "")
      return(!all_genus_disagree)
    }
    assigned_taxonomy <- NULL
    if (nrow(high_hits) == 1) {
      if (genus_check(high_hits[1,], high_hits[-1,])) {
        assigned_taxonomy <- as.list(high_hits[1, taxonomic_levels])
        assigned_taxonomy$qseqid <- otu
      }
    }
    if (is.null(assigned_taxonomy) && nrow(high_hits) > 1) {
      pidents <- as.numeric(high_hits$pident)
      best_idx <- which.max(pidents)
      pdiff <- pidents[best_idx] - pidents[-best_idx]
      if (any(pdiff >= 1)) {
        if (genus_check(high_hits[best_idx, ], high_hits[-best_idx, ])) {
          assigned_taxonomy <- as.list(high_hits[best_idx, taxonomic_levels])
          assigned_taxonomy$qseqid <- otu
        }
      }
    }
    if (is.null(assigned_taxonomy) && nrow(high_hits) > 1) {
      unique_genus <- unique(high_hits$genus)
      unique_genus <- unique_genus[!is.na(unique_genus) & unique_genus != "Unclassified" & unique_genus != ""]
      if (length(unique_genus) == 1) {
        assigned_taxonomy <- as.list(high_hits[1, taxonomic_levels[1:6]])
        assigned_taxonomy$species <- "Unclassified"
        assigned_taxonomy$qseqid <- otu
      }
    }
    if (!is.null(assigned_taxonomy)) {
      for (lvl in taxonomic_levels[1:6]) {
        if (assigned_taxonomy[[lvl]] == "Unclassified" || is.na(assigned_taxonomy[[lvl]]) || assigned_taxonomy[[lvl]] == "") {
          check_hits <- otu_hits[as.numeric(otu_hits$pident) >= 90, ]
          defined_vals <- check_hits[[lvl]][!is.na(check_hits[[lvl]]) & check_hits[[lvl]] != "" & check_hits[[lvl]] != "Unclassified"]
          if (length(defined_vals) > 0 && length(unique(defined_vals)) == 1) {
            assigned_taxonomy[[lvl]] <- unique(defined_vals)
          }
        }
      }
      assigned_otus[[otu]] <- assigned_taxonomy
      next
    }
    # If genus does not agree, skip for now
  }
  
  # -- Build assigned_otus_df, guaranteeing all columns exist and are ordered --
  assigned_otus_df <- do.call(rbind, lapply(assigned_otus, function(x) {
    vals <- setNames(rep(NA, length(expected_columns)), expected_columns)
    for (name in names(x)) {
      if (name %in% expected_columns) {
        vals[[name]] <- x[[name]]
      }
    }
    as.data.frame(as.list(vals), stringsAsFactors = FALSE)
  }))
  rownames(assigned_otus_df) <- NULL  # remove rownames
  
  still_remaining <- setdiff(all_otus, assigned_otus_df$qseqid)
  remaining_otus_df <- blast_filtered[blast_filtered$qseqid %in% still_remaining, ]
  
  # === FINAL KINGDOM CHECK ===
  bad_otus <- c()
  for (otu in assigned_otus_df$qseqid) {
    otu_hits <- blast_filtered[blast_filtered$qseqid == otu, ]
    n_hits <- nrow(otu_hits)
    n_fungal <- sum(tolower(otu_hits$kingdom) == "fungi")
    if (n_hits == 0) next
    if (n_fungal / n_hits <= 0.5) {
      bad_otus <- c(bad_otus, otu)
    }
  }
  if (length(bad_otus) > 0) {
    assigned_otus_df <- assigned_otus_df[!assigned_otus_df$qseqid %in% bad_otus, ]
    still_remaining <- union(still_remaining, bad_otus)
    remaining_otus_df <- blast_filtered[blast_filtered$qseqid %in% still_remaining, ]
  }
  
  # === FINAL PHYLUM, CLASS, ORDER DOUBLE CHECK ===
  bad_ranks <- c()
  ranks_to_check <- c("phylum", "class", "order")
  for (otu in assigned_otus_df$qseqid) {
    otu_hits <- blast_filtered[blast_filtered$qseqid == otu, ]
    for (rank in ranks_to_check) {
      assigned_val <- assigned_otus_df[assigned_otus_df$qseqid == otu, rank]
      # Is assigned value unclassified/blank/NA?
      if (is.na(assigned_val) || assigned_val == "" || assigned_val == "Unclassified") {
        # Does any BLAST hit have a defined value at this rank?
        blast_vals <- otu_hits[[rank]]
        blast_defined <- blast_vals[!is.na(blast_vals) & blast_vals != "" & blast_vals != "Unclassified"]
        if (length(blast_defined) > 0) {
          bad_ranks <- c(bad_ranks, otu)
          break # move on to next OTU after marking this one
        }
        # If all BLAST hits at this rank are unclassified/blank/NA, keep assignment
      }
    }
  }
  if (length(bad_ranks) > 0) {
    assigned_otus_df <- assigned_otus_df[!assigned_otus_df$qseqid %in% bad_ranks, ]
    still_remaining <- union(still_remaining, bad_ranks)
    remaining_otus_df <- blast_filtered[blast_filtered$qseqid %in% still_remaining, ]
  }
  
  # Final: guarantee column order
  assigned_otus_df <- assigned_otus_df[, expected_columns, drop=FALSE]
  
  return(list(
    assigned_otus_df = assigned_otus_df,
    remaining_otus_df = remaining_otus_df
  ))
}