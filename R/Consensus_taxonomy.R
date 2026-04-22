#' Per-rank consensus filter for taxonomy assignment
#'
#' Only confirms or demotes, never promotes Unclassified.
#' FINAL hierarchy check: if any rank is Unclassified, all lower ranks are forced to Unclassified.
#'
#' @param final_table Data frame of taxonomic assignments.
#' @param blast_qc Data frame of filtered BLAST hits for each OTU.
#' @return Data frame of consensus assignments (same structure as input).
#' @export
consensus_taxonomy_assignment <- function(final_table, blast_qc) {
  tax_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  final_table_consensus <- final_table
  
  norm_val <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | x == "" | tolower(x) == "unclassified"] <- "Unclassified"
    x
  }
  
  for (i in seq_len(nrow(final_table))) {
    otu_id <- final_table$qseqid[i]
    otu_hits <- blast_qc[blast_qc$qseqid == otu_id, , drop = FALSE]
    if (nrow(otu_hits) == 0) next
    
    # ensure numeric e-value ordering for "top 2" rule
    otu_hits$evalue_num <- suppressWarnings(as.numeric(otu_hits$evalue))
    otu_hits <- otu_hits[order(otu_hits$evalue_num, na.last = TRUE), , drop = FALSE]
    
    # ---- Gate: compare BEST vs SECOND; only if disagreement do consensus ----
    run_consensus <- FALSE
    if (nrow(otu_hits) >= 2) {
      best_hit   <- otu_hits[1, , drop = FALSE]
      second_hit <- otu_hits[2, , drop = FALSE]
      
      assigned_vals <- norm_val(final_table[i, tax_ranks, drop = TRUE])
      classified_ranks <- tax_ranks[assigned_vals != "Unclassified"]
      
      if (length(classified_ranks) > 0) {
        for (rank in classified_ranks) {
          if (norm_val(best_hit[[rank]]) != norm_val(second_hit[[rank]])) {
            run_consensus <- TRUE
            break
          }
        }
      }
    }
    
    # if no disagreement (or no 2nd hit), keep original row unchanged
    if (!run_consensus) next
    
    # ---- Per-rank consensus demotion (top2 agreement OR majority vote) ----
    for (rank in tax_ranks) {
      orig_val <- norm_val(final_table[[rank]][i])
      
      if (orig_val == "Unclassified") {
        final_table_consensus[[rank]][i] <- "Unclassified"
        next
      }
      
      ranks_vector <- norm_val(otu_hits[[rank]])
      ranks_vector <- ranks_vector[ranks_vector != "Unclassified"]
      
      if (length(ranks_vector) == 0) {
        final_table_consensus[[rank]][i] <- "Unclassified"
        next
      }
      
      # top2 agree at this rank AND match original
      if (nrow(otu_hits) >= 2) {
        top2_vals <- norm_val(otu_hits[1:2, rank, drop = TRUE])
        if (length(top2_vals) == 2 && length(unique(top2_vals)) == 1 && top2_vals[1] == orig_val) {
          final_table_consensus[[rank]][i] <- orig_val
          next
        }
      }
      
      # majority (>50%) among defined hits matches original
      val_tab <- table(ranks_vector)
      max_ratio <- max(val_tab) / sum(val_tab)
      majority_val <- names(val_tab)[which.max(val_tab)]
      
      if (max_ratio > 0.5 && majority_val == orig_val) {
        final_table_consensus[[rank]][i] <- orig_val
      } else {
        final_table_consensus[[rank]][i] <- "Unclassified"
      }
    }
  }
  
  # === FINAL hierarchy cleanup: no classification allowed below an Unclassified above ===
  for (i in seq_len(nrow(final_table_consensus))) {
    seen_uncl <- FALSE
    for (rank in tax_ranks) {
      v <- norm_val(final_table_consensus[[rank]][i])
      if (seen_uncl || v == "Unclassified") {
        seen_uncl <- TRUE
        final_table_consensus[[rank]][i] <- "Unclassified"
      } else {
        final_table_consensus[[rank]][i] <- v
      }
    }
  }
  
  final_table_consensus
}


#' Safely rbinds list of data frames, ensuring columns match
#'
#' @param dfs List of data frames
#' @param all_cols Vector of required columns
#' @return Combined data frame
#' @export
safe_rbind_list <- function(dfs, all_cols = NULL) {
  dfs <- Filter(function(x) !is.null(x) && nrow(x) > 0, dfs)
  if (length(dfs) == 0) {
    if (!is.null(all_cols)) {
      empty <- as.data.frame(
        matrix(NA_character_, nrow = 0, ncol = length(all_cols)),
        stringsAsFactors = FALSE
      )
      colnames(empty) <- all_cols
      return(empty)
    } else {
      return(data.frame())
    }
  }
  if (!is.null(all_cols)) {
    dfs <- lapply(dfs, function(x) ensure_cols(x, all_cols))
  }
  do.call(rbind, dfs)
}

#' Ensure data frame has all required columns (as character)
#'
#' @param df Data frame to fix
#' @param all_cols Vector of required columns
#' @return Fixed data frame (in correct order, with all columns present)
#' @export
ensure_cols <- function(df, all_cols) {
  for (col in all_cols) {
    if (!col %in% colnames(df)) df[[col]] <- NA_character_
    df[[col]] <- as.character(df[[col]])
  }
  df <- df[, all_cols, drop = FALSE]
  rownames(df) <- NULL
  return(df)
}

#' Create and write the initial assignments table including drops at all steps
#'
#' @param easy_df Data frame of easy-assigned OTUs
#' @param consensus_df Data frame of consensus-assigned OTUs (hard ones)
#' @param rep_seqs DNAStringSet or named character vector of rep seqs
#' @param blast Data frame of all BLAST results
#' @param blast_filtered Data frame of filtered BLAST results
#' @param file Path for output CSV. If NULL (default), no file is written.
#' @param verbose Logical; if TRUE emit a message when a file is written. Default FALSE.
#' @return Data frame of assignments (written if file is not NULL)
#' @export
write_initial_assignments <- function(
    easy_df, consensus_df, rep_seqs, blast, blast_filtered,
    file = NULL,
    verbose = FALSE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  tax_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  all_cols <- c("qseqid", tax_cols, "notes")
  
  # 1. Successful assignments: add note
  easy_df <- ensure_cols(easy_df, all_cols)
  easy_df$notes <- "ITS pipeline completed"
  consensus_df <- ensure_cols(consensus_df, all_cols)
  consensus_df$notes <- "ITS pipeline completed"
  assigned <- rbind(easy_df, consensus_df)
  assigned_otus <- unique(assigned$qseqid)
  
  # 2. Failed BLAST
  rep_otus <- names(rep_seqs)
  blast_otus <- unique(blast$qseqid)
  blast_filtered_otus <- unique(blast_filtered$qseqid)
  otus_no_blast <- setdiff(rep_otus, blast_otus)
  df_no_blast <- if (length(otus_no_blast) > 0) {
    data.frame(
      qseqid = otus_no_blast,
      kingdom = NA_character_,
      phylum  = NA_character_,
      class   = NA_character_,
      order   = NA_character_,
      family  = NA_character_,
      genus   = NA_character_,
      species = NA_character_,
      notes   = "Failed QC: no BLAST result passed Quality Control step (no BLAST result, too many N, alignment too short, low quality results)",
      stringsAsFactors = FALSE
    )
  } else NULL
  
  # 3. Failed QC
  otus_failed_qc <- setdiff(blast_otus, blast_filtered_otus)
  otus_failed_qc <- setdiff(otus_failed_qc, assigned_otus)
  otus_failed_qc <- setdiff(otus_failed_qc, otus_no_blast)
  df_failed_qc <- if (length(otus_failed_qc) > 0) {
    data.frame(
      qseqid = otus_failed_qc,
      kingdom = NA_character_,
      phylum  = NA_character_,
      class   = NA_character_,
      order   = NA_character_,
      family  = NA_character_,
      genus   = NA_character_,
      species = NA_character_,
      notes   = "Failed QC: no BLAST result passed Quality Control step (no BLAST result, too many N, alignment too short, low quality results)",
      stringsAsFactors = FALSE
    )
  } else NULL
  
  # 4. Combine all results and clean up
  all_results <- rbind(
    assigned[, all_cols, drop = FALSE],
    if (!is.null(df_no_blast)) df_no_blast[, all_cols, drop = FALSE] else NULL,
    if (!is.null(df_failed_qc)) df_failed_qc[, all_cols, drop = FALSE] else NULL
  )
  all_results <- all_results[!duplicated(all_results$qseqid), ]
  all_results <- all_results[all_results$qseqid != "qseqid", , drop = FALSE]
  all_results <- all_results[order(all_results$qseqid), ]
  rownames(all_results) <- NULL
  
  # Guarantee NA for failed rows in taxonomy columns
  is_fail <- grepl("Failed", all_results$notes)
  for (col in tax_cols) {
    all_results[[col]][is_fail] <- NA_character_
  }
  
  # Write output only if requested
  if (!is.null(file)) {
    dir_out <- dirname(file)
    if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
    data.table::fwrite(all_results, file, sep = ",", na = "NA")
    if (isTRUE(verbose)) message("Taxonomy assignments exported to ", file)
  }
  
  all_results
}