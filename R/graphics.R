```r
#' Save taxonomy summary charts and tables to multi-page PDF
#'
#' @param all_results Combined assignments table from write_initial_assignments
#' @param hist_plot ggplot2 object for histogram
#' @param pdf_file Output path for multi-page PDF. If NULL (default), no file is written.
#' @param caption_texts Vector of captions for PDF pages (optional)
#' @param rank_names Vector of rank names (default: c("Phylum",...))
#' @param verbose Logical; if TRUE, emit a message when a PDF is written. Default FALSE.
#' @return List with plots/tables; includes pdf_file when written.
#' @importFrom grDevices pdf dev.off
#' @export
save_taxonomy_graphics <- function(
    all_results,
    hist_plot,
    pdf_file = NULL,
    caption_texts = NULL,
    rank_names = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
    verbose = FALSE
) {
  pkgs <- c("ggplot2", "gridExtra", "reshape2", "RColorBrewer", "grid")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package %s needed for graphics.", pkg))
    }
  }
  
  if (is.null(caption_texts)) {
    caption_texts <- c(
      "",
      "Summary of each major assingment step",
      "",
      ".",
      "Table of unique fungal taxa per rank"
    )
  }
  
  # Step summary table
  step_labels <- c(
    "# of OTUs in Representative Sequence",
    "# of OTUs pass quality check",
    "OTUs assigned to fungal kingdom",
    "OTUs assigned to fungal phylum",
    "OTUs assigned to fungal class",
    "OTUs assigned to fungal order",
    "OTUs assigned to fungal family",
    "OTUs assigned to fungal genus",
    "OTUs assigned to fungal species"
  )
  
  total_repseqs <- nrow(all_results)
  not_failed <- !grepl("failed", all_results$notes, ignore.case = TRUE)
  count_qc <- sum(not_failed)
  
  kingdom_fungi <- not_failed & all_results$kingdom == "Fungi"
  count_kingdom <- sum(kingdom_fungi)
  
  tax_levels <- c("phylum", "class", "order", "family", "genus", "species")
  assigned_counts <- numeric(length = length(tax_levels))
  unassigned_counts <- numeric(length = length(tax_levels))
  for (i in seq_along(tax_levels)) {
    lvl <- tax_levels[i]
    assigned_counts[i] <- sum(
      kingdom_fungi &
        !is.na(all_results[[lvl]]) &
        all_results[[lvl]] != "" &
        all_results[[lvl]] != "Unclassified"
    )
    unassigned_counts[i] <- count_kingdom - assigned_counts[i]
  }
  
  all_counts <- c(
    total_repseqs,
    count_qc,
    count_kingdom,
    assigned_counts
  )
  percents <- round(100 * all_counts / total_repseqs, 1)
  step_table <- data.frame(
    Step = step_labels,
    Count = paste0(all_counts, " (", percents, "%)"),
    stringsAsFactors = FALSE
  )
  step_table_grob <- gridExtra::tableGrob(step_table, rows = NULL)
  
  # Assignment bar chart: Assigned and Unassigned per rank (within kingdom=Fungi)
  df_bar <- data.frame(
    Rank = factor(rank_names, levels = rank_names),
    Assigned = assigned_counts,
    Unassigned = unassigned_counts,
    PercentAssigned = round(100 * assigned_counts / count_kingdom, 1),
    PercentUnassigned = round(100 * unassigned_counts / count_kingdom, 1)
  )
  bar_df <- reshape2::melt(df_bar, id.vars = "Rank", measure.vars = c("Assigned", "Unassigned"))
  bar_df$Percent <- c(df_bar$PercentAssigned, df_bar$PercentUnassigned)
  
  assignment_bar_plot <- ggplot2::ggplot(bar_df, ggplot2::aes(x = Rank, y = Percent, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(Percent, "%\n(", value, ")")),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 4, color = "black"
    ) +
    ggplot2::scale_fill_manual(values = c("Assigned" = "skyblue", "Unassigned" = "gray80")) +
    ggplot2::ylim(0, 100) +
    ggplot2::labs(
      title = NULL,
      x = "Taxonomic Rank",
      y = "Percentage of kingdom=Fungi OTUs",
      fill = ""
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_blank()
    )
  
  # Phylum stacked bar chart
  phylum_vals <- all_results$phylum[kingdom_fungi]
  phylum_vals <- phylum_vals[!is.na(phylum_vals) & phylum_vals != "" & phylum_vals != "Unclassified"]
  phylum_tbl <- table(phylum_vals)
  phylum_count <- as.data.frame(phylum_tbl)
  colnames(phylum_count) <- c("Fungal Phylum", "OTU Count")
  phylum_count$Sample <- "All OTUs"
  phylum_plot <- ggplot2::ggplot(phylum_count, ggplot2::aes(x = Sample, y = `OTU Count`, fill = `Fungal Phylum`)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_brewer(palette = "Set3") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = "Fungal Phyla",
      y = "OTU Count",
      fill = "Fungal Phylum"
    ) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  
  # Unique taxa counts table (only for kingdom = "Fungi")
  kingdom_fungi <- !grepl("failed", all_results$notes, ignore.case = TRUE) & all_results$kingdom == "Fungi"
  
  unique_tax_counts <- sapply(tax_levels, function(lvl) {
    vals <- all_results[[lvl]][kingdom_fungi]
    vals <- vals[!is.na(vals) & vals != "" & vals != "Unclassified"]
    length(unique(vals))
  })
  final_count_tbl <- data.frame(
    "Rank" = rank_names,
    "Unique Count" = unique_tax_counts,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  final_count_grob <- gridExtra::tableGrob(final_count_tbl, rows = NULL)
  
  # Caption helper
  add_caption <- function(caption, y = 0.025, gp = grid::gpar(fontsize = 10, col = "gray30")) {
    if (!is.null(caption) && !is.na(caption) && nchar(caption) > 0) {
      grid::grid.text(label = caption, x = 0.5, y = y, just = "center", gp = gp)
    }
  }
  
  # PDF export sequence (only if requested)
  if (!is.null(pdf_file)) {
    dir_out <- dirname(pdf_file)
    if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    
    grDevices::pdf(pdf_file, width = 7, height = 5)
    on.exit(grDevices::dev.off(), add = TRUE)
    
    page_counter <- 1
    
    # 1. Alignment length histogram (FIRST PAGE)
    grid::grid.newpage()
    grid::grid.draw(ggplot2::ggplotGrob(hist_plot))
    if (length(caption_texts) >= page_counter) add_caption(caption_texts[page_counter])
    page_counter <- page_counter + 1
    
    # 2. Step summary table
    grid::grid.newpage()
    grid::grid.draw(step_table_grob)
    if (length(caption_texts) >= page_counter) add_caption(caption_texts[page_counter])
    page_counter <- page_counter + 1
    
    # 3. Assignment bar chart
    grid::grid.newpage()
    grid::grid.draw(ggplot2::ggplotGrob(assignment_bar_plot))
    if (length(caption_texts) >= page_counter) add_caption(caption_texts[page_counter])
    page_counter <- page_counter + 1
    
    # 4. Phylum stacked bar chart
    grid::grid.newpage()
    grid::grid.draw(ggplot2::ggplotGrob(phylum_plot))
    if (length(caption_texts) >= page_counter) add_caption(caption_texts[page_counter])
    page_counter <- page_counter + 1
    
    # 5. Unique taxa count table
    grid::grid.newpage()
    grid::grid.draw(final_count_grob)
    if (length(caption_texts) >= page_counter) add_caption(caption_texts[page_counter])
    
    if (isTRUE(verbose)) {
      message("Saved taxonomy graphics and summary tables to: ", pdf_file)
    }
  }
  
  invisible(list(
    pdf_file = pdf_file,
    hist_plot = hist_plot,
    step_table = step_table,
    assignment_bar_plot = assignment_bar_plot,
    phylum_plot = phylum_plot,
    final_count_tbl = final_count_tbl
  ))
}
```