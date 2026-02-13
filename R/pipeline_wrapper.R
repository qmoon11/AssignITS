#' Complete Fungal Assignment Pipeline
#'
#' Runs all steps: QC, filtering, plotting, assignments, writes outputs.
#' @param blast_file Path to BLAST results TSV file
#' @param rep_fasta Path to representative sequences FASTA file
#' @param cutoffs_file Path to taxonomy cutoffs CSV file (optional; defaults to package example if omitted)
#' @param cutoff_fraction Numeric, fraction of median rep-seq length for BLAST filtering (default: 0.6)
#' @param n_cutoff Numeric, N base percentage cutoff (default: 1)
#' @param outdir Output directory for results (default: "outputs")
#' @return Named list of results and output file paths
#' @export
ITS_assignment <- function(
    blast_file,
    rep_fasta,
    cutoffs_file = NULL,
    cutoff_fraction = 0.6,
    n_cutoff = 1,
    outdir = "outputs"
) {
  # --- Robust path resolution for cutoffs_file ---
  if (is.null(cutoffs_file) || !file.exists(cutoffs_file)) {
    cutoffs_file <- system.file("extdata", "taxonomy_cutoffs.csv", package = "ClassifyITS")
    # Fallback for dev mode: try local package path if system.file fails (returns "")
    if (!file.exists(cutoffs_file) || nchar(cutoffs_file) == 0) {
      cutoffs_file <- file.path("inst", "extdata", "taxonomy_cutoffs.csv")
    }
  }
  if (!file.exists(cutoffs_file) || nchar(cutoffs_file) == 0) {
    stop("Could not locate taxonomy_cutoffs.csv. Please reinstall the package or supply a custom file.")
  }
  
  # Set up output paths
  hist_dir <- file.path(outdir, "individual_plots")
  if (!dir.exists(hist_dir)) dir.create(hist_dir, recursive = TRUE)
  assignments_file <- file.path(outdir, "initial_assignments.csv")
  pdf_file <- file.path(hist_dir, "combined_taxonomy_graphics.pdf")
  
  # Step 1: QC load
  ldat <- load_and_check(blast_file, rep_fasta)
  rep_seqs <- ldat$rep_seqs
  blast_all <- ldat$blast
  
  # Step 2: Trim BLAST for alignment length
  blast_filtered <- trim_alignments(blast_all, rep_seqs, fraction = cutoff_fraction)
  
  # Step 3: Check N-bases
  N_check <- check_N(rep_seqs, cutoff = n_cutoff)
  n_fail_otus <- N_check$qseqid[N_check$N_flag]
  
  # Step 4: Warn if any OTUs failed QC
  total_otus <- length(names(rep_seqs))
  qc_pass <- blast_all$qseqid
  failed_otus <- setdiff(names(rep_seqs), qc_pass)
  num_failed <- length(failed_otus)
  if (num_failed > 0) {
    warning(sprintf(
      "Warning: %d of your %d FASTA sequences failed QC and could not be classified using this pipeline due to missing or poor BLAST results.",
      num_failed, total_otus
    ))
  }
  
  # Step 5: Create alignment histogram plot object (for graphics)
  hist_plot <- plot_alignment_hist(blast_all, rep_seqs, cutoff_fraction = cutoff_fraction)
  
  # Step 6: Easy assignments
  easy_results <- easy_assignments(blast_filtered, cutoffs_file = cutoffs_file)
  easy_assignments_df <- easy_results$assigned_otus_df
  remaining_otus_df <- easy_results$remaining_otus_df
  
  # Step 7: Consensus assignment for remaining
  consensus_assignments_df <- consensus_taxonomy_assignment(remaining_otus_df, blast_filtered)
  
  # Step 8: Final output CSV table: Save and return results object
  all_results <- write_initial_assignments(
    easy_df = easy_assignments_df,
    consensus_df = consensus_assignments_df,
    rep_seqs = rep_seqs,
    blast = blast_all,
    blast_filtered = blast_filtered,
    file = assignments_file
  )
  
  # Step 9 - Save taxonomy summary graphics as PDF using in-memory plot objects
  save_taxonomy_graphics(
    all_results = all_results,
    hist_plot = hist_plot,
    pdf_file = pdf_file
  )
  
  # Return top-level results and file paths
  invisible(list(
    rep_seqs            = rep_seqs,
    blast_all           = blast_all,
    blast_filtered      = blast_filtered,
    N_check             = N_check,
    n_fail_otus         = n_fail_otus,
    easy_assignments_df = easy_assignments_df,
    consensus_assignments_df = consensus_assignments_df,
    all_results         = all_results,
    assignments_file    = assignments_file,
    pdf_file            = pdf_file
  ))
}