setwd("/Users/quinnmoon/Downloads/age/ClassifyITS/")

# Source R scripts
source("R/QC.R")
source("R/Easy_Assingments.R")
source("R/Best_Assingment.R")
source("R/Concensus_taxonomy.R")
source("R/Graphics.R")   # Make sure your updated save_taxonomy_graphics is in this file

# Define file paths
blast_file   <- "../megablast_ITS2.tsv"
rep_fasta    <- "../dna-sequences_ITS.fasta"
cutoffs_file <- "/Users/quinnmoon/Downloads/age/ClassifyITS/inst/extdata/taxonomy_cutoffs.csv"

# Run QC and assignments
ldat <- load_and_check(blast_file, rep_fasta)
rep_seqs <- ldat$rep_seqs
blast_all <- ldat$blast

blast_filtered <- trim_alignments(blast_all, rep_seqs, fraction = 0.6)
N_check <- check_N(rep_seqs, cutoff = 1)
n_fail_otus <- N_check$qseqid[N_check$N_flag]

easy_results <- easy_assignments(blast_filtered, cutoffs_file = cutoffs_file)
easy_assignments <- easy_results$assigned_otus_df
remaining_otus_df <- easy_results$remaining_otus_df

consensus_assignments <- consensus_taxonomy_assignment(remaining_otus_df, blast_filtered)

# Now, capture the combined results for downstream graphics/reporting
all_results <- write_initial_assignments(
  easy_df = easy_assignments,
  consensus_df = consensus_assignments,
  rep_seqs = rep_seqs,
  blast = blast_all,      # Unfiltered BLAST dataframe, consistent with earlier steps
  blast_filtered = blast_filtered,
  file = "outputs/initial_assignments.csv"
)

# Generate summary graphics/tables from all_results
hist_plot <- plot_alignment_hist(blast_all, rep_seqs, cutoff_fraction = 0.6)
save_taxonomy_graphics(all_results, hist_plot)



ITS_assignment(
  blast_file = "path/to/blast.tsv",
  rep_fasta = "path/to/sequences.fasta",
  cutoffs_file = "path/to/taxonomy_cutoffs.csv",   # optional
  cutoff_fraction = 0.6,
  n_cutoff = 1,
  outdir = "outputs"
)


install.packages("devtools")
library(devtools)
devtools::document()
devtools::install()
list.files()

library(ClassifyITS)
?check_N 



results <- ITS_assignment(
  blast_file = "age/megablast_ITS2.tsv",
  rep_fasta = "age/dna-sequences_ITS.fasta"
)

getwd()
