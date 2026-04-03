## Quick Start

Once you have installed the package, and prepared your BLAST results as described in Data Preparation, you can run the assignment pipeline with just a few lines:

    ITS_taxonomy <- ITS_assignment(
        blast_file = "input/blast_results.tsv",            # Path to BLAST .TSV file
        rep_fasta = "input/representative_seqs.fasta"      # Path to FASTA file containg representative sequences used to generate the BLAST results
    )

By default, ITS_assignment() does not write any files. If you’d like the pipeline to also export the assignments table (CSV) and summary graphics (PDF), supply an output directory:

    ITS_taxonomy <- ITS_assignment(
      blast_file = "input/blast_results.tsv",
      rep_fasta  = "input/representative_seqs.fasta",
      outdir     = "outputs"
    )


While not necessary for a successful run, the assignment process can easily be customized with additional parameters:

    ITS_taxonomy <- ITS_assignment(
      blast_file       = "input/blast_results.tsv",
      rep_fasta        = "input/representative_seqs.fasta",
      cutoffs_file     = NULL,       # Optional: custom taxonomy cutoffs CSV file
      cutoff_fraction  = 0.6,        # Optional: fraction for alignment length QC
      n_cutoff         = 1,          # Optional: percent N cutoff
      outdir           = "outputs",  # Optional: output directory (writes CSV/PDF when provided)
     verbose          = FALSE       # Optional: print progress messages
    )

- `cutoffs_file`: Path to a custom taxonomy cutoffs CSV file (see extdata in package).
- `cutoff_fraction`: Fraction of median alignment length for BLAST QC filter.
- `n_cutoff`: Maximum percent N bases (ambiguous bases) allowed in any sequence.
- `outdir`: Output directory for results. If NULL, nothing is written; if provided, the pipeline writes initial_assignments.csv and combined_taxonomy_graphics.pdf.
- `verbose`: If TRUE, emits progress messages.


The function returns a list including the main taxonomy assignment table, QC information, and summary figures. See documentation for advanced workflow options.