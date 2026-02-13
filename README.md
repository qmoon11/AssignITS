# ClassifyITS

**ClassifyITS** is an R package for quality control, taxonomy assignment, and visualization of fungal OTU representative sequences based on BLAST results.

---

## Features

- Performs QC on user-imported representative sequence FASTA and BLAST results
- Assigns taxonomy using kingdom through species taxon-specific e-values and sequence similarity thresholds
- Outputs an initial CSV taxonomy table for manual inspection, as well as summary graphics and statistics

---
## Installation

You can download and install ClassifyITS directly from GitHub using the `devtools` package in R:

```r
# Install devtools if you don't have it yet
install.packages("devtools")

# Install ClassifyITS from GitHub
devtools::install_github("qmoon11/ClassifyITS")

# Load the package
library(ClassifyITS)
```
---

## Prerequisites

This package depends on [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), a Bioconductor package for manipulating biological sequences in R.

**Important:**  
Biostrings is _not_ installed automatically via CRAN or `install.packages()`.  
You must install BiocManager and then Biostrings before installing or using ClassifyITS.

```r
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Biostrings
BiocManager::install("Biostrings")

library(Biostrings)
```


## Documentation & User Guides

See [Data Preparation](data-preparation.md) for detailed BLAST setup and input requirements, including:
- How to use NCBI BLAST+ with the UNITE eukaryote database
- Proper FASTA and BLAST result formatting
- Example shell commands for HPC or local runs

See [Quick Start](quick-start.md) for a step-by-step usage guide, including:
- Example R code for running the pipeline with `ITS_assignment`
- Minimal installation and first run walkthrough

See [Custom Cutoffs](custom-cutoffs.md) to fine-tune assignments, including:
- How to provide your own taxonomy cutoffs file
- Example cutoff file formats
- Using the `cutoffs_file`, `cutoff_fraction`, `n_cutoff`, and `outdir` arguments

See [Troubleshooting](troubleshooting.md) for help with common errors, including:
- Installation tips for Biostrings and dependencies
- BLAST output format issues
- FAQ and support contact information

See [Check Contamination](check-contamination.md) for descriptions of how to check for common contaminants in your dataset, including:
- Common human associated contaminants and how to identify them in your BLAST results

See [Citations](citations.md) for a list of relevant literature
- please cite ClassifyITS in your publications if you use it for your research!
