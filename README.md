# ClassifyITS

**ClassifyITS** is an R package for quality control, taxonomy assignment, and visualization of fungal OTU representative sequences based on user provided BLAST results.

Fungi are ubiquitous in Earth's wonderfully diverse ecosystems. The AssignITS package aids in the taxonomic classification of environmental ITS short-read barcoding data.
Unlike previous methods, it employs taxon-specific e-value and percent identity cutoffs at each taxonomic rank from kingdom to species. 
The package takes a conservative approach and outputs both graphics and user-friendly files to help users manually inspect fungal OTUs that fail classification at relevant levels (e.g., Phylum). 
AssignITS is based on taxonomic cutoff criteria from "The Global Soil Mycobiome consortium dataset for boosting fungal diversity research" (Fungal Diversity, Tedersoo et al., 2021, doi:10.1007/s13225-021-00493-7) and "Best practices in metabarcoding of fungi: From experimental design to results" (Molecular Ecology, Tedersoo et al., 2022, doi:10.1111/mec.16460).

---

## Features

- Performs QC on user-imported representative sequence FASTA and BLAST results
- Assigns taxonomy using kingdom through species taxon-specific e-values and sequence similarity thresholds
- Outputs an initial CSV taxonomy table for manual inspection, as well as summary graphics and statistics
- [View vignette](vignettes/ClassifyITS.html) for detailed usage instructions and examples
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

## Inspecting taxonomic assingments

It is recommended to carefully examine 1) OTUs that failed the pipeline 2) any OTUs that failed to assign at kingdom level (see list in outputs/intial.assingments.csv) after classifying) as these are likely sequencing errors.
See [Inspection](docs/inspection.md) for information on inspecting results, including:
- How to interpret the `initial_assignments.csv` output file
- Common reasons for failed assignments and how to identify them
- Tips for manual curation of taxonomy assignments (e.g. manually checking BLAST results for fungal OTUs unassigned at phylum or class level)

## Documentation & User Guides

See [Data Preparation](docs/data-preparation.md) for detailed BLAST setup and input requirements, including:
- How to use NCBI BLAST+ with the UNITE eukaryote database
- Proper FASTA and BLAST result formatting
- Example shell commands for HPC or local runs

See [Quick Start](docs/quick-start.md) for a step-by-step usage guide, including:
- Example R code for running the pipeline with `ITS_assignment`
- Minimal installation and first run walkthrough

See [Custom Cutoffs](docs/custom-cutoffs.md) to fine-tune assignments, including:
- How to provide your own taxonomy cutoffs file
- Example cutoff file formats
- Using the `cutoffs_file`, `cutoff_fraction`, `n_cutoff`, and `outdir` arguments

See [Troubleshooting](docs/troubleshooting.md) for help with common errors, including:
- Installation tips for Biostrings and dependencies
- BLAST output format issues
- FAQ and support contact information

See [Check Contamination](docs/check-contamination.md) for descriptions of how to check for common contaminants in your dataset, including:
- Common human associated contaminants and how to identify them in your BLAST results

See [Citations](docs/citations.md) for a list of relevant literature
- please cite ClassifyITS in your publications if you use it for your research!

---

## About Taxonomic Cutoffs

ClassifyITS uses taxon-specific cutoffs for kingdom through species. At minimum, default cutoffs are provided for all fungal phyla, common classes and orders found in ITS data, and a subset of well studied families, genera, and species. The current list can be seen in `inst/extdata/taxonomy_cutoffs.csv` in the package source code.

You can provide your own custom cutoffs file to fine-tune assignments for your dataset or research question.

If you are knowledgeable in a specific taxonomic group and think ClassifyITS could benefit from modified cutoffs for a particular group, please reach out with proposed cutoffs and supporting citations. Community contributions are welcomed!

> The package is designed so cutoffs can be easily updated as new fungal groups are discovered and taxonomy is refined.

---

## Cutoff Logic

- **E-values** are favored at higher taxonomic levels (kingdom, phylum) due to variability in sequence similarity and database effects.
- **Percent identity** is favored at lower taxonomic levels (genus, species), where research has established thresholds for ITS similarity within species and genera.

Metabarcoding pipelines are inherently reliant on the idea that percent identity can differentiate species or OTUs. ClassifyITS provides a conservative starting point for taxonomy assignment, with flexibility for manual curation and extension.

---

## Future Directions and Community Ideas

- Automatic OTU binning with taxonomic-informed cutoffs remains an exciting area for environemental microbiology.
- For example: updating SWARM software to integrate BLAST+taxonomic cutoffs for more accurate OTU clustering could be valuable. I hope someone who is very good at computer can do this!

---