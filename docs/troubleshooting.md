# Troubleshooting

This document addresses common problems users may encounter with the ClassifyITS package.

---

## 1. Biostrings Installation and Errors

**Error:**
Error in library(Biostrings): there is no package called 'Biostrings'

**Solution:**
Biostrings is not installed automatically via CRAN. Before using ClassifyITS, install Biostrings with BiocManager:

    # Install BiocManager if needed
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("Biostrings")

Be sure to restart the R session after installation.

---

## 2. BLAST Output Format Errors

**Error:**
Error: BLAST file is missing required columns  
Error reading BLAST results: column 'stitle' not found

**Solution:**
ClassifyITS requires tab-delimited BLAST output with columns exactly as below:

    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

Make sure to use the following BLAST command:

    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

---

## 3. File Not Found Errors

**Error:**
Error: FASTA file not found.  
Error: BLAST results file not found.

**Solution:**
Check that the file paths are correct and that the files exist in the specified location.

---

## 4. "Too many N bases" or QC Fails

**Error:**
Warning: X of the Y FASTA sequences failed QC ... too many N ...

**Solution:**
By default, ClassifyITS will flag sequences where >1% of bases are ambiguous (N).
Representative sequences should not have many (or any) N bases.

- It is recommended to run DADA2 denoise and chimera removal.
- If issues persist, consider running LULU to further filter out low-quality sequences.
- If using Illumina data, consider filtering out polyG tails (for example, using fastp).

Without careful examination, as many as (~1%) of representative sequences may fail to produce BLAST results or be unassigned at kingdom level—not because they are novel taxa, but because they are sequence errors.

---

## 5. BLAST Must Be Run Externally

**Note:**
ClassifyITS does NOT generate BLAST databases or run BLAST itself.
Users must generate BLAST output externally (see Data Preparation).
BLAST is well implemented and meant to be run on HPC, not within R.
It is a bad idea to try to run BLAST within R on large datasets (like UNITE), and it is not recommended to use the rBLAST package for this purpose.
Theoretically, you could run the R package rBLAST locally on a small dataset, but it is not recommended for large datasets like UNITE.

---

## 6. Still Stuck? Need Help?

For technical issues, bugs, or feedback:

Quinn Moon  
Email: qmoon@umich.edu

---