# Inspecting Failed OTUs and Assignments

It is recommended to carefully examine:
1. **OTUs that failed the pipeline**
2. **Any OTUs that failed to assign at kingdom level**  
   (see list in `outputs/initial_assignments.csv` after classification)

These are likely indicators of sequencing errors or problematic input data.

---

# Notes on Taxonomic Assignment Using ClassifyITS

This pipeline is designed to assign taxonomy to fungal ITS sequences, as ITS is the universal fungal barcode and fungal-specific primers do not classify non-fungal taxa (such as plants).

**It is normal to observe some OTUs assigned to non-fungal kingdoms** (e.g., Plantae, Protista, etc.); these assignments are made based on kingdom-specific cutoffs. Such OTUs are likely non-fungal and can typically be discarded.

---

## Unclassified OTUs at Phylum or Class Level

A small percentage of taxa (up to approximately 10% in novel habitats) may be unclassified at the phylum or class level.  
This does **not automatically indicate discovery of a new order or phylum**.  
It is strongly recommended to manually check the BLAST results for these OTUs:

- If they have strong BLAST hits to known fungal taxa, they may represent novel taxa not well-represented in the database.
- Occasionally, errors in the UNITE database (e.g., an entry classified as `fungi_sp` at the phylum level) will prevent assignment, even if other BLAST results indicate a confident placement (such as Ascomycota).

```r
# Read in BLAST results and assignments
blast <- read.table("input/blast_results.tsv", sep = "\t", header = TRUE) #user provided BLAST results
assignments <- read.csv("outputs/initial_assignments.csv", stringsAsFactors = FALSE) #preliminary taxonomy assignments from ClassifyITS


# Subset assignments to fungal OTUs unclassified at class level
unclassified_class_otus <- assignments[
  assignments$kingdom == "Fungi" &
  (is.na(assignments$class) | assignments$class == "" | assignments$class == "Unclassified"),
  "qseqid"
]

# subset blast results to these OTUs
blast_sub <- blast[blast$qseqid %in% unclassified_class_otus, ]

# View the BLAST results for these OTUs
View(blast_sub)

# You can also export this subset for manual inspection in Excel
write.table(
  blast_sub, 
  "outputs/unclassified_class_blast_results.tsv", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)
```

---

## Manual Inspection

Depending on your research question and sample size, it is feasible to **subset the BLAST results to only those fungal taxa unclassified at phylum or class**, and manually inspect and assign taxonomy:

- This is the recommended approach, as ClassifyITS is designed to conservatively assign taxonomy and allow for further inspection at relevant taxonomic levels.
- Manual curation is especially valuable for ambiguous or novel taxa.

---

*Careful review and manual curation enhance the reliability of your fungal taxonomy assignments and support robust downstream analyses.*