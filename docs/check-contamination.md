# Screening for Potential Contaminants

This guide describes how to identify possible contaminants in a fungal ITS amplicon dataset using the provided `potential-contamination.fasta` file and BLAST+.
This is not implemented in the package but I find it useful to check. 
---

## Step 1: Download the Contaminant FASTA File

Download the contaminant reference file from the github repository:

    wget https://github.com/qmoon/ClassifyITS/inst/extdata/potential-contamination.fasta


---

## Step 2: Build a BLAST Database from the Contaminant FASTA

Using BLAST+ (already installed from [Data Preparation](data-preparation.md)), turn the FASTA into a database:
*Again this will be fastest on HPC

    makeblastdb \
        -in potential-contamination.fasta \
        -dbtype nucl \
        -out contam_db

---

## Step 3: BLAST Rep-Seqs Against the Contaminant Database

Run BLASTn (megablast recommended for high similarity):

    blastn \
        -task megablast \
        -query rep-seqs/dna-sequences.fasta \
        -db contam_db \
        -out blast_contam.tsv \
        -word_size 28 \
        -reward 1 \
        -penalty -2 \
        -gapopen 0 \
        -gapextend 2 \
        -max_target_seqs 10 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
        -num_threads 8

---

## Step 4: Filter Results for High Identity and Sufficient Alignment Length

To find matches that are likely contaminants:

- Filter for percent identity (`pident`) ≥ 98%  
- Filter for alignment length (`length`) > 100 bp

Do this with `awk` or Excel, but here’s how with `awk`:

    awk '$3 >= 98 && $4 > 100' blast_contam.tsv > filtered_contaminants.tsv

Where:
- `$3` is the `pident` column (98 or higher = likely contaminant)
- `$4` is the `length` column (>100 bp = strong match)

---

## Step 5: Use Results

- Review `filtered_contaminants.tsv` for sequences potentially matching contaminants in the data.
- Remove or flag these sequences for downstream analysis.

---

## Notes

- Adjust the percent identity (98%) or alignment length (100bp) as needed for your project.
- The contaminant screening step can be completed before or after taxonomic classification.
- The provided `potential-contamination.fasta` is not exhaustive, so consider adding other known contaminants relevant to your environment or sequencing method. 
It currently contains ITS sequence from the type specimen of the obligatory human associated taxa.

| Species                      | Type Strain/Source            | GenBank Accession |
|------------------------------|-------------------------------|-------------------|
| Candida albicans             | CBS 562                       | MW284503.1        |
| Trichophyton rubrum          | CBS 392.58                    | NR_131330.1       |
| Trichophyton violaceum       | CBS 374.92                    | NR_144901.1       |
| Trichophyton soudanense      | IHEM 19751                    | NR_155948.1       |
| Trichophyton concentricum    | CBS 196.26                    | NR_144903.1       |
| Trichophyton schoenleinii    | CBS 458.59                    | NR_153261.1       |
| Trichophyton tonsurans       | CBS 496.48                    | NR_144891.1       |
| Saccharomyces cerevisiae     | CBS 1171                      | NR_111007.1       |
| Malassezia pachydermatis     | CBS 1879                      | NR_126114.1       |
| Malassezia globosa           | ATCC MYA-4612                 | NR_111475.1       |
| Malassezia sympodialis       | CBS 7222                      | NR_103583.1       |
| Malassezia restricta         | CBS 7877                      | NR_103585.1       |

---