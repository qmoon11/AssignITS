## Optional: Custom Taxonomy Cutoffs

Optionally, you can provide a custom taxonomy cutoffs file via the `cutoffs_file` parameter in the pipeline.

Format:
The file should be a CSV or TSV table with columns for each rank, name and its cutoffs (either e-value or percent identity), for example:

taxonomic_rank,evalue_cutoff,pident_cutoff
kingdom,1e-5,80
phylum,1e-10,85

How to use:
- Specify the path to your custom cutoffs file using `cutoffs_file` in your pipeline call:

    ITS_assignment(
        blast_file = "input/blast_results.tsv",
        rep_fasta = "input/representative_seqs.fasta",
        cutoffs_file = "input/taxonomy_cutoffs.csv"
    )

- If omitted, the package default cutoffs will be used.

Example file:
A sample file is provided at `inst/extdata/taxonomy_cutoffs.csv` in package source code.

Reference example (in CSV format):

    rank,Kingdom,Phylum,Class,Order,Family,Genus,Species,e-value-kingdom,e-value-phylum,e-value-class,e-value-order,e-value-family,e-value-genus,per-ident-genus,per-ident-species
    kingdom,Metazoa,,,,,,1e-45
    kingdom,Fungi,,,,,1e-50
    order,Fungi,,,,,,,1e-80 
    phylum,Fungi,Ascomycota,,,,,1e-70
    species,Fungi,Ascomycota,Saccharomycetes,Saccharomycetales,Saccharomycetaceae,Saccharomyces,,,,,,,,98

Tip:
Check out `inst/extdata/taxonomy_cutoffs.csv` in the package folder for the full list used in ClassifyITS example and recommended formatting.

ClassifyITS is built to look for specific taxonmic groups (Leotiomycetes) if it does not find that level, it will use the default for the level above (Ascomycota default class cutoffs). Do not worry common Asco classes are very much present.