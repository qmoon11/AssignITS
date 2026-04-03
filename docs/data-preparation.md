## Data Preparation

ClassifyITS does not run NCBI BLAST searches or build reference databases. These steps are already well implemented and should be run on HPC (high performance computing cluster). You must generate BLAST results externally and provide them as input. Below is a typical workflow to prepare your inputs (run on HPC):

1. Download BLAST+ (Linux example):

    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
    tar -xzf ncbi-blast-2.17.0+-x64-linux.tar.gz
    export PATH=/your/path/to/ncbi-blast-2.17.0+/bin:$PATH

2. Download the latest UNITE eukaryotic ITS database from https://unite.ut.ee/
It is recommended to use the latest UNITE eukaryotic ITS database for best results. Download the "General FASTA release" that matches your needs. For example:  
   
    UNITE/euk/sh_general_release_dynamic_s_all_19.02.2025.fasta

3. Build the BLAST database:

    makeblastdb \
      -in UNITE/euk/sh_general_release_dynamic_s_all_19.02.2025.fasta \
      -dbtype nucl \
      -out BLAST/uniteITSEuk_db

4. Run MegaBLAST search (example parameters):

    blastn \
      -task megablast \
      -query rep-seqs/dna-sequences.fasta \
      -db BLAST/uniteITSEuk_db \
      -out BLAST/megablast_ITS2.tsv \
      -word_size 28 \
      -reward 1 \
      -penalty -2 \
      -gapopen 0 \
      -gapextend 2 \
      -max_target_seqs 10 \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
      -num_threads 8

This will generate a tab-delimited BLAST output file (for example, BLAST/megablast_ITS2.tsv) that can be used as input for ClassifyITS. It will contain the top 10 hits for each query sequence, along with relevant information such as percent identity, alignment length, e-value, and taxonomic information from the UNITE database.

**Required BLAST output format:**  
Make sure you use:  
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

Do not change the column order—this is required for ClassifyITS to work!