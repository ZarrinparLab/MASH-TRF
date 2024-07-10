#!/bin/bash

# File containing SRR accession numbers (one per line)
SRR_FILE="/mnt/zarrinpar/scratch/sfloresr/NASH_KF/transcriptomics/list_liver2_files_clean.txt"

# Output directory for downloaded FASTQ files
OUTPUT_DIR="/mnt/zarrinpar/archive/Notebooks/sfloresr/20240109_RNAseq_Panda_multiomics/"

# Loop through each SRR accession in the file and run fastq-dump
while IFS= read -r SRR; do
    fasterq-dump "$SRR" -p --outdir "$OUTPUT_DIR" --skip-technical --split-3
    gzip "$OUTPUT_DIR/${SRR}_1.fastq"
    gzip "$OUTPUT_DIR/${SRR}_2.fastq"
done < "$SRR_FILE"
