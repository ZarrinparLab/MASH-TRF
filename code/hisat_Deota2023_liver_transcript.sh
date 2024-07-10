#!/bin/bash
list_file="/mnt/zarrinpar/scratch/sfloresr/NASH_KF/transcriptomics/Panda_multiorgan_TRF/SRR_Acc_List-liver-2.txt"

path=/mnt/archive/Notebooks/sfloresr/20240109_RNAseq_Panda_multiomics
output_path=/mnt/zarrinpar/archive/Notebooks/sfloresr/20240109_RNAseq_Panda_multiomics/hisat/

while IFS= read -r item; do
    /mnt/zarrinpar/Pynchon/Notebooks/rrichter/scripts/hisat_stringtie.sh -s "$item" -o "$output_path" "$path/${item}_1.fastq.gz" "$path/${item}_2.fastq.gz"
done < "$list_file"
