#!/bin/bash
#SBATCH --chdir=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/
#SBATCH --partition=short
#SBATCH --mail-user="sfloresr@ucsd.edu"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --array=1-20

pwd; hostname; date

set -e

source ~/anaconda3/bin/activate birdman

echo Chunk $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX


TABLEID="04.taxonomy_filtered.asv.counts.g1pg2p"
TABLE="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/data/human_analysis/caussy_nafld_16s/"$TABLEID".biom"
SLURMS="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/results/human_analysis/caussy_nafld_16s/slurm_out/"$TABLEID
OUTDIR="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/results/human_analysis/caussy_nafld_16s/inferences/"$TABLEID
LOGDIR="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/MASHomics/results/human_analysis/caussy_nafld_16s/logs/"$TABLEID
mkdir -p $SLURMS
mkdir -p $OUTDIR
mkdir -p $LOGDIR

echo Starting Python script...
time python /home/sfloresr/projects/NALFD_ml/nafld_birdman/nafld_birdman_chunked.py \
    --table-path $TABLE \
    --inference-dir $OUTDIR \
    --num-chunks $SLURM_ARRAY_TASK_MAX \
    --chunk-num $SLURM_ARRAY_TASK_ID \
    --logfile "${LOGDIR}/chunk_${SLURM_ARRAY_TASK_ID}.log" && echo Finished Python script!

