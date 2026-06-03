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
#SBATCH --array=1-50

pwd; hostname; date

set -e

source ~/anaconda3/bin/activate birdman

echo Chunk $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX


TABLEID="NASH_NASH_FT_taxonomy_filtered.gg2.asv.counts"
TABLE="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/mashstool16s_preprocessed_20211020_ID_13785_gg2/"$TABLEID".biom"
SLURMS="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/birdman_outputs/slurm_out/"$TABLEID
OUTDIR="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/birdman_outputs/inferences/"$TABLEID
LOGDIR="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/dat/stool_mtb/birdman_outputs/logs/"$TABLEID
mkdir -p $SLURMS
mkdir -p $OUTDIR
mkdir -p $LOGDIR

echo Starting Python script...
time python /mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/code/stool_mtb/03-nash_birdman/nash_birdman_chunked.py \
    --table-path $TABLE \
    --inference-dir $OUTDIR \
    --num-chunks $SLURM_ARRAY_TASK_MAX \
    --chunk-num $SLURM_ARRAY_TASK_ID \
    --logfile "${LOGDIR}/chunk_${SLURM_ARRAY_TASK_ID}.log" && echo Finished Python script!

x
