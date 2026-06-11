#!/bin/bash
#SBATCH --job-name=gatsv-classify
#SBATCH --output=/data/beroukhim1/oumayma/logs/gatsv_%A_%a.out
#SBATCH --error=/data/beroukhim1/oumayma/logs/gatsv_%A_%a.err
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=1-885

module purge
module load R

# Get the specific patient prefix for this array task
PREFIX_LIST="/data/beroukhim1/oumayma/PCAWG_metadata/patient_prefixes.txt"
PATIENT_PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PREFIX_LIST)

# Run the R script, passing the prefix as an argument
Rscript /data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/delly_process_classify.R "$PATIENT_PREFIX"
