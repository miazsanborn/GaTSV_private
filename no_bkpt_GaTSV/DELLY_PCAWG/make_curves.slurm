#!/bin/bash

#SBATCH --job-name=gatsv-stats
#SBATCH --output=/data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/gatsv_%A_%a.out
#SBATCH --error=/data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/gatsv_%A_%a.err
#SBATCH --partition=normal
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --nodes=1

# Load R
module purge
module load R

# Move to working directory
cd /data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG

# Run R script with that line as an argument
Rscript /data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/DELLY_curves.R
