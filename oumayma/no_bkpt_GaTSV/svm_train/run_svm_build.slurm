#!/bin/bash

#SBATCH --job-name=train-svm
#SBATCH --output=/data/beroukhim1/oumayma/no_bkpt_GaTSV/svm_train/logs/gatsv_%A_%a.out
#SBATCH --error=/data/beroukhim1/oumayma/no_bkpt_GaTSV/svm_train/logs/gatsv_%A_%a.err
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --nodes=1

# Load R
module load conda/latest
source activate r44

# Move to working directory
cd /data/beroukhim1/oumayma/no_bkpt_GaTSV/svm_train

# Run R script with that line as an argument
Rscript /data/beroukhim1/oumayma/no_bkpt_GaTSV/svm_train/svm_build.R

