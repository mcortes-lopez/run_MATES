#!/bin/bash
#SBATCH --mail-user=mcorteslopez@nygenome.org
#SBATCH --mail-type=ALL
#SBATCH --job-name=MATES_long_all_files_inclusive
#SBATCH --mem=48G
#SBATCH --partition=cpu
#SBATCH --output=logs/%x.log
#SBATCH --error=logs/%x.e
#SBATCH --cpus-per-task=12



CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate mates_env

python mates_lr.py