#!/bin/bash
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=246G
#SBATCH --mail-type=END,FAIL
# (We omit variables that are set dynamically by R)

# 1. Load Modules
module load gfbf/2023b
module load R/4.4.1-gfbf-2023b
module load pandoc/2.19.2
module load libdeflate/1.19-GCCcore-13.2.0

# 2. Navigate to the project directory
# This ensures Rscript finds .Rprofile and activates renv automatically.
cd $HOME/github/ARscape

# 3. Run the R Script
# $1 = The path to slurm_payload.R
# $2 = The input fc path
# $3 = The input hfc path
# $4 = The output .rds path
Rscript "$1" "$2" "$3" "$4"
