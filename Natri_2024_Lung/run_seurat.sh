#!/bin/bash
#SBATCH --job-name=lung_sc               # Job name
#SBATCH --ntasks=1                       # Number of tasks (processes)
#SBATCH --cpus-per-task=1                # Number of CPU cores per task
#SBATCH --mem=64G                        # Memory per node
#SBATCH --time=24:00:00                  # Time limit 

# Run your R script
R -e 'source("explore_seurat.R")'