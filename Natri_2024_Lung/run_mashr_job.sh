#!/bin/bash
#SBATCH --job-name=mashr_Natri           # Job name
#SBATCH --ntasks=1                       # Number of tasks (processes)
#SBATCH --cpus-per-task=4                # Number of CPU cores per task
#SBATCH --mem=64G                        # Memory per node
#SBATCH --time=48:00:00                  # Time limit 

# Run your R script
R -e 'source("4_run_mashr.R")'