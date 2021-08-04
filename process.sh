#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --array=1-125
#SBATCH --cpus-per-task=32
#SBATCH --output=Rexample1.out
#SBATCH --mail-user=ssule25@umd.edu
#SBATCH --mail-type=ALL

# module load R/3.6.1
# firstarg="$1"
# Rscript trial.R $firstarg
 
module load R/3.6.1
# Rscript trial.R $firstarg
Rscript treecomputer.R $SLURM_ARRAY_TASK_ID
# Rscript treecomputer.R $firstarg

