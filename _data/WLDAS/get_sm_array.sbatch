#!/bin/bash
#SBATCH --job-name=sm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --array=2017
#SBATCH --output=get_sm_%a.log

## load modules
ml R/4.1.2
ml physics
ml proj gdal geos udunits

## execute script(s)
srun Rscript get_sm_array.R $SLURM_ARRAY_TASK_ID
