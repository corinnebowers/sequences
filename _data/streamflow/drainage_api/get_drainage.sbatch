#!/bin/bash
#SBATCH --job-name=drainage
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --array=1-12
#SBATCH --output=get_drainage_%a.log

## load modules
ml R/4.1.2
ml openssl curl
ml physics
ml proj gdal geos udunits

## execute script(s)
srun Rscript get_drainage.R $SLURM_ARRAY_TASK_ID
