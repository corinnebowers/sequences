#!/bin/bash
#SBATCH --job-name=biascorr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --array=3,5,6
#SBATCH --output=biascorr%a.log

## load modules
ml R/4.1.2
ml physics
ml proj gdal geos udunits

## execute script(s)
Rscript biascorr_array.R $SLURM_ARRAY_TASK_ID