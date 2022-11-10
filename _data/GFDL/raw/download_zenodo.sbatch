#!/bin/bash
#SBATCH --job-name=zenodo
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --array=1-16
#SBATCH --output=zenodo_%a.log

## load modules
ml R/4.1.2
ml physics proj gdal geos udunits

## Zenodo: execute R script
srun Rscript download_zenodo.R $SLURM_ARRAY_TASK_ID