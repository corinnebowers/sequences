#!/bin/bash
#SBATCH --job-name=sm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --output=get_sm.log

## load modules
ml R/4.1.2
ml physics
ml proj gdal geos udunits

## execute script(s)
Rscript get_sm.R
