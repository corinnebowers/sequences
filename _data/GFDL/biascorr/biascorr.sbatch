#!/bin/bash
#SBATCH --job-name=biascorr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --output=biascorr.log

## load modules
ml R/4.1.2
ml physics
ml proj gdal geos udunits

## execute script(s)
Rscript biascorr.R