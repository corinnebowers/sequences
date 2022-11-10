#!/bin/bash
#SBATCH --job-name=bc_setup
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --output=biascorr_setup.log

## load modules
ml R/4.1.2
ml physics
ml proj gdal geos udunits

## execute script(s)
Rscript biascorr_setup.R