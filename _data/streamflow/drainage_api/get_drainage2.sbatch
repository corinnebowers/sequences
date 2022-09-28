#!/bin/bash
#SBATCH --job-name=drain2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --output=get_drainage_retry.log

## load modules
ml R/4.1.2
ml openssl curl
ml physics
ml proj gdal geos udunits

## execute script(s)
Rscript get_drainage2.R