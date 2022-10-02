#!/bin/bash
#SBATCH --job-name=freqbycat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --output=ar_freq_by_cat.log

## load modules
ml R/4.1.2
ml openssl curl
ml physics
ml proj gdal geos udunits

## execute script(s)
Rscript ar_freq_by_cat.R
