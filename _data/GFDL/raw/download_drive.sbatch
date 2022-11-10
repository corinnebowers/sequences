#!/bin/bash
#SBATCH --job-name=drive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbowers@stanford.edu
#SBATCH --partition=cee
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --output=drive.log

## load modules
ml system rclone

## Google Drive: execute command
# https://rclone.org/googlecloudstorage/
# rclone copy -vv bowerscorinne:GFDL/LE_data . 
rclone copy -vv bowerscorinne:GFDL/filler .