
###################################################################################################
## download_zenodo.R: This script downloads SPEAR IVT, AR, and precipitation files associated with 
## the paper "When Will Humanity Notice Its Influence on Atmospheric Rivers?" (Tseng et al., 2022)
## Created by Corinne Bowers 8/22/22
###################################################################################################

#### setup ########################################################################################

## load packages
library(inborutils)

## load helper functions
toNumber <- function(x) as.numeric(paste(x))

## define loop number
loop <- toNumber(commandArgs(trailingOnly=TRUE)[1])

#### zenodo data ##################################################################################
# https://ict.ipbes.net/ipbes-ict-guide/data-management/technical-guidelines/Zenodo

doi <- c(
  '10.5281/zenodo.6589996',
  '10.5281/zenodo.6590721', 
  '10.5281/zenodo.6591231', 
  '10.5281/zenodo.6591690', 
  '10.5281/zenodo.6591940', 
  '10.5281/zenodo.6600704', 
  '10.5281/zenodo.6600880', 
  '10.5281/zenodo.6600982', 
  '10.5281/zenodo.6609371', 
  '10.5281/zenodo.6612644', 
  '10.5281/zenodo.6612648', 
  '10.5281/zenodo.6613931', 
  '10.5281/zenodo.6629478', 
  '10.5281/zenodo.6629488', 
  '10.5281/zenodo.6633832', 
  '10.5281/zenodo.6633830')
download_zenodo(doi[loop], '/scratch/users/cbowers/GFDL')