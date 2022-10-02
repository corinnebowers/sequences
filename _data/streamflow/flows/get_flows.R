
###################################################################################################
## get_flows.R: This script downloads daily streamflow timeseries from NWIS.
## Created by Corinne Bowers 8/9/22
## Last updated 9/28/22
###################################################################################################

#### setup ########################################################################################
cat('loading setup information...\n')

## set working directory
setwd('/scratch/users/cbowers/sequences/')

## load functions & packages 
source('_data/setup.R')

## load gages (filtered to only those with valid drainage geometries)
load('_data/streamflow/drainage_0928.Rdata')


#### get flows ####################################################################################

## define loop number
i <- toNumber(commandArgs(trailingOnly=TRUE)[1])
index <- ((i-1)*100+1):(100*i)
index <- index[index <= nrow(gages)]
cat(paste0('getting data for gages ', min(index), '-', max(index), '...\n'))

## download streamflow timeseries
flows <- 
  readNWISdata(
    sites = paste(gages$site_no[index]),
    parameterCd = unique(paste(gages$param_cd[index])),
    startDate = max(c(ymd('1980-01-01'), min(gages$begin_date[index]))), 
    endDate = min(c(today(), max(gages$end_date[index])))) %>% 
  renameNWISColumns


#### save out #####################################################################################
save(flows, file = paste0('_data/streamflow/flows/flow', i, '.Rdata'))
cat('done!\n\n') 
