
###################################################################################################
## get_drainage.R: This script downloads drainage area shapefiles from the StreamStats API.
## Created by Corinne Bowers 8/9/22
## Last updated 9/26/22
###################################################################################################

#### setup ########################################################################################
cat('loading setup information...\n')

## set working directory
setwd('/scratch/users/cbowers/sequences/')

## load functions & packages 
source('_data/setup.R')

## define loop number
loop <- toNumber(commandArgs(trailingOnly=TRUE)[1])
dx <- 100


#### load gages ###################################################################################
cat('loading gages...\n')
load('_data/streamflow/gages_0922.Rdata')


#### get drainage area for each gage ##############################################################
cat('calculating drainage area...\n')

### initialize loop
drainage <- list()
index <- ((loop-1)*dx+1):(dx*loop)
index <- index[index <= nrow(gages)]

### start loop
pb <- txtProgressBar(min = min(index), max = max(index), style = 3)
for (i in index) try({
  ## download json list
  h <- handle_setopt(new_handle())
  api_call <-
    paste0(
      'https://streamstats.usgs.gov/streamstatsservices/watershed.geojson?',
      'rcode=CA&',
      'xlocation=', st_coordinates(gages)[i,1], '&',
      'ylocation=', st_coordinates(gages)[i,2] ,'&',
      'crs=4269&',
      'includeparameters=false&',
      'includeflowtypes=false&',
      'includefeatures=true&',
      'simplify=false')
  api <- api_call %>%
    curl_download(tempfile(), handle = h) %>%
    fromJSON(simplifyVector = FALSE)

  ## convert to sf polygon
  coords <-
    api$featurecollection[[2]]$feature$features[[1]]$geometry$coordinates[[1]]
  drainage[[i]] <- 
    cbind(
      x = coords %>% lapply(function(x) x[[1]]) %>% unlist,
      y = coords %>% lapply(function(x) x[[2]]) %>% unlist) %>%
    list %>% st_polygon %>%
    st_sfc(crs = 4269) %>%
    st_sf %>%
    mutate(site_no = gages$site_no[i]) %>%
    select(site_no, geometry)
  
  ## report progress
  setTxtProgressBar(pb, i)
})


#### save out #####################################################################################
save(drainage, file = paste0('_data/streamflow/drainage_api/drainage', loop, '.Rdata'))
cat('\ndone!\n\n')
