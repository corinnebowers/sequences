
###################################################################################################
## get_drainage2.R: This script identifies failed API calls from get_drainage.R and retries them. 
## Created by Corinne Bowers 8/9/22
## Last updated 9/26/22
###################################################################################################

#### setup ########################################################################################
cat('loading setup information...\n')

## set working directory
setwd('/scratch/users/cbowers/sequences/')

## load functions & packages 
source('_data/setup.R')

## load loop information from first round
dx <- 100
n <- 12


#### load gages ###################################################################################
cat('loading gages...\n')
load('_data/streamflow/gages_0922.Rdata')


#### combine first-attempt files ##################################################################
cat('loading first-attempt drainage files...\n')
drainage_full <- 
  foreach (i = 1:n, .combine = 'c') %do% {
    load(paste0('_data/streamflow/drainage_api/drainage', i, '.Rdata'))
    index <- ((i-1)*dx+1):(dx*i)
    index <- index[index <= nrow(gages)] 
    drainage[index]
  }


#### identify failed attempts #####################################################################
cat('identifying failed attempts...\n')
bad <- drainage_full %>% lapply(is.null) %>% unlist %>% which
write.table(
  bad, 
  file = '_data/streamflow/drainage_api/bad.txt', 
  row.names = FALSE, col.names = FALSE)


#### retry failed attempts ########################################################################
cat('retrying failed attempts...\n')
pb <- txtProgressBar(min = 0, max = nrow(gages), style = 3)
for (i in bad) try({
  ## download json list
  h <- handle_setopt(new_handle())
  api_call <-
    paste0(
      'https://streamstats.usgs.gov/streamstatsservices/watershed.geojson?',
      'rcode=CA&',
      'xlocation=', st_coordinates(gages)[i,1], '&',
      'ylocation=', st_coordinates(gages)[i,2] ,'&',
      'crs=4269&',
      'includeparameters=false&includeflowtypes=false&',
      'includefeatures=true&simplify=false')
  api <- api_call %>%
    curl_download(tempfile(), handle = h) %>%
    fromJSON(simplifyVector = FALSE)
      
  ## convert to sf polygon
  coords <-
    api$featurecollection[[2]]$feature$features[[1]]$geometry$coordinates[[1]]
  drainage_full[[i]] <- 
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
save(drainage_full, file = '_data/streamflow/drainage_api/drainage_full.Rdata')
cat('\ndone!\n\n')    
