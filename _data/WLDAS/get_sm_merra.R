
###################################################################################################
## get_sm.R: This script retrieves soil moisture data from NASA WLDAS.
## Created by Corinne Bowers 10/21/22
## Last updated 10/27/22
###################################################################################################

#### setup ########################################################################################
cat('setting up...\n') 

## set working directory
setwd('/scratch/users/cbowers/sequences/')

## load packages & functions
source('_data/setup.R')

## set up parallel backend
cores <- parallel::detectCores()-2
cat(paste0('  cores: ', cores, '\n'))

## define WLDAS metadata
lat_wldas <- seq(25.065, 52.925, 0.01)
lon_wldas <- seq(-124.925, -89.025, 0.01)

## define download URL
url <- 'https://portal.nccs.nasa.gov/datashare/WLDAS/wldas_domain/'

## define function to produce soil moisture raster 
wldas_to_raster <- function(x) {
  ## transform to raster
  sm <- x %>% 
    aperm(c(2,1,3)) %>% 
    brick(
      xmn = min(lon_wldas)-0.01/2, xmx = max(lon_wldas)+0.01/2,
      ymn = min(lat_wldas)-0.01/2, ymx = max(lat_wldas)+0.01/2,
      crs = projection(wus)) %>%
    raster::flip('y')
  
  ## convert to mm/m
  sm <- 1e3 * (sm[[1]]*0.1 + sm[[2]]*0.3 + sm[[3]]*0.6)
  
  ## resample to MERRA-2 grid
  sm <- sm %>% raster::resample(grid_ca)
  return(sm)
}


#### get soil moisture ############################################################################
cat('downloading soil moisture...\n') 

start <- Sys.time()
pb <- txtProgressBar(min = 0, max = 42*12, style = 3)
cl <- parallel::makeCluster(cores)
registerDoSNOW(cl)
sm <- 
  foreach (yr = 2010:2021, .inorder = FALSE) %:%   ### 1980:2021
  foreach (
    mo = 1:12, 
    .combine = 'c',
    .inorder = FALSE,
    .packages = c('raster', 'ncdf4', 'foreach', 'rvest', 'lubridate', 'tidyverse'),
    .export = 'wldas_to_raster',
    .options.snow = opts) %dopar% {
      files <- paste0(url, yr, str_pad(mo,2,pad ='0')) %>% 
        read_html %>% html_elements('a') %>% html_text()
      files <- files[grepl('.nc', files)]
    
      options(timeout = 3600)
      foreach (file = files) %do% {
        ## download file
        filepath <- paste0('_data/WLDAS/files/', file)
        try({
          if (!file.exists(filepath)) {
            download.file(
              url = paste0(url, yr, str_pad(mo,2,pad ='0'), '/', file),
              destfile = filepath)
          }     

          sm_name <- paste0('_data/WLDAS/files_sm/', str_remove(file, '\\.nc'), '.tif')  
          if (!file.exists(sm_name)) {
            nc <- nc_open(filepath)
            sm <- ncvar_get(nc, 'SoilMoist_tavg') %>% wldas_to_raster(.)
            nc_close(nc)
            writeRaster(sm, filename = sm_name)
          }
        })
        return(0)
      }
    }
stopCluster(cl)
cat('\n')
Sys.time() - start


###################################################################################################
cat('done!\n\n')


