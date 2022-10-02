
#### setup ########################################################################################

## set working directory
setwd('/scratch/users/cbowers/sequences/')

## load functions & packages
source('_data/setup.R')
source('_scripts/create_df_functions.R')

## set up parallel backend
cores <- parallel::detectCores()-2

## define custom function to convert .nc array output to raster bricks
array_to_raster <- function(x) {
  x %>% 
    aperm(c(2,1,3)) %>% 
    brick(
      xmn = min(lon_subset)-xres/2, xmx = max(lon_subset)+xres/2,
      ymn = min(lat_subset)-yres/2, ymx = max(lat_subset)+yres/2,
      crs = projection(wus)) %>%
    raster::flip('y')
}

## load AR & IVT data
load('_data/MERRA/ivt/ivt_merra_0928.Rdata')
load('_data/MERRA/ivt/ar_merra_0928.Rdata')


#### calculate AR frequency by category ###########################################################

## set up dataframe
ar_freq <- 
  matrix(ncol = 8, nrow = ncell(ivt_merra)) %>% 
  as.data.frame %>% 
  setNames(c('x', 'y', paste0('cat', 0:5)))

## loop over all grid cells
start <- Sys.time()
pb <- txtProgressBar(min = 0, max = ncell(ivt_merra), style = 3)
for (j in 1:ncell(ivt_merra)) {
  ## create AR event catalog
  index <- rowColFromCell(ar_merra, j)
  catalog <- 
    data.frame(
      ts = ts_merra,
      ar = c(ar_merra[index[1],index[2],])==1,
      ivt = c(ivt_merra[index[1],index[2],])) %>%
    create_catalog(., name = 'ar', interval = 3, cat = TRUE)
  
  ## summarize annual occurrence frequency of AR intensity categories
  freq <- catalog %>%
    group_by(cat) %>%
    summarize(freq = length(ivt_max) / length(years)) %>% 
    full_join(data.frame(cat = 0:5), by = 'cat') %>% 
    arrange(cat) %>% 
    pull(freq)
  
  ## save out
  ar_freq[j,] <- c(coordinates(ar_merra)[j,], freq)
  setTxtProgressBar(pb,j)
}
Sys.time() - start


#### save out #####################################################################################

save(ar_freq, file = '_data/MERRA/ivt/ar_freq_0929.Rdata')

