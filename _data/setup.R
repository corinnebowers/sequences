
#### packages #####################################################################################

suppressPackageStartupMessages({
  require(sf)
  require(terra)
  require(raster)
  require(tigris) 
  require(censusapi)
  require(lubridate)
  require(ncdf4)
  require(rnoaa)
  require(mvtnorm)
  require(evd)
  require(caret)
  require(pracma)
  require(dataRetrieval)
  require(exactextractr)
  require(fitdistrplus)
  require(foreach)
  require(doSNOW)
  require(doRNG)
  require(parallel)
  require(tidyverse)
  require(ncdf4)
  require(abind)
  require(curl)
  require(httr)
  require(jsonlite)
  require(rvest)
  require(ggplot2)
  require(GGally)
  require(scales) 
  require(ggspatial)
  require(cowplot)
  require(ggridges)
  require(leaflet)
  require(mapboxapi)
  require(scico)
  require(extrafont)
  require(mapview)
  require(grid)
  require(gridExtra)
  require(ggpubr)
  require(gt)
  require(ggnewscale)
})

options(tigris_use_cache = TRUE) #tigris
Sys.setenv(CENSUS_KEY = 'f2e090156b02ced027d4ed756f82c9a3a1aa38c9') #censusapi
rnoaa_options(cache_messages = FALSE) #rnoaa


#### constants & utility functions ################################################################

## set parallel backend
cores <- 5
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

## useful helper functions

add_index <- function(n) (1:n)/(n+1)

toNumber <- function(x) as.numeric(paste(x))
strip <- function(x) unname(unlist(x))

Mean <- function(x) ifelse(all(is.na(x)), NA, mean(x, na.rm = TRUE))
Sum <- function(x) ifelse(all(is.na(x)), NA, sum(x, na.rm = TRUE))
Min <- function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE))
Max <- function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))
mode <- function(x) names(which.max(table(x)))[1]

setNA <- function(x, new) ifelse(is.na(x), new, x)
sum.na <- function(x) sum(is.na(x))

wateryear <- function(d) year(d) + ifelse(month(d) %in% 10:12, 1, 0)

positive <- function(x) {
  x[x < 0] <- 0
  x
}

## useful constants

mft <- 3.28084


#### themes & axes ################################################################################

scale_x_origin <- function(...) {
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,NA), ...) }
scale_y_origin <- function(...) {
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,NA), ...) }
geom_parity <- function() geom_abline(slope = 1, intercept = 0, linetype = 'dashed')

theme_set(
  theme_classic() + theme(
    text = element_text(family = 'Segoe UI'),
    axis.line = element_line(size = 0.35),
    axis.ticks = element_line(size = 0.35, color = 'black'),
    legend.key.size = unit(0.35, 'cm')))
theme_bw_custom <- 
  function() {
    theme_bw() + theme(
      text = element_text(family = 'Segoe UI'),
      panel.border = element_rect(fill = NA, color = 'black', size = 0.35),
      axis.ticks = element_line(size = 0.35, color = 'black'),
      legend.key.size = unit(0.35, 'cm'))}


#### colors #######################################################################################

baker <- c()
baker[1] <- rgb(56, 95, 150, maxColorValue = 255)
baker[2] <- rgb(207, 89, 33, maxColorValue = 255)
baker[3] <- rgb(158, 184, 219, maxColorValue = 255)
baker[4] <- rgb(231, 184, 0, maxColorValue = 255)
baker[5] <- rgb(128, 0, 0, maxColorValue = 255)

ggcolor <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

roma.colors <- scico(n = 5, palette = 'roma')


#### geospatial ###################################################################################

## function to plot rasters in ggplot
raster.df <- function(x) x %>% as.data.frame(xy = TRUE) %>% setNames(c('x', 'y', 'value'))

## define coordinate systems
wgs84 <- 4326
nad83 <- 4269
projected <- 6417

## define useful spatial geometries
wus <- states() %>%
  filter(STUSPS %in% c('OR','WA','CA','NV')) %>%
  st_transform(wgs84)
california <- counties(state = 'CA')
conus <- states() %>%
  filter(!(STUSPS %in% c('AK','HI','MP','GU','PR','AS','VI'))) %>%
  st_transform(wgs84)
sonoma <- tracts(state = 'CA', county = 'Sonoma') %>%
  st_union %>% st_sf %>%
  st_transform(wgs84)

## define global raster grid 
xres <- 0.625
yres <- 0.5
grid <- 
  raster(
    xmn = -180, xmx = 180, ymn = -90, ymx = 90, 
    resolution = c(xres,yres), crs = projection(wus))
lon <- raster.df(grid) %>% pull(x) %>% unique %>% sort
lon360 <- lon + 180
lat <- raster.df(grid) %>% pull(y) %>% unique %>% sort

## subset global grid to western US 
lon_subset <- lon[lon >= -130 & lon <= -110]
lat_subset <- lat[lat >= 30 & lat <= 50]

## identify grid cells that cover california
grid_ca <- grid %>% 
  raster::crop(california, snap = 'out') %>% 
  setValues(1:ncell(.))
grid_ca[] <- ifelse(rasterize(california, grid_ca, getCover = TRUE)[]>=0.05, grid_ca[], NA)
index_ca <- grid_ca %>% raster.df %>% filter(!is.na(value)) %>% pull(value)

## find out which county covers the majority of each grid cell
fact <- 100
grid_county <- grid %>% 
  crop(california, snap = 'out') %>% 
  disaggregate(fact) %>% 
  rasterize(california %>% mutate(county = toNumber(COUNTYFP)), ., field = 'county') %>% 
  aggregate(fact, 'modal')


#### watersheds ###################################################################################

## load HUC4 watersheds
huc4 <- st_read('D:/2-sequences/_data/WBD/WBDHU4.shp', quiet = TRUE)
huc4 <- huc4 %>% 
  mutate(huc4 = toNumber(huc4)) %>% 
  arrange(huc4)
colors.huc <- scico(11, palette = 'roma')[-6]

## calculate representative cells
grid_hucrep <- grid_ca
grid_hucrep[] <- 0
for (huc in 1:10) {
  ## choose based on HUC4 centroid
  id <- grid_ca %>% 
    raster::extract(huc4[huc,] %>% st_centroid %>% st_transform(crs(grid_ca)))
  grid_hucrep[id] <- huc+1800
}

## find out which HUC4 covers the majority of each cell
fact <- 100
grid_huc4 <- grid %>% 
  crop(california, snap = 'out') %>% 
  disaggregate(fact) %>% 
  rasterize(huc4 %>% mutate(huc4 = toNumber(huc4)), ., field = 'huc4') %>%
  aggregate(fact, 'modal')



#### functions: AR catalog ########################################################################

calculate_duration <- function(df, int = 3) {
  ## calculates AR duration
  #' df: dataframe with the columns ar, duration, & count
  #' int: integer measuring the time interval of the data (hours)
  
  counter <- 0
  for (ii in 1:nrow(df)) {
    if (df$ar[ii]==1) {
      if (ii==1) {
        counter <- counter+1
        df$duration[ii] <- int
      } else if (df$ar[ii-1]==0) {
        counter <- counter+1
        df$duration[ii] <- int
      } else {
        df$duration[ii] <- df$duration[ii-1]+int
      }
      df$count[ii] <- counter
    }
  }
  return(df)
}

assign_AR_cat <- function(IVT, duration) {
  ## calculates AR intensity category
  if (duration > 48) {
    return(case_when(IVT>=1000 ~ 5, IVT>=750 ~ 4, IVT>=500 ~ 3, TRUE ~ 2))
  } else if (duration > 24) {
    return(case_when(IVT>=1250 ~ 5, IVT>=1000 ~ 4, IVT>=750 ~ 3, IVT>=500 ~ 2, TRUE ~ 1))
  } else {
    return(case_when(IVT>=1250 ~ 4, IVT>=1000 ~ 3, IVT>=750 ~ 2, IVT>=500 ~ 1, TRUE ~ 0))
  }
}

lag <- function(x, agg, fun, align = 'right') {
  ## performs a running calculation across a certain number of rows
  #' x: vector of numbers 
  #' agg: integer specifying the number of rows to include in the running calculation
  #' fun: string specifying what calculation to perform (accepts sum, mean, min, max, events)
  #' align: should the calculated number be saved to the left, right, or center 
  #'   of the running interval?

  if (align == 'left') {
    offset <- agg
  } else if (align == 'center') {
    offset <- agg/2
  } else if (align == 'right') {
    offset <- 0
  } else stop('Unknown input for "align."')
  
  y <- rep(NA, length(x))
  if (fun == 'sum') {
    for (i in (agg+1):length(x)) y[i-offset] <- sum(x[(i-agg):(i-1)])
  } else if (fun == 'mean') {
    for (i in (agg+1):length(x)) y[i-offset] <- mean(x[(i-agg):(i-1)])
  } else if (fun == 'Mean') {
    for (i in (agg+1):length(x)) y[i-offset] <- Mean(x[(i-agg):(i-1)])
  } else if (fun == 'max') {
    for (i in (agg+1):length(x)) y[i-offset] <- max(x[(i-agg):(i-1)])
  } else if (fun == 'min') { 
    for (i in (agg+1):length(x)) y[i-offset] <- min(x[(i-agg):(i-1)])
  } else if (fun == 'events') {
    for (i in (agg+1):length(x)) y[i-offset] <- sum(unique(x[(i-agg):(i-1)])>0)
  } else stop('Unknown input for "agg."')
  y
}


###################################################################################################
