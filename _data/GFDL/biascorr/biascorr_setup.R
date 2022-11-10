
###################################################################################################
## biascorr_setup.R: This script creates .Rdata files to support the file biascorr.Rmd.
## Created by Corinne Bowers 11/3/22
## Last updated 11/4/22
###################################################################################################

#### setup ########################################################################################
# cat('setting up...\n')

## set working directory
# setwd('/scratch/users/cbowers/sequences/')
# setwd('D:/2-sequences')

## get functions & packages
source('_data/setup.R')
source('_data/GFDL/biascorr/biascorr_functions.R')

## set up parallel backend
cores <- parallel::detectCores()-2

## load additional packages
suppressPackageStartupMessages(require(POT))


#### load metadata ################################################################################
cat('loading metadata...\n')

## MERRA
ts_merra <- 
  seq(ymd_hms('1980-1-1 00:00:00'), ymd_hms('2021-12-31 21:00:00'), by = '3 hours')
dates_merra <- ts_merra %>% as.Date %>% unique

## GFDL
load('_data/GFDL/gfdl_metadata.Rdata')


#### define constants ############################################################################
cat('defining constants...\n')

## define future decades for distribution fits
decades <- lapply(0:6, function(x) (2021:2030) + 10*x)
tridecades <- lapply(decades, function(x) (min(x)-10):(max(x)+10))

## define interpolation domains
ivt.domain <- c(1e-1, 5e3, 1e-1)
precip.domain <- c(1e-1, 500, 1e-2)

## define nonzero precipitation threshold
nonzero <- 0.1

## define length of data used for fitting
len <- data.frame(date = ts_hist) %>%
  filter(month(date) %in% c(10:12,1:3) & year(date) %in% 1981:2010) %>% nrow


#### load IVT #####################################################################################
cat('loading IVT from file...\n')

# start <- Sys.time()
# 
# ## MERRA
# load('_data/MERRA/ivt/ivt_merra_0928.Rdata')  # ivt_merra, ts_merra
# ivt_merra <- ivt_merra %>% crop(grid_ca)
# 
# ## GFDL historic
# ivt_gfdl <- foreach (i = 1:5) %do% {
#   load(paste0('_data/GFDL/raw/files/ivt_hist_ens0', i, '.Rdata'))
#   if (i==1) ivt_hist_ens01 else ivt_hist
# } %>% lapply(function(x) crop(x, grid_ca))
# 
# ## GFDL SSP 2-4.5
# ivt_ssp245 <- foreach (i = 1:5) %do% {
#   load(paste0('_data/GFDL/raw/files/ivt_ssp245_ens0', i, '.Rdata'))
#   if (i==1) ivt_ssp245_ens01 else ivt_ssp245
# } %>% lapply(function(x) crop(x, grid_ca))
# 
# ## GFDL SSP 5-8.5
# ivt_ssp585 <- foreach (i = 1:5) %do% {
#   load(paste0('_data/GFDL/raw/files/ivt_ssp585_ens0', i, '.Rdata'))
#   if (i==1) ivt_ssp585_ens01 else ivt_ssp585
# } %>% lapply(function(x) crop(x, grid_ca))
# 
# ## clean up
# rm(ivt_hist, ivt_hist_ens01, ivt_ssp245_ens01, ivt_ssp585_ens01)
# Sys.time() - start

## checkpoint 
# save(ivt_merra, ivt_gfdl, ivt_ssp245, ivt_ssp585, file = '_data/GFDL/raw/ivt_all.Rdata')
load('_data/GFDL/raw/ivt_all.Rdata')


#### load precipitation ###########################################################################
cat('loading precipitation from file...\n')

# ## MERRA
# load('_data/MERRA/precip/precip_24hr_0928.Rdata')
# precip_merra <- prcp %>% crop(grid_ca)
# 
# ## GFDL historic
# precip_gfdl <- foreach (i = 1:5) %do% {
#   load(paste0('_data/GFDL/raw/files/prcp_hist_ens0', i, '.Rdata'))
#   if (i==1) prcp_hist_ens01 else prcp_hist
# } %>% lapply(function(x) crop(x, grid_ca))
# 
# ## GFDL SSP 2-4.5
# precip_ssp245 <- foreach (i = 1:5) %do% {
#   load(paste0('_data/GFDL/raw/files/prcp_ssp245_ens0', i, '.Rdata'))
#   if (i==1) prcp_ssp245_ens01 else prcp_ssp245
# } %>% lapply(function(x) crop(x, grid_ca))
# 
# ## GFDL SSP 5-8.5
# precip_ssp585 <- foreach (i = 1:5) %do% {
#   load(paste0('_data/GFDL/raw/files/prcp_ssp585_ens0', i, '.Rdata'))
#   if (i==1) prcp_ssp585_ens01 else prcp_ssp585
# } %>% lapply(function(x) crop(x, grid_ca))
# 
# ## clean up 
# rm(prcp, prcp_hist, prcp_hist_ens01, prcp_ssp245, prcp_ssp245_ens01, prcp_ssp585, prcp_ssp585_ens01)

## checkpoint
# save(precip_merra, precip_gfdl, precip_ssp245, precip_ssp585, file = '_data/GFDL/raw/precip_all.Rdata')
load('_data/GFDL/raw/precip_all.Rdata')


#### mixture model: MERRA IVT #####################################################################
cat('calculating mixture model parameters for MERRA-2 IVT...\n')

# start <- Sys.time()
# pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
# ivt_merra_param <-
#   foreach (
#     i = 1:ncell(grid_ca),
#     .export =
#       c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt', 'calculate.rp100', 'calculate.aic'),
#     .packages = c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## calculate mixture model parameters for MERRA
#         values_merra <- ivt_merra %>%
#           raster::extract(i) %>% c %>%
#           data.frame(ts = ts_merra, ivt = .) %>%
#           group_by(date = as.Date(ts)) %>%
#           summarize(values = sample(ivt, 1)) %>%
#           filter(month(date) %in% c(10:12,1:3) & year(date) %in% 1981:2010) %>%
#           mutate(values = cbind(values,ivt.domain[1]) %>% apply(1,max)) %>%
#           pull(values)
#         par_merra <- tryCatch(
#           optimize.aic(values_merra, ivt.domain),
#           error = function(e) c(rep(0,6),e))
#         ## return results
#         list(values_merra, par_merra[1:6], par_merra[7])
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start

# ## troubleshoot
# data.frame(
#   converge = ivt_merra_param %>%
#     lapply(function(x) if (is.null(x)) NA else sum(x[[2]])!=0) %>% reduce(c),
#   ks.pass = ivt_merra_param %>%
#     lapply(function(x) if (is.null(x)) NA else x[[2]][7]) %>% reduce(c))

## checkpoint
# save(ivt_merra_param, file = '_data/GFDL/biascorr/params/ivt_merra_param.Rdata')
load('_data/GFDL/biascorr/params/ivt_merra_param.Rdata')


#### mixture model: historic IVT ##################################################################
cat('calculating mixture model parameters for GFDL historic IVT...\n')

# start <- Sys.time()
# pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
# ivt_hist_param <-
#   foreach (
#     i = 1:ncell(grid_ca),
#     .export =
#       c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt', 'calculate.rp100', 'calculate.aic'),
#     .packages = c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## calculate mixture model parameters for GFDL historic
#         values_hist <- matrix(nrow = len, ncol = 5)
#         par_hist <- matrix(nrow = 7, ncol = 5)
#         for (ens in 1:5) {
#           values_hist[,ens] <- ivt_gfdl[[ens]] %>%
#             raster::extract(i) %>% c %>%
#             data.frame(date = ts_hist, values = .) %>%
#             filter(month(date) %in% c(10:12,1:3) & year(date) %in% 1981:2010) %>%
#             mutate(values = cbind(values,ivt.domain[1]) %>% apply(1,max)) %>%
#             pull(values)
#           par_hist[,ens] <- tryCatch(
#             optimize.aic(values_hist[,ens], ivt.domain),
#             error = function(e) rep(0,7))
#         }
#         ## return results
#         list(values_hist, par_hist[1:6,], par_hist[7,])
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start

# ## troubleshoot
# ivt_hist_param %>%
#   lapply(function(x) if (is.null(x)) rep(NA,5) else apply(x[[2]], 2, function(x) sum(x)==0)) %>%
#   reduce(rbind) %>%
#   apply(1, Sum) %>% unname
# ivt_hist_param %>%
#   lapply(function(x) if (is.null(x)) rep(NA,5) else x[[3]]) %>%
#   reduce(rbind)

## checkpoint
# save(ivt_hist_param, file = '_data/GFDL/biascorr/params/ivt_hist_param.Rdata')
load('_data/GFDL/biascorr/params/ivt_hist_param.Rdata')


#### mixture model: MERRA precipitation ###########################################################
cat('calculating mixture model parameters for MERRA-2 precipitation...\n')

# start <- Sys.time()
# pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
# precip_merra_param <-
#   foreach (
#     i = 1:ncell(grid_ca),
#     .export =
#       c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt',
#         'calculate.rp100', 'calculate.aic', 'add_index'),
#     .packages =
#       c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## calculate mixture model parameters for MERRA
#         values_merra <- precip_merra %>%
#           raster::extract(i) %>% c %>%
#           data.frame(date = dates_merra, values = .) %>%
#           filter(month(date) %in% c(10:12,1:3) & year(date) %in% 1981:2010) %>%
#           mutate(values = ifelse(values>nonzero, values, NA)) %>%
#           pull(values)
#         par_merra <- tryCatch(
#           optimize.aic(values_merra[!is.na(values_merra)], precip.domain),
#           error = function(e) rep(0,7))
#         ## return results
#         list(values_merra, par_merra[1:6], par_merra[7])
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start

# ## troubleshoot
# cbind(
#   converge = precip_merra_param %>% 
#     lapply(function(x) if (is.null(x)) NA else sum(x[[2]])!=0) %>% reduce(c),
#   kspass = precip_merra_param %>% 
#     lapply(function(x) if (is.null(x)) NA else x[[3]]) %>% reduce(c))

## checkpoint
# save(precip_merra_param, file = '_data/GFDL/biascorr/params/precip_merra_param.Rdata')
load('_data/GFDL/biascorr/params/precip_merra_param.Rdata')


#### mixture model: historic precipitation ########################################################
cat('calculating mixture model parameters for GFDL historic precipitation...\n')

# start <- Sys.time()
# pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
# precip_hist_param <-
#   foreach (
#     i = 1:ncell(grid_ca),
#     .export =
#       c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt',
#         'calculate.rp100', 'calculate.aic', 'add_index'),
#     .packages =
#       c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## calculate mixture model parameters for GFDL historic
#         values_hist <- matrix(nrow = len, ncol = 5)
#         par_hist <- matrix(nrow = 7, ncol = 5)
#         for (ens in 1:5) {
#           values_hist[,ens] <- precip_gfdl[[ens]] %>%
#             raster::extract(i) %>% c %>%
#             data.frame(date = ts_hist, values = .) %>%
#             filter(month(date) %in% c(10:12,1:3) & year(date) %in% 1981:2010) %>%
#             mutate(values = ifelse(values>nonzero, values, NA)) %>%
#             pull(values)
#           par_hist[,ens] <- tryCatch(
#             optimize.aic(values_hist[!is.na(values_hist[,ens]),ens], precip.domain),
#             error = function(e) rep(0,7))
#         }
#         ## return results
#         list(values_hist, par_hist[1:6,], par_hist[7,])
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start

# ## troubleshoot
# precip_hist_param %>%
#   lapply(function(x) ifelse(is.null(x), NA, sum(apply(x[[2]], 2, sum)==0))) %>%
#   reduce(c)

## checkpoint
# save(precip_hist_param, file = '_data/GFDL/biascorr/params/precip_hist_param.Rdata')
load('_data/GFDL/biascorr/params/precip_hist_param.Rdata')


###################################################################################################
cat('done!\n\n')
