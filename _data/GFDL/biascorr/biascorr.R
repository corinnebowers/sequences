
###################################################################################################

#### setup ########################################################################################
cat('setting up...\n')

## set working directory
# setwd('D:/2-sequences/')
setwd('/scratch/users/cbowers/sequences/')

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


#### mixture model: MERRA IVT ###############################################################################
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
#     .packages =
#       c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## calculate mixture model parameters for MERRA
#         values_merra <-
#           ivt_merra %>%
#           raster::extract(i) %>% c %>%
#           data.frame(ts = ts_merra, ivt = .) %>%
#           group_by(date = as.Date(ts)) %>%
#           summarize(values = sample(ivt, 1)) %>%
#           filter(month(date) %in% c(10:12,1:3)) %>%
#           filter(year(date) %in% 1981:2010) %>%
#           mutate(values = cbind(values,ivt.domain[1]) %>% apply(1,max)) %>%
#           pull(values)
#         par_merra <- tryCatch(
#           optimize.aic(values_merra, ivt.domain)$par, 
#           error = function(e) rep(0,6))
# 
#         ## return results
#         list(values_merra, par_merra)
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start
# 
# ## checkpoint
# save(ivt_merra_param, file = '_data/GFDL/biascorr/params/ivt_merra_param.Rdata')
load('_data/GFDL/biascorr/params/ivt_merra_param.Rdata')


#### mixture model: historic IVT ############################################################################
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
#     .packages =
#       c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## calculate mixture model parameters for GFDL historic
#         len <- data.frame(date = ts_hist) %>% 
#           filter(month(date) %in% c(10:12,1:3)) %>% 
#           filter(year(date) %in% 1981:2010) %>% 
#           nrow
#         values_hist <- matrix(nrow = len, ncol = 5)
#         par_hist <- matrix(nrow = 6, ncol = 5)
#         for (ens in 1:5) {
#           values_hist[,ens] <- 
#             ivt_gfdl[[ens]] %>% 
#             raster::extract(i) %>% c %>%
#             data.frame(date = ts_hist, values = .) %>% 
#             filter(month(date) %in% c(10:12,1:3)) %>% 
#             filter(year(date) %in% 1981:2010) %>% 
#             mutate(values = cbind(values,1) %>% apply(1,max)) %>% 
#             pull(values)
#           par_hist[,ens] <- tryCatch(
#             optimize.aic(values_hist, ivt.domain)$par, 
#             error = function(e) rep(0,6))
#         }
#         
#         ## return results
#         list(values_hist, par_hist)
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start
# 
# ## checkpoint
# save(ivt_hist_param, file = '_data/GFDL/biascorr/params/ivt_hist_param.Rdata')
load('_data/GFDL/biascorr/params/ivt_hist_param.Rdata')


#### mixture model: MERRA precipitation #####################################################################
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
#         values_merra <- 
#           precip_merra %>% 
#           raster::extract(i) %>% c %>%
#           data.frame(date = dates_merra, values = .) %>% 
#           filter(month(date) %in% c(10:12,1:3)) %>% 
#           filter(year(date) %in% 1981:2010) %>% 
#           mutate(values = cbind(values,precip.domain[1]) %>% apply(1,max)) %>%
#           pull(values)
#         par_merra <- tryCatch(
#           optimize.aic(values_merra, precip.domain)$par, 
#           error = function(e) rep(0,6))
# 
#         ## return results
#         list(values_merra, par_merra)
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start
# 
# ## checkpoint
# save(precip_merra_param, file = '_data/GFDL/biascorr/params/precip_merra_param.Rdata')
load('_data/GFDL/biascorr/params/precip_merra_param.Rdata')


#### mixture model: historic precipitation ##################################################################
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
#         len <- data.frame(date = ts_hist) %>% 
#           filter(month(date) %in% c(10:12,1:3)) %>% 
#           filter(year(date) %in% 1981:2010) %>% 
#           nrow
#         values_hist <- matrix(nrow = len, ncol = 5)
#         par_hist <- matrix(nrow = 6, ncol = 5)
#         for (ens in 1:5) {
#           values_hist[,ens] <- 
#             precip_gfdl[[ens]] %>% 
#             raster::extract(i) %>% c %>%
#             data.frame(date = ts_hist, values = .) %>% 
#             filter(month(date) %in% c(10:12,1:3)) %>% 
#             filter(year(date) %in% 1981:2010) %>% 
#             pull(values)
#           par_hist[,ens] <- tryCatch(
#             optimize.aic(values_hist[,ens], precip.domain)$par, 
#             error = function(e) rep(0,6))
#         }
#         
#         ## return results
#         list(values_hist, par_hist)
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start
# 
# ## checkpoint
# save(precip_hist_param, file = '_data/GFDL/biascorr/params/precip_hist_param.Rdata')
load('_data/GFDL/biascorr/params/precip_hist_param.Rdata')


#### bias correction: historic IVT ######################################################################
cat('bias correcting IVT...\n')
cat('...historic\n')

# start <- Sys.time()
# pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
# ivt_hist_sdm <-
#   foreach (
#     i = 1:ncell(grid_ca),
#     .export =
#       c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt', 'calculate.rp100', 'calculate.aic'),
#     .packages =
#       c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## load MERRA mixture model parameters
#         values_merra <- ivt_merra_param[[i]][[1]]
#         par_merra <- ivt_merra_param[[i]][[2]]
# 
#         ## load GFDL historic mixture model parameters
#         values_hist <- ivt_hist_param[[i]][[1]]
#         par_hist <- ivt_hist_param[[i]][[2]]
# 
#         ## apply SDM to bias-correct GFDL historic
#         unbiased <-
#           foreach (ens = 1:5) %do% {
#             df_hist <-
#               ivt_gfdl[[ens]] %>%
#               raster::extract(i) %>% c %>%
#               data.frame(date = ts_hist, values = .) %>%
#               mutate(values = cbind(values,1) %>% apply(1,max))
#             if (sum(par_merra) > 0 & sum(par_hist[,ens]) > 0) {
#               sdm2(
#                 par_future = par_hist[,ens], data = df_hist,
#                 par_merra = par_merra, values_merra = values_merra,
#                 par_hist = par_hist[,ens], values_hist = values_hist[,ens],
#                 domain = ivt.domain) %>%
#                 select(date, unbiased)
#             } else {
#               df_hist %>% select(date) %>% mutate(empty = NA)
#             }
#           } %>% reduce(full_join, by = 'date') %>%
#           arrange(date) %>% select(-date) %>%
#           setNames(paste0('ens0', 1:5))
#         return(unbiased)
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start
# 
# ## checkpoint
# save(ivt_hist_sdm, file = '_data/GFDL/biascorr/files/ivt_hist_sdm.Rdata')


#### bias correction: SSP2-4.5 IVT ######################################################################
cat('...SSP 2-4.5\n')

start <- Sys.time()
pb <- txtProgressBar(min = 0, max = length(decades)*ncell(grid_ca), style = 3)
cl <- parallel::makeCluster(cores)
registerDoSNOW(cl)
ivt_ssp245_sdm <-
  foreach (
    i = 1:ncell(grid_ca),
    .export =
      c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt', 'calculate.rp100', 'calculate.aic'),
    .packages =
      c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
    .options.snow = opts) %:%
  foreach (dec = length(decades), .combine = 'rbind') %dopar% {
    if (i %in% index_ca) {
      ## load MERRA mixture model parameters
      values_merra <- ivt_merra_param[[i]][[1]]
      par_merra <- ivt_merra_param[[i]][[2]]

      ## load GFDL historic mixture model parameters
      values_hist <- ivt_hist_param[[i]][[1]]
      par_hist <- ivt_hist_param[[i]][[2]]

      ## apply SDM to bias-correct GFDL future
      unbiased_ssp245 <-
        foreach (ens = 1:5) %do% {
          df_ssp245 <-
            ivt_ssp245[[ens]] %>%
            raster::extract(i) %>% c %>%
            data.frame(date = ts_future, values = .) 
          values_ssp245 <- df_ssp245 %>%
            filter(month(date) %in% c(10:12,1:3)) %>%
            filter(year(date) %in% tridecades[[dec]]) %>%
            mutate(values = cbind(values,1) %>% apply(1,max)) %>% 
            pull(values)
          par_ssp245 <- tryCatch(
            optimize.aic(values_ssp245, ivt.domain)$par,
            error = function(e) rep(0,6))

          if (sum(par_merra) > 0 & sum(par_hist[,ens]) > 0 & sum(par_ssp245) > 0) {
            sdm2(
              par_future = par_ssp245,
              data = df_ssp245 %>% filter(year(date) %in% decades[[dec]]),
              par_merra = par_merra, values_merra = values_merra,
              par_hist = par_hist[,ens], values_hist = values_hist[,ens],
              domain = ivt.domain) %>%
              select(date, unbiased)
          } else {
            df_ssp245 %>% select(date) %>% mutate(empty = NA)
          }
        } %>% reduce(full_join, by = 'date') %>%
        arrange(date) %>% #select(-date) %>%
        setNames(c('date', paste0('ens0', 1:5)))
      return(unbiased_ssp245)
    } else NULL
  }
stopCluster(cl)
Sys.time() - start

## checkpoint
save(ivt_ssp245_sdm, file = '_data/GFDL/biascorr/files/ivt_ssp245_sdm.Rdata')


#### bias corr: SSP5-8.5 IVT ######################################################################
cat('...SSP 5-8.5\n')

# start <- Sys.time()
# pb <- txtProgressBar(min = 0, max = length(decades)*ncell(grid_ca), style = 3)
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
# ivt_ssp585_sdm <- 
#   foreach (dec = 1:length(decades), .combine = 'rbind') %:%
#   foreach (
#     i = 1:ncell(grid_ca),
#     .export = 
#       c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt', 'calculate.rp100', 'calculate.aic'),
#     .packages = 
#       c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## load MERRA mixture model parameters
#         values_merra <- ivt_merra_param[[i]][[1]]
#         par_merra <- ivt_merra_param[[i]][[2]]
#         
#         ## load GFDL historic mixture model parameters
#         values_hist <- ivt_hist_param[[i]][[1]]
#         par_hist <- ivt_hist_param[[i]][[2]]
#         
#         ## apply SDM to bias-correct GFDL future
#         unbiased_ssp585 <- 
#           foreach (ens = 1:5) %do% {
#             df_ssp585 <- 
#               ivt_ssp585[[ens]] %>% 
#               raster::extract(i) %>% c %>%
#               data.frame(date = ts_future, values = .) %>% 
#               filter(month(date) %in% c(10:12,1:3)) %>% 
#               filter(year(date) %in% tridecades[[dec]]) %>% 
#               mutate(values = cbind(values,1) %>% apply(1,max))
#             par_ssp585 <- tryCatch(
#               optimize.aic(df_ssp585$values, ivt.domain)$par,
#               error = function(e) rep(0,6))
#             
#             if (sum(par_merra) > 0 & sum(par_hist[,ens]) > 0 & sum(par_ssp585) > 0) {
#               sdm2(
#                 par_future = par_ssp585, 
#                 data = df_ssp585 %>% filter(year(date) %in% decades[[dec]]),
#                 par_merra = par_merra, values_merra = values_merra, 
#                 par_hist = par_hist[,ens], values_hist = values_hist[,ens],
#                 domain = ivt.domain) %>%
#                 select(date, unbiased)
#             } else {
#               df_ssp585 %>% select(date) %>% mutate(empty = NA)
#             }
#           } %>% reduce(full_join, by = 'date') %>% 
#           arrange(date) %>% select(-date) %>% 
#           setNames(paste0('ens0', 1:5))
#       return(unbiased_ssp585)
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start
# 
# ## checkpoint
# save(ivt_ssp585_sdm, file = '_data/GFDL/biascorr/files/ivt_ssp585_sdm.Rdata')


#### bias corr: historic precipitation ############################################################
cat('bias correcting precipitation...\n')
cat('...historic\n')

# start <- Sys.time()
# pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
# precip_hist_sdm <- 
#   foreach (
#     i = 1:ncell(grid_ca),
#     .export = 
#       c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt', 
#         'calculate.rp100', 'calculate.aic', 'add_index'),
#     .packages = 
#       c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## load MERRA mixture model parameters
#         values_merra <- precip_merra_param[[i]][[1]]
#         par_merra <- precip_merra_param[[i]][[2]]
#         rd_merra <- sum(values_merra > nonzero)/length(values_merra)
# 
#         ## load GFDL historic mixture model parameters
#         values_hist <- precip_hist_param[[i]][[1]]
#         par_hist <- precip_hist_param[[i]][[2]]
#         rd_hist <- apply(values_hist, 2, function(x) sum(x > nonzero))/nrow(values_hist)
#         
#         ## apply SDM to bias-correct GFDL historic
#         unbiased <- 
#           foreach (ens = 1:5) %do% {
#             df_hist <- 
#               precip_gfdl[[ens]] %>% 
#               raster::extract(i) %>% c %>%
#               data.frame(date = ts_hist, values = .) %>% 
#               filter(month(date) %in% c(10:12,1:3)) %>% 
#               filter(year(date) %in% 1981:2010)
#             rd_hist <- sum(values_hist[,ens] > nonzero)/nrow(values_hist)
#             
#             ## apply rainy-day correction
#             if (sum(par_merra) > 0 & sum(par_hist[,ens]) > 0) {
#               temp <- sdm2(
#                 par_future = par_hist[,ens], data = df_hist,
#                 par_merra = par_merra, values_merra = values_merra, 
#                 par_hist = par_hist[,ens], values_hist = values_hist[,ens],
#                 domain = precip.domain) %>%
#                 select(date, unbiased)
#               rd <- sum(temp$unbiased > nonzero)/nrow(temp)
#               rd_scaled <- rd * rd_merra / rd_hist
#               temp %>% 
#                 arrange(desc(unbiased)) %>% 
#                 mutate(
#                   position = add_index(nrow(.)),
#                   unbiased_rd = ifelse(position > rd_scaled, 0, unbiased)) %>% 
#                 select(date, unbiased_rd)
#             } else {
#               df_hist %>% select(date) %>% mutate(empty = NA)
#             }
#           } %>% reduce(full_join, by = 'date') %>% 
#           arrange(date) %>% select(-date) %>% 
#           setNames(paste0('ens0', 1:5))
#         return(unbiased) 
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start
# 
# ## checkpoint
# save(precip_hist_sdm, file = '_data/GFDL/biascorr/files/precip_hist_sdm.Rdata')


#### bias corr: SSP2-4.5 precipitation ############################################################
cat('...SSP 2-4.5\n')

# start <- Sys.time()
# pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
# precip_ssp245_sdm <- 
#   foreach (dec = length(decades), .combine = 'rbind') %:% 
#   foreach (
#     i = 1:ncell(grid_ca),
#     .export = 
#       c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt', 
#         'calculate.rp100', 'calculate.aic', 'add_index'),
#     .packages = 
#       c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## load MERRA mixture model parameters
#         values_merra <- precip_merra_param[[i]][[1]]
#         par_merra <- precip_merra_param[[i]][[2]]
#         rd_merra <- sum(values_merra > nonzero)/length(values_merra)
# 
#         ## load GFDL historic mixture model parameters
#         values_hist <- precip_hist_param[[i]][[1]]
#         par_hist <- precip_hist_param[[i]][[2]]
#         rd_hist <- apply(values_hist, 2, function(x) sum(x > nonzero))/nrow(values_hist)
#         
#         unbiased_ssp245 <- 
#           foreach (ens = 1:5) %do% {
#             ## calculate mixture model parameters for GFDL future
#             df_ssp245 <- 
#               precip_ssp245[[ens]] %>% 
#               raster::extract(i) %>% c %>%
#               .[1:length(ts_future)] %>% 
#               data.frame(date = ts_future, values = .) %>% 
#               filter(month(date) %in% c(10:12,1:3)) %>% 
#               filter(year(date) %in% tridecades[[dec]])
#             par_ssp245 <- tryCatch(
#               optimize.aic(df_ssp245$values, precip.domain)$par,
#               error = function(e) rep(0,6))
#             
#             if (sum(par_merra) > 0 & sum(par_hist[,ens]) > 0 & sum(par_ssp245) > 0) {
#               ## apply SDM to bias-correct GFDL future
#               temp <- sdm2(
#                 par_future = par_ssp245, 
#                 data = df_ssp245 %>% filter(year(date) %in% decades[[dec]]),
#                 par_merra = par_merra, values_merra = values_merra, 
#                 par_hist = par_hist[,ens], values_hist = values_hist[,ens],
#                 domain = ivt.domain) %>%
#                 select(date, unbiased)
#           
#               ## apply rainy-day correction
#               rd <- sum(temp$unbiased > nonzero)/nrow(temp)
#               rd_scaled <- rd * rd_merra / rd_hist
#               temp %>% 
#                 arrange(desc(unbiased)) %>% 
#                 mutate(
#                   position = add_index(nrow(.)),
#                   unbiased_rd = ifelse(position > rd_scaled, 0, unbiased)) %>% 
#                 select(date, unbiased_rd)
#             } else {
#               df_ssp245 %>% select(date) %>% mutate(empty = NA)
#             }
#           } %>% reduce(full_join, by = 'date') %>% 
#           arrange(date) %>% select(-date) %>% 
#           setNames(paste0('ens0', 1:5))
#         return(unbiased_ssp245) 
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start
# 
# ## checkpoint
# save(precip_ssp245_sdm, file = '_data/GFDL/biascorr/files/precip_ssp245_sdm.Rdata')


#### bias corr: SSP5-8.5 precipitation ############################################################
cat('...SSP 5-8.5\n')

# start <- Sys.time()
# pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
# precip_ssp585_sdm <- 
#   foreach (dec = length(decades), .combine = 'rbind') %:% 
#   foreach (
#     i = 1:ncell(grid_ca),
#     .export = 
#       c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt', 
#         'calculate.rp100', 'calculate.aic', 'add_index'),
#     .packages = 
#       c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
#     .options.snow = opts) %dopar% {
#       if (i %in% index_ca) {
#         ## load MERRA mixture model parameters
#         values_merra <- precip_merra_param[[i]][[1]]
#         par_merra <- precip_merra_param[[i]][[2]]
#         rd_merra <- sum(values_merra > nonzero)/length(values_merra)
#         
#         ## calculate mixture model parameters for GFDL historic
#         values_hist <- matrix(nrow = length(values_merra), ncol = 5)
#         par_hist <- matrix(nrow = 6, ncol = 5)
#         rd_hist <- apply(values_hist, 2, function(x) sum(x > nonzero))/nrow(values_hist)
# 
#         unbiased_ssp585 <- 
#           foreach (ens = 1:5) %do% {
#             ## calculate mixture model parameters for GFDL future
#             df_ssp585 <- 
#               precip_ssp585[[ens]] %>% 
#               raster::extract(i) %>% c %>%
#               .[1:length(ts_future)] %>% 
#               data.frame(date = ts_future, values = .) %>% 
#               filter(month(date) %in% c(10:12,1:3)) %>% 
#               filter(year(date) %in% tridecades[[dec]])
#             par_ssp585 <- tryCatch(
#               optimize.aic(df_ssp585$values, precip.domain)$par,
#               error = function(e) rep(0,6))
#             
#             if (sum(par_merra) > 0 & sum(par_hist[,ens]) > 0 & sum(par_ssp245) > 0) {
#               ## apply SDM to bias-correct GFDL future
#               temp <- sdm2(
#                 par_future = par_ssp585, 
#                 data = df_ssp585 %>% filter(year(date) %in% decades[[dec]]),
#                 par_merra = par_merra, values_merra = values_merra, 
#                 par_hist = par_hist[,ens], values_hist = values_hist[,ens],
#                 domain = ivt.domain) %>%
#                 select(date, unbiased)
#           
#               ## apply rainy-day correction
#               rd <- sum(temp$unbiased > nonzero)/nrow(temp)
#               rd_scaled <- rd * rd_merra / rd_hist
#               temp %>% 
#                 arrange(desc(unbiased)) %>% 
#                 mutate(
#                   position = add_index(nrow(.)),
#                   unbiased_rd = ifelse(position > rd_scaled, 0, unbiased)) %>% 
#                 select(date, unbiased_rd)
#             } else {
#               df_ssp585 %>% select(date) %>% mutate(empty = NA)
#             }
#           } %>% reduce(full_join, by = 'date') %>% 
#           arrange(date) %>% select(-date) %>% 
#           setNames(paste0('ens0', 1:5))
#         return(unbiased_ssp585) 
#       } else NULL
#     }
# stopCluster(cl)
# Sys.time() - start
# 
# ## checkpoint
# save(precip_ssp585_sdm, file = '_data/GFDL/biascorr/files/precip_ssp585_sdm.Rdata')


###################################################################################################
cat('done!\n\n')


