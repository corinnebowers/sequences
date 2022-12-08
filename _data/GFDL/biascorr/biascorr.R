
###################################################################################################
## biascorr.R: This script creates generates .Rdata files of bias-corrected IVT and precipitation
##   (historic, SSP 2-4.5, and SSP 5-8.5) for all raster grid cells in California.
## Created by Corinne Bowers 9/20/22
## Last updated 11/3/22
###################################################################################################

#### setup ########################################################################################
cat('setting up...\n')

## set working directory
# setwd('D:/2-sequences/')
setwd('/scratch/users/cbowers/sequences/')

## load setup file
source('_data/GFDL/biascorr/biascorr_setup.R')


#### historic IVT #################################################################################
cat('bias correcting IVT...\n')
cat('...historic\n')

start <- Sys.time()
pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)
cl <- parallel::makeCluster(cores)
registerDoSNOW(cl)
ivt_hist_sdm <-
  foreach (
    i = 1:ncell(grid_ca),
    .export =
      c('mixture', 'qmixture', 'dmixture', 'pmixture', 'wt', 'calculate.rp100', 'calculate.aic'),
    .packages =
      c('POT', 'fitdistrplus', 'foreach', 'pracma', 'lubridate', 'tidyverse'),
    .options.snow = opts) %dopar% {
      if (i %in% index_ca) {
        ## load MERRA mixture model parameters
        values_merra <- ivt_merra_param[[i]][[1]]
        par_merra <- ivt_merra_param[[i]][[2]]
        ks_merra <- ivt_merra_param[[i]][[3]]

        ## load GFDL historic mixture model parameters
        values_hist <- ivt_hist_param[[i]][[1]]
        par_hist <- ivt_hist_param[[i]][[2]]
        ks_hist <- ivt_hist_param[[i]][[3]]

        ## apply SDM to bias-correct GFDL historic
        unbiased <-
          foreach (ens = 1:5) %do% {
            df_hist <- ivt_gfdl[[ens]] %>%
              raster::extract(i) %>% c %>%
              data.frame(date = ts_hist, values = .) %>%
              mutate(values = cbind(values,ivt.domain[1]) %>% apply(1,max))
            if (sum(par_merra) > 0 & sum(par_hist[,ens]) > 0 &
                ks_merra > 0.05 & ks_hist[ens] > 0.05) {
              sdm2(
                par_future = par_hist[,ens], data = df_hist,
                par_merra = par_merra, values_merra = values_merra,
                par_hist = par_hist[,ens], values_hist = values_hist[,ens],
                domain = ivt.domain) %>%
                select(date, unbiased)
            } else {
              df_hist %>% select(date) %>% mutate(unbiased = NA)
            }
          } %>% reduce(full_join, by = 'date') %>%
          arrange(date) %>% select(-date) %>%
          setNames(paste0('ens0', 1:5))
        return(unbiased)
      } else NULL
    }
# stopCluster(cl)
Sys.time() - start

## checkpoint
save(ivt_hist_sdm, file = '_data/GFDL/biascorr/files/ivt_hist_sdm.Rdata')


#### SSP2-4.5 IVT #################################################################################
cat('...SSP 2-4.5\n')

ivt_ssp245_sdm <- list()
pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)

start <- Sys.time()
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
for (i in 1:ncell(grid_ca)) {
  if (i %in% index_ca) {
    ## load MERRA mixture model parameters
    values_merra <- ivt_merra_param[[i]][[1]]
    par_merra <- ivt_merra_param[[i]][[2]]
    ks_merra <- ivt_merra_param[[i]][[3]]
    
    ## load GFDL historic mixture model parameters
    values_hist <- ivt_hist_param[[i]][[1]]
    par_hist <- ivt_hist_param[[i]][[2]]
    ks_hist <- ivt_hist_param[[i]][[3]]
    
    ## bias-correct GFDL future 
    ivt_ssp245_sdm[[i]] <- 
      foreach(
        ens = 1:5,
        .export =
          c('mixture', 'qmixture', 'dmixture', 'pmixture', 
            'wt', 'calculate.rp100', 'calculate.aic'),
        .packages = c('POT', 'fitdistrplus', 'pracma', 'lubridate', 'tidyverse')) %:% 
      foreach(dec = 1:length(decades), .combine = 'rbind') %dopar% {
        ## calculate mixture model parameters
        df_ssp245 <- ivt_ssp245[[ens]] %>%
          raster::extract(i) %>% c %>%
          data.frame(date = ts_future, values = .) %>%
          mutate(values = cbind(values,ivt.domain[1]) %>% apply(1,max))
        values_ssp245 <- df_ssp245 %>%
          filter(month(date) %in% c(10:12,1:3) & year(date) %in% tridecades[[dec]]) %>% 
          pull(values)
        par_ssp245 <- tryCatch(
          optimize.aic(values_ssp245, ivt.domain),
          error = function(e) rep(0,6))
        
        ## apply SDM to bias-correct one decade+ensemble
        if (sum(par_merra) > 0 & sum(par_hist[,ens]) > 0 & sum(par_ssp245[-7]) > 0 & 
            ks_merra > 0.05 & ks_hist[ens] > 0.05 & par_ssp245[7] > 0.05) {
          sdm2(
            par_future = par_ssp245[-7],
            data = df_ssp245 %>% filter(year(date) %in% decades[[dec]]),
            par_merra = par_merra, values_merra = values_merra,
            par_hist = par_hist[,ens], values_hist = values_hist[,ens],
            domain = ivt.domain) %>%
            select(date, unbiased)
        } else {
          df_ssp245 %>% 
            filter(year(date) %in% decades[[dec]]) %>% 
            select(date) %>% mutate(unbiased = NA)
        }
      } %>% reduce(full_join, by = 'date') %>%
      arrange(date) %>% select(-date) %>%
      setNames(paste0('ens0', 1:5))
  } else ivt_ssp245_sdm[[i]] <- NULL
  
  setTxtProgressBar(pb,i)
}
# stopCluster(cl)
Sys.time() - start

## checkpoint
save(ivt_ssp245_sdm, file = '_data/GFDL/biascorr/files/ivt_ssp245_sdm.Rdata')


#### SSP5-8.5 IVT #################################################################################
cat('...SSP 5-8.5\n')

ivt_ssp585_sdm <- list()
pb <- txtProgressBar(min = 0, max = ncell(grid_ca), style = 3)

start <- Sys.time()
# cl <- parallel::makeCluster(cores)
# registerDoSNOW(cl)
for (i in 1:ncell(grid_ca)) {
  if (i %in% index_ca) {
    ## load MERRA mixture model parameters
    values_merra <- ivt_merra_param[[i]][[1]]
    par_merra <- ivt_merra_param[[i]][[2]]
    ks_merra <- ivt_merra_param[[i]][[3]]
    
    ## load GFDL historic mixture model parameters
    values_hist <- ivt_hist_param[[i]][[1]]
    par_hist <- ivt_hist_param[[i]][[2]]
    ks_hist <- ivt_hist_param[[i]][[3]]
    
    ## bias-correct GFDL future 
    ivt_ssp585_sdm[[i]] <- 
      foreach(
        ens = 1:5,
        .export =
          c('mixture', 'qmixture', 'dmixture', 'pmixture', 
            'wt', 'calculate.rp100', 'calculate.aic'),
        .packages = c('POT', 'fitdistrplus', 'pracma', 'lubridate', 'tidyverse')) %:% 
      foreach(dec = 1:length(decades), .combine = 'rbind') %dopar% {
        ## calculate mixture model parameters
        df_ssp585 <- ivt_ssp585[[ens]] %>%
          raster::extract(i) %>% c %>%
          data.frame(date = ts_future, values = .) %>%
          mutate(values = cbind(values,ivt.domain[1]) %>% apply(1,max))
        values_ssp585 <- df_ssp585 %>%
          filter(month(date) %in% c(10:12,1:3) & year(date) %in% tridecades[[dec]]) %>% 
          pull(values)
        par_ssp585 <- tryCatch(
          optimize.aic(values_ssp585, ivt.domain),
          error = function(e) rep(0,6))
        
        ## apply SDM to bias-correct one decade+ensemble
        if (sum(par_merra) > 0 & sum(par_hist[,ens]) > 0 & sum(par_ssp585[-7]) > 0 & 
            ks_merra > 0.05 & ks_hist[ens] > 0.05 & par_ssp585[7] > 0.05) {
          sdm2(
            par_future = par_ssp585[-7],
            data = df_ssp585 %>% filter(year(date) %in% decades[[dec]]),
            par_merra = par_merra, values_merra = values_merra,
            par_hist = par_hist[,ens], values_hist = values_hist[,ens],
            domain = ivt.domain) %>%
            select(date, unbiased)
        } else {
          df_ssp585 %>%
            filter(year(date) %in% decades[[dec]]) %>% 
            select(date) %>% mutate(unbiased = NA)
        }
      } %>% reduce(full_join, by = 'date') %>%
      arrange(date) %>% select(-date) %>%
      setNames(paste0('ens0', 1:5))
  } else ivt_ssp585_sdm[[i]] <- NULL
  
  setTxtProgressBar(pb,i)
}
stopCluster(cl)
Sys.time() - start

## checkpoint
save(ivt_ssp585_sdm, file = '_data/GFDL/biascorr/files/ivt_ssp585_sdm.Rdata')


#### historic precipitation #######################################################################
# cat('bias correcting precipitation...\n')
# cat('...historic\n')

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


#### SSP2-4.5 precipitation #######################################################################
# cat('...SSP 2-4.5\n')

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


#### SSP5-8.5 precipitation #######################################################################
# cat('...SSP 5-8.5\n')

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


