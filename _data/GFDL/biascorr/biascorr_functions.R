
###################################################################################################
## biascorrection_functions.R: This script contains functions to support the file biascorrection.Rmd.
## Created by Corinne Bowers 9/20/22
## Last updated 10/19/22
###################################################################################################

#### mixture model functions ######################################################################

wt <- function(x,m,t) { 1/2 + 1/pi * atan((x-m)/t) }

mixture <- 
  function(par, qmin, qmax, dq, fun, val) {
    ## define distribution parameters
    gamma.shape <- par[1]
    gamma.rate <- par[2]
    gpd.scale <- par[3]
    gpd.shape <- par[4]
    m <- par[5]
    t <- par[6]
    
    ## generate interpolation lookup
    q <- seq(qmin, qmax, dq)
    pdf <- 
      (1-wt(q,m,t)) * dgamma(q, shape = gamma.shape, rate = gamma.rate) +
      wt(q,m,t) * dgpd(q, scale = gpd.scale, shape = gpd.shape)
    pdf[is.infinite(pdf)] <- 1e9
    cdf <- cumsum(c(0,pdf)/sum(pdf))
    pdf <- diff(cdf)
    
    ## return values
    if (fun == 'd') {
      suppressWarnings(interp1(x = q, y = pdf, xi = val))
    } else if (fun == 'p') {
      suppressWarnings(interp1(x = q, y = cdf[-length(cdf)], xi = val))
    } else if (fun == 'q') {
      suppressWarnings(interp1(x = cdf[-length(cdf)], y = q, xi = val))
    } else {
      stop('Not a valid function label.')
    }
  }

dmixture <- 
  function(par, x, domain) {
    mixture(par, domain[1], domain[2], domain[3], 'd', x)
  }
pmixture <- 
  function(par, q, domain) {
    mixture(par, domain[1], domain[2], domain[3], 'p', q)
  }
qmixture <- 
  function(par, p, domain) {
    mixture(par, domain[1], domain[2], domain[3], 'q', p)
  }


#### distribution optimization functions ##########################################################

calculate.rp100 <- 
  function(par, val, threshold, nyr) {
    ## define distribution parameters for GPD
    gpd.loc <- threshold
    gpd.scale <- par[3]
    gpd.shape <- par[4]
    m <- 100*nyr
    pr <- sum(val > gpd.loc)/length(val)
    
    ## calculate 95% confidence interval
    x_m <- gpd.loc + (gpd.scale/gpd.shape) * ((m*pr)^gpd.shape - 1)
    sigma <- ((x_m - gpd.loc)*gpd.shape) / ((m*pr)^gpd.shape - 1)
    rp100_ci <- qnorm(c(0.025,0.975), mean = x_m, sd = sigma)
    return(rp100_ci)  
  }

calculate.aic <- 
  function(par, val, rp100_ci, domain) {
    # print(par) #to troubleshoot
    
    ## reject invalid parameter choices
    if (par[1] < 0 | par[2] < 0) return(1e10)
    
    ## calculate AIC
    loglik <- sum(log(dmixture(par, val, domain)))
    aic <- -2*loglik + 2*length(par)
    
    ## penalize AIC based on 100-year RP criterion
    rp100 <- qmixture(par, p = 1 - (0.01/182), domain)
    if (rp100 < rp100_ci[1] | rp100 > rp100_ci[2]) return(aic + 1e6) else return(aic)
  }

optimize.aic <- function(values, domain) {
  ## filter bad domain input
  if (domain[1] == 0) stop('Domain minimum must be greater than zero.')
  if (length(seq(domain[1], domain[2], domain[3])) < 10) {
    stop('Please define at least 10 values for the interpolation domain.')
  }
  
  ## create initial parameter guesses
  par <- rep(NA, 6)
  par[1:2] <- fitdist(values, 'gamma')$estimate
  par[3:4] <- fitgpd(values, threshold = quantile(values, 0.95))$fitted.values
  par[5] <- quantile(values, 0.95)
  par[6] <- par[5]/2
  
  ## calculate GPD-only 100-year IVT value
  rp100_ci <- calculate.rp100(par, values, threshold = par[5], nyr = 30)
  
  ## run parameter optimization
  bestfit <- optim(
    par = par, 
    fn = calculate.aic, 
    val = values, rp100_ci = rp100_ci, domain = domain,
    control = list(maxit = 1e6))
  
  ## check convergence
  if (bestfit$convergence > 0) stop('Optimization did not converge.')  
  
  ## check K-S statistic
  ks <- ks.test(values, qmixture(bestfit$par, p = add_index(length(values)), domain))$p.value
  if (ks < 0.05) warning('Mixture model distribution may not be well fit.')
  
  ## return best-fit parameters
  return(bestfit)
}

## plot mixture models
plot.mixture <- function(par, values, qq = TRUE, pp = TRUE) {
  df <- data.frame(q = sort(values), p = add_index(length(values))) %>% 
    mutate(
      fitted.q = qmixture(par, p, domain),
      fitted.p = pmixture(par, q, domain))
  val.max <- df %>% select(q, fitted.q) %>% apply(2, max) %>% max
  ks.pass <- ks.test(df$q, df$fitted.q)$p.value > 0.05
  ks.col <- ifelse(ks.pass, 'black', 'grey70')
  
  g <- ggplot(df) + 
    geom_parity() + 
    labs(x = 'Observed', y = 'Fitted') +
    theme(panel.grid.major = element_line())
  if (pp) {
    g1 <- g + 
      geom_point(aes(x = p, y = fitted.p), color = ks.col) + 
      scale_x_origin() + scale_y_origin() + coord_fixed() + 
      ggtitle('Probability Plot')
    print(g1)
  }
  if (qq) {
    g2 <- g + 
      geom_point(aes(x = q, y = fitted.q), color = ks.col) + 
      scale_x_origin(breaks = 250*(-5:10)) + scale_y_origin(breaks = 250*(-5:10)) + 
      coord_fixed(xlim = c(0,val.max), ylim = c(0,val.max)) + 
      ggtitle('Quantile Plot')
    print(g2)
  }
}

#### SDM functions ################################################################################

# sdm <- 
#   function(
#     par_future, data, par_merra, values_merra, par_hist, values_hist, domain) {
#     # note: data is a dataframe with the columns date, values, paramid
#     
#     ## calculate return intervals for observed (MERRA) data
#     ri_merra <- 1 / (1 - pmixture(par_merra, values_merra, domain))
#     
#     ## calculate return intervals for simulated historic data
#     ri_hist <- 1 / (1 - pmixture(par_hist, values_hist, domain))
#     
#     ## calculate return intervals for simulated future data & scale to match 
#     datadec <- 
#       foreach (i = 1:nrow(par_future), .combine = 'rbind') %do% {
#         data %>% 
#           filter(paramid == i) %>% 
#           arrange(values) %>% 
#           mutate(index = 1:nrow(.)) %>% 
#           mutate(
#             p = pmixture(par = par_future[i,], q = values, domain),
#             ri = 1 / (1-p)) %>% 
#           mutate(
#             ri_merra_interp = interp1(
#               x = (1:length(values_merra))/length(values_merra)*nrow(.), 
#               y = sort(ri_merra), 
#               xi = index),
#             ri_hist_interp = interp1(
#               x = (1:length(values_hist))/length(values_hist)*nrow(.),
#               y = sort(ri_hist), 
#               xi = index))
#       }
#     
#     ## calculate bias-corrected data
#     datadec <- datadec %>% 
#       mutate(sf = values/qmixture(par = par_hist, p = p, domain)) %>% 
#       mutate(
#         ri_scaled = (ri*ri_merra_interp/ri_hist_interp) %>% cbind(1) %>% apply(1,max),
#         p_scaled = 1 - (1/ri_scaled)) %>% 
#       mutate(unbiased = qmixture(par = par_merra, p = p_scaled, domain) * sf)
#     
#     ## re-insert bias-corrected data into correct locations
#     datadec <- datadec %>% 
#       group_by(paramid) %>% 
#       mutate(unbiased = sort(unbiased)) %>% 
#       ungroup
#     return(data %>% left_join(datadec %>% select(date, unbiased), by = 'date'))
#   }

sdm2 <- 
  function(
    par_future, data, par_merra, values_merra, par_hist, values_hist, domain) {
    # note: data is a dataframe with the columns date, values
    
    ## calculate return intervals for observed (MERRA) data
    ri_merra <- 1 / (1 - pmixture(par_merra, values_merra, domain))
    
    ## calculate return intervals for simulated historic data
    ri_hist <- 1 / (1 - pmixture(par_hist, values_hist, domain))
    
    ## calculate return intervals for simulated future data & scale to match 
    x_merra <- seq(0, length(values_merra), length.out = length(values_merra))
    x_hist <- seq(0, length(values_hist), length.out = length(values_hist))
    df <- data %>% 
      arrange(values) %>% 
      mutate(
        index = 1:nrow(.),
        p = pmixture(par = par_future, q = values, domain),
        ri = 1 / (1-p)) %>% 
      mutate(
        ri_merra_interp = interp1(
          x = x_merra/(length(values_merra))*nrow(.), 
          y = sort(ri_merra), 
          xi = index),
        ri_hist_interp = interp1(
          x = x_hist/(length(values_hist))*nrow(.),
          y = sort(ri_hist), 
          xi = index)) %>%
      mutate(sf = values/qmixture(par = par_hist, p = p, domain)) %>% 
      mutate(
        ri_scaled = (ri*ri_merra_interp/ri_hist_interp) %>% cbind(1) %>% apply(1,max),
        p_scaled = 1 - (1/ri_scaled)) %>% 
      mutate(unbiased = sort(qmixture(par = par_merra, p = p_scaled, domain) * sf)) 
      
    return(data %>% left_join(df %>% select(date, unbiased), by = 'date'))
  }

rainy <- function(data) {
  data %>% 
    mutate(rainy = ifelse(values > 0.1, TRUE, FALSE)) %>% 
    count(rainy) %>% 
    mutate(pct = prop.table(n)) %>% 
    filter(rainy) %>% pull(pct)
}


#### wrapper function #############################################################################

biascorrect <- function(merra, gfdl_hist, gfdl_future, type) {
  #' Outputs a bias-corrected timeseries dataframe of the same length as gfdl_future.
  #' 
  #' @param merra dataframe with columns "date" and "values"
  #' @param gfdl_hist dataframe with columns "date" and "values"
  #' @param gfdl_future dataframe with columns "date" and "values"
  #' @param type string specifying whether to bias-correct IVT or precip
  #' 
  #' @return 
  #' 

  start <- Sys.time()
  
  ## define constants based on analysis type
  cat('setting up variables...\n')
  if (type == 'IVT' | type == 'ivt') {
    precip <- FALSE
    domain <- c(1e-1, 5e3, 1e-1)
  } else if (type == 'precip' | type == 'precipitation') {
    precip <- TRUE
    domain <- c(1e-1, 500, 1e-2)
  } else stop('Not a valid keyword for "type".')

  ## format provided datasets
  data_merra <- merra %>% 
    filter(month(date) %in% c(10:12,1:3)) %>% 
    filter(year(date) %in% 1981:2010)
  values_merra <- data_merra %>% pull(values)
  
  data_hist <- gfdl_hist %>% 
    filter(month(date) %in% c(10:12,1:3)) %>% 
    filter(year(date) %in% 1981:2010)
  values_hist <- data_hist %>% pull(values)
  
  data_future <- gfdl_future %>% 
    filter(month(date) %in% c(10:12,1:3)) %>% 
    mutate(paramid = case_when(
      year(date) %in% 2021:2030 ~ 1,
      year(date) %in% 2031:2040 ~ 2,
      year(date) %in% 2041:2050 ~ 3,
      year(date) %in% 2051:2060 ~ 4,
      year(date) %in% 2061:2070 ~ 5,
      year(date) %in% 2071:2080 ~ 6,
      year(date) %in% 2081:2090 ~ 7))
  
  ## get distribution parameters
  cat('fitting mixture model distribution parameters...\n')
  cat('...observed\n')
  # start <- Sys.time()
  par_merra <- optimize.aic(values_merra, domain)$par
  Sys.time() - start

  cat('...simulated historic\n')
  # start <- Sys.time()
  par_hist <- optimize.aic(values_hist, domain)$par
  Sys.time() - start
  
  cat('...simulated future\n')
  # start <- Sys.time()
  decades <- 
    list(2011:2040, 2021:2050, 2031:2060, 2041:2070, 2051:2080, 2061:2090, 2071:2100)
  par_future <- 
    foreach(i = 1:length(decades), .combine = 'rbind') %do% { 
      data_future %>% 
        filter(year(date) %in% decades[[i]]) %>% 
        pull(values) %>% 
        optimize.aic(., domain) %>% .$par
    }
  Sys.time() - start
  
  ## perform SDM
  cat('performing SDM bias correction...\n')
  # start <- Sys.time()
  data_bc <- 
    sdm(
      par_future = par_future, data = data_future,
      par_merra = par_merra, values_merra = values_merra, 
      par_hist = par_hist, values_hist = values_hist,
      domain = domain) %>% 
    arrange(date)
  Sys.time() - start
  
  ## for precip only: perform dry-day frequency scaling
  if (precip) {
    cat('performing dry-day frequency correction...\n')
    # start <- Sys.time()
    rd_merra <- rainy(data_merra)
    rd_hist <- rainy(data_hist)
    data_bcrd <- 
      foreach (i = 1:length(decades), .combine = 'rbind') %do% {
        ## calculate rainy day correction factor
        rd_future <- data_future %>%
          filter(paramid == i) %>% 
          rainy(.)
        rd_scaled <- rd_future * rd_merra / rd_hist
        
        ## set drizzle events to zero
        data_bc %>% 
          filter(paramid==i) %>% 
          arrange(unbiased) %>% 
          mutate(rank = 1 - add_index(nrow(.))) %>% 
          mutate(unbiased = ifelse(rank > rd_scaled, 0, unbiased))
      }
    data_bc <- data_bc %>% select(date, values) %>% 
      left_join(data_bcrd %>% select(date, unbiased), by = 'date')
    Sys.time() - start
  } else {
    data_bc <- data_bc %>% select(date, values, unbiased)
  }
  
  ## return bias-corrected timeseries
  cat('done!\n\n')
  return(data_bc)
}

###################################################################################################
