
###################################################################################################
## biascorrection_functions.R: This script contains functions to support the file biascorrection.Rmd.
## Created by Corinne Bowers 9/20/22
## Last updated 11/3/22
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
    cdf <- if(sum(pdf)==0) c(0,pdf) else cumsum(c(0,pdf)/sum(pdf))
    pdf <- diff(cdf)
    cdf <- cdf[-length(cdf)]/max(cdf[-length(cdf)])
    
    ## return values
    if (fun == 'd') {
      suppressWarnings(interp1(x = q, y = pdf, xi = val))
    } else if (fun == 'p') {
      suppressWarnings(interp1(x = q, y = cdf, xi = val))
    } else if (fun == 'q') {
      suppressWarnings(interp1(x = cdf, y = q, xi = val))
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


#### return period function #######################################################################

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


#### distribution optimization functions ##########################################################

calculate.aic <- 
  function(par, val, rp100_ci, domain) {
    # print(par) #to troubleshoot
    
    ## reject invalid parameter choices
    if (par[1] < 0 | par[2] < 0 | par[3] < 0) return(1e10)
    
    ## calculate AIC
    loglik <- sum(log(dmixture(par, val, domain)))
    if (is.infinite(loglik)) return(1e10)
    aic <- -2*loglik + 2*length(par)
    
    ## penalize AIC based on 100-year RP criterion
    rp100 <- qmixture(par, p = 1 - (0.01/182), domain)
    penalty <- (rp100-mean(rp100_ci))^2
    if (rp100 < rp100_ci[1] | rp100 > rp100_ci[2]) return(aic + penalty) else return(aic)
  }

optimize.aic <- function(values, domain, par = NA) {
  ## filter bad domain input
  if (domain[1] == 0) stop('Domain minimum must be greater than zero.')
  if (length(seq(domain[1], domain[2], domain[3])) < 10) {
    stop('Please define at least 10 values for the interpolation domain.')
  }
  
  ## create initial parameter guesses, if not provided
  if (any(is.na(par)) | sum(par)==0) {
    par <- rep(NA, 6)
    par[1:2] <- fitdist(values, 'gamma')$estimate
    par[3:4] <- fitgpd(values, threshold = quantile(values, 0.95))$fitted.values
    par[5] <- quantile(values, 0.95)
    par[6] <- par[5]/2
  }

  ## calculate GPD-only 100-year IVT value
  rp100_ci <- calculate.rp100(par, values, threshold = par[5], nyr = 30)
  
  ## run parameter optimization
  bestfit <- optim(
    par = par, fn = calculate.aic, 
    val = values, rp100_ci = rp100_ci, domain = domain,
    method = 'BFGS',
    control = list(maxit = 1e6))
  
  ## check convergence
  if (bestfit$convergence > 0) stop('Optimization did not converge.')  
  
  ## check K-S statistic and return best-fit parameters
  ks <- ks.test(values, qmixture(bestfit$par, p = add_index(length(values)), domain))$p.value
  if (ks < 0.05) {
    bestfit2 <- optim(
      par = par, fn = calculate.aic, 
      val = values, rp100_ci = rp100_ci, domain = domain,
      control = list(maxit = 1e6))
    ks2 <- ks.test(values, qmixture(bestfit2$par, p = add_index(length(values)), domain))$p.value
    if (ks2 > ks) return(c(bestfit2$par, ks2))
  } 
  return(c(bestfit$par, ks))
}


#### SDM functions ################################################################################

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

rainy <- function(vals, nonzero) sum(vals > nonzero)/length(vals)


###################################################################################################
