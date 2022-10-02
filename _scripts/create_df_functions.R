
###################################################################################################

#' add_counter2 <- function(ar) {
#'   #' Creates a vector of ordered indices to distinguish sequences of true values.
#'   #' @param ar vector of true/false values
#'   #' @return vector of numbers (where ar=true) and NA (where ar=false)
#' 
#'   counter <- 0
#'   foreach (i = 1:length(ar), .combine = 'c') %do% {
#'     if (ar[i]) {
#'       if (i==1 | !ar[i-1]) counter <- counter+1
#'       counter
#'     } else NA
#'   }
#' }

add_counter <- function(ar) {
  #' Creates a vector of ordered indices to distinguish sequences of true values.
  #' @param ar vector of true/false values
  #' @return vector of numbers (where ar=true) and NA (where ar=false)
  
  ## identify when the vector changes from true to false
  index <- which(ar)
  change <- rep(0, length(ar))
  change[index][which(c(0,diff(index))>1)] <- 1
  
  ## convert changepoints to event indices
  count <- cumsum(change)+1
  count[!ar] <- NA
  return(count)
  }

# all(setNA(add_counter(ar),0) == setNA(add_counter2(ar),0))


###################################################################################################

create_catalog <- function(df, name, cat = TRUE, interval = 3) {
  #' Creates a catalog of events.
  #' @param df dataframe with column ts (datetime), ar (logical), and ivt (double)
  #' @param name name of the column to base events on 
  #' @param cat logical indicating whether to calculate AR intensity categories
  #' @param interval temporal resolution of the dataframe, in hours
  #' @return dataframe of AR events with columns ar (integer), start (datetime), end (datetime), 
  #'   ivt_max (double), duration (double), and cat (integer)
  
  catalog <- df %>% 
    mutate(count = add_counter(get(name))) %>% 
    filter(get(name)) %>% 
    group_by(count) %>% 
    summarize(
      start = min(ts), 
      end = max(ts),
      maxivt = max(ivt), 
      duration = length(ivt)*interval) 
  if (cat) {
    return(
      catalog %>% 
        mutate(cat = map2_dbl(.x = maxivt, .y = duration, .f = ~assign_AR_cat(.x, .y))))
  } else return(catalog)
}


###################################################################################################

attach_impacts <- function(daily, i) {
  #' Attaches all impacts data sources to daily dataframe.
  #' @param daily daily dataframe
  #' @param i index of grid cell under consideration
  #' @return daily dataframe with additional columns
  
  ## calculate start-end daily intervals of disaster declarations
  intervals <-
    disasters %>%
    filter(fips == 6e3+grid_county[i]) %>%
    map2(.x = ymd(.$start_day), .y = ymd(.$end_day), .f = ~interval(.x, .y)) %>%
    do.call('c', .)
  
  ## attach information to df_merra
  daily <- daily %>%
    ## disaster declarations
    left_join(
      data.frame(
        date = dates_merra,
        disaster = map_dbl(
          .x = dates_merra,
          .f = ~any(.x %within% intervals)) %>% as.logical),
      by = 'date') %>%
    ## public assistance
    left_join(
      pa_grid %>% filter(cell == i) %>% select(-cell),
      by = 'date') %>%
    mutate(pa_amt = case_when(date >= min(pa$start_day) ~ setNA(pa_amt,0))) %>%
    ## NCEI
    left_join(
      ncei_grid %>% filter(cell == i) %>% select(-cell),
      by = 'date') %>%
    mutate(ncei_damage = setNA(ncei_damage,0)) %>%
    ## claims
    left_join(
      claims %>% filter(cell == i) %>% select(-cell),
      by = 'date') %>%
    mutate(
      claims_num = setNA(claims_num,0),
      claims_value = setNA(claims_value,0),
      claims_coverage = setNA(claims_coverage,0),
      claims_dr = setNA(claims_dr,0)) %>%
    ## policies
    mutate(wy = wateryear(date)) %>%
    left_join(
      policies %>% filter(cell == i) %>% select(-cell),
      by = c('wy' = 'coverage_wy')) %>%
    mutate(
      policies_num =
        case_when(wy >= min(policies$coverage_wy) ~ setNA(policies_num,0)),
      policies_value_adj =
        case_when(wy >= min(policies$coverage_wy) ~ setNA(policies_value_adj,0))) %>%
    rename(policies_value = policies_value_adj) %>% 
    select(-wy) %>%
    ## population
    mutate(pop = case_when(
      year(date) <= 2000 ~ pop2000_raster[i],
      year(date) <= 2010 ~ pop2010_raster[i],
      year(date) <= 2015 ~ pop2015_raster[i],
      TRUE ~ pop2020_raster[i])) %>%
    ## WWA
    left_join(
      wwa_grid %>% filter(cell == i) %>% select(-cell),
      by = 'date') %>%
    mutate(
      wwa_phenom = case_when(date >= min(wwa_grid$date) ~ setNA(wwa_phenom, 'X')),
      wwa_vtec = case_when(date >= min(wwa_grid$date) ~ setNA(wwa_vtec, 'none'))) %>%
    arrange(date)
  
  ## return dataframe
  return(daily)
}
