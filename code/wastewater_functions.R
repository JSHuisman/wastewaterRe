###########################################################
# Function to estimate Re from WW
# Author: JS Huisman
###########################################################
# These are mostly wrapper functions around the 
# deconvolution and Re estimation functions from 
# https://github.com/covid-19-Re/shiny-dailyRe

###########################################################
# Delay/Shedding Load Distributions #####

# find gamma parameters from mean/sd of distribution
getGammaParams <- function(meanParam, sdParam){
  shapeParam <- meanParam^2 / (sdParam^2)
  scaleParam <- (sdParam^2) / meanParam
  return(list(shape = shapeParam, scale = scaleParam))
}

#find mean/sd from gamma shape/scale parameters
getInvGammaParams <- function(shapeParam, scaleParam){
  meanParam <- scaleParam * shapeParam
  sdParam <- sqrt(scaleParam^2 * shapeParam)
  return(list(mean = meanParam, sd = sdParam))
}

# specific distributions we use in the paper
getCountParams <- function(obs_type){
  switch(obs_type,
         incubation = getGammaParams(5.3, 3.2),
         zero = list(shape = 0, scale = 0),
         #confirmed = getGammaParams(5.5, 3.8),
         confirmed_zh = getGammaParams(2.83, 2.96),
         confirmed_cali = getGammaParams(4.51, 3.16),
         death = getGammaParams(15.0, 6.9),
         han = getGammaParams(4.7, 1.7),
         wolfel = getGammaParams(8.6, 0.9),
         benefield = list(shape = 0.929639, scale = 7.241397))
}

###########################################################
# Deconvolve #####

# To achieve an input format required by the 
# get_infection_incidence_by_deconvolution function
addUselessColumns <- function(df, inc_var = 'n1'){
  if (!'region' %in% colnames(df)){
    df <- df %>%
      mutate(region = 'CHE')
  }
  
  observation_df <- df %>%
    dplyr::select(date, region, value = inc_var) %>%
    mutate(data_type = inc_var,
           source = 'ETH',
           variable = 'incidence',
           country = 'Switzerland',
           date_type = 'report',
           local_infection = TRUE)
  
  observation_df <- observation_df %>%
    filter(!is.na(value))
  
  return(observation_df)
}

# wrapper around the actual deconvolution function
deconvolveIncidence <- function(df, incidence_var = 'n1',
                                IncubationParams, OnsetToCountParams,
                                smooth_param = FALSE, n_boot = 50){
  infection_df <- addUselessColumns(df, inc_var = incidence_var)
  
  constant_delay_distributions <- list("Simulated" = get_vector_constant_waiting_time_distr(
    IncubationParams$shape, IncubationParams$scale,
    OnsetToCountParams$shape, OnsetToCountParams$scale),
    "Symptoms" = get_vector_constant_waiting_time_distr(
      IncubationParams$shape, IncubationParams$scale,
      0, 0))
  
  estimatedInfections <- get_infection_incidence_by_deconvolution(
    infection_df,
    is_local_cases = T,
    constant_delay_distribution = constant_delay_distributions[['Simulated']],
    constant_delay_distribution_incubation = constant_delay_distributions[["Symptoms"]],
    max_iterations = 100,
    smooth_incidence = smooth_param,
    empirical_delays = tibble(),
    n_bootstrap = n_boot,
    verbose = FALSE)
  
  return(estimatedInfections)
}

###########################################################
# Estimate Re ####

# wrapper around the Re estimation function
getReBootstrap <- function(deconvoluted_data){
  
  all_delays <- lapply(unique(deconvoluted_data$data_type), function(x){ c(Cori = 0)})
  names(all_delays) <- unique(deconvoluted_data$data_type)
  
  truncations <- list(left = c(Cori = 5),
                      right = c(Cori = 0))
  
  rawReEstimates <- suppressWarnings(doAllReEstimations(
    deconvoluted_data,
    slidingWindow = 3,
    methods = c("Cori"),
    variationTypes = c("slidingWindow"),
    all_delays,
    truncations,
    interval_ends = list()) )
  
  cleanEstimates <- cleanCountryReEstimate(rawReEstimates, method = 'bootstrap',
                                           rename_types = F, alpha=0.95)
  
  return(cleanEstimates)
}

###########################################################
compareTraces <- function(Re_i, Re_j){
  compare_df = Re_i %>%
    left_join(Re_j, by = 'date', suffix = c(".i", ".j")) %>%
    mutate(se = (median_R_mean.i - median_R_mean.j)^2,
           rele = abs((median_R_mean.j - median_R_mean.i)/median_R_mean.j),
           coverage = (median_R_mean.i > median_R_lowHPD.j) & (median_R_mean.i < median_R_highHPD.j) ) 
  
  se = compare_df %>% pull(se)
  rele = compare_df %>% pull(rele)
  coverage = compare_df %>% pull(coverage) %>% sum(na.rm = T) /length(Re_i$date)
  
  rmse = sqrt(sum(se, na.rm = T)/length(Re_i$date))
  mape = sum(rele, na.rm = T)/length(Re_i$date)
  
  return(c(rmse, coverage, mape))
}

