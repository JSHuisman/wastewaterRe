###########################################################
# Analysis of Swiss wastewater data
# Author: J.S. Huisman
###########################################################
# This script first loads Zurich wastewater data (public)
# and catchment specific case data (not public)
# From these two sources of data, we then estimate Re
# This is summarized in Fig 1 from the paper
#
# Second, we scan across SLDs to find the best fit between
# both traces, here Rww is re-estimated each time (Fig 3)
#
# Third, we investigate the influence of different SLDs
# on Rww estimates (shown in supplement)
#
# Lastly, we subsample the wastewater data by different
# weekdays, and compare the resulting Rww's
###########################################################
library(tidyverse)
library(lubridate)
library(patchwork)
library(viridis)
library(EpiEstim)
# the code further requires zoo for data imputation

app_location = '../../covid-19-re-shiny-app'
source(paste0(app_location,'/app/otherScripts/2_utils_getInfectionIncidence.R'))
source(paste0(app_location,'/app/otherScripts/3_utils_doReEstimates.R'))

source('wastewater_functions.R')

plot_dir = '../figures'

theme_set(theme_minimal() +
            theme(
              strip.text = element_text(size=20),
              axis.text= element_text(size=17),
              axis.title =  element_text(size=20),
              legend.text= element_text(size=17),
              legend.title= element_text(size=20)
            ))

###########################################################
### ZURICH ####
ZH_flow_url = "http://parsivel-eawag.ch/sarscov2/__data__/ARA%20Werdhoelzli_flow_cases.csv"
ZH_genes_url = "http://parsivel-eawag.ch/sarscov2/__data__/ARA%20Werdhoelzli_genes.csv"

raw_flow_data_ZH <- read_delim(ZH_flow_url, delim = ';',
                               col_names = c('date', 'cases', 'cases_smooth', 
                                             'flow', 'n1_smooth', 'n2_smooth'),
                               col_types = cols(date = col_date(format = '')),
                               skip = 1) 

raw_gene_data_ZH <- read_delim(ZH_genes_url, delim = ';',
                               col_names = c('date', 'n1', 'n2'),
                               col_types = cols(date = col_date(format = '')),
                               skip = 1) 

# We filter after Sept 1 because the data prior is often below LOQ
# We select data up to January 20th, because the sampling changed after
# 29th of October is removed because of quality control
# Missing data is imputed with linear interpolation
raw_data_ZH <- raw_flow_data_ZH %>%
  left_join(raw_gene_data_ZH, c('date')) %>%
  filter(!is.na(n1),
         date >= as_date("2020-09-01"),
         date <= as_date("2021-01-20"),
         date != as_date("2020-10-29")) %>%
  mutate(orig_data = TRUE) %>%
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>%
  mutate(region = 'ZH')

# Plot `raw' WW data #####

ggplot(raw_data_ZH, aes(x=date, y = n2)) +
  geom_point(colour = 'blue') +
  geom_line(colour = 'black', linetype = 'dashed') +
  geom_smooth(method = 'loess', colour = 'black',
              method.args = list(span = 0.05, degree = 1) ) +
  labs(x = 'Date' , y='N2 load')

###########################################################
##### Catchment Level: Clinical Case Data  #####
## Extract case data ####
orig_cases <- read_csv('../data/catchment/incidence_data.csv') %>%
  filter(date >= as_date("2020-08-15"),
         date <= as_date("2021-01-20")) %>%
  mutate(across(c(-date), ~ ifelse(is.na(.), 0, .)))  %>%
  select(date, confirmed)

# Plot Case incidence data #####

plot_orig_cases <- orig_cases %>% 
  pivot_longer(cols = c(-date))

case_data_plot <- ggplot(plot_orig_cases) +
  geom_bar(aes(x=date, y= value, fill = name), alpha = 0.5,
             position = 'identity', stat = 'identity', show.legend = F) +
  labs(x = 'Date' , y='Cases per day') +
  scale_x_date(limits = c(as_date('2020-08-15'), as_date('2021-01-20')) ) +
  scale_fill_manual(values = c(viridis(4)[1:2])) + 
  labs(colour = 'Variable') 

case_data_plot
ggsave(plot = case_data_plot, paste0(plot_dir, '/', 'Case_data_raw_catchment.pdf'), height = 11, width = 9)

### Deconvolution and Re estimation for Case data ####
deconv_cases <- data.frame()
Re_cases <- data.frame()
for(inc_var in c('confirmed')){
  new_deconv_data = deconvolveIncidence(orig_cases, 
                                        incidence_var = inc_var,
                                        getCountParams('incubation'), 
                                        getCountParams(paste0(inc_var, '_zh')),
                                        smooth_param = TRUE, n_boot = 50)
  
  new_deconv_data <- new_deconv_data 

  new_Re = getReBootstrap(new_deconv_data) 

  deconv_cases <- bind_rows(deconv_cases, new_deconv_data)
  Re_cases = bind_rows(Re_cases, new_Re)
}

## Clean up a bit for later ####
plot_deconv_cases <- deconv_cases %>%
  mutate(region = 'ZH') %>%
  mutate(data_type = recode_factor(data_type,
                            infection_confirmed = "Confirmed cases",
                            infection_hospitalised = "Hospitalized patients",
                            infection_death = "Deaths")
         ) %>%
  #filter(data_type == 'Confirmed cases') %>%
  dplyr::select(-local_infection, -country, -source) %>%
  group_by(date, region, data_type) %>%
  summarise(sd = sd(value),
            value = mean(value),
            .groups = 'drop')

Restimates <- Re_cases %>%
  mutate(region = 'ZH') %>%
  mutate(data_type = recode_factor(data_type,
                            infection_confirmed = "Confirmed cases",
                            infection_hospitalised = "Hospitalized patients",
                            infection_death = "Deaths") ) %>%
  mutate(countryIso3 = 'CHE')

###########################################################
## Normalisation of WW data ####

# At LOQ we see 1 case # Min observed
norm_min <- min(raw_data_ZH$n1, raw_data_ZH$n2)

### ALL Normalised WW DATA ####
ww_data = bind_rows(raw_data_ZH)  %>%
  mutate(norm_n1 = n1/norm_min,
         norm_n2 = n2/norm_min)

# Plot normalised data #####
plot_raw_ww_data <- ww_data %>%
  dplyr::select(date, orig_data, region, N1 = n1, N2 = n2) %>%
  pivot_longer(cols = c(N1, N2)) %>%
  mutate(name_orig = ifelse(!is.na(orig_data), name, 'Imputed'))

ww_data_plot <- ggplot() +
  geom_point(data = plot_raw_ww_data, aes(x=date, y= value, colour = name_orig),
             size = 2, show.legend = F) +
  geom_line(data = plot_raw_ww_data %>% filter(orig_data), 
            aes(x=date, y= value,colour = name), linetype = 'dotted', show.legend = F) +
  labs(x = 'Date' , y='Gene copies per day') +
  scale_x_date(limits = c(as_date('2020-08-15'), as_date('2021-01-20')) ) +
  scale_colour_manual(values = c(viridis(4)[3:4], 'lightgrey'), 
                      labels = c('N1', 'N2', 'Imputed'),
                      breaks = c('N1', 'N2', 'Imputed'),
                      name = 'Variable') + 
  labs(colour = 'Variable') 

ww_data_plot

ggsave(plot = ww_data_plot, paste0(plot_dir, '/', 'Wastewater_data_raw_CH.pdf'), height = 11, width = 9)


###########################################################
##### Deconvolve and Estimate WW Re #####

config_df = expand.grid("region" = c('ZH'),  
                        'incidence_var' = c('norm_n1', 'norm_n2'),
                        'GammaParams' = list(c('incubation', 'benefield')) )


deconv_ww_data <- data.frame()
Re_ww <- data.frame()

for(row_i in 1:nrow(config_df)){
  new_deconv_data = deconvolveIncidence(ww_data %>% filter(region == config_df[row_i, 'region']), 
                                        incidence_var = config_df[row_i, 'incidence_var'],
                                        getCountParams(unlist(config_df[row_i, 'GammaParams'])[1]), 
                                        getCountParams(unlist(config_df[row_i, 'GammaParams'])[2]),
                                        smooth_param = TRUE, n_boot = 50)
  
  new_deconv_data <- new_deconv_data %>%
    mutate(incidence_var = config_df[row_i, 'incidence_var'],
           incubationParams = unlist(config_df[row_i, 'GammaParams'])[1], 
           onsetToCountParams = unlist(config_df[row_i, 'GammaParams'])[2],
           GammaParams = paste0(incubationParams, '_', onsetToCountParams),
           source = GammaParams)
  
  ##### Get Re #####
  new_Re_ww = getReBootstrap(new_deconv_data)
  new_Re_ww <- new_Re_ww %>%
    mutate(variable = config_df[row_i, 'incidence_var'],
           incubationParams = unlist(config_df[row_i, 'GammaParams'])[1], 
           onsetToCountParams = unlist(config_df[row_i, 'GammaParams'])[2],
           GammaParams = paste0(incubationParams, '_', onsetToCountParams),
           region = config_df[row_i, 'region'])
  
  deconv_ww_data <- bind_rows(deconv_ww_data, new_deconv_data)
  Re_ww = bind_rows(Re_ww, new_Re_ww)
}

##### Plot Deconvolved Cases #####
all_deconv_results <- bind_rows(deconv_ww_data, deconv_cases) %>%
  mutate(plot_row = ifelse(data_type == 'infection_confirmed', 'Cases', 'Wastewater'),
         plot_row = factor(plot_row, levels = c('Cases', 'Wastewater')))

mean_deconv_data <- all_deconv_results %>%
  group_by(date, region, country, source, data_type, incidence_var, 
           incubationParams, onsetToCountParams, GammaParams, plot_row) %>%
  summarise(sd = sd(value),
            value = mean(value),
            .groups = 'drop')

deconv_plot <- ggplot() +
  geom_errorbar(data = mean_deconv_data, 
                aes(x=date, ymin = value -sd,  ymax = value +sd, colour = data_type),
                show.legend = T) +
  facet_wrap(vars(plot_row), ncol = 1, scale = 'free_y') +
  labs(x = 'Date' , y='Infections per day') +
  scale_x_date(limits = c(as_date('2020-08-15'), as_date('2021-01-20')) ) +
  scale_colour_manual(values = viridis(4), 
                      breaks = c("infection_confirmed", "x",
                                 "infection_norm_n1", "infection_norm_n2"),
                      labels = c("Confirmed cases", "x", "N1", "N2"),
                      name = 'Variable') + 
  guides(color = guide_legend(override.aes = list(size=5))) +
  theme(
    legend.position = 'bottom',
    panel.spacing.y = unit(2, "lines")
  )

deconv_plot

ggsave(plot = deconv_plot, paste0(plot_dir, '/', 'Wastewater_data_deconv_CH.png'), height = 11, width = 16)

###########################################################
##### Plot Re ####

# Matching Re estimates from report data ####

date_ranges <- Re_ww %>%
  group_by(region) %>%
  summarise(min_date = min(date),
            max_date = max(date)) 

plotData <- Restimates %>%
  filter(region %in% c('ZH', 'VD'),
         estimate_type == 'Cori_slidingWindow',
         date >= date_ranges[date_ranges$region == 'ZH', ]$min_date,
         date <= date_ranges[date_ranges$region == 'ZH',]$max_date ) %>%
  dplyr::select(-estimate_type, -countryIso3, -country, -source) %>%
  mutate(source = 'Case data',
         plot_row = case_when(data_type == 'Confirmed cases' ~ 'top',
                              data_type == 'Hospitalized patients' ~ 'middle',
                              data_type == 'Deaths' ~ 'bot'),
         plot_row = factor(plot_row, levels = c('top', 'middle', 'bottom')))

## Plot Re ####
Re_ww <- Re_ww %>%
  mutate(variable = recode(variable,
                           'norm_n1' = 'N1',
                           'norm_n2' = 'N2'))

Re_plot <- ggplot() +
  geom_ribbon(data = plotData, aes(x = date, ymin = median_R_lowHPD,
                                   ymax = median_R_highHPD, fill = data_type),
              alpha = 0.2, show.legend = F) +
  geom_line(data = plotData,
            aes(x = date, y = median_R_mean, colour = data_type), alpha = 0.7, show.legend = F) +
  geom_ribbon(data = Re_ww, aes(x = date, ymin = median_R_lowHPD,
                                ymax = median_R_highHPD,  fill = variable), alpha = 0.2, show.legend = F) +
  geom_line(data = Re_ww, 
            aes(x = date, y = median_R_mean, colour = variable),alpha = 0.7, show.legend = F) +
  geom_hline(yintercept = 1) +
  facet_grid(cols = vars(region), rows = vars(plot_row), scale = 'free_x') +
  coord_cartesian(ylim = c(0, 2.0)) +
  scale_colour_viridis(discrete = T) + 
  scale_fill_viridis(discrete = T) + 
  scale_x_date(limits = c(as_date('2020-08-15'), as_date('2021-01-20'))
  ) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Observation Type', fill = 'Observation Type') +
  theme(
    panel.spacing.y = unit(2, "lines"),
    strip.text= element_blank()
    )

Re_plot

ggsave(plot = Re_plot, paste0(plot_dir, '/', 'Wastewater_Re_CH.png'), height = 12, width = 14)

##############
# Fig 1 for the paper

(ww_data_plot | case_data_plot ) / ( deconv_plot | Re_plot ) + 
  plot_layout(heights=c(1, 1), guides = "collect") + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 30),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.key.size = unit(3,"line") )

ggsave(paste0(plot_dir, '/', 'Fig1.png'), height = 16, width = 14)

###########################################################
## Compare RMSE of wastewater and case report traces ####

## Compare Rww and Rcc in all combinations ####
all_Re <- as_tibble(Re_ww) %>%
  select(region, data_type, date, median_R_mean, 
         median_R_highHPD, median_R_lowHPD, source) %>%
  mutate(data_type = recode(data_type,
                            'infection_norm_n1' = 'N1',
                            'infection_norm_n2' = 'N2')) %>%
  bind_rows(plotData) %>%
  filter(region == 'ZH') %>%
  select(-region, -source)


data_types = unique(all_Re$data_type)
rmse_matrix = matrix(data = NA, nrow = length(data_types), ncol = length(data_types),
                     dimnames = list(data_types, data_types))
coverage_matrix = matrix(data = NA, nrow = length(data_types), ncol = length(data_types),
                         dimnames = list(data_types, data_types))
mape_matrix = matrix(data = NA, nrow = length(data_types), ncol = length(data_types),
                     dimnames = list(data_types, data_types))

for(ind_i in seq_along(data_types)){
  for (ind_j in seq_along(data_types)){
    if (ind_i == ind_j){next}
    
    subRe_i = all_Re %>% 
      filter(data_type == data_types[ind_i])
    subRe_j = all_Re %>% 
      filter(data_type == data_types[ind_j])
    
    result = compareTraces(subRe_i, subRe_j)
    
    rmse_matrix[ind_i, ind_j] = result[1]
    coverage_matrix[ind_i, ind_j] = result[2]
    mape_matrix[ind_i, ind_j] = result[3]
  }
}

rmse_matrix

######################################################################################################################
#Scan across deconvolution parameters ####

# This takes about 2 hrs per scan
for (incubationParam in c('zero', 'incubation')){
  for (incidence_var in c('norm_n1')){
    for (canton in c('ZH')){
      
      proc_data <- ww_data %>% filter(region == canton)
      plotData_subset <- plotData %>% filter(region == canton)
      
      # the scan values
      meanOpts = seq(0.5, 15, 0.5)
      sdOpts = seq(0.5, 10, 0.5)
      
      deconv_results = cbind(expand_grid(meanOpts, sdOpts), 
                             'rmse_cc' = NA, 'coverage_cc' = NA, 'mape_cc' = NA,
                             'rmse_h' = NA, 'coverage_h' = NA, 'mape_h' = NA)
      
      for (row_id in 1:nrow(deconv_results)){
        deconv_config = try(deconvolveIncidence(proc_data, incidence_var,
                                                getCountParams(incubationParam), 
                                                getGammaParams(deconv_results[row_id, 'meanOpts'],
                                                               deconv_results[row_id, 'sdOpts']),
                                                smooth_param = TRUE, n_boot = 50))
        
        if('try-error' %in% class(deconv_config)){
          deconv_results[row_id, c('rmse_cc', 'coverage_cc', 'mape_cc')] = c(Inf, 0, Inf)
          deconv_results[row_id, c('rmse_h', 'coverage_h', 'mape_h')] = c(Inf, 0, Inf)
          next
        }
        
        Re_config = getReBootstrap(deconv_config)
        
        deconv_results[row_id, c('rmse_cc', 'coverage_cc', 'mape_cc')] = compareTraces(Re_config, plotData_subset %>%
                                                                                         filter(data_type == 'Confirmed cases'))
      }
      
      write_csv(deconv_results, paste0('../scan/deconv_', incubationParam, '_', 
                                       canton, '_', incidence_var, '.csv'))
      
    }
  }
}

###############
# Results of the scan

scan_options <- expand.grid(incubationParam = c('zero', 'incubation'),
                            incidence_var = c('norm_n1'),
                            canton = c('ZH'))

for (row_i in 1:nrow(scan_options)){
  deconv_results <- read_csv(paste0('../scan/deconv_', 
                                    scan_options[row_i, 'incubationParam'], '_', 
                                    scan_options[row_i, 'canton'], '_', 
                                    scan_options[row_i, 'incidence_var'], '.csv'),
                             col_types = cols(
                               meanOpts = col_double(),
                               sdOpts = col_double(),
                               rmse_cc = col_double(),
                               coverage_cc = col_double(),
                               mape_cc = col_double(),
                               rmse_h = col_double(),
                               coverage_h = col_double(),
                               mape_h = col_double()
                             ) )
  
  ######## Plots
  methods = setNames(c("RMSE", "Coverage", "MAPE",
                       "RMSE", "Coverage", "MAPE"),
                     c("rmse_cc", "coverage_cc", "mape_cc",
                       "rmse_h", "coverage_h", "mape_h"))
  
  for (method in names(methods)){
    plot <- ggplot(deconv_results) +
      geom_tile(aes(x = meanOpts, y = sdOpts, fill = get(method) )) +
      labs(x = 'Mean', y = 'Standard deviation', fill = methods[method])
    
    if (method %in% c('coverage_cc', 'coverage_h')){
      plot +
        scale_fill_viridis()
    } else {
      plot +
        scale_fill_viridis(direction = -1)
    }
    
    ggsave(paste0(plot_dir, '/', 'Scan_', method, '_',
                  scan_options[row_i, 'canton'], '_',
                  scan_options[row_i, 'incidence_var'], '_',
                  scan_options[row_i, 'incubationParam'], '.png'), height = 10, width = 14)
  }
  
  
  print('##################')
  print(paste(scan_options[row_i, 'incubationParam'], scan_options[row_i, 'canton'],
              scan_options[row_i, 'incidence_var']))
  print(deconv_results[which.min(deconv_results$rmse_cc), ])
  print(deconv_results[which.max(deconv_results$coverage_cc), ])
  print(deconv_results[which.min(deconv_results$mape_cc), ])
}

## Looking at the marginal distributions ####
deconv_results <- read_csv(paste0('../scan/deconv_', 
                                  'incubation_', 
                                  'ZH_', 
                                  'norm_n1.csv'),
                           col_types = cols(
                             meanOpts = col_double(),
                             sdOpts = col_double(),
                             rmse_cc = col_double(),
                             coverage_cc = col_double(),
                             mape_cc = col_double(),
                             rmse_h = col_double(),
                             coverage_h = col_double(),
                             mape_h = col_double()
                           ) )

# Mean
deconv_results %>%
  select(-c("rmse_h", "coverage_h", "mape_h")) %>%
  group_by(meanOpts) %>%
  summarise(across(c(-sdOpts), .fns = list(mean = mean)),
            .groups = 'drop') %>%
  mutate(across(c(-meanOpts, -coverage_cc_mean),
                ~1-.)) %>%
  summarise(across(c(-meanOpts), .fns = list(max = ~ meanOpts[which.max(.)], 
                                             low_ci = ~ meanOpts[which(. == (.[. >= quantile(., probs =c(0.95)) ])[1] )],
                                             high_ci = ~ meanOpts[which(. == (.[. >= quantile(., probs =c(0.95)) ])[-1] )]
                                             ) )) %>%
  pivot_longer(cols = everything())

# SD
deconv_results %>%
  select(-c("rmse_h", "coverage_h", "mape_h")) %>%
  group_by(sdOpts) %>%
  summarise(across(c(-meanOpts), .fns = list(mean = mean)),
            .groups = 'drop') %>%
  mutate(across(c(-sdOpts, -coverage_cc_mean),
                ~1-.)) %>%
  summarise(across(c(-sdOpts), .fns = list(max = ~ meanOpts[which.max(.)],
                                           low_ci = ~ meanOpts[which(. == min(.[. >= quantile(., probs =c(0.95)) ]) )],
                                           high_ci = ~ meanOpts[which(. == max(.[. >= quantile(., probs =c(0.95)) ]) )] 
                                           )) ) %>%
  pivot_longer(cols = everything())


######################################################################################################################
# Results for a number of SLDs - Fig SX ####

# We select data up to January 20th, because the sampling changed after
ww_data = raw_data_ZH %>%
  dplyr::select(date, n1, orig_data, region) %>%
  filter(date <= as_date("2021-01-20")) %>%
  mutate(norm_n1 = n1/5e10,
         norm_n1_high = n1/1e12,
         norm_n1_low = n1/1e5)


config_df = expand.grid("region" = c('ZH'), 'incidence_var' = c('n1','norm_n1', 'norm_n1_high', 'norm_n1_low'),
                        'GammaParams' = list(c('incubation', 'zero'),
                          c('incubation', 'death'), 
                          c('incubation', 'han') ,
                          c('incubation', 'benefield')
                        ) )


deconv_ww_data <- data.frame()
Re_ww <- data.frame()

for(row_i in 1:nrow(config_df)){
  new_deconv_data = deconvolveIncidence(ww_data %>% filter(region == config_df[row_i, 'region']), 
                                        incidence_var = config_df[row_i, 'incidence_var'],
                                        getCountParams(unlist(config_df[row_i, 'GammaParams'])[1]), 
                                        getCountParams(unlist(config_df[row_i, 'GammaParams'])[2]),
                                        smooth_param = TRUE, n_boot = 50)
  
  new_deconv_data <- new_deconv_data %>%
    mutate(incidence_var = config_df[row_i, 'incidence_var'],
           incubationParams = unlist(config_df[row_i, 'GammaParams'])[1], 
           onsetToCountParams = unlist(config_df[row_i, 'GammaParams'])[2],
           GammaParams = paste0(incubationParams, '_', onsetToCountParams),
           source = GammaParams)
  
  ##### Get Re #####
  new_Re_ww = getReBootstrap(new_deconv_data)
  new_Re_ww <- new_Re_ww %>%
    mutate(variable = config_df[row_i, 'incidence_var'],
           incubationParams = unlist(config_df[row_i, 'GammaParams'])[1], 
           onsetToCountParams = unlist(config_df[row_i, 'GammaParams'])[2],
           GammaParams = paste0(incubationParams, '_', onsetToCountParams),
           region = config_df[row_i, 'region'])
  
  deconv_ww_data <- bind_rows(deconv_ww_data, new_deconv_data)
  Re_ww = bind_rows(Re_ww, new_Re_ww)
}

## Deconvolution Plot ####
mean_deconv_data <- deconv_ww_data %>%
  group_by(date, region, country, source, data_type, incidence_var, 
           incubationParams, onsetToCountParams, GammaParams) %>%
  summarise(sd = sd(value),
            value = mean(value),
            .groups = 'drop')

deconv_plot <- ggplot() +
  geom_errorbar(data = mean_deconv_data %>% filter(incidence_var == 'norm_n1'), 
                aes(x=date, ymin = value -sd,  ymax = value +sd, colour = GammaParams),
                show.legend = F) +
  #facet_wrap(vars(incidence_var), scale = 'free_y') +
  labs(x = 'Date' , y='Infections per day') +
  scale_x_date(date_breaks = "4 weeks", 
               date_labels = '%b\n%d',
               limits = c(as_date('2020-08-15'), as_date('2021-01-20'))
  ) +
  scale_colour_manual(values = viridis(4),
                      breaks = c("incubation_zero","incubation_han",
                                 "incubation_benefield", "incubation_death"),
                      labels = c("Incubation only", "Incubation + Han", "Incubation + Benefield", "Incubation + Death"),
                      name = 'Variable') +
  theme(
    strip.text.x= element_blank(),
    axis.title.x =  element_blank()
  )

deconv_plot

## Re Plot ####

Re_plot <- ggplot(Re_ww %>% filter(data_type == 'infection_norm_n1')) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                                ymax = median_R_highHPD,  fill = GammaParams), alpha = 0.2, show.legend = T) +
  geom_line(aes(x = date, y = median_R_mean, colour = GammaParams), alpha = 0.7, show.legend = F) +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim = c(0.5, 2)) +
  scale_x_date(date_breaks = "4 weeks", 
               date_labels = '%b\n%d',
               limits = c(as_date('2020-08-15'), as_date('2021-01-20'))
  ) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Observation Type', fill = 'Observation Type') +
  scale_fill_manual(values = viridis(4),
                      breaks = c("incubation_zero","incubation_han",
                                 "incubation_benefield", "incubation_death"),
                      labels = c("Incubation only", "Incubation + Han", "Incubation + Benefield", "Incubation + Death"),
                      name = 'Variable') +
  scale_colour_manual(values = viridis(4),
                      breaks = c("incubation_zero","incubation_han",
                                 "incubation_benefield", "incubation_death"),
                      labels = c("Incubation only", "Incubation + Han", "Incubation + Benefield", "Incubation + Death"),
                      name = 'Variable') +
  theme(
    panel.spacing.y = unit(2, "lines"),
    #strip.text= element_blank(),
    legend.position = 'bottom'
  )

Re_plot

## Combine plots ####

deconv_plot + Re_plot + 
  plot_layout(ncol=1, heights=c(1, 2)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
ggsave(paste0(plot_dir, '/', 'Supp_Fig_SLDs.png'), height = 16, width = 14)

######################################################################################################################
# Subsample across weekdays ####

daily_data_ZH <- raw_gene_data_ZH %>%
  filter(date > as_date("2020-11-22") & date < as_date("2021-01-11")) %>%
  mutate(day = wday(date),
         region = 'ZH') %>% # 1 is Monday
  mutate(norm_n1 = n1/min(raw_data_ZH$n1, raw_data_ZH$n2),
         norm_n2 = n2/min(raw_data_ZH$n1, raw_data_ZH$n2))

sample_methods = list('MTWTFSS' = c(1, 2, 3, 4, 5, 6, 7),
                      'MTWTF' = c(1, 2, 3, 4, 5), 'MWF' = c(1, 3, 5),
                      'TTS' = c(2, 4, 6), 'MT' = c(1, 4), 'TF' = c(2, 5), 
                      'WS' = c(3, 6),
                      'M' = c(1), 'W' = c(3), 'F' = c(5))

sampled_Restimates <- data.frame()

for (s_ind in 1:length(sample_methods)){
  subsampled_data <- daily_data_ZH %>%
    filter(day %in% sample_methods[[s_ind]]) %>%
    mutate(orig_data = TRUE) %>%
    complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
    mutate(across(c(norm_n1, norm_n2), ~ zoo::na.approx(.x, na.rm = F) )) %>%
    mutate(region = 'ZH') 
  
  for (incidence_var in c('norm_n1')){
    deconv_config = try(deconvolveIncidence(subsampled_data, incidence_var,
                                            getCountParams('incubation'), 
                                            getCountParams('benefield'),
                                            smooth_param = TRUE, n_boot = 50))
    
    Re_config = getReBootstrap(deconv_config) %>%
      mutate(sampling = names(sample_methods)[s_ind])
    
    write_csv(Re_config, paste0('../subsampling/', names(sample_methods)[s_ind], '_',
                                incidence_var, '.csv'))
    
    sampled_Restimates <- bind_rows(sampled_Restimates, Re_config)
  }
}

ZH_sampled_Restimates <- sampled_Restimates %>%
  mutate(category = str_length(sampling),
         sampling = factor(sampling, levels = c('M', 'W', 'F', 'MT', 'TF', 'WS', 
                                                'MWF', 'TTS', 'MTWTF', 'MTWTFSS'))) %>% 
  filter(data_type == 'infection_norm_n1')


ZH_subsample_plot <- ggplot(ZH_sampled_Restimates) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, fill = sampling),
              alpha = 0.4, show.legend = F) +
  geom_line(aes(x = date, y = median_R_mean, colour = sampling), 
            alpha = 1, show.legend = T) +
  geom_hline(yintercept = 1) +
  facet_grid(rows = vars(category)) +
  coord_cartesian(ylim = c(0.5, 1.5)) +
  scale_colour_manual(values = viridis(10)[c(1,4,9,2,5,10,3,8,7,6)]) + 
  scale_fill_manual(values = viridis(10)[c(1,4,9,2,5,10,3,8,7,6)]) +
  scale_x_date(date_breaks = "1 weeks", 
               date_labels = '%b\n%d'#,
               #limits = c(as_date("2020-09-01"), as_date("2021-01-20"))
  ) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Sampling Type', fill = 'Sampling Type') +
  guides(color = guide_legend(override.aes = list(size=5))) 
ZH_subsample_plot

ggsave(paste0(plot_dir, '/', 'Re_subsampling_ZH_N1.png'), height = 12, width = 14)

## Comparison of RMSE ####

rmse_matrix = matrix(data = NA, nrow = length(sample_methods), ncol = length(sample_methods),
                     dimnames = list(names(sample_methods), names(sample_methods)))
coverage_matrix = matrix(data = NA, nrow = length(sample_methods), ncol = length(sample_methods),
                     dimnames = list(names(sample_methods), names(sample_methods)))
mape_matrix = matrix(data = NA, nrow = length(sample_methods), ncol = length(sample_methods),
                     dimnames = list(names(sample_methods), names(sample_methods)))

for(ind_i in seq_along(sample_methods)){
  for (ind_j in seq_along(sample_methods)){
    if (ind_i >= ind_j){next}
    
    subRe_i = sampled_Restimates %>% 
      filter(sampling == names(sample_methods)[ind_i])
    subRe_j = sampled_Restimates %>% 
      filter(sampling == names(sample_methods)[ind_j])
    
    result = compareTraces(subRe_i, subRe_j)
    
    rmse_matrix[ind_i, ind_j] = result[1]
    coverage_matrix[ind_i, ind_j] = result[2]
    mape_matrix[ind_i, ind_j] = result[3]
  }
}

rmse_matrix

formatC(rmse_matrix, digits = 2, format = "f")
formatC(coverage_matrix, digits = 2, format = "f")
formatC(mape_matrix, digits = 2, format = "f")
