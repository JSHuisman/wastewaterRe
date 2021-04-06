###########################################################
# Californian wastewater data
# Author: J.S. Huisman
###########################################################
# This script first loads Rcc estimates from CalCat (public)
# and case data (public) to estimate Rcc
# Then, we load San Jose wastewater data (not public)
# From these two sources of data, we estimate Re (Fig 2)
#
# Second, we scan across SLDs to find the best fit between
# both traces, here Rww is re-estimated each time (Fig 3)
#
# Lastly, we subsample the wastewater data by different
# weekdays, and compare the resulting Rww's
###########################################################
library(tidyverse)
library(lubridate)
library(viridis)
library(patchwork)
library(EpiEstim)

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
###### Rt data #####

R_data_raw <- read_csv('../data/CalCAT_Santa_Clara.csv', skip = 2,
                       col_types = cols(
                         Date = col_date(format = ""),
                         CovidActNow = col_double(),
                         `covid19-projections.com` = col_double(),
                         UCLA = col_double(),
                         LEMMA = col_double(),
                         Stanford = col_double(),
                         UCSF = col_double(),
                         Harvard = col_double(),
                         Ensemble = col_double()
                       ))

# We exclude SEIR models, and the ensemble because it is not
# properly combined by date of infection
# JHU is not documented on the CalCat page
methods_to_exclude = c("LEMMA", "Stanford", "UCLA", "Ensemble", "JHU")

R_data <- R_data_raw %>%
  filter(!is.na(Date)) %>%
  pivot_longer(cols = c(-Date), names_to = 'Method', values_to = 'R_estimate') %>%
  filter(!Method %in% methods_to_exclude)

###### Case data #####

urlfile <- "https://data.chhs.ca.gov/dataset/f333528b-4d38-4814-bebb-12db1f10f535/resource/046cdd2b-31e5-4d34-9ed3-b48cdbc4be7a/download/covid19cases_test.csv"

cali_case_data <- read_csv(urlfile)
county_list = "Santa Clara"

ggplot(cali_case_data %>% filter(area %in% county_list, date >= as_date("2020-11-01")) )  +
  geom_bar(aes(x=date, y= abs(reported_cases)), stat = "identity", fill = 'red', size = 1.5) +
  labs(x = 'Date' , y='COVID cases and deaths') +
  scale_x_date(date_breaks = "4 weeks",
               date_labels = '%b\n%d'
  )

## Clean data 
first_clean_cali_case_data <- cali_case_data %>%
  filter(area %in% county_list,
         !is.na(date)) %>%
  mutate(confirmed = ifelse(reported_cases <0, 0, reported_cases),
         deaths = ifelse(reported_deaths <0, 0, reported_deaths),
         county = area) %>%
  mutate(orig_data = TRUE) %>%
  select(date, county, confirmed, deaths, orig_data) %>%
  filter(date >= as_date("2020-11-01"))

clean_cali_case_data <- first_clean_cali_case_data %>%
  mutate(location = recode(county,
                           "Santa Clara" = "SanJose") )


###### Plot #######################

plot_orig_cases <- clean_cali_case_data %>% 
  pivot_longer(cols = c('confirmed'))

case_data_plot <- ggplot(plot_orig_cases) +
  geom_bar(aes(x=date, y= value, fill = name), alpha = 0.5,
           position = 'identity', stat = 'identity', show.legend = F) +
  labs(x = 'Date' , y='Cases per day') +
  scale_x_date(limits = c(as_date('2020-11-15'), as_date('2021-03-18')) ) +
  scale_fill_manual(values = c(viridis(4)[1])) + 
  labs(colour = 'Variable') 

case_data_plot
ggsave(plot = case_data_plot, paste0(plot_dir, '/', 'Case_data_Cali.pdf'), height = 11, width = 9)



## Confirmed case Rt #####

deconv_result_c <- data.frame()
result_c <- data.frame()
for (location_i in unique(clean_cali_case_data$county)){
  for (incidence_var_i in c('confirmed')){
    select_data <- clean_cali_case_data %>%
      filter(county == location_i,
             !is.na(get(incidence_var_i)) ) %>%
      mutate(region = county) #%>%
    
    if (incidence_var_i == 'confirmed'){
      delayParams = getCountParams('confirmed_cali')
    } else if (incidence_var_i == 'deaths'){
      delayParams = getCountParams('death')
    }
    #delayParams = getCountParams('zero')
    
    new_deconv_result = deconvolveIncidence(select_data,
                                            incidence_var = incidence_var_i,
                                            getCountParams('incubation'),
                                            delayParams,
                                            smooth_param = TRUE, n_boot = 50) %>%
      mutate(county = location_i,
             data_type = incidence_var_i)
    
    new_result = getReBootstrap(new_deconv_result)
    new_result = new_result %>%
      mutate(county = location_i,
             data_type = incidence_var_i)
    new_result['variable'] = incidence_var_i
    
    deconv_result_c = bind_rows(deconv_result_c, new_deconv_result)
    result_c = bind_rows(result_c, new_result)
    
  }
}

deconv_result_c['location'] = 'SanJose'

###########################################################
## BOEHM - SCAN Pilot Project ####

xls_path = '../data/DATA_SCAN_pilot_project.xlsx'

raw_data <- readxl::read_excel(xls_path, sheet = "SanJose", na = c("", "ND") ) %>%
    mutate(date = as_date(Date),
           location = "SanJose") %>%
    select(date, location, n_gene = `N Gene gc/g dry weight`, 
           n_top = `N Gene gc/g dry weight UCI`,
           n_bot = `N Gene gc/g dry weight LCI`, 
           ORF1a = `ORF1a gc/g dry weight`,
           s_gene = `S Gene gc/g dry weight`, 
           s_top = `S Gene gc/g dry weight UCI`, s_bot = `S Gene gc/g dry weight LCI`,
           bcov = `BCoV Recovery`, pmmov = `PMMoV gc/g dry weight`) %>%
    mutate(n_orig = n_gene,
           s_orig = s_gene,
           orf1a_orig = ORF1a) %>%
    filter(!is.na(date) ) 

## Normalisation ####
norm_min <- raw_data %>%
  filter(n_gene != 0) %>%
  pull(n_gene) %>% min()

## Final Cleaning ####
#  gaps in the data prior to 15 Nov.
# stop on 18.3.2021
# 03 January, 18 February do not meet quality control
clean_data <- raw_data %>%
    mutate(orig_data = TRUE) %>%
    filter(date >= as_date('2020-11-15'),
         date <= as_date('2021-03-18'),
         date != as_date('2021-01-03'),
         date != as_date('2021-02-18') ) %>%
    complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
    mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>%
    mutate(across(c(n_gene, s_gene, ORF1a), ~ ./norm_min ) ) %>%
  mutate(county = 'Santa Clara',
         region = county)

#### Plot wastewater data ####
plot_raw_ww_data <- clean_data %>%
  #select(date, orig_data, location, N = n_gene, S = s_gene, ORF1a) %>%
  select(date, orig_data, location, N = n_orig, S = s_orig, ORF1a = orf1a_orig) %>%
  pivot_longer(cols = c(N, S, ORF1a)) %>%
  mutate(name_orig = ifelse(!is.na(orig_data), name, 'Imputed'))

SCANpilot_data_plot <- ggplot() +
  geom_point(data = plot_raw_ww_data, aes(x=date, y= value, colour = name_orig),
             size = 2, show.legend = F) +
  geom_line(data = plot_raw_ww_data %>% filter(orig_data), 
            aes(x=date, y= value, colour = name), linetype = 'dotted', show.legend = F) +
  labs(x = 'Date' , y='Gene copies per gr dry weight') +
  scale_x_date(limits = c(as_date('2020-11-15'), as_date('2021-03-18')) 
  ) +
  scale_colour_manual(values = c(viridis(4)[2:4], 'lightgrey'), 
                      labels = c('N', 'S', 'ORF1a', 'Imputed'),
                      breaks = c('N', 'S', 'ORF1a', 'Imputed'),
                      name = 'Variable') + 
  labs(colour = 'Variable') +
  theme(
    axis.title.x =  element_blank(),
  )
SCANpilot_data_plot

ggsave(plot = SCANpilot_data_plot, paste0(plot_dir, '/', 'Wastewater_data_raw_Cali.png'), height = 11, width = 9)

###### Deconvolve #####

deconv_result <- data.frame()
result <- data.frame()
for (incidence_var_i in c('n_gene', 's_gene', 'ORF1a')){
  new_deconv_result = deconvolveIncidence(clean_data, 
                                          incidence_var = incidence_var_i,
                                          getCountParams('incubation'), 
                                          getCountParams('benefield'),
                                          smooth_param = TRUE, n_boot = 50) %>%
    mutate(data_type = incidence_var_i)
  
  new_result = getReBootstrap(new_deconv_result)
  new_result = new_result %>%
    mutate(data_type = incidence_var_i)
  new_result['variable'] = incidence_var_i
  
  deconv_result = bind_rows(deconv_result, new_deconv_result)
  result = bind_rows(result, new_result)
  
}


##### Plot #####
all_deconv_results <- bind_rows(deconv_result, deconv_result_c) %>%
  filter(date >= as_date("2020-11-15"),
         date <= as_date("2021-03-08")) %>%
  mutate(plot_row = ifelse(data_type == 'confirmed', 'Cases', 'Wastewater'),
         plot_row = factor(plot_row, levels = c('Cases', 'Wastewater')))

mean_deconv_data <- all_deconv_results %>%
  group_by(date, region, country, data_type, plot_row, location) %>%
  summarise(sd = sd(value),
            value = mean(value),
            .groups = 'drop')

SCANpilot_deconv_plot <- ggplot() +
  geom_errorbar(data = mean_deconv_data, 
                aes(x=date, ymin = value -sd,  ymax = value +sd, colour = data_type),
                show.legend = F) +
  labs(x = 'Date' , y='Infections per day') +
  scale_x_date(limits = c(as_date('2020-11-15'), as_date('2021-03-18')) ) +
  facet_wrap(vars(plot_row), ncol = 1, scale = 'free_y') +
  #coord_cartesian(ylim = c(0, 500)) +
  scale_colour_manual(values = viridis(4), 
                      breaks = c("confirmed",
                                 "n_gene", "s_gene", 'ORF1a'),
                      labels = c("Confirmed cases", "N", "S", 'ORF1a'),
                      name = 'Variable')
SCANpilot_deconv_plot

ggsave(plot = SCANpilot_deconv_plot,
       paste0(plot_dir, '/', 'Deconv_SCANpilot.png'), height = 10, width = 14)
ggsave(paste0(plot_dir, '/', 'Deconvolution_Boehm.pdf'), height = 10, width = 15)


#### R estimate Plots ########
all_results <- result_c %>%
  mutate(location = 'SanJose')%>%
  filter(date >= as_date("2020-11-15"),
         date <= as_date("2021-03-08")) %>%
  bind_rows(result)
# the deconvolution is 10 days prior to last data

SCANpilot_Re_plot <- ggplot(all_results) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                  ymax = median_R_highHPD,  fill = variable), alpha = 0.3, show.legend = F) +
  geom_line(aes(x = date, y = median_R_mean, colour = variable), show.legend = T) +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim = c(0, 2)) +
  scale_x_date(limits = c(as_date('2020-11-15'), as_date('2021-03-18')) ) +
  scale_colour_manual(values = viridis(4), 
                      breaks = c("confirmed",
                                 "n_gene", "s_gene", 'ORF1a'),
                      labels = c("Confirmed cases", "N", "S", 'ORF1a'),
                      name = 'Variable') +
  scale_fill_manual(values = viridis(4), 
                    breaks = c("confirmed",
                               "n_gene", "s_gene", 'ORF1a'),
                    labels = c("Confirmed cases", "N", "S", 'ORF1a'),
                    name = 'Variable') +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Observation Type', fill = 'Observation Type') +
  guides(color = guide_legend(override.aes = list(size=5))) +
  theme(
    legend.position = 'bottom' 
  )
SCANpilot_Re_plot 

ggsave(plot = SCANpilot_Re_plot, paste0(plot_dir, '/', 'Wastewater_Re_Boehm.pdf'), height = 14, width = 12)
ggsave(plot = SCANpilot_Re_plot, paste0(plot_dir, '/', 'Wastewater_Re_SCANpilot.png'), height = 14, width = 12)

###########################################################
# Fig 2 for the paper


(SCANpilot_data_plot | case_data_plot ) / (SCANpilot_deconv_plot | SCANpilot_Re_plot ) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 30),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.key.size = unit(3,"line") )
ggsave(paste0(plot_dir, '/', 'Fig_Cali_SCANpilot.png'), height = 16, width = 14)

#ggsave(paste0(plot_dir, '/', 'Fig_Cali_SCANpilot_SCdata.png'), height = 16, width = 14)

###########################################################
## With all official R ####

ggplot(result %>% filter(variable == 'n_gene')) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                  ymax = median_R_highHPD), fill = viridis(4)[2], alpha = 0.2, show.legend = F) +
  geom_line(aes(x = date, y = median_R_mean), colour = viridis(4)[2], show.legend = F) +
  geom_line(data = R_data, aes(x = Date, y = R_estimate, colour = Method), size = 1, show.legend = F )+
  geom_ribbon(data = result_c %>% mutate(Method = 'ETH'), aes(x = date, ymin = median_R_lowHPD,
                                                              ymax = median_R_highHPD), fill = viridis(4)[1], alpha = 0.5, show.legend = F) +
  geom_line(data = result_c %>% mutate(Method = 'ETH'), aes(x = date, y = median_R_mean), colour = viridis(4)[1], show.legend = F) +
  #geom_ribbon(data = R_data,aes(x = Date, ymin = min, ymax = max), alpha = 0.2, fill = 'black') +
  facet_wrap(vars(Method), ncol = 2) +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim = c(0.5, 1.5)) +
  scale_colour_viridis(discrete = T, begin = 0.3) + 
  #scale_fill_viridis(discrete = T) + 
  scale_x_date(limits = c(as_date('2020-11-18'), as_date('2021-03-15')) ) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Observation Type', fill = 'Observation Type') +
  guides(color = guide_legend(override.aes = list(size=5))) +
  theme(
    legend.position = 'bottom' #c(0.9, 0.1)
  )

ggsave(paste0(plot_dir, '/', 'Cali_Re_WW_CalCat.png'), height = 12, width = 10)

######################################################################################################################
##### Scan ##################
for(incubationParam in c('zero', 'incubation')){
  for (incidence_var in c('s_gene')){
    
    meanOpts = seq(0.5, 15, 0.5)
    sdOpts = seq(0.5, 10, 0.5)
    
    deconv_results = cbind(expand_grid(meanOpts, sdOpts), 
                           'rmse' = NA, 'coverage' = NA, 'mape' = NA)
    
    for (row_id in 1:nrow(deconv_results)){
      deconv_config = try(deconvolveIncidence(clean_data, incidence_var,
                                              getCountParams(incubationParam), 
                                              getGammaParams(deconv_results[row_id, 'meanOpts'],
                                                             deconv_results[row_id, 'sdOpts']),
                                              smooth_param = TRUE, n_boot = 50))
      
      if('try-error' %in% class(deconv_config)){
        deconv_results[row_id, c('rmse', 'coverage', 'mape')] = c(Inf, 0, Inf)
        next
      }
      
      Re_config = getReBootstrap(deconv_config)
      
      deconv_results[row_id, c('rmse', 'coverage', 'mape')] = compareTraces(Re_config, result_c)
    }
    
    write_csv(deconv_results, paste0('../scan/deconv_Cali_', incubationParam, '_', 
                                     incidence_var, '_SCdata.csv'))
    
  }
}

## Results of scan ####
scan_options <- expand.grid(incubationParam = c('incubation', 'zero'),
                            incidence_var = c('s_gene'))

for (row_i in 1:nrow(scan_options)){
  deconv_results <- read_csv(paste0('../scan/deconv_Cali_', 
                                    scan_options[row_i, 'incubationParam'], '_', 
                                    scan_options[row_i, 'incidence_var'], '_SCdata.csv'),
                             col_types = cols(
                               meanOpts = col_double(),
                               sdOpts = col_double(),
                               rmse = col_double(),
                               coverage = col_double(),
                               mape = col_double()
                             ) )
  
  ######## Plots
  methods = setNames(c("RMSE", "Coverage", "MAPE"),
                     c("rmse", "coverage", "mape"))
  
  for (method in names(methods)){
    plot <- ggplot(deconv_results) +
      geom_tile(aes(x = meanOpts, y = sdOpts, fill = get(method) )) +
      labs(x = 'Mean', y = 'Standard deviation', fill = methods[method])
    
    if (method %in% c('coverage')){
      plot +
        scale_fill_viridis()
    } else {
      plot +
        scale_fill_viridis(direction = -1)
    }
    
    ggsave(paste0(plot_dir, '/', 'Scan_Cali_', method, '_',
                  scan_options[row_i, 'incidence_var'], '_',
                  scan_options[row_i, 'incubationParam'], '_SCdata.png'), height = 10, width = 14)
  }
  
  
  print('##################')
  print(paste(scan_options[row_i, 'incidence_var'], ' ', scan_options[row_i, 'incubationParam']))
  print(deconv_results[which.min(deconv_results$rmse), ])
  print(deconv_results[which.max(deconv_results$coverage), ])
  print(deconv_results[which.min(deconv_results$mape), ])
}

###########################################################
# Plot scan together with the one for Zurich, Fig 3 paper ####

ZH_deconv_results <- read_csv(paste0('../scan/deconv_', 
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
  
Cali_deconv_results <- read_csv(paste0('../scan/deconv_Cali_', 
                                  'incubation_', 
                                  's_gene.csv'),
                           col_types = cols(
                             meanOpts = col_double(),
                             sdOpts = col_double(),
                             rmse = col_double(),
                             coverage = col_double(),
                             mape = col_double()
                           ) ) %>% 
  mutate(city = 'San Jose')

all_deconv_results <- ZH_deconv_results %>%
  mutate(city = 'Zurich') %>%
  select(city, meanOpts, sdOpts, rmse = rmse_cc, coverage = coverage_cc, mape = mape_cc) %>%
  bind_rows(Cali_deconv_results) %>%
  mutate(city = factor(city, levels = c('Zurich', 'San Jose')))

combined_plot <- ggplot(all_deconv_results) +
    geom_tile(aes(x = meanOpts, y = sdOpts, fill = rmse )) +
  facet_wrap(vars(city)) +
    labs(x = 'Mean', y = 'Standard deviation', fill = 'RMSE') +
    scale_fill_viridis(direction = -1)
combined_plot

#ggsave(paste0(plot_dir, '/', 'Scan_Cali_ZH.png'), height = 7, width = 20)
ggsave(paste0(plot_dir, '/', 'Scan_Cali_ZH_incubation.png'), height = 7, width = 20)

###########################################################
## Looking at the marginal distributions ####
deconv_results <- read_csv(paste0('../scan/deconv_Cali_', 
                                  'incubation_', 
                                  's_gene.csv'),
                           col_types = cols(
                             meanOpts = col_double(),
                             sdOpts = col_double(),
                             rmse = col_double(),
                             coverage = col_double(),
                             mape = col_double()
                           ) )

# Mean
deconv_results %>%
  group_by(meanOpts) %>%
  summarise(across(c(-sdOpts), .fns = list(mean = mean)),
            .groups = 'drop') %>%
  mutate(across(c(-meanOpts, -coverage_mean),
                ~1-.)) %>%
  summarise(across(c(-meanOpts), .fns = list(max = ~ meanOpts[which.max(.)],
                                             #low_ci = ~ meanOpts[which.min(abs(. - quantile(., probs =c(0.05)) )) ],
                                             #high_ci = ~ meanOpts[which.min(abs(. - quantile(., probs =c(0.95)) )) ] 
                                             low_ci = ~ meanOpts[which(. == (.[. >= quantile(., probs =c(0.95)) ])[1] )],
                                             high_ci = ~ meanOpts[which(. == (.[. >= quantile(., probs =c(0.95)) ])[-1] )] 
                                             #top = ~ meanOpts[which.min(abs(. - (max(.) + 1.96*sd(.)) ))],
                                             #bot = ~ meanOpts[which.min(abs(. - (max(.) - 1.96*sd(.)) ))]
  ) )) %>%
  pivot_longer(cols = everything())

# SD
deconv_results %>%
  group_by(sdOpts) %>%
  summarise(across(c(-meanOpts), .fns = list(mean = mean)),
            .groups = 'drop') %>%
  mutate(across(c(-sdOpts, -coverage_mean),
                ~1-.)) %>%
  summarise(across(c(-sdOpts), .fns = list(max = ~ meanOpts[which.max(.)],
                                           #low_ci = ~ meanOpts[which.min(abs(. - quantile(., probs =c(0.05)) )) ],
                                           #high_ci = ~ meanOpts[which.min(abs(. - quantile(., probs =c(0.95)) )) ] 
                                           low_ci = ~ meanOpts[which(. == min(.[. >= quantile(., probs =c(0.95)) ]) )],
                                           high_ci = ~ meanOpts[which(. == max(.[. >= quantile(., probs =c(0.95)) ]) )] 
                                           #top = ~ meanOpts[which.min(abs(. - (max(.) + 1.96*sd(.)) ))],
                                           #bot = ~ meanOpts[which.min(abs(. - (max(.) - 1.96*sd(.)) ))]
  )) ) %>%
  pivot_longer(cols = everything())

######################################################################################################################
# Subsample across weekdays ####

daily_data_ww <- clean_data %>%
  mutate(day = wday(date)) # 1 is Monday

sample_methods = list('MTWTFSS' = c(1, 2, 3, 4, 5, 6, 7),
                      'MTWTF' = c(1, 2, 3, 4, 5), 'MWF' = c(1, 3, 5),
                      'TTS' = c(2, 4, 6), 'MT' = c(1, 4), 'TF' = c(2, 5), 
                      'WS' = c(3, 6),
                      'M' = c(1), 'W' = c(3), 'F' = c(5))

sampled_Restimates <- data.frame()

for (s_ind in 1:length(sample_methods)){
  subsampled_data <- daily_data_ww %>%
    filter(day %in% sample_methods[[s_ind]]) %>%
    mutate(orig_data = TRUE) %>%
    complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
    mutate(across(c(n_gene, s_gene, ORF1a), ~ zoo::na.approx(.x, na.rm = F) ))
  
  for (incidence_var in c('s_gene')){
    deconv_config = try(deconvolveIncidence(subsampled_data, incidence_var,
                                            getCountParams('incubation'), 
                                            getCountParams('benefield'),
                                            smooth_param = TRUE, n_boot = 50))
    
    Re_config = getReBootstrap(deconv_config) %>%
      mutate(sampling = names(sample_methods)[s_ind])
    
    write_csv(Re_config, paste0('../subsampling/Cali_', names(sample_methods)[s_ind], '_',
                                incidence_var, '.csv'))
    
    sampled_Restimates <- bind_rows(sampled_Restimates, Re_config)
  }
}

Cali_sampled_Restimates <- sampled_Restimates %>%
  mutate(category = str_length(sampling),
         sampling = factor(sampling, levels = c('M', 'W', 'F', 'MT', 'TF', 'WS', 
                                                'MWF', 'TTS', 'MTWTF', 'MTWTFSS')))  %>% 
  filter(data_type == 'infection_s_gene')

Cali_subsample_plot <- ggplot(Cali_sampled_Restimates) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, fill = sampling),
              alpha = 0.4, show.legend = F) +
  geom_line(aes(x = date, y = median_R_mean, colour = sampling), 
            alpha = 1, show.legend = T) +
  geom_hline(yintercept = 1) +
  facet_grid(rows = vars(category)) +
  coord_cartesian(ylim = c(0, 2)) +  
  scale_colour_manual(values = viridis(10)[c(1,4,9,2,5,10,3,8,7,6)]) + 
  scale_fill_manual(values = viridis(10)[c(1,4,9,2,5,10,3,8,7,6)]) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Sampling Type', fill = 'Sampling Type') +
  guides(color = guide_legend(override.aes = list(size=5))) 
Cali_subsample_plot

ggsave(paste0(plot_dir, '/', 'Re_subsampling_Cali_s.png'), height = 12, width = 14)

### Add Zurich-based subsampling ######
# A bit hacky - requires opening the CHE script and running the plot there

# ZH_subsample_plot + Cali_subsample_plot +
#   plot_layout(guides = 'collect') +
#   plot_annotation(tag_levels = 'A') &
#   theme(legend.position = 'bottom')

all_sampled_Restimates = bind_rows(ZH_sampled_Restimates, Cali_sampled_Restimates) %>%
  mutate(region = factor(recode(region,
                         "Santa Clara" = "San Jose",
                         "ZH" = "Zurich"), levels = c("Zurich", "San Jose")) )

ggplot(all_sampled_Restimates) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, fill = sampling),
              alpha = 0.4, show.legend = F) +
  geom_line(aes(x = date, y = median_R_mean, colour = sampling), 
            alpha = 1, show.legend = T) +
  geom_hline(yintercept = 1) +
  facet_grid(rows = vars(category), cols = vars(region), scale = 'free_x') +
  coord_cartesian(ylim = c(0, 2)) +  
  #scale_x_date(date_breaks = "4 weeks") +
  scale_colour_manual(values = viridis(10)[c(1,4,9,2,5,10,3,8,7,6)]) + 
  scale_fill_manual(values = viridis(10)[c(1,4,9,2,5,10,3,8,7,6)]) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Sampling Type', fill = 'Sampling Type') +
  guides(color = guide_legend(override.aes = list(size=5))) + 
  theme(legend.position = 'bottom')

ggsave(paste0(plot_dir, '/', 'Re_subsampling_Cali_ZH.png'), height = 14, width = 12)

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








