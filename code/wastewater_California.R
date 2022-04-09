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
#urlfile <- "https://data.chhs.ca.gov/dataset/f333528b-4d38-4814-bebb-12db1f10f535/resource/046cdd2b-31e5-4d34-9ed3-b48cdbc4be7a/download/covid19cases_test.csv"

#cali_case_data <- read_csv(urlfile)
cali_case_data <- read_csv('../data/covid19cases_test.csv')
county_list = "Santa Clara"

## Clean data 
first_clean_cali_case_data <- cali_case_data %>%
  filter(area %in% county_list,
         !is.na(date)) %>%
  mutate(confirmed = ifelse(reported_cases <0, 0, reported_cases),
  #mutate(confirmed = ifelse(reported_cases <0, 1000, reported_cases), # For Fig. S5
         county = area) %>%
  mutate(orig_data = TRUE) %>%
  select(date, county, confirmed, deaths, orig_data) %>%
  filter(date >= as_date("2020-11-01"),
         date <= as_date("2021-04-15"))

clean_cali_case_data <- first_clean_cali_case_data %>%
  mutate(location = recode(county,
                           "Santa Clara" = "SanJose") )


## Add testing-adjusted cases (test positivity) ####
mean_tests <- first_clean_cali_case_data %>%
  pull(total_tests) %>%
  mean(., na.rm = T)
# 14960.25

# Note: positive and total tests are by estimated
# testing date, not reporting date!
cali_testpos_data <- first_clean_cali_case_data %>%
  mutate(test_pos = (positive_tests / total_tests)*mean_tests) 

sub_testpos <- cali_testpos_data %>%
  select(date, test_pos)

clean_cali_case_data <- clean_cali_case_data %>%
  left_join(sub_testpos, by = "date") 

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
#ggsave(plot = case_data_plot, paste0(plot_dir, '/', 'SJ_case_data.pdf'), height = 11, width = 9)

test_case_plot <- ggplot(clean_cali_case_data %>% 
                           pivot_longer(cols = c('test_pos'))) +
  geom_bar(aes(x=date, y= value, fill = name), alpha = 0.5,
           position = 'identity', stat = 'identity', show.legend = F) +
  labs(x = 'Date' , y='Testing-adjusted \ncases per day') +
  scale_x_date(limits = c(as_date('2020-11-15'), as_date('2021-03-18')) ) +
  scale_fill_manual(values = "#A72B02") + 
  labs(colour = 'Variable') 

test_case_plot

## Confirmed case Rt #####

deconv_result_c <- data.frame()
result_c <- data.frame()
for (location_i in unique(clean_cali_case_data$county)){
  for (incidence_var_i in c('confirmed', 'test_pos')){
    select_data <- clean_cali_case_data %>%
      filter(county == location_i,
             !is.na(get(incidence_var_i)) ) %>%
      mutate(region = county) #%>%
    
    if (incidence_var_i  == 'confirmed'){
      delayParams = getCountParams('confirmed_cali')
    } else if (incidence_var_i == 'test_pos'){
      delayParams = getCountParams('confirmed_cali')
      #delayParams = getCountParams('zero') # For Fig. S4
    } else if (incidence_var_i == 'deaths'){
      delayParams = getCountParams('death')
    }
    
    new_deconv_result = deconvolveIncidence(select_data,
                                            incidence_var = incidence_var_i,
                                            getCountParams('incubation'),
                                            delayParams,
                                            smooth_param = TRUE, n_boot = 50) %>% #n_boot = 1000 for the paper
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
## Wastewater data - SCAN Pilot Project ####

csv_path = '../data/CA_WW_data.csv'

raw_data <- read.csv(csv_path) %>%
    mutate(n_orig = n_gene,
           s_orig = s_gene,
           orf1a_orig = ORF1a)

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
  select(date, orig_data, location, N = n_orig, S = s_orig, ORF1a = orf1a_orig) %>%
  pivot_longer(cols = c(N, S, ORF1a)) %>%
  mutate(name_orig = ifelse(!is.na(orig_data), name, 'Imputed'))

SCANpilot_data_plot <- ggplot() +
  geom_point(data = plot_raw_ww_data, aes(x=date, y= value, colour = name_orig, shape = name_orig),
             size = 2, show.legend = F) +
  geom_line(data = plot_raw_ww_data %>% filter(orig_data), 
            aes(x=date, y= value, colour = name), linetype = 'dotted', show.legend = F) +
  labs(x = 'Date' , y='Gene copies per \ngr dry weight') +
  scale_x_date(limits = c(as_date('2020-11-15'), as_date('2021-03-18')) 
  ) +
  scale_colour_manual(values = c(viridis(4)[2:4], 'lightgrey'), 
                      labels = c('N', 'S', 'ORF1a', 'Imputed'),
                      breaks = c('N', 'S', 'ORF1a', 'Imputed'),
                      name = 'Variable') + 
  labs(colour = 'Variable')

SCANpilot_data_plot
#ggsave(plot = SCANpilot_data_plot, paste0(plot_dir, '/', 'SJ_wastewater_data.pdf'), height = 11, width = 9)

###### Deconvolve #####

deconv_result <- data.frame()
result <- data.frame()
for (incidence_var_i in c('n_gene', 's_gene', 'ORF1a')){
  new_deconv_result = deconvolveIncidence(clean_data, 
                                          incidence_var = incidence_var_i,
                                          getCountParams('incubation'), 
                                          getCountParams('benefield'),
                                          smooth_param = TRUE, n_boot = 50) %>% #n_boot = 1000 in paper
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
         date <= as_date("2021-03-08"),
         #data_type != 'test_pos'
  ) %>%
  mutate(plot_row = ifelse(data_type %in% c('confirmed', 'test_pos'), data_type, 'Wastewater'),
         plot_row = factor(plot_row, levels = c('Wastewater', 'confirmed', 'test_pos')))

mean_deconv_data <- all_deconv_results %>%
  group_by(date, region, country, data_type, plot_row, location) %>%
  summarise(sd = sd(value),
            value = mean(value),
            .groups = 'drop')

SCANpilot_deconv_plot <- ggplot() +
  geom_ribbon(data = mean_deconv_data %>% filter(plot_row != 'Wastewater'), 
              aes(x=date, ymin = value -sd,  ymax = value +sd, fill = data_type),
              alpha = 0.5, show.legend = F) +
  ggpattern::geom_ribbon_pattern(data = mean_deconv_data %>% filter(plot_row == 'Wastewater'), 
                                 aes(x = date, ymin = value -sd,  ymax = value +sd, 
                                     colour = data_type, fill = data_type,
                                     pattern_colour = data_type, 
                                     pattern = 'stripe', pattern_angle = data_type), 
                                 pattern_spacing = 0.1, pattern_density = 0.01, 
                                 alpha = 0.2, show.legend = F) +
  labs(x = 'Date' , y='Estimated infection incidence per day') +
  scale_x_date(limits = c(as_date('2020-11-15'), as_date('2021-03-18')) ) +
  facet_wrap(vars(plot_row), ncol = 1, scale = 'free_y',
             labeller = labeller(plot_row = setNames(c('Cases', 'Testing-adjusted cases', 'Wastewater'),
                                                     c('confirmed', 'test_pos', 'Wastewater')) )) +
  scale_fill_manual(values = c(viridis(4), "#A72B02"),
                    breaks = c("confirmed",
                               "n_gene", "s_gene", 'ORF1a', 'test_pos'),
                    labels = c("Confirmed cases", "N", "S", 'ORF1a', "Testing-adjusted cases"),
                    name = 'Variable') +
  scale_colour_manual(values = c(viridis(4), "#A72B02"),
                      breaks = c("confirmed",
                                 "n_gene", "s_gene", 'ORF1a', 'test_pos'),
                      labels = c("Confirmed cases", "N", "S", 'ORF1a', "Testing-adjusted cases"),
                      name = 'Variable') +
  ggpattern::scale_pattern_angle_manual(values = c(n_gene = 180, s_gene = 135, ORF1a = 45)) +
  ggpattern::scale_pattern_colour_manual(values = c(n_gene = viridis(4)[2], 
                                                    s_gene = viridis(4)[3], ORF1a = viridis(4)[4]))

SCANpilot_deconv_plot
#ggsave(paste0(plot_dir, '/', 'SJ_deconv_data.pdf'), height = 10, width = 15)


#### R estimate Plots ########
all_results <- result_c %>%
  mutate(location = 'SanJose')%>%
  filter(date >= as_date("2020-11-15"),
         date <= as_date("2021-03-08"),
         #data_type == 'confirmed'
  ) %>%
  bind_rows(result) %>%
  mutate(plot_row = ifelse(data_type %in% c('confirmed', 'test_pos'), data_type, 'confirmed')) %>%
  bind_rows(result) %>%
  mutate(plot_row = ifelse(is.na(plot_row), 'test_pos', plot_row),
         plot_row = factor(plot_row, levels = c('confirmed', 'test_pos')))
# the deconvolution is 10 days prior to last data

SCANpilot_Re_plot <- ggplot(all_results) +
  geom_ribbon(data = all_results %>% filter(variable %in% c('confirmed', 'test_pos')),
              aes(x = date, ymin = median_R_lowHPD,
                  ymax = median_R_highHPD,  fill = variable), alpha = 0.5, show.legend = F) +
  ggpattern::geom_ribbon_pattern(data = all_results %>% filter(!variable %in% c('confirmed', 'test_pos')),
                                 aes(x = date, ymin = median_R_lowHPD,
                                     ymax = median_R_highHPD, colour = variable, 
                                     fill = variable,
                                     pattern_colour = variable, 
                                     pattern = 'stripe', pattern_angle = variable), 
                                 pattern_spacing = 0.1, pattern_density = 0.01, 
                                 alpha = 0.2, show.legend = F) +
  geom_line(aes(x = date, y = median_R_mean, colour = variable), size = 1.2, show.legend = T) +
  facet_wrap(vars(plot_row), ncol = 1, scale = 'free_y',
             labeller = labeller(plot_row = setNames(c('Cases', 'Testing-adjusted cases', 'Wastewater'),
                                                     c('confirmed', 'test_pos', 'Wastewater')) )) +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim = c(0.25, 1.75)) +
  scale_x_date(limits = c(as_date('2020-11-15'), as_date('2021-03-18')) ) +
  scale_colour_manual(values = c(viridis(4), "#A72B02"),
                      breaks = c("confirmed",
                                 "n_gene", "s_gene", 'ORF1a', "test_pos"),
                      labels = c("Confirmed cases", "N", "S", 'ORF1a', "Testing-adjusted cases"),
                      name = 'Variable') +
  scale_fill_manual(values = c(viridis(4), "#A72B02"),
                    breaks = c("confirmed",
                               "n_gene", "s_gene", 'ORF1a', 'test_pos'),
                    labels = c("Confirmed cases", "N", "S", 'ORF1a', "Testing-adjusted cases"),
                    name = 'Variable') +
  ggpattern::scale_pattern_angle_manual(values = c(n_gene = 180, s_gene = 135, ORF1a = 45)) +
  ggpattern::scale_pattern_colour_manual(values = c(n_gene = viridis(4)[2], 
                                                    s_gene = viridis(4)[3], ORF1a = viridis(4)[4])) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Observation Type', fill = 'Observation Type') +
  guides(color = guide_legend(override.aes = list(size=5)), pattern_angle = "none", pattern_colour = "none") +
  theme(
    legend.position = 'bottom',
    strip.text = element_blank()
  )

SCANpilot_Re_plot 
#ggsave(plot = SCANpilot_Re_plot, paste0(plot_dir, '/', 'SJ_wastewater_Re.pdf'), height = 14, width = 12)

###########################################################
# Fig 2 for the paper

((SCANpilot_data_plot / case_data_plot / test_case_plot ) | SCANpilot_deconv_plot) / SCANpilot_Re_plot + 
  #plot_layout(guides = "collect") + 
  plot_layout(height = c(3,2), guides = "collect") + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 30),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.key.size = unit(3,"line") )
ggsave(paste0(plot_dir, '/', 'Fig_Cali_SCANpilot.png'), height = 18, width = 14)
ggsave(paste0(plot_dir, '/', 'Fig2.pdf'), height = 18, width = 14)

## with Dec 30 filled with 1500
#ggsave(paste0(plot_dir, '/', 'Fig_Cali_SCANpilot_filled.png'), height = 18, width = 14)

## with symptom-onset to test delay = 0
#ggsave(paste0(plot_dir, '/', 'Fig_Cali_SCANpilot_notestdelay.png'), height = 18, width = 14)

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
# Takes long to compute, results included in scan folder

# for(incubationParam in c('zero', 'incubation')){
#   for (incidence_var in c('s_gene')){
#     
#     meanOpts = seq(0.5, 15, 0.5)
#     sdOpts = seq(0.5, 10, 0.5)
#     
#     deconv_results = cbind(expand_grid(meanOpts, sdOpts), 
#                            'rmse' = NA, 'coverage' = NA, 'mape' = NA)
#     
#     for (row_id in 1:nrow(deconv_results)){
#       deconv_config = try(deconvolveIncidence(clean_data, incidence_var,
#                                               getCountParams(incubationParam), 
#                                               getGammaParams(deconv_results[row_id, 'meanOpts'],
#                                                              deconv_results[row_id, 'sdOpts']),
#                                               smooth_param = TRUE, n_boot = 50))
#       
#       if('try-error' %in% class(deconv_config)){
#         deconv_results[row_id, c('rmse', 'coverage', 'mape')] = c(Inf, 0, Inf)
#         next
#       }
#       
#       Re_config = getReBootstrap(deconv_config)
#       
#       deconv_results[row_id, c('rmse', 'coverage', 'mape')] = compareTraces(Re_config, result_c)
#     }
#     
#     write_csv(deconv_results, paste0('../scan/deconv_Cali_', incubationParam, '_', 
#                                      incidence_var, '_SCdata.csv'))
#     
#   }
# }

## Results of scan ####
scan_options <- expand.grid(incubationParam = c('incubation', 'zero'),
                            incidence_var = c('s_gene'))

for (row_i in 1:nrow(scan_options)){
  deconv_results <- read_csv(paste0('../scan/deconv_Cali_', 
                                    scan_options[row_i, 'incubationParam'], '_', 
                                    scan_options[row_i, 'incidence_var'], '.csv'),
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
                  scan_options[row_i, 'incubationParam'], '.png'), height = 10, width = 14)
  }
  
  
  print('##################')
  print(paste(scan_options[row_i, 'incidence_var'], ' ', scan_options[row_i, 'incubationParam']))
  print(deconv_results[which.min(deconv_results$rmse), ])
  print(deconv_results[which.max(deconv_results$coverage), ])
  print(deconv_results[which.min(deconv_results$mape), ])
}

###########################################################
# Plot scan together with the one for Zurich, Fig 3 paper ####

#scanmode = 'incubation'
scanmode = 'zero'

ZH_deconv_results <- read_csv(paste0('../scan/deconv_', 
                                     scanmode, '_', 
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
                                       scanmode, '_', 
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


## Set method here
combined_plot <- ggplot(all_deconv_results) +
  geom_tile(aes(x = meanOpts, y = sdOpts, fill = rmse )) +
  geom_contour(aes(x = meanOpts, y = sdOpts, z = rmse ), colour = 'black', 
               binwidth = 0.1*min(all_deconv_results$rmse)) +
  facet_wrap(vars(city)) +
  labs(x = 'Mean', y = 'Standard deviation', fill = 'RMSE') +
  scale_fill_viridis(direction = -1)
combined_plot

ggsave(paste0(plot_dir, '/', 'Scan_Cali_ZH_RMSE_', scanmode, '.png'), height = 7, width = 20)


combined_plot <- ggplot(all_deconv_results) +
  geom_tile(aes(x = meanOpts, y = sdOpts, fill = coverage )) +
  geom_contour(aes(x = meanOpts, y = sdOpts, z = coverage ), colour = 'black', 
               binwidth = 0.1*max(all_deconv_results$coverage)) +
  facet_wrap(vars(city)) +
  labs(x = 'Mean', y = 'Standard deviation', fill = 'Coverage') +
  scale_fill_viridis(direction = 1)
combined_plot

ggsave(paste0(plot_dir, '/', 'Scan_Cali_ZH_COV_', scanmode, '.png'), height = 7, width = 20)


combined_plot <- ggplot(all_deconv_results) +
  geom_tile(aes(x = meanOpts, y = sdOpts, fill = mape )) +
  geom_contour(aes(x = meanOpts, y = sdOpts, z = mape ), colour = 'black', 
               binwidth = 0.1*min(all_deconv_results$mape)) +
  facet_wrap(vars(city)) +
  labs(x = 'Mean', y = 'Standard deviation', fill = 'MAPE') +
  scale_fill_viridis(direction = -1)
combined_plot

ggsave(paste0(plot_dir, '/', 'Scan_Cali_ZH_MAPE_', scanmode, '.png'), height = 7, width = 20)

###########

min_method = min(all_deconv_results$mape)
#min_method = max(all_deconv_results$coverage)
error_tol = min_method*0.1

all_deconv_results %>%
  group_by(city) %>%
  mutate(inv_coverage = - coverage) %>%
  summarise(across(c(inv_coverage), .fns = list(max_mean = ~ meanOpts[which.min(.)],
                                                max_sd = ~ sdOpts[which.min(.)],
                                                bot_mean = ~ range(meanOpts[(. <= min(.)+error_tol)])[1],
                                                top_mean = ~ range(meanOpts[(. <= min(.)+error_tol)])[2],
                                                bot_sd = ~ range(sdOpts[(. <= min(.)+error_tol)])[1],
                                                top_sd = ~ range(sdOpts[(. <= min(.)+error_tol)])[2]
  ) ))


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

##### New revised plot with separate legends #####

plot1 <- ggplot(all_sampled_Restimates %>% filter(category %in% c(1,7))) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, fill = sampling),
              alpha = 0.4, show.legend = F) +
  geom_line(aes(x = date, y = median_R_mean, colour = sampling), 
            alpha = 1, show.legend = T) +
  geom_hline(yintercept = 1) +
  facet_wrap(vars(region), scale = 'free_x') +
  coord_cartesian(ylim = c(0, 2)) +  
  scale_colour_manual(values = viridis(10)[c(1,4,9,6)]) + 
  scale_fill_manual(values = viridis(10)[c(1,4,9,6)]) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Sampling \nScheme', fill = 'Sampling \nScheme') +
  guides(color = guide_legend(override.aes = list(size=5))) + 
  theme(legend.position = 'right',
        panel.spacing = unit(2, 'lines'),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

plot2 <- ggplot(all_sampled_Restimates %>% filter(category %in% c(2,7))) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, fill = sampling),
              alpha = 0.4, show.legend = F) +
  geom_line(aes(x = date, y = median_R_mean, colour = sampling), 
            alpha = 1, show.legend = T) +
  geom_hline(yintercept = 1) +
  facet_wrap(vars(region), scale = 'free_x') +
  coord_cartesian(ylim = c(0, 2)) +  
  scale_colour_manual(values = viridis(10)[c(2,5,10,6)]) + 
  scale_fill_manual(values = viridis(10)[c(2,5,10,6)]) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Sampling Type', fill = 'Sampling Type') +
  guides(color = guide_legend(override.aes = list(size=5))) + 
  theme(legend.position = 'right',
        panel.spacing = unit(2, 'lines'),
        strip.text = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

plot3 <- ggplot(all_sampled_Restimates %>% filter(category %in% c(3,7))) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, fill = sampling),
              alpha = 0.4, show.legend = F) +
  geom_line(aes(x = date, y = median_R_mean, colour = sampling), 
            alpha = 1, show.legend = T) +
  geom_hline(yintercept = 1) +
  facet_wrap(vars(region), scale = 'free_x') +
  coord_cartesian(ylim = c(0, 2)) +  
  scale_colour_manual(values = viridis(10)[c(3,8,6)]) + 
  scale_fill_manual(values = viridis(10)[c(3,8,6)]) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Sampling Type', fill = 'Sampling Type') +
  guides(color = guide_legend(override.aes = list(size=5))) + 
  theme(legend.position = 'right',
        panel.spacing = unit(2, 'lines'),
        strip.text = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

plot5 <- ggplot(all_sampled_Restimates %>% filter(category %in% c(5,7))) +
  geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                  ymax = median_R_highHPD, fill = sampling),
              alpha = 0.4, show.legend = F) +
  geom_line(aes(x = date, y = median_R_mean, colour = sampling), 
            alpha = 1, show.legend = T) +
  geom_hline(yintercept = 1) +
  facet_wrap(vars(region), scale = 'free_x') +
  coord_cartesian(ylim = c(0, 2)) +  
  scale_colour_manual(values = c("#A72B02", viridis(10)[6]) ) + 
  scale_fill_manual(values = c("#A72B02", viridis(10)[6]) ) +
  labs( x = 'Date', y = 'Estimated Re', 
        colour = 'Sampling Type', fill = 'Sampling Type') +
  guides(color = guide_legend(override.aes = list(size=5))) + 
  theme(legend.position = 'right',
        panel.spacing = unit(2, 'lines'),
        strip.text = element_blank(),
        legend.title = element_blank())

plot1 + plot2 + plot3 + plot5 +
  plot_layout(ncol = 1)

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








