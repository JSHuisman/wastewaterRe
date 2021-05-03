# Wastewater-based estimation of the effective reproductive number of SARS-CoV-2 

This repository contains the code for the manuscript "Wastewater-based estimation of the effective reproductive number of SARS-CoV-2" by [Huisman et al.](https://www.medrxiv.org/content/10.1101/2021.04.29.21255961v1).

The scripts in this repository assume that the code of the Re estimation pipeline of [Huisman, Scire et al.](https://www.medrxiv.org/content/10.1101/2020.11.26.20239368v1) is available in a folder "covid-19-re-shiny-app" at the same level as this repository. The git-repo for the pipeline can be found [here](https://github.com/covid-19-Re/shiny-dailyRe).

## Code
The code is split into 3 parts: 
 - [general functions](./code/wastewater_functions.R)
 - [the analysis of the Zurich data](./code/wastewater_Zurich.R)
 - [the analysis of the San Jose data](./code/wastewater_California.R).

When run from within this repository, the code recreates the [results](./results) and [figures](./figures) for Zurich shown in the manuscript. The settled solids measurements for California will be made accessible soon, whereafter the code can be used to reproduce those figures as well.

## Data
Includes case incidence data in [Zurich](./data/ZH_case_incidence_data.csv) and [California](./data/COVID-19_case_counts_by_date.csv).

## Results
The scan across SLD parameter pairs takes quite long to compute, so we have added the results (shown in the paper) to the folder ["scan"](./scan).

The folder ["subsampling"](./subsampling) contains the Rww estimates for the data that was subsampled by weekday/sampling regime (this takes less time to compute but is included for convenience).
