################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-25
################################################################################
setwd("~/Documents/GitHub/coronavirus_immunity/")
# Overall setup script

####### LOAD THE INITIAL REQUIREMENTS ######
# load the relevant packages
library(data.table)
library(MASS)
library(deSolve)
library(rootSolve)
library(Rcpp)
library(tmvtnorm)
library(ggplot2)
library(gridExtra)
library(here)
library(parallel)
library(doParallel)
library(coda)

## This is the package that contains the actual model code. It's in a package 
# for convenience when running on a cluster.
#install.packages("RcppCoronaImmunitty_0.1.0.tar.gz",
#               repos = NULL, type = "source")
library(RcppCoronaImmunitty)

# load in the functions and parameters for fitting
source(here("fitting_seasonal/fit", "BM_PT_functions.R"))
source(here("fitting_seasonal/fit", "BM_PT.R"))
source(here("fitting_seasonal/fit", "BM_PT_parameters.R"))

###### RUN THE FIT ###### 

#run the fit (from the below script. can alter preferences etc.)
#Note: this takes a long time! The end traces are already saved so you can
#carry on to the analysis if wanted
file.edit(here("fitting_seasonal/fit", "BM_PT_run.R"))

##### ANALYSE THE FIT ######

# load in the plotting functions
source(here("fitting_seasonal/analysis", "analyse_traces_funcs.R"))

#run the analysis and create the figures
# if you want to change the traces you want to analyze go into this file and change
source(here("fitting_seasonal/analysis", "analyse_traces.R"))

######## COVID SIMULATIONS ######

source(here("simulating_covid", "covid_sims_data.R"))
source(here("simulating_covid", "covid_sims_functions.R"))

# take the samples and make figures
# if any fixed parameters have changed need to change upper limit transmission
# rate by running last line of this script and then putting in as limit
source(here("simulating_covid", "Covid_sims_main.R"))


####### COVID FOREWARD PROJECTIONS ######

# source any required functions
source(here("projections", "projections_functions.R"))
# run the simulations and create graphs
source(here("projections", "projections_run.R"))

###### EXTRAS #####
# extra plots
source(here("figures", "plotting_data.R"))
source(here("figures", "extra_plots.RDS"))
# compare fit with Kissler parameters
file.edit(here("kissler_comparison","kissler_comparison.R"))








