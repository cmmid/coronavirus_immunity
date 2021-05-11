# Coronavirus_immunity
 
This repository contains the code for the paper: How Immunity from and interaction with seasonal coronaviruses can shape future SARS-CoV-2 epidemiology

All the fitting, analysis and figure generation can be run from "master_script.R". This loads in all the dependencies and functions. It has to be run in the order specified, and opens any files where manual edits can be made (e.g. number of simulations to run). All the code and data (which is publically available) are in the repo, and can be run directly from the master script. However please note that some elements (particularly the fitting) is time-consuming. The fitted traces are included in the repo, so that further analysis can be followed without having to run the whole fitting mechanism. 

Note that the main model code is stored within a package (RcppCoronaImmunity_0.1.0.tar.gz) for convienience of use on HPCs. The package will be loaded through a call in master_script.R, so that the functions are available within R. The model scripts for reference are in the reference_models folder.

Authors: Naomi R Waterlow, Edwin Van Leewen, Nick Davies, Stefan Flasche, Rosalind M Eggo
See LICENSE file for licensing details.

Corresponding author: Naomi R Waterlow, naomi.waterlow1@lshtm.ac.uk

Input/Output files - see paper for references: 
- all_data_ages.csv : seasonal coronavirus data
- contacts*.csv : mixing patterns based of POLYMOD
- theta_bounds.csv : parameter limits for fitting
- covmat_3.csv: covariance matrix used in fitting
- SEIR*.Rdata : output traces from fitting
- trace_to_sample.Rdata: combined posterior trace, used as input for simulating_covid
- ONS_daily_deaths.csv : Covid-19 deaths
- ONS_location_death.csv : Covid-19 deaths by location - not currently used
- uk_google_mobility.csv : google mobility data used to simulate reduction in contacts

Package versions are: 

To install the correct versions of the packages, you can use the install_packages.R script.

coda_0.19-3
doParallel_1.0.15
iterators_1.0.12
foreach_1.5.0
here_0.1
gridExtra_2.3
ggplot2_3.3.2
tmvtnorm_1.4-10
gmm_1.6-5
sandwich_2.5-1 
Matrix_1.2-18
mvtnorm_1.1-1 
Rcpp_1.0.6
rootSolve_1.8.2.1
deSolve_1.28
MASS_7.3-51.6
data.table_1.13.4
viridis_0.5.1 
forcats_0.5.0


All analysis except for the fitting was done in R version 4.0.0 on macOS Catalina 10.15.7. The fitting was done on R 3.4 on AWS ec2 machines, running a linux AMI. On AWS the package mvtnorm was version 1.0.8.

