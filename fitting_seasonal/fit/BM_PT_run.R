################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################

# runnning script for paralell tempering 
########## Load in params and functions and initialise paralisation

# Specify the characteristics of the cluster
cores_to_use <- detectCores()
cl <- parallel::makeCluster(cores_to_use, setup_strategy = "sequential")
registerDoParallel(cl)
array_num <- 1
# setup the initial conditions

lower_bounds <- c(waning_day = 100, 
                    seasonal_R0 = 1,
                    seasonal_reported_1 = -Inf, 
                    seasonal_reported_2 = -Inf, 
                    seasonal_reported_3 = -Inf, 
                    seasonal_reported_5 = -Inf, 
                    seasonal_amplitude = 0, 
                    phi = -364)

upper_bounds <-  c(waning_day = 3000, 
                   seasonal_R0 = 8.5,
                   seasonal_reported_1 = Inf, 
                   seasonal_reported_2 = Inf, 
                   seasonal_reported_3 = Inf, 
                   seasonal_reported_5 = Inf, 
                   seasonal_amplitude = 2.5, 
                   phi = 364)

init_theta <- c(waning_day = 200,
                seasonal_R0 = 2,
                seasonal_reported_1 = -8,
                seasonal_reported_2 = -8,
                seasonal_reported_3 = -8,
                seasonal_reported_5 = -8,
                seasonal_amplitude = 1.1,
                phi = 20)


names(init_theta) <- c("waning_day", "seasonal_R0", "seasonal_reported_1", 
                       "seasonal_reported_2", "seasonal_reported_3", 
                       "seasonal_reported_5", "seasonal_amplitude", 
                       "phi")

Temperatures <- c(1,2,3,6.5,12,25,50,70,100,200,300,400,500,600,800,1000)
n_chains <- 16


covmat <- as.matrix(read.csv(here::here("fitting_seasonal/fit","covmat_3.csv"), row.names = 1))/2
rownames(covmat) <- colnames(covmat) <- names(init_theta)
proposal_sd <- diag(covmat)

###### Can include this in order to start from the end of a previous chain
#   load("SEIR_PT_newlims_5_2021-02-16.Rdata")
#  new_start <- lapply(1:n_chains, function (x)
#    as.list(tail(total_trace[[x]],1)[,2:(length(proposal_sd)+1)]))
#  for ( section in 1:n_chains){
#    names(new_start[[section]]) <- names(proposal_sd)}
# init_theta <- new_start


# Export parameters to the c
clusterExport(cl, list("n_chains","fixed_parameters",
                        "age_groups", "pop_numbers", "pop_params_base", "Temperatures",
                       "covmat", "proposal_sd", "init_theta",
                       "lower_bounds", "upper_bounds", "contacts_all", "rho",
                       "length_to_run", "lp_15", "all_oths", "seasonal_dates_15", 
                       "seasonal_15_20", "init_state_temp"))
# Export required functions to the cores
clusterExport(cl, list("get_log_posterior_no_covid",
                       "calc_beta_SEIPRR", "optim_beta_SEIPRR", "calc_R0_SEIPRR",
                       "calculate_init_state_reported2", "run_model_seasonal", "calculate_init_state_reported2",
                       "naming_states","summary_stats_reported_seasonal","reject_seasonal_Tech",
                       "calc_pop_per_age", "summary_groups","calc_lik_seasonal_ages_binomial",
                       "temp_func", "step_fn", "swapping_rates", "convert_sig",
                       "adjust_unsymmetrical", "infected_year_age", "swapping_props", "swapping_acceps", 
                       "run_model_get_logliks_seasonalonly", "create_parameters", "add_on", "get_llprior"
                       ))



clusterEvalQ(cl, library("Rcpp"))
clusterEvalQ(cl, library("data.table"))
clusterEvalQ(cl, library("deSolve"))
clusterEvalQ(cl, library("rootSolve"))
clusterEvalQ(cl, library("MASS"))
clusterEvalQ(cl, library("tmvtnorm"))
clusterEvalQ(cl, library("RcppSims"))


time_start <- Sys.time()
name <- "SEIR_PT_test"

total_trace <- trace_wrapper(length_run = 10, # length of total run
                             each_run = 10, # lenght of each subset
                             n_chains = n_chains, # number of chians
                             init_theta = init_theta, # initial parameter guesses
                             indep = 5, # number steps to run chains independetly
                             covmat = covmat, # covariance matrix
                             adapt_rate = 0.55, # temperature adaption cooling
                             adaption_starter = 1000000, # time point to start temp adapt
                             proposal_sd = proposal_sd, # proposal s.d.
                             Temperatures = Temperatures,  # Temperatures
                             name = name, # name to save file
                             virus = "seasonal",  #only works for seasonal
                             model_type = "SEIR", # currently only works for SEIR
                             bb = "binomial", # currently only works for binomial
                             cont = F # whether or not to continue from previous chain
)

time_end <- Sys.time()
print(time_end - time_start)
save(total_trace, file=here("fitting_seasonal/fit",paste0(name,"_", array_num,"_", Sys.Date() , ".Rdata")))
stopCluster(cl)
