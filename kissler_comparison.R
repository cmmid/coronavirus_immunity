################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################

# see what the model output is with Kissler estimates

# update the fixed parameters
fixed_parameters["waning_day"] <- 45*7
fixed_parameters["seasonal_R0"] <- 2
fixed_parameters["seasonal_amplitude"] <- 0.231
fixed_parameters["phi"] <- 2*7
fixed_parameters["seasonal_incubation"] <- 3

#overwrite the parameter creating function so that just
# reporting rates and timing are estimated
source("kissler_create_params_func.R")

# start dates
params_to_run <-  c(seasonal_reported_1 = -7,
                    seasonal_reported_2 = -8,
                    seasonal_reported_3 = -8,
                    seasonal_reported_5 = -9, 
                    phi = 14)
# pop it into an optim function
op_out <- optim(
  par = params_to_run, 
  upper = c(0,0,0,0,300), 
  lower = c(-30, -30,-30,-30,0),
  method = "L-BFGS-B",
  control=list(fnscale=-1), # maximise instead of minimise
  fn = get_log_posterior_no_covid, 
  bb= "binomial",
  model_type = "SEIR"
)

# run the model
params_to_run <- data.table(
  waning_duration = 45*7,
  seasonal_R0 = 2, 
  reporting_rate_1 = exp((op_out$par[1])/(1 + exp(op_out$par[1]))),
  reporting_rate_2 = exp((op_out$par[2])/(1 + exp(op_out$par[2]))),
  reporting_rate_3 = exp((op_out$par[3])/(1 + exp(op_out$par[3]))),
  reporting_rate_4 = exp((op_out$par[4])/(1 + exp(op_out$par[4]))),
  seasonal_amplitude = 0.231,
  seasonal_timing = op_out$par[5]
)
# plot the output
IND_RUN <- plot_rbinom(samples = 1, 
                       trace_period = 1,
                       trace_dt = params_to_run,
                       model_type = "SEIR")

tiff(here("figures","rbinom_kissler_params.tiff"), height = 2000, width = 3200, res = 300)
IND_RUN
dev.off()

# rename parameters for ll check
params_to_run <- data.table(
  waning_day = 45*7,
  seasonal_R0 = 2, 
  seasonal_reported_1 = op_out$par[1],
  seasonal_reported_2 = op_out$par[2],
  seasonal_reported_3 = op_out$par[3],
  seasonal_reported_5 = op_out$par[4],
  seasonal_amplitude = 0.231,
  phi = op_out$par[5]
)

#check ll
get_log_posterior_no_covid(params_to_run, bb="binomial",model_type = "SEIR")



