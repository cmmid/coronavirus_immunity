# R sensitivity of seasonal timing

sourceCpp("figures", "SEIR_365.cpp")
sourceCpp("figures","SEIR_lockdown_seasonal_365.cpp")

RBINOM <- plot_rbinom(samples = 100, 
                      trace_period = c(1:dim(trace_to_sample)[1]),
                      trace_dt = trace_to_sample,
                      model_type = "SEIR")


tiff(here("figures","rbinom_365.tiff"), height = 2000, width = 3200, res = 300)
RBINOM

# optimise the phi parameter


optim_ll_seasonal <- function(input_guess){
  init_theta <- c(seasonal_sample["waning_day"],
                  seasonal_sample["seasonal_R0"],
                  seasonal_sample["seasonal_reported_1"],
                  seasonal_sample["seasonal_reported_2"],
                  seasonal_sample["seasonal_reported_3"],
                  seasonal_sample["seasonal_reported_5"],
                  seasonal_sample["seasonal_amplitude"],
                  phi = input_guess)
  
  ll <- get_log_posterior_no_covid(init_theta, model_type = "SEIR",bb = "binomial")
  
  return(-ll)}


optim(-15, 
      fn = optim_ll_seasonal, 
      method = "Brent", 
      lower = -364,
      upper = 364)


##### Then use in the extra plots and projections ####
