# create parameters extra function
create_parameters <- function(parameter_guesses){
  num_groups <- length(pop_params_base[["age_groups"]])
  
  rep_f_1 <- exp(as.numeric(parameter_guesses["seasonal_reported_1"]))/
    (1+ exp(as.numeric(parameter_guesses["seasonal_reported_1"])))
  rep_f_2 <- exp(as.numeric(parameter_guesses["seasonal_reported_2"]))/
    (1+ exp(as.numeric(parameter_guesses["seasonal_reported_2"])))
  rep_f_3 <- exp(as.numeric(parameter_guesses["seasonal_reported_3"]))/
    (1+ exp(as.numeric(parameter_guesses["seasonal_reported_3"])))
  rep_f_4 <- exp(as.numeric(parameter_guesses["seasonal_reported_2"]))/
    (1+ exp(as.numeric(parameter_guesses["seasonal_reported_2"])))
  rep_f_5 <- exp(as.numeric(parameter_guesses["seasonal_reported_5"]))/
    (1+ exp(as.numeric(parameter_guesses["seasonal_reported_5"])))
  
  parameters <- list(
    beta_covid_0 = 0,
    beta_covid_1 = 0,
    beta_other = 0,
    waning_covid = 0,
    waning_other = 1/as.numeric(fixed_parameters["waning_day"]),
    gamma_covid = 0,
    gamma_other = 1/as.numeric(fixed_parameters["seasonal_gamma"]),
    incubation_covid = 0,
    incubation_other = 1/as.numeric(fixed_parameters["seasonal_incubation"]),
    num_grps = num_groups,
    births = pop_params_base[["births"]],
    deaths =  pop_params_base[["deaths"]],
    red_susc_covid = rep(0,num_groups),
    contacts_all = pop_params_base$contacts_all,
    contacts_hh = pop_params_base$contacts_hh,
    age_rates = pop_params_base[["ageing"]],
    covid_imm = rep(0,num_groups), 
    reporting_delay_covid = 0,
    reporting_delay_other = (1/as.numeric(fixed_parameters["s_symp_to_report"])), 
    inf_to_symp_covid = 0,
    inf_to_symp_other = (1/as.numeric(fixed_parameters["s_inf_to_symp"])),
    contacts_other = pop_params_base$contacts_other, 
    contacts_school = pop_params_base$contacts_school*1,
    covid_intros = rep(0,num_groups), 
    covid_deaths = rep(0,5), 
    # seasonal_reported = c(as.numeric(parameter_guesses["seasonal_reported_1"]),
    #                       as.numeric(parameter_guesses["seasonal_reported_2"]),
    #                       as.numeric(parameter_guesses["seasonal_reported_3"]),
    #                       as.numeric(parameter_guesses["seasonal_reported_2"]),
    #                       as.numeric(parameter_guesses["seasonal_reported_5"])), 
    seasonal_reported = c(rep_f_1,
                          rep_f_2,
                          rep_f_3,
                          rep_f_2,
                          rep_f_5), 
    start_dying_off = 14, 
    child_extra_reduc = 1
  )
  
  parameters <- c(parameters,
                  sig_R = convert_sig(as.numeric(fixed_parameters["sig_R"])),
                  sig_I = convert_sig(0),
                  sig_solid_R = as.numeric(fixed_parameters["sig_solid_R"]),
                  sig_solid_I = 0,
                  rho = rho,
                  covid_date=unname(as.character(as.Date("3000-01-01"))),
                  covid_time=Inf,
                  covid_change_time = Inf,
                  amplitude = NA, 
                  phi = as.numeric(parameter_guesses["phi"]), 
                  mobility_start = Inf, 
                  contacts_reduc = unname((data.frame(pop_params_base$contacts_reduc))), 
                  social_distancing = 0)
  
  beta_val_seasonal <- calc_beta_SEIPRR(start_val = 0.025,
                                        requiredR0 = fixed_parameters[["seasonal_R0"]],
                                        parameters = parameters,
                                        pop_params = pop_params_base, 
                                        virus = "other",
                                        timestep = 1)
  
  parameters$beta_other <- beta_val_seasonal$par
  
  parameters$amplitude <- as.numeric(fixed_parameters["seasonal_amplitude"])/
    fixed_parameters[["seasonal_R0"]]*beta_val_seasonal$par
  
  
  return(parameters)
}
