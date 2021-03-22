################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################


# BM_PT functions

get_log_posterior_no_covid <- function(theta, bb = "no", model_type){
  
  ll <- run_model_get_logliks_seasonalonly(parameter_guesses = theta, bb, model_type)
  print(ll)
#   prior <- get_llprior(parameter_guesses = theta) 
  log_posterior <- ll #+ prior
  return(log_posterior)
}

run_model_get_logliks_seasonalonly <- function(parameter_guesses, bb = "no", 
                                               model_type){
  #combine the parameters

  parameters <- create_parameters(unlist(parameter_guesses))

  # print(parameters)
  # Run the model - seasonal corona only. For 20 years + a bit, to find low point (lp)
  output_s <- run_model_seasonal(parameters, model_type)
  colnames(output_s) <- naming_states(model_type)
  # Get the incidence of reporting
  reporting <- summary_stats_reported_seasonal(output_s, type = model_type)
  
  # reporting_selection2 <- reporting[required_set, c("time", ..all_others)
  # ][ , rowSums(.SD), by= time]
  
  # reject based on Techumseh
  likelihood_test <-0# reject_seasonal_Tech(output_s, model_type)
  
  if(likelihood_test ==0){ 
    
    # summaries
    reportin_2020 <- summary_stats_reported_seasonal(output_s, type = model_type)
    reportin_2020_daily <- summary_groups(reportin_2020)
    # calculate likelihood of the data: monthly seasonal age 
    # if(bb == "no"){
    #   lik_seasonal_ages <- calc_lik_seasonal_ages(reportin_2020_daily, parameters)
    # } else if (bb == "yes"){
    #   lik_seasonal_ages <- calc_lik_seasonal_ages_bb(reportin_2020_daily, parameters)
    # }
    if(bb == "binomial"){
      lik_seasonal_ages <- calc_lik_seasonal_ages_binomial(reportin_2020_daily, parameters)
    }
    
    likelihood_data <- lik_seasonal_ages   
  } else {likelihood_data <- likelihood_test
  
  lik_seasonal_ages <- NA ;lik_covid_deaths <- NA }

  likelihood_total <- likelihood_data #+ get_llprior
  return(c(likelihood_total) )
}

get_llprior <- function(parameter_guesses){
  
 prior <- dgamma(x = as.numeric(parameter_guesses["seasonal_R0"]), 
                                shape = 2, 
                                rate =1, 
                                log=T)
 
 return(prior)
  
}


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
    waning_other = 1/as.numeric(parameter_guesses["waning_day"]),
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
                                        requiredR0 = parameter_guesses[["seasonal_R0"]],
                                        parameters = parameters,
                                        pop_params = pop_params_base, 
                                        virus = "other",
                                        timestep = 1)

  parameters$beta_other <- beta_val_seasonal$par
  
  parameters$amplitude <- as.numeric(parameter_guesses["seasonal_amplitude"])/
    parameter_guesses[["seasonal_R0"]]*beta_val_seasonal$par
  
  
  return(parameters)
}


run_model_seasonal <- function(parameters, model_type){
  # Calculate the initial state
  init_state <- init_state_temp
  # init_state <- calculate_init_state_reported2(pop_sizes = pop_numbers,
  #                                              theta = parameters,
  #                                              covid = rep(0,length(age_groups)),
  #                                              other = rep(1,length(age_groups)),
  #                                              type = model_type)
  # times to run model
  times <- c(1:length_to_run)
  # run the model
  if(model_type == "SEIR"){
    func_to_use <- SEIR_2virus_cons
  } else if(model_type =="SEIPRR") {
    func_to_use <- SEIPRR_2virus_interact_cons
  } else if(model_type == "SEIRR"){
    func_to_use <- SEIRR_2virus_cons
  } else {
    message("Invalid model type! Stopped.")
    stop()
  }
  outall <- as.data.table(ode(
    y = init_state,
    t = times,
    func = func_to_use,
    parms = parameters,
    method = "rk4"))
  
  
  return(outall)
  
}

calculate_init_state_reported2 <- function(pop_sizes, theta,
                                           covid = 0, other = 0, 
                                           type){
  init_state_season <- c()
  
  if(type == "SIPRR"){
    # for each age group
    for(a in c(1:length(pop_sizes))) {
      
      cov_imm <- theta[["covid_imm"]][a]
      
      temp = c(
        (pop_sizes[a] -
           (pop_sizes[a] * cov_imm)-covid[a]-other[a]),#SS
        covid[a],0, #ISPS
        (pop_sizes[a] *  cov_imm), #,RS,
        other[a],#,SI
        0,0, 0,0,0,0,0, #II, PI, RI, SP, IP, PP, RP,
        0, #SR
        0,0,0,# IR, PR,RR, 
        0,0,0,0,0,0,0,0,0, # all the extra R2s
        0,0, 0,0,0,0) #Rcases I cases reported
      
      init_state_season <- c(init_state_season, temp)
    }
  } else if(type == "SEIPRR"){
    for(a in c(1:length(pop_sizes))) {
      
      cov_imm <- theta[["covid_imm"]][a]
      
      temp = c(
        (pop_sizes[a] -
           (pop_sizes[a] * cov_imm)-covid[a]-other[a]),#SS
        covid[a],0, #IS PS
        (pop_sizes[a] * cov_imm), #,RS,
        other[a],#,SI
        0,0, 0,0,0,0,0, #II, PI, RI, SP, IP, PP, RP,
        0, #SR
        0,0,0,# IR, PR,RR, 
        0,0,0,0,0,0,0,0,0, # all the extra R2s, 
        0,0,0,0,0,0,0,0,0,0,0, # all the extra Es.
        0,0,0,0,0,0) #Rcases I cases reported
      
      init_state_season <- c(init_state_season, temp)
    }}else if(type == "SEIR"){
      for(a in c(1:length(pop_sizes))) {
        
        cov_imm <- theta[["covid_imm"]][a]
        
        temp = c(
          (pop_sizes[a] -
             (pop_sizes[a] * cov_imm)-covid[a]-other[a]),#SS
          0,covid[a], #ES, IS,
          0,0,0,0, # SE, EE, IE, RE
          (pop_sizes[a] * cov_imm), #,RS,
          other[a],0,0,0,#,SI, EI, II, RI,
          0,0, 0,0, #SR, ER, IR, RR
          0,0,0,0,0,0) #Reporting compartments
        
        init_state_season <- c(init_state_season, temp)
      }}else if(type == "SEIRR"){
        for(a in c(1:length(pop_sizes))) {
          
          cov_imm <- theta[["covid_imm"]][a]
          
          temp = c(
            (pop_sizes[a] -
               (pop_sizes[a] * cov_imm)-covid[a]-other[a]),#SS
            0,covid[a], #ES, IS,
            0,0,0,0, # SE, EE, IE, RE
            (pop_sizes[a] * cov_imm), #,RS,
            other[a],0,0,0,#,SI, EI, II, RI,
            0,0, 0,0, #SR, ER, IR, RR
            0,0,0,0,0,0,0,0,0, # The extra Rs
            0,0,0,0,0,0) #Reporting compartments
          
          init_state_season <- c(init_state_season, temp)
        }
      } else {
        print("type not acceptable")
      }
  
  return(init_state_season)
}

# function to create state names
naming_states <- function(type){
  
  if (type == "SEIPRR"){
    # Specify state names
    Names <- c("time")
    for(i in 1:length(age_groups)){
      x<-c(paste0("SS", i), paste0("IS", i),paste0("PS", i),paste0("RS", i),
           paste0("SI", i),paste0("II", i),paste0("PI", i),paste0("RI", i),
           paste0("SP", i), paste0("IP", i), paste0("PP", i), paste0("RP", i),
           paste0("SR", i),paste0("IR", i), paste0("PR", i),paste0("RR", i),
           paste0("R2S", i),paste0("R2I", i), paste0("R2P", i),paste0("R2R", i),
           paste0("R2R2", i),paste0("RR2", i), paste0("PR2", i),paste0("IR2", i),
           paste0("SR2", i),
           paste0("ES", i),paste0("SE", i), paste0("EE", i),paste0("IE", i),
           paste0("PE", i),paste0("R1E", i), paste0("R2E", i),paste0("EI", i),
           paste0("EP", i),paste0("ER1", i),paste0("ER2", i),
           paste0("Covid_cases", i),
           paste0("Other_cases", i),
           paste0("Covid_symptom", i),
           paste0("Other_symptom", i),
           paste0("Covid_reported", i),
           paste0("Other_reported", i)
      )
      Names <- c(Names, x)
    }
  } else if (type == "SIPRR"){
    # Specify state names
    Names <- c("time")
    for(i in 1:length(age_groups)){
      x<-c(paste0("SS", i), paste0("IS", i),paste0("PS", i),paste0("RS", i),
           paste0("SI", i),paste0("II", i),paste0("PI", i),paste0("RI", i),
           paste0("SP", i), paste0("IP", i), paste0("PP", i), paste0("RP", i),
           paste0("SR", i),paste0("IR", i), paste0("PR", i),paste0("RR", i),
           paste0("R2S", i),paste0("R2I", i), paste0("R2P", i),paste0("R2R", i),
           paste0("R2R2", i),paste0("RR2", i), paste0("PR2", i),paste0("IR2", i),
           paste0("SR2", i),
           paste0("Covid_cases", i),
           paste0("Other_cases", i),
           paste0("Covid_symptom", i),
           paste0("Other_symptom", i),
           paste0("Covid_reported", i),
           paste0("Other_reported", i))
      Names <- c(Names, x)
    }
  } else if (type == "SEIR"){
    # Specify state names
    Names <- c("time")
    for(i in 1:length(age_groups)){
      x<-c(paste0("SS", i), paste0("ES", i),paste0("IS", i),paste0("RS", i),
           paste0("SE", i),paste0("EE", i),paste0("IE", i),paste0("RE", i),
           paste0("SI", i), paste0("EI", i), paste0("II", i), paste0("RI", i),
           paste0("SR", i),paste0("ER", i), paste0("IR", i),paste0("RR", i),
           paste0("Covid_cases", i),
           paste0("Other_cases", i),
           paste0("Covid_symptom", i),
           paste0("Other_symptom", i),
           paste0("Covid_reported", i),
           paste0("Other_reported", i))
      Names <- c(Names, x)
    }
  } else if (type == "SEIRR"){
    # Specify state names
    Names <- c("time")
    for(i in 1:length(age_groups)){
      x<-c(paste0("SS", i), paste0("ES", i),paste0("IS", i),paste0("RS", i),
           paste0("SE", i),paste0("EE", i),paste0("IE", i),paste0("RE", i),
           paste0("SI", i), paste0("EI", i), paste0("II", i), paste0("RI", i),
           paste0("SR", i),paste0("ER", i), paste0("IR", i),paste0("RR", i),
           paste0("SR2", i), paste0("ER2", i),paste0("IR2", i),paste0("RR2", i),
           paste0("R2S", i), paste0("R2E", i),paste0("R2I", i),paste0("R2R", i),
           paste0("R2R2", i),
           paste0("Covid_cases", i),
           paste0("Other_cases", i),
           paste0("Covid_symptom", i),
           paste0("Other_symptom", i),
           paste0("Covid_reported", i),
           paste0("Other_reported", i))
      Names <- c(Names, x)
    }
  } else {print("not acceptable type")}
  return(Names)
}

# summary of seasonal only
summary_stats_reported_seasonal <- function(outall, type) {
  # Specify state names
  #Calculate the incidence
  colnames(outall) <- naming_states(type)
  
  outall[, OTHER_0 := Other_reported1 - shift(Other_reported1, 1L, type = "lag")]
  outall[, OTHER_5 := Other_reported2 - shift(Other_reported2, 1L, type = "lag")]
  outall[, OTHER_10 := Other_reported3 - shift(Other_reported3, 1L, type = "lag")]
  outall[, OTHER_15 := Other_reported4 - shift(Other_reported4, 1L, type = "lag")]
  outall[, OTHER_20 := Other_reported5 - shift(Other_reported5, 1L, type = "lag")]
  outall[, OTHER_25 := Other_reported6 - shift(Other_reported6, 1L, type = "lag")]
  outall[, OTHER_30 := Other_reported7 - shift(Other_reported7, 1L, type = "lag")]
  outall[, OTHER_35 := Other_reported8 - shift(Other_reported8, 1L, type = "lag")]
  outall[, OTHER_40 := Other_reported9 - shift(Other_reported9, 1L, type = "lag")]
  outall[, OTHER_45 := Other_reported10 - shift(Other_reported10, 1L, type = "lag")]
  outall[, OTHER_50 := Other_reported11 - shift(Other_reported11, 1L, type = "lag")]
  outall[, OTHER_55 := Other_reported12 - shift(Other_reported12, 1L, type = "lag")]
  outall[, OTHER_60 := Other_reported13 - shift(Other_reported13, 1L, type = "lag")]
  outall[, OTHER_65 := Other_reported14 - shift(Other_reported14, 1L, type = "lag")]
  outall[, OTHER_70 := Other_reported15 - shift(Other_reported15, 1L, type = "lag")]
  outall[, OTHER_75 := Other_reported16 - shift(Other_reported16, 1L, type = "lag")]
  
  outall <- outall[-1,]
  
  return(outall)
}
#  reject if any agegroup annual infection is over 30%
reject_seasonal_Tech <- function(output_s, model_type){
  population <- calc_pop_per_age(output_s,timepoint=(length_to_run -364), model_type)
  population$inf_year_age <- c(infected_year_age(output_s, model_type,
                                                 year_start  = length_to_run-364))
  population[,percent := as.numeric(inf_year_age)/V1*100]
  if(any(population$percent > 30, na.rm = T)){
    likelihood = -1000000
  } else {likelihood =0}
  return(likelihood)
}

# calc population per age group and summarise into PHE groups
calc_pop_per_age <- function(outall, timepoint, model_type){
  
  if(model_type == "SEIR"){
    cols_to_exclude <- c(1)
    for(i in 1:length(age_groups)){
      cols_to_exclude <- c(cols_to_exclude, c((((i-1)*22)+18):(((i-1)*22)+23)))
    } } else if(model_type == "SEIRR"){
      cols_to_exclude <- c(1)
      for(i in 1:length(age_groups)){
        cols_to_exclude <- c(cols_to_exclude, c((((i-1)*31)+27):(((i-1)*31)+32))) # SEIRR
      } }else if (model_type == "SEIPRR"){
        cols_to_exclude <- c(1)
        for(i in 1:length(age_groups)){
          cols_to_exclude <- c(cols_to_exclude, c((((i-1)*42)+38):(((i-1)*42)+43))) # SEIRR
          
        }}
  
  
  
  named_cols_to_exclude <- naming_states(model_type)[cols_to_exclude]
  colnames(outall) <- c(naming_states(model_type))
  # calculate the population over time (have to remove tracking compartments
  suppressWarnings({
    population <- melt.data.table(outall[timepoint,!cols_to_exclude, with=F], id.vars =NULL)
  })
  for(i in 1:length(age_groups)){
    population[ grep(paste0(i,"$"), variable),agegroup := i]
  }
  
  # relabel the incorrectly labeled ones
  population[variable == "ER11", agegroup:= 1]
  population[variable == "ER12", agegroup:= 2]
  population[variable == "ER13", agegroup:= 3]
  population[variable == "ER14", agegroup:= 4]
  population[variable == "ER15", agegroup:= 5]
  population[variable == "ER16", agegroup:= 6]
  
  # calculate the sum by age group
  pop2 <- population[,sum(value),by = agegroup]
  
  pop2[agegroup == 1, PHE := "Under5"]
  pop2[agegroup == 2 | agegroup==3, PHE := "5-14 Y"]
  pop2[agegroup > 3 & agegroup <=9, PHE := "15-44 Y"]
  pop2[agegroup > 9 & agegroup <=13, PHE := "45-64 Y"]
  pop2[agegroup > 13 , PHE := "65 Y+"]
  
  return(pop2[,sum(V1), by=PHE])
}

# calculate the betabinomial likelihood for seasonal age cases monthly
calc_lik_seasonal_ages_binomial <- function(reportin_2020_daily, parameters, 
                                            covid_run = F){
  # reported seasonal cases by age over time
  reportin_seasonal <- melt.data.table(reportin_2020_daily[,c("time", ..all_oths)],
                                       id= "time")
  # get the right values
  if(covid_run == T){
    relevant_dates <- seasonal_dates
    real_data <- seasonal_19_20
  } else if (covid_run == F){
    relevant_dates <- seasonal_dates_15
    real_data <- seasonal_15_20
  }
  reportin_seasonal[, date := as.Date(time, origin = lp_15)]
  #change from daily to monthly
  for(i in 1:(length(relevant_dates)-1)){
    reportin_seasonal[ date >= relevant_dates[i] &
                         date < relevant_dates[i+1],
                       year_week := relevant_dates[i]]
  }
  
  to_match <- na.omit(reportin_seasonal[,sum(value, na.rm=T), by=c("year_week",
                                                                   "variable")
  ])
  # add the real data
  to_match[real_data, on = c("year_week", "variable"), true_value:= i.value ]
  
  
  for( stepper in 1:dim(to_match)[1]){
    if(to_match[stepper,"variable"] == "OTHER_p0"){age_set = 1
    } else if(to_match[stepper,"variable"] == "OTHER_p5"){age_set = 2
    } else if(to_match[stepper,"variable"] == "OTHER_p15"){age_set = 3
    } else if(to_match[stepper,"variable"] == "OTHER_p45"){age_set = 4
    } else if(to_match[stepper,"variable"] == "OTHER_p65"){age_set = 5}
    # work out quantile intervals
    
    to_match[stepper,"likelihood"] <- dbinom(x=as.numeric(to_match[stepper, "true_value"]),
                                             size = round(as.numeric(to_match[stepper,"V1"])),
                                             prob = parameters$seasonal_reported[age_set]
                                             , log=T)
    
  }

  # weight the off_seasons by half
  lik_summary <- to_match[, sum(likelihood, na.rm = T)]
  
  lik_out <- lik_summary
  
  if(any(is.nan(to_match$likelihood))){lik_out <- -10000}
  return(lik_out)
}

# function to get desired beta from R0 using optimisation algorithm.
calc_beta_SEIPRR <- function(start_val, requiredR0, parameters,
                             pop_params, timestep, virus){
  
  best_fit <- optim(par = start_val, 
                    fn = optim_beta_SEIPRR,
                    Ro_wanted = requiredR0, 
                    parameters = parameters,
                    pop_params = pop_params,
                    timestep = timestep,
                    virus = virus,
                    method="Brent", 
                    lower = 0, 
                    upper = 0.2)

  return(best_fit)
}

#function used to optimise the beta parameter based on a desired Ro. SEIPRR
optim_beta_SEIPRR <- function(beta_input, Ro_wanted, parameters,
                              pop_params, timestep, virus){
  R0 <- calc_R0_SEIPRR(beta_input = beta_input, parameters = parameters,
                       pop_params = pop_params, timestep = timestep,
                       virus = virus)
  difference = (abs(R0 - Ro_wanted)^2)

    return(difference)
}  

#calculate the R0 value, assuming no interaction, SEIPRR, from parameter values
calc_R0_SEIPRR <- function(beta_input, parameters, pop_params,
                           beta_or_factor = "beta", virus, timestep){
  # if (beta_or_factor == "beta"){
  #   beta  = beta_input
  # } else if(beta_or_factor == "factor"){
  #   beta = parameters["beta_covid"]
  # }

  # specify parameters
  age_rates = pop_params[["ageing"]]
  pop_numbers = pop_params[["pop_numbers"]]
  
  #calc deaths for this timestep (last two age groups)
  deaths <- (pop_numbers[15:16]/sum(pop_numbers[15:16])*pop_params[["deaths"]])
  deaths <- deaths/ pop_numbers[15:16]
  # add dying to the older age groups
  age_rates[15] <- age_rates[15] + deaths[1]
  age_rates[16] <- age_rates[16] + deaths[2]
  
  beta <- beta_input
  
  # select the relevant beta
  if(virus == "covid"){
    gamma = parameters$gamma_covid
    incubation = parameters$incubation_covid
    sus_rates = parameters$red_susc_covid
    immune = parameters$covid_imm
    contactmatrix <- calc_contacts_timestep(timestep,parameters,"polymod")
    
    
  }else{
    gamma = parameters$gamma_other
    incubation = parameters$incubation_other
    sus_rates = rep(1, length(age_groups))
    immune = rep(0,length(age_groups))
    # seasonal forcing to the beta

    if (timestep ==1){
      beta = beta
    }else {
      beta = parameters$amplitude * cos(2*pi*(timestep-parameters$phi)/364) + beta
      }
    contactmatrix = parameters$contacts_all
      #parameters$contacts_hh +  parameters$contacts_school +  parameters$contacts_other
  }
  
  
  Transmission <- matrix(nrow=length(pop_numbers)*2, ncol=length(pop_numbers)*2, 0)
  for (i in 1:length(pop_numbers)){
    for (j in 1:length(pop_numbers)) {
      Transmission[((i*2)-1),((j*2))] <- beta*contactmatrix[i,j]*pop_numbers[j]*
        sus_rates[i]* (1-immune[i])
    }
  }
  
  #Transition matrix
  Transition <- matrix(0,nrow = length(pop_numbers)*2, ncol = length(pop_numbers)*2 )

  for(i in 1:length(pop_numbers)){
    for (j in 1:length(pop_numbers)){
      # create diagonal subsections
      if(i ==j){
        # top left
        Transition[((i*2)-1),((j*2)-1)] <- -(incubation + age_rates[i])
        # top right
        Transition[((i*2)-1),((j*2))] <- 0
        # bottom left 
        Transition[((i*2)),((j*2)-1)] <- incubation
        # bottom right
        Transition[((i*2)),((j*2))] <- -(gamma + age_rates[i])
      }
      if (i == (j+1) ){
        Transition[((i*2)-1),((j*2)-1)] <- age_rates[j]
        # top right
        Transition[((i*2)-1),((j*2))] <- 0
        # bottom left 
        Transition[((i*2)),((j*2)-1)] <- 0
        # bottom right
        Transition[((i*2)),((j*2))] <- age_rates[j]
      }
    }
  }
  #Inverse 
  
  Transition_inverse<- ginv(Transition)
  NGM <- -Transmission%*%Transition_inverse
  
  Eigen<- unlist((eigen(NGM)[1]))
  R0 <-  max(Eigen)

  return(R0)
}  

# summatisr into the PHE groups
summary_groups <- function(dt, ages_wanted= "PHE"){
  
  # reduce to just incidence measures
  cols_wanted <- colnames(dt)[grep("COVID_", x=colnames(dt))]
  cols_wanted <- c("time",cols_wanted,colnames(dt)[grep("OTHER_", x=colnames(dt))] )
  dt <- dt[,..cols_wanted]
  
  if(ages_wanted == "PHE"){

    
    dt[, OTHER_p0 := OTHER_0]
    dt[, OTHER_p5 := OTHER_5 + OTHER_10]
    dt[, OTHER_p15 := OTHER_15 + OTHER_20 + OTHER_25 + OTHER_30+ OTHER_35 + 
         OTHER_40]
    dt[, OTHER_p45 := OTHER_45 + OTHER_50 + OTHER_55 + OTHER_60]
    dt[, OTHER_p65 := OTHER_65 + OTHER_70 + OTHER_75]
    
    keepers <- c("time",colnames(dt)[grep("p", colnames(dt))])
    dt <- dt[,..keepers]
  }
  
  
  return(dt)
}


temp_func <- function(Temperatures, i){
  beta_i <- 1/Temperatures[i]
  return(beta_i)
}

step_fn <- function(parameters,  covmat,
                    lower_bounds, upper_bounds,
                    theta_names){
  #print(theta_names)  


  parameters_new <- c(rtmvnorm(1,
                              mean = unlist(parameters[theta_names], use.names=F),
                              sigma = covmat[theta_names,
                                             theta_names], 
                              lower = lower_bounds, 
                              upper = upper_bounds))
  
  names(parameters_new) <- theta_names
  
  for(param_name in names(parameters_new)){
    parameters[param_name] <- parameters_new[param_name]}
  
  return(parameters)
}

swapping_rates <- function(trace_in, n_pars, prev_swap_proposed, prev_swap_accepted){
  
  swap_tot <- sum(trace_in[,n_pars+4], na.rm = T) # total swaps proposed
  swap_tot <- swap_tot + prev_swap_proposed
  swap_accep <- length(which(trace_in[,n_pars+3]!= 0)) # total swaps accepted
  swap_accep <- swap_accep + prev_swap_accepted
  
  if(swap_tot ==0){
    swap_rate <- 0} else {
      swap_rate <- swap_accep/swap_tot
      
    }
  return(swap_rate)
}


adjust_unsymmetrical <- function(log_acceptance, parameters_old,
                                 parameters_new, covmat, lower_bounds, upper_bounds,
                                 theta_names){
  log_acceptance <- log_acceptance + dtmvnorm(x = unlist(parameters_old[theta_names]),
                                              mean = unlist(parameters_new[theta_names]),
                                              sigma = covmat[theta_names,
                                                             theta_names],
                                              lower = lower_bounds[theta_names],
                                              upper = upper_bounds[theta_names],
                                              log = TRUE)
  
  log_acceptance <- log_acceptance - dtmvnorm(x = unlist(parameters_new[theta_names]),
                                              mean = unlist(parameters_old[theta_names]),
                                              sigma = covmat[theta_names,
                                                             theta_names],
                                              lower = lower_bounds[theta_names],
                                              upper = upper_bounds[theta_names],
                                              log = TRUE)
  
  return(log_acceptance)
}

convert_sig <- function(value){
  # Convert to model sigma value
  # value <- exp(value)/ (1+exp(value))
  
  if (value == 0 ){
    value <- 1
  } else {
    value <- (1 - value)
  }
  
  return(value)
}

# calculate the percent infected per age group (reporeted)
infected_year_age <- function(outall, type, year_start){
  
  Names <-naming_states(type = type)
  #Calculate the incidence
  colnames(outall) <- Names
  # calculate the number infected over 364 days by age group
  year_change <- matrix(nrow=1, ncol = length(age_groups))
  colnames(year_change) = paste0("infected_year_",age_groups)
  
  
  year_change[,1] <- outall[year_start+364, Other_reported1]- outall[year_start, Other_reported1]
  year_change[,2] <- outall[year_start+364, Other_reported2]- outall[year_start, Other_reported2]
  year_change[,3] <- outall[year_start+364, Other_reported3]- outall[year_start, Other_reported3]
  year_change[,4] <- outall[year_start+364, Other_reported4]- outall[year_start, Other_reported4]
  year_change[,5] <- outall[year_start+364, Other_reported5]- outall[year_start, Other_reported5]
  year_change[,6] <- outall[year_start+364, Other_reported6]- outall[year_start, Other_reported6]
  year_change[,7] <- outall[year_start+364, Other_reported7]- outall[year_start, Other_reported7]
  year_change[,8] <- outall[year_start+364, Other_reported8]- outall[year_start, Other_reported8]
  year_change[,9] <- outall[year_start+364, Other_reported9]- outall[year_start, Other_reported9]
  year_change[,10] <- outall[year_start+364, Other_reported10]- outall[year_start, Other_reported10]
  year_change[,11] <- outall[year_start+364, Other_reported11]- outall[year_start, Other_reported11]
  year_change[,12] <- outall[year_start+364, Other_reported12]- outall[year_start, Other_reported12]
  year_change[,13] <- outall[year_start+364, Other_reported13]- outall[year_start, Other_reported13]
  year_change[,14] <- outall[year_start+364, Other_reported14]- outall[year_start, Other_reported14]
  year_change[,15] <- outall[year_start+364, Other_reported15]- outall[year_start, Other_reported15]
  year_change[,16] <- outall[year_start+364, Other_reported16]- outall[year_start, Other_reported16]
  # for(test in 1:length(age_groups)){
  #   year_change[,test] <- outall[year_start+364, get(paste0("Other_reported",test))] - 
  #     outall[year_start, get(paste0("Other_reported",test))]
  # }
  # 
  year_change <- data.table(year_change)
  year_change[, infected_year_p0 := infected_year_0]
  year_change[, infected_year_p5 := infected_year_5 + infected_year_10]
  year_change[, infected_year_p15 := infected_year_15 + infected_year_20 + 
                infected_year_25 + infected_year_30 + infected_year_35 + infected_year_40]
  year_change[, infected_year_p45 := infected_year_45 + infected_year_50 + 
                infected_year_55 + infected_year_60]
  year_change[, infected_year_p65 := infected_year_65 + infected_year_70 + 
                infected_year_75 ]
  
  colnames_wanted <- colnames(year_change)[grep("year_p", colnames(year_change))]
  
  return(year_change[,..colnames_wanted])
  
}

# number swaps proposed(takes into account those from previous wrap)
swapping_props <- function(trace_in, n_pars, prev_swap_proposed){
  
  swap_tot <- sum(trace_in[,n_pars+4], na.rm = T) # total swaps proposed
  swap_tot <- swap_tot + prev_swap_proposed
  
  return(swap_tot)
}

# number swaps accepted (takes into account those from previous wrap)
swapping_acceps <- function(trace_in, n_pars, prev_swap_accepted){
  
  swap_accep <- length(which(trace_in[,n_pars+3]!= 0)) # total swaps accepted
  swap_accep <- swap_accep + prev_swap_accepted
  
  return(swap_accep)
}

#add on to previous
add_on <- function(tt, t){
  tout <- rbind(tt, t)
  return(tout)
  
}