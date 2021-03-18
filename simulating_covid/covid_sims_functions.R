################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################

# Functions for covid simulations
library(RcppSims)
library(ggplot2)
#function that takes a sample as input and outputs the modelled daily deaths +
# the modelled serology data and the Reffective over time
run_model_on_sample <- function(test_params, parameters_2020,
                                sigma, init_state_2020){
  # combine the parameters
  parameters <- update_parameters_specific(test_params = test_params,
                                           parameters = parameters_2020,
                                           sigma = sigma)
  # run the model 

  output_2020 <- run_model_2020(parameters = parameters,
                                init_state_2020 = init_state_2020)
  # calculate the deaths from model otuput
  model_deaths <- calculate_deaths(model_output = output_2020,
                                   parameters = parameters)

  # calculate the serology data
  serology <- calc_youngest_ages(output_2020, parameters)
 
  # calculate Reffective over time
  r_effective <- calc_Reff_time(parameters = parameters,
                                outall = output_2020, 
                                timespan = c(1,dim(output_2020)[1]))
  r_effective <- data.table(r_eff = r_effective, 
                            timestep = c(1:dim(output_2020)[1]))

  #calculate the R0
  r_0 <- calc_R0_SEIR(parameters$beta_covid_0, parameters, timestep=1)
  # return the modelled outputs as a list
  return(list(model_deaths= model_deaths[,c("date", "total")], 
              r_eff = r_effective, 
              r_0 = r_0, 
              serology = serology))
}

#FUNCTION: Take a set of seasonal parameters and generate 2020 initial state
run_seasonal_for_init <- function(parameter_guesses, model_type){
  # put into parameter list
  parameters <- create_parameters(parameter_guesses = parameter_guesses)
  #run the model
  output_s <- run_model_seasonal(parameters, model_type = model_type)
  colnames(output_s) <- naming_states(model_type)
  # change step to date and save the states @ required start date (run_start_2)
  output_s[, date := as.Date(time, origin = lp_15)]
  init_state_2020 <- output_s[date == run_start_2,2:(ncol(output_s)-1)]
  #return the new initial state
  return(list(init_state_2020 =init_state_2020, parameters_2019 = parameters))
}

optimise_ll_deaths <- function(test_params, parameters_2020,
                               sigma, init_state_2020){
  # combine the parametes
  parameters <- update_parameters_specific(test_params = test_params,
                                           parameters = parameters_2020,
                                           sigma = sigma)
  # run the model 
  output_2020 <- run_model_2020(parameters = parameters,
                                init_state_2020 = init_state_2020)
  # calculate the deaths from model otuput
  model_deaths <- calculate_deaths(model_output = output_2020,
                                   parameters = parameters)
  # calculate the log likelihood
  ll_covid_deaths <- calc_ll_covid_deaths(covid_deaths_daily = model_deaths)
  # return the log likelihood
  
  return(ll_covid_deaths)
}

#FUNCTION: Calculate the log likelihood of deaths, using a poisson distribution
calc_ll_covid_deaths <- function(covid_deaths_daily){
  # total number of deaths
  deaths_model <- unique(covid_deaths_daily[,c("date", "total")])
  colnames(deaths_model) <- c("date", "model_deaths")
  # combine wiht real data
  deaths_model[deaths_covid, on = "date", actual_deaths := i.actual_deaths]
  # calculate with poisson
  deaths_model[, likelihood := dpois(x = as.numeric(actual_deaths), 
                                     lambda = as.numeric(model_deaths), log=T) ]
  #return the total likelihood
  if(any(is.nan(deaths_model$likelihood))){browser()}
  return(-sum(deaths_model$likelihood, na.rm = T))
}

#FUNCTION: take the model output and use it to calculate daily modelled deaths
calculate_deaths <- function(model_output, parameters){
  # create summaries
  reportin_2020 <- summary_stats_reported_both(model_output, type = "SEIR")
  daily_infections <- summary_groups_both(reportin_2020)
  # extract relevant dates
  daily_infections[, date := as.Date(time, origin = run_start_2)]
  daily_infections <- daily_infections[date > run_start_2 & 
                                        date < run_end,]
# calculate the deaths using the death rate by age group
  daily_infections[, COVID_p0 := COVID_p0 * parameters$covid_deaths[1]]
  daily_infections[, COVID_p5 := COVID_p5 * parameters$covid_deaths[2]]
  daily_infections[, COVID_p15 := COVID_p15 * parameters$covid_deaths[3]]
  daily_infections[, COVID_p45 := COVID_p45 * parameters$covid_deaths[4]]
  daily_infections[, COVID_p65 := COVID_p65 * parameters$covid_deaths[5]]
# combine to get total daily infections  
  daily_infections[, total := sum(COVID_p0 +COVID_p5 + COVID_p15 + 
                                    COVID_p45 + COVID_p65), by=date][
                                      , c("COVID_p0","COVID_p5","COVID_p15",
                                          "COVID_p45","COVID_p65") := NULL
                                    ]
#return the total daily deaths
  return(daily_infections)
}

#FUNCTION: run the model for the second set up to including 2019
run_model_2020 <- function(parameters, init_state_2020){
  #times over which to run the model
  times_20 <- c(1:length_to_run_2)
  # run the model
  outall <- as.data.table(ode(
    y = unlist(init_state_2020),
    t = times_20,
    func = SEIR_2virus_cons, #SEIR_2virus_cons_ld
    parms = parameters,
    method = "rk4"))
  
  return(outall)
  
}

update_parameters <- function(parameters){
  parameters$covid_intros <- c(rep(0,4),
                               rep(4,9),
                               rep(0,3))

  parameters$sig_solid_R <- 0
  parameters$social_distancing <- 0.33
  parameters$red_susc_covid <- rep(1,16)

  parameters$covid_time <-  (as.Date(covid_date) - as.Date(run_start_2))[[1]]
  parameters$covid_change_time <- (as.Date(covid_change_time)- as.Date(covid_date))[[1]]
  parameters$mobility_start <-  (as.Date(mobility_date) - run_start_2)[[1]]
  
  parameters$contacts_reduc <-(((google_mobility_use$rolling_mean/100)) +1)
  parameters$waning_covid <- parameters$waning_other
  parameters$gamma_covid <- 1/4
  parameters$incubation_covid <- 1/3
  parameters$inf_to_symp_covid = 1/(22/2)
  parameters$reporting_delay_covid = 1/(22/2)
  parameters$covid_deaths[1] =  as.numeric(IFRs_actual[1,"proportion"]) #0.00004,
  parameters$covid_deaths[2] =  as.numeric(IFRs_actual[2,"proportion"]) #0.00004,
  parameters$covid_deaths[3] =  as.numeric(IFRs_actual[3,"proportion"])  #0.000253,
  parameters$covid_deaths[4] =  as.numeric(IFRs_actual[4,"proportion"]) #0.0049,
  parameters$covid_deaths[5] =  as.numeric(IFRs_actual[5,"proportion"]) #0.131,
  

  return(parameters)
}

#Update the transmission and interaction parameters
update_parameters_specific <- function(test_params, parameters, sigma){
 
  parameters$sig_R <- convert_sig(sigma)
  parameters$sig_I <- convert_sig(sigma*0.5)
  parameters$beta_covid_0 <- test_params[1]
  parameters$beta_covid_1 <- test_params[1]
  parameters$covid_intros <- c(rep(0,4),
                               rep(test_params[2],9),
                               rep(0,3))

  return(parameters)
}


# summarise both - specifically coded for speed
summary_stats_reported_both <- function(outall, type) {

  #Calculate the incidence
  colnames(outall) <- naming_states(type = type)

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
  
  outall[, COVID_0 := Covid_reported1 - shift(Covid_reported1, 1L, type = "lag")]
  outall[, COVID_5 := Covid_reported2 - shift(Covid_reported2, 1L, type = "lag")]
  outall[, COVID_10 := Covid_reported3 - shift(Covid_reported3, 1L, type = "lag")]
  outall[, COVID_15 := Covid_reported4 - shift(Covid_reported4, 1L, type = "lag")]
  outall[, COVID_20 := Covid_reported5 - shift(Covid_reported5, 1L, type = "lag")]
  outall[, COVID_25 := Covid_reported6 - shift(Covid_reported6, 1L, type = "lag")]
  outall[, COVID_30 := Covid_reported7 - shift(Covid_reported7, 1L, type = "lag")]
  outall[, COVID_35 := Covid_reported8 - shift(Covid_reported8, 1L, type = "lag")]
  outall[, COVID_40 := Covid_reported9 - shift(Covid_reported9, 1L, type = "lag")]
  outall[, COVID_45 := Covid_reported10 - shift(Covid_reported10, 1L, type = "lag")]
  outall[, COVID_50 := Covid_reported11 - shift(Covid_reported11, 1L, type = "lag")]
  outall[, COVID_55 := Covid_reported12 - shift(Covid_reported12, 1L, type = "lag")]
  outall[, COVID_60 := Covid_reported13 - shift(Covid_reported13, 1L, type = "lag")]
  outall[, COVID_65 := Covid_reported14 - shift(Covid_reported14, 1L, type = "lag")]
  outall[, COVID_70 := Covid_reported15 - shift(Covid_reported15, 1L, type = "lag")]
  outall[, COVID_75 := Covid_reported16 - shift(Covid_reported16, 1L, type = "lag")]
  
 outall <- outall[-1,]
  
  
  return(outall)
}

summary_groups_both <- function(dt){
  
  # reduce to just incidence measures
  cols_wanted <- colnames(dt)[grep("COVID_", x=colnames(dt))]
  cols_wanted <- c("time",cols_wanted,colnames(dt)[grep("OTHER_", x=colnames(dt))] )
  dt <- dt[,..cols_wanted]
    
    dt[, COVID_p0 := COVID_0]
    dt[, COVID_p5 := COVID_5 + COVID_10]
    dt[, COVID_p15 := COVID_15 + COVID_20 + COVID_25 + COVID_30+ COVID_35 + 
         COVID_40]
    dt[, COVID_p45 := COVID_45 + COVID_50 + COVID_55 + COVID_60]
    dt[, COVID_p65 := COVID_65 + COVID_70 + COVID_75]
    
    # dt[, OTHER_p0 := OTHER_0]
    # dt[, OTHER_p5 := OTHER_5 + OTHER_10]
    # dt[, OTHER_p15 := OTHER_15 + OTHER_20 + OTHER_25 + OTHER_30+ OTHER_35 + 
    #      OTHER_40]
    # dt[, OTHER_p45 := OTHER_45 + OTHER_50 + OTHER_55 + OTHER_60]
    # dt[, OTHER_p65 := OTHER_65 + OTHER_70 + OTHER_75]
    
    keepers <- c("time",colnames(dt)[grep("p", colnames(dt))])
    dt <- dt[,..keepers]
  
  return(dt)
}


calc_R0_SEIR <- function(beta_input, parameters, timestep){
  # specify parameters
  age_rates = parameters$age_rates
  pop_numbers = pop_numbers
  
  #calc deaths for this timestep (last two age groups)
  deaths <- (pop_numbers[15:16]/sum(pop_numbers[15:16])*parameters$deaths)
  deaths <- deaths/ pop_numbers[15:16]
  # add dying to the older age groups
  age_rates[15] <- age_rates[15] + deaths[1]
  age_rates[16] <- age_rates[16] + deaths[2]
  
  beta <- beta_input
  gamma = parameters$gamma_covid
  incubation = parameters$incubation_covid
  sus_rates = parameters$red_susc_covid
  immune = parameters$covid_imm
  contactmatrix <- calc_contacts_timestep(timestep=1,parameters,"polymod")
    

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


#function to get desired beta from R0 using optimisation algorithm.
calc_beta_SEIR <- function(start_val, requiredR0, parameters,
                            timestep=1){
  
  best_fit <- optim(par = start_val, 
                    fn = optim_beta_SEIR,
                    Ro_wanted = requiredR0, 
                    parameters = parameters,
                    timestep = timestep,
                    method="Brent", 
                    lower = 0, 
                    upper = 0.5)
  
  return(best_fit)
}

# function used to optimise the beta parameter based on a desired Ro. SEIPRR
optim_beta_SEIR <- function(beta_input, Ro_wanted, parameters,
                             timestep){
  
  R0 <- calc_R0_SEIR(beta_input = beta_input, parameters = parameters,
                       timestep = timestep)
  difference = (abs(R0 - Ro_wanted)^2)
  return(difference)
}  

#calculate the Reffective over a specified time period
calc_Reff_time <- function(parameters, outall, timespan, thinning = 1, excl_sus = "no"){
  
  colnames(outall) <- naming_states("SEIR")
  # Calculate proportion susceptible
  outall_temp <- calc_sus_for_NGM(outall, parameters)
  outall_temp <- outall_temp[,(dim(outall_temp)[2]-length(age_groups)+1): dim(outall_temp)[2]]
  
  R0_store <- c()
  # for each timestep
  for(timestep in seq(from = timespan[1], to = timespan[2], by=thinning)){
    #calculate the R0
    NGM_to_use <- calc_NGM_change(outall = outall, 
                                 parameters = parameters,
                                  timestep=timestep)
    
    # Times by the proportion susceptible 
    if(excl_sus == "no"){
      temp <- unlist(outall_temp[timestep,])
      vector_sus <- rep(temp,each=2)
      
      NGM_to_use <- sweep(NGM_to_use, MARGIN = 1, vector_sus,"*")
    }
    
    # eigen values
    Eigen<- unlist((eigen(NGM_to_use)[1]))
    # max eigenvalue
    R0 <-  max(Eigen)
    #store!
    R0_store <- c(R0_store, R0)
    
  }
  return(R0_store)}


calc_sus_for_NGM <- function(outall, parameters){
  
      #Susceptible for each age group
        for(y in 1:length(age_groups)){
          
          agegroup <- age_groups[y]
          
          outall[,paste0("Susceptibility_", agegroup) :=
                   ((get(paste0("SS",y)) + get(paste0("SE",y))) + 
                      parameters$sig_R*(get(paste0("SI",y))  + get(paste0("SR",y))))/
                   (get(paste0("SS",y)) + get(paste0("IS",y))+get(paste0("RS",y))+
                      get(paste0("SI",y)) + get(paste0("II",y)) +get(paste0("RI",y))+
                      get(paste0("SR",y))+get(paste0("IR",y))+get(paste0("RR",y))+
                      get(paste0("ES",y))+get(paste0("SE",y))+get(paste0("EE",y))+get(paste0("IE",y))+
                      get(paste0("RE",y))+get(paste0("EI",y))+
                      get(paste0("ER",y)))]
  
  }
  return(outall)
}





# calculate NGM for either covid or seasonal
calc_NGM_change <- function(outall,parameters,timestep){
  
  # specify parameters
  age_rates = parameters$age_rates
  pop_numbers = calc_pop_per_age_all(outall,timestep, model_type="SEIR")
  
  #calc deaths for this timestep (last two age groups)
  deaths <- (pop_numbers[15:16]/sum(pop_numbers[15:16])*parameters$deaths)
  deaths <- deaths/ pop_numbers[15:16]
  # add dying to the older age groups
  age_rates[15] <- age_rates[15] + deaths[1]
  age_rates[16] <- age_rates[16] + deaths[2]
  
  # select the relevant beta
  beta <- parameters$beta_covid_0
  gamma = parameters$gamma_covid
  incubation = parameters$incubation_covid
  sus_rates = parameters$red_susc_covid
  immune = parameters$covid_imm
  contactmatrix <- calc_contacts_timestep(timestep,parameters,"polymod")
  
  #Work out the transmssion matrix (as per page 6 of notebook)
  Transmission <- matrix(nrow=length(pop_numbers)*2, ncol=length(pop_numbers)*2, 0)
  for (i in 1:length(pop_numbers)){
    for (j in 1:length(pop_numbers)) {
      
      Transmission[((i*2)-1),((j*2))] <- beta*contactmatrix[i,j]*pop_numbers[j]*
        sus_rates[i]* (1-immune[i])
    }
  }
  
  #Transition matrix (as per page 7 of notebook)
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
  return(NGM)
}


calc_contacts_timestep <- function(timestep, parameters, type){
  
  if(type == "polymod"){
    if(timestep<parameters$mobility_start){
      contactmatrix = parameters$contacts_hh +  parameters$contacts_school +  parameters$contacts_other;
    } else if(timestep>=(parameters$mobility_start)){
      contactmatrix = parameters$contacts_hh +  parameters$contacts_school + 
        (parameters$contacts_other*parameters$contacts_reduc[timestep - as.numeric(parameters$mobility_start)+1]);
      if(timestep>= parameters$covid_time + parameters$covid_change_time ){
        contactmatrix = contactmatrix - parameters$contacts_school
        contactmatrix = contactmatrix * parameters$social_distancing
        contactmatrix[1:3, 1:16] <- contactmatrix[1:3, 1:16]*parameters$child_extra_reduc
        contactmatrix[4:16, 1:3] <- contactmatrix[4:16, 1:3]*parameters$child_extra_reduc
      }}
  } else  if(type == "bbc" ){
    
    if(timestep< (parameters$covid_time + parameters$covid_change_time)){
      contactmatrix <- parameters$contacts_other
    } else if(t>=(parameters$covid_time + parameters$covid_change_time)){
      contactmatrix = parameters$contacts_hh * parameters$social_distancing
    }
    
  }
  
  return(contactmatrix)
}


# calc population per age group and return all age groups
calc_pop_per_age_all <- function(outall, timepoint, model_type){
  
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
  # calculate the population over time (have to remove tracking compartments)
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
  
  
  return(pop2$V1)
}


calc_youngest_ages <- function(input, parameters,type="SEIR"){

  # need the proportion in the positive compartment
  incidences <- summary_stats_reported_both(input, type = type )
  timepoint_needed_end <- (serology_time - run_start_2)[[1]]
  subset_input <- incidences
      for(y in 1:length(age_groups)){ 
        agegroup <- age_groups[y]
        
        subset_input[,paste0("pop_size", agegroup) :=
                       (get(paste0("SS",y)) + get(paste0("IS",y))+get(paste0("RS",y))+
                          get(paste0("SI",y)) + get(paste0("II",y)) +get(paste0("RI",y))+
                          get(paste0("SR",y))+get(paste0("IR",y))+get(paste0("RR",y))+
                          get(paste0("ES",y))+get(paste0("SE",y))+get(paste0("EE",y))+get(paste0("IE",y))+
                          get(paste0("RE",y))+get(paste0("EI",y))+
                          get(paste0("ER",y)))]
      }

  proportions <- c(sum(subset_input[,"COVID_0"])/as.numeric(subset_input[timepoint_needed_end, "pop_size0"]),
                   sum(subset_input[,"COVID_5"])/as.numeric(subset_input[timepoint_needed_end, "pop_size5"]),
                   sum(subset_input[,"COVID_10"])/as.numeric(subset_input[timepoint_needed_end, "pop_size10"]),
                   sum(subset_input[,"COVID_15"])/as.numeric(subset_input[timepoint_needed_end, "pop_size15"]),
                   sum(subset_input[,"COVID_20"])/as.numeric(subset_input[timepoint_needed_end, "pop_size20"]),
                   sum(subset_input[,"COVID_25"])/as.numeric(subset_input[timepoint_needed_end, "pop_size25"]),
                   sum(subset_input[,"COVID_30"])/as.numeric(subset_input[timepoint_needed_end, "pop_size30"]),
                   sum(subset_input[,"COVID_35"])/as.numeric(subset_input[timepoint_needed_end, "pop_size35"]),
                   sum(subset_input[,"COVID_40"])/as.numeric(subset_input[timepoint_needed_end, "pop_size40"]),
                   sum(subset_input[,"COVID_45"])/as.numeric(subset_input[timepoint_needed_end, "pop_size45"]),
                   sum(subset_input[,"COVID_50"])/as.numeric(subset_input[timepoint_needed_end, "pop_size50"]),
                   sum(subset_input[,"COVID_55"])/as.numeric(subset_input[timepoint_needed_end, "pop_size55"]),
                   sum(subset_input[,"COVID_60"])/as.numeric(subset_input[timepoint_needed_end, "pop_size60"]),
                   sum(subset_input[,"COVID_65"])/as.numeric(subset_input[timepoint_needed_end, "pop_size65"]),
                   sum(subset_input[,"COVID_70"])/as.numeric(subset_input[timepoint_needed_end, "pop_size70"]),
                   sum(subset_input[,"COVID_75"])/as.numeric(subset_input[timepoint_needed_end, "pop_size75"])
  )
  

  to_plot <- data.table(
    model = proportions,
    data = c(0.007,0.038,0.027,0.03,0.077, rep(NA,11)),
    ages = factor(c("0-4", "5-9","10-14", "15-19", "20-24", "25-29", "30-34", "35-39", 
                    "40-44", "45-49", "50-54", "55-59", "60-64" , "65-69", "70-74", "75+"), 
                  levels = c("0-4", "5-9","10-14", "15-19", "20-24", "25-29", "30-34", "35-39", 
                             "40-44", "45-49", "50-54", "55-59", "60-64",  "65-69", "70-74", "75+")),
    lower = c(0.000,0.002,0.00,0.001,0.018,rep(NA,11)), 
    upper = c(0.058,0.101,0.077,0.084,0.159,rep(NA,11))
  )
  return(to_plot)
}





  