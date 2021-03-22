################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################

# projections simulations functions

#FUNCTION: run the model for the second set up to including 2019
run_model_project <- function(parameters, init_state_2020){
  #times over which to run the model
  times_20 <- c(1:length_to_run_3)
  # run the model
  outall <- as.data.table(ode(
    y = unlist(init_state_2020),
    t = times_20,
    func = SEIR_2virus_cons_season,
    parms = parameters,
    method = "rk4"))
  
  return(outall)
  
}

run_model_to_project <- function(test_params, parameters_2020,
                                sigma, init_state_2020){

  # combine the parameters

  parameters <- update_parameters_specific(test_params = test_params,
                                           parameters = parameters_2020,
                                           sigma = sigma, 
                                           other_way = T)
  
if(seasonal_factor_covid == "yes"){
  # run the model
  parameters$amplitude_covid <- parameters$beta_covid_0 * (parameters$amplitude/
                                                             parameters$beta_other)
} else {parameters$amplitude_covid <- 0}
  output_project <- run_model_project(parameters = parameters,
                                init_state_2020 = init_state_2020)
 
  # create summaries
  reportin_2020 <- summary_stats_reported_proj(output_project, type = "SEIR")
  #out <- summary_groups_both(reportin_2020)
  # deaths <- calculate_deaths(output_project, parameters)
  return(reportin_2020)
  
}
  
summary_stats_reported_proj <- function(outall, type) {
  
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
  
  outall[,other_all := OTHER_0 +  OTHER_5 +  OTHER_10 +  OTHER_15 +
           OTHER_20 +  OTHER_25 +  OTHER_30 +  OTHER_35 +
           OTHER_40 +  OTHER_45 +  OTHER_50 +  OTHER_55 +
           OTHER_60 +  OTHER_65 +  OTHER_70 +  OTHER_75 ]
  
  outall[,covid_all := COVID_0 +  COVID_5 +  COVID_10 +  COVID_15 +
           COVID_20 +  COVID_25 +  COVID_30 +  COVID_35 +
           COVID_40 +  COVID_45 +  COVID_50 +  COVID_55 +
           COVID_60 +  COVID_65 +  COVID_70 +  COVID_75 ]
outall <- outall[,c("time", "other_all", "covid_all")]
  
  return(outall)
}

  
  
summary_groups_proj <- function(dt){
  
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
  
  dt[, OTHER_p0 := OTHER_0]
  dt[, OTHER_p5 := OTHER_5 + OTHER_10]
  dt[, OTHER_p15 := OTHER_15 + OTHER_20 + OTHER_25 + OTHER_30+ OTHER_35 +
       OTHER_40]
  dt[, OTHER_p45 := OTHER_45 + OTHER_50 + OTHER_55 + OTHER_60]
  dt[, OTHER_p65 := OTHER_65 + OTHER_70 + OTHER_75]
  
  keepers <- c("time",colnames(dt)[grep("p", colnames(dt))])
  dt <- dt[,..keepers]
  
  return(dt)
}
  