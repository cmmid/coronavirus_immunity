################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################



label_trace <- function(trace_in){
 
    colnames(trace_in) <- c("log_liklihood", 
                            "waning_duration", "seasonal_R0", 
                            "reporting_rate_1", "reporting_rate_2", "reporting_rate_3", 
                            "reporting_rate_4", "seasonal_amplitude", "seasonal_timing",
                            "Step_accepted","Swap_accepted_to",
                            "swap_attempted")
  
  return(trace_in)
}


format_trace <- function(trace_in, thin = 1, burnin=1, keep_all = F){
  
  trace_dt <- data.table(trace_in)
  trace_dt[,"time_steps"] = timesteps
  
  trace_dt[,reporting_rate_1 := exp(reporting_rate_1)/(1+exp(reporting_rate_1))]
  trace_dt[,reporting_rate_2 := exp(reporting_rate_2)/(1+exp(reporting_rate_2))]
  trace_dt[,reporting_rate_3 := exp(reporting_rate_3)/(1+exp(reporting_rate_3))]
  trace_dt[,reporting_rate_4 := exp(reporting_rate_4)/(1+exp(reporting_rate_4))]
  
  
  trace_dt <- trace_dt[burnin:dim(trace_dt)[1],]
  trace_dt <- trace_dt[seq(1,dim(trace_dt)[1], by=thin),]
  if(keep_all == F){
    trace_dt[,c("Step_accepted", "Swap_accepted_to", "swap_attempted") := NULL]}
  
  return(trace_dt)
}

plot_trace <- function(trace_in, multi=T){
  if (multi == T){
    TRACE_PLOTS <- ggplot(as.data.frame(trace_in)) +
      geom_line(aes(y=value, x=time_steps, colour=factor(temp))) +
      facet_wrap(variable~., scales="free", ncol=2 ) +
      theme(
        legend.title = element_text(size = 10),
        legend.text  = element_text(size = 10),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10, margin = margin(l=0,r=100,b=0,t=0)),
        title = element_text(size=10),
        legend.position = "none")
  } else {
    TRACE_PLOTS <- ggplot(as.data.frame(trace_in)) +
      geom_line(aes(y=value, x=time_steps)) +
      facet_wrap(variable~., scales="free", ncol=2) +
      theme(
        legend.title = element_text(size = 10),
        legend.text  = element_text(size = 10),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10, margin = margin(l=0,r=100,b=0,t=0)),
        title = element_text(size=10))
  }
  
  return(TRACE_PLOTS)
}

swaps_accepted <- function(trace_in, subset_group){
  
  num <- sum(trace_in[subset_group,"Swap_accepted_to"]!=0)
  denom <- sum(trace_in[subset_group,"swap_attempted"])

  rate <- num / denom
  return(rate)
}

plot_density <- function(trace_in, multi = T){
  if (multi == T){
    DENS_PLOTS <- ggplot(as.data.frame(trace_combined)) +
      geom_density(aes(x=value, colour=factor(temp))) +
      facet_wrap(variable~., scales="free")
  } else {
    DENS_PLOTS <- ggplot(as.data.frame(trace_combined)) +
      geom_density(aes(x=value)) +
      facet_wrap(variable~., scales="free")
  }
  return(DENS_PLOTS)
}

plot_rbinom <- function(samples, 
                        trace_period,
                        trace_dt, 
                        model_type){
  samples_to_take <- samples
  reporting_store <- data.frame()
  
  for(i in 1:samples_to_take){
    # or manually write them
    sample_num <- sample(trace_period, size=1)
    params_in <- unlist(trace_dt[sample_num,])
    init_theta <- c()
    init_theta["waning_day"] <- params_in["waning_duration"]
    init_theta["seasonal_R0"] <- params_in["seasonal_R0"]
    init_theta["seasonal_reported_1"] <- log(params_in["reporting_rate_1"]/(1-params_in["reporting_rate_1"]))
    init_theta["seasonal_reported_2"] <-  log(params_in["reporting_rate_2"]/(1-params_in["reporting_rate_2"]))
    init_theta["seasonal_reported_3"] <-  log(params_in["reporting_rate_3"]/(1-params_in["reporting_rate_3"]))
    init_theta["seasonal_reported_5"] <-  log(params_in["reporting_rate_4"]/(1-params_in["reporting_rate_4"]))
    init_theta["seasonal_amplitude"]<- params_in["seasonal_amplitude"]
    init_theta["phi"] <- params_in["seasonal_timing"]

    #run the model up unitl likelihood point
    parameters <- create_parameters(unlist(init_theta))
    # Run the model - seasonal corona only. For 20 years + a bit, to find low point (lp)
    output_s <- run_model_seasonal(parameters, model_type)
    colnames(output_s) <- naming_states(model_type)
    # summaries
    reportin_2020 <- summary_stats_reported_seasonal(output_s, type = model_type)
    reportin_2020_daily <- summary_groups(reportin_2020)
    reportin_seasonal <- melt.data.table(reportin_2020_daily[,c("time", ..all_oths)],
                                         id= "time")
    reportin_seasonal[, date := as.Date(time, origin = lp_15)]
    #change from daily to monthly
    for(weeker in 1:(length(seasonal_dates_15)-1)){
   
      reportin_seasonal[ date >= seasonal_dates_15[weeker] &
                           date < seasonal_dates_15[weeker+1],
                         year_week := seasonal_dates_15[weeker]]
    }
    reportin_seasonal <- na.omit(reportin_seasonal[,sum(value, na.rm=T), by=c("year_week",
                                                                              "variable")
    ])

    for( stepper in 1:dim(reportin_seasonal)[1]){
      if(reportin_seasonal[stepper,"variable"] == "OTHER_p0"){age_set = 1
      } else if(reportin_seasonal[stepper,"variable"] == "OTHER_p5"){age_set = 2
      } else if(reportin_seasonal[stepper,"variable"] == "OTHER_p15"){age_set = 3
      } else if(reportin_seasonal[stepper,"variable"] == "OTHER_p45"){age_set = 4
      } else if(reportin_seasonal[stepper,"variable"] == "OTHER_p65"){age_set = 5}
      # work out quantile intervals
      reportin_seasonal[stepper,"rbb"] <- rbinom(n=1,
                                                 size = as.numeric(round(reportin_seasonal[stepper,"V1"])),
                                                 prob = as.numeric(parameters$seasonal_reported[age_set])
      )

    }
    
    reportin_seasonal[,sample_num := sample_num]
    
    reporting_store <- rbind(reporting_store, reportin_seasonal)
    print(i)
  }
  reporting_store[seasonal_15_20, on = c("year_week", "variable"), true_value:= i.value ]

  reporting_store <- rbind(reporting_store,data.frame("2017-04-03", "OTHER_p5", NA, NA, NA ,NA), use.names=F)
  reporting_store <- rbind(reporting_store,data.frame("2017-04-03", "OTHER_p15", NA, NA, NA ,NA), use.names=F)
  reporting_store <- rbind(reporting_store,data.frame("2017-04-03", "OTHER_p45", NA, NA, NA ,NA), use.names=F)
  reporting_store <- rbind(reporting_store,data.frame("2017-04-03", "OTHER_p65", NA, NA, NA ,NA), use.names=F)
  reporting_store <- rbind(reporting_store,data.frame("2017-04-03", "OTHER_p0", NA, NA, NA ,NA), use.names=F)
  reporting_store$variance_indicator<- as.factor(reporting_store$variance_indicator)
  
  reporting_store[variable == "OTHER_p0", nice_label := "Age 0 - 4"]
  reporting_store[variable == "OTHER_p5", nice_label := "Age 5 - 14"]
  reporting_store[variable == "OTHER_p15", nice_label := "Age 15 - 44"]
  reporting_store[variable == "OTHER_p45", nice_label := "Age 45 - 64"]
  reporting_store[variable == "OTHER_p65", nice_label := "Age 65 +"]

  reporting_store$year_week <- as.Date(reporting_store$year_week)
  RBINOM <- ggplot(reporting_store, aes(x = year_week, y = rbb),) + 
    geom_point(alpha = 0.4, colour = "deepskyblue2") +
    facet_wrap(nice_label~., scales = "free_y", ncol=1) + 
    theme_linedraw()+
    scale_x_date(date_breaks = "1 year")+
    labs(x = "Month", y = "Number of infections reported") +
    geom_point(aes(y = true_value), colour = "black",) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    theme(legend.position = "none", strip.text.y = element_text(angle = 0), 
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(colour = 'black',hjust = 0))
  return(RBINOM)}



attack_rate<- function(samples, 
                       trace_period,
                       trace_dt, 
                       model_type){
  samples_to_take <- samples
  reporting_store <- data.frame()
  
  for(i in 1:samples_to_take){
    for(j in c(1:5)){
      
      # or manually write them
      sample_num <- sample(trace_period, size=1)
      params_in <- unlist(trace_dt[sample_num,])
      init_theta <- params_in[2:9]
      #run the model up unitl likelihood point
      parameters <- create_parameters(unlist(init_theta))
      # Run the model - seasonal corona only. For 20 years + a bit, to find low point (lp)
      output_s <- run_model_seasonal(parameters, model_type)
      colnames(output_s) <- naming_states(model_type)
      
      population <- calc_pop_per_age(output_s,timepoint=(length_to_run -(j*364)), model_type)
      population$inf_year_age <- c(infected_year_age(output_s, model_type,
                                                     year_start  = length_to_run-364))
      population[,percent := as.numeric(inf_year_age)/V1*100]
      population$percent
      reporting_store <- rbind(  reporting_store, c(population$percent,sample_num,j))
    }}
  reporting_store <- data.table(reporting_store)
  colnames(reporting_store)<- c("Age1", "Age2", "Age3", "Age4", "Age5", "sample", "year")
  reporting_store_m <- melt.data.table(reporting_store, id.vars= c("sample"))
  reporting_store_m[,mean(value), by = variable]
  
  return(reporting_store_m[,mean(value), by = variable])
  
}


