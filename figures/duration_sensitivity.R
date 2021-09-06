# Duration of Immunity sensitivity 


immunity_durations <- seq(365, 3300, by = 365)

sero_sens <- data.frame()
for(j in immunity_durations){
for(i in 1:nrow(transmission)){
    # need to run a loop or an apply over the number of samples required.
    # take sample of seasonal parameters
    sample_num <- as.numeric(transmission[i, "sample"])
    seasonal_sample <- trace_to_sample[sample_num,]
    
    # run the seaonal model in order to get init start for the covid model
    seasonal_out <- run_seasonal_for_init(seasonal_sample, model_type = "SEIR")
    parameters_2020 <- update_parameters(seasonal_out[["parameters_2019"]])
    
    # loop here over different interaction parameters
    interaction_param <- as.numeric(transmission[i,"interaction"])
  
    parameters_2020$waning_covid <- parameters_2020$waning_other <- 1/j
    
    # combine the parameters
    parameters <- update_parameters_specific(test_params =  c(as.numeric(transmission[i,"transmission_val"]), 
                                                              as.numeric(introductions[i, "intro_val"])),
                                             parameters = parameters_2020,
                                             sigma = interaction_param)
    # run the model 
    output_2020 <- run_model_2020(parameters = parameters,
                                  init_state_2020 = unlist(seasonal_out["init_state_2020"]))
    # calculate the deaths from model otuput
    model_deaths <- calculate_deaths(model_output = output_2020,
                                     parameters = parameters)
    
    # calculate the serology data
    sero_sens_temp<-calc_youngest_ages(output_2020, parameters)
    sero_sens_temp$sample <- sample_num
    sero_sens_temp$interaction <- interaction_param
    sero_sens_temp$duration <- j

    sero_sens <- rbind(sero_sens,sero_sens_temp)
    print(i)
}}


#create the serology plot  
sero_sens$interaction <- as.factor(sero_sens$interaction)
sero_sens$ages <- forcats::fct_rev(sero_sens$ages)
# sero[, source := 2]
# sero[ages == "0-4" | ages == "5-9" | ages == "10-14" |
#        ages == "15-19" | ages == "20-24", source := 1]
# sero$source <- as.factor(sero$source)
# dup_point <- data.frame(x = "20-24", y = 0.102)


SERO_PLOT <- ggplot(sero_sens, aes(x= ages)) + 
  # geom_errorbar(aes(y=data, ymin=lower, ymax=upper)) +
  scale_colour_manual(values=cc) +
  labs(y = "Proportion positive", title = "Impact of different durations of immunity (days) on age-susceptibility", x = "Age group", colour = "Cross-protection",
       shape = "Source")+ 
  theme_linedraw() + 
  #  lims(y = c(0,0.2)) +
  geom_point(aes(y=model, colour = interaction, group = sample), alpha=0.8) + 
  theme(legend.position = "bottom") + 
  guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 1.5)))+
  facet_wrap(duration~. ) +
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8))+ 
  theme( legend.position = "bottom",
         legend.key.height = unit(0.3, "cm"), 
         axis.text.x = element_text(angle =95)) +
  geom_line(aes(x = ages, y = model, group = interaction(interaction,sample), colour = interaction))  + 
  geom_pointrange(data = sero_data, aes(x = age_group, y = data, ymin= lower, ymax = upper, shape = type), 
                  position = position_dodge(width=0.3))

tiff(here("figures","duration_sensitivity.tiff"), height = 2000, width = 4000, res = 300)
SERO_PLOT
dev.off()
