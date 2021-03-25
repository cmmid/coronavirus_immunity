################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################

# run the projections - will need to loop over the different cross-protections 
# and loop over different posterior samples

# overhead parameters that need changing
length_to_run_3 <- 4000
covid_change_time <- "3000-01-01"
mobility_date <- "3000-01-01"
year_times <- seq(from = 41, by = 365, length.out = 10)

set.seed(985)
storage <- data.frame()
storage2 <- data.frame()
######## first set #######

seasonal_factor_covid <- "yes" # "yes" or "no
# 1. get the sample 
for(i in 1:1#n_samples
){
projections_both <- data.frame()

for(sigma_otherway in c(0,"same")){
  projections_all <- data.frame()

    sample_num <- samples_to_take[i]
    seasonal_sample <- unlist(trace_to_sample[sample_num,])
    
    # run the seasonal model in order to get init start for the covid model
    seasonal_out <- run_seasonal_for_init(seasonal_sample, model_type = "SEIR")
    parameters_2020 <- update_parameters(seasonal_out[["parameters_2019"]])
    
    # loop here over different interaction parameters
    for(sig in sig_list){
      interaction_param <- sig
      # 2. get the fitted value
      transmission_val <- transmission[sample == sample_num & 
                                         interaction == interaction_param,
                                       transmission_val]
      introductions_val <- introductions[sample == sample_num & 
                                           interaction == interaction_param,
                                         intro_val]
      
      projections_out <- run_model_to_project(test_params  =  c(transmission_val, introductions_val) ,
                                          parameters_2020 = parameters_2020,
                                          sigma = interaction_param, 
                                          init_state_2020 = seasonal_out["init_state_2020"])
      
      for(year_time in year_times){
      age_inf <- infected_year_age(projections_out$all_output, type = "SEIR", 
                              year_start = year_time, time_frame = 364)
      tot <-calc_pop_per_age(projections_out$all_output, timepoint = year_time, "SEIR")
      storage_temp1 <- data.frame(age_inf/tot$V1)
      storage_temp1$time_start = year_time
      storage_temp1$sigma = sig
      storage_temp1$type = sigma_otherway
      storage <- rbind(storage, storage_temp1)
      if (year_time == 406){print(age_inf)}
    
      age_inf <- infected_year_age(projections_out$all_output, type = "SEIR",
                                   year_start = year_time, time_frame = 364)
      tot <- sum(age_inf)
      storage_temp2 <- data.frame(age_inf/tot)
      storage_temp2$time_start = year_time
      storage_temp2$sigma = sig
      storage_temp2$type = sigma_otherway
      storage2 <- rbind(storage2, storage_temp2)
      }
      
     projections <- projections_out$reported
      projections$sample <- sample_num
      projections$inter <- sig
      
      projections_all <- rbind(projections_all, projections)
    }

    if(sigma_otherway == "same"){
      projections_all$type = "Bi-directional Cross-protection"
    } else if(sigma_otherway == "0"){projections_all$type ="One-way Cross-protection"}
    
    projections_both <- rbind(projections_both, projections_all)
  }

projections_both_m <- melt(projections_both, id.vars= c("time","sample", "inter", "type"))

# plot some projections
projections_both_m$sample <- as.factor(projections_both_m$sample)
projections_both_m$inter<- as.factor(projections_both_m$inter)
projections_both_m <- as.data.frame(projections_both_m)
projections_both_m$date <- as.Date(projections_both_m$time, origin = run_start_2)
projections_both_m <- data.table(projections_both_m)
projections_both_m$inter <- plyr::revalue(projections_both_m$inter, c("0" = "Cross-protection: 0", 
                                                                "0.2" = "Cross-protection: 0.2",
                                                                "0.4" = "Cross-protection: 0.4",
                                                                "0.6" = "Cross-protection: 0.6",
                                                                "0.8" = "Cross-protection: 0.8",
                                                                "1" = "Cross-protection: 1"))
projections_both_m[variable == "other_all", variable := "Seasonal HCoVs"]
projections_both_m[variable == "covid_all", variable := "SARS-CoV-2"]
projections_both_m$type = forcats::fct_rev(projections_both_m$type)

# colours
cc_two <- c("navyblue","orangered2")

comp_coef <- 5
projections_both_m[date >= as.Date("2020-08-01"), scale := 2]
projections_both_m[date < as.Date("2020-08-01"), scale := 1]


PROJ_all <- ggplot(projections_both_m[scale ==1], aes(x = date, y=value,
                                                   colour = variable,
                                                   group=interaction(sample, inter, variable))) +
  geom_line()+
  geom_line(data = projections_both_m[scale ==2],aes(x = date, y=value*comp_coef,
                                                    colour = variable,
                                                    group=interaction(sample, inter, variable))) + 
  facet_grid(inter~type, switch= "y",) +
  theme_linedraw() +
  scale_y_continuous( name = "Infections (pre 2020-08-01)", 
                      sec.axis = sec_axis(~.*comp_coef, name = "Infections (post 2020-08-01)")) +
  geom_vline(xintercept =as.Date("2020-08-01"), linetype="twodash" ) +
  scale_colour_manual(values=cc_two) +
  theme( strip.text.y.left = element_blank(),
         strip.background = element_rect(colour="white", fill="white"),
         strip.text = element_text(colour = 'black',hjust = 0),
         strip.placement = "outside",
         legend.position = "bottom")+
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
  labs(colour = "Virus") +
  geom_vline(xintercept = as.Date(year_times, origin = run_start_2), colour = "black")
PROJ_all

dummy_table <- data.frame(inter = c(0,0.2,0.4,0.6,0.8,1), 
                          R0 = c(as.numeric(r_0s[sample == sample_num &
                                      interaction == 0.0, "r_0"]),
                                      as.numeric(r_0s[sample == sample_num &
                                        interaction == 0.2, "r_0"]),
                                        as.numeric(r_0s[sample == sample_num &
                                        interaction == 0.4, "r_0"]),
                                        as.numeric(r_0s[sample == sample_num &
                                        interaction == 0.6, "r_0"]),
                                        as.numeric(r_0s[sample == sample_num &
                                        interaction == 0.8, "r_0"]),
                                        as.numeric(r_0s[sample == sample_num &
                                        interaction == 1, "r_0"]))
                          )

LABELS <- ggplot(dummy_table, aes()) + 
  geom_text(aes(y = 1-inter, x =1, label = paste0("Cross-protection: ",inter))) +
  geom_text(aes(y = (1-inter)-0.05, x=1, label = paste0("R0 is: ", round(R0,1))))+
theme_void() +
  theme(plot.margin = margin(1, 0, 3, 0, "cm"))


tiff(here("figures", paste0("projections_",sample_num,".tiff")), height = 2000, width = 3200, res = 300)
grid.arrange(LABELS, PROJ_all, ncol=2, widths = c(1,4))
dev.off()
}


colnames(storage) <- c("under_5", "5_to_14", "15_to_44", "45_to_64", "65+", "time_taken", "interaction", "type")
storage <- data.table(storage)
storage_m <- melt.data.table(storage, id.vars = c("time_taken", "interaction", "type"))
storage_m$interaction <- as.factor(storage_m$interaction)
storage_m <- storage_m[type == "0"]
ONEWAY_POP <- ggplot(storage_m, aes(x = variable, y = value, fill = variable)) + 
  geom_bar(position = "dodge", stat="identity") + 
  facet_grid(interaction~time_taken) + labs(y = "PERCENT OF AGE GROUP", 
                                            title = "ONEWAY INTERACTION")

storage_m <- melt.data.table(storage, id.vars = c("time_taken", "interaction", "type"))
storage_m$interaction <- as.factor(storage_m$interaction)
storage_m <- storage_m[type == "same"]
BI_POP <- ggplot(storage_m, aes(x = variable, y = value, fill = variable)) + 
  geom_bar(position = "dodge", stat="identity") + 
  facet_grid(interaction~time_taken) + labs(y = "PERCENT OF AGE GROUP", 
                                            title = "BIDIRECTIONAL INTERACTION")


colnames(storage2) <- c("under_5", "5_to_14", "15_to_44", "45_to_64", "65+", "time_taken", "interaction", "type")
storage2 <- data.table(storage2)
storage_m <- melt.data.table(storage2, id.vars = c("time_taken", "interaction", "type"))
storage_m$interaction <- as.factor(storage_m$interaction)
storage_m <- storage_m[type == "0"]
ONEWAY_PERC <- ggplot(storage_m, aes(x = variable, y = value, fill = variable)) + 
  geom_bar(position = "dodge", stat="identity") + 
  facet_grid(interaction~time_taken) + labs(y = "PERCENT OF CASES", 
                                            title = "ONEWAY INTERACTION")

colnames(storage2) <- c("under_5", "5_to_14", "15_to_44", "45_to_64", "65+", "time_taken", "interaction", "type")
storage2 <- data.table(storage2)
storage_m <- melt.data.table(storage2, id.vars = c("time_taken", "interaction", "type"))
storage_m$interaction <- as.factor(storage_m$interaction)
storage_m <- storage_m[type == "same"]
BI_PERC <- ggplot(storage_m, aes(x = variable, y = value, fill = variable)) + 
  geom_bar(position = "dodge", stat="identity") + 
  facet_grid(interaction~time_taken) + labs(y = "PERCENT OF CASES", 
                                            title = "BIDIRECTIONAL INTERACTION")
