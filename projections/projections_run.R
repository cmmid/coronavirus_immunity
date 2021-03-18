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

projections_all <- data.frame()
# 1. get the sample 
for(i in 1:n_samples){
sample_num <- samples_to_take[i]
seasonal_sample <- trace_to_sample[sample_num,]

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

projections <- run_model_to_project(test_params  =  c(transmission_val,introductions_val) ,
                                   parameters_2020 = parameters_2020,
                                   sigma = interaction_param, 
                                   init_state_2020 = seasonal_out["init_state_2020"])

projections$sample <- sample_num
projections$inter <- sig

projections_all <- rbind(projections_all, projections)
}
print(i)}

projections_all_m <- melt(projections_all, id.vars= c("time","sample", "inter"))

# plot some projections
projections_all_m$sample <- as.factor(projections_all_m$sample)
projections_all_m$inter<- as.factor(projections_all_m$inter)
projections_all_m <- as.data.frame(projections_all_m)
projections_all_m$date <- as.Date(projections_all_m$time, origin = run_start_2)
projections_all_m <- data.table(projections_all_m)
projections_all_1 <- projections_all_m[date<"2020-08-01"]
projections_all_2 <- projections_all_m[date>="2020-08-01"]

PROJ_1 <- ggplot(projections_all_1, aes(x = date,y=value, colour = variable,
                                        group=interaction(sample, inter, variable))) + 
  geom_line() + 
  facet_grid(inter~.) +
  theme_linedraw() + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())


PROJ_2 <- ggplot(projections_all_2, aes(x = date,y=value, colour = variable,
                              group=interaction(sample, inter, variable))) + 
  geom_line() + 
  facet_grid(inter~.) +
   theme_linedraw() 

grid.arrange(PROJ_1, PROJ_2, widths = c(1,2))




