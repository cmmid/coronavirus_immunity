################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################

##### SETUP #######
run_start_2 <- as.Date("2020-01-01") # When to start the smulation
run_end <- as.Date("2020-06-01") # When to end teh simulation
length_to_run_2 <- (run_end - run_start_2)[[1]] #total days to run
mobility_date <- "2020-02-21" # date teh google mobility data starts
covid_date <- as.Date("2020-02-10") # date that covid is first introduced
covid_change_time <- as.Date("2020-03-23") #date of lockdown (i.e. schools closing)
serology_time <- as.Date("2020-05-31") # time point of seology. CHECK
#trace from which to take samples. Should already be formated appropiately. 
load(here("fitting_seasonal/analysis", "trace_to_sample.Rdata"))
# change to correct names and reporting as log odds
colnames(trace_to_sample) <- c("ll","waning_day", "seasonal_R0", "seasonal_reported_1", 
                                 "seasonal_reported_2", "seasonal_reported_3", 
                                 "seasonal_reported_5", "seasonal_amplitude", 
                                 "phi", "step")
trace_to_sample$seasonal_reported_1 <- log(trace_to_sample$seasonal_reported_1/(1-trace_to_sample$seasonal_reported_1))
trace_to_sample$seasonal_reported_2 <- log(trace_to_sample$seasonal_reported_2/(1-trace_to_sample$seasonal_reported_2))
trace_to_sample$seasonal_reported_3 <- log(trace_to_sample$seasonal_reported_3/(1-trace_to_sample$seasonal_reported_3))
trace_to_sample$seasonal_reported_5 <- log(trace_to_sample$seasonal_reported_5/(1-trace_to_sample$seasonal_reported_5))
trace_to_sample <- data.frame(trace_to_sample)



# storage
model_deaths <- data.table()
r_0s <- data.table()
r_effs <- data.table()
sero <- data.table()

sig_list <- c(0,0.2,0.4,0.6,0.8,1)

for(i in 1:100){
  for(sig in sig_list){
    # need to run a loop or an apply over the number of samples required.
    # take sample of seasonal parameters
    sample_num <- sample(x=1:dim(trace_to_sample)[1],size=1,replace = T)
    seasonal_sample <- trace_to_sample[sample_num,]
    
    # run the seaonal model in order to get init start for the covid model
    seasonal_out <- run_seasonal_for_init(seasonal_sample, model_type = "SEIR")
    parameters_2020 <- update_parameters(seasonal_out[["parameters_2019"]])
    
    # loop here over different interaction parameters
    interaction_param <- sig
    
    transmission_fit <- optim(par = parameters_2020$beta_other, # use as the transmission parameter (convert to R0 after) 
                              fn = optimise_ll_deaths, 
                              method = "Brent",
                              lower = 0.01, 
                              upper =  0.3631551, 
                              parameters_2020 = parameters_2020, 
                              sigma = interaction_param,
                              init_state_2020 = seasonal_out["init_state_2020"])
    # run maximum likelihood estimation 
    # run covid model
    # ll_deaths
    # run the model with the ML estimate parameters

    output_list <- run_model_on_sample(transmission_rate =  transmission_fit$par,
                                       parameters_2020 = parameters_2020,
                                       sigma = interaction_param, 
                                       init_state_2020 = seasonal_out["init_state_2020"])
    
   
    # store the death output with the output from other samples
    temp_deaths <- output_list$model_deaths
    temp_deaths$sample <- sample_num
    temp_deaths$interaction <- interaction_param
    
    model_deaths <- rbind(model_deaths, temp_deaths)
    
    r_0s_temp <- data.table(r_0 = output_list$r_0)
    r_0s_temp$sample <- sample_num
    r_0s_temp$interaction <- interaction_param
    
    r_0s <- rbind(r_0s,r_0s_temp)
    
    r_effs_temp <- data.table(r_eff = output_list$r_eff)
    r_effs_temp$sample <- sample_num
    r_effs_temp$interaction <- interaction_param
    
    r_effs <- rbind(r_effs,r_effs_temp)
    
    sero_temp <- data.table(output_list$serology)
    sero_temp$sample <- sample_num
    sero_temp$interaction <- interaction_param
    
    sero <- rbind(sero,sero_temp)
    
    
  }
  print(i)
}
# store the serology output with the output from other samples

# store the Reffective output with the output from other samples




model_deaths$sample <- factor(model_deaths$sample)

 FIT_PLOT <- ggplot() + 
  geom_point(data=model_deaths, aes(x = date, y = total, colour = interaction ),alpha=0.5) +
  facet_grid(interaction~.)+
   scale_color_gradient(low = "deepskyblue", high = "royalblue4") +
  geom_point(data = deaths_covid, aes(x = date, y = actual_deaths), size = 1) + 
  geom_vline(xintercept = covid_date, colour = "red" ) + 
  labs(x = "Date", y = "Number of deaths", colour = "strength of protection") + 
   lims(x = c(as.Date("2020-02-01"), as.Date("2020-05-31")))  + 
   theme_linedraw()
 
colnames(r_effs) <- c("r_eff", "timestep", "sample", "interaction")

r_effs$sample <- factor(r_effs$sample)
r_effs <- as.data.frame(r_effs)
r_effs$date <- as.Date(r_effs$timestep, origin = run_start_2)

tiff(here("figures","reffs.tiff"), height = 2000, width = 3200, res = 300)

ggplot(r_effs, aes(x=date, y = r_eff, group=sample, colour = interaction )) + 
  geom_line(aes(group=interaction(sample, interaction))) + 
  geom_hline(yintercept = 2.25) + 
  geom_hline(yintercept = 3.75) + 
  geom_vline(xintercept = as.Date((covid_date-run_start_2)[[1]], origin = run_start_2), colour = "red") + 
  labs(x = "Date", y = "R-effective", colour = "Strength of protection") + 
  theme_linedraw()

dev.off()

r_0s_summary <- r_0s[,mean(r_0), by=interaction]
colnames(r_0s_summary) <- c("Protection", "Mean_R0")
r_0s_summary$Mean_R0 <- round(r_0s_summary$Mean_R0,2)
r_0s_summary$Minimum_R0 <- round(r_0s[,min(r_0), by=interaction][,2],2)
r_0s_summary$Maximum_R0 <- round(r_0s[,max(r_0), by=interaction][,2],2)

tiff(here("figures","fit.tiff"), height = 2000, width = 3200, res = 300)

grid.arrange(FIT_PLOT, tableGrob(r_0s_summary, rows = NULL, theme =  gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.75)),
  colhead = list(fg_params=list(cex = 0.75)),
  rowhead = list(fg_params=list(cex = 0.75)))), layout_matrix = rbind(1,1,1,2))
dev.off()

SERO_PLOT <- ggplot(sero, aes(x= ages)) + 
  geom_pointrange(aes(y=data, ymin=lower, ymax=upper)) +
  geom_point(aes(y=model, colour = interaction, group = sample), alpha=0.8) + 
  scale_color_gradient(low = "deepskyblue", high = "royalblue4") +
   labs(y = "Proportion positive", x = "Age group", colour = "Strength of protection")+ 
  theme_linedraw() + 
  lims(y = c(0,0.2))

tiff(here("figures","sero.tiff"), height = 2000, width = 3200, res = 300)
SERO_PLOT
dev.off()

  
# Used to calculate the upper limit of R0
calc_beta_SEIR(0.2, 18, parameters_2020, 1)
  
