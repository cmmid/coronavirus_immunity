################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################
set.seed(19930625)
##### SETUP #######
run_start_2 <- as.Date("2020-01-01") # When to start the smulation
run_end <- as.Date("2020-06-01") # When to end teh simulation
length_to_run_2 <- (run_end - run_start_2)[[1]] #total days to run
mobility_date <- "2020-02-21" # date teh google mobility data starts
covid_date <- as.Date("2020-02-10") # date that covid is first introduced
covid_change_time <- as.Date("2020-03-23") #date of lockdown (i.e. schools closing)
serology_time <- as.Date("2020-05-31") # time point of seology.
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
transmission <- data.table()
introductions <- data.table()
ll_values <- data.table()

sig_list <- c(0,0.2,0.4,0.6,0.8,1)

n_samples <- 50
samples_to_take <- sample(x=1:dim(trace_to_sample)[1],size=n_samples,replace = T)
i <- 1
while(i <= n_samples){
 
  for(sig in sig_list){
    
    # need to run a loop or an apply over the number of samples required.
    # take sample of seasonal parameters
    sample_num <- samples_to_take[i]
    seasonal_sample <- trace_to_sample[sample_num,]
    
    # run the seaonal model in order to get init start for the covid model
    seasonal_out <- run_seasonal_for_init(seasonal_sample, model_type = "SEIR")
    parameters_2020 <- update_parameters(seasonal_out[["parameters_2019"]])
    
    # loop here over different interaction parameters
    interaction_param <- sig
    
    transmission_fit <- optim(par = c(parameters_2020$beta_other,1), # use as the transmission parameter (convert to R0 after) 
                              fn = optimise_ll_deaths, 
                              method = "L-BFGS-B",
                              lower = c(0.001,0.1), 
                              upper = c(0.3,4.5),
                              control=list(parscale=c(0.0001,0.1)),
                              parameters_2020 = parameters_2020, 
                              sigma = interaction_param,
                              init_state_2020 = seasonal_out["init_state_2020"])
    # run maximum likelihood estimation 
    # run covid model
    # ll_deaths
    # run the model with the ML estimate parameters

    output_list <- run_model_on_sample(test_params  =  transmission_fit$par,
                                       parameters_2020 = parameters_2020,
                                       sigma = interaction_param, 
                                       init_state_2020 = seasonal_out["init_state_2020"])

    if(any(is.na(output_list$r_0)) |
           any(is.na(output_list$r_eff)) |
               any(is.nan(output_list$model_deaths4total))
           ){
      rerun_new_sample <- TRUE
      samples_to_take[i] <- samples_to_take[i] +1
      print("about to jump")
      break}
    rerun_new_sample <- FALSE
    # store the death output with the output from other samples
    temp_ll_values <- data.frame(transmission_fit$value, sample_num, interaction_param)
    
    ll_values <- rbind(ll_values, temp_ll_values)
    
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
    
    transmission_temp <- data.table(transmission_val = transmission_fit$par[1])
    transmission_temp$sample <- sample_num
    transmission_temp$interaction <- interaction_param

    transmission <- rbind(transmission,transmission_temp)
    
    introductions_temp <- data.table(intro_val = transmission_fit$par[2])
    introductions_temp$sample <- sample_num
    introductions_temp$interaction <- interaction_param
    
    introductions <- rbind(introductions,introductions_temp)
    print(sig)
  }
    if(rerun_new_sample == F){
      i <- i + 1
    }
  print(i)
}

#colours for all
cc <- scales::seq_gradient_pal("darkmagenta",  "darkorange", "Lab")(seq(0,1,length.out=6))
saveRDS(transmission, file ="transmission.RDS")
saveRDS(introductions, file ="introductions.RDS")
saveRDS(sero, file ="sero.RDS")
saveRDS(r_effs, file ="r_effs.RDS")
saveRDS(r_0s, file ="r_0s.RDS")
saveRDS(model_deaths, file ="model_deaths.RDS")
saveRDS(samples_to_take, file ="samples_to_take.RDS")


# create the fit plot
model_deaths$sample <- factor(model_deaths$sample)
model_deaths$interaction <- factor(model_deaths$interaction)
model_deaths$interaction <- plyr::revalue(model_deaths$interaction, c("0" = "Cross-protection: 0", 
                                    "0.2" = "Cross-protection: 0.2",
                                    "0.4" = "Cross-protection: 0.4",
                                    "0.6" = "Cross-protection: 0.6",
                                    "0.8" = "Cross-protection: 0.8",
                                    "1" = "Cross-protection: 1"))

  FIT_PLOT <- ggplot() + 
  geom_point(data=model_deaths, aes(x = date, y = total, colour = interaction )) +
  facet_wrap(interaction~., ncol=1)+
scale_colour_manual(values=cc) +
   geom_point(data = deaths_covid, aes(x = date, y = actual_deaths), size = 1) + 
  geom_vline(xintercept = covid_date, colour = "red" ) + 
  labs(x = "Date", title = "", y = "Number of deaths", colour = "strength of protection") + 
   lims(x = c(as.Date("2020-02-01"), as.Date("2020-05-31")))  + 
   theme_linedraw() +
   theme(strip.text.y = element_text(angle = 0), 
         strip.background = element_rect(colour="white", fill="white"),
         strip.text = element_text(colour = 'black',hjust = 0))
 
# create the R0 plot
r_0s_summary <- r_0s[,mean(r_0), by=interaction]
colnames(r_0s_summary) <- c("Protection", "Mean_R0")
r_0s_summary$Mean_R0 <- round(r_0s_summary$Mean_R0,2)
r_0s_summary$Minimum_R0 <- round(r_0s[,min(r_0), by=interaction][,2],2)
r_0s_summary$Maximum_R0 <- round(r_0s[,max(r_0), by=interaction][,2],2)
r_0s_summary$Protection <- as.factor(r_0s_summary$Protection)

# R0 <- ggplot(r_0s_summary, aes(x = Protection, y = Mean_R0, fill = Protection)) + 
#   geom_bar(stat="identity") + 
#   theme_linedraw() + 
#  scale_fill_manual(values=cc)+
#   geom_errorbar(aes(ymin = Minimum_R0, ymax = Maximum_R0))+
#   labs(x = "Cross-protection", y = "Mean R0", title = "C") + 
#   theme(legend.position = "none")
# create the r_effective plot
colnames(r_effs) <- c("r_eff", "timestep", "sample", "interaction")

r_effs$sample <- factor(r_effs$sample)
r_effs <- as.data.frame(r_effs)
r_effs$date <- as.Date(r_effs$timestep, origin = run_start_2)
r_effs$interaction <- factor(r_effs$interaction)

r_0s$interaction <- forcats::fct_rev(as.factor(r_0s$interaction))
R0 <- ggplot(r_0s, aes(x = interaction, y = r_0, colour = interaction)) + 
  geom_jitter() + 
  theme_linedraw() + 
  scale_colour_manual(values=rev(cc))+
  scale_y_continuous(breaks = seq(from=0,30,by=5), 
                     limits = c(0,30))+
  labs(x = "Cross-protection", y = "R0 (points), Reff (lines)", title = "A") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size =12), 
        plot.margin = margin(t = 3,r=2,b= 0, l=2, unit = "pt")) + 
  coord_flip()  + 
  geom_line(data = r_effs,aes(x = interaction, y = r_eff, colour = interaction))

R0

# create the serology plot  
sero$interaction <- as.factor(sero$interaction)
sero$ages <- forcats::fct_rev(sero$ages)
# sero[, source := 2]
# sero[ages == "0-4" | ages == "5-9" | ages == "10-14" |
#        ages == "15-19" | ages == "20-24", source := 1]
# sero$source <- as.factor(sero$source)
# dup_point <- data.frame(x = "20-24", y = 0.102)
#sero$model <- sero$model * unlist(sero[5,"model"])

SERO_PLOT <- ggplot(sero, aes(x= ages)) + 
  # geom_errorbar(aes(y=data, ymin=lower, ymax=upper)) +
  scale_colour_manual(values=cc) +
  labs(y = "Proportion positive", title = "B", x = "Age group", colour = "Cross-protection",
       shape = "Source")+ 
  theme_linedraw() + 
#  lims(y = c(0,0.2)) +
  geom_point(aes(y=model, colour = interaction, group = sample), alpha=0.8) + 
  theme(legend.position = "bottom") + 
  guides(shape = guide_legend(override.aes = list(size = 0.5)))+
guides(color = guide_legend(override.aes = list(size = 1.5)))+
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8))+ 
   theme( legend.margin=margin(0.5,0.5,0.5,0.5),
          legend.box.margin=margin(-10,-10,-10,-10),
          legend.box.just = "center", 
          legend.box = "horizontal",
          legend.position = c(0.8, 0.75),
          legend.key.height = unit(0.3, "cm")) +
  geom_line(aes(x = ages, y = model, group = interaction(interaction,sample), colour = interaction))  + 
  geom_pointrange(data = sero_data, aes(x = age_group, y = data, ymin= lower, ymax = upper, shape = type), 
                  position = position_dodge(width=0.3))
  # geom_pointrange(data = sero_out, aes(x = age_group, y = mean,ymin = lower, ymax=upper, 
  #                                              group = source, shape = source),
  #                 position=position_dodge(width=c(0.4))) +
  # scale_shape_manual(values=1:7)
 # geom_point(aes(y=data,shape = source), size=3)# +
#  geom_point(data = dup_point, aes(x = x, y = y), shape = "triangle", size = 3)
  
# plot the combined plots
tiff(here("figures","covid_sims.tiff"), height = 2000, width = 3200, res = 300)

grid.arrange( R0,SERO_PLOT, layout_matrix = rbind(c(2,2,2), 
                                                            c(2,2,2),
                                                            c(2,2,2),
                                                            c(3,3,3), 
                                                            c(3,3,3), 
                                                            c(3,3,3),
                                                            c(3,3,3)))
dev.off()

tiff(here("figures","fit_to_deaths.tiff"), height = 2000, width = 3000, res = 300)

FIT_PLOT

dev.off()

tiff(here("figures","reffs.tiff"), height = 2000, width = 3200, res = 300)

ggplot(r_effs, aes(x=date, y = r_eff, group=sample, colour = interaction )) + 
  geom_line(aes(group=interaction(sample, interaction))) + 
  geom_hline(yintercept = 2.25) + 
  geom_hline(yintercept = 3.75) + 
  geom_vline(xintercept = as.Date((covid_date-run_start_2)[[1]], origin = run_start_2), colour = "red") + 
  labs(x = "Date", y = "R-effective", colour = "Cross-protection") + 
  theme_linedraw() + 
  scale_colour_manual(values=cc) + 
  scale_y_continuous(breaks = seq(0,5.5, by = 0.5), 
                     limits = c(0,5.5))+
  theme(strip.text.y = element_text(angle = 0))

dev.off()
  
