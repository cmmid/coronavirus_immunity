# Analyse the traces and create figures for paper. 

######## LOAD THE TRACES ######
n_chains <- 16

# load the traces
load(here("fitting_seasonal/analysis","SEIR_PT_newlims5_5_2021-02-28.Rdata"))
trace_1_all <- total_trace
trace_temp <- lapply(total_trace[1:n_chains], function(i) label_trace(i))
timesteps <- 1:dim(total_trace[[1]])[1]
trace_1 <- lapply(trace_temp, function(i) format_trace(i, thin=, burnin=12000,
                                                       keep_all = F))

load(here("fitting_seasonal/analysis","SEIR_AWS_tracelims2_2020_03.Rdata"))
trace_2_all <- total_trace
trace_temp <- lapply(total_trace[1:n_chains], function(i) label_trace(i))
timesteps <- 1:dim(total_trace[[1]])[1]
trace_2 <- lapply(trace_temp, function(i) format_trace(i, thin=, burnin=12000,
                                                       keep_all = F))

######### CONVERGENCE DIAGNOSTIC #######

compare_mcmc<-mcmc.list(mcmc(trace_1[[1]][1:22000,2:9]), 
                        mcmc(trace_2[[1]][1:22000,2:9]))

gelman_rubin <-gelman.diag(compare_mcmc)

print(paste0("The Gelman Rubin statistic is ", gelman_rubin[2]))

####### PLOT ONE MULTI-TRACE ######

trace_temps <- trace_1_all[[length(trace_1_all)-2]]
temperatures <- 1/trace_temps[dim(trace_temps)[1],]
timesteps <- c(1:dim(trace_1[[1]])[1])
# thinning 0 here as already thnned above
trace_edited <- lapply(trace_1, function(i) format_trace(i, thin=1, burnin=1, keep_all = T))
trace_edited <- Map(cbind, trace_1, temp = temperatures)

trace_combined <- do.call(rbind.data.frame, trace_edited)
trace_combined <- melt(trace_combined, id=c("time_steps", "temp"))

TRACE_PLOTS <- plot_trace(trace_combined)
DENS_PLOTS <- plot_density(trace_combined)

tiff(here("figures","Multi_trace_plot.tiff"), height = 2000, width = 3200, res = 300)

TRACE_PLOTS

dev.off()

####### COMBINE TRACES ####### (and save)

trace_to_sample <- rbind(trace_1[[1]][], trace_2[[1]])
save(trace_to_sample, file = here("fitting_seasonal/analysis","trace_to_sample.Rdata"))

####### PLOT POSTERIOR ######
trace_using <- data.table(trace_to_sample)

trace_using[,step:= 1:dim(trace_using)[1]]
ggtrace <- melt.data.table(trace_using, id.vars= c("step"))
colnames(ggtrace) <- c("step", "parameter", "value")
ggtrace <- ggtrace[parameter != "time_steps"]

tiff(here("figures","posterior.tiff"), height = 2000, width = 3200, res = 300)

ggplot(ggtrace) + 
  geom_density(aes(x = value)) + 
  facet_wrap(parameter~., scales = "free") + 
  theme_linedraw() + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  labs(x = "Parameter value", y = "Density")

dev.off()
######## RBINOM SAMPLING ######

RBINOM <- plot_rbinom(samples = 100, 
                      trace_period = c(1:dim(trace_to_sample)[1]),
                      trace_dt = trace_to_sample,
                      model_type = "SEIR")

tiff(here("figures","rbinom.tiff"), height = 2000, width = 3200, res = 300)
RBINOM
dev.off()

######## QUANTILES ########
print("R0 quantiles")
print(quantile(trace_to_sample$waning_duration, probs = c(0.025, 0.5, 0.975))/364)
print("waning quantiles")
print(quantile(trace_to_sample$seasonal_R0, probs = c(0.025, 0.5, 0.975)))
