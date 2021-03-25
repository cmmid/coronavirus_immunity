################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################
# extra plots

# values to rest
R0_trials <- seq(1,10, by = 0.2)
waning_trials <- seq(0,10, by = 0.2)*364
store_ll <- data.frame()
# choose sample from trace
sample_num <- samples_to_take[i] #in paper is sample 74904
seasonal_sample <- trace_to_sample[sample_num,]

# loop pver different parameters
for(required_r0 in R0_trials){
  for(required_waning in waning_trials){
    
    parameter_guesses <- c(waning_day = required_waning,
                           seasonal_R0 = required_r0,
                           seasonal_sample["seasonal_reported_1"],
                           seasonal_sample["seasonal_reported_2"],
                           seasonal_sample["seasonal_reported_3"],
                           seasonal_sample["seasonal_reported_5"],
                           seasonal_sample["seasonal_amplitude"],
                           seasonal_sample["phi"])
    #calculate the likelihood
    ll <- run_model_get_logliks_seasonalonly(parameter_guesses = parameter_guesses, 
                                             bb = "binomial", model_type ="SEIR")
    #store the likelihood
    temp <-c(r0 = required_r0, waning = required_waning, ll = ll)
    store_ll <- rbind(store_ll, temp)
    
  }
} 

#format
colnames(store_ll) <- c('r0', 'waning', 'll')
# choose zooming scales
zoom1 <- c(4,7,1800,3500)
zoom2 <- c(5.2,6,2300,2900)

#make the plots
 A <- ggplot(store_ll, aes(x = r0, y = waning, fill = ll)) + 
  geom_tile() + 
  scale_fill_viridis(option = "magma", na.value = "lightgrey") + 
  geom_point(aes(x = seasonal_sample$seasonal_R0, y = seasonal_sample$waning_day), 
             shape = 3) +#+ facet_zoom(xlim=c(4,8), ylim =c(2000,3500))
   geom_rect( mapping=aes(xmin=zoom1[1], xmax=zoom1[2], ymin=zoom1[3], ymax=zoom1[4]), color="black",
              alpha=0.5, fill = NA) + 
   labs( x = "R0", y = "",fill = "LL", title = "A: full range") + 
 theme_linedraw() #+ 

 
store_ll <- data.table(store_ll)
store_ll_sub <- store_ll[r0>zoom1[1] & r0 <zoom1[2] &waning<zoom1[4] & waning>zoom1[3]]
B <- ggplot(store_ll_sub, aes(x = r0, y = waning, fill = ll)) + 
  geom_tile() + 
  scale_fill_viridis(option = "magma", na.value = "lightgrey") + 
  geom_point(aes(x = seasonal_sample$seasonal_R0, y = seasonal_sample$waning_day), 
             shape = 3) +
  geom_rect( mapping=aes(xmin=zoom2[1], xmax=zoom2[2], ymin=zoom2[3], ymax=zoom2[4]), color="black",
             alpha=0.5, fill = NA) + 
  labs(y = "Waning duration (days)", x = "R0", fill = "LL", title = "B: zoom one")+
  theme_linedraw()

store_ll_sub2 <- store_ll[r0>zoom2[1] & r0 <zoom2[2] &waning<zoom2[4] & waning>zoom2[3]]
C <- ggplot(store_ll_sub2, aes(x = r0, y = waning, fill = ll)) + 
  geom_tile() + 
  scale_fill_viridis(option = "magma", na.value = "lightgrey") + 
  geom_point(aes(x = seasonal_sample$seasonal_R0, y = seasonal_sample$waning_day), 
             shape = 3) + 
  labs( x = "R0", fill = "LL", title = "C: zoom two", y  ="") +
  theme_linedraw()


tiff(here("figures","heatmap.tiff"), height = 3200, width = 2000, res = 300)

grid.arrange(A,B,C)

dev.off()

