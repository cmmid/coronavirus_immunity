################################################################################
# Coronavirus cross-protection
# Author: Naomi R Waterlow
# Date: 2021-04-08
################################################################################

empty_plot <- data.frame(x = seq(from = as.Date("2014-06-09"), 
                                 to = as.Date("2020-05-31"), 
                                 by =1), 
                         y = seq(0,1400, length.out = 2184))

# Format seasonal corona data for plotting
seasonal_15_20$year_week <- as.Date(seasonal_15_20$year_week)
seasonal_hold <- seasonal_15_20[,c("variable", "year_week", "value")]
seasonal_hold <- dcast.data.table(seasonal_hold, formula = year_week ~ variable )
seasonal_hold[, All := OTHER_p0 + OTHER_p5 + OTHER_p15 + OTHER_p45 + OTHER_p45 + OTHER_p65]
seasonal_to_plot <- seasonal_hold[,c("year_week", "All")]
seasonal_to_plot$type <- "Monthly seasonal HCoV reports"

# what ratio between the two different scales
comp_coef <- 1

combine_deaths <- data.frame(
  year_week = deaths_covid$date,
  All = deaths_covid$actual_deaths /comp_coef, 
  type = "Daily Covid19 deaths"
)

# combine the data
together <- rbind(as.data.frame(seasonal_to_plot), combine_deaths)
together <- na.omit(together)

# create the time series plot
TIMESERIES <- ggplot(empty_plot, aes(x=x)) + 
  geom_line(data = together, aes(x = year_week, y = All, colour = type)) + 
  theme_linedraw() + 
  scale_colour_manual(values = c("orangered2","navyblue" ))+
   scale_y_continuous(breaks = seq(0,1400,by = 200))+
  theme(
    axis.title.y = element_text(color = "black", size=13),
    axis.title.y.right = element_text(color = "orangered2", size=13),
    legend.position = c(0.2, 0.75))  + 
  labs(x = "Date", title = "A", colour = "Time Series", y = "Count")

# calculate proportion in each age group
tot <- seasonal_hold[,-1]
tot_sum <- apply(na.omit(tot), 2, sum)
to_plot <- data.frame(ages = c("0-4", "5-14", "15-44","45-64", "65+"), 
                      percent_cases = c(tot_sum["OTHER_p0"]/tot_sum["All"],
                                        tot_sum["OTHER_p5"]/tot_sum["All"],
                                        tot_sum["OTHER_p15"]/tot_sum["All"],
                                        tot_sum["OTHER_p45"]/tot_sum["All"],
                                        tot_sum["OTHER_p65"]/tot_sum["All"]))
to_plot$ages <- factor(to_plot$ages, levels =c("0-4", "5-14", "15-44","45-64", "65+") )
# create the age hcov plot
HCOV_AGES <- ggplot(to_plot, aes(x = ages, y = percent_cases )) + 
  geom_bar(stat="identity", fill= "navyblue") +
  theme_linedraw() + 
  labs(x = "Age group", y = "Proportion of seasonal HCoV cases", title = "B")

# add the death data
to_plot_covid <- data.table(
  data = c(0.007,0.038,0.027,0.03,0.077, rep(NA,11)),
  ages = factor(c("0-4", "5-9","10-14", "15-19", "20-24", "25-29", "30-34", "35-39", 
                  "40-44", "45-49", "50-54", "55-59", "60-64" , "65-69", "70-74", "75+"), 
                levels = c("0-4", "5-9","10-14", "15-19", "20-24", "25-29", "30-34", "35-39", 
                           "40-44", "45-49", "50-54", "55-59", "60-64",  "65-69", "70-74", "75+")),
  lower = c(0.000,0.002,0.00,0.001,0.018,rep(NA,11)), 
  upper = c(0.058,0.101,0.077,0.084,0.159,rep(NA,11))
)

# # make the death covid plot
# COVID_AGES <- ggplot(na.omit(to_plot_covid), aes(x = ages, y= data) ) + 
#   geom_pointrange(aes(ymin = lower, ymax = upper), colour="orangered2") +
#   theme_linedraw() + 
#   labs(y = "Proportion positive for SARS-CoV-2 Antibodies", x = "Age groups",
#        title = "D")

#make the empty plot (for adding in model)
EMPTY <- ggplot() +
  theme(panel.background = element_blank())

# save the combined plot
tiff(here("figures","intro.tiff"), height = 2000, width = 3200, res = 300)

grid.arrange(TIMESERIES, HCOV_AGES, EMPTY, layout_matrix = rbind(c(1,1,1), 
                                                                      c(2,3,3)))
dev.off()



seasonal_hold[,All := NULL]
seasonal_hold <- na.omit(seasonal_hold)
sh_m <- melt.data.table(seasonal_hold, id.vars = c("year_week"))
sh_m[variable == "OTHER_p0", variable:= "Age 0 - 4"]
sh_m[variable == "OTHER_p5", variable:= "Age 5 - 14"]
sh_m[variable == "OTHER_p15", variable:= "Age 15 - 44"]
sh_m[variable == "OTHER_p45", variable:= "Age 45 - 64"]
sh_m[variable == "OTHER_p65", variable:= "Age Age 65+"]

tiff(here("figures","S_data.tiff"), height = 2000, width = 3200, res = 300)

ggplot(sh_m, aes( x = year_week, y = value)) + 
  labs(x = "Date", y = "Reported Seasonal Coronavirus Infections") +
  geom_line() + 
  facet_grid(variable~.) +
  theme_linedraw() +
  scale_x_date(breaks =c(as.Date("2015-01-01"),as.Date("2016-01-01"),
                         as.Date("2017-01-01"),as.Date("2018-01-01"),
                         as.Date("2019-01-01"),as.Date("2020-01-01")),
               labels = c("2015","2016","2017","2018","2019","2020"))
             
dev.off()
