trace_to_sample <- data.table(trace_to_sample)

trace_to_sample[(seasonal_R0 >5) &
                  (seasonal_R0 <= 6.5), section := "5_to_6.5"]
trace_to_sample[(seasonal_R0 >6.5) &
                  (seasonal_R0 <= 8), section := "6.5_to_8"]
#trace_to_sample[(seasonal_R0 >6) &
#                  (seasonal_R0 <= 6.5), section := "6_to_6.5"]

trace_to_sample_m <- melt.data.table(trace_to_sample)
trace_to_sample_m$type <- "excluding_2014"

trace_to_sample_m[variable == "time_steps",] <- NA
trace_to_sample_m <- na.omit(trace_to_sample_m)


estimates$X1 <- as.numeric(estimates$X1)
estimates$X2 <- as.numeric(estimates$X2)
estimates$X3 <- as.numeric(estimates$X3)
colnames(estimates) <- c("X1", "X2", "X3", "variable")
library(ggstance)
ggplot(trace_to_sample_m, aes(x = value, fill = "navyblue")) + 
  geom_density() + 
  facet_wrap(variable~., scales = "free")

ALL <- ggplot(trace_to_sample_m, aes(x = value)) + 
  geom_histogram(bins =500, position = "stack", fill="deepskyblue") + 
  facet_wrap(variable~., scales = "free") + 
  geom_pointrangeh(data = estimates, aes(xmin = X1, xmax = X3, x = X2, y = 1)) + 
  labs(x = "Parameter value", y = "Count", title = "A: Posterior estimates") + 
  theme_linedraw(
  )


key_subset <- trace_to_sample_m[variable=="log_liklihood" |
                                  variable == "seasonal_R0", ]

estimates_subset <- estimates[c(2),]

SUSBET <- ggplot(key_subset, aes(x = value)) + 
  geom_histogram(bins =500, position = "stack", aes(fill = section)) + 
  facet_wrap(variable~., scales = "free", ncol = 1) + 
  geom_pointrangeh(data = estimates_subset, aes(xmin = X1, xmax = X3, x = X2, y = 1)) + 
  labs(x = "Parameter value", y = "Count", title = "B: Posteriors by mode", 
       fill = "R0 subset") + 
  theme_linedraw()


COR <- ggplot(trace_to_sample, aes(x = seasonal_R0, y = waning_duration, colour = log_liklihood)) + 
  geom_point() + 
  theme_linedraw() + 
  labs(x = "Seasonal R0", y = "Waning duration (days)", 
       colour = "log likelihood", title= "C: Correlation")




tiff(here("figures","betasesitivity.tiff"), height = 2000, width = 3200, res = 300)
grid.arrange(ALL, SUSBET, COR, layout_matrix =rbind(c(1,1,2),
                                                    c(1,1,2),
                                                    c(1,1,3)))
dev.off()

tiff(here("figures","2014sesitivity.tiff"), height = 2000, width = 3200, res = 300)
ALL
dev.off()

