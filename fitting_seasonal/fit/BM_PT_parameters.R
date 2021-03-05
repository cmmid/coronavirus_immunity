#BM_PT_parameters etc.
age_groups <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)

###### load seasonal data and format ######

# Load in the age_specific data
all_data_ages <- data.table(read.csv(here("fitting_seasonal/fit","all_data_ages.csv"),
                                     stringsAsFactors =F))
# melt
all_data_ages_m <- melt.data.table(all_data_ages, id.vars=c("Year", "startweek",
                                                            "Age", "year_week", 
                                                            "timestep"))
all_data_ages_m$value <- as.numeric(all_data_ages_m$value)

# if want to convert to under 5 run these two lines
all_data_ages_m[Age == "Under 1", Age := "Under5"]
all_data_ages_m[Age == "1-4 Y", Age := "Under5"]
# subset and combine the Under 5 age groups
all_data_ages_m <- all_data_ages_m[,c("Age", "year_week", "variable", "value")]
all_data_ages_m<-all_data_ages_m[, lapply(.SD, function(x) sum(x)), by = 
                                   .(Age, year_week,variable)]
all_data_ages_m$Age <- factor(all_data_ages_m$Age, 
                              levels = c("Under5", "5-14 Y", "15-44 Y",
                                         "45-64 Y", "65 Y+"))

corona_15_20<- all_data_ages_m[variable =="Coronavirus" & 
                                 year_week < as.Date("2020-03-01")] 

# subset the 2018/2019/2020 from lowest point, in the first year
corona_19_20 <- all_data_ages_m[variable =="Coronavirus" ]
# use this to calculate the lowest point in the season
corona_19_20[,sum(value), by = year_week]
corona_15_20[,sum(value), by = year_week]
#lp_15 = as.Date("2014-09-01 ")
lp = as.Date("2018-09-03")
lp_15 = as.Date("2002-01-01")

#corona_19_20 <- corona_19_20[year_week > as.Date("2019-05-01")]
seasonal_19_20 <- corona_19_20[ year_week >= lp_15,]
seasonal_19_20[ Age == "Under5",variable := "OTHER_p0"]
seasonal_19_20[ Age == "5-14 Y",variable := "OTHER_p5"]
seasonal_19_20[ Age == "15-44 Y",variable := "OTHER_p15"]
seasonal_19_20[ Age == "45-64 Y",variable := "OTHER_p45"]
seasonal_19_20[ Age == "65 Y+",variable := "OTHER_p65"]
seasonal_19_20<-seasonal_19_20[year_week < as.Date("2020-06-08")]

seasonal_15_20 <- corona_15_20[ year_week >= lp_15,]
seasonal_15_20[ Age == "Under5",variable := "OTHER_p0"]
seasonal_15_20[ Age == "5-14 Y",variable := "OTHER_p5"]
seasonal_15_20[ Age == "15-44 Y",variable := "OTHER_p15"]
seasonal_15_20[ Age == "45-64 Y",variable := "OTHER_p45"]
seasonal_15_20[ Age == "65 Y+",variable := "OTHER_p65"]

to_exclude1 <- which(seasonal_15_20[,"year_week"] == "2017-04-03")
to_exclude2 <- which(seasonal_19_20[,"year_week"] == "2017-04-03")
seasonal_dates_15 <- unique(seasonal_15_20$year_week)

seasonal_15_20[c(to_exclude1),"value"] <- NA
seasonal_19_20[c(to_exclude2), "value"] <- NA



###### Other parameters #####
run_end <- as.Date("2020-07-01")
length_to_run <-  (run_end - lp_15 )[[1]]

# Population numbers: source - iukpppsumpop18 ons pop 2019 mid year estimate
pop_numbers <- c(3465179.00, 3721990.00, 3535065.00, 3262613.00, 3690265.00,
                 4009669.00, 4000908.00, 3918562.00, 3583853.00, 3914884.00,
                 4127941.00, 3888131.00, 3304688.00, 2978882.00, 2958612.00,
                 (2067471.00+  1529682.00+ 933656.00+ 547789.00))
pop_numbers_5year <- c(3465179.00, 3721990.00, 3535065.00, 3262613.00, 3690265.00,
                       4009669.00, 4000908.00, 3918562.00, 3583853.00, 3914884.00,
                       4127941.00, 3888131.00, 3304688.00, 2978882.00, 2958612.00,
                       2067471.00,  1529682.00, 933656.00, 547789.00)
# Births: 1.79 births per woman per lifetime. so in model times by population
births <- 640370/(52*7)
deaths <- births
# Ageing: in years
ageing <- c(rep(1/5,length(pop_numbers)-1)/(52*7),0)


fixed_parameters <- c(
  covid_incubation = 3, # 3.5 fixed
  covid_gamma = 5, # fixed
  seasonal_incubation = 2.5, #fixed 2.5
  seasonal_gamma = 5, # fixed
  covid_inf_to_symp = 21/2,
  covid_symp_to_death = 21/2,
  s_inf_to_symp = 2, 
  s_symp_to_report= 3,
  death_rate_1 =  0, #0.00004,
  death_rate_2 =  0, #0.00004,
  death_rate_3 = 0,  #0.000253,
  death_rate_4 =  0, #0.0049,
  death_rate_5 = 0, #0.131,
  mobility_start = 10000000,
  covid_date = 45 ,# where 1 is "2020-01-01"
  # covid_intro_adult = 6, 
  covid_R0 = 0, 
  sig_R = 0,
  sig_solid_R = 0,
  social_distancing = 0
)

# create the UK contact matrix
# Load polymod data
contacts_all <- unname(as.matrix(read.csv2(here("fitting_seasonal/fit","contacts_all.csv"), stringsAsFactors = F)))[,2:17]
contacts_hh <- unname(as.matrix(read.csv2(here("fitting_seasonal/fit","contacts_hh.csv"), stringsAsFactors = F)))[,2:17]
contacts_school <- unname(as.matrix(read.csv2(here("fitting_seasonal/fit","contacts_school.csv"), stringsAsFactors = F)))[,2:17]
contacts_other <- unname(as.matrix(read.csv2(here("fitting_seasonal/fit","contacts_other.csv"), stringsAsFactors = F)))[,2:17]


pop_params_base <- list(
  age_groups = age_groups,
  pop_numbers = pop_numbers,
  births = births,
  deaths = deaths,
  ageing = ageing, 
  contacts_all = contacts_all,
  contacts_other = contacts_other, 
  contacts_school = contacts_school, 
  contacts_hh = contacts_hh,
  contacts_reduc = 0)


rho <- 1 # temporary interaction compartment
covid_time <- 100000
mobility_date <- "2020-02-21"

#colnames other reported states
other_reported_states <- c("Other_reported1",  "Other_reported2", "Other_reported3",
                           "Other_reported4", "Other_reported5",  "Other_reported6",
                           "Other_reported7",  "Other_reported8", "Other_reported9",
                           "Other_reported10", "Other_reported11", "Other_reported12",
                           "Other_reported13", "Other_reported14", "Other_reported15",
                           "Other_reported16")

all_others <- c("OTHER_0", "OTHER_5",   "OTHER_10",  "OTHER_15",  "OTHER_20",
                "OTHER_25",  "OTHER_30", "OTHER_35",  "OTHER_40",  "OTHER_45", 
                "OTHER_50",  "OTHER_55",  "OTHER_60",  "OTHER_65", "OTHER_70",
                "OTHER_75")

#TODO these need fixing
all_covs <- c("COVID_p0", "COVID_p5", "COVID_p15", "COVID_p45", "COVID_p65")
all_oths <- c("OTHER_p0" , "OTHER_p5",  "OTHER_p15", "OTHER_p45", "OTHER_p65")


init_state_temp <- c(
9.754441e+05,     0.000000e+00,     0.000000e+00,     0.000000e+00,     1.077689e+03, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     2.350224e+03,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     2.229414e+06,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     8.913192e+02,     0.000000e+00,     1.409270e+03, 
0.000000e+00,     1.509306e+07,     3.486487e+05,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     6.729822e+02,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
1.457567e+03,     0.000000e+00,     0.000000e+00,     0.000000e+00,     2.887670e+06, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     0.000000e+00,     5.541048e+02, 
0.000000e+00,     8.701346e+02,     0.000000e+00,     1.180942e+07,     3.259776e+05, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     5.130331e+02,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     1.115641e+03,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     2.973897e+06,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     4.232164e+02,     0.000000e+00,     6.665336e+02,     0.000000e+00, 
9.720894e+06,     3.559257e+05,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
5.147016e+02,     0.000000e+00,     0.000000e+00,     0.000000e+00,     1.122411e+03, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     3.018398e+06,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     4.251990e+02,     0.000000e+00, 
6.711055e+02,     0.000000e+00,     9.122837e+06,     5.060726e+05,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     4.794893e+02,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     1.052552e+03,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
2.926998e+06,     0.000000e+00,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
3.970725e+02,     0.000000e+00,     6.290323e+02,     0.000000e+00,     8.837782e+06, 
5.134229e+05,     0.000000e+00,     0.000000e+00,     0.000000e+00,     5.284126e+02, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     1.158071e+03,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     2.968699e+06,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     4.373936e+02,     0.000000e+00,     6.924428e+02, 
0.000000e+00,     9.804763e+06,     5.225209e+05,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     5.370261e+02,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
1.176280e+03,     0.000000e+00,     0.000000e+00,     0.000000e+00,     3.027664e+06, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     0.000000e+00,     4.443876e+02, 
0.000000e+00,     7.031912e+02,     0.000000e+00,     1.002337e+07,     5.128891e+05, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     5.594808e+02,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     1.224302e+03,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     3.132963e+06,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     4.628144e+02,     0.000000e+00,     7.319785e+02,     0.000000e+00, 
1.021714e+07,     5.520182e+05,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
5.582162e+02,     0.000000e+00,     0.000000e+00,     0.000000e+00,     1.223188e+03, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     3.188034e+06,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     4.619854e+02,     0.000000e+00, 
7.311901e+02,     0.000000e+00,     9.627853e+06,     6.057192e+05,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     5.600158e+02,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     1.229463e+03,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
3.194441e+06,     0.000000e+00,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
4.638140e+02,     0.000000e+00,     7.348977e+02,     0.000000e+00,     9.777033e+06, 
6.743860e+05,     0.000000e+00,     0.000000e+00,     0.000000e+00,     5.540832e+02, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     1.218332e+03,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     3.155897e+06,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     4.591651e+02,     0.000000e+00,     7.281668e+02, 
0.000000e+00,     9.843428e+06,     7.417126e+05,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     5.509990e+02,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
1.212909e+03,     0.000000e+00,     0.000000e+00,     0.000000e+00,     3.112879e+06, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     0.000000e+00,     4.568065e+02, 
0.000000e+00,     7.249003e+02,     0.000000e+00,     9.486276e+06,     7.532212e+05, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     5.667658e+02,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     1.246823e+03,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     3.114512e+06,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     4.697803e+02,     0.000000e+00,     7.452541e+02,     0.000000e+00, 
8.896737e+06,     8.184531e+05,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
5.295006e+02,     0.000000e+00,     0.000000e+00,     0.000000e+00,     1.165918e+03, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     3.010569e+06,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     4.389972e+02,     0.000000e+00, 
6.966741e+02,     0.000000e+00,     7.883138e+06,     6.647327e+05,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     3.593526e+02,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     7.929791e+02,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
2.099166e+06,     0.000000e+00,     0.000000e+00,     0.000000e+00,     0.000000e+00, 
2.981594e+02,     0.000000e+00,     4.737156e+02,     0.000000e+00,     6.052353e+06, 
1.854464e+06,     0.000000e+00,     0.000000e+00,     0.000000e+00,     8.382373e+02, 
0.000000e+00,     0.000000e+00,     0.000000e+00,     1.850774e+03,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     4.643035e+06,     0.000000e+00,     0.000000e+00, 
0.000000e+00,     0.000000e+00,     6.956788e+02,     0.000000e+00,     1.105732e+03, 
0.000000e+00,     1.191610e+07)
