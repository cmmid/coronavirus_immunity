#Covid sims load in extra data

####### COVID DATA #######

deaths_covid <- data.table(read.csv(here("simulating_covid","ONS_daily_deaths.csv")))
deaths_covid[, All := Other + Home + Care.home + Hospital..acute.or.community..not.psychiatric.]
deaths_covid <- deaths_covid[, c("Date","All") ]
deaths_covid[, date := paste0(Date, "-2020")]
deaths_covid[, date := as.Date(date, format = "%d-%b-%Y")]
deaths_covid <- deaths_covid[ date < as.Date("2020-06-01"), ]
colnames(deaths_covid) <- c("Date_char", "actual_deaths", "date")

# hospital and other care home deaths only (excluding care homes) - if want to
deaths_covid_ch <- data.table(read.csv(here("simulating_covid","ONS_location_deaths.csv")))
deaths_covid_ch[, deaths_wanted := Hospital + Elsewhere + Care.home]
deaths_covid_ch[, c("Care.home", "Hospital", "Elsewhere") := NULL]
colnames(deaths_covid_ch) <- c("date", "resident_deaths")
deaths_covid_ch[, date := as.Date(date, format = "%d/%m/%Y")]
deaths_covid[deaths_covid_ch, on = c("date"), resident_deaths := i.resident_deaths]

# deaths_covid[, non_resident := actual_deaths - resident_deaths]
deaths_covid[, non_resident := actual_deaths ]
deaths_covid <- na.omit(deaths_covid)
# deaths_covid[, week := week(as.Date(date))]
# deaths_covid[nchar(week)==1, week_end :=  as.Date(paste0("2020-0",week,"-1"), 
#                                                   format = "%Y-%W-%w")]
# deaths_covid[nchar(week)==2, week_end :=  as.Date(paste0("2020-",week,"-1"), 
#                                                   format = "%Y-%W-%w")]
# deaths_covid <- deaths_covid[,sum(actual_deaths), by=week_end]

deaths_covid <- deaths_covid[,c("date", "non_resident")]

colnames(deaths_covid) <- c("date", "actual_deaths")
# remove the last point as it's not a full week
deaths_covid[date == "2020-06-15", actual_deaths := NA]
deaths_covid <- na.omit(deaths_covid)

# remove points after intervention time
deaths_covid <- deaths_covid[date <"2020-06-01",]

####### MOBILITY DATA ######
google_mobility <- data.table(read.csv(here("simulating_covid","uk_google_mobility.csv")))
google_mobility$date <- as.Date(google_mobility$date, format = "%d/%m/%Y")

google_mobility[,average := ( retail_and_recreation_percent_change_from_baseline+
                                workplaces_percent_change_from_baseline +
                                grocery_and_pharmacy_percent_change_from_baseline +
                                transit_stations_percent_change_from_baseline)  /4]
google_mobility$rolling_mean <- frollmean(google_mobility$average, n=7)

google_mobility <- na.omit(google_mobility)
google_mobility_use <- google_mobility[,c("date", "rolling_mean")]

# Infection fatality rates - from Levin paper
IFRs <- data.table(
  lower = c(0,35,45,55,65,75,85), 
  upper = c(34,44,54,64,74,84,Inf), 
  proportion = c(0.004,0.064,0.21,0.7,2.3,7.6,22.3)
)
#IFRs population weigthed in correct age group
IFRs_actual <- data.table(
  lower = c(0,5,15,45,65), 
  upper = c(4,14,44,64,Inf), 
  percent = as.numeric(c(IFRs[lower==0,"proportion"],
                         IFRs[lower==0,"proportion"],
                         ((IFRs[lower==0,"proportion"]*sum(pop_numbers[4:7]) + 
                             IFRs[lower==35,"proportion"]*sum(pop_numbers[8:9]))/
                            sum(pop_numbers[4:9])), 
                         ((IFRs[lower==45,"proportion"]*sum(pop_numbers[10:11]) + 
                             IFRs[lower==55,"proportion"]*sum(pop_numbers[12:13]))/
                            sum(pop_numbers[10:13])),
                         ((IFRs[lower==65,"proportion"]*sum(pop_numbers[14:15]) + 
                             IFRs[lower==75,"proportion"]*sum(pop_numbers_5year[16:17]) + 
                             IFRs[lower==85,"proportion"]*sum(pop_numbers_5year[18:19]))/
                            sum(pop_numbers_5year[14:19])))
  )
)

IFRs_actual[, proportion := percent/100 ]
IFRs_actual[, log_odds := log(proportion/(1-proportion))]
