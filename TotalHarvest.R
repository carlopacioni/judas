library(readxl)  # I think this is better than XLconnect and it doesn't require java !
library(ggplot2)
library(data.table)
library(sf, quietly=T)
library(dplyr)

data.path <- "../Data"
dir.create(file.path(data.path, "Analysis"))

opp.shoot <- data.table(read_excel(
  path=file.path(data.path, "Tracking_History_Opportunistic_2.xlsx"), 
  sheet="OPP_Tracking_Data_All Species"))
nms.o <- c(names(opp.shoot))

# Add Region col
opp.shoot[, Region_ID:=ifelse(Shire == "PB", "PB", "KM")]
setcolorder(opp.shoot, c("Region_ID", nms.o))
setnames(opp.shoot, "N_Donkeys", "N_Ferals")
nms.o <- c(names(opp.shoot))

judas.master <- data.table(read_excel(
  path=file.path(data.path, "Tracking_History_Judas_with Habitat_Soil_and_CC_2.xlsx"), 
  sheet="Judas_Tracking_History"))
nms.jm <- c(names(judas.master))

# Pull in Kimberly shire LGAs
kimb_LGA<- st_read("../Data/Shapefile/.","Kimberley_LGA")
kimb_LGA<- kimb_LGA %>% st_transform(4326) # convert to Lat/Long

# Insert a col with location (shire) where the animals have been tracked, e.g. 'TrackedShire'
# Judas.master
judas.master<- judas.master %>% filter(!is.na(LONG))
judas.master_sf<- st_as_sf(judas.master, coords=c("LONG", "LAT"), crs=4326, agr="identity")
curr_LGA<- st_intersects(judas.master_sf, kimb_LGA, sparse=F)
judas.master<- judas.master %>% mutate(TrackedShire=case_when(
  curr_LGA[,1] ~ "HALL",
  curr_LGA[,2] ~ "PHED",
  curr_LGA[,3] ~ "EPIL",
  curr_LGA[,4] ~ "BROOME",
  curr_LGA[,5] ~ "DERB",
  curr_LGA[,6] ~ "WYND",
  TRUE ~ "other"))

# opp.shoot
opp.shoot_sf <- st_as_sf(opp.shoot, coords=c("Longitude", "Latitude"), 
                           crs=4326, agr="identity")
curr_LGA <- st_intersects(opp.shoot_sf, kimb_LGA, sparse=F)
opp.shoot <- opp.shoot %>% mutate(TrackedShire=case_when(
  curr_LGA[,1] ~ "HALL",
  curr_LGA[,2] ~ "PHED",
  curr_LGA[,3] ~ "EPIL",
  curr_LGA[,4] ~ "BROOME",
  curr_LGA[,5] ~ "DERB",
  curr_LGA[,6] ~ "WYND",
  TRUE ~ "other"))

# Convert back to data.table
judas.master <- as.data.table(judas.master)
opp.shoot <- as.data.table(opp.shoot)

# Group actual shire in four SHIRE
judas.master[, TrackedShire:=
               ifelse(TrackedShire == "PHED" | TrackedShire == "EPIL", "PB",
                      ifelse(TrackedShire == "BROOME" | TrackedShire == "DERB", "WK", 
                             ifelse(TrackedShire == "HALL", "HC", 
                                    ifelse(TrackedShire == "WYND", "EK", "other"))))]
judas.master[LAT > -17 & LONG > 128.5 & TrackedShire == "other", TrackedShire:="EK"]
judas.master[LAT < -17 & LONG > 128.5& TrackedShire == "other", TrackedShire:="HC"]
judas.master[LAT < -19.98 & LONG < 120.25 & TrackedShire == "other", TrackedShire:="PB"]

opp.shoot[, TrackedShire:=
               ifelse(TrackedShire == "PHED" | TrackedShire == "EPIL", "PB",
                      ifelse(TrackedShire == "BROOME" | TrackedShire == "DERB", "WK", 
                             ifelse(TrackedShire == "HALL", "HC", 
                                    ifelse(TrackedShire == "WYND", "EK", "other"))))]
opp.shoot[Latitude > -17 & Longitude > 128.5 & TrackedShire == "other", TrackedShire:="EK"]
opp.shoot[Latitude < -17 & Longitude > 128.5 & TrackedShire == "other", TrackedShire:="HC"]
opp.shoot[Latitude < -19.98 & Longitude < 120.25 & TrackedShire == "other", TrackedShire:="PB"]

opp.shoot[TrackedShire == "other", .N]
judas.master[TrackedShire == "other", .N]
judas.master[TrackedShire == "other", ]

# match headers
setnames(judas.master, "REGION", "Region")
setnames(judas.master, nms.jm[c(3:4, 6:8, 12)], nms.o[c(2:3, 5:8)])

# Add Region and Method to opp.shoot
opp.shoot[, Region := ifelse(Region_ID == "PB", "PILBARA", "KIMBERLEY")]

# Create year and month columns
judas.master[, Year:=year(Date)]
judas.master[, Month:=month(Date)]

opp.shoot[, Year:=year(Date)]
opp.shoot[, Month:=month(Date)]

#### Effort ####


# Count instances
runs <- judas.master[, .(N_JudasTracked=.N, N_Ferals=sum(N_Ferals)), 
                     by=c("Date", "Year", "Month", "Region", "SurvShire")]
# effort in secs over a hr (3600 secs) is per individual. Longer time is needed 
        # for larger groups (except camels) because animals separate
runs[, EffortShootingFeralsWithJudas:=N_Ferals * ifelse(N_Ferals<10, 10/3600, 60/3600)]
ncollared <- judas.master[ACTION == "COLLARED", .(Ncollared=.N),
                          by=c("Date", "Year", "Month", "Region", "SurvShire")]
nshot <- judas.master[ACTION == "SHOT", .(N_JudasShot=.N),
                          by=c("Date", "Year", "Month", "Region", "SurvShire")]

opp.shoot[, EffortShootingHorses:=N_Horses * 60/3600]
opp.shoot[, EffortShootingCamel:=N_Camels * ifelse(N_Camels<15, 10/3600, 20/3600)]
opp.shoot_month <- opp.shoot[, .(EffortShootingHorses=sum(EffortShootingHorses),
                                 EffortShootingCamel=sum(EffortShootingCamel)), 
                             by=c("Year", "Month", "Region", "SurvShire")]
opp.shoot_month[, lapply(.SD, summary), .SDcols=c("EffortShootingHorses", "EffortShootingCamel")]

ggplot(opp.shoot_month) + geom_histogram(aes(EffortShootingHorses), col="blue") +
  geom_histogram(aes(EffortShootingCamel), col="red")

# Merge together
runs <- merge(runs, ncollared, all.x=T, 
              on=c("Date", "Year", "Month", "Region", "SurvShire"))
runs <- merge(runs, nshot, all.x=T, 
              on=c("Date", "Year", "Month", "Region", "SurvShire"))

runs_month <- runs[, .(NdayTrips=.N, N_JUdasTracked=sum(N_JUdasTracked, na.rm=T),
                       Ncollared=sum(Ncollared, na.rm=T),
                       N_JudasShot=sum(N_JudasShot, na.rm=T),
                       N_Ferals=sum(N_Ferals, na.rm=T)), 
                   by=c("Year", "Month", "Region", "SurvShire")]

runs_month <- merge(runs_month, opp.shoot_month, all=TRUE, 
                    on=c("Year", "Month", "Region", "SurvShire"))
####
runs_month <- runs_month[TrackedShire != "other",]
runs_month <- runs_month[Year < 2018,]

runs_month[, summary(N_JUdasTracked)]

runs_month <- runs_month[, 
                         Effort_opp:=N * 6.5 - Ncollared * 0.5 - N_shot * 0.08 - Effort_judas]
surveyed_months <- runs_month[, .N, by=c("Region", "TrackedShire", "Month")]

runs_month[, summary(Effort_opp)]
runs_month[, summary(Effort_judas)]
runs_month[Effort_opp<0,]
runs_year <- runs_month[, .(Effort_opp=sum(Effort_opp), 
                            Effort_judas=sum(Effort_judas)), 
                      by=c("Year", "Region", "TrackedShire")]

# Number of year of the program per region
runs_month[, .(Min=min(Year), Max=max(Year), nYear=max(Year) - min(Year)), 
           by=Region]

ggplot(surveyed_months, aes(factor(Month), N)) + geom_histogram(stat="identity") + 
  facet_grid(TrackedShire~Region) + geom_hline(yintercept=19, col="red") + 
  geom_hline(yintercept=23, col="blue")

ggplot(data=runs_month, aes(x=factor(Month), y=Year)) +
  geom_raster(aes(fill=N), interpolate=F) +
  scale_fill_gradient(low="yellow", high="red") + facet_grid(TrackedShire~Region)

write.csv(runs_month, file.path(data.path, "Effort.csv"), row.names=F)  
write.csv(judas.master, file.path(data.path, "judas.master.Effort.csv"), row.names=F)  

#### Total Harvest ####

# Add col Method
opp.shoot[, Method := "Opportunistic"]

# Cross check dates
opp.dates <- opp.shoot[, .(NC=sum(N_Camels)), by=c("Date", "TrackedShire")]
judas.dates <- judas.master[, .(NF=sum(N_Ferals)), by=c("Date", "TrackedShire")]
test_dates <- merge(opp.dates[, .(NC=sum(NC)), by=Date], 
                         judas.dates[, .(NF=sum(NF)), by=Date], 
                         all.x=T, on="Date")
test_dates_Shire <- merge(opp.dates, judas.dates, all.x=T, on=c("TrackedShire", "Date"))
test_dates[is.na(NF), ]
test_dates_Shire[is.na(NF), ]
nrow(test_dates_Shire[is.na(NF), ])

sum(opp.shoot[, Date] %in% judas.master[, Date])
sum(judas.master[, Date] %in% opp.shoot[, Date])

# combine opportunistic and judas culling
setkey(judas.master, EVENT_ID)
judas.culling <- judas.master[unique(EVENT_ID), 
                              names(judas.master)[c(1:4, 6:8, 12, 18)], 
                              with=FALSE, mult="first"]
judas.culling <- judas.culling[N_Ferals > 0,]
judas.culling[, Method := "Judas"]
judas.dead <- judas.master[ACTION == "DEAD" | ACTION == "SHOT", 
                           names(judas.master)[c(1:5, 6:8, 12, 14, 18)], with=FALSE]

# cross check that judas don't appear twice
nrow(judas.dead) == length(
  jd_ID<-(judas.master[ACTION == "DEAD" | ACTION == "SHOT", unique(JUDAS_ID)]))

# find out which ones appear twice. I'm removing the dead entries not to have duplicates. 
# I have to confirm with Magda what happened here
jd <- judas.dead[, .N, by=JUDAS_ID]
judas.dead <- judas.dead[!(JUDAS_ID == jd[N>1, JUDAS_ID] & ACTION == "DEAD"),]

# Set N_Ferals=1 and rm cols that are not needed
judas.dead[, N_Ferals := 1]
judas.dead[, Method := "Judas"]
judas.dead[, ':='(JUDAS_ID=NULL, ACTION=NULL)]

tot.harvest <- rbindlist(list(judas.culling, opp.shoot, judas.dead), 
                         use.names=TRUE, fill=TRUE)

# Add months and years
tot.harvest[, ':='(Month=month(Date), Year=year(Date))]

totals <- tot.harvest[, .(N=sum(N_Ferals)), 
                      by=c("Year", "Region", "TrackedShire", "Method")]
totals_Region <- totals[, .(N=sum(N)), by=c("Year", "Region", "Method")]
totals_Region[, Total := sum(N), by=Year]


ggplot(totals_Region, aes(Year, N, colour=Region, shape=Method)) + geom_point() + 
  geom_line(aes(Year, Total), colour="black")

ggsave(filename=file.path(data.path, "Analysis", "TotalHarvest.pdf"))
write.csv(tot.harvest, file=file.path(data.path, "tot.harvest.csv"), row.names=F)
write.csv(totals, file=file.path(data.path, "totals.csv"), row.names=F)

#### Combine in the same table ###
#setnames(runs_year, "REGION", "Region")
TotalHarvest_Effort <- merge(totals, runs_year, all.x=T, 
                             on=c("Year", "Region", "TrackedShire"))
write.csv(TotalHarvest_Effort, file=file.path(data.path, "totals_with_effort.csv"), 
          row.names=F)

TotalHarvest_Effort[, summary(Effort_opp)]
TotalHarvest_Effort[, summary(log(Effort_opp))]
TotalHarvest_Effort[, summary(log(Effort_judas+1))]


p_EffOpp <- ggplot(TotalHarvest_Effort[Method == "Opportunistic" & TrackedShire != "other",], 
       aes(Year, log(Effort_opp))) + geom_point(aes(size=N, col=N)) +
  scale_colour_gradient(low="yellow", high="red") + facet_grid(TrackedShire~.)

p_EffJudas <- ggplot(TotalHarvest_Effort[Method == "Judas" & TrackedShire != "other",], 
                   aes(Year, log(Effort_judas + 1))) + geom_point(aes(size=N, col=N)) +
  scale_colour_gradient(low="yellow", high="red") + facet_grid(TrackedShire~.)

ggsave(filename = file.path(data.path, "Plot_Harvest_effort_Opportunistic.pdf"), plot = p_EffOpp)
ggsave(filename = file.path(data.path, "Plot_Harvest_effort_judas.pdf"), plot = p_EffJudas)










































