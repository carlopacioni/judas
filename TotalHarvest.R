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

write.csv(judas.master, file.path(data.path, "Judas.master.TrackedShire.csv"))

#### Effort ####

# Effort in secs over a hr (3600 secs) is per individual. Longer time is needed 
# for larger groups (except camels) because animals separate
judas.master[, EffortShootingFeralsWithJudas:=N_Ferals * 
       ifelse(N_Ferals < 10, 10/3600, 60/3600)]
opp.shoot[, EffortShootingHorses:=N_Horses * 60/3600]
opp.shoot[, EffortShootingCamel:=N_Camels * ifelse(N_Camels<15, 10/3600, 20/3600)]
opp.shoot[, EffortShootingFerals:=N_Ferals * ifelse(N_Ferals < 10, 10/3600, 60/3600)]
opp.shoot.days <- opp.shoot[, .(N_Ferals=sum(N_Ferals), N_Horses=sum(N_Horses), 
                                N_Camels=sum(N_Camels), 
                                EffortShootingHorses=sum(EffortShootingHorses),
                                EffortShootingCamel=sum(EffortShootingCamel),
                                EffortShootingFerals=sum(EffortShootingFerals)),
                            by=c("Date", "Year", "Month", "Region")]

# Count instances
#feralsWithJudas <- judas.master[SEARCHED == "YES", .(N_Ferals, EffortShootingFeralsWithJudas), 
#                     by=c("EVENT_ID", "Date", "Year", "Month", "Region")]
runs <- judas.master[SEARCHED == "YES", .(N_FeralsWithJudas=sum(N_Ferals),
                         EffortShootingFeralsWithJudas=sum(EffortShootingFeralsWithJudas)), 
                     by=c("Date", "Year", "Month", "Region")]
ncollared <- judas.master[ACTION == "COLLARED", .(Ncollared=.N),
                          by=c("Date", "Year", "Month", "Region")]
nshot <- judas.master[ACTION == "SHOT", .(N_JudasShot=.N),
                          by=c("Date", "Year", "Month", "Region")]
njudasTracked <- judas.master[SEARCHED == "YES", .(N_JudasTracked=.N), 
                              by=c("Date", "Year", "Month", "Region")]
judas.dead <- judas.master[ACTION == "DEAD", .(N_JudasDead=.N), 
                           by=c("Date", "Year", "Month", "Region")]

# Merge together
runs <- merge(runs, ncollared, all.x=T, on=c("Date", "Year", "Month", "Region"))
runs <- merge(runs, nshot, all.x=T, on=c("Date", "Year", "Month", "Region"))
runs <- merge(runs, njudasTracked, all.x=T, on=c("Date", "Year", "Month", "Region"))
runs <- merge(runs, judas.dead, all.x=T, on=c("Date", "Year", "Month", "Region"))

# Note all=T here. That is because there are some dates that are not in the judas data
runs <- merge(runs, opp.shoot.days, all=TRUE, on=c("Date", "Year", "Month", "Region"))

runs <- runs[Year < 2018,]

# Replace NA 
replace.na <- function(DT, col_names) {
  for (j in col_names)
    set(DT,which(is.na(DT[[j]])),j,0)
}

col_names <- c("N_FeralsWithJudas", "EffortShootingFeralsWithJudas", "Ncollared",
               "N_JudasShot", "N_Ferals", "N_Horses", "N_Camels", "N_JudasDead",
               "EffortShootingHorses", "EffortShootingCamel", "EffortShootingFerals")

replace.na(runs, col_names)

# Update judasTracked removing judas collared
runs[!is.na(N_JudasTracked), N_JudasTracked:=N_JudasTracked - Ncollared]

runs[Ncollared==0, summary(N_JudasTracked), by=Region]
runs[Ncollared==0, quantile(N_JudasTracked, probs = 0.90, na.rm = T), by=Region]
runs[, summary(N_JudasTracked), by=Region]
runs[N_JudasTracked>100,]
runs[is.na(N_JudasTracked),]
runs[N_JudasTracked>30, .N]

ggplot(runs[Ncollared==0,], aes(N_JudasTracked)) + geom_histogram(bins = 185) +facet_grid(Region~.)
ggsave(filename = file.path(data.path, "Analysis", "N_JudasTracked.pdf"))

runs[, N_JudasTracked := ifelse(is.na(N_JudasTracked), 15, N_JudasTracked)]

runs[, AdjNdayTrips:=ifelse(N_JudasTracked/30>1, ceiling(N_JudasTracked/30), 1)]

runs_years <- runs[, .(NdayTrips=.N, AdjNdayTrips=sum(AdjNdayTrips), 
                       N_JudasTracked=sum(N_JudasTracked, na.rm=T),
                       Ncollared=sum(Ncollared, na.rm=T),
                       N_JudasShot=sum(N_JudasShot, na.rm=T),
                       N_JudasDead=sum(N_JudasDead),
                       N_FeralsWithJudas=sum(N_FeralsWithJudas, na.rm = T),
                       N_Ferals=sum(N_Ferals, na.rm=T),
                       N_Horses=sum(N_Horses, na.rm = T),
                       N_Camels=sum(N_Camels, na.rm = T),
                       EffortShootingFeralsWithJudas=sum(EffortShootingFeralsWithJudas),
                       EffortShootingFerals=sum(EffortShootingFerals),
                       EffortShootingHorses=sum(EffortShootingHorses),
                       EffortShootingCamel=sum(EffortShootingCamel)), 
                   by=c("Year", "Region")]

# Totals harveted 
runs_years[, Judas.culling:=N_FeralsWithJudas + N_JudasShot + N_JudasDead] 
runs_years[, Total_Harvest:=sum(Judas.culling) + sum(N_Ferals), by=Year]

runs_years[, Effort_judas:=
             N_JudasTracked * 6.5/30 + # Baseline effort: time allocated to tracked a judas when no shoot 
             # occurs (obtained by dividing the flying time (6.5 hrs) by the 
             # max number of judas taht can be tracked in a day
             EffortShootingFeralsWithJudas +
             N_JudasShot * 10/3600]

runs_years <- runs_years[, Effort_opp:=AdjNdayTrips * 6.5 - # Flying time
                           Ncollared * 0.5 - # half hr to collar
                           N_JudasShot * 300/3600 - # 5' to recover the collars
                           EffortShootingFeralsWithJudas -
                           EffortShootingHorses -
                           EffortShootingCamel]

ggplot(runs_years) + geom_histogram(aes(EffortShootingHorses), alpha=0.5, fill="blue") +
  geom_histogram(aes(EffortShootingCamel), fill="yellow", alpha=0.5)

#surveyed_months <- runs_years[, .N, by=c("Region", "TrackedShire", "Month")]

runs_years[, summary(Effort_opp)]
runs_years[, summary(Effort_judas)]

# Long format
runs_years_long <- melt(runs_years, id.vars = c("Year", "Region"),
                        measure.vars = c("Judas.culling", "N_Ferals"),
                        variable.name = "Method", value.name = "N_Donkeys")

Effort_years_long <- melt(runs_years, id.vars = c("Year", "Region"),
                        measure.vars = c("Effort_judas", "Effort_opp"),
                        variable.name = "Method", value.name = "Effort")

runs_years_long[, Effort:=Effort_years_long[, Effort]]
runs_years_long <- merge(runs_years_long, 
                         runs_years[, .(Year, Region, N_Horses, N_Camels, Total_Harvest)], 
                         all.x = T,
                         by=c("Year", "Region"))

runs_years_long[, Method:=ifelse(Method == "N_Ferals", "Opportunistic", "Judas")]

ggplot(runs_years_long, aes(Year, N_Donkeys, colour=Region, shape=Method)) + geom_point() + 
  geom_line(aes(Year, Total_Harvest), colour="black")
ggsave(filename=file.path(data.path, "Analysis", "TotalHarvest.pdf"))

write.csv(runs_years, file=file.path(data.path, "tot.harvest.wide.csv"), row.names=F)
write.csv(runs_years_long, file=file.path(data.path, "tot.harvest.long.csv"), row.names=F)

ggplot(runs_years_long, aes(Year, Effort, fill=Method)) + 
  geom_histogram(stat="identity") +
  facet_grid(Region~.) + ylab("Effort")
ggsave(filename = file.path(data.path, "Analysis", "Effort_Year.pdf"))

# Number of year of the program per region
runs_years[, .(Min=min(Year), Max=max(Year), nYear=max(Year) - min(Year)), 
           by=Region]


p_Effort_harvest <- ggplot(runs_years_long, aes(Year, Effort)) + 
  geom_point(aes(size=N_Donkeys, col=N_Donkeys)) +
  scale_colour_gradient(low="yellow", high="red") + facet_grid(Region~Method)
p_Effort_harvest

ggsave(filename=file.path(data.path, "Plot_Harvest_effort.pdf"), plot=p_Effort_harvest)

CumEffort_harvest <- runs_years_long[, .(CumEff=cumsum(Effort), CumHarvest=cumsum(N_Donkeys)), 
                by=c("Region", "Method")]

p_CumEffort_harvest <- ggplot(CumEffort_harvest, aes(CumEff, CumHarvest)) + 
  geom_point(shape=1) + facet_grid(Region~Method) 
p_CumEffort_harvest

