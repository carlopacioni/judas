library(readxl)  # I think this is better than XLconnect and it doesn't require java !
library(ggplot2)
library(data.table)
library(sf, quietly=T)

data.path <- "../Data"
dir.create(file.path(data.path, "Analysis"))

opp.shoot <- data.table(read_excel(
  path=file.path(data.path, "Tracking_History_Opportunistic.xlsx"), 
  sheet="OPP_Tracking_Data"))
nms.o <- names(opp.shoot)

judas.master <- data.table(read_excel(
  path=file.path(data.path, "Tracking_History_Judas_with habitat_MZ.xlsx"), 
  sheet="Judas_Tracking_History"))
nms.jm <- names(judas.master)

# Pull in Kimberly shire LGAs
kimb_LGA<- st_read("../Data/Shapefile/.","Kimberley_LGA")
kimb_LGA<- kimb_LGA %>% st_transform(4326) # convert to Lat/Long

# Insert a col with location (shire) where the animals have been tracked, e.g. 'TrackedShire'
# Judas.master
judas.master<- judas.master %>% filter(!is.na(LONG))
judas.master_sf<- st_as_sf(judas.master, coords = c("LONG", "LAT"), crs = 4326, agr = "identity")
curr_LGA<- st_intersects(judas.master_sf, kimb_LGA, sparse=F)
judas.master<- judas.master %>% mutate(TrackedShire = case_when(
  curr_LGA[,1] ~ "HALL",
  curr_LGA[,2] ~ "PHED",
  curr_LGA[,3] ~ "EPIL",
  curr_LGA[,4] ~ "BROOME",
  curr_LGA[,5] ~ "DERB",
  curr_LGA[,6] ~ "WYND",
  TRUE ~ "other"))

# opp.shoot
#judas.master<- judas.master %>% filter(!is.na(LONG))
opp.shoot_sf <- st_as_sf(opp.shoot, coords = c("Longitude", "Latitude"), 
                           crs = 4326, agr = "identity")
curr_LGA <- st_intersects(opp.shoot_sf, kimb_LGA, sparse=F)
opp.shoot <- opp.shoot %>% mutate(TrackedShire = case_when(
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
                                    ifelse(TrackedShire == "WYND", "EW", "other"))))]

opp.shoot[, TrackedShire:=
               ifelse(TrackedShire == "PHED" | TrackedShire == "EPIL", "PB",
                      ifelse(TrackedShire == "BROOME" | TrackedShire == "DERB", "WK", 
                             ifelse(TrackedShire == "HALL", "HC", 
                                    ifelse(TrackedShire == "WYND", "EW", "other"))))]

# Rm useless spaces in REGION
judas.master[, REGION := sub(pattern=" +", replacement="", x=REGION)]

# Create year and month columns
judas.master[, Year:=year(EVENT_DATE)]
judas.master[, Month:=month(EVENT_DATE)]

#### Effort ####
# Count instances
runs <- judas.master[, .N, by=c("EVENT_DATE", "Year", "Month", "REGION", "TrackedShire")]
ncollared <- judas.master[ACTION == "COLLARED", .(Ncollared=.N),
                          by=c("EVENT_DATE", "Year", "Month", "REGION", "TrackedShire")]
# Meerge together
runs <- merge(runs, ncollared, all.x=T, 
              on=c("EVENT_DATE", "Year", "Month", "REGION", "TrackedShire"))

runs_month <- runs[, .(.N, Ncollared=sum(Ncollared, na.rm = T)), 
                   by=c("Year", "Month", "REGION", "TrackedShire")]

runs_month <- runs_month[, Effort := N * 7.5 - Ncollared * 0.5]
surveyed_months <- runs_month[, .N, by=c("REGION", "TrackedShire", "Month")]

# Number of year of the program per region
runs_month[, .(Min=min(Year), Max=max(Year), nYear=max(Year) - min(Year)), 
           by=REGION]

ggplot(surveyed_months, aes(factor(Month), N)) + geom_histogram(stat="identity") + 
  facet_grid(TrackedShire~REGION) + geom_hline(yintercept = 19, col="red") + 
  geom_hline(yintercept = 23, col="blue")

ggplot(data = runs_month, aes(x=factor(Month), y=Year)) +
  geom_raster(aes(fill = N), interpolate = F) +
  scale_fill_gradient(low = "yellow", high = "red") + facet_grid(TrackedShire~REGION)

write.csv(runs_month, file.path(data.path, "Effort.csv"), row.names = F)  
write.csv(judas.master, file.path(data.path, "judas.master.Effort.csv"), row.names = F)  

#### Total Harvest ####

# match headers
setnames(opp.shoot, "Region", "Region_ID")
setnames(judas.master, "REGION", "Region")
setnames(judas.master, nms.jm[c(3:4, 6:8, 12)], nms.o[c(2:3, 5:8)])

# Add Region to opp.shoot
opp.shoot[, Region := ifelse(Region_ID == "PB", "PILBARA", "KIMBERLEY")]
opp.shoot[, Method := "Opportunistic"]

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

ggsave(filename = file.path(data.path, "Analysis", "TotalHarvest.pdf"))
write.csv(tot.harvest, file = file.path(data.path, "tot.harvest.csv"))
write.csv(totals, file = file.path(data.path, "totals.csv"))
