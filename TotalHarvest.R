options(java.parameters="-Xmx8024m")
library(XLConnect)
library(ggplot2)
library(data.table)

data.path <- "../Data"
dir.create(file.path(data.path, "Analysis"))

opp.shoot <- data.table(readWorksheetFromFile(
  file=file.path(data.path, "Tracking_History_Opportunistic.xlsx"), 
  sheet="OPP_Tracking_Data"))
nms.o <- names(opp.shoot)

judas.master <- data.table(readWorksheetFromFile(
  file=file.path(data.path, "Tracking_History_Judas_with habitat_MZ.xlsx"), 
  sheet="Judas_Tracking_History"))
nms.jm <- names(judas.master)

# Rm useless spaces in REGION
judas.master[, REGION := sub(pattern=" +", replacement="", x=REGION)]

# match headers
setnames(opp.shoot, "Region", "Region_ID")
setnames(judas.master, "REGION", "Region")
setnames(judas.master, nms.jm[c(3:4, 6:8, 12)], nms.o[c(2:3, 5:8)])

# Add Region to opp.shoot
opp.shoot[, Region := ifelse(Region_ID == "PB", "PILBARA", "KIMBERLEY")]
opp.shoot[, Method := "Opportunistic"]

# combine opportunistic and judas culling
setkey(judas.master, EVENT_ID)
judas.culling <- judas.master[unique(EVENT_ID), nms.jm[c(1:4, 6:8, 12)], 
                              with=FALSE, mult="first"]
judas.culling <- judas.culling[N_Ferals > 0,]
judas.culling[, Method := "Judas"]
judas.dead <- judas.master[ACTION == "DEAD" | ACTION == "SHOT", 
                           nms.jm[c(1:5, 6:8, 12, 14)], with=FALSE]

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

totals <- tot.harvest[, .(N=sum(N_Ferals)), by=c("Year", "Region", "Method")]
totals[, Total := sum(N), by=Year]


ggplot(totals, aes(Year, N, colour=Region, shape=Method)) + geom_point() + 
  geom_line(aes(Year, Total), colour="black")

ggsave(filename = file.path(data.path, "Analysis", "TotalHarvest.pdf"))
write.csv(tot.harvest, file = file.path(data.path, "tot.harvest.csv"))
write.csv(totals, file = file.path(data.path, "totals.csv"))
