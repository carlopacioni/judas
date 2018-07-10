library(ggplot2, quietly = T)
library(data.table, quietly = T)
library(geosphere)
library(RcppAlgos)

calc.latlong.dist<- function(xy1,xy2)
{
  # uses spherical law of cosines to calculate distance between two lat/long
  # coordinates in decimal degrees
  R <- 6371 # Earths radius
  xy1 <- (pi * xy1)/180 # radians
  xy2 <- (pi * xy2)/180 
  D <- acos(sin(xy1[,1])*sin(xy2[,1]) + cos(xy1[,1])*cos(xy2[,1])*cos(xy2[,2]-xy1[,2]))  
  return(R*D)
}
#--------------------------------------------------

data.path <- "../Data/"
load(file = file.path(data.path, "judas.cleaned.rda"))

# Home Range centres
judas.cleaned[, ':='(HRlat=mean(LAT), HRlong=mean(LONG)), by=JUDAS_ID]

names(judas.cleaned)
setkey(judas.cleaned, JUDAS_ID)

# List all info re each judas
list.IDs <- judas.cleaned[, unique(JUDAS_ID)]
judas.info <- judas.cleaned[list.IDs, .SD, 
                            .SDcols=c("JUDAS_ID","REGION", "start.date", "end.date", 
                                      "Time.deployment", "Time.dep.years", 
                                      "HRlong", "HRlat"),
                            mult="first"]
setkey(judas.info, JUDAS_ID)

# Find all possible associations
# RcppAlgos::comboGeneral is much faster than combn
combs <- RcppAlgos::comboGeneral(judas.cleaned[, unique(JUDAS_ID)], 2)

# Collate info together
judas.Dist<- data.table(combs)
setnames(judas.Dist, c("JUDAS_ID", "id2"))
setkey(judas.Dist, JUDAS_ID)

judas.Dist <- merge(judas.Dist, judas.info, all.x = T)
nms <- names(judas.Dist)[c(-1, -2)]
setnames(judas.Dist, c("JUDAS_ID", "id2", nms), 
         c("ID.1", "JUDAS_ID", sapply(nms, paste0, ".1")))

length.nms <- length(names(judas.Dist))
setkey(judas.Dist, JUDAS_ID)
judas.Dist <- merge(judas.Dist, judas.info, all.x = T)
ncols <- length(judas.Dist)
setnames(judas.Dist, names(judas.Dist)[(length.nms + 1):ncols], 
         sapply(names(judas.Dist)[(length.nms + 1):ncols], paste0, ".2"))
setnames(judas.Dist, "JUDAS_ID", "ID.2")

# Calculate distance between HR centres
judas.Dist[, Dist := calc.latlong.dist(judas.Dist[, .(HRlat.1, HRlong.1)], 
                                    judas.Dist[, .(HRlat.2, HRlong.2)])]
key(judas.cleaned) # Confirm that Judas_ID is the key

judas.Dist[, Concurrent := 
             ifelse(end.date.1 < start.date.2 | end.date.2 < start.date.1, 'no', "yes")]
setkey(judas.Dist, Concurrent)
judas.Dist <- judas.Dist["yes", ] # get rid of judas that were not present at same time

# Prepare columns
# event is 0 when there is no encounter, and 1 when there is
# time is the length of time when both were 'deployed' (the beginning of the 
   # period when concurrently present) and when they were observed together if event=1
   # otherwise the max time when both were 'deployed'  
# maxtime is the max time when both were 'deployed'  
judas.Dist[, ':='(SHIRE=as.character(NA), AREA=as.character(NA), EVENT_DATE=as.Date(NA), 
                  event=as.numeric(NA), time=as.numeric(NA), maxtime=as.numeric(NA), 
                  EVENT_ID=as.numeric(NA), EVENT_CODE=as.character(NA), 
                  ACTION=as.character(NA), 
                  LAT=as.numeric(NA), LONG=as.numeric(NA), Habitat.Type=as.character(NA))]
jd.dates <- c("start.date.1", "end.date.1", "start.date.2", "end.date.2")
judas.Dist[, (jd.dates) := lapply(.SD, as.Date), .SDcols=jd.dates]

judas.cleaned[, EVENT_DATE := as.Date(EVENT_DATE)] 


rns <- nrow(judas.Dist) # Number of judas pairs
# There must be a clever way to do this, but couldn't think any...
system.time(
for(rn in seq_len(rns)) {
  data.ID1 <- judas.cleaned[judas.Dist[rn, ID.1], ]
  data.ID2 <- judas.cleaned[judas.Dist[rn, ID.2], ]
  matching.Events <- data.ID1[, EVENT_ID] %in% data.ID2[, EVENT_ID]
  data.Evs <- data.ID1[matching.Events, ]
  if(nrow(data.Evs) > 0) {
    line.no <- which.min(data.Evs[, EVENT_DATE])
    judas.Dist[rn, ':='(SHIRE=data.Evs[line.no, SHIRE], AREA=data.Evs[line.no, AREA], 
                        EVENT_DATE=data.Evs[line.no, EVENT_DATE], event=1, 
                        EVENT_ID=data.Evs[line.no, EVENT_ID], 
                        EVENT_CODE=data.Evs[line.no, EVENT], 
                        ACTION=data.Evs[line.no, ACTION], 
                        LAT=data.Evs[line.no, LAT], LONG=data.Evs[line.no, LONG], 
                        Habitat.Type=data.Evs[line.no, Habitat.Type])]
  } else {
    judas.Dist[rn, ':='(EVENT_DATE=min(end.date.1, end.date.2), event=0)]
  }
  judas.Dist[rn, ':='(time=EVENT_DATE - max(start.date.1, start.date.2), 
                      maxtime=min(end.date.1, end.date.2) - max(start.date.1, start.date.2))]
}
)

write.csv(judas.Dist, file = file.path(data.path, "judas.Dist.csv"), row.names = F)
save(judas.Dist, file = file.path(data.path, "judas.Dist.RDA"))

