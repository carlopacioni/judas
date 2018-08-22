library(ggplot2, quietly = T)
library(data.table, quietly = T)
library(RcppAlgos)

#### Helper functions ####
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
#------------------------------------------------------------------------------#

# Start from judas.master dataset generated with TotalHarvest.R
data.path <- "../Data/"
judas.master <- fread(file.path(data.path, "Judas.master.TrackedShire.csv"))

# Rm found dead because they are out of the program
judas.cleaned <- judas.master[ACTION != "DEAD", ]

# Rm not found because they do not contribute 
judas.cleaned <- judas.cleaned[ACTION != "NONE", ]

# Set start and end date
judas.cleaned[, start.date := min(Date), by=JUDAS_ID]
judas.cleaned[, end.date := max(Date), by=JUDAS_ID]

# Cross check start.date matches collared
judas.cleaned[ACTION == "COLLARED", date.coll := Date, by=JUDAS_ID]
judas.cleaned[ACTION == "COLLARED", sum(start.date != date.coll, na.rm = T)]
judas.cleaned[start.date != date.coll, ]
# Some animals have been collared after their start date. This is a collar replacement
# Rm date.coll
judas.cleaned[, date.coll := NULL]

# Length in the program
judas.cleaned[, Time.deployment := difftime(end.date, start.date, units="weeks")]
judas.cleaned[, Time.dep.years := round(as.numeric(Time.deployment) / 52, 2)]

# Keep only animals that were tracked for less than 13 yrs. Longer time may 
  # reflect collars moved to other animals while not changing animal IDs
judas.cleaned <- judas.cleaned[Time.dep.years < 13, ]

# Rm animals that were collared on start.date and searched and not found (Year==0)
judas.cleaned <- judas.cleaned[Time.deployment>0, ]

# Rm animals that have coordinate outside Australia
judas.cleaned[TrackedShire == "other", .N]
judas.cleaned <- judas.cleaned[TrackedShire != "other",]

# subset judas.clean to 1995-2007 for KM and 1999-2001&2007-2010 for PB (PB runthe analysis with 2 month timewindow)
judas.cleaned[, Date := as.Date(Date)] 
judas.cleaned.surv <- judas.cleaned[(Region == "KIMBERLEY" & 
                                      Year >= 1995 & Year <= 2007) |
                                     (Region == "PILBARA" & 
                                        (Year >= 2007 & Year <= 2010)), ]

# Rm judas whose only entry is when they were collared
ntrack.events.surv <- judas.cleaned.surv[, .N, by=JUDAS_ID]
ntrack.events.surv[, summary(N)]
judas.cleaned.surv.one <- judas.cleaned.surv[JUDAS_ID %in% ntrack.events.surv[N==1, JUDAS_ID],]

judas.cleaned.surv <- judas.cleaned.surv[!JUDAS_ID %in% 
                                           judas.cleaned.surv.one[ACTION == "COLLARED", JUDAS_ID], ]

# Retain only IDs that are going to be used for surv analysis
# judas.cleaned.HR has more location s per animals because it retains also the 
    # location outside the year interval selected for survival analysis
judas.cleaned.HR <- judas.cleaned[JUDAS_ID %in% judas.cleaned.surv[, unique(JUDAS_ID)],]

# Home Range centres
judas.cleaned.HR[, ':='(HRlat=mean(Latitude), HRlong=mean(Longitude)), by=JUDAS_ID]

# Calculate deviations from HRcentres 
judas.cleaned.HR[, xdev:=calc.latlong.dist(judas.cleaned.HR[, .(Latitude, HRlong)],
                                        judas.cleaned.HR[, .(HRlat, HRlong)])]
judas.cleaned.HR[, ydev:=calc.latlong.dist(judas.cleaned.HR[, .(HRlat, Longitude)],
                                        judas.cleaned.HR[, .(HRlat, HRlong)])]
judas.cleaned.HR[, summary(xdev)]
judas.cleaned.HR[, summary(ydev)]

ggplot(judas.cleaned.HR) + geom_density(aes(xdev), col="blue") + 
  geom_density(aes(ydev), col="red") + xlim(c(0, 60))

# rm animals that didn't move (hence possibly wrong coords)
dev.checks <- judas.cleaned.HR[, .(sumxdev=sum(xdev), sumydev=sum(ydev)), by=JUDAS_ID]
dev.checks[sumxdev == 0 | sumydev == 0,]
setkey(judas.cleaned.HR, JUDAS_ID)
judas.cleaned.HR <- judas.cleaned.HR[dev.checks[sumxdev > 0 & sumydev > 0, JUDAS_ID], ]

# Check whether there are judas with < minloc data points and rm
locs <- judas.cleaned.HR[, .N, by=JUDAS_ID]
minloc <- 10
locs[, sum(N < minloc)]
judas.cleaned.HR <- judas.cleaned.HR[!JUDAS_ID %in% locs[N < minloc, JUDAS_ID], ]

judas.cleaned.HR[, summary(xdev)]
judas.cleaned.HR[, summary(ydev)]

ggplot(judas.cleaned.HR) + geom_density(aes(xdev), col="blue") + 
  geom_density(aes(ydev), col="red") + xlim(c(0, 60))


# Retain only IDs that were kept for HR analysis
judas.cleaned.surv <- judas.cleaned.surv[JUDAS_ID %in% judas.cleaned.HR[, JUDAS_ID]]

ntrack.events.HR <- judas.cleaned.HR[, .N, by=JUDAS_ID]
ntrack.events.HR[, summary(N)]

ntrack.events.surv <- judas.cleaned.surv[, .N, by=JUDAS_ID]
ntrack.events.surv[, summary(N)]
###########################################################################

# Adjust start.date and end.date to stay within the year interval considered
judas.cleaned.surv[, start.date := ifelse(test = Region == "KIMBERLEY" & year(start.date) < 1995,
                            "1995-03-01", start.date)]

judas.cleaned.surv[, end.date := ifelse(test = Region == "KIMBERLEY" & year(end.date) > 2007,
                                          "2007-11-30", end.date)]

judas.cleaned.surv[, start.date := ifelse(test = Region == "PILBARA" & year(start.date) < 2007,
                                          "2007-04-01", start.date)]

judas.cleaned.surv[, end.date := ifelse(test = Region == "PILBARA" & year(end.date) > 2010,
                                        "2010-12-31", end.date)]

names(judas.cleaned.surv)
setkey(judas.cleaned.surv, JUDAS_ID)
setkey(judas.cleaned.HR, JUDAS_ID)

# List all info re each judas
list.IDs <- judas.cleaned.surv[, unique(JUDAS_ID)]
judas.info <- judas.cleaned.surv[list.IDs, .SD, 
                            .SDcols=c("JUDAS_ID", "Region", "start.date", "end.date"),
                            mult="first"]
setkey(judas.info, JUDAS_ID)
HR.info <- judas.cleaned.HR[list.IDs, .SD, 
                                 .SDcols=c("JUDAS_ID", "HRlat", "HRlong"  ),
                                 mult="first"]
judas.info <- merge(judas.info, HR.info, all.x = T)

# Find all possible associations
# RcppAlgos::comboGeneral is much faster than combn
combs <- RcppAlgos::comboGeneral(judas.cleaned.surv[, unique(JUDAS_ID)], 2)

# Collate info together
judas.Dist<- data.table(combs)
setnames(judas.Dist, c("JUDAS_ID", "id2"))
setkey(judas.Dist, JUDAS_ID)

judas.Dist <- merge(judas.Dist, judas.info, all.x = T)
nms <- c(names(judas.Dist)[c(-1, -2)])
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
key(judas.cleaned.surv) # Confirm that Judas_ID is the key

judas.Dist[, Concurrent := 
             ifelse(end.date.1 < start.date.2 | end.date.2 < start.date.1, 'no', "yes")]
setkey(judas.Dist, Concurrent)
judas.Dist <- judas.Dist["yes", ] # get rid of judas pairs that were not present at same time

# Keep only pairs from the same region
judas.Dist[, RegCheck:=Region.1 == Region.2]
judas.Dist <- judas.Dist[(RegCheck),]

# Prepare columns
# event is 0 when there is no encounter, and 1 when there is
# time is the length of time when both were 'deployed' (the beginning of the 
   # period when concurrently present) and when they were observed together if event=1
   # otherwise the max time when both were 'deployed'  
# maxtime is the max time when both were 'deployed'  
judas.Dist[, ':='(TrackedLoc=as.character(NA), EVENT_DATE=as.Date(NA), 
                  event=as.numeric(NA), time=as.numeric(NA), maxtime=as.numeric(NA), 
                  EVENT_ID=as.numeric(NA), EVENT_CODE=as.character(NA), 
                  ACTION=as.character(NA), 
                  LAT=as.numeric(NA), LONG=as.numeric(NA), Habitat.Type=as.character(NA),
                  Soil.Type=as.character(NA), 
                  Capacity=as.character(NA))]
jd.dates <- c("start.date.1", "end.date.1", "start.date.2", "end.date.2")
judas.Dist[, (jd.dates) := lapply(.SD, as.Date), .SDcols=jd.dates]

# summary stats
judas.Dist[, summary(Dist)]

# Apply a cut off of 100 km
judas.Dist <- judas.Dist[Dist < 100, ]

rns <- nrow(judas.Dist) # Number of judas pairs

# There must be a clever way to do this, but couldn't think any...
# This takes about 1.5 mins, an lapply and parallel execution would be faster, 
# but I didn't bother. Worth considering for larger datasets though.
system.time(
for(rn in seq_len(rns)) {
  data.ID1 <- judas.cleaned.surv[judas.Dist[rn, ID.1], ]
  data.ID2 <- judas.cleaned.surv[judas.Dist[rn, ID.2], ]
  matching.Events <- data.ID1[, EVENT_ID] %in% data.ID2[, EVENT_ID]
  data.Evs <- data.ID1[matching.Events, ]
  if(nrow(data.Evs) > 0) {
    line.no <- which.min(data.Evs[, Date])
    judas.Dist[rn, ':='(TrackedLoc=data.Evs[line.no, TrackedShire], 
                        EVENT_DATE=data.Evs[line.no, Date], event=1, 
                        EVENT_ID=data.Evs[line.no, EVENT_ID], 
                        EVENT_CODE=data.Evs[line.no, EVENT], 
                        ACTION=data.Evs[line.no, ACTION], 
                        LAT=data.Evs[line.no, Latitude], LONG=data.Evs[line.no, Longitude], 
                        Habitat.Type=data.Evs[line.no, `Habitat Type`],
                        Soil.Type=data.Evs[line.no, `Soil Type`],
                        Capacity=data.Evs[line.no, Capacity])]
  } else {
    judas.Dist[rn, ':='(EVENT_DATE=min(end.date.1, end.date.2), event=0)]
  }
  judas.Dist[rn, ':='(time=EVENT_DATE - max(start.date.1, start.date.2), 
                      maxtime=min(end.date.1, end.date.2) - max(start.date.1, start.date.2))]
}
)

#### Balance id1 and id2 to similar counts ####
# Compute descriptive stats of distance when event=1
judas.Dist[event == 1, summary(Dist)]
ggplot(judas.Dist[event == 1, ], aes(Dist)) + geom_density() + xlim(c(0, 100))

max.Dist <- judas.Dist[event == 1, max(Dist)]
judas.Dist.Bal <- judas.Dist[Dist <= max.Dist,]
judas.Dist.Bal[, Rn:=seq_len(nrow(judas.Dist[Dist <= max.Dist,]))]

JUDAS_ID <- unique(c(judas.Dist.Bal[, ID.1], judas.Dist.Bal[, ID.2]))
ID.1cnt.table <- judas.Dist.Bal[, .(ID.1cnt=.N), by=ID.1]
ID.2cnt.table <- judas.Dist.Bal[, .(ID.2cnt=.N), by=ID.2]
setnames(ID.1cnt.table, "ID.1", "JUDAS_ID")
setnames(ID.2cnt.table, "ID.2", "JUDAS_ID")

counts <- data.table(JUDAS_ID)
counts <- merge(counts, ID.1cnt.table, on=JUDAS_ID, all.x=TRUE)
counts <- merge(counts, ID.2cnt.table, on=JUDAS_ID, all.x=TRUE)
counts[is.na(ID.1cnt), ID.1cnt:=0]
counts[is.na(ID.2cnt), ID.2cnt:=0]
counts[, Tot:=ID.1cnt + ID.2cnt]
counts[, prop:=ID.1cnt/Tot]
counts[, summary(prop)]

nrefinement <- 5
rf <- 1

# ~14 secs
system.time(
while(rf <= nrefinement & 
      (counts[, quantile(prop, probs = 0.025)] < 0.4 & 
       counts[, quantile(prop, probs = 0.975)] > 0.6)) {
  for(rn in nrow(counts):1) {
    if(counts[rn, ID.1cnt] - counts[rn, ID.2cnt] < -3) {
      num <- floor(counts[rn, Tot]/2) - counts[rn, ID.1cnt]
      inds <- judas.Dist.Bal[ID.2 %in% counts[rn, JUDAS_ID], Rn]
      ID.1s <- judas.Dist.Bal[inds, ID.1]
      keep <- ID.1s %in% counts[prop > 0.5, JUDAS_ID]
      if(sum(keep) == 0) next
      sinds <- inds[keep]
      if(length(sinds) > num) sinds <- sample(inds,size=num,replace=F)
      move2ID.2 <- judas.Dist.Bal[sinds, ID.1]
      move2ID.1 <- judas.Dist.Bal[sinds, ID.2]
      judas.Dist.Bal[sinds, ID.1:=move2ID.1]
      judas.Dist.Bal[sinds, ID.2:=move2ID.2]
    }
  }
  
  for(rn in 1:nrow(counts)) {
    if(counts[rn, ID.1cnt] - counts[rn, ID.2cnt] > 3) {
      num <- floor(counts[rn, Tot]/2) - counts[rn, ID.2cnt]
      inds <- judas.Dist.Bal[ID.1 %in% counts[rn, JUDAS_ID], Rn]
      ID.2s <- judas.Dist.Bal[inds, ID.2]
      keep <- ID.2s %in% counts[prop < 0.5, JUDAS_ID]
      if(sum(keep) == 0) next
      sinds <- inds[keep]
      if(length(sinds) > num) sinds <- sample(inds,size=num,replace=F)
      move2ID.2 <- judas.Dist.Bal[sinds, ID.1]
      move2ID.1 <- judas.Dist.Bal[sinds, ID.2]
      judas.Dist.Bal[sinds, ID.1:=move2ID.1]
      judas.Dist.Bal[sinds, ID.2:=move2ID.2]
    }
  }
  
  ID.1cnt.table <- judas.Dist.Bal[, .(ID.1cnt=.N), by=ID.1]
  ID.2cnt.table <- judas.Dist.Bal[, .(ID.2cnt=.N), by=ID.2]
  setnames(ID.1cnt.table, "ID.1", "JUDAS_ID")
  setnames(ID.2cnt.table, "ID.2", "JUDAS_ID")
  
  counts <- data.table(JUDAS_ID)
  counts <- merge(counts, ID.1cnt.table, on=JUDAS_ID, all.x=TRUE)
  counts <- merge(counts, ID.2cnt.table, on=JUDAS_ID, all.x=TRUE)
  counts[is.na(ID.1cnt), ID.1cnt:=0]
  counts[is.na(ID.2cnt), ID.2cnt:=0]
  counts[, Tot:=ID.1cnt + ID.2cnt]
  counts[, prop:=ID.1cnt/Tot]
  message(paste("completed refinement", rf))
  rf <- rf + 1
}
)

counts[, quantile(prop, probs = 0.025)] 
counts[, quantile(prop, probs = 0.975)]
counts[, summary(prop)]

# Subset judas for HR analysis that were kept in survival analysis as ID1 in Dist.Bal
judas.cleaned.HR <- judas.cleaned.HR[JUDAS_ID %in% judas.Dist.Bal[, unique(ID.1)],]

# Number of judas retained
judas.Dist.Bal[, length(unique(ID.1))]

write.csv(judas.Dist, file = file.path(data.path, "judas.Dist.csv"), row.names = F)
save(judas.Dist, file = file.path(data.path, "judas.Dist.rda"))
save(judas.cleaned.surv, file = file.path(data.path, "judas.cleaned.surv.rda"))
save(judas.cleaned.HR, file = file.path(data.path, "judas.cleaned.HR.rda"))
# Note that data regarding ID.1 and ID.2 may be swapped now
# To track these it should be possible to add a col "swab" and add 1 each time a 
# line is swabbed so that it is then possible to adjust other cols
save(judas.Dist.Bal, file = file.path(data.path, "judas.Dist.Bal.rda"))
