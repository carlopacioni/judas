options(java.parameters = "-Xmx8024m")
library(XLConnect, quietly = T)
library(ggplot2, quietly = T)
library(data.table, quietly = T)

data.path <- "../Data/"
judas.master <- data.table(readWorksheetFromFile(
  file = file.path(data.path, "Tracking_History_Judas_with habitat_MZ.xlsx"), 
  sheet="Judas_Tracking_History"))
names(judas.master)

# Rm useless spaces in REGION
judas.master[, REGION := sub(pattern = " +", replacement = "", x = REGION)]

# Create year and month columns
judas.master[, Year:=year(EVENT_DATE)]
judas.master[, Month:=month(EVENT_DATE)]

# Insert a col with location (shire) where the animals have been tracked, e.g. 'TrackedShire'
judas.master[, TrackedShire := *** your intercept function here ***]

# Search and replace all SHIRE instances for TrackedShire and things should work

# Count instances
runs <- judas.master[, .N, by=c("EVENT_DATE", "Year", "Month", "REGION", "SHIRE")]
ncollared <- judas.master[ACTION == "COLLARED", .(Ncollared=.N),
                          by=c("EVENT_DATE", "Year", "Month", "REGION", "SHIRE")]
# Meerge together
runs <- merge(runs, ncollared, all.x=T, 
              on=c("EVENT_DATE", "Year", "Month", "REGION", "SHIRE"))

runs_month <- runs[, .(.N, Ncollared=sum(Ncollared, na.rm = T)), 
                   by=c("Year", "Month", "REGION", "SHIRE")]

runs_month <- runs_month[, Effort := N * 7.5 - Ncollared * 0.5]
surveyed_months <- runs_month[, .N, by=c("REGION", "SHIRE", "Month")]

# Number of year of the program per region
runs_month[, .(Min=min(Year), Max=max(Year), nYear=max(Year) - min(Year)), 
           by=REGION]

ggplot(surveyed_months, aes(factor(Month), N)) + geom_histogram(stat="identity") + 
  facet_grid(SHIRE~REGION) + geom_hline(yintercept = 19, col="red") + 
  geom_hline(yintercept = 23, col="blue")

ggplot(data = runs_month, aes(x=factor(Month), y=Year)) +
  geom_raster(aes(fill = N), interpolate = F) +
  scale_fill_gradient(low = "yellow", high = "red") + facet_grid(SHIRE~REGION)

write.csv(runs_month, file.path(data.path, "Effort.csv"), row.names = F)  
