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

# Count instances
runs <- judas.master[, .N, by=c("EVENT_DATE", "Year", "Month", "REGION", "SHIRE")]
runs_month <- runs[, .N, by=c("Year", "Month", "REGION", "SHIRE")]

ggplot(runs_month, aes(factor(Month), N)) + geom_histogram(stat="identity") + 
  facet_grid(Year~REGION)

write.csv(runs_month, file.path(data.path, "Effort.csv"), row.names = F)  
