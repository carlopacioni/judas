library(dplyr)
library(rgdal)
library(rgeos)
library(sf)
library(lubridate)

wa<-st_read("../Data/Shapefile/.","AUS_adm1")
wa<- wa %>% filter(NAME_1=="Western Australia")
wa<- st_simplify(wa, dTolerance=0.009)

kimb_LGA<- st_read("../Data/Shapefile/.","Kimberley_LGA")
kimb_LGA<- kimb_LGA %>% st_transform(4326) # convert to Lat/Long

Pil_LGA<- st_read("../Data/Shapefile/.","Pilbara_LGA")
Pil_LGA<- Pil_LGA %>% st_transform(4326) # convert to Lat/Long

win.graph(10,10)
plot(st_geometry(wa))
plot(st_geometry(kimb_LGA),add=T)
plot(st_geometry(Pil_LGA), add=T)
text(st_coordinates(st_centroid(kimb_LGA)), labels=kimb_LGA$LGA_CODE, cex=0.6)
text(st_coordinates(st_centroid(Pil_LGA)), labels=Pil_LGA$LGA_CODE, cex=0.6)
#plot(st_geometry(tracked_sf[tracked_sf$Shire=="HC",]), add=T, cex=0.5)
plot(st_geometry(tracked_sf), add=T, cex=0.5)
# Remove locations that are outside the ROI

