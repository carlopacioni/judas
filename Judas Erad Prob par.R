library(simecol)  # for neighborhood function
library(rgdal)
library(raster)
library(parallel)
library(data.table)
library(ggplot2)
source("./Spatial.inference.par.functions.R")

# load shape file (region 6 --> East Kimberley)
region<- readOGR("../Data/Shapefile/Kimberley_LGA.shp")
tmpzone<- region[6,]
tmpzone<- spTransform(tmpzone, CRS("+proj=utm +zone=51 +south +ellps=GRS80 +units=km +no_defs")) #km

# Get params of hierarchical dists
analysis.path <- "../Data/Analysis"
load(file = file.path(analysis.path, "fit1.KM2.rda"))
load(file.path(analysis.path, "fit2.KM.rda"))

dparms<- fit1.KM2$mean
hparms<- fit2.KM$mean

# Sim params
cellsizes <- c(1, 5, 10)  # resolution of raster in km
nsims<- 8
njudas<- c(10,50,100,500)
CE=0.9
Pus=c(1, 5, 10)
nperiod=6

raw_results<- vector("list", length = length(njudas) * length(cellsizes) * length(Pus))
for(cellsize in cellsizes) {
  cat(paste("doing cellsize ", cellsize, sep=""),"\n")
  # set up raster
  rr<- raster(tmpzone, resolution=cellsize)
  rr<- rasterize(tmpzone, rr, 0)
  for(Pu in Pus) {
    cat(paste("doing Pu ", Pu, sep=""),"\n")
    for(i in 1:length(njudas)) {
      
      cat(paste("doing sample size ",njudas[i],sep=""),"\n")
      sample_size <- njudas[i]
     
      # set up cluster
      ncores <- if(detectCores() > nsims) nsims else detectCores() 
      cl <- makeCluster(ncores)
      clusterSetRNGStream(cl, iseed = NULL)
      clusterEvalQ(cl, library("simecol"))
      clusterEvalQ(cl, library("rgdal"))
      clusterEvalQ(cl, library("raster"))
      clusterEvalQ(cl, source("./Spatial.inference.functions.R"))
      
      clusterExport(cl, 
                    varlist=c("rr", "tmpzone", "hparms", "dparms",
                              "sample_size", "cellsize", "nperiod", "CE", "Pu"), 
                    envir=.GlobalEnv) 
      st <-  system.time(
        tmpres <- 
          parLapply(cl, seq_len(nsims), spat.det.prob,rr=rr,
                    tmpzone=tmpzone,  hparms=hparms, dparms=dparms,
                    sample_size=sample_size, cellsize=cellsize, nperiod=nperiod, 
                    CE=CE, Pu=Pu) 
      )
      stopCluster(cl)
      print("Time elapsed")
      print(st)
      raw_results[[(which(cellsizes %in% cellsize) -1) * length(njudas) * length(Pus) +
                     (which(Pus %in% Pu) -1) * length(njudas) + i]] <-
        data.table(do.call(rbind, tmpres))
    }
  }
}
dt.raw_results <- rbindlist(raw_results)
results <- dt.raw_results[, c(Seu=mean(seu_avg), 
                              sapply(quantile(seu_avg,c(0.025, 0.5, 0.975)), list),
                              Cov=mean(Cov), 
                              sapply(quantile(Cov,c(0.025, 0.5, 0.975)), list), 
                              SSe=mean(SSe), 
                              sapply(quantile(SSe,c(0.025, 0.5, 0.975)), list)),
                          by=c("NJudas", "CellSize", "CE", "Pu")] 

pos <- grep("2.5%", x = names(results))
nnms <- sapply(c("Seu", "Cov", "SSe"), paste, c("lower", "median", "upper"), sep="_")
for(i in 0:2) {
setnames(results, pos + i, nnms[1 + i,])  
}
analysis.path <- "../Data/Analysis"
write.csv(dt.raw_results, file.path(analysis.path, "BootstrapRawResults.csv"), row.names = F)
write.csv(results, file.path(analysis.path, "BootstrapResults.csv"), row.names = F)

ggplot(results, aes(x=NJudas)) + geom_point(aes(y=SSe), col="black") + 
  geom_errorbar(aes(ymin=SSe_lower, ymax=SSe_upper), col="black") + 
  geom_point(aes(y=Cov), col="red") + 
  geom_errorbar(aes(ymin=Cov_lower, ymax=Cov_upper), col="red") +
  facet_grid(CellSize~Pu)
ggsave(file.path(analysis.path, "SimResPlots.tiff"), dpi = "print")
