
library(simecol)  # for neighborhood function
library(rgdal)
library(raster)
library(progress)
source("./Spatial.inference.functions.R")
# Do Kimberly region 6

region<- readOGR("../Data/Shapefile/Kimberley_LGA.shp")
tmpzone<- region[6,]
tmpzone<- spTransform(tmpzone, CRS("+proj=utm +zone=51 +south +ellps=GRS80 +units=km +no_defs")) #km

# Get params of hierarchical dists
dparms<- fit1.KM2$mean
hparms<- fit2.KM$mean

cellsize <- 1  # resolution of raster in km
nsims<- 10
njudas<- c(10,50,100,500)
results<- matrix(NA,nrow=length(njudas),ncol=9)

rr<- raster(tmpzone, resolution=cellsize)
rr<- rasterize(tmpzone, rr, 0)
plot(rr)
start_time <- Sys.time()
for(i in 1:length(njudas)) {
  tmpsims<- matrix(NA,nsims,3)
  cat(paste("doing sample size ",njudas[i],sep=""),"\n")
  pb<- progress_bar$new(total=nsims, width=60, show_after = 0)
  pb$tick(0)
  for(j in 1:nsims) {
    pb$tick()
    cells<- sampleRandom(rr,size=njudas[i],xy=TRUE) # random locations of HR centres
    dat<- data.frame(xbar=cells[,1],ybar=cells[,2])
    n<- nrow(dat)
    
    sx<- exp(rnorm(n, hparms$mu.sigmax, hparms$sigx))
    sx[sx < cellsize]<- cellsize
    sy<- exp(rnorm(n, hparms$mu.sigmay, hparms$sigy))
    sy[sy < cellsize]<- cellsize
    ro<- rbeta(n, hparms$a, hparms$b)
    rho<- 2 * (ro - 0.5)
    b1<- rnorm(n, dparms$mu.b1, dparms$sigma.b1)
    b2<- rnorm(n, dparms$mu.b2, dparms$sigma.b2)
    
    params<- data.frame(sigmax=sx,sigmay=sy,rho=rho,b1=b1,b2=b2)
    
    bvn.surf<- make.surface(dat,params,tmpzone, nperiod=12, cellsize=cellsize, CE=0.95, Pu=1)
    
    tmpsims[j,]<- unlist(bvn.surf$Table)
  }
  results[i,1:3]<- c(mean(tmpsims[,1]),quantile(tmpsims[,1],c(0.025,0.975)))
  results[i,4:6]<- c(mean(tmpsims[,2]),quantile(tmpsims[,2],c(0.025,0.975)))
  results[i,7:9]<- c(mean(tmpsims[,3]),quantile(tmpsims[,3],c(0.025,0.975)))
  pb$terminate()
}
end_time <- Sys.time()
end_time - start_time

results<- data.frame(results)
names(results)<- c("Seu","Seu.l","Seu.u","Cov","Cov.l","Cov.u","SSe","SSe.l","SSe.u")
analysis.path <- "../Data/Analysis"
write.csv(results, file.path(analysis.path, "ProbDetResults.csv"), row.names = F)

#------------------------------------------
#
# Single realisation for plotting
#

njudas<- c(10,50,100,500)
cellsize<- 1

rast<- list()
centers<- list()
Pd<- list()

rr<- raster(tmpzone, resolution=cellsize)
rr<- rasterize(tmpzone, rr, 0)

for(i in 1:length(njudas)) {
  cat(paste("doing sample size ",njudas[i],sep=""),"\n")
  
  cells<- sampleRandom(rr,size=njudas[i],xy=TRUE) 
  dat<- data.frame(xbar=cells[,1],ybar=cells[,2])
  centers[[i]]<- dat
  n<- nrow(dat)
  
  sx<- exp(rnorm(n, hparms$mu.sigmax, hparms$sigx))
  sx[sx < cellsize]<- cellsize
  sy<- exp(rnorm(n, hparms$mu.sigmay, hparms$sigy))
  sy[sy < cellsize]<- cellsize
  ro<- rbeta(n, hparms$a, hparms$b)
  rho<- 2 * (ro - 0.5)
  b1<- rnorm(n, dparms$mu.b1, dparms$sigma.b1)
  b2<- rnorm(n, dparms$mu.b2, dparms$sigma.b2)
  
  params<- data.frame(sigmax=sx,sigmay=sy,rho=rho,b1=b1,b2=b2)
  
  bvn.surf<- make.surface(dat,params,tmpzone, nperiod=12, cellsize=cellsize, CE=0.99, Pu=1)
  
  rast[[i]]<- bvn.surf$Raster
  Pd[[i]]<- bvn.surf$Table
}

Pd<- do.call('rbind',Pd)
#---------------------------------------------------------------------------------
win.graph(12,12)
pal.1=colorRampPalette(c("blue", "lightblue", "orange","red"), space="rgb")
min.rast<- min(unlist(lapply(rast,cellStats,'min')))
max.rast<- max(unlist(lapply(rast,cellStats,'max')))
par(mfrow=c(2,2),mar=c(2,2,2,2))
for(i in 1:4) {
  plot(rast[[i]],col=pal.1(100))
  plot(tmpzone,add=T)
  points(centers[[i]]$xbar,centers[[i]]$ybar,cex=0.6,pch=16,col="white")
  mtext(paste("n=",njudas[i],sep=""),side=3,cex=1)
  box()
}

