
library(simecol)  # for neighborhood function
library(rgdal)
library(raster)
library(progress)

# Do Kimberly region 6

region<- readOGR("../Data/Shapefile/Kimberley_LGA.shp")
tmpzone<- region[6,]
tmpzone<- spTransform(tmpzone, CRS("+proj=utm +zone=51 +south +ellps=GRS80 +units=km +no_defs")) #km

# Get params of hierarchical dists
dparms<- fit1.KM$mean
hparms<- fit2.KM$mean

cellsize <- 1  # resolution of raster in km
nsims<- 10
njudas<- c(10,25,50,100)
results<- matrix(NA,nrow=length(njudas),ncol=9)

rr<- raster(tmpzone, resolution=cellsize)
rr<- rasterize(tmpzone, rr, 0)
plot(rr)

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
  results[i,1:3]<- c(mean(tmpsims[,1]),quantile(tmpsims[,1],c(0.25,0.975)))
  results[i,4:6]<- c(mean(tmpsims[,2]),quantile(tmpsims[,2],c(0.25,0.975)))
  results[i,7:9]<- c(mean(tmpsims[,3]),quantile(tmpsims[,3],c(0.25,0.975)))
  pb$terminate()
}
results<- data.frame(results)
names(results)<- c("Seu","Seu.l","Seu.u","Cov","Cov.l","Cov.u","SSe","SSe.l","SSe.u")


#------------------------------------------
#
# Single realisation for plotting
#

njudas<- c(10,25,50,100)
cellsize<- 1

rast<- list()
centers<- list()
Pd<- list()

rr<- raster(tmpzone, resolution=cellsize)
rr<- rasterize(tmpzone, rr, 0)

for(i in 1:length(njudas)) {
  cat(paste("doing sample size ",np[i],sep=""),"\n")
  
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
  
  bvn.surf<- make.surface(dat,params,tmpzone, nperiod=1, cellsize=cellsize, CE=0.95, Pu=1)
  
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
  mtext(paste("n=",np[i],sep=""),side=3,cex=1)
  box()
}

#===============================================================================
# Functions
#================================================================================
kernel2D<- function(parms, maxdim, nperiod, eps, volume=0.95) {
  # 2D kernel density function fully vectorised
  require(mvtnorm)
  calc.BN<- function(x,y, parm){
    covar<- parm$sigmax * parm$sigmay * parm$rho
    covmat<- matrix(c(parm$sigmax^2,covar,covar,parm$sigmay^2),2,2)  
    dmvnorm(cbind(x,y), c(0,0), covmat)
  }  
  calc.dist<- function(x,y){
    return(sqrt((x^2 + y^2))) 
  }  
  max.disp<- max(c(parms$sigmax,parms$sigmay)) * 4
  mx<-  ceiling(max.disp/eps)
  while(length(-mx:mx) >= maxdim) {mx<- mx - 1}  #make sure kernel size is compatible with raster
  gxy <- -mx:mx * eps
  zprob <- outer(gxy, gxy, FUN=calc.BN, parm=parms)
  zhr<- getvol(zprob, volume=volume) # truncate kernel
  dmat<- outer(gxy, gxy, FUN=calc.dist)
  lam<- exp(parms$b1 + parms$b2 * dmat) # rate
  lam<- lam * nperiod  # total rate is rate per period * number of periods (np)
  return(lam * zhr) # mask lam to kernel
}

#------------------------------------------------
getvol<- function(X,volume=0.95){
  # find isopleth of kernel that gives the specified volume
  ff<- function(x,lalpha){
    if(x >= lalpha & !is.na(x)) 1 else 0
  }
  dims<- dim(X)
  vec<- sort(as.vector(X),decreasing=T)
  px<- vec/sum(vec)
  cum.px<- cumsum(px) #cdf
  ind<- length(cum.px[cum.px<=volume])
  tmpvol<- cum.px[ind]
  tmpvol2<- cum.px[ind+1]
  interp<- (volume-tmpvol)/(tmpvol2-tmpvol)
  lalpha<- vec[ind]*interp + vec[ind+1]*(1-interp)
  tmp<- matrix(sapply(as.vector(X),ff,lalpha),dims[1],dims[2])
  return(tmp)
}

#====================================

make.surface<- function(dat, parms, shape, nperiod, cellsize, CE=0.95, Pu=1, verbose=F) {
  n<- nrow(dat)
  rast<- raster(shape, resolution=cellsize)
  rast<- rasterize(shape, rast, 0, background = -1)
  rast.mat<- getValues(rast, format="matrix")
  dims<- dim(rast.mat)
  for(i in 1:n){
    if(verbose) cat(paste("doing judas ",i," of ",n,sep=""),"\n")
      tmprast<- rast
      cells<- cellFromXY(tmprast, cbind(dat$xbar[i],dat$ybar[i]))
      tmprast[cells]<- 1
      wkern<- kernel2D(parms=parms[i,], maxdim=min(dims), nperiod=nperiod, eps=cellsize, volume=CE)
      tmpmat<- getValues(tmprast, format="matrix")
      tmpmat<- simecol::neighbors(tmpmat, state=1, wdist=wkern)
      rast.mat<- rast.mat + tmpmat
  }
  den<- raster(rast.mat, template=rast)
  den <- mask(den, shape)
  
  N<- length(Which(den >= 0, cells=T))
  n<- length(Which(den > 1e6, cells=T)) # sampled cells
  den<- calc(den, function(x){1-exp(-x)}) #Probability scale
  seu_avg<- cellStats(den, 'mean', na.rm=TRUE)
  SSe<- 1 - (1 - seu_avg * n/N)^(Pu/N * N)
  result<- list(seu_avg=round(seu_avg,3),Cov=round(n/N,3),SSe=round(SSe,3))
  list(Table=result, Raster=den)
}

