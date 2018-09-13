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
# params data.frame(sigmax=sx,sigmay=sy,rho=rho,b1=b1,b2=b2)
# nperiod number of discrete time unit the judas are deployed for
# cellsize the size of the cell (resolution of raster) in km
# CE=0.95 proportion of the Judas home range coverage
# Pu minimum number of (occupied) cells that are considered a detection 
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
  n<- length(Which(den > 1e-6, cells=T)) # sampled cells
  den<- calc(den, function(x){1-exp(-x)}) #Probability scale
  seu_avg<- cellStats(den, 'mean', na.rm=TRUE)
  SSe<- 1 - (1 - seu_avg * n/N)^(Pu/N * N)
  result<- list(seu_avg=round(seu_avg,3),Cov=round(n/N,3),SSe=round(SSe,3))
  list(Table=result, Raster=den)
}

