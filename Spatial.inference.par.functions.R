

spat.det.prob <- function(i, rr=rr, tmpzone,  hparms, dparms,
                          sample_size, cellsize, nperiod, 
                          CE=0.80, Pu=1) {
    
    cells<- sampleRandom(rr, size=sample_size, xy=TRUE) # random locations of HR centres
    dat<- data.frame(xbar=cells[,1], ybar=cells[,2])
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
    
    bvn.surf<- make.surface(dat,params,tmpzone, nperiod=nperiod, cellsize=cellsize, CE=CE, Pu=Pu)
    
    tmpsims<- unlist(bvn.surf$Table)
    tmpsims<- c(Sim=i, NJudas=sample_size, CellSize=cellsize, CE=CE, Pu=Pu, tmpsims)
    
    return(tmpsims)
  }
