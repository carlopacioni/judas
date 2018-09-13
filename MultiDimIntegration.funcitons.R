library(cubature) # The package pracma may be faster

# 2D
integrand2d1<- function(xy, b1, b2, p) {
  # ellipse
  require(mvtnorm)
  d2<- (xy[1]^2 + xy[2]^2)
  d1<- sqrt(d2)
  covar<- prod(p)
  covmat<- matrix(c(p[1]^2,covar,covar,p[2]^2),2,2)  
  e<- dmvnorm(xy,c(0,0),covmat) 
  int<- e * (1-exp(-exp(b1 + b2*d1)))
  return(int)
}

integrand2d2<- function(xy, b1, b2, p, centre, Window) {
  
  testx<- centre[1]+xy[1]
  testy<- centre[2]+xy[2]
  if(inside.owin(testx,testy,Window)) {
    d2<- (xy[1]^2 + xy[2]^2)
    d1<- sqrt(d2)
    xdev<- xy[1]/p[1]
    ydev<- xy[2]/p[2]
    covar<-  2 * p[3] * xdev * ydev
    e<- (1/(2*pi*p[1]*p[2]*sqrt((1-p[3]^2))))*exp(-(xdev^2 + ydev^2 - covar)/(2*(1-p[3]^2)))
    int<- e * (1-exp(-exp(b1 + b2*d1)))
  }
  else int<- 0
  return(int)
}

integrand2d3<- function(xy, b1, b2, p) {
  
  d2<- (xy[1]^2 + xy[2]^2)
  d1<- sqrt(d2)
  xdev<- xy[1]/p[1]
  ydev<- xy[2]/p[2]
  covar<-  2 * p[3] * xdev * ydev
  e<- (1/(2*pi*p[1]*p[2]*sqrt((1-p[3]^2))))*exp(-(xdev^2 + ydev^2 - covar)/(2*(1-p[3]^2)))
  int<- e * (1-exp(-exp(b1 + b2*d1)))
  
  return(int)
}

integrand2d3_bis<- function(xy, b1, b2, sigmax, sigmay, rho) {
  
  d2<- (xy[1]^2 + xy[2]^2)
  d1<- sqrt(d2)
  xdev<- xy[1]/sigmax
  ydev<- xy[2]/sigmay
  covar<-  2 * rho * xdev * ydev
  e<- (1/(2*pi*sigmax*sigmay*sqrt((1-rho^2))))*exp(-(xdev^2 + ydev^2 - covar)/(2*(1-rho^2)))
  int<- e * (1-exp(-exp(b1 + b2*d1)))
  return(int)
}

# wrap integrand2d3_bis with mapply
mapp_int <- function(b1, b2, sigmax, sigmay, rho,f=integrand2d3_bis, 
                     lowerLimit=c(-50,-50), upperLimit=c(50,50), tol = 1e-05, 
                     fDim = 2, maxEval = 5e05, absError = 1e-10, vectorInterface=FALSE,
                     SIMPLIFY=FALSE) {
  mapp_list <- mapply(FUN = hcubature, b1=b1, b2=b2, sigmax=sigmax, sigmay=sigmay, rho=rho, 
                      MoreArgs=list(f=f,lowerLimit=lowerLimit, upperLimit=upperLimit, 
                                    tol=tol, fDim=fDim, maxEval=maxEval, 
                                    absError=absError, vectorInterface=vectorInterface), 
                      SIMPLIFY=SIMPLIFY)
  
  p.det<- data.frame(DProb=t(sapply(mapp_list, "[[", 1))[, 1],
                     Abs_error=t(sapply(mapp_list, "[[", 1))[, 2])
  return(p.det)
}

# Tried a vectorised fun passing a matrix, which is supposed to be way faster, but it doesn't work
integrand2d3_v<- function(xy, b1, b2, sigmax, sigmay, rho) {
  
  d2<- (xy[, 1]^2 + xy[, 2]^2)
  d1<- sqrt(d2)
  xdev<- xy[, 1]/sigmax
  ydev<- xy[, 2]/sigmay
  covar<-  2 * rho * xdev * ydev
  e<- (1/(2*pi*sigmax*sigmay*sqrt((1-rho^2))))*exp(-(xdev^2 + ydev^2 - covar)/(2*(1-rho^2)))
  int<- e * (1-exp(-exp(b1 + b2*d1)))
  return(int)
}

vec_int <- function(b1, b2, sigmax=sigmax, sigmay=sigmay, rho=rho, limit=50, 
                    tol = 1e-05, fDim = 2, maxEval = 5e05, absError = 1e-10) {
  njudas <- length(b1)
  int <- hcubature(f=integrand2d3_v, lowerLimit=matrix(-limit, njudas, 2), 
                   upperLimit=matrix(limit, njudas, 2), 
                   b1=b1, b2=b2, sigmax=sigmax, sigmay=sigmay, rho=rho,
                   tol=tol, fDim=fDim, maxEval=maxEval, absError=absError, 
                   vectorInterface = TRUE)
  
  p.det<- data.frame(DProb=int$integral,
                     Abs_error=int$error)
  return(p.det)
}

