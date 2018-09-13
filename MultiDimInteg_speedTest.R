analysis.path <- "../Data/Analysis"
load(file = file.path(analysis.path, "fit1.PB2.rda"))
load(file.path(analysis.path, "fit2.PB.rda"))
b1<- fit1.PB2$mean$b1
b2<- fit1.PB2$mean$b2
sigmax<- fit2.PB$mean$sigmax 
sigmay<- fit2.PB$mean$sigmay
rho <- fit2.PB$mean$rho
p.det <- matrix(NA, njudas, 2)

# limit njudas to test speed
njudas<-5
t_int <- vector("list", length = 6)
t_int[[1]] <- system.time(
  for(i in 1:njudas) {
    message(paste("doing judas No. ", i, sep=""))
    int <- hcubature(f=integrand2d1, lowerLimit=c(-50,-50), upperLimit=c(50,50), 
                     b1=b1[i], b2=b2[i], p=c(sigmax[i], sigmay[i], rho[i]),
                     tol = 1e-05, fDim = 2, maxEval = 5e05, absError = 1e-10)
    p.det[i, 1] <- int$integral[1]
    p.det[i, 2] <- int$error[1]
  } 
)

t_int[[2]] <- system.time(
  for(i in 1:njudas) {
    message(paste("doing judas No. ", i, sep=""))
    int <- pcubature(f=integrand2d1, lowerLimit=c(-50,-50), upperLimit=c(50,50), 
                     b1=b1[i], b2=b2[i], p=c(sigmax[i], sigmay[i], rho[i]),
                     tol = 1e-05, fDim = 2, maxEval = 5e05, absError = 1e-10)
    p.det[i, 1] <- int$integral[1]
    p.det[i, 2] <- int$error[1]
  } 
)

t_int[[3]] <- system.time(
  for(i in 1:njudas) {
    message(paste("doing judas No. ", i, sep=""))
    int <- hcubature(f=integrand2d3, lowerLimit=c(-50,-50), upperLimit=c(50,50), 
                     b1=b1[i], b2=b2[i], p=c(sigmax[i], sigmay[i], rho[i]),
                     tol = 1e-05, fDim = 2, maxEval = 5e05, absError = 1e-10)
    p.det[i, 1] <- int$integral[1]
    p.det[i, 2] <- int$error[1]
  } 
)

t_int[[4]] <- system.time(
  for(i in 1:njudas) {
    message(paste("doing judas No. ", i, sep=""))
    int <- pcubature(f=integrand2d3, lowerLimit=c(-50,-50), upperLimit=c(50,50), 
                     b1=b1[i], b2=b2[i], p=c(sigmax[i], sigmay[i], rho[i]),
                     tol = 1e-05, fDim = 2, maxEval = 5e05, absError = 1e-10)
    p.det[i, 1] <- int$integral[1]
    p.det[i, 2] <- int$error[1]
  } 
)

t_int[[5]] <- system.time(
  pdet <- mapp_int(b1[1:njudas], b2[1:njudas], sigmax[1:njudas], 
                   sigmay[1:njudas], rho[1:njudas])
)

# This doesn't work :-(
t_int[[6]] <- system.time(
  p.det <- vec_int(b1=b1[1:njudas], b2=b2[1:njudas], sigmax=sigmax[1:njudas], 
                   sigmay=sigmay[1:njudas], rho=rho[1:njudas],
                   tol = 1e-05, fDim = 2, maxEval = 5e05, absError = 1e-10)
)
names(t_int) <- c("For_d1", "For_d1_pcub", "For_d3", "For_d3_pcub", "mapply_d3", "vec_d3")
sapply(t_int, "[[", 3)
