library(data.table)
library(ggplot2)
library(jagsUI)

#### Helper functions ####

# b1, b2 vectors of coefficients
# Dist vector of distances of same length as b1 & b2
compute.hDist <- function(b1, b2, Dist) {
  h_i <- 1-exp(-exp(b1 + b2 * Dist))
  return(h_i)
}

# M number of judas
# fit the model output 
compute.hJudas <- function (fit, M) {
  lhJudas <- list()
  for(i in seq_len(M)) {
    h <- 1-exp(-exp(fit$sims.list$b1[,i]))
    Judas <- rep(i, length(h))
    lhJudas[[i]] <- data.frame(Judas, h) 
  }
  hJudas <- do.call('rbind', lhJudas)
  return(hJudas)
}

#------------------------------------------------------------------------------#


data.path <- "../Data/"
load(file = file.path(data.path, "judas.cleaned.sub.rda"))
load(file = file.path(data.path, "judas.Dist.Bal.rda"))

#### Judas association probability - Discrete survival analysis ####
# Compute descriptive stats of distance when event=1
judas.Dist.Bal[event == 1, summary(Dist)]
judas.Dist.Bal[, summary(Dist)]
ggplot(judas.Dist.Bal[event == 1, ], aes(Dist)) + geom_density() + xlim(c(0, 100))
ggplot(judas.Dist.Bal, aes(Dist)) + geom_density() + xlim(c(0, 100))
hist(judas.Dist.Bal[, time])

#### Fit model to PILBARA donkeys ####
judas.cleaned.PB <- judas.cleaned.sub[Region == "PILBARA",]
judas.Dist.Bal.PB <- judas.Dist.Bal[Region.1 == "PILBARA",]
judas.Dist.Bal.PB[event == 1, summary(Dist)]
judas.Dist.Bal.PB[event == 0, summary(Dist)]  
hist(judas.Dist.Bal.PB[, time])

# Number of all possible judas pairs concurrently deployed 
npairs.PB <- nrow(judas.Dist.Bal.PB)

# Number of all retained observation of each judas
N2.PB <- nrow(judas.cleaned.PB)
JUDAS_ID.PB <- judas.Dist.Bal.PB[, unique(ID.1)] 
judas_id_cleaned.PB <- judas.cleaned.PB[, unique(JUDAS_ID)]

# y is the matrix of times and events taking d=1 for an event and d=0 for censoring;
# N1 is the number of observations (pairs); 
# M is the number of individuals
# time is in (integer) months

N1 <- npairs.PB
M <- length(JUDAS_ID.PB)
id1 <- judas.Dist.Bal.PB[, as.numeric(unclass(as.factor(ID.1)))]
time <- judas.Dist.Bal.PB[, floor(time / (30.4 * 2)) + 1]
d <- judas.Dist.Bal.PB[, event]
distance <- judas.Dist.Bal.PB[, Dist]

y <- matrix(rep(NA, N1 * max(time)), nrow=N1)
for (i in 1:N1) {
  y[i, time[i]] <- d[i]   
  for (j in seq_len(time[i] - 1)) {
    y[i,j] <- 0
  }
}

t <- apply(y, 1, sum, na.rm=TRUE)
sum(t)
sum(d)
dim(y)


data.PB <- list(y=y, N1=N1, M=M, id1=id1, time=time, distance=distance)
sapply(data.PB, class)
hist(data.PB$time)


# b1.init <- rnorm(length(JUDAS_ID.PB), sd = 0.5)
# b2.init <- rnorm(length(JUDAS_ID.PB), sd = 0.5)

inits <- function(){list(mu.b1= -1, mu.b2=0, sigma.b1=1, sigma.b2=1) 
                        }

params <- c("mu.b1","mu.b2","sigma.b1","sigma.b2","b1","b2")

ni <- 20000
nb <- 10000
nthin <- 1
nc <- 3
np <- 8 # Number of CPUs

fit1.PB = jags(data.PB, inits, params,  model.file="./Models/SurvDist2.txt", 
                n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin,  
                parallel=ifelse(nc>1, TRUE, FALSE), 
               n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit1.PB)
print(fit1.PB,digits=3) 
sapply(fit1.PB$n.eff, summary)
sapply(fit1.PB$Rhat, summary)
plot(fit1.PB)

# fit1.PB$sims.list

analysis.path <- "../Data/Analysis"
save(fit1.PB, file = file.path(analysis.path, "fit1.PB.rda"))
load(file = file.path(analysis.path, "fit1.PB.rda"))

# Prob profile
hDist <- mapply(compute.hDist, b1=fit1.PB$mean$b1, b2=fit1.PB$mean$b2, 
             Dist=seq(0, 50, length.out=length(fit1.PB$mean$b1)))

pProf <- data.frame(h=hDist, Distance=seq(0, 50, length.out=length(fit1.PB$mean$b1)))
ggplot(pProf, aes(Distance, h)) + geom_point() + geom_smooth()

mu.b1.PB.backtrans <- 1-exp(-exp(fit1.PB$mean$mu.b1))
mu.b2.PB.backtrans <- 1-exp(-exp(fit1.PB$mean$mu.b2))
#------------------------------------------------------------------------------#

#### Fit model to KM donkeys ####
judas.cleaned.KM <- judas.cleaned.sub[Region == "KIMBERLEY",]
judas.Dist.Bal.KM <- judas.Dist.Bal[Region.1 == "KIMBERLEY",]
judas.Dist.Bal.KM[event == 1, summary(Dist)]
judas.Dist.Bal.KM[event == 0, summary(Dist)]  
# judas.Dist.Bal.KM <- judas.Dist.Bal.KM[Dist < 50,]

# Number of all possible judas pairs concurrently deployed 
npairs.KM <- nrow(judas.Dist.Bal.KM)

# Number judas
JUDAS_ID.KM <- judas.Dist.Bal.KM[, unique(ID.1)] 

# y is the matrix of times and events taking d=1 for an event and d=0 for censoring;
# N1 is the number of observations (pairs); M is the number of individuals
# time is in (integer) months
N1 <- npairs.KM
M <- length(JUDAS_ID.KM)
id1 <- judas.Dist.Bal.KM[, as.numeric(unclass(as.factor(ID.1)))]
time <- judas.Dist.Bal.KM[, floor(time / (30.4 * 2)) + 1]
d <- judas.Dist.Bal.KM[, event]
distance <- judas.Dist.Bal.KM[, Dist]

y <- matrix(rep(NA, N1 * max(time)), nrow=N1)
for (i in 1:N1) {
  y[i, time[i]] <- d[i]   
  for (j in seq_len(time[i] - 1)) {
    y[i,j] <- 0
  }
}

t <- apply(y, 1, sum, na.rm=TRUE)
sum(t)
sum(d)
dim(y)

data.KM <- list(y=y, N1=N1, M=M, id1=id1, time=time, distance=distance)

sapply(data.KM, class)
hist(data.KM$time)
apply(y, 1, class)

# b1.init <- rnorm(length(JUDAS_ID.KM), sd = 0.5)
# b2.init <- rnorm(length(JUDAS_ID.KM), sd = 0.5)

inits <- function(){list(mu.b1= -1, mu.b2=0, sigma.b1=1, sigma.b2=1)}

params <- c("mu.b1","mu.b2","sigma.b1","sigma.b2", "b1","b2")


ni <- 20000
nb <- 10000
nthin <- 1
nc <- 3
np <- 8 # Number of CPUs
fit1.KM = jags(data.KM, inits, params,  model.file="./Models/SurvDist2.txt", 
               n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin,  
               parallel=ifelse(nc>1, TRUE, FALSE), 
               n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit1.KM)
print(fit1.KM,digits=3) 
sapply(fit1.KM$n.eff, summary)
sapply(fit1.KM$Rhat, summary)
plot(fit1.KM)

# fit1.KM$sims.list

analysis.path <- "../Data/Analysis"
save(fit1.KM, file = file.path(analysis.path, "fit1.KM.rda"))

# Prob profile
hDist <- mapply(compute.hDist, b1=fit1.KM$mean$b1, b2=fit1.KM$mean$b2, 
                Dist=seq(0, 50, length.out=length(fit1.KM$mean$b1)))

pProf <- data.frame(h=hDist, Distance=seq(0, 50, length.out=length(fit1.KM$mean$b1)))
ggplot(pProf, aes(Distance, h)) + geom_point() + geom_smooth()

mu.b1.KM.backtrans <- 1-exp(-exp(fit1.KM$mean$mu.b1))
mu.b2.KM.backtrans <- 1-exp(-exp(fit1.KM$mean$mu.b2))
#------------------------------------------------------------------------------#

#### Fit model to whole dataset ####

# Number of all possible judas pairs concurrently deployed 
npairs <- nrow(judas.Dist.Bal)

# Number of all retained observation of each judas
N2 <- nrow(judas.cleaned.sub)
JUDAS_ID <- judas.Dist.Bal[, unique(ID.1)] 

# y is the matrix of times and events taking d=1 for an event and d=0 for censoring;
# N1 is the number of observations (pairs); M is the number of individuals
# time is in (integer) months
N1 <- npairs
M <- length(JUDAS_ID)
id1 <- judas.Dist.Bal[, as.numeric(unclass(as.factor(ID.1)))]
time <- judas.Dist.Bal[, floor(time / (30.4 * 2)) + 1]
d <- judas.Dist.Bal[, event]
distance <- judas.Dist.Bal[, Dist]

y <- matrix(rep(NA, N1 * max(time)), nrow=N1)
for (i in 1:N1) {
  y[i, time[i]] <- d[i]   
  for (j in seq_len(time[i] - 1)) {
    y[i,j] <- 0
  }
}

t <- apply(y, 1, sum, na.rm=TRUE)
sum(t)
sum(d)
dim(y)

data <- list(y=y, N1=N1, M=M, id1=id1, time=time, distance=distance)

sapply(data, class)

# b1.init <- rnorm(length(JUDAS_ID), sd = 0.5)
# b2.init <- rnorm(length(JUDAS_ID), sd = 0.5)

inits <- function(){list(mu.b1= -1, mu.b2=0, sigma.b1=1, sigma.b2=1, 
                         rho.b=runif(1, -1, 1), 
                         B.hat=cbind(rnorm(length(JUDAS_ID), sd = 0.5), 
                                     rnorm(length(JUDAS_ID), sd = 0.5)))}

params <- c("mu.b1","mu.b2","sigma.b1","sigma.b2","rho.b", "b1","b2")


ni <- 20000
nb <- 10000
nthin <- 10
nc <- 3
np <- 8 # Number of CPUs
fit1 = jags(data, inits, params,  model.file="./Models/SurvDist.txt", 
               n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, n.adapt=5000, 
               parallel=ifelse(nc>1, TRUE, FALSE), 
               n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit1)
print(fit1,digits=3) 
sapply(fit1$n.eff, summary)
sapply(fit1$Rhat, summary)

# fit1$sims.list

analysis.path <- "../Data/Analysis"
save(fit1, file = file.path(analysis.path, "fit1.rda"))
load(file.path(analysis.path, "fit1.rda"))

# Prob profile
hDist <- mapply(compute.hDist, b1=fit1$mean$b1, b2=fit1$mean$b2, 
                Dist=seq(0, 50, length.out=length(fit1$mean$b1)))

pProf <- data.frame(h=hDist, Distance=seq(0, 50, length.out=length(fit1$mean$b1)))
ggplot(pProf, aes(Distance, h)) + geom_point() + geom_smooth()

mu.b1.backtrans <- 1-exp(-exp(fit1$mean$mu.b1))
mu.b2.backtrans <- 1-exp(-exp(fit1$mean$mu.b2))


#=============================================================================#

#### Bivariate HR model ####

# Fit to PB data
# dev are the deviates in the X and Y  dimensions; 
# N2 is the total number of observations (locations) for each judas; 
# M is the number of individuals; mean[] is given as data c(0,0)
# id2 is the index of individuals of the retained observation of judas
n <- length(judas.cleaned.PB[, unique(JUDAS_ID)])
data <- list(N2=N2.PB, M=n, mean=c(0,0), 
             id2=judas.cleaned.PB[, as.numeric(unclass(as.factor(JUDAS_ID)))], 
            dev=as.matrix(judas.cleaned.PB[, .(xdev, ydev)]))

inits <- function(){list(sigmax=runif(n, 0.1, 10), 
                         sigmay=runif(n, 0.1, 10), 
                         rho=runif(n, -1, 1))}

params<- c("rho","sigmax","sigmay")

# fit model to data using WinBUGS code
ni <- 2000
nb <- 1000
nthin <- 1
nc <- 3
np <- 8 # Number of CPUs
fit2.PB = jags(data, inits, params,  model.file="./Models/HRmodel.txt", 
            n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin,  
            parallel=ifelse(nc>1, TRUE, FALSE), 
            n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit2.PB)
print(fit2.PB,digits=3) 
sapply(fit2.PB$n.eff, summary)
sapply(fit2.PB$Rhat, summary)
plot(fit2.PB)
traceplot(fit2.PB, parameters="rho")
traceplot(fit2.PB, parameters="sigmay")
save(fit2.PB, file = file.path(analysis.path, "fit2.PB.rda"))
load(file.path(analysis.path, "fit2.PB.rda"))
#------------------------------------------------------------------------------#

# Fit to KM data
# dev are the deviates in the X and Y  dimensions; 
# N2 is the total number of observations (locations) for each judas; 
# M is the number of individuals; mean[] is given as data c(0,0)
# id2 is the index of individuals of the retained observation of judas

# Number of all retained observation of each judas
nobs<- judas.cleaned.KM %>% group_by(JUDAS_ID) %>% summarise(n=n())
nobs<- filter(nobs, n>=10)
judas.tmp.KM<- filter(judas.cleaned.KM, JUDAS_ID %in% nobs$JUDAS_ID)
judas.tmp.KM<- data.table(judas.tmp.KM)

N2.KM <- nrow(judas.tmp.KM)

n <- length(judas.tmp.KM[, unique(JUDAS_ID)])
data <- list(N2=N2.KM, M=n, mean=c(0,0), 
             id2=judas.tmp.KM[, as.numeric(unclass(as.factor(JUDAS_ID)))], 
             dev=as.matrix(judas.tmp.KM[, .(xdev, ydev)]))

inits <- function(){list(sigmax=runif(n, 0.1, 10), 
                         sigmay=runif(n, 0.1, 10), 
                         rho=runif(n, -1, 1))}

params<- c("sigmax","sigmay","rho")

# fit model to data using WinBUGS code
ni <- 200
nb <- 100
nthin <- 1
nc <- 3
np <- 8 # Number of CPUs
fit2.KM = jags(data, inits, params,  model.file="./Models/HRmodel.txt", 
               n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin,  
               parallel=ifelse(nc>1, TRUE, FALSE), 
               n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit2.KM)
print(fit2.KM,digits=3) 
sapply(fit2.KM$n.eff, summary)
sapply(fit2.KM$Rhat, summary)
plot(fit2.KM)
traceplot(fit2.KM, parameters="rho")
traceplot(fit2.KM, parameters="sigmay")
save(fit2.KM, file = file.path(analysis.path, "fit2.KM.rda"))
load(file.path(analysis.path, "fit2.KM.rda"))
#------------------------------------------------------------------------------#

# Fit to whole dataset
# dev are the deviates in the X and Y  dimensions; 
# N2 is the total number of observations (locations) for each judas; 
# M is the number of individuals; mean[] is given as data c(0,0)
# id2 is the index of individuals of the retained observation of judas

# Number of all retained observation of each judas
N2 <- nrow(judas.cleaned.sub)

n <- length(judas.cleaned.sub[, unique(JUDAS_ID)])
data <- list(N2=N2, M=n, mean=c(0,0), 
             id2=judas.cleaned.sub[, as.numeric(unclass(as.factor(JUDAS_ID)))], 
             dev=as.matrix(judas.cleaned.sub[, .(xdev, ydev)]))

inits <- function(){list(sigmax=runif(n, 0.1, 10), 
                         sigmay=runif(n, 0.1, 10), 
                         rho=runif(n, -1, 1))}

params<- c("sigmax","sigmay","rho")

# fit model to data using WinBUGS code
ni <- 20000
nb <- 1000
nthin <- 1
nc <- 1
np <- 8 # Number of CPUs
fit2 = jags(data, inits, params,  model.file="./Models/HRmodel.txt", 
               n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, #n.adapt = 2000, 
               parallel=ifelse(nc>1, TRUE, FALSE), 
               n.cores=ifelse(floor(nc/np) < np, nc, np))

fit2 = jags(data, inits, params,  model.file="./Models/HRmodel.vcov.txt", 
               n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, #n.adapt = 2000, 
               parallel=ifelse(nc>1, TRUE, FALSE), 
               n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit2)
print(fit2,digits=3) 
sapply(fit2$n.eff, summary)
sapply(fit2$Rhat, summary)
plot(fit2)
traceplot(fit2, parameters="rho")
traceplot(fit2, parameters="sigmay")
save(fit2, file = file.path(analysis.path, "fit2.KM.rda"))
load(file.path(analysis.path, "fit2.rda"))
#------------------------------------------------------------------------------#

#### multidimensional integration ####

# R2Cuba package for multidimensional integration
library(R2Cuba)
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
#------------------------------------------------------------------------------#

b1<- fit1$mean$b1
b2<- fit1$mean$b2
sigmax<- fit2$mean$sigmax 
sigmay<- fit2$mean$sigmay
rho <- fit2$mean$rho

njudas <- length(b1)

p.det <- matrix(NA, njudas, 4)

for(i in 1:njudas) {
  message(paste("doing judas No. ", i, sep=""))
  int <- cuhre(2, 1, integrand2d1, lower=c(-50,-50), upper=c(50,50),
              b1=b1[i], b2=b2[i], p=c(sigmax[i], sigmay[i], rho[i]),
              flags=list(verbose=0))
  p.det[i, 1] <- int$value
  p.det[i, 2] <- int$abs.error
  p.det[i, 3] <- int$prob
  p.det[i, 4] <- ifelse(int$prob>0.1,1,0)
} 

p.det<- data.frame(p.det)
names(p.det)<- c("DProb","Abs_error","chiprob","flag")
#-------------
# Use individual samples to incorporate uncertainty
#
# CAUTION - LONG RUN-TIME:  Best use parallel version

nsims<- 10
sims1<- fit1$sims.list
sims2<- fit2$sims.list

n.samples<- fit1$mcmc.info$n.samples
ss<- sample(1:n.samples,size=nsims,replace=F)

b1<- sims1$b1[ss,]
b2<- sims1$b2[ss,]
sigmax<- sims2$sigmax[ss,] 
sigmay<- sims2$sigmay[ss,]
rho<- sims2$rho[ss,]

dd<- dim(sigmax)

p.det<- matrix(NA,dd[1],dd[2])

for(i in 1:dd[1]) {
  message(paste("doing judas No. ",i,sep=""))
  for(j in 1:dd[2]) {
    int<- cuhre(2,1,integrand2d1, lower=c(-50,-50),upper=c(50,50),
                b1=b1[i,j],b2=b2[i,j],p=c(sigmax[i,j],sigmay[i,j],rho[i,j]),
                flags=list(verbose=0))
    p.det[i,j]<- int$value
  } }

#-----------------------------------------------------------
# Parallel implementation
################################################################################
####### Need to check code doSMP not available on CRAN #########################
################################################################################

workers <- startWorkers(12) # 10 cores
registerDoSMP(workers)
start.time<- Sys.time()

p.det2 <- foreach(i=1:dd[1], .combine=rbind, 
                 .packages=c("spatstat","R2Cuba")) %:% foreach(j=1:dd[2], .combine=c) %dopar% {
  centre <- c(goats$xbar[goats$ID == id[i]],goats$ybar[goats$ID == id[i]])
  cuhre(2,1,integrand2d2, lower=c(-50,-50),upper=c(50,50),b1=b1[i,j],b2=b2[i,j],p=c(sigmax[i,j],sigmay[i,j],rho[i,j]),centre=centre,Window=tmpzone,flags=list(verbose=0))$value
  
}
end.time = Sys.time()
elapsed.time = difftime(end.time, start.time, units='mins')
elapsed.time
stopWorkers(workers) 



#------------------------------------------------------------
require(doSMP)
workers <- startWorkers(12) # 10 cores
registerDoSMP(workers)
start.time<- Sys.time()

p.det2<- foreach(i=1:dd[1], .combine=rbind, .packages=c("spatstat","R2Cuba")) %:% foreach(j=1:dd[2], .combine=c) %dopar% {
  centre<- c(goats$xbar[goats$ID == id[i]],goats$ybar[goats$ID == id[i]])
  cuhre(2,1,integrand2d2, lower=c(-50,-50),upper=c(50,50),b1=b1[i,j],b2=b2[i,j],p=c(sigmax[i,j],sigmay[i,j],rho[i,j]),centre=centre,Window=tmpzone,flags=list(verbose=0))$value
  
}
end.time = Sys.time()
elapsed.time = difftime(end.time, start.time, units='mins')
elapsed.time
stopWorkers(workers) 
#=====================================================

require(doSMP)
workers <- startWorkers(12) # 10 cores
registerDoSMP(workers)
start.time<- Sys.time()

p.det<- foreach(i=1:dd[1], .combine=rbind, .packages=c("mvtnorm","R2Cuba")) %:% foreach(j=1:dd[2], .combine=c) %dopar% {
  
  cuhre(2,1,integrand2d3, lower=c(-50,-50),upper=c(50,50),b1=b1[i,j],b2=b2[i,j],p=c(sigmax[i,j],sigmay[i,j],rho[i,j]),flags=list(verbose=0))$value
  
}

end.time = Sys.time()
elapsed.time = difftime(end.time, start.time, units='mins')
elapsed.time
stopWorkers(workers) 
#rmSessions(all.names=TRUE)
#----------------------------------------------

par(mfrow=c(2,1))

n<- dim(nontrans)[2]
boxplot(nontrans[,1:n], notch=F, outline=F,ylim=c(0,1),col="grey80",las=1,ylab="Probability of association",xlab="Goat ID",axes=F)
axis(1, at=1:n, labels=as.character(jud.tmp$id2[which(jud.tmp$count==0)]),cex.axis=0.5,las=2)
axis(2, at=seq(0,1,0.2),las=1)
mtext("a",line=1)

n<- dim(trans)[2]
boxplot(trans[,1:n], notch=F, outline=F,ylim=c(0,1),col="grey80",las=1,ylab="Probability of association",xlab="Goat ID",axes=F)
axis(1, at=1:n, labels=as.character(jud.tmp$id2[which(jud.tmp$count==1)]),cex.axis=0.5,las=2)
axis(2, at=seq(0,1,0.2),las=1)
mtext("b",line=1) 


mean(apply(nontrans,2,mean))
mean(apply(trans,2,mean))


#============================================================================

