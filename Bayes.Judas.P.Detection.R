library(data.table)
library(ggplot2)
library(jagsUI)
library(dplyr)

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
# Number of all retained observation of each judas
nobs<- judas.cleaned.PB %>% group_by(JUDAS_ID) %>% summarise(n=n())
nobs<- filter(nobs, n>=10)
judas.tmp.PB<- filter(judas.cleaned.PB, JUDAS_ID %in% nobs$JUDAS_ID)
judas.tmp.PB<- data.table(judas.tmp.PB)

n <- length(judas.tmp.PB[, unique(JUDAS_ID)])
data <- list(N2=nrow(judas.tmp.PB), M=n, mean=c(0,0), 
             id2=judas.tmp.PB[, as.numeric(unclass(as.factor(JUDAS_ID)))], 
            dev=as.matrix(judas.tmp.PB[, .(xdev, ydev)]))

inits <- function(){list(mu.sigmax=0,r=rep(0.5,n),
                         mu.sigmay=0,mu.rho=0.5,sigx=1,sigy=1,
                         sig.rho=1)}

params<- c("mu.sigmax","mu.sigmay","a","sigx","sigy","b","rho","sigmax","sigmay")

# fit model to data using WinBUGS code
ni <- 2000
nb <- 1000
nthin <- 1
nc <- 3
np <- 8 # Number of CPUs
fit2.PB = jags(data, inits, params,  model.file="./Models/HRmodel3.txt", 
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
#   Fit model with overarching priors
inits <- function(){list(mu.sigmax.pop=runif(1, 0, 4),
                         tau.sigmax.pop=runif(1, 0, 4),
                         mu.sigmay.pop=runif(1, 0, 4),
                         tau.sigmay.pop=runif(1, 0, 4), 
                         mu.rho.pop=runif(1, -1, 1),
                         sigma.rho.pop=runif(1, 0, 0.5))}

params<- c("mu.sigmax.pop", "tau.sigmax.pop", "mu.sigmay.pop", "tau.sigmay.pop",
           "mu.rho.pop", "sigma.rho.pop", "rho", "sigmax","sigmay")

# fit model to data using WinBUGS code
ni <- 3000
nb <- 1000
nthin <- 1
nc <- 3
np <- 8 # Number of CPUs
fit2.PB2 <- jags(data, inits, params,  model.file="./Models/HRmodel.vcov.hyperpriors.txt", 
                n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, #n.adapt = 2000, 
                parallel=ifelse(nc>1, TRUE, FALSE), 
                n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit2.PB2)
print(fit2.PB2, digits=3) 
sapply(fit2.PB2$n.eff, summary)
sapply(fit2.PB2$Rhat, summary)
plot(fit2.PB2)
traceplot(fit2.PB2, parameters="rho")
traceplot(fit2.PB2, parameters="sigmay")
save(fit2.PB2, file = file.path(analysis.path, "fit2.PB.rda"))
load(file.path(analysis.path, "fit2.PB2.rda"))

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

n <- length(judas.tmp.KM[, unique(JUDAS_ID)])
data <- list(N2=nrow(judas.tmp.KM), M=n, mean=c(0,0), 
             id2=judas.tmp.KM[, as.numeric(unclass(as.factor(JUDAS_ID)))], 
             dev=as.matrix(judas.tmp.KM[, .(xdev, ydev)]))

inits <- function(){list(mu.sigmax=0,
                         mu.sigmay=0,mu.rho=0.5,sigx=1,sigy=1,
                         sig.rho=1)}

params<- c("mu.sigmax","mu.sigmay","a","sigx","sigy","b","rho","sigmax","sigmay")

# fit model to data using WinBUGS code
ni <- 2000
nb <- 1000
nthin <- 1
nc <- 3
np <- 8 # Number of CPUs
fit2.KM = jags(data, inits, params,  model.file="./Models/HRmodel3.txt", 
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
#   Fit model with overarching priors
inits <- function(){list(mu.sigmax.pop=runif(1, 0, 4),
                         tau.sigmax.pop=runif(1, 0, 4),
                         mu.sigmay.pop=runif(1, 0, 4),
                         tau.sigmay.pop=runif(1, 0, 4), 
                         mu.rho.pop=runif(1, -1, 1),
                         sigma.rho.pop=runif(1, 0, 0.5))}

params<- c("mu.sigmax.pop", "tau.sigmax.pop", "mu.sigmay.pop", "tau.sigmay.pop",
           "mu.rho.pop", "sigma.rho.pop", "rho", "sigmax","sigmay")

# fit model to data using WinBUGS code
ni <- 5000
nb <- 1000
nthin <- 1
nc <- 4
np <- 8 # Number of CPUs
fit2.KM2 <- jags(data, inits, params,  model.file="./Models/HRmodel.vcov.hyperpriors.txt", 
                 n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, #n.adapt = 2000, 
                 parallel=ifelse(nc>1, TRUE, FALSE), 
                 n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit2.KM2)
print(fit2.KM2, digits=3) 
sapply(fit2.KM2$n.eff, summary)
sapply(fit2.KM2$Rhat, summary)
plot(fit2.KM2)
traceplot(fit2.KM2, parameters="rho")
traceplot(fit2.KM2, parameters="sigmay")
save(fit2.KM2, file = file.path(analysis.path, "fit2.KM2.rda"))
load(file.path(analysis.path, "fit2.KM2.rda"))


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
#   Fit model with overarching priors
inits <- function(){list(mu.sigmax.pop=runif(1, 0, 4),
                         tau.sigmax.pop=runif(1, 0, 4),
                         mu.sigmay.pop=runif(1, 0, 4),
                         tau.sigmay.pop=runif(1, 0, 4), 
                         mu.rho.pop=runif(1, -1, 1),
                         sigma.rho.pop=runif(1, 0, 0.5))}

params<- c("mu.sigmax.pop", "tau.sigmax.pop", "mu.sigmay.pop", "tau.sigmay.pop",
           "mu.rho.pop", "sigma.rho.pop", "rho", "sigmax","sigmay")

# fit model to data using WinBUGS code
ni <- 6000
nb <- 1000
nthin <- 1
nc <- 4
np <- 8 # Number of CPUs
fit2.2 <- jags(data, inits, params,  model.file="./Models/HRmodel.vcov.hyperpriors.txt", 
                 n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, #n.adapt = 2000, 
                 parallel=ifelse(nc>1, TRUE, FALSE), 
                 n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit2.2)
print(fit2.2, digits=3) 
sapply(fit2.2$n.eff, summary)
sapply(fit2.2$Rhat, summary)
plot(fit2.2)
traceplot(fit2.2, parameters="rho")
traceplot(fit2.2, parameters="sigmay")
save(fit2.2, file = file.path(analysis.path, "fit2.2.rda"))
load(file.path(analysis.path, "fit2.2.rda"))

#------------------------------------------------------------------------------#

#### multidimensional integration ####


library(cubature) # The package pracma may be faster

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

#------------------------------------------------------------------------------#
#### speed test ####
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

sapply(t_int, "[[", 3)
#> 3.64 324.21   0.45  41.64   0.44

#### Integration for PB ####
b1<- fit1.PB2$mean$b1
b2<- fit1.PB2$mean$b2
sigmax<- fit2.PB$mean$sigmax 
sigmay<- fit2.PB$mean$sigmay
rho <- fit2.PB$mean$rho


system.time(
  pdet.PB <- mapp_int(b1, b2, sigmax, sigmay, rho)
)
ggplot(pdet.PB, aes(DProb)) + geom_density() 
ggsave(filename = file.path(analysis.path, "DProb.density.PB.pdf"))
quantile(pdet.PB[, 1], probs = c(0.025, 0.975))
# 2.5%       97.5% 
#  0.007055694 0.070258962 

# Using fit2.PB2
analysis.path <- "../Data/Analysis"
load(file = file.path(analysis.path, "fit1.PB2.rda"))
load(file.path(analysis.path, "fit2.PB2.rda"))

b1<- fit1.PB2$mean$b1
b2<- fit1.PB2$mean$b2
sigmax<- fit2.PB2$mean$sigmax 
sigmay<- fit2.PB2$mean$sigmay
rho <- fit2.PB2$mean$rho

system.time(
  pdet.PB2 <- mapp_int(b1, b2, sigmax, sigmay, rho)
)

ggplot(pdet.PB2, aes(DProb)) + geom_density() 
ggsave(filename = file.path(analysis.path, "DProb.density.PB2.pdf"))
quantile(pdet.PB2[, 1], probs = c(0.025, 0.975))
# 2.5%       97.5% 
#  0.006474863 0.058781619 

pdet.PB2.mu <- mapp_int(b1=fit1.PB2$mean$mu.b1, b2=fit1.PB2$mean$mu.b2, 
                 sigmax=fit2.PB2$mean$mu.sigmax.pop, 
                 sigmay=fit2.PB2$mean$mu.sigmay.pop, 
                 rho=fit2.PB2$mean$mu.rho.pop) 
pdet.PB2.mu
# DProb     Abs_error
#1 0.02106019 2.065035e-310

#### Integration for KM ####
analysis.path <- "../Data/Analysis"
load(file = file.path(analysis.path, "fit1.KM2.rda"))
load(file.path(analysis.path, "fit2.KM.rda"))

b1<- fit1.KM2$mean$b1
b2<- fit1.KM2$mean$b2
sigmax<- fit2.KM$mean$sigmax 
sigmay<- fit2.KM$mean$sigmay
rho <- fit2.KM$mean$rho

system.time(
  pdet.KM <- mapp_int(b1, b2, sigmax, sigmay, rho)
)

ggplot(pdet.KM, aes(DProb)) + geom_density() 
ggsave(filename = file.path(analysis.path, "DProb.density.KM.pdf"))



# Using fit2.KM2
analysis.path <- "../Data/Analysis"
load(file = file.path(analysis.path, "fit1.KM2.rda"))
load(file.path(analysis.path, "fit2.KM2.rda"))

b1<- fit1.KM2$mean$b1
b2<- fit1.KM2$mean$b2
sigmax<- fit2.KM2$mean$sigmax 
sigmay<- fit2.KM2$mean$sigmay
rho <- fit2.KM2$mean$rho

system.time(
  pdet.KM2 <- mapp_int(b1, b2, sigmax, sigmay, rho)
)

ggplot(pdet.KM2, aes(DProb)) + geom_density() 
ggsave(filename = file.path(analysis.path, "DProb.density.KM2.pdf"))
quantile(pdet.KM2[, 1], probs = c(0.025, 0.975))
#        2.5%       97.5% 
# 0.008297998 0.063008072 

pdet.KM2.mu <- mapp_int(b1=fit1.KM2$mean$mu.b1, b2=fit1.KM2$mean$mu.b2, 
                        sigmax=fit2.KM2$mean$mu.sigmax.pop, 
                        sigmay=fit2.KM2$mean$mu.sigmay.pop, 
                        rho=fit2.KM2$mean$mu.rho.pop) 
pdet.KM2.mu
#        DProb     Abs_error
# 0.03877888 2.300748e-310

#### Integration for whole dataset ####
analysis.path <- "../Data/Analysis"
load(file = file.path(analysis.path, "fit1.rda"))
load(file.path(analysis.path, "fit2.rda"))


#============================================================================

quantile(pdet[, 1], probs = c(0.025, 0.975))
#   2.5%       97.5% 
#  0.008305563 0.069610217 

# Using fit2.2
analysis.path <- "../Data/Analysis"
load(file = file.path(analysis.path, "fit1.rda"))
load(file.path(analysis.path, "fit2.2.rda"))

b1<- fit1$mean$b1
b2<- fit1$mean$b2
sigmax<- fit2.2$mean$sigmax 
sigmay<- fit2.2$mean$sigmay
rho <- fit2.2$mean$rho

system.time(
  pdet.2 <- mapp_int(b1, b2, sigmax, sigmay, rho)
)

ggplot(pdet.2, aes(DProb)) + geom_density() 
ggsave(filename = file.path(analysis.path, "DProb.density.2.pdf"))
quantile(pdet.2[, 1], probs = c(0.025, 0.975))
# 2.5%       97.5% 
#  0.008043299 0.061871118  

pdet.2.mu <- mapp_int(b1=fit1$mean$mu.b1, b2=fit1$mean$mu.b2, 
                        sigmax=fit2.2$mean$mu.sigmax.pop, 
                        sigmay=fit2.2$mean$mu.sigmay.pop, 
                        rho=fit2.2$mean$mu.rho.pop) 
pdet.2.mu
# DProb     Abs_error
# 0.03653901 2.563781e-310
