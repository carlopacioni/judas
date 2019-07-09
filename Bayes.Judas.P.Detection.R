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
load(file = file.path(data.path, "judas.Dist.Bal.rda"))

#### Judas association probability - Discrete survival analysis ####
# Compute descriptive stats of distance when event=1
judas.Dist.Bal[event == 1, summary(Dist)]
judas.Dist.Bal[, summary(Dist)]
ggplot(judas.Dist.Bal[event == 1, ], aes(Dist)) + geom_density() + xlim(c(0, 100))
ggplot(judas.Dist.Bal, aes(Dist)) + geom_density() + xlim(c(0, 100))
hist(judas.Dist.Bal[, time])

#### Fit model to PILBARA donkeys ####
judas.Dist.Bal.PB <- judas.Dist.Bal[Region.1 == "PILBARA",]
judas.Dist.Bal.PB[event == 1, summary(Dist)]
judas.Dist.Bal.PB[event == 0, summary(Dist)]  
hist(judas.Dist.Bal.PB[, time])

# y is the matrix of times and events taking d=1 for an event and d=0 for censoring;
# N1 is the number of observations (pairs); 
# M is the number of individuals
# time is in (integer) months

N1 <- nrow(judas.Dist.Bal.PB)
M <- length(judas.Dist.Bal.PB[, unique(ID.1)])
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

ni <- 40000
nb <- 10000
nthin <- 1
nc <- 4
np <- 8 # Number of CPUs

fit1.PB2 = jags(data.PB, inits, params,  model.file="./Models/SurvDist2.txt", 
                n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin,  
                parallel=ifelse(nc>1, TRUE, FALSE), 
               n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit1.PB2)
print(fit1.PB2,digits=3) 
sapply(fit1.PB2$n.eff, summary)
sapply(fit1.PB2$Rhat, summary)
plot(fit1.PB2)
fit1.PB2$mean[1:4]

hist(fit1.PB2$mean$b1)
hist(fit1.PB2$mean$b2)
summary(fit1.PB2$mean$b2)

ggplot(data.frame(b1=fit1.PB2$mean$b1, b2=fit1.PB2$mean$b2), aes(b1, b2)) +
  geom_point() + geom_smooth()
analysis.path <- "../Data/Analysis"
save(fit1.PB2, file = file.path(analysis.path, "fit1.PB2.rda"))

# Prob profile
hDist <- mapply(compute.hDist, b1=fit1.PB2$mean$b1, b2=fit1.PB2$mean$b2, 
                Dist=seq(0, 50, length.out=length(fit1.PB2$mean$b1)))

pProf <- data.frame(h=hDist, Distance=seq(0, 50, length.out=length(fit1.PB2$mean$b1)))
ggplot(pProf, aes(Distance, h)) + geom_point() + geom_smooth()
ggsave(file.path(analysis.path, "Density.dist.Dprob_Dist.PB.pdf"))

mu.b1.PB.backtrans <- 1-exp(-exp(fit1.PB2$mean$mu.b1))
mu.b2.PB.backtrans <- 1-exp(-exp(fit1.PB2$mean$mu.b2))
#------------------------------------------------------------------------------#

#### Fit model to KM donkeys ####
judas.Dist.Bal.KM <- judas.Dist.Bal[Region.1 == "KIMBERLEY",]
judas.Dist.Bal.KM[event == 1, summary(Dist)]
judas.Dist.Bal.KM[event == 0, summary(Dist)]  
# judas.Dist.Bal.KM <- judas.Dist.Bal.KM[Dist < 50,]


# y is the matrix of times and events taking d=1 for an event and d=0 for censoring;
# N1 is the number of observations (pairs); M is the number of individuals
# time is in (integer) 2-month periods
N1 <- nrow(judas.Dist.Bal.KM)
M <- length(judas.Dist.Bal.KM[, unique(ID.1)])
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

# b1.init <- rnorm(length(JUDAS_ID.KM), sd = 0.5)
# b2.init <- rnorm(length(JUDAS_ID.KM), sd = 0.5)

inits <- function(){list(mu.b1= -1, mu.b2=0, sigma.b1=1, sigma.b2=1)}

params <- c("mu.b1","mu.b2","sigma.b1","sigma.b2", "b1","b2")


ni <- 400000
nb <- 10000
nthin <- 10
nc <- 4
np <- 8 # Number of CPUs
fit1.KM2 = jags(data.KM, inits, params,  model.file="./Models/SurvDist2.txt", 
                n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, #n.adapt=100000, 
                parallel=ifelse(nc>1, TRUE, FALSE), 
                n.cores=ifelse(floor(nc/np) < np, nc, np))

summary(fit1.KM2)
print(fit1.KM2,digits=3) 
sapply(fit1.KM2$n.eff, summary)
sapply(fit1.KM2$Rhat, summary)
plot(fit1.KM2)
fit1.KM2$mean[1:4]

ggplot(data.frame(b1=fit1.KM2$mean$b1, b2=fit1.KM2$mean$b2), aes(b1, b2)) +
  geom_point() + geom_smooth()

hist(fit1.KM2$mean$b1)
hist(fit1.KM2$mean$b2)
summary(fit1.KM2$mean$b2)

analysis.path <- "../Data/Analysis"
save(fit1.KM2, file = file.path(analysis.path, "fit1.KM2.rda"))

# Prob profile
hDist <- mapply(compute.hDist, b1=fit1.KM2$mean$b1, b2=fit1.KM2$mean$b2, 
                Dist=seq(0, 50, length.out=length(fit1.KM2$mean$b1)))

pProf <- data.frame(h=hDist, Distance=seq(0, 50, length.out=length(fit1.KM2$mean$b1)))
ggplot(pProf, aes(Distance, h)) + geom_point() + geom_smooth()
ggsave(file.path(analysis.path, "Density.dist.Dprob_Dist.KM.pdf"))

mu.b1.KM.backtrans <- 1-exp(-exp(fit1.KM2$mean$mu.b1))
mu.b2.KM.backtrans <- 1-exp(-exp(fit1.KM2$mean$mu.b2))
#------------------------------------------------------------------------------#

#### Fit model to whole dataset ####

# y is the matrix of times and events taking d=1 for an event and d=0 for censoring;
# N1 is the number of observations (pairs); M is the number of individuals
# time is in (integer) months
N1 <- nrow(judas.Dist.Bal)
M <- length(judas.Dist.Bal[, unique(ID.1)])
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

inits <- function(){list(mu.b1= -1, mu.b2=0, sigma.b1=1, sigma.b2=1)}

params <- c("mu.b1","mu.b2","sigma.b1","sigma.b2", "b1","b2")

ni <- 400000
nb <- 10000
nthin <- 20
nc <- 4
np <- 8 # Number of CPUs
fit1 = jags(data, inits, params,  model.file="./Models/SurvDist2.txt", 
            n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, #n.adapt=5000, 
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
ggplot(pProf, aes(Distance, h)) + geom_point() + geom_smooth() + theme_classic() +
  ylab("Detection probability")
ggsave(file.path(analysis.path, "Density.dist.Dprob_Dist.CA.pdf"))
ggsave(file.path(analysis.path, "Density.dist.Dprob_Dist.CA.jpeg"))

mu.b1.backtrans <- 1-exp(-exp(fit1$mean$mu.b1))
mu.b2.backtrans <- 1-exp(-exp(fit1$mean$mu.b2))
#=============================================================================#

#### Bivariate HR model ####
data.path <- "../Data/"
load(file = file.path(data.path, "judas.cleaned.HR.rda"))

#### Fit to PB data ####
judas.cleaned.HR.PB <- judas.cleaned.HR[Region == "PILBARA",]
# dev are the deviates in the X and Y  dimensions; 
# N2 is the total number of observations (locations) for each judas; 
# M is the number of individuals; mean[] is given as data c(0,0)
# id2 is the index of individuals of the retained observation of judas
# Confirm min number of observations
nobs.PB <- judas.cleaned.HR.PB[, .N, by=JUDAS_ID]
nobs.PB[, summary(N)]
# setkey(judas.cleaned.HR.PB, JUDAS_ID)
# judas.cleaned.HR.PB <- judas.cleaned.HR.PB[nobs.PB[N>=5, JUDAS_ID], ]

N2.PB <- nrow(judas.cleaned.HR.PB)
n <- length(judas.cleaned.HR.PB[, unique(JUDAS_ID)])
data <- list(N2=N2.PB, M=n, mean=c(0,0), 
             id2=judas.cleaned.HR.PB[, as.numeric(unclass(as.factor(JUDAS_ID)))], 
             dev=as.matrix(judas.cleaned.HR.PB[, .(xdev, ydev)]))

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
analysis.path <- "../Data/Analysis"
save(fit2.PB, file = file.path(analysis.path, "fit2.PB.rda"))
load(file.path(analysis.path, "fit2.PB.rda"))
#------------------------------------------------------------------------------#

#### Fit to KM data ####
# dev are the deviates in the X and Y  dimensions; 
# N2 is the total number of observations (locations) for each judas; 
# M is the number of individuals; mean[] is given as data c(0,0)
# id2 is the index of individuals of the retained observation of judas

judas.cleaned.HR.KM <- judas.cleaned.HR[Region == "KIMBERLEY",]

# Confirm min number of observations
nobs.KM <- judas.cleaned.HR.KM[, .N, by=JUDAS_ID]
nobs.KM[, summary(N)]

N2.KM <- nrow(judas.cleaned.HR.KM)
n <- length(judas.cleaned.HR.KM[, unique(JUDAS_ID)])
data <- list(N2=N2.KM, M=n, mean=c(0,0), 
             id2=judas.cleaned.HR.KM[, as.numeric(unclass(as.factor(JUDAS_ID)))], 
             dev=as.matrix(judas.cleaned.HR.KM[, .(xdev, ydev)]))

inits <- function(){list(mu.sigmax=0,
                         mu.sigmay=0,mu.rho=0.5,sigx=1,sigy=1,
                         sig.rho=1)}

params<- c("mu.sigmax","mu.sigmay","a","sigx","sigy","b","rho","sigmax","sigmay")

# fit model to data using WinBUGS code
ni <- 2000
nb <- 1000
nthin <- 1
nc <- 4
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
analysis.path <- "../Data/Analysis"
save(fit2.KM, file = file.path(analysis.path, "fit2.KM.rda"))
load(file.path(analysis.path, "fit2.KM.rda"))

#------------------------------------------------------------------------------#

#### Fit to whole dataset ####
# dev are the deviates in the X and Y  dimensions; 
# N2 is the total number of observations (locations) for each judas; 
# M is the number of individuals; mean[] is given as data c(0,0)
# id2 is the index of individuals of the retained observation of judas

# Number of all retained observation of each judas
N2 <- nrow(judas.cleaned.HR)

n <- length(judas.cleaned.HR[, unique(JUDAS_ID)])
data <- list(N2=N2, M=n, mean=c(0,0), 
             id2=judas.cleaned.HR[, as.numeric(unclass(as.factor(JUDAS_ID)))], 
             dev=as.matrix(judas.cleaned.HR[, .(xdev, ydev)]))

inits <- function(){list(mu.sigmax=0,
                         mu.sigmay=0,mu.rho=0.5,sigx=1,sigy=1,
                         sig.rho=1)}

params<- c("mu.sigmax","mu.sigmay","a","sigx","sigy","b","rho","sigmax","sigmay")

# fit model to data using WinBUGS code
ni <- 3000
nb <- 1000
nthin <- 1
nc <- 4
np <- 8 # Number of CPUs
fit2 = jags(data, inits, params,  model.file="./Models/HRmodel3.txt", 
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
save(fit2, file = file.path(analysis.path, "fit2.rda"))
load(file.path(analysis.path, "fit2.rda"))
#------------------------------------------------------------------------------#

#### multidimensional integration ####
source("./MultiDimIntegration.funcitons.R")

#------------------------------------------------------------------------------#
#### speed test ####
source("./MultiDimInteg_speedTest.R")
#> 3.64 324.21   0.45  41.64   0.44

#### Integration for PB ####
analysis.path <- "../Data/Analysis"
load(file = file.path(analysis.path, "fit1.PB2.rda"))
load(file.path(analysis.path, "fit2.PB.rda"))

# Individual Judas
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
quantile(pdet.PB[, 1], probs = c(0.025, 0.5, 0.975))
#       2.5%         50%       97.5% 
# 0.006285564 0.019556476 0.058474979 

# Mean h
pdet.PB2.mu <- mapp_int(b1=fit1.PB2$mean$mu.b1, b2=fit1.PB2$mean$mu.b2, 
                 sigmax=exp(qnorm(0.5, fit2.PB$mean$mu.sigmax, fit2.PB$mean$sigx)),
                 sigmay=exp(qnorm(0.5, fit2.PB$mean$mu.sigmay, fit2.PB$mean$sigy)), 
                 rho=2 * (qbeta(0.5, fit2.PB$mean$a, fit2.PB$mean$b) - 0.5)) 
pdet.PB2.mu
# DProb     Abs_error
#1 0.02109411 1.048961e-310

# Compute variablility around mean unconditional prob of detection
LL.fit1 <- length(fit1.PB2$sims.list$mu.b1) # the MCMC was run for longer in fit1 than it is in fit2

ss <- sample(1:LL.fit1, size = length(fit2.PB$sims.list$mu.sigmax), replace = FALSE)
system.time(
pdet.PB2.mu.mcmc <- mapp_int(b1=fit1.PB2$sims.list$mu.b1[ss], 
                        b2=fit1.PB2$sims.list$mu.b2[ss], 
                 sigmax=exp(qnorm(0.5, fit2.PB$sims.list$mu.sigmax, fit2.PB$sims.list$sigx)),
                 sigmay=exp(qnorm(0.5, fit2.PB$sims.list$mu.sigmay, fit2.PB$sims.list$sigy)), 
                 rho=2 * (qbeta(0.5, fit2.PB$sims.list$a, fit2.PB$sims.list$b) - 0.5)) 
)
quantile(pdet.PB2.mu.mcmc[, 1], probs = c(0.025, 0.5, 0.975))
#      2.5%        50%      97.5% 
# 0.01095804 0.02171224 0.03542346 
      
#### Integration for KM ####
analysis.path <- "../Data/Analysis"
load(file = file.path(analysis.path, "fit1.KM2.rda"))
load(file.path(analysis.path, "fit2.KM.rda"))

# Individual Judas
b1<- fit1.KM2$mean$b1
b2<- fit1.KM2$mean$b2
sigmax<- fit2.KM$mean$sigmax 
sigmay<- fit2.KM$mean$sigmay
rho <- fit2.KM$mean$rho

system.time(
  pdet.KM <- mapp_int(b1, b2, sigmax, sigmay, rho)
)

ggplot(pdet.KM, aes(DProb)) + geom_density() +xlim(c(0, 0.1))
ggsave(filename = file.path(analysis.path, "DProb.density.KM.pdf"))
quantile(pdet.KM[, 1], probs = c(0.025, 0.975))
# 2.5%       97.5% 
# 0.008268773 0.062367318 

# Mean h
pdet.KM.mu <- mapp_int(b1=fit1.KM2$mean$mu.b1, b2=fit1.KM2$mean$mu.b2, 
                 sigmax=exp(qnorm(0.5, fit2.KM$mean$mu.sigmax, fit2.KM$mean$sigx)),
                 sigmay=exp(qnorm(0.5, fit2.KM$mean$mu.sigmay, fit2.KM$mean$sigy)), 
                 rho=2 * (qbeta(0.5, fit2.KM$mean$a, fit2.KM$mean$b) - 0.5)) 
pdet.KM.mu
# DProb     Abs_error
#1 0.0366868 1.270618e-310

# Compute variablility around mean unconditional prob of detection
LL.fit1 <- length(fit1.KM2$sims.list$mu.b1) # the MCMC was run for longer in fit1 than it is in fit2

ss <- sample(1:LL.fit1, size = length(fit2.KM$sims.list$mu.sigmax), replace = FALSE)
system.time(
pdet.KM.mu.mcmc <- mapp_int(b1=fit1.KM2$sims.list$mu.b1[ss], 
                        b2=fit1.KM2$sims.list$mu.b2[ss], 
                 sigmax=exp(qnorm(0.5, fit2.KM$sims.list$mu.sigmax, fit2.KM$sims.list$sigx)),
                 sigmay=exp(qnorm(0.5, fit2.KM$sims.list$mu.sigmay, fit2.KM$sims.list$sigy)), 
                 rho=2 * (qbeta(0.5, fit2.KM$sims.list$a, fit2.KM$sims.list$b) - 0.5)) 
)
quantile(pdet.KM.mu.mcmc[, 1], probs = c(0.025, 0.5, 0.975))
#      2.5%        50%      97.5% 
# 0.03128371 0.03681183 0.04241456  

#### Integration for whole dataset ####
analysis.path <- "../Data/Analysis"
load(file = file.path(analysis.path, "fit1.rda"))
load(file.path(analysis.path, "fit2.rda"))

b1<- fit1$mean$b1
b2<- fit1$mean$b2
sigmax<- fit2$mean$sigmax 
sigmay<- fit2$mean$sigmay
rho <- fit2$mean$rho

system.time(
  pdet.CA <- mapp_int(b1, b2, sigmax, sigmay, rho)
)

ggplot(pdet.CA, aes(DProb)) + geom_density() 
ggsave(filename = file.path(analysis.path, "DProb.density.CA.pdf"))
quantile(pdet.CA[, 1], probs = c(0.025, 0.5, 0.975))
#  2.5%         97.5% 
#  0.007890469 0.034008355 0.061669921  

pdet.CA.mu <- mapp_int(b1=fit1$mean$mu.b1, b2=fit1$mean$mu.b2, 
                         sigmax=exp(qnorm(0.5, fit2$mean$mu.sigmax, fit2$mean$sigx)),
                 sigmay=exp(qnorm(0.5, fit2$mean$mu.sigmay, fit2$mean$sigy)), 
                 rho=2 * (qbeta(0.5, fit2$mean$a, fit2$mean$b) - 0.5))
pdet.CA.mu
# DProb     Abs_error
# 0.03478861 1.863031e-310

# Compute variablility around mean unconditional prob of detection
LL.fit1 <- length(fit1$sims.list$mu.b1) # the MCMC was run for longer in fit1 than it is in fit2

ss <- sample(1:LL.fit1, size = length(fit2$sims.list$mu.sigmax), replace = FALSE)
system.time(
pdet.CA.mu.mcmc <- mapp_int(b1=fit1$sims.list$mu.b1[ss], 
                        b2=fit1$sims.list$mu.b2[ss], 
                 sigmax=exp(qnorm(0.5, fit2$sims.list$mu.sigmax, fit2$sims.list$sigx)),
                 sigmay=exp(qnorm(0.5, fit2$sims.list$mu.sigmay, fit2$sims.list$sigy)), 
                 rho=2 * (qbeta(0.5, fit2$sims.list$a, fit2$sims.list$b) - 0.5)) 
)
quantile(pdet.CA.mu.mcmc[, 1], probs = c(0.025, 0.5, 0.975))
#      2.5%        50%      97.5% 
# 0.02975004 0.03483036 0.04010203 