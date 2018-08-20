library(tidyverse)
library(lubridate)
library(jagsUI)

#---------------------------------------------------
# summarise harvest and effort totals for Donkey data

harvest_tots<- tot.harvest %>% group_by(TrackedShire, Year, Month) %>% summarise(Rem=sum(N_Ferals))
effort_tots<- runs_month %>% group_by(TrackedShire, Year, Month) %>% summarise(Effort=sum(Effort)) 

harvest<- inner_join(harvest_tots, effort_tots)

dat<- harvest %>% filter(TrackedShire=="PB" & Year < 2017)

cumdat<- dat %>% group_by(Year) %>% mutate(cumR=cumsum(Rem),cumE=cumsum(Effort))

# Plot cumulative catch and effort by year
win.graph(10,10)
cumdat %>% ggplot(aes(cumE, cumR)) +
  geom_line() +
  facet_wrap(~ Year, nrow=4)

#---------------------------------------------------
# try simulated data
source("Simulation Functions.r")

sim<- sim.remove(5000, 20, 5, log.phi= -5.5, lam=1.2, lam.sigma=0.05)
dat<- do.call('rbind',sim$removed)

cumdat<- dat %>% group_by(Year) %>% mutate(cumR=cumsum(Rem),cumE=cumsum(Effort))

win.graph(10,10)
cumdat %>% ggplot(aes(cumE, cumR)) +
        geom_line() +
        facet_wrap(~ Year, nrow=4)

#---------------------------------------------------
# model specification in WinBUGS
modelFilename = 'donkeyremoval.txt'
cat('
    data {
      for(i in 1:nyears) {
        for(j in os[i]:(os[i+1]-1)) {
        cumx[j]<- sum(y[os[i]:j]) - y[j]
        }
      }
    }
    
    model {
    # priors
    #roi ~ dunif(1, 1.3)
    roi<- 1.2
    alpha ~ dnorm(0, 0.1)
    
    N[1] ~ dpois(lambda[1])
    log(lambda[1]) <- phi
    phi ~ dnorm(0, 0.1)

    # Abundance model
    for (j in 2:nyears) {
      N[j] ~ dpois(lambda[j])
      lambda[j]<- (lambda[j-1] - sumx[j-1]) * roi
    }
    
    # Data model for removals
    for (i in 1:nyears) {
      for(j in os[i]:(os[i+1]-1)) {
        n[j] <- N[i] - cumx[j] 
        cloglog(p[j]) <- alpha + log(eff[j]) 
        y[j] ~ dbin(p[j],n[j])
      }
    }
}', fill=TRUE, file=modelFilename)

#---------------------------------------------------
# Assemble data

tmp<- table(dat$Year)
os<- c(1,cumsum(tmp)+1)  # Offset variable
nyears<- length(tmp)
nrem<- as.vector(dat$Rem)
eff<- as.vector(dat$Effort)
sumx<- sapply(sim$removed, function(x) sum(x$Rem))


data = list(y=nrem, eff=eff, nyears=nyears,os=os, sumx=sumx)

N.start<- tapply(dat$Rem, dat$Year, sum)

params = c("alpha","N","roi")

inits = function() {
  list(alpha= -5, phi=log(6000), N=N.start/0.2)
}

#---------------------------------------------------
# call to JAGS

ni<- 500000
nb<- 100000
nt<- 50
nc<- 3


fit = jags(data, inits, params, model.file=modelFilename,n.chains=nc, n.iter=ni, n.burnin=nb, 
           n.thin=nt, parallel = TRUE)

print(fit, digits=3)

# check fit

win.graph(8,8)
plot(sim$pop,fit$mean$N)
segments(sim$pop, fit$q2.5$N, sim$pop, fit$q97.5$N)
abline(a=0,b=1)


#############################################
#
# Try model using total removals per year
#
#############################################

modelFilename = 'donkeyremoval.txt'
cat('
    model {
    
    # priors
    N0 ~ dpois(lam)
    log(lam) <- phi 
    phi ~ dnorm(0, 0.1)
    
    alpha ~ dnorm(0, 0.1)
    roi ~ dunif(1, 1.3)
    #roi<- 1.2
    # Model
    pop[1] ~ dpois(lam)
    N[1]<- pop[1] - n[1]
    n[1] ~ dbin(p[1], pop[1])
    
    
    for(i in 2:K) {
    n[i] ~ dbin(p[i], pop[i]) 
    pop[i] ~ dpois(N[i-1]*roi)
    N[i]<- pop[i] - n[i]
    
    }
    
    for(j in 1:K){
    cloglog(p[j])<- alpha + log(eff[j])
    }
    
    }', fill=TRUE, file=modelFilename)

#---------------------------------------------------
# Simulate yearly removal data

sim<- sim.remove2(10000, 20, log.phi = -4, lam=1.2, lam.sigma=0.05)

nrem<- as.vector(sim$rem[-1])
eff<- as.vector(sim$eff[-1])
K<- length(nrem)

win.graph(8,8)
plot(cumsum(sim$eff),cumsum(sim$pop))

data = list(K=K, n=nrem, eff=eff)


# arguments for bugs()

params = c("alpha","roi","N0","N")

inits = function() {
  list(alpha= -5, phi=runif(1,0,15), pop=round(nrem/0.2))
}


# call to JAGS

library(jagsUI)

ni<- 200000
nb<- 100000
nt<- 10
nc<- 3


fit = jags(data, inits, params, model.file=modelFilename,n.chains=nc, n.iter=ni, n.burnin=nb, 
           n.thin=nt, parallel = TRUE)

print(fit, digits=3)

# check fit

win.graph(8,8)
plot(sim$pop,c(fit$mean$N0, fit$mean$N))
abline(a=0,b=1)


