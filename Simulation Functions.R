
#########################################################################
#
# Removal Model Simulation functions
#
#######################################################################
sim.remove<- function(N.init, nyears, mpy, lam=1.2, lam.sigma=0.05, log.phi= -5) {
  # Simulated multi-year removal data with
  # multiple removal periods per year (mpy)
  # Population subject to annual changes following removal
  # lam - mean rate of increase (following removal)
  # log.phi - removal efficiency parameter (log scale)
  #
  removed<- list()
  pop<- rep(NA, nyears)
  ROI<- rep(NA, nyears-1)
  pop[1]<- N.init
  
  for(i in 1:(nyears-1)) {
    remyear<- c(0, rep(NA, mpy))
    popyear<- c(pop[i],rep(NA, mpy))
    effyear<- c(0,rnorm(mpy, 25, 2)) # 0 effort initially
    
    for(j in 2:(mpy+1)) {
      p<- 1-exp(-exp(log.phi + log(effyear[j])))  # Poisson model
      remyear[j]<- rbinom(1, popyear[j-1], p)
      popyear[j]<- popyear[j-1] - remyear[j]
    }
    lam.hat<- rnorm(1, lam, lam.sigma) # no decrease
    removed[[i]]<- data.frame(Year=i, N=popyear[-1], Rem=remyear[-1], Effort=effyear[-1])
    pop[i+1]<- round(popyear[mpy+1] * lam.hat)
    ROI[i]<- round(lam.hat, 3)
  }
  #last year
  remyear<- c(0, rep(NA, mpy))
  popyear<- c(pop[nyears],rep(NA, mpy))
  effyear<- c(0,rnorm(mpy, 25, 2)) # 0 effort initially
  
  for(j in 2:(mpy+1)) {
    p<- 1-exp(-exp(log.phi + log(effyear[j])))
    remyear[j]<- rbinom(1, popyear[j-1], p)
    popyear[j]<- popyear[j-1] - remyear[j]
  }
  removed[[nyears]]<- data.frame(Year=nyears, N=popyear[-1], Rem=remyear[-1], Effort=effyear[-1])
  
  list(pop=pop,ROI=ROI,removed=removed)
}

#---------------------------------------------------------------------

sim.remove2<- function(N.init, nyears, lam=1.2, lam.sigma=0.05, log.phi=-5) {
  # Simulated multi-year removal data with
  # single removal per year
  # Population subject to annual changes following removal
  # roi - instantaneous rate of increase (following removal)
  # log.phi - removal efficiency parameter (log scale)
  #
  pop<- c(N.init, rep(NA, nyears))
  rem<- c(0, rep(NA, nyears))
  eff<- c(0, rep(NA, nyears))
  ROI<- c(1, rep(NA, nyears))
  n<- nyears+1
  for(i in 2:n) {
    eff[i]<- rnorm(1, 25, 2)
    p<- 1-exp(-exp(log.phi + log(eff[i])))  # Poisson model
    rem[i]<- rbinom(1, pop[i-1], p)
    lam.hat<- rnorm(1, lam, lam.sigma)
    pop[i]<- round((pop[i-1] - rem[i]) * lam.hat)
    ROI[i]<- round(lam.hat, 3)
  }
  data.frame(year=0:nyears, pop=pop, rem=rem, eff=eff, ROI=ROI)
}  

