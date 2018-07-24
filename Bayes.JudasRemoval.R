library(jagsUI)

# Prep data and fit model
data.path <- "../Data"
analysis.path <- "../Data/Analysis"
load(file.path(data.path, "runs.rda"))

# log removal rate from data
runs_years_long[, .(minlog=log(min(N_Donkeys/Effort)), 
                    meanlog=log(mean(N_Donkeys/Effort)),
                                maxlog=log(mean(N_Donkeys/Effort))), 
                    by=c("Region", "Method")]

# arguments for bugs()
params = c("alpha","roi","N0","N", "pop", "p")


#### Fit model to KM data ####
runs_years_KM <- runs_years[Region == "KIMBERLEY",]
data <- list(nyears=nrow(runs_years_KM), n=runs_years_KM[, .(Judas.culling, N_Ferals)], 
             eff=runs_years_KM[, .(Effort_judas, Effort_opp)], nmethods=2)

inits <- function() {
  #list(a=rnorm(5,0,1), b=rnorm(5,0,1),u=10, prior.p= runif(1),log.eta=10)
  list(alpha=c(-5, -5), phi=runif(1, 8, 12), 
       pop=round(runs_years_KM[, Judas.culling + N_Ferals]/0.2))
}

# call to JAGS
ni<- 2000000
nb<- 1000000
nt<- 100
nc<- 4
np <- 8 # Number of CPUs

fit.KM = jags(data, inits, params, model.file="./Models/JudasRm_multipleMethods.txt", 
           n.chains=nc, n.iter=ni, n.burnin=nb, 
           n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
           n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fit.KM, digits=3)
sapply(fit.KM$n.eff, summary)
sapply(fit.KM$Rhat, summary)
plot(fit.KM)

# Plot fit
n.est <- fit.KM$mean$p * cbind(fit.KM$mean$pop, fit.KM$mean$pop)
check_fit.KM <- data.frame(runs_years_KM[, .(Year, Judas.culling, N_Ferals)], judas.est=n.est[, 1], opp.est=n.est[,2])

ggplot(check_fit.KM) + geom_point(aes(Year, Judas.culling), col="blue", shape=16) + 
  geom_point(aes(Year, judas.est), col="red", shape=16) +
  geom_point(aes(Year, N_Ferals), col="blue", shape=7) +
  geom_point(aes(Year, opp.est), col="red", shape=7) +
  ylab("Number of animals removed")

ggsave(filename=file.path(analysis.path, "Judas_Rm_Mod_fit_KM.pdf"))
  
# Mean difference and sd  of judas and opportunistic
mean(check_fit.KM$Judas.culling - check_fit.KM$judas.est)
sd(check_fit.KM$Judas.culling - check_fit.KM$judas.est)

mean(check_fit.KM$N_Ferals - check_fit.KM$opp.est)
sd(check_fit.KM$N_Ferals - check_fit.KM$opp.est)

print(fit.KM, digits=3)

#### Fit model to PB data ####
runs_years_PB <- runs_years[Region == "PILBARA",]
data <- list(nyears=nrow(runs_years_PB), n=runs_years_PB[, .(Judas.culling, N_Ferals)], 
             eff=runs_years_PB[, .(Effort_judas, Effort_opp)], nmethods=2)

inits <- function() {
  #list(a=rnorm(5,0,1), b=rnorm(5,0,1),u=10, prior.p= runif(1),log.eta=10)
  list(alpha=c(-5, -5), phi=runif(1, 6, 8), 
       pop=round(runs_years_PB[, Judas.culling + N_Ferals]/0.2))
}

# call to JAGS
ni<- 200000
nb<- 100000
nt<- 10
nc<- 4
np <- 8 # Number of CPUs

fit.PB = jags(data, inits, params, model.file="./Models/JudasRm_multipleMethods.txt", 
           n.chains=nc, n.iter=ni, n.burnin=nb, 
           n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
           n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fit.PB, digits=3)
sapply(fit.PB$n.eff, summary)
sapply(fit.PB$Rhat, summary)

# Plot fit
n.est <- fit.PB$mean$p * cbind(fit.PB$mean$pop, fit.PB$mean$pop)
check_fit.PB <- data.frame(runs_years_PB[, .(Year, Judas.culling, N_Ferals)], 
                           judas.est=n.est[, 1], opp.est=n.est[,2])

ggplot(check_fit.PB) + geom_point(aes(Year, Judas.culling), col="blue", shape=16) + 
  geom_point(aes(Year, judas.est), col="red", shape=16) +
  geom_point(aes(Year, N_Ferals), col="blue", shape=7) +
  geom_point(aes(Year, opp.est), col="red", shape=7) +
  ylab("Number of animals removed")

ggsave(filename=file.path(analysis.path, "Judas_Rm_Mod_fit_PB.pdf"))

# Mean difference and sd  of judas and opportunistic
mean(check_fit.PB$Judas.culling - check_fit.PB$judas.est)
sd(check_fit.PB$Judas.culling - check_fit.PB$judas.est)

mean(check_fit.PB$N_Ferals - check_fit.PB$opp.est)
sd(check_fit.PB$N_Ferals - check_fit.PB$opp.est)

#------------------------------------------------------------------------------#

# Fit with fixed roi
fit.PB.roi = jags(data, inits, params, model.file="./Models/JudasRm_multipleMethods_Fixed_roi.txt", 
              n.chains=nc, n.iter=ni, n.burnin=nb, 
              n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
              n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fit.PB.roi, digits=3)
sapply(fit.PB.roi$n.eff, summary)
sapply(fit.PB.roi$Rhat, summary)

# Plot fit
n.est <- fit.PB.roi$mean$p * cbind(fit.PB.roi$mean$pop, fit.PB.roi$mean$pop)
check_fit.PB.roi <- data.frame(runs_years_PB[, .(Year, Judas.culling, N_Ferals)], 
                           judas.est=n.est[, 1], opp.est=n.est[,2])

ggplot(check_fit.PB.roi) + geom_point(aes(Year, Judas.culling), col="blue", shape=16) + 
  geom_point(aes(Year, judas.est), col="red", shape=16) +
  geom_point(aes(Year, N_Ferals), col="blue", shape=7) +
  geom_point(aes(Year, opp.est), col="red", shape=7) +
  ylab("Number of animals removed")

ggsave(filename=file.path(analysis.path, "Judas_Rm_Mod_fit_PB_roi.pdf"))

# Mean difference and sd  of judas and opportunistic
mean(check_fit.PB.roi$Judas.culling - check_fit.PB.roi$judas.est)
sd(check_fit.PB.roi$Judas.culling - check_fit.PB.roi$judas.est)

mean(check_fit.PB.roi$N_Ferals - check_fit.PB.roi$opp.est)
sd(check_fit.PB.roi$N_Ferals - check_fit.PB.roi$opp.est)



