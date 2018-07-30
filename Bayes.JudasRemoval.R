library(data.table)
library(ggplot2)
library(jagsUI)

# Prep data and fit model
data.path <- "../Data"
dir.create(file.path(data.path, "Analysis"))
analysis.path <- "../Data/Analysis"
load(file.path(data.path, "runs.rda"))

# log removal rate from data
runs_years_long[, .(minlog=log(min(N_Donkeys/Effort)), 
                    meanlog=log(mean(N_Donkeys/Effort)),
                                maxlog=log(mean(N_Donkeys/Effort))), 
                    by=c("Region", "Method")]

# arguments for bugs()
params = c("alpha","roi","N0","N", "pop", "p", "fit", "fit.sim", "n.est", "n.sim")

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
check_fit.KM <- data.frame(runs_years_KM[, .(Year, Judas.culling, N_Ferals)], 
                           judas.est=n.est[, 1], opp.est=n.est[,2])

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

fit.KM$mean$fit/fit.KM$mean$fit.sim
mean(fit.KM$sims.list$fit.sim > fit.KM$sims.list$fit)

ggplot(data.frame(fit=fit.KM$sims.list$fit, fit.sim=fit.KM$sims.list$fit.sim),
       aes(fit, fit.sim)) + geom_point() + xlim(c(0, 15000)) + ylim(c(0, 15000)) +
  geom_abline(aes(slope=1, intercept=0), col="red")
#------------------------------------------------------------------------------#

#### Fit model to PB data ####
runs_years_PB <- runs_years[Region == "PILBARA",]
data <- list(nyears=nrow(runs_years_PB), n=runs_years_PB[, .(Judas.culling, N_Ferals)], 
             eff=runs_years_PB[, .(Effort_judas, Effort_opp)], nmethods=2)

params = c("alpha","roi","N0","N", "pop", "p", "fit", "fit.sim", "n.est", "n.sim")
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

fit.PB$mean$fit/fit.PB$mean$fit.sim
mean(fit.PB$sims.list$fit.sim > fit.PB$sims.list$fit)

ggplot(data.frame(fit=fit.PB$sims.list$fit, fit.sim=fit.PB$sims.list$fit.sim),
       aes(fit, fit.sim)) + geom_point() + xlim(c(0, 15000)) + ylim(c(0, 15000)) +
  geom_abline(aes(slope=1, intercept=0), col="red")

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

fit.PB.roi$mean$fit/fit.PB.roi$mean$fit.sim
mean(fit.PB.roi$sims.list$fit.sim > fit.PB.roi$sims.list$fit)

ggplot(data.frame(fit=fit.PB.roi$sims.list$fit, fit.sim=fit.PB.roi$sims.list$fit.sim),
       aes(fit, fit.sim)) + geom_point() + xlim(c(0, 15000)) + ylim(c(0, 15000)) +
  geom_abline(aes(slope=1, intercept=0), col="red")
#------------------------------------------------------------------------------#

#### Fit model to both pops' data ####
yrs_diff <- nrow(runs_years_KM) - nrow(runs_years_PB)
data <- list(nyears=c(nrow(runs_years_KM), nrow(runs_years_PB)), 
             n=array(c(runs_years_KM[, Judas.culling], runs_years_KM[, N_Ferals], 
                         runs_years_PB[, Judas.culling], rep(NA, yrs_diff),
                           runs_years_PB[, N_Ferals], rep(NA, yrs_diff)), c(24, 2, 2)), 
             eff=array(c(runs_years_KM[, Effort_judas], runs_years_KM[, Effort_opp],
                         runs_years_PB[, Effort_judas], rep(NA, yrs_diff),
                         runs_years_PB[, Effort_opp], rep(NA, yrs_diff)),
                       c(24, 2, 2)),
                         nmethods=2, npops=2)

inits <- function() {
  #list(a=rnorm(5,0,1), b=rnorm(5,0,1),u=10, prior.p= runif(1),log.eta=10)
  list(alpha=runif(2, -10, -4), phi=c(runif(1, 8, 12), runif(1, 6, 8)),
       pop=matrix(round(c(runs_years_KM[, Judas.culling + N_Ferals]/0.2, 
                          runs_years_PB[, Judas.culling + N_Ferals]/0.2, 
                          rep(NA, yrs_diff))), ncol=2))
}

# call to JAGS
ni<- 2000000
nb<- 1000000
nt<- 100
nc<- 4
np <- 8 # Number of CPUs

params = c("alpha","roi","N0", "fit", "fit.sim", "N", "pop", "n.est", "n.sim")

fit.both = jags(data, inits, params, model.file="./Models/JudasRm_multipleMethods_multPops.txt", 
              n.chains=nc, n.iter=ni, n.burnin=nb, 
              n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
              n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fit.both, digits=3)
sapply(fit.both$n.eff, summary)
sapply(fit.both$Rhat, summary)
plot(fit.both)

# Plot fit
ns <- rbind(fit.both$mean$n.est[, , 1], fit.both$mean$n.est[1:20, , 2])
ns <- as.data.frame(cbind(ns, c(runs_years_KM[, Judas.culling], runs_years_PB[, Judas.culling]),
                         c(runs_years_KM[, N_Ferals], runs_years_PB[, N_Ferals])))

ns <- cbind(Year=c(1994:2017, 1998:2017), Pop=c(rep("KM", 24), rep("PB", 20)), ns)
names(ns)[3:6] <- c("judas.est", "opp.est", "Judas.culling", "N_Ferals")

ggplot(ns) + geom_point(aes(Year, Judas.culling), col="blue", shape=16) + 
  geom_point(aes(Year, judas.est), col="red", shape=16) +
  geom_point(aes(Year, N_Ferals), col="blue", shape=7) +
  geom_point(aes(Year, opp.est), col="red", shape=7) +
  ylab("Number of animals removed") + facet_grid(Pop~.)

ggsave(filename=file.path(analysis.path, "Judas_Rm_Mod_fit_both.pdf"))

# Mean difference and sd  of judas and opportunistic
mean(ns$Judas.culling - ns$judas.est)
sd(ns$Judas.culling - ns$judas.est)

mean(ns$N_Ferals - ns$opp.est)
sd(ns$N_Ferals - ns$opp.est)

fit.both$mean$fit/fit.both$mean$fit.sim
mean(fit.both$sims.list$fit.sim > fit.both$sims.list$fit)

ggplot(data.frame(fit=fit.both$sims.list$fit, fit.sim=fit.both$sims.list$fit.sim),
       aes(fit, fit.sim)) + geom_point() + xlim(c(0, 15000)) + ylim(c(0, 15000)) +
  geom_abline(aes(slope=1, intercept=0), col="red")
#------------------------------------------------------------------------------#

# Fit model to both pops' data with two alphas
yrs_diff <- nrow(runs_years_KM) - nrow(runs_years_PB)
data <- list(nyears=c(nrow(runs_years_KM), nrow(runs_years_PB)), 
             n=array(c(runs_years_KM[, Judas.culling], runs_years_KM[, N_Ferals], 
                       runs_years_PB[, Judas.culling], rep(NA, yrs_diff),
                       runs_years_PB[, N_Ferals], rep(NA, yrs_diff)), c(24, 2, 2)), 
             eff=array(c(runs_years_KM[, Effort_judas], runs_years_KM[, Effort_opp],
                         runs_years_PB[, Effort_judas], rep(NA, yrs_diff),
                         runs_years_PB[, Effort_opp], rep(NA, yrs_diff)),
                       c(24, 2, 2)),
             nmethods=2, npops=2)

inits <- function() {
  #list(a=rnorm(5,0,1), b=rnorm(5,0,1),u=10, prior.p= runif(1),log.eta=10)
  list(alpha=matrix(runif(4, -10, -4), nrow=2), phi=c(runif(1, 8, 12), runif(1, 6, 8)),
       pop=matrix(round(c(runs_years_KM[, Judas.culling + N_Ferals]/0.2, 
                          runs_years_PB[, Judas.culling + N_Ferals]/0.2, 
                          rep(NA, yrs_diff))), ncol=2))
}

# call to JAGS
ni<- 2000000
nb<- 1000000
nt<- 100
nc<- 4
np <- 8 # Number of CPUs

fit.both.2a = jags(data, inits, params, model.file="./Models/JudasRm_multipleMethods_multPops_sepAlpha.txt", 
                n.chains=nc, n.iter=ni, n.burnin=nb, 
                n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fit.both.2a, digits=3)
sapply(fit.both.2a$n.eff, summary)
sapply(fit.both.2a$Rhat, summary)
plot(fit.both.2a)

# Plot fit
ns <- rbind(fit.both.2a$mean$n.est[, , 1], fit.both.2a$mean$n.est[1:20, , 2])
ns <- as.data.frame(cbind(ns, c(runs_years_KM[, Judas.culling], runs_years_PB[, Judas.culling]),
                          c(runs_years_KM[, N_Ferals], runs_years_PB[, N_Ferals])))

ns <- cbind(Year=c(1994:2017, 1998:2017), Pop=c(rep("KM", 24), rep("PB", 20)), ns)
names(ns)[3:6] <- c("judas.est", "opp.est", "Judas.culling", "N_Ferals")

ggplot(ns) + geom_point(aes(Year, Judas.culling), col="blue", shape=16) + 
  geom_point(aes(Year, judas.est), col="red", shape=16) +
  geom_point(aes(Year, N_Ferals), col="blue", shape=7) +
  geom_point(aes(Year, opp.est), col="red", shape=7) +
  ylab("Number of animals removed") + facet_grid(Pop~.)

ggsave(filename=file.path(analysis.path, "Judas_Rm_Mod_fit_both.2a.pdf"))

# Mean difference and sd  of judas and opportunistic
mean(ns$Judas.culling - ns$judas.est)
sd(ns$Judas.culling - ns$judas.est)

mean(ns$N_Ferals - ns$opp.est)
sd(ns$N_Ferals - ns$opp.est)

fit.both.2a$mean$fit/fit.both.2a$mean$fit.sim
mean(fit.both.2a$sims.list$fit.sim > fit.both.2a$sims.list$fit)

ggplot(data.frame(fit=fit.both.2a$sims.list$fit, fit.sim=fit.both.2a$sims.list$fit.sim),
       aes(fit, fit.sim)) + geom_point() + xlim(c(0, 10000)) + ylim(c(0, 10000)) +
  geom_abline(aes(slope=1, intercept=0), col="red")

ggsave(filename=file.path(analysis.path, "Judas_Rm_Mod_discrep_fit_both.2a.pdf"))
#------------------------------------------------------------------------------#

# Fit model to both pops' data with two alphas and fixed roi
yrs_diff <- nrow(runs_years_KM) - nrow(runs_years_PB)
data <- list(nyears=c(nrow(runs_years_KM), nrow(runs_years_PB)), 
             n=array(c(runs_years_KM[, Judas.culling], runs_years_KM[, N_Ferals], 
                       runs_years_PB[, Judas.culling], rep(NA, yrs_diff),
                       runs_years_PB[, N_Ferals], rep(NA, yrs_diff)), c(24, 2, 2)), 
             eff=array(c(runs_years_KM[, Effort_judas], runs_years_KM[, Effort_opp],
                         runs_years_PB[, Effort_judas], rep(NA, yrs_diff),
                         runs_years_PB[, Effort_opp], rep(NA, yrs_diff)),
                       c(24, 2, 2)),
             nmethods=2, npops=2)

inits <- function() {
  #list(a=rnorm(5,0,1), b=rnorm(5,0,1),u=10, prior.p= runif(1),log.eta=10)
  list(alpha=matrix(runif(4, -10, -4), nrow=2), phi=c(runif(1, 8, 12), runif(1, 6, 8)),
       pop=matrix(round(c(runs_years_KM[, Judas.culling + N_Ferals]/0.2, 
                          runs_years_PB[, Judas.culling + N_Ferals]/0.2, 
                          rep(NA, yrs_diff))), ncol=2))
}

# call to JAGS
ni<- 200000
nb<- 100000
nt<- 10
nc<- 4
np <- 8 # Number of CPUs

fit.both.2a.Froi = jags(data, inits, params, model.file="./Models/JudasRm_multMtds_multPops_sepAlpha_Fixed_roi.txt", 
                   n.chains=nc, n.iter=ni, n.burnin=nb, 
                   n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                   n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fit.both.2a.Froi, digits=3)
sapply(fit.both.2a.Froi$n.eff, summary)
sapply(fit.both.2a.Froi$Rhat, summary)
plot(fit.both.2a.Froi)

# Plot fit
ns <- rbind(fit.both.2a.Froi$mean$n.est[, , 1], fit.both.2a.Froi$mean$n.est[1:20, , 2])
ns <- as.data.frame(cbind(ns, c(runs_years_KM[, Judas.culling], runs_years_PB[, Judas.culling]),
                          c(runs_years_KM[, N_Ferals], runs_years_PB[, N_Ferals])))

ns <- cbind(Year=c(1994:2017, 1998:2017), Pop=c(rep("KM", 24), rep("PB", 20)), ns)
names(ns)[3:6] <- c("judas.est", "opp.est", "Judas.culling", "N_Ferals")

ggplot(ns) + geom_point(aes(Year, Judas.culling), col="blue", shape=16) + 
  geom_point(aes(Year, judas.est), col="red", shape=16) +
  geom_point(aes(Year, N_Ferals), col="blue", shape=7) +
  geom_point(aes(Year, opp.est), col="red", shape=7) +
  ylab("Number of animals removed") + facet_grid(Pop~.)

ggsave(filename=file.path(analysis.path, "Judas_Rm_Mod_fit_both.2a.Froi.pdf"))

# Mean difference and sd  of judas and opportunistic
mean(ns$Judas.culling - ns$judas.est)
sd(ns$Judas.culling - ns$judas.est)

mean(ns$N_Ferals - ns$opp.est)
sd(ns$N_Ferals - ns$opp.est)

fit.both.2a.Froi$mean$fit/fit.both.2a.Froi$mean$fit.sim
mean(fit.both.2a.Froi$sims.list$fit.sim > fit.both.2a.Froi$sims.list$fit)

ggplot(data.frame(fit=fit.both.2a.Froi$sims.list$fit, fit.sim=fit.both.2a.Froi$sims.list$fit.sim),
       aes(fit, fit.sim)) + geom_point() + xlim(c(0, 10000)) + ylim(c(0, 10000)) +
  geom_abline(aes(slope=1, intercept=0), col="red")

#------------------------------------------------------------------------------#
save(list= c("fit.KM", "fit.PB", "fit.PB.roi", "fit.both", "fit.both.2a", "fit.both.2a.Froi"), 
     file = file.path(analysis.path, "FittedMods.rda"))
#------------------------------------------------------------------------------#

# compare the model results
fit.both$DIC - fit.both.2a$DIC
fit.both.2a$DIC - fit.both.2a.Froi$DIC


extr <- function(mod_fit, param, j=1) {
  m <- mod_fit$mean[[param]][j]
  low <- mod_fit$q2.5[[param]][j]
  upp <- mod_fit$q97.5[[param]][j]
  return(c(low, m, upp))
}

extr.bypop <- function(mod_fit, param, pn=1, j=1) {
  m <- mod_fit$mean[[param]][pn,j]
  low <- mod_fit$q2.5[[param]][pn, j]
  upp <- mod_fit$q97.5[[param]][pn, j]
  return(c(low, m, upp))
}

lmods <- list(fit.KM, fit.PB, fit.PB.roi, fit.both)
Alphaj <- lapply(lmods, extr, "alpha")
Alphaopp <- lapply(lmods, extr, "alpha", j=2)
RoI <- lapply(lmods, extr, "roi")

lmods.sep <- list(fit.both.2a, fit.both.2a.Froi)
AlphajKM <- lapply(lmods.sep, extr.bypop, "alpha")
AlphajPM <- lapply(lmods.sep, extr.bypop, "alpha", pn=2)
AlphaoppKM <- lapply(lmods.sep, extr.bypop, "alpha", j=2)
AlphaoppPB <- lapply(lmods.sep, extr.bypop, "alpha", pn=2, j=2)
RoI.CAsep <- lapply(lmods.sep, extr, "roi")

DIC <- lapply(c(lmods, lmods.sep), "[[", "DIC")

Analysis <- c("KM", "PB", "PB.roi", "CA", "CA.sepKM", "CA.sepKM.roi", "CA.sepPB", "CA.sepPB.roi")
Params <- c("Lower.Alpha.Judas", "Mean.Alpha.Judas", "Upper.Alpha.Judas", 
            "Lower.Alpha.opp", "Mean.Alpha.opp", "Upper.Alpha.opp", 
            "Lower.RoI", "Mean.RoI", "Upper.RoI", "DIC")

comp_a_roi <- as.data.frame(rbind(
                           do.call("cbind", args=c(Alphaj, AlphajKM, AlphajPB)),
                           do.call("cbind", args=c(Alphaopp, AlphaoppKM, AlphaoppPB)),
                           do.call("cbind", args=c(RoI, RoI.CAsep, RoI.CAsep)),
                           do.call(c, args=c(DIC, DIC[5:6]))))

names(comp_a_roi) <- Analysis
comp_a_roi <- cbind(Params, comp_a_roi)
comp_a_roi[, -1] <- round(comp_a_roi[,-1], 2)
write.csv(comp_a_roi, 
          file=file.path(analysis.path, "Comparison_alphas_roi.csv"), row.names=F)

# Plots
lmods <- list(fit.PB, fit.PB.roi)
PB.N <- lapply(lmods, extr, "N", j=1:20)
PB.N <- lapply(PB.N, matrix, nrow=20)

lmods.both <- list(fit.both, fit.both.2a, fit.both.2a.Froi)
PB.N.both <- lapply(lmods.both, extr.bypop, "N", pn=1:20, j=2)
PB.N.both <- lapply(PB.N.both, matrix, nrow=20)

PB.Ns <- do.call(rbind, c(PB.N, PB.N.both))

PB.Ns <- data.frame(Year=rep(1998:2017, 5), 
                    Analysis=rep(c("PB", "PB.roi", "CA", "CA.sep", "CA.sep.roi"), each=20),
                    PB.Ns)
names(PB.Ns)[3:5] <- c("Lower", "Mean", "Upper")

ggplot(PB.Ns, aes(Year, log(Mean, 10), col=Analysis, shape=Analysis)) + 
  geom_point(alpha=0.5) + 
  geom_errorbar(aes(ymin=log(Lower, 10), ymax=log(Upper, 10)), alpha=0.5)

ggsave(filename=file.path(analysis.path, "Popsize.est_comp_PB.pdf"))


KM.N <- lapply(list(fit.KM), extr, "N", j=1:24)
KM.N <- lapply(KM.N, matrix, nrow=24)

KM.N.both <- lapply(lmods.both, extr.bypop, "N", pn=1:24, j=1)
KM.N.both <- lapply(KM.N.both, matrix, nrow=24)

KM.Ns <- do.call(rbind, c(KM.N, KM.N.both))

KM.Ns <- data.frame(Year=rep(1994:2017, 4), 
                    Analysis=rep(c("KM", "CA", "CA.sep", "CA.sep.roi"), each=24),
                    KM.Ns)
names(KM.Ns)[3:5] <- c("Lower", "Mean", "Upper")

ggplot(KM.Ns, aes(Year, log(Mean, 10), col=Analysis, shape=Analysis)) + 
  geom_point(alpha=0.5) + 
  geom_errorbar(aes(ymin=log(Lower, 10), ymax=log(Upper, 10)), alpha=0.5)

ggsave(filename=file.path(analysis.path, "Popsize.est_comp_KM.pdf"))


































