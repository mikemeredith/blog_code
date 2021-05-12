# R code for the blog post at
#  https://mmeredith.net/blog/2020/GOF_2.htm


library(jagsUI)
library(HDInterval)
library(makeJAGSmask)
library(mcmcOutput)

# Load the data
# =============

load("SCRdata.RData")
ls()
# [1] "JAGSmask" "nOcc"     "sex"      "yObs" 

# Data augmentation
dim(yObs)  # 22 animals detected x 44 traps
nInd <- nrow(yObs)
nTraps <- ncol(yObs)
M <- 80  # Number of individuals after augmentation
Yaug <- matrix(0, M, nTraps)
Yaug[1:nInd, ] <- yObs
sexAug <- c(sex, rep(NA, M - length(sex)))

# Put together data for JAGS
JAGSdata0 <- list(Y = Yaug,
                  M = M,
                  nTraps=nTraps,
                  sex=sexAug - 1, # Convert factor to 0/1 vector
                  nOcc=nOcc,
                  OK=rep(1, M))
JAGSdata <- c(JAGSdata0, JAGSmask)
str(JAGSdata)  # check

# Starting values for AC locations:
#   for animals caught, use centroid of capture locations
ACcapt <- matrix(NA, nInd, 2)
for(i in 1:nInd) {
  captTraps <- which(yObs[i, ] > 0) # Which traps caught the animal
  captLocs <- JAGSmask$trapMat[captTraps, , drop=FALSE] # Locations of the traps
  ACcapt[i, ] <- colMeans(captLocs)
  # Check it's in good habitat (might not be if trap is on edge):
  stopifnot(JAGSmask$habMat[ACcapt[i, , drop=FALSE]] == 1)
}
head(ACcapt, 7) # Check
inits <- function() list(w = rep(1, M), AC = randomPoints(M, JAGSmask, ACcapt))
str(inits())  # check

# Estimates to extract:
wanted <- c("N", "p0pc", "sigmakm", "pi", "D", "omega", "mean.p0", "sig.p0", "lp0.sex",
            "T1obs", "T1sim", "T2obs", "T2sim", "T3obs", "T3sim",
            # "AC", "w", "sex", "Ynew")
           "Ysim")

# p0(sex) sigma(sex) model
# ========================

out1 <- jags(JAGSdata, inits, wanted, "SCRsex.jags", DIC=FALSE,
    n.chains=3, n.adapt=1000, n.iter=10000, n.thin=10, parallel=TRUE)
diagPlot(out1)

SL1 <- out1$sims.list

# Ratio of sigmas
ratio <- SL1$sigmakm[,2]/SL1$sigmakm[,1]
mean(ratio)
hdi(ratio)
mean(ratio > 1)  # >99% certain

# Ratio of p0's
ratio <- SL1$p0pc[,2]/SL1$p0pc[,1]
mean(ratio)
hdi(ratio)
mean(ratio > 1)  # >99% certain

# Check data augmentation
M
max(SL1$N[,1] + SL1$N[,2])  # << M


# Look at Bayesian p-value results
# --------------------------------

pp.check(out1, "T1obs", "T1sim", main="T1: animals x detectors")
pp.check(out1, "T2obs", "T2sim", main="T2: animals")
pp.check(out1, "T3obs", "T3sim", main="T3: detectors")

# Compare data with posterior predictive (simulated) data sets
# ------------------------------------------------------------
# Load the post. pred. simulations
Ysim <- out1$sims.list$Ysim
str(Ysim)

# Note that these capture histories have many all-zero rows,
#  keep in mind when calculating summary stats, esp. means

# jpeg(filename = "p0(sex)sigma(sex)hists_%03d.jpg")

# 1. total detections
( obs <- sum(yObs) )
sims <- apply(Ysim, 1, sum)
hist(sims, xlim=range(sims, obs), main="Total detections")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# 2. animals detected
getAnimals <- function(x) sum(rowSums(x) > 0)
( obs <- getAnimals(yObs) )
sims <- apply(Ysim, 1, getAnimals)
hist(sims, xlim=range(sims, obs), main="Animals detected")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# 3. total detectors visited
getTrapsVisited <- function(x) {
  visitsPerTrap <- colSums(x)
  sum(visitsPerTrap > 0) }
( obs <- getTrapsVisited(yObs) )
sims <- apply(Ysim, 1, getTrapsVisited)
hist(sims, xlim=range(sims, obs), main="Detectors visited")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# 4. max visits per detector
getMaxVisitsperTrap <- function(x) {
  visitsPerTrap <- colSums(x)
  max(visitsPerTrap) }
( obs <- getMaxVisitsperTrap(yObs) )
sims <- apply(Ysim, 1, getMaxVisitsperTrap)
hist(sims, xlim=range(sims, obs), main="Max visits per detector")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# 5. max detections for 1 animal at 1 detector
( obs <- max(yObs) )
sims <- apply(Ysim, 1, max)
hist(sims, xlim=range(sims, obs), main="Max detections for 1 animal at 1 detector")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# 6. max detections for 1 animal overall
getMaxCaps <- function(x) {
  capsPerAnimal <- rowSums(x)
  max(capsPerAnimal)}
( obs <- getMaxCaps(yObs) )
sims <- apply(Ysim, 1, getMaxCaps)
hist(sims, xlim=range(sims, obs), main="Max detections per animal")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# dev.off()

# Model with detector heterogeneity, p0(sex+det) sigma(sex) model
# ===============================================================

out2 <- jags(JAGSdata, inits, wanted, "SCRsex_trapHetero.jags", DIC=FALSE,
    n.chains=3, n.adapt=1000, n.iter=10000, n.thin=10, parallel=TRUE)
diagPlot(out2)

SL2 <- out2$sims.list

# Ratio of sigmas
ratio <- SL2$sigmakm[,2]/SL2$sigmakm[,1]
mean(ratio)
hdi(ratio)
mean(ratio > 1)  # >99% certain

# Check data augmentation
M
max(SL2$N[,1] + SL2$N[,2])  # << M


# Look at Bayesian p-value results
# --------------------------------

pp.check(out2, "T1obs", "T1sim", main="T1: animals x detectors")
pp.check(out2, "T2obs", "T2sim", main="T2: animals")
pp.check(out2, "T3obs", "T3sim", main="T3: detectors")

# Compare data with posterior predictive (simulated) data sets
# ------------------------------------------------------------
# Load the post. pred. simulations
Ysim <- out2$sims.list$Ysim
str(Ysim)

# jpeg(filename = "p0(sex)sigma(sex)hists_%03d.jpg")

# 1. total detections
( obs <- sum(yObs) )
sims <- apply(Ysim, 1, sum)
hist(sims, xlim=range(sims, obs), main="Total detections")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# 2. animals detected
getAnimals <- function(x) sum(rowSums(x) > 0)
( obs <- getAnimals(yObs) )
sims <- apply(Ysim, 1, getAnimals)
hist(sims, xlim=range(sims, obs), main="Animals detected")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# 3. total detectors visited
getTrapsVisited <- function(x) {
  visitsPerTrap <- colSums(x)
  sum(visitsPerTrap > 0) }
( obs <- getTrapsVisited(yObs) )
sims <- apply(Ysim, 1, getTrapsVisited)
hist(sims, xlim=range(sims, obs), main="Detectors visited")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# 4. max visits per detector
getMaxVisitsperTrap <- function(x) {
  visitsPerTrap <- colSums(x)
  max(visitsPerTrap) }
( obs <- getMaxVisitsperTrap(yObs) )
sims <- apply(Ysim, 1, getMaxVisitsperTrap)
hist(sims, xlim=range(sims, obs), main="Max visits per detector")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# 5. max detections for 1 animal at 1 detector
( obs <- max(yObs) )
sims <- apply(Ysim, 1, max)
hist(sims, xlim=range(sims, obs), main="Max detections for 1 animal at 1 detector")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# 6. max detections for 1 animal overall
getMaxCaps <- function(x) {
  capsPerAnimal <- rowSums(x)
  max(capsPerAnimal)}
( obs <- getMaxCaps(yObs) )
sims <- apply(Ysim, 1, getMaxCaps)
hist(sims, xlim=range(sims, obs), main="Max detections per animal")
abline(v=obs, col='red', lwd=3)
legend('topright', legend=paste("p =", round(mean(sims > obs),2)), bty='n')

# dev.off()
