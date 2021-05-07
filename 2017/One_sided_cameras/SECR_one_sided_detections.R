
# Code for the blog post at
#  https://mmeredith.net/blog/2017/SECR_unpaired_cameras.htm

# SECR with irregular block of habitat with parallel processing

# Several possible scenarios:
#
# A: single cameras, ca. 10 clusters of 8, 
#  analyse only data for flank with largest number of captures = toss out 1 flank
# B: single cameras, ca. 10 clusters of 8, 
#  analyse data for both flanks as separate sessions = two sessions
# C: paired cameras, ca. 10 clusters of 8 (so double the number of cameras) = many pairs
# D: paired cameras, ca. 10 clusters of 4 (so same number of cameras) = few pairs
# 

library(secr)
library(spatstat)
library(beeswarm) # For the beeswarm plot

## PRELIMINARIES
## =============
# Create the habitat patch polygon:
source("patch.R")

# A function to calculate the number of captures and the number of relocations,
#   ie, captures of an animal in a new location 
getCapsRelocs <- function(CH) {
  # function to do 1 session
  getCapsRelocs1 <- function(ch) {
    nCap <- nrow(ch)
    if(nCap > 0) {
      tmp <- apply(ch, c(1,3), sum)
      reLoc <- sum(rowSums(tmp > 0) > 1)  # number of relocations
    } else {
      reLoc <- NA
    }
    return(c(nCap=nCap, reLoc=reLoc))
  }
  
  if(is.list(CH)) {
    nSes <- length(CH)
    out <- matrix(NA_real_, nSes, 2)
    colnames(out) <- c("nCap", "reLoc")
    for(ses in 1:nSes)
      out[ses, ] <- getCapsRelocs1(CH[[ses]])
  } else {
    out <- getCapsRelocs1(CH)
  }
  return(out)
}

## WORK THROUGH ONE SIMULATION
## ===========================

# Check the habitat patch
head(patch)
MASS::eqscplot(patch, type='l') # Small "L", not "1" one
# Convert to 'owin' object
win <- owin(poly=patch[115:1, ])
plot(win)  # check
( area <- area.owin(win) / 1e4 )  # ha
area / 100  # km2

# Data generation
# ---------------
# Biological model
# ''''''''''''''''
# Generate activity centres
N <- 20
( D <- N / area )  # animals / ha
( DD <- round(D *100*100, 3) )  # animals / 100 sq km
AC <- rSSI(3000,  N, win)      # ACs at least 3km apart
points(AC, pch=19, col='red')  # add to the plot

## Observation model
## '''''''''''''''''
# Lots of trap clusters spread out
trapSpacing <- 2500
# for A, B, and C, each cluster 3 x 3 but hollow 
cluster <- make.grid(3, 3, spacing=trapSpacing, hollow = TRUE,
  detector='proximity')
traps <- make.systematic(c(5, 4), cluster, patch)
plot(traps, add=TRUE)
nrow(traps)  # number of traps
# For scenario D
clusterD <- make.grid(2, 2, spacing=trapSpacing, detector='proximity')
trapsD <- make.systematic(c(5, 4), clusterD, patch)
points(trapsD, col= 'blue', pch=4)

# Generate capture histories
g0 <- 0.02
sigma <- 2000
nOcc <- 90  # number of trapping occasions
popn <- data.frame(x=AC$x, y=AC$y)
class(popn) <- c("popn", "data.frame")
# Both Left and Right captured (scenario C):
simCH_C <- sim.capthist(traps, popn=popn,
    detectpar=list(g0=g0, sigma=sigma), noccasions=nOcc)
# Split randomly into 2
is.L <- array(rbinom(prod(dim(simCH_C)), 1, 0.5), dim=dim(simCH_C)) # random 1/0 array
simCHL <- simCH_C * is.L
simCHL <- subset(simCHL, dropnullCH = TRUE)  # Keep only flanks captured
simCHR <- simCH_C * (1 - is.L)
simCHR <- subset(simCHR, dropnullCH = TRUE) 
sum(simCHL, simCHR) == sum(simCH_C)  # check total number of capture events
( nCapL <- nrow(simCHL) )
( nCapR <- nrow(simCHR) )

# For scenario A, use whichever flank has most captures:
simCH_A <- if (nCapR > nCapL) simCHR else simCHL

# For scenario B, combine into a multisession capture history:
simCH_B <- MS.capthist(L = simCHL, R = simCHR)

# Scenario D, use the smaller trap array
simCH_D <- sim.capthist(trapsD, popn=popn,
    detectpar=list(g0=g0, sigma=sigma),
    noccasions=nOcc)
nrow(simCH_D)

# Data analysis
# -------------
# Create a fairly sparse mask to reduce run time.
mask <- make.mask(traps, buffer=25000, nx=32, ny=32, poly=patch)
plot(mask, add=TRUE)
# Starting values may be necessary for some simulated data sets
secrStart <- c(D=log(D), g0=qlogis(g0), sigma=log(sigma))

fitA <- secr.fit(simCH_A, mask=mask, start=secrStart)  # A = toss out one flank
predict(fitA)
predict(fitA)[1, -1]*100*100 # Convert to animals per 100km2
DD

fitB <- secr.fit(simCH_B, mask=mask, start=secrStart)  # B = two sessions
predict(fitB)
predict(fitB)[[1]][1, -1]*100*100 # Convert to animals per 100km2
DD

fitC <- secr.fit(simCH_C, mask=mask, start=secrStart)  # C = many pairs
predict(fitC)
predict(fitC)[1, -1]*100*100 # Convert to animals per 100km2
DD

fitD <- secr.fit(simCH_D, mask=mask, start=secrStart)  # D = few pairs
predict(fitD)
predict(fitD)[1, -1]*100*100 # Convert to animals per 100km2
DD


# DO HUNDREDS
# ===========
# Set up the parallel stuff:
library(foreach)
library(parallel)  # parallel is part of the R core distribution
library(doParallel)

( ncore <- detectCores() - 1 )	# I like to keep 1 thread free
( cl <- makeCluster(ncore) )

clusterSetRNGStream(cl = cl) # set up random number stream
registerDoParallel(cl)  # register as backend for 'foreach'

# The bits needed to run the loop:
# be sure to run patch and getCapsRelocs
win <- owin(poly=patch[115:1, ])
area <- area.owin(win) / 1e4  # ha
trapSpacing <- 2500
cluster8 <- make.grid(3, 3, spacing=trapSpacing, hollow = TRUE, detector='proximity')
cluster4 <- make.grid(2, 2, spacing=trapSpacing, detector='proximity')
N <- 20
D <- N / area  # animals / ha
g0 <- 0.02
sigma <- 2000
nOcc <- 90
trapsTmp <- make.systematic(c(5, 4), cluster8, patch) # only for mask
mask <- make.mask(trapsTmp, buffer=25000, nx=32, ny=32, poly=patch)
secrStart <- c(D=log(D), g0=qlogis(g0), sigma=log(sigma))

# seeds <- 1:301  # I was using 7 cores, so 43 per core = 301
seeds <- 1:6  # for testing
system.time(
result <- foreach(x=seeds, .combine=rbind, .packages=c("secr", "spatstat")) %dopar% {
  set.seed(x)

  # Trap layouts
  traps <- make.systematic(c(5, 4),	cluster8, patch) # scenarios A, B, C
  trapsD <- make.systematic(c(5, 4),	cluster4, patch) # scenario D
  outTraps <- c(nTrapsABC = nrow(traps), nTrapsD = nrow(trapsD))
  
  # Population, ie, activity centre locations (same for all scenarios)
  AC <- rSSI(3000,  N, win)      # ACs at least 3km apart
  popn <- data.frame(x=AC$x, y=AC$y)
  class(popn) <- c("popn", "data.frame")
  
  # List to hold capture histories (allows us to use a loop for the model fitting)
  CHlist <- vector("list", 4)
  names(CHlist) <- c("A", "B", "C", "D")
  
  # Capture histories for scenario C, "many pairs"
  CHlist$C <- simCH_C <- sim.capthist(traps, popn=popn,
      detectpar=list(g0=g0, sigma=sigma), noccasions=nOcc)
  # Split into 2 for use with scenarios A and B
  is.L <- array(rbinom(prod(dim(simCH_C)), 1, 0.5), dim=dim(simCH_C)) # random 1/0 array
  simCHL <- simCH_C * is.L
  simCHL <- subset(simCHL, dropnullCH = TRUE)  # Keep only flanks captured
  simCHR <- simCH_C * (1 - is.L)
  simCHR <- subset(simCHR, dropnullCH = TRUE) 
  nCapL <- nrow(simCHL)
  nCapR <- nrow(simCHR)
  CHlist$A <- if (nCapR > nCapL) simCHR else simCHL  # scenario A, "toss out"
  CHlist$B <- MS.capthist(simCHL, simCHR)            # scenario B, "2 sessions"
  
  # capture histories for scenario D, "few pairs"
  CHlist$D <- sim.capthist(trapsD, popn=popn,
    detectpar=list(g0=g0, sigma=sigma), noccasions=nOcc)

  # Analyse each of the capture histories
  outMat <- matrix(NA_real_, 6, 4)
  # rownames(outMat) <- c("nCaps", "nRelocs", "Dest", "Dse", "Dlow", "Dupp")
  # colnames(outMat) <- names(CHlist)  # names are stripped out before returning anyway
  for(ch in 1:4) {
    # Check the number of animals captured, and relocations:
    capRel <- getCapsRelocs(CHlist[[ch]])
    if(is.matrix(capRel))
      capRel <- colSums(capRel)
    outMat[1:2, ch] <- capRel
    if(capRel[1] > 0 && capRel[2] > 0) {  # proceed with analysis
      fit <- secr.fit(CHlist[[ch]], mask=mask, start=secrStart)
      pred <- predict(fit)
      if(!is.data.frame(pred))  # multisession fit gives a list of data frames
        pred <- pred[[1]]
      outMat[3:6,  ch] <- as.matrix(pred[-1])[1, ]*100*100  # extract the density estimate
    }
  }
  c(outTraps, as.vector(outMat))
} ) # Took 1 hr

# save(result, file="OneSide_sims_20170205.Rda")
dim(result)  # check
# Give the columns informative names
colnames(result) <- c("nTrapsABC", "nTrapsD", 
  t(outer(c("A", "B", "C", "D"), c("nCaps", "nRelocs", "Dest", "Dse", "Dlow", "Dupp"),
    paste, sep="_")))

colMeans(result, na.rm=TRUE)
( DD <- D*100*100)  # true density, animals per 100 sq km
colSums(is.na(result))  # How many failed?
mean(result[, 1]) # average numbers of traps
mean(result[, 2])
density <- result[, c(5, 11,17,23)] # per 100km2
colMeans(density, na.rm=TRUE)

# Calculate bias and RMSE
( bias <- colMeans(density, na.rm=TRUE) - DD )
getRmse <- function(x, true)
  sqrt(mean((x - true)^2, na.rm=TRUE))
( RMSE <- apply(density, 2, getRmse, true = DD) )

# Organise data for plotting
density0 <- density
density0[is.na(density0)] <- -0.1  # Replace NAs with small neg value for plotting
densityDF <- as.data.frame(density)
densityDF0 <- as.data.frame(density0)

beeswarm(densityDF0, method='hex', pch=19,	cex=0.5, las=1,
  pwcol=c('blue', 'red')[(density0 < 0) +1], ylim=c(-0.2, 5),
  labels=c("A:Toss out", "B:Two sessions", "C:Many pairs", "D:Few pairs"),
  main="Estimates of density", ylab="Density (animals per 100 sq km)")
boxplot(densityDF, yaxt='n', xaxt='n', pch=NA, add=TRUE)
abline(h=DD, col='red', lwd=2, lty=3)
abline(h=0)
text(1:4, 3.5, paste("bias:", round(bias, 2)), cex=0.8, pos=4)
text(1:4, 3.3, paste("RMSE:", round(RMSE, 2)), cex=0.8, pos=4)

# Check confidence interval coverage:
colnames(result)
mean(result[, 'A_Dlow'] < DD & result[, 'A_Dupp'] > DD, na.rm=TRUE)
mean(result[, 'B_Dlow'] < DD & result[, 'B_Dupp'] > DD, na.rm=TRUE)
mean(result[, 'C_Dlow'] < DD & result[, 'C_Dupp'] > DD, na.rm=TRUE)
mean(result[, 'D_Dlow'] < DD & result[, 'D_Dupp'] > DD, na.rm=TRUE)
