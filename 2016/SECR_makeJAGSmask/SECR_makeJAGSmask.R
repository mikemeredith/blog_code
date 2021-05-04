
# Code for the blog post at
#  https://mmeredith.net/blog/2016/SECR_patchy_habitat_makeJAGSmask.htm

library(secr)
library(jagsUI)
library(makeJAGSmask)

remotes::install_github("mikemeredith/makeJAGSmask")
packageVersion("makeJAGSmask") # Should be 0.1.0 or later.

# You may need this:
# install.packages(c("abind", "plotrix", "raster", "secr", "rgeos"))

# Using secr::make.mask
# ---------------------

# Get the example data and plot them:
# Get the habitat polygon and traps data and plot them:
data(simSCR)
str(simSCR, 1)
attach(simSCR)
MASS::eqscplot(patchDF, type='l')
points(traps, pch=3, col='red')

# Generate a 'secr' mask and plot:
secrmask <- secr::make.mask(traps, spacing=1000, type='polygon', poly=patchDF)
str(secrmask)
plot(secrmask)
plot(traps, add=TRUE)

# Convert to a 'JAGSmask':
mymask <- convertMask(secrmask, traps)
str(mymask)

# Generate a (sparse) population of Activity Centres
RNGversion("3.3") # use the old (buggy) RNG
set.seed(2)
popn <- sim.popn (D=3e-05, core=traps, buffer = 1e5, poly = patchDF)
nrow(popn) # True number of animals in the polygon
MASS::eqscplot(patchDF, type='l')
plot(traps, add=TRUE)
points(popn, pch=17)

# Simulate capture histories
simCH <- sim.capthist(traps, popn, detectfn = 0,
  detectpar = list(g0 = 0.02, sigma = 2500), noccasions = 90)
plot(secrmask)
plot(simCH, add=TRUE, tracks=TRUE, rad=500, icolours=palette())
plot(traps, add=TRUE)
points(popn, pch=17)

# JAGS code for the model vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
JAGScode <- "
model {
  sigma ~ dunif(0, 10)     # set up the priors for the parameters
  alpha <- 1/(2*sigma^2)
  p0 ~ dbeta(1, 1)
  omega ~ dbeta(0.001, 1)

  for (i in 1:M){             # loop through the augmented population
    w[i] ~ dbern(omega)       # state of individual i (real or imaginary)
    S[i, 1] ~ dunif(1, upperLimit[1]) # priors for the activity centers for each individual
    S[i, 2] ~ dunif(1, upperLimit[2]) # lower x and y coordinates = 1
    pOK[i] <- habMat[trunc(S[i,1]), trunc(S[i,2])] # habitat check
    OK[i] ~ dbern(pOK[i])     # OK[i] = 1, the ones trick
    for(j in 1:nTraps) {      # loop through the camera trap locations
      Dsq[i,j] <- (S[i,1]-trapMat[j,1])^2 + (S[i,2]-trapMat[j,2])^2
      p[i,j] <- p0 * exp(-Dsq[i,j] * alpha)
      y[i,j] ~ dbin(p[i,j] * w[i], nOcc)
    }
  }
  N <- sum(w)   # derive number (check that N < M)
}"
writeLines(JAGScode, "patch.jags")
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Organise capture histories:
dim(simCH)
# simCH is an array, 5 animals x 90 occasions x 89 traps
# Add up number of detections across occasions
yObs <- apply(simCH, c(1,3), sum)
dim(yObs)

# Data augmentation: add all-zero rows
M <- 25
nTraps <- nrow(traps)
y <- rbind(yObs, matrix(0, M - nrow(yObs), nTraps))

# organise data for JAGS
dataList <- c(mymask, list(y = y, M = M, nTraps = nTraps, OK = rep(1, M), nOcc = 90))
str(dataList)

# Initial values - several options:
# (1) Quick and dirty, but works as {15, 15} is in good habitat:
inits <- function() list(w = rep(1, M), S = matrix(15, M, 2))

# (2) A sounder approach is to start all ACs at random locations in good habitat:
inits <- function() list(w = rep(1, M), S = randomPoints(M, mymask))
# ... but all-random AC locations will produce errors if the starting point for a captured
# animal is too far from the place it was captured.

# (3) Use starting values close to traps where the animal was caught.
Scapt <- matrix(NA, nrow(yObs), 2)
for(i in 1:nrow(yObs)) {
  captTraps <- which(yObs[i, ] > 0) # Which traps caught the animal
  captLocs <- mymask$trapMat[captTraps, , drop=FALSE] # Locations of the traps
  Scapt[i, ] <- colMeans(captLocs)
  # Check it's in good habitat (might not be if trap is on edge):
  stopifnot(mymask$habMat[Scapt[i, , drop=FALSE]] == 1)
}
head(Scapt, 7) # Check
inits <- function() list(w = rep(1, M), S = randomPoints(M, mymask, Scapt))
head(inits()$S, 7)
head(inits()$S, 7) # First 5 unchanged, last 2 different.

# Estimates to extract:
wanted <- c("N", "omega", "p0", "sigma", "S", "w")

# Run it, takes about 5 mins
result <- jags(dataList, inits, wanted, "patch.jags", DIC=FALSE,
  n.chains=3, n.adapt=1000, n.iter=5000, n.burnin=500, n.thin=10,
  parallel=TRUE)
result # n.eff very small, but ok for this test.
mcmcOutput::diagPlot(result)
# Check that data augmentation was adequate
max(result$sims.list$N) # Must be less than M
M

# Extract the chains for X and Y coordinates and w for 10 animals,
#   first 5 captured, next 5 not captured or phantoms.
S <- result$sims.list$S[, 1:10, ]
w <- result$sims.list$w[, 1:10]
# Remove phantom animals
S[w == 0] <- NA # Yes, this works fine, even though S has an extra dim:
dim(S)
dim(w)

# Convert to the original coordinate system
AC <- convertOutput(S, JAGSmask = mymask)
str(AC)

# Do plots with the original coordinates
MASS::eqscplot(patchDF, type='l')
for(i in 1:5)
  points(AC[, i, ], col=i)
points(AC[, 6:10, 1], AC[, 6:10, 2], col='grey')
points(traps, pch=3, col='red')
points(popn, pch=21, bg='white', cex=2) # true AC locations.

detach(simSCR) # clean up
