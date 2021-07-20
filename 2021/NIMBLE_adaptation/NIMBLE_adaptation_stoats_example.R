
# Code for the blog post at
#  https://mmeredith.net/blog/2021/NIMBLE_adaptation.htm

library(nimble)
library(mcmcOutput)

# Helper function for plotting
# calculate the y coordinate for adding text to a plot
ypos <- function(ht) {
  usr <- par("usr")
  usr[3] + ht*(usr[4] - usr[3])
}


# Get the data from the 'secr' package
# ------------------------------------
library(secr)
?stoatDNA
data(stoatDNA)
str(stoatCH)
plot(stoatCH, rad=50, tracks=TRUE)

### keep distances in metres for the first runs ###

# Extract the detector locations
traps <- traps(stoatCH)
str(traps)
# Convert to a matrix
Detlocs <- as.matrix(traps)
( nDetlocs <- nrow(Detlocs) )  # 94 detectors

# Decide on state space
# ---------------------
# allow a buffer around the detector array
buffer <- 1000  # metres
# outer edges of the state space
( xmin <- min(Detlocs[,1]) - buffer )
( xmax <- max(Detlocs[,1]) + buffer )
( ymin <- min(Detlocs[,2]) - buffer )
( ymax <- max(Detlocs[,2]) + buffer )

# Plot it:
plot(Detlocs, xlim=c(xmin, xmax), ylim=c(ymin, ymax),
    bty='n', pch=3, col='red')
rect(xmin, ymin, xmax, ymax, border='red', lwd=3)

# Reformat the detection data as a matrix
# ---------------------------------------
# Aggregate to an animals x detectors matrix with number of occasions detected
table(stoatCH)  # 0/1
dim(stoatCH)    # animals x occasions x detectors
y <- apply(stoatCH, c(1,3), sum)
dim(y)
table(y)
sum(y) == sum(stoatCH)

# All the detectors were active for all 7 days
nOcc <- 7

# Data augmentation
# -----------------
( nCaps <- nrow(y) )  # 20 animals detected
nAug <- 200           # number of rows to add
yAug <- rbind(y, matrix(0, nAug, nDetlocs))
( M <- nrow(yAug) )   # 220

w <- c(rep(1, nCaps), rep(NA, nAug)) # first 20 known to be present

# Bundle data and constants for NIMBLE
# ------------------------------------
ndata <- list(y = yAug, w = w)
str(ndata)

nconst <- list(M = M, nOcc = nOcc, sigmaMax = 3000, # sigmaMax in metres
    Detlocs = Detlocs, nDetlocs = nDetlocs,
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
str(nconst)

# To work in kilometres instead of metres, change 'nconst' to this:
# nconst <- list(M = M, nOcc = nOcc, sigmaMax = 3000/1000, # sigmaMax in km
    # Detlocs = Detlocs/1000, nDetlocs = nDetlocs,
    # xmin = xmin/1000, xmax = xmax/1000, ymin = ymin/1000, ymax = ymax/1000)

# Specify model in NIMBLE dialect of BUGS language
# ------------------------------------------------

SCR.code <- nimbleCode({
  # Likelihood
  for(i in 1:M){                   # Loop over all M individuals
    w[i] ~ dbern(omega)            # w = 1 if animal is real/present
    AC[i, 1] ~ dunif(xmin, xmax)   # x-coord of activity centre
    AC[i, 2] ~ dunif(ymin, ymax)   # y coord of activity centre

    for(j in 1:nDetlocs){          # Loop over all detectors
      d2[i,j] <- (AC[i,1] - Detlocs[j,1])^2 +
        (AC[i,2] - Detlocs[j,2])^2              # d2 = distance^2
      # p[i,j] <- p0 * exp(- alpha * d2[i,j])   # Detection prob
      p[i,j] <- p0 * exp(- d2[i,j]/(2*sigma^2)) # Detection prob
      y[i,j] ~ dbin(p[i,j] * w[i], nOcc)        # The observed data
    }
  }

  # Priors
  p0 ~ dbeta(1, 1)            # Baseline detection probability
  sigma ~ dunif(0, sigmaMax)  # Half-normal scale
  omega ~ dbeta(1, 1)         # Data augmentation parameter

  # Derived quantities
  N <- sum(w[])                       # Population size
})

# Initial values
# --------------
# Running with no initial values specified (usually) works fine.

# To demonstrate the problem with bad initial values use this:
inits <- function()list(sigma = runif(1, 2, 5))

# For the same thing in kilometers, use this
# inits <- function()list(sigma = runif(1, 2, 5)/1000)

# Parameters monitored
# --------------------
wanted <- c("p0", "sigma", "omega", "N", "AC", "w")

# Run with the nimbleMCMC wrapper
# ===============================

set.seed(99)  # to get the values shown in the blog
system.time(
out <- nimbleMCMC(SCR.code, data=ndata, constants=nconst,
    # inits=inits,   # optional, include to see effect of bad initial values
    monitors=wanted, niter=6000, nburnin=1500, nchains=3,
    samplesAsCodaMCMC=TRUE) )  # 10 mins

( mco <- mcmcOutput(out) )
diagPlot(mco, params=wanted[1:4])
plot(mco, wanted[1:4])
View(summary(mco))

# check that N << M
max(mco$N)
M

# Do it bit by bit
# ================

# set options to make adaptation values accessible
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)


# Build the model
# ---------------
# 1. the version with a sensible initial value for sigma
Rmodel <- nimbleModel(SCR.code, name="SCR",
   data=ndata, constants=nconst, inits=list(sigma=255))

# 2. the version with a crazy starting value
# Rmodel <- nimbleModel(SCR.code, name="SCR",
    # data=ndata, constants=nconst, inits=list(sigma=3))

# 3. when distances are in km, and starting value still crazy
# Rmodel <- nimbleModel(SCR.code, name="SCR",
    # data=ndata, constants=nconst, inits=list(sigma=3/1000))

# Set up MCMC configuration (default values)
# -------------------------
mcmcConf <- configureMCMC(Rmodel, print=TRUE)
# ===== Monitors =====
# thin = 1: AC, p0, sigma, omega
# ===== Samplers =====
# RW sampler (442)
# - AC[]  (440 elements)
# - p0
# - sigma
# conjugate sampler (1)
# - omega
# binary sampler (200)
# - w[]  (200 elements)

# Build and compile the MCMC (and compile the model too)
# --------------------------
Rmcmc <- buildMCMC(mcmcConf)
Cmcmc <- compileNimble(Rmodel, Rmcmc)

# Generate MCMC output
# --------------------
# Run 1 chain with 2400 draws, keep all 2400
set.seed(1)
Cmcmc$Rmcmc$run(niter=2400, nburnin=0)
outList <- as.list(Cmcmc$Rmcmc$mvSamples)
str(outList)
tail(outList$sigma)  # if chain is stuck, these will all be the same

# Look at the sampler tuning phase
# --------------------------------
# Which sampler is used for sigma?
mcmcConf$printSamplers("sigma")
# [442] RW sampler: sigma  # note the index number

# Get the sampler scale and acceptance rate for the target sampler
scales <- Cmcmc$Rmcmc$samplerFunctions[[442]]$getScaleHistory()
accept <- Cmcmc$Rmcmc$samplerFunctions[[442]]$getAcceptanceHistory()
length(scales)  # 12

windows(12, 3)
plot(outList$sigma[,1], type='l', ylab="sigma")
abline(v=0:12*200, col='red')
text(x=0:11*200 + 100, y=ypos(0.95), round(scales, 2), col='red')
text(x=0:11*200 + 100, y=ypos(0.85), round(accept, 2), col='blue')

