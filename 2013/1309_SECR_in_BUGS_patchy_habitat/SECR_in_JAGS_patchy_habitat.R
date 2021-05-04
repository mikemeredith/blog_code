#
# Analysis of SECR data for an irregular patch of habitat using JAGS
#
# For details see the blog post at
# https://mmeredith.net/blog/2013/1309_SECR_in_JAGS_patchy_habitat.htm

# With current packages, the simulated capture histories are  different, even
#  when using the same seed. Only animals in the west of the patch are captured.
# This is likely due to changes in the simulation functions in the 'secr' package
#  since 2013.

# The original code used 'rjags' functions, here I use jagsUI, so the MCMC results
#  are different.

library(jagsUI)
library(spatstat) # Functions to manipulate the polygon
library(secr)     # For simulation of trapping data
library(MASS)     # For equi-scale plot

source("patch.R")
eqscplot(patch, type='l')
# This is the patch of good habitat, surrounded by unsuitable terrain.
#   Units are metres. Total area = 2000 sq km.

# SIMULATE THE TRAPPING EXERCISE
# ==============================

# Put a block of traps in the middle of the polygon
# (this isn't the optimal layout but will do for this!)
traps <- make.grid(nx=16, ny=6, spacing=2500,
    originxy=c(135100, 783100), detector='proximity')
eqscplot(patch, type='l')
plot(traps, add=TRUE)

# Generate a (sparse) population of Activity Centres
set.seed(1)
popn <- sim.popn (D=3e-05, core=traps, buffer = 1e5, poly = patch)
nrow(popn)  # Number of animals in the polygon
points(popn)

# Simulate capture histories
simCH <- sim.capthist(traps, popn, detectfn = 0,
      detectpar = list(g0 = 0.02, sigma = 2500), noccasions = 90)
plot(simCH, add=TRUE, tracks=TRUE, rad=500, icolours=palette())
# Only animals in the west of the patch have been captured.
# In 2013, one of the animals in the SE was captured.

# CREATE THE HABITAT MATRIX
# =========================

# This is the original code doing it 'by hand',
#   using 'makeJAGSmask' is much much simpler!

# Find the SW corner of a box bounding the polygon:
( swCorner <- sapply(patch, min) )
# Sweep this from the data frame:
patch2 <- sweep(patch, 2, swCorner, "-")
sapply(patch2, min) == c(0, 0)  # just checking

# Find divisor needed to rescale to max dimension of (say) 50,
#  so maximum size of habitat matrix is 50x50.
( div <- max(patch2) / 50 )
patch3 <- patch2 / div
eqscplot(patch3, type='l')
sapply(patch3, max)
  # 50, 45.9, so habitat matrix will be 50x46

# Turn the polygon into an 'owin' object.
# The x, y matrix needs to go anticlockwise and
#   the first and last points must NOT be the same.
tmp <- nrow(patch3)
win <- owin(poly=patch3[tmp:2, ])
plot(win)

# Generate habitat data frame (one row per pixel)
( matdims <- round(sapply(patch3, max)) )
hab <- data.frame(x=rep(1:matdims[1], each=matdims[2]), y=1:matdims[2]) - 0.5
head(hab)         # Coordinates of centres of pixels
plot(win, lwd=3)
points(hab)
# See which points (pixel centres) are inside 'win'
#   This means that a pixel is classed as 'good habitat' if its centre
#   is in the patch.
hab$hab <- inside.owin(hab$x, hab$y, win)
plot(win, lwd=3)
points(hab$x, hab$y, col=c('grey', 'red')[hab$hab+1], cex=0.7)
  # Red if inside, grey if outside

# Turn the data frame into a matrix
# hab$hab is logical; the " * 1" coerces it to numeric
habmat <- t(matrix(hab$hab, nrow=matdims[2], ncol=matdims[1]) * 1)
image(z=habmat)


# DO THE JAGS ANALYSIS
# ====================

# Convert the trap locations to 'pixel' units to match habmat:
traps1 <- sweep(traps, 2, swCorner, "-") / div
plot(win)
points(traps1) # plot to check
trapmat <- as.matrix(traps1)  # Need to convert from data frame to matrix

# Add up number of detections across occasions
dim(simCH)
# simCH is an array, 4 animals x 90 occasions x 96 traps
yObs <- apply(simCH, c(1,3), sum)
dim(yObs)

# Data augmentation: add all-zero rows
M <- 50
J <- nrow(traps)
y <- rbind(yObs, matrix(0, M - nrow(yObs), J))

# organise data for JAGS
dataList <- list(
      y = y,
      M = M,
      trapmat = trapmat,
      J = J,
      xu = nrow(habmat), yu = ncol(habmat),
      habmat = habmat,
      OK = rep(1, M),  # all the ACs are known to be in good habitat
      A = 2000, # sq km
      K = 90  )

# initial values:
# Starting points for ACs in middle of the patch
inits <- function() list(z = rep(1, M),
              SX = rep(25, M), SY = rep(23, M) )

# estimates to extract:
wanted <- c("N", "D", "psi", "lam0", "sigma", "SX", "SY", "z")

# Run the model (takes 2 mins on my laptop)
( out <- jags(dataList, inits, wanted, "patch_model.jags", DIC=FALSE,
    n.chains=3, n.adapt=500, n.iter=2000,
    codaOnly = c("SX", "SY", "z"), parallel=TRUE) )
# Diagnostics are ok.

# result <- as.matrix(mcmcList)
# Check largest value of N << M
M
max(out$sims.list$N)  # 25, OK.

# In this simulation there were 8 animals in 2000 sq km, so density =
8/2000  # per sq km
# The estimate is very good

# Capture histories were simulated with g0 = 0.02 (g0 is roughly equivalent to
#   lam0 for small values) and sigma = 2500.
# Here we have underestimated sigma and overestimated lam0.

# CHECK ACTIVITY CENTRE LOCATIONS
# ===============================

( niter <- out$mcmc.info$n.samples )
# Too many, only need to plot 50 points
toplot <- round(seq.int(1, niter, length=50))

# Plot ACs of animals captured
( nCap <- nrow(simCH) ) # 4 animals captured in this simulation
palette("R3")
plot(win, main="")
for(i in 1:nCap)
  points(x=out$sims.list$SX[toplot, i],
         y=out$sims.list$SY[toplot, i], pch=19, cex=0.3, col=i)
points(trapmat, pch='+', col='grey')
title(sub="50 simulated ACs\nfor each animal captured")

# Plot ACs for all REAL animals not captured
# Replace coordinates for phantoms with NA
SXreal <- out$sims.list$SX
SXreal[out$sims.list$z == 0] <- NA
SYreal <- out$sims.list$SY
SYreal[out$sims.list$z == 0] <- NA
# Keep only uncaptured animals and 50 draws
nocapx <- SXreal[toplot, -(1:nCap)]
nocapy <- SYreal[toplot, -(1:nCap)]
plot(win, main="")
# Uncomment these lines if you want to display the habitat matrix pixels:
# for(x in 1:50)   for(y in 1:46)    if(habmat[x, y])
#       points(x-0.5, y-0.5, pch=22, col='grey', bg='grey', cex=1.3)
# plot(win, lwd=3, add=TRUE)
points(nocapx, nocapy, pch=19, cex=0.3, col='blue')
points(trapmat, pch='+', col='grey')
title(sub="Simulated ACs\nfor real animals not captured")
# ACs don't go outside the area defined by the habitat matrix, which is slightly
#   larger than the habitat polygon.
# And they do appear to be random (with reassuring 'faces-in-the-clouds').
