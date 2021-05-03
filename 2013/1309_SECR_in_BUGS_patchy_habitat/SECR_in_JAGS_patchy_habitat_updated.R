# 
# Analysis of SECR data for an irregular patch of habitat using JAGS
#
# For details see the blog post at
# https://www.mmeredith.net/blog/1309_SECR_in_BUGS_patchy_habitat.htm

library(rjags)
library(spatstat) # Functions to manipulate the polygon
library(secr)     # For simulation of trapping data
library(MASS)     # For equi-scale plot
library(BEST)     # For the 'hdi' (Highest Density Interval) function


# JAGS code for the model as given on the web page
# Run everything from vvvvv to ^^^^ in one go.
# (This will not work as-is in WinBUGS or OpenBUGS) vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
JAGScode <- "
model {
	sigma2 ~ dunif(0, 20)     # set up the priors for the parameters
  sigma <- sqrt(sigma2 / 2)
	lam0 ~ dgamma(0.1,0.1)
  psi ~ dunif(0, 1)
     
  for (i in 1:M){           # loop through the augmented population
    z[i] ~ dbern(psi)       # state of individual i (real or imaginary)
    SX[i] ~ dunif(0,xu)     # priors for the activity centers for each individual
    SY[i] ~ dunif(0,yu)     # lower x coordinate = 0, xu is the upper x value
    pOK[i] <- habmat[trunc(SX[i]+1), trunc(SY[i]+1)] # habitat check
    OK[i] ~ dbern(pOK[i])   # OK[i] = 1, the ones trick
    for(j in 1:J) {         # loop through the J camera trap locations
      Dsq[i,j] <- pow(SX[i]-trapmat[j,1], 2) + pow(SY[i]-trapmat[j,2],2)
      g[i,j] <- exp(-Dsq[i,j]/sigma2)
      lambda[i,j] <- K * g[i,j] * lam0 * z[i]
      y[i,j] ~ dpois(lambda[i,j])
    }
  }
  N <- sum(z[1:M])   # derive number (check that N << M)
  D <- N / A         # derive density

}"
writeLines(JAGScode, "patch_model.txt")
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# This is the patch of good habitat, surrounded by unsuitable terrain.
#   Units are metres. Total area = 2000 sq km.

patch <- data.frame(x = c(130455, 132789, 134476, 136032, 138497, 
139794, 144204, 144594, 144983, 145761, 147577, 150820, 155749, 
156398, 159122, 160289, 162494, 165218, 165478, 167553, 167812, 
168980, 170537, 172482, 174687, 176763, 178060, 179617, 181562, 
183638, 183767, 184935, 185194, 185065, 179876, 178968, 179876, 
181303, 184416, 186881, 189086, 189215, 187659, 187270, 186362, 
183897, 182730, 182211, 182470, 183249, 183508, 183249, 182859, 
181173, 179227, 178190, 178319, 176114, 175855, 175336, 174428, 
173650, 173779, 172482, 171185, 168331, 167683, 168331, 168980, 
168202, 166645, 164829, 162883, 160159, 159640, 158732, 158603, 
157824, 154971, 154063, 151987, 151209, 146669, 143686, 142648, 
141221, 139924, 137848, 136162, 136292, 134995, 133697, 131881, 
130844, 129676, 128509, 126044, 124488, 124358, 124747, 125525, 
126174, 125525, 124617, 125007, 124488, 125655, 127731, 128249, 
129287, 129287, 128898, 130065, 130325, 129936, 130455),
    y = c(802999, 
803388, 802220, 802350, 805852, 806371, 805593, 804426, 798588, 
798588, 797940, 798329, 804815, 805334, 806242, 805982, 809874, 
812209, 814154, 816878, 819991, 821937, 821548, 820381, 818824, 
818175, 815841, 812857, 809744, 806631, 804426, 802350, 799626, 
798070, 795605, 793011, 790935, 788990, 787303, 786395, 786136, 
785358, 784190, 782504, 782115, 782504, 781855, 780688, 778353, 
775240, 771737, 765641, 763306, 762398, 762917, 763825, 765771, 
770570, 772775, 773683, 773553, 771867, 766289, 765641, 766030, 
767976, 769403, 770051, 771867, 773683, 774980, 775629, 777704, 
780039, 779780, 779131, 777575, 776018, 775110, 773553, 773553, 
772905, 772775, 773683, 775759, 777185, 776667, 775369, 773294, 
771997, 771478, 771608, 771348, 772775, 774202, 774461, 774202, 
774721, 777185, 783671, 784450, 785358, 786266, 787433, 789249, 
790416, 791195, 790935, 792103, 792751, 794178, 795475, 797032, 
800794, 801702, 802999))

eqscplot(patch, type='l')

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


# CREATE THE HABITAT MATRIX
# =========================

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
      OK = rep(1, M),
      A = 2000, # sq km
      K = 90  )

# initial values:
inits <- function() list(z = rep(1, M),
              SX = rep(25, M), SY = rep(23, M) )

# estimates to extract:
wanted <- c("N", "D", "psi", "lam0", "sigma")

# Initialize the model (takes a while, 3 mins on my old laptop)
system.time(
jagsModel <- jags.model("patch_model.txt", data=dataList, inits=inits, 
                  n.chains=3, n.adapt=500)
)

mcmcList <- coda.samples(jagsModel, variable.names=wanted, n.iter=500)

# Do some diagnostic checks:
varnames(mcmcList)
gelman.diag(mcmcList[, -1]) # D and N have correlation = 1, so can't include both. 
effectiveSize(mcmcList)

result <- as.matrix(mcmcList)
# Check largest value of N << M
M
max(result[, 'N'])  # 27, OK.

# Get means and 95% HDI
colMeans(result)
hdi(result)
# In this simulation there were 8 animals in 2000 sq km =
8/2000  # per sq km
# Capture histories were simulated with g0 = 0.02 (g0 is roughly equivalent to 
#   lam0 for small values) and sigma = 2500.
# Here we have underestimated sigma and (consequently) overestimated D.


# CHECK ACTIVITY CENTRE LOCATIONS
# ===============================

# To check the AC locations, we'll get 51 values (17 per chain)
#   after thinning by 10
wantedXY <- c("SX", "SY", "z")
mcmcListXY <- coda.samples(jagsModel, variable.names=wantedXY, n.iter=170, thin=10)

res <- as.matrix(mcmcListXY)[-1, ] # Drop the first row to make 50.
dim(res)
colnames(res)

# Plot ACs of animals captured
nrow(simCH) # 4 animals captured
plot(win, main="")
for(i in 1:4) 
  points(x=res[, i], y=res[, i+50], pch=19, cex=0.3, col=i)
points(trapmat, pch='+', col='grey')
title(sub="50 simulated ACs\nfor each animal captured")

# Plot ACs for all REAL animals not captured
real <- as.vector(res[, 105:150])
nocapx <- as.vector(res[, 5:50])[real == 1]
nocapy <- as.vector(res[, 55:100])[real == 1]
plot(win, main="")
# Uncomment these lines if you want to display the habitat matrix pixels:
# for(x in 1:50)   for(y in 1:46)    if(habmat[x, y])
#       points(x-0.5, y-0.5, pch=22, col='grey', bg='grey', cex=1.3)
# plot(win, lwd=3, add=TRUE)
points(nocapx, nocapy, pch=19, cex=0.3, col='blue')
points(trapmat, pch='+', col='grey')
title(sub="Simulated ACs\nfor real animals not captured")
# ACs don't go outside the area defined by the habitat matrix.
# And they do appear to be random (with reassuring 'faces-in-the-clouds').
