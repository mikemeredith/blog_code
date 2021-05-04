
# Animation of Gibbs sampler with salamander data

# Code for the web page at
#  https://mmeredith.net/blog/2015/Gibbs_sampler.htm

# The salamanders data:
# ---------------------
S <- 39    # Number of sites
K <- 5     # Number of visits to each site
xObs <- 18 # Number of sites where salamanders were detected
d <- 30    # total number of detections
# This is all we need, we don't need the site-by-site data.

# Priors: independent beta(1, 1) priors
# -------------------------------------
priPsi <- c(1, 1)
priPi <- c(1, 1)

# Run the Gibbs sampler
# ---------------------
nIter <- 10100               # number of iterations
psi <- pi <- numeric(nIter) # objects to hold results
psi[1] <-  pi[1] <- 0.5     # starting values

for(i in 2:nIter) {
  # 1. Calculate prob(occupied | not detected):
  psi0 <- (psi[i-1] * (1 - pi[i-1])^K) / (psi[i-1] * (1 - pi[i-1])^K + (1 - psi[i-1]))
  # ...and draw number of additional sites occupied
  xAdd <- rbinom(1, S - xObs, psi0)
  x <- xObs + xAdd
  # 2a. Draw new psi from beta(occupied+prior, unoccupied+prior)
  psi[i] <- rbeta(1, x + priPsi[1], S - x + priPsi[2])
  # 2b. Draw new pi from beta(detected+prior, undetected+prior).
  pi[i] <- rbeta(1, d + priPi[1], x * K - d + priPi[2])
}

# Plot the results
# ----------------
plot(psi, pi, xlim=c(0.2, 1), ylim=c(0, 0.8), type='l', col='grey',
  xlab="occupancy, psi", ylab="detection, pi")
points(psi, pi, pch=16, cex=0.1)
points(psi[1], pi[1], pch=16, col='red')  # starting point

# Remove first 100 iterations (burn-in)
# -------------------------------------
psi <- psi[101:10100]
pi <- pi[101:10100]

# Plot histograms with mean and SD
# -----------------------------------
par(mfrow=1:2)
if(require(mcmcOutput)) {
  postPlot(psi)
  postPlot(pi)
} else {
  hist(psi)
  abline(v=mean(psi), col='red')
  abline(v=quantile(psi, c(0.025, 0.975)), col='red', lty=3)

  hist(pi)
  abline(v=mean(pi), col='red')
  abline(v=quantile(pi, c(0.025, 0.975)), col='red', lty=3)
}
par(mfrow=c(1,1))

# Relationship between psi and pi
# -------------------------------
cor(psi, pi)

if(require(emdbook)) {
  plot(psi, pi, pch=19, cex=0.1)
  HPDregionplot(cbind(psi, pi), col="red", lwd=3, add=TRUE)
}

# Trace plots:
# ------------
par(mfrow=1:2)
plot(psi[1:200], type='l')
abline(h=mean(psi), col='red')
plot(pi[1:200], type='l')
abline(h=mean(pi), col='red')
par(mfrow=c(1, 1))

head(cbind(psi, pi))

# Autocorrelation plots and effective size
# ----------------------------------------
par(mfrow=1:2)
acf(psi)
acf(pi)
par(mfrow=c(1, 1))

library(coda)
effectiveSize(psi)
effectiveSize(pi)
