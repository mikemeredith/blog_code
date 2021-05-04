
# Animation of random walk sampler with salamander data

# Code for the web page at
#   https://mmeredith.net/blog/2015/RandomWalk_MCMC.htm

# The salamander data:
# ---------------------
y <- c(2, 1, 1, 4, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
n <- 5

# Function to calculate likelihood:
get_like <- function(psi, pi, y, n) {
  prod(psi * dbinom(y, n, pi) + ifelse(y==0, 1 - psi, 0))
}

# Priors: independent beta(1, 1) priors
# -------------------------------------
priorPsi <- function(psi, shape1=1, shape2=1) dbeta(psi, shape1, shape2)
priorPi <- function(pi, shape1=1, shape2=1) dbeta(psi, shape1, shape2)
# Change shape1 and shape2 to match your prior distribution.

# Run the random walk sampler
# ---------------------------
nIter <- 10100               # number of iterations
# stepmaxPsi <- 0.1   # maximum step size (too small)
# stepmaxPi <- 0.1
# stepmaxPsi <- 0.6   # maximum step size (too big)
# stepmaxPi <- 0.6
stepmaxPsi <- 0.3   # maximum step size (optimum)
stepmaxPi <- 0.2
psi <- pi <- numeric(nIter) # objects to hold results
psi[1] <-  pi[1] <- 0.5     # starting values
likelihood <- get_like(psi[1], pi[1], y, n)  # initial likelihood

for(i in 2:nIter) {
  # 1. Update psi:
  # Generate candidate value
  cand <- psi[i-1] + runif(1, -stepmaxPsi, stepmaxPsi)
  # If candidate value is outside [0,1], we definately reject
  if (cand < 0 | cand > 1) {
    A <- 0   # acceptance probability
  } else {
    # Calculate likelihood at candidate value
    candLlh <- get_like(cand, pi[i-1], y, n)
    # Calculate likelihood * prior at old value and candidate value
    jointOld <- likelihood * priorPsi(psi[i-1])
    jointCand <- candLlh * priorPsi(cand)
    # Acceptance probability is ratio, or 1 if ratio > 1
    A <- min(1, jointCand / jointOld)
  }
  # Decide whether to accept or not
  if(A > runif(1)) {   # if accepted
    psi[i] <- cand
    likelihood <- candLlh
  } else {
    psi[i] <- psi[i-1]
  }

  # 1. Update pi:
  # Generate candidate value
  cand <- pi[i-1] + runif(1,  -stepmaxPi, stepmaxPi)
  # If candidate value is outside [0,1], definately reject
  if (cand < 0 | cand > 1) {
    A <- 0
  } else {
    # Calculate likelihood at candidate value
    candLlh <- get_like(psi[i], cand, y, n)
    # Calculate likelihood * prior at old value and candidate value
    jointOld <- likelihood * priorPi(pi[i-1])
    jointCand <- candLlh * priorPi(cand)
    # Acceptance probability is ratio, or 1 if ratio > 1
    A <- min(1, jointCand / jointOld)
  }
  # Decide whether to accept or not
  if(A > runif(1)) {   # if accepted
    pi[i] <- cand
    likelihood <- candLlh
  } else {
    pi[i] <- pi[i-1]
  }
}

# Plot the results
# ----------------
plot(psi, pi, xlim=c(0.2, 1), ylim=c(0, 0.8), type='l', col='grey',
  xlab="occupancy, psi", ylab="detection, pi")
points(psi, pi, pch=16, cex=0.1)
points(psi[1], pi[1], pch=16, col='red')  # starting point


# Remove first 100 iterations (burn-in)
# -------------------------------------
psi <- psi[-(1:100)]
pi <- pi[-(1:100)]

# Trace plots:
# ------------
par(mfrow=1:2)
plot(psi[1:200], type='l')
abline(h=mean(psi), col='red')
plot(pi[1:200], type='l')
abline(h=mean(pi), col='red')
par(mfrow=c(1, 1))

# Effective size, autocorrelation plots, and acceptance rate
# ----------------------------------------------------------
library(coda)
effectiveSize(psi)
effectiveSize(pi)

par(mfrow=1:2)
acf(psi)
acf(pi)
par(mfrow=c(1, 1))

# Acceptance rate:
mean(diff(psi) != 0)
mean(diff(pi) != 0)


# Plot histograms with mean and credible interval
# -----------------------------------------------
# (Do this after checking that the effective sample soze is big enough.)
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
