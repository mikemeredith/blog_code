
# Bivariate normal priors in JAGS
# ===============================

# See blog post at
# https://mmeredith.net/blog/2020/Correlated_priors.htm

library(rjags)
library(mcmcOutput)

# 1. Linked univariate normal priors
# ==================================

wanted <- c("sigma", "rho", "eta")

jm <- jags.model("BV1_link.jags", data=NULL)
out <- coda.samples(jm, wanted, 1e3)

( mc <- mcmcOutput(out) )
# Object of class 'mcmcOutput'; approx. size 8.1 MB 
# MCMC values from mcmc.list object ‘out’ 
# The output has 1 chains each with 1000 draws.
# It has 3 parameters with 1003 nodes monitored:
#       nodes
# eta    1000
# rho       1
# sigma     2

# Calculate the means, SDs and correlation coefficients for the realised priors
eta <- mc$eta
niter <- nrow(mc)
meanx <- sdx <- matrix(NA, niter, 2)
corx <- numeric(niter)
for(n in 1:niter) {
  meanx[n, ] <- colMeans(eta[n,,])
  sdx[n, ] <- apply(eta[n,,], 2, sd)
  corx[n] <- cor(eta[n,,])[1,2]
}

# Compare the realised values of SD and correlation with sigma and rho:
plot(mc$sigma[,1], sdx[,1])
plot(mc$sigma[,2], sdx[,2])
plot(mc$rho, corx)

hist(sdx[,1])  # uniform(0, 5)
hist(sdx[,2])  # uniform(0, 5)
hist(corx)     # uniform (-1, 1)


# 2. dmvnorm with a Wishart (dwish) prior
# =======================================

( R <- diag(c(5,1)) )
#      [,1] [,2]
# [1,]    5    0
# [2,]    0    1
( sigma.guess <- sqrt(diag(R)/3) )
# [1] 1.2909944 0.5773503

jdat <- list(R = R, k = 3)

wanted <- c("sigma", "rho", "eta")

jm <- jags.model("BV2_Wishart.jags", data=jdat)
out <- coda.samples(jm, wanted, 1e3)

( mc <- mcmcOutput(out) )

# Calculate the means, SDs and correlation coefficients for the realised priors
eta <- mc$eta
niter <- nrow(mc)
meanx <- sdx <- matrix(NA, niter, 2)
corx <- numeric(niter)
for(n in 1:niter) {
  meanx[n, ] <- colMeans(eta[n,,])
  sdx[n, ] <- apply(eta[n,,], 2, sd)
  corx[n] <- cor(eta[n,,])[1,2]
}

# Check shape of realised distribution of SD
#    add 'prior guess' and 90th percentile
op <- par(mfrow=c(1,2))
hist(sdx[,1], main="sd[1]", breaks=50)
abline(v=sigma.guess[1], col='red', lwd=2)
abline(v=quantile(sdx[,1], 0.9), col='blue', lwd=2, lty=2)
hist(sdx[,2], main="sd[2]", breaks=50)
abline(v=sigma.guess[2], col='red', lwd=2)
abline(v=quantile(sdx[,2], 0.9), col='blue', lwd=2, lty=2)
par(op)

# 3. dmvnorm with a scaled Wishart (dscaled.wishart) prior
# ========================================================

# Plots of half-t
library(wiqid)
curve(2*dt2(x, 0, 1, 2), 0, 20, ylab="Density", las=1, lwd=2,
    main="t-distribution with df = 2")
abline(v=qt2(0.95, 0, 1, 2))
curve(2*dt2(x, 0, 2, 2), col=2, add=TRUE, lwd=2)
abline(v=qt2(0.95, 0, 2, 2), col=2)
curve(2*dt2(x, 0, 3, 2), col=3, add=TRUE, lwd=2)
abline(v=qt2(0.95, 0, 3, 2), col=3)
curve(2*dt2(x, 0, 5, 2), col=4, add=TRUE, lwd=2)
abline(v=qt2(0.95, 0, 5, 2), col=4)
legend('topright', paste("scale =", c(1,2,3,5)), col=1:4,
  lwd=2, bty='n')
legend('right', "Vertical lines show 90th percentile", bty='n')


jdata <- list(s = c(5,1), df=2)

wanted <- c("sigma", "rho", "eta")

load.module("glm")
jm <- jags.model("BV3_scaled.jags", data=jdata)
out <- coda.samples(jm, wanted, 1e3)

# Calculate the means, SDs and correlation coefficients for the realised priors
( mc <- mcmcOutput(out) )
eta <- mc$eta
niter <- nrow(mc)
meanx <- sdx <- matrix(NA, niter, 2)
corx <- numeric(niter)
for(n in 1:niter) {
  meanx[n, ] <- colMeans(eta[n,,])
  sdx[n, ] <- apply(eta[n,,], 2, sd)
  corx[n] <- cor(eta[n,,])[1,2]
}

# Check shape of realised distribution of SD
op <- par(mfrow=c(1,2))
hist(sdx[,1], main="sd[1]", breaks=50)
abline(v=quantile(sdx[,1], 0.9), col='blue', lwd=2, lty=2)
hist(sdx[,2], main="sd[2]", breaks=50)
abline(v=quantile(sdx[,2], 0.9), col='blue', lwd=2, lty=2)

#   add t2 curve and 90th percentile 
sd <- sdx[,1][sdx[,1] < 50]
hist(sd, main="sd[1]", breaks=30, freq=FALSE, las=1)
curve(2*dt2(x, 0, 5, 2), col='red', lwd=2, add=TRUE)
abline(v=quantile(sdx[,1], 0.9), col='blue', lwd=2, lty=2)

sd <- sdx[,2][sdx[,2] < 10]
hist(sd, main="sd[2]", breaks=30, freq=FALSE, las=1)
curve(2*dt2(x, 0, 1, 2), col='red', lwd=2, add=TRUE)
abline(v=quantile(sdx[,2], 0.9), col='blue', lwd=2, lty=2)
par(op)


# 4. dmvnorm.vcov with a variance-covariance matrix prior
# -------------------------------------------------------

wanted <- c("sigma", "rho", "eta")

load.module("glm")
jm <- jags.model("BV4_VCOV.jags", data=NULL)
out <- coda.samples(jm, wanted, 1e3)

# Calculate the means, SDs and correlation coefficients for the realised priors
( mc <- mcmcOutput(out) )
eta <- mc$eta
niter <- nrow(mc)
meanx <- sdx <- matrix(NA, niter, 2)
corx <- numeric(niter)
for(n in 1:niter) {
  meanx[n, ] <- colMeans(eta[n,,])
  sdx[n, ] <- apply(eta[n,,], 2, sd)
  corx[n] <- cor(eta[n,,])[1,2]
}

# Check shape of realised distribution of SD
op <- par(mfrow=c(1,2))
hist(sdx[,1], main="sd[1]")
hist(sdx[,2], main="sd[2]")
par(op)


