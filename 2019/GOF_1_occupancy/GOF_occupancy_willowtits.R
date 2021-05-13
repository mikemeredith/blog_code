
# R code for the blog post at
#  https://mmeredith.net/blog/2019/GOF_1.htm

# SBBS willow tit occupancy : model with covariates
# Includes code for Goodness-Of-Fit and for plotting trends.

library(jagsUI)
# library(wiqid)

# Read in data
wt <- read.csv("http://mmeredith.net/data/willowtits.csv", comment="#")
summary(wt)
head(wt)

# Aggregate the detection data:
( y <- rowSums(wt[, 1:3], na.rm=TRUE) )
( n <- rowSums(!is.na(wt[, 1:3])) )

# standardise site covariates and create model matrix
covs <- cbind(inter=1, scale(wt[8:7]))
str(covs)
covs <- cbind(covs, elev2=covs[, 'elev']^2)
head(covs)

# Freeman-Tukey discrepancy
# -------------------------
# (see the file "site_covs.jags")

# Organise the data etc
# ... null model (intercept only)
jagsData <- list(y = y, n = n, nSites = length(n), covs=covs[, 1, drop=FALSE], nCovs=1)
# ... full model
# jagsData <- list(y = y, n = n, nSites = length(n), covs=covs, nCovs=ncol(covs))
str(jagsData)

inits <- function() list(z = rep(1, length(n)))

wanted <- c("p", "beta", "N", "Tobs", "Tsim")

# Run the model (6 secs)
jagsOut <- jags(jagsData, inits, wanted, "site_covs.jags", DIC=FALSE,
    n.chains=3, n.iter=5000, n.adapt=1000, n.thin=5, parallel=TRUE)
plot(jagsOut)

pp.check(jagsOut, "Tobs", "Tsim", main="Freeman-Tukey discrepancy\nnull model")
# pp.check(jagsOut, "Tobs", "Tsim", main="Freeman-Tukey discrepancy\nfull model")

# FT with resampling of z[i]
# --------------------------
# Uses the jags file "site_covs_zSim.jags"
# The rest of the code is the same.


# Deviance
# --------
# Uses the jags file "site_covs_deviance.jags"

wanted <- c("p", "beta", "N", "Dobs", "Dsim")

jagsOutD <- jags(jagsData, inits, wanted, "site_covs_deviance.jags", DIC=FALSE,
    n.chains=3, n.iter=5000, n.adapt=1000, n.thin=5, parallel=TRUE)
plot(jagsOutD)

pp.check(jagsOutD, "Dobs", "Dsim", main="Deviance\nnull model")
# pp.check(jagsOutD, "Dobs", "Dsim", main="Deviance\nfull model")


# Calculating GOF after the JAGS run
# ----------------------------------

# 'wanted' must include 'z'
wanted <- c("p", "beta", "N", "z")

jagsOut <- jags(jagsData, inits, wanted, "site_covs.jags", DIC=FALSE,
    n.chains=3, n.iter=5000, n.adapt=1000, n.thin=5, parallel=TRUE)

attach(jagsOut$sims.list)
nIter <- length(p)
Tobs <- Tsim <- numeric(nIter)
for(iter in 1:nIter) {
  psi <- plogis(beta[iter])               # for the null model
  # psi <- plogis(covs %*% beta[iter, ])  # for the full model
  Tobs[iter] <- sum((sqrt(y) - sqrt(p[iter]*z[iter, ]*n))^2)
  ySim <- rbinom(237, n, p[iter]*z[iter, ])
  Tsim[iter] <- sum((sqrt(ySim) - sqrt(p[iter]*z[iter, ]*n))^2)
}
detach(jagsOut$sims.list)

MASS::eqscplot(Tobs, Tsim, xlim=range(Tobs, Tsim), ylim=range(Tobs, Tsim),
  xlab="Observed data", ylab="Simulated data")
abline(0, 1, lwd=2, col='red')
mean(Tsim > Tobs)
