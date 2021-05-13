
# R code for blog post at
#  https://mmeredith.net/blog/2021/parallel_chains_foreach.htm

# Blue ridge salamander occupancy analysis

library(nimble)
library(mcmcOutput)

# Create vector of number of detections
( y <- rep(4:0, c(1,4,1,12,21)) )
n <- 5                               # No. of surveys at each site
( nSites <- length(y) )

code <- nimbleCode({
  # Likelihood
  for(i in 1:nSites) {
     # Ecological model
     z[i] ~ dbern(psi)  # z=1 if occupied, z=0 if not occupied
     # Observation model
     y[i] ~ dbin(p * z[i], n)
     # pz[i] <- p * z[i]
     # y[i] ~ dbin(pz[i], n)
  }

  # Priors
  p ~ dbeta(1, 1) # Uninformative prior
  psi ~ dbeta(1, 1)

  # Derived values
  N <- sum(z[])
})

data <- list(y = y)
const <- list(n = n, nSites = nSites)
inits <- list(z = rep(1, nSites))
monitor <- c("p", "psi", "N")

ML1 <- nimbleMCMC(code, data=data, constants=const, inits=inits, monitors=monitor,
    niter=11000, nburnin=1000, nchains=3, samplesAsCodaMCMC=TRUE)
class(ML1)  # mcmc.list


# Parallel processing with nimble and foreach
# ===========================================

library(doParallel) # also attaches 'foreach' and 'parallel'
detectCores()  # Number of cores available on this machine

ncore <- 3     # Number to use
cl <- makeCluster(ncore)
registerDoParallel(cl)

seeds <- 1:ncore

result <- foreach(x = seeds, .packages="nimble") %dopar% {
  set.seed(x)
  nimbleMCMC(code, data=data, constants=const, inits=inits, monitors=monitor,
      niter=11000, nburnin=1000, nchains=1, samplesAsCodaMCMC=TRUE)
}
str(result) # list of `mcmc` objects

ML2 <- coda::mcmc.list(result)

library(mcmcOutput)
( mco <- mcmcOutput(ML2) )
summary(mco)
diagPlot(mco)
plot(mco)
