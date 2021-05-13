
# R code for blog post at
#  https://mmeredith.net/blog/2021/parallel_chains_foreach.htm

# Blue ridge salamander occupancy analysis

library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14"
library(mcmcOutput)

# Create vector of number of detections
( y <- rep(4:0, c(1,4,1,12,21)) )
n <- 5                               # No. of surveys at each site
( nSites <- length(y) )

cat(file="occu.bugs", "
model{
  # Likelihood
  for(i in 1:nSites) {
     # Ecological model
     z[i] ~ dbern(psi)  # z=1 if occupied, z=0 if not occupied
     # Observation model
     pz[i] <- p * z[i]
     y[i] ~ dbin(pz[i], n)
  }

  # Priors
  p ~ dbeta(1, 1) # Uninformative prior
  psi ~ dbeta(1, 1)

  # Derived values
  N <- sum(z[])
}")

data <- list(y = y, n = n, nSites = nSites)
inits <- function() list(z = rep(1, nSites))
wanted <- c("p", "psi", "N")

( WB1 <- bugs(data=data, inits=inits, parameters=wanted, model="occu.bugs",
    n.iter=11000, n.burnin=1000, n.chains=3, DIC=FALSE,
    bugs.dir=bugs.dir, debug=TRUE) )
diagPlot(WB1)

# Parallel processing with WinBUGS and foreach
# ============================================

library(doParallel) # also attaches 'foreach' and 'parallel'
detectCores()  # Number of cores available on this machine

ncore <- 3     # Number to use
cl <- makeCluster(ncore)
registerDoParallel(cl)

seeds <- 1:ncore

result <- foreach(x = seeds, .combine=rbind, .packages="R2WinBUGS") %dopar% {
  set.seed(x)
  wb <- bugs(data=data, inits=inits, parameters=wanted, model="occu.bugs",
      n.iter=11000, n.burnin=1000, n.chains=1, DIC=FALSE,
      bugs.dir=bugs.dir, debug=FALSE)
  wb$sims.matrix
}
str(result) # one big matrix

library(mcmcOutput)
( mco <- mcmcOutput(result, nChains=ncore) )
summary(mco)
diagPlot(mco)
plot(mco)
