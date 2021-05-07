# Royle-Nichols model in JAGS with simulated data
#  plus comparison of RN model with psi(.)p(.) model

library(jagsUI)
library(loo)

library(parallel)
library(doParallel)
detectCores()
ncore <- 7
cl <- makeCluster(ncore)
clusterSetRNGStream(cl = cl)
registerDoParallel(cl)

# Parameters for simulated data:
nSite <- 200   # number of sites
nRep <- 10     # number of replicate surveys
lambda <- 1.5  # mean number of individuals at each site
r <- 0.1       # probability of detection of one individual

# Define models vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
modelRN <- "
model {
  # Likelihood
  for (i in 1:nSite) {
    # Process model
    N[i] ~ dpois(lambda)
    z[i] <- step(N[i] - 1) # z=1 if N>0, ie. site is occupied

    # Observation model
    p[i] <- 1 - pow(1-r, N[i])
    y[i] ~ dbin(p[i], nRep)

    # Loglikelihood for LOO and WAIC calculation
    loglik[i] <- log(dbin(y[i], p[i], nRep) * psi +
            step(-y[i]) * (1 - psi))

  }
  # Priors
  lambda ~ dunif(0, 10)
  r ~ dunif(0, 1)

  # Derived quantities
  psi.sample <- mean(z[])   # Proportion of occupied sites in the sample
  psi <- 1 - exp( -lambda)  # Prob(N > 0 | lambda), probability that a random site is occupied.
} "
writeLines(modelRN, con="modelRN.jags")

model0 <- "
model {
  # Likelihood
  for (i in 1:nSite) {
    # Process model
    z[i] ~ dbern(psi) # z=1 if N>0, ie. site is occupied

    # Observation model
    y[i] ~ dbin(p * z[i], nRep)

    # Loglikelihood for LOO and WAIC calculation
    loglik[i] <- log(dbin(y[i], p, nRep) * psi +
            step(-y[i]) * (1 - psi))
  }
  # Priors
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)

  # Derived quantities
  psi.sample <- mean(z[])   # Proportion of occupied sites in the sample
} "
writeLines(model0, con="model0.jags")
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

## Do the loop
## ===========

# nSims <- 1001 # 143 per core * 7
nSims <- ncore * 2 # for testing, 8 mins

system.time(
result <- foreach(x=1:nSims, .combine=rbind, .packages=c("jagsUI", "loo"),
    .errorhandling='pass') %dopar% {
  out <- numeric(13)

  # Generate simulated data
   out[1] <- x
  set.seed(x)
  # Ecological process:
  N <- rpois(nSite, lambda)  # Latent occurrence state
  out[2] <- mean(N > 0)  # True occupancy
  # Observation process:
  p <- 1 - (1-r)^N
  y <- rbinom(nSite, nRep, p)
  out[3] <- mean(y > 0)  # Naive estimate

  # Run JAGS - simple model
  JAGSdata <- list(y = y, nRep = nRep, nSite = nSite)
  zst <- as.numeric(y > 0)  # z must be 1 if y > 0, so start at 1
  inits <- function() list(z = zst)
  params <- c("p", "psi", "psi.sample", "loglik")
  JAGSout0 <- jags(JAGSdata, inits, params, "model0.jags", DIC=FALSE,
    n.chains = 3, n.adapt=1000, n.burnin=0, n.iter=10000, n.thin=1,
    codaOnly="loglik", parallel=FALSE)
  out[4] <- JAGSout0$mean$psi
  loglik <- JAGSout0$sims.list$loglik
  loo <- loo(loglik, cores=1)
  out[5] <- loo$looic
  out[6] <- sum(loo$pareto_k > 1)
  waic <- waic(loglik)
  out[7] <- waic(loglik)$waic
  out[8] <- sum(waic$pointwise[, 'p_waic'] > 0.4)

  # Run JAGS - RN model
  # JAGSdata <- list(y = y, nRep = nRep, nSite = nSite)
  Nst <- as.numeric(y > 0)  # N must be >0 if y > 0, so start at 1
  inits <- function() list(N = Nst, lambda=runif(1,0,3))
  params <- c("psi", "psi.sample", "loglik")
  JAGSoutRN <- jags(JAGSdata, inits, params, "modelRN.jags", DIC=FALSE,
    n.chains = 3, n.adapt=1000, n.burnin=0, n.iter=50000, n.thin=10,
    codaOnly="loglik", parallel=FALSE)
  out[9] <- JAGSoutRN$mean$psi
  loglik <- JAGSoutRN$sims.list$loglik
  loo <- loo(loglik, cores=1)
  out[10] <- loo$looic
  out[11] <- sum(loo$pareto_k > 1)
  waic <- waic(loglik)
  out[12] <- waic(loglik)$waic
  out[13] <- sum(waic$pointwise[, 'p_waic'] > 0.4)

  out
} )  # 143 * 7 took 5.4 hrs

colnames(result) <- c("seed", "psiSample", "psiNaive",
  "psi0", "loo0", "loo0bad", "waic0", "waic0bad",
  "psiRN", "looRN", "looRNbad", "waicRN", "waicRNbad")
write.csv(result, "RN_vs_null_RNtrue_JAGS.csv")
dim(result)
head(result)

colSums(is.na(result))
colMeans(result, na.rm=TRUE)

# Bias and RMSE with data-generating value
( psiTrue <- 1 - exp(-lambda) )
mean(result[, 'psi0']) - psiTrue     # -0.0766
mean(result[, 'psiRN']) - psiTrue    # 0.006022
sqrt(mean((result[, 'psi0'] - psiTrue)^2))   # 0.0868
sqrt(mean((result[, 'psiRN'] - psiTrue)^2))  # 0.0533
# What's "RMSE" for sample proportion
sqrt(mean((result[, 'psiSample'] - psiTrue)^2)) # 0.029

# Bias and RMSE with sample-specific value
err0 <- result[, 'psi0'] - result[, 'psiSample']
errRN <- result[, 'psiRN'] - result[, 'psiSample']
mean(err0)   # -0.076
mean(errRN)  # 0.0066
sqrt(mean(err0^2))   # 0.0827
sqrt(mean(errRN^2))  # 0.0439

# LOOaic and WAIC
# How many 'bad's?
colMeans(result[, c("loo0bad", "waic0bad", "looRNbad", "waicRNbad")] > 0)
  # No bads for looaic; waic null model 1/1001, RN model waic 100% bad.
# How do they compare
mean(result[, "looRN"] < result[, "loo0"])
mean(result[, "waicRN"] < result[, "waic0"])
