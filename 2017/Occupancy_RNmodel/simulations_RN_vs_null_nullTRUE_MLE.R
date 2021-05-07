
# Comparison of Royle-Nichols model vs psi(.)p(.) model with MLE

# This time the data are generated from the psi(.)p(.) model

library(wiqid)

library(parallel)
library(doParallel)
ncore <- 7 # on my machine
cl <- makeCluster(ncore)
clusterSetRNGStream(cl = cl)
registerDoParallel(cl)

# Parameters for simulated data:
nSite <- 200   # number of sites
nRep <- 10     # number of replicate surveys
n <- rep(nRep, 200)
lambda <- 1.5  # mean number of individuals at each site
r <- 0.1       # probability of detection of one individual

tmp <- rpois(1e7, lambda)  # simulate N for 1e7 sites
tmp <- tmp[tmp > 0]        # remove sites where species is absent
p.tmp <- 1 - (1 - r)^tmp   # p for each site
( p <- mean(p.tmp) )       # mean p.detect | present for current lambda and r; 0.179


## Do the loop
## ===========

nSims <- 1001 # 205
system.time( 
result <- foreach(x=1:nSims, .combine=rbind, .packages="wiqid", 
    .errorhandling='pass') %dopar% {
  out <- numeric(7)
  
  # Generate simulated data
  out[1] <- x
  set.seed(x)
  # Ecological process:
  N <- rpois(nSite, lambda)  # Latent occurrence state
  z <- as.numeric(N > 0)
  out[2] <- mean(z)  # True occupancy
  # Observation process:
  # p <- 1 - (1-r)^N
  y <- rbinom(nSite, nRep, p * z)
  out[3] <- mean(y > 0)  # Naive estimate
  
  # Run wiqid - simple model
  null <- occSS0(y, n)
  out[4] <- null$real[1,1]
  out[5] <- AIC(null)

  # Run wiqid - RN model
  rn <- occSSrn0(y, n)
  out[6] <- rn$real[1,1]
  out[7] <- AIC(rn)
  
  out
} )  # 1000 took 80 secs

colnames(result) <- c("seed", "psiSample", "psiNaive",
  "psi0", "AIC0", "psiRN", "AICRN")
write.csv(result, "simulation_RN_vs_null_nullTRUE_MLE.csv")
dim(result)
head(result)

colSums(is.na(result))
colMeans(result, na.rm=TRUE)

# Bias and RMSE with data-generating value
( psiTrue <- 1 - exp(-lambda) )
mean(result[, 'psi0']) - psiTrue     # 0.002
mean(result[, 'psiRN']) - psiTrue    # 0.122
sqrt(mean((result[, 'psi0'] - psiTrue)^2))   # 0.0421
sqrt(mean((result[, 'psiRN'] - psiTrue)^2))  # 0.1337
# What's "RMSE" for sample proportion
sqrt(mean((result[, 'psiSample'] - psiTrue)^2)) # 0.029

# Bias and RMSE with sample-specific value
err0 <- result[, 'psi0'] - result[, 'psiSample']
errRN <- result[, 'psiRN'] - result[, 'psiSample']
mean(err0)   # 0.0019
mean(errRN)  # 0.1227
sqrt(mean(err0^2))   # 0.0311
sqrt(mean(errRN^2))  # 0.1297

# AIC
mean(result[, "AICRN"] < result[, "AIC0"])





