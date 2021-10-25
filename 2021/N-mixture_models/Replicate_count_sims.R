
# Simulate replicate count data for multiple sites with imperfect detection
#  and analyse with different models.

# Analysis models
#  Poisson
#  Zero-inflated Poisson (ZIP)
#  Negative binomial (NB)
#  Zero-inflated negative binomial (ZINB)
ananames <- c("Pois", "ZIP", "NB", "ZINB")
nana <- length(ananames)
modfiles <- c("Nmix_Pois.jags", "Nmix_ZIP.jags", "Nmix_NB.jags", "Nmix_ZINB.jags")

# Data-generating scenarios for true N:
#  Poisson with mean = variance = 3.5 (SD = 1.87)
mu <- 3.5
#  ZIP1 with zero-inflation parameter 0.5
zinf1 <- 0.5
#  ZIP2 with zero-inflation parameter 0.7
zinf2 <- 0.7
#  NB1 with dispersion parameter 0.45 (SD = 3.00)
disp1 <- 0.45
#  NB2 with dispersion parameter 1.0 (SD = 3.97)
disp2 <- 1
#  ZINB with zero-inflation parameter 0.5 and dispersion parameter 0.45
gennames <- c("Pois", "ZIP1", "ZIP2", "NB1", "NB2", "ZINB")
ngen <- length(gennames)

library(jagsUI)

# Do lots of simulations
# ======================
library(foreach)
library(parallel)
library(doParallel)

detectCores()  # number of cores on your machine
ncore <- 20     # <-- adjust this for your machine

cl <- makeCluster(ncore)
registerDoParallel(cl)

# True parameters
nSites <- 300      # Spatial replicates
nOcc <- 3          # number of replicate counts
pdet <- 0.6        # probability of detection of individual at one visit

# Stuff we need inside the loop
trueN <- matrix(NA, nSites, ngen)
colnames(trueN) <- gennames
C <- array(NA, dim=c(nSites, nOcc, ngen))

# All models include the following parameters
# mu     : expected number of animals per site (excluding structural 0s for ZI models)
# zinf   : proportion of sites with structural zeros for ZIP/ZINB models, 0 otherwise
# disp   : overdispersion in NB/ZINB models, 0 otherwise
# p      : probability of detection for each individual per visit
# Ntotal : estimated total number of individuals, summed across all sites
# FT.obs : Freeman-Tukey discrepancy for observed counts
# FT.sim : Freeman-Tukey discrepancy for simulated counts
wanted <- c("mu", "zinf", "disp", "p", "Ntotal", "FT.obs", "FT.sim")
npars <- length(wanted)

# The output also includes
# trueNtotal : true total in the simulated population (independent of analysis method)
# GOF.p      : the proportion of discrepancies for simulated data that exceed the
#              discrepancy for the observed data.

result <- array(NA, dim=c(nana, ngen, npars))
dimnames(result) <- list(analysis=ananames, generation=gennames,
    param=c(wanted[1:5], "trueNtotal", "GOF.p"))


# seeds <- 1:20
# seeds <- 21:200
# seeds <- 201:500
seeds <- 21:500
system.time(
res5 <- foreach(x = seeds, .combine=rbind, .packages=c("jagsUI"),
  .errorhandling='remove') %dopar% {
  set.seed(x)
  # Generate populations
  trueN[,1] <- rpois(nSites, mu)               # Poisson
  trueN[,2] <- rpois(nSites, mu)
  trueN[,2][sample(nSites, nSites*zinf1)] <- 0 # ZIP1
  trueN[,3] <- rpois(nSites, mu)
  trueN[,3][sample(nSites, nSites*zinf2)] <- 0 # ZIP2
  trueN[,4] <- rnbinom(nSites, 1/disp1, mu=mu) # NB1
  trueN[,5] <- rnbinom(nSites, 1/disp2, mu=mu) # NB2
  trueN[,6] <- rnbinom(nSites, 1/disp1, mu=mu)
  trueN[,6][sample(nSites, nSites*zinf1)] <- 0 # ZINB

  # Get counts
  for(i in 1:ngen)
    C[,,i] <- rbinom(nSites*nOcc, trueN[,i], pdet)

  # Do the analysis
  for(ana in 1:nana) {
    for(gen in 1:ngen) {
      Cmax <- apply(C[,,gen], 1, max)
      jdata <- list(nSites = nSites, C = C[,,gen], nOcc = nOcc,
          z = ifelse(Cmax > 0, 1, NA))
      inits <- function() list(N = Cmax)
      out <- jags(jdata, inits, wanted, modfiles[ana], DIC=FALSE,
          n.chains=3, n.adapt=5000, n.iter=20000, parallel=FALSE)
      P <- mean(out$sims.list$FT.sim > out$sims.list$FT.obs)
      tmp <- unlist(out$mean)
      Rhat <- unlist(out$Rhat)
      tmp[Rhat > 1.1] <- NA
      result[ana, gen, ] <- c(tmp[1:5], sum(trueN[, gen]), round(P, 3))
    }
  }

  # Return results as a vector
  c(result)
} )  # 20 took 85 mins, 180 took 13 hrs
save(res5, file="output_res5.RData")

# stopCluster(cl)

# Combine the different runs
# --------------------------
resfst <- rbind(res1, res2, res3)
niter <- nrow(resfst)
dim(resfst)
dim(resfst) <- c(niter, 4, 5, 7)

resmore <- rbind(res4, res5)
niter == nrow(resmore)
dim(resmore)
dim(resmore) <- c(niter, 4, 6, 7)
dimnames(resmore) <- c(list(iter=NULL), dimnames(result))
str(resmore)
resall <- resmore
resall[,,1:5,] <- resfst
str(resall)
resall[,,,'GOF.p'] <- round(resall[,,,'GOF.p'], 3)
save(resall, file="resall_500.RData")
# Post-processing
# ===============

# Check for NAs
sum(is.na(resall))  # 2
which(is.na(resall), arr.ind = TRUE)  # both ZINB analysis of ZIP2 data: mu and Ntotal

( means <- apply(resall, 2:4, mean, na.rm=TRUE) )
dim(means)
round(means[,,1:4], 2)
round(means[,,5:6])

# Check GOF
GOF <- resall[,,,'GOF.p']
GOFok <- GOF > 0.1 & GOF < 0.9
apply(GOFok, 2:3, mean)

# estimates of bias and error
error <- (resall[,,,'Ntotal'] - resall[,,,'trueNtotal'])
bias <- error / resall[,,,'trueNtotal']
mean_bias <- apply(bias, 2:3, mean, na.rm=TRUE)*100 # percent
round(mean_bias, 2)

RMSE <- sqrt(apply(error^2, 2:3, mean, na.rm=TRUE))
RMSEpc <- RMSE / means[,,'trueNtotal'] * 100 # percent
round(RMSEpc, 2)

dim(error)
par(mfcol=c(4,6), mar=c(2,2,5,1))#, oma = c(2,2,4,0))
for(scen in 1:6) {
  for(fit in 1:4) {
    hist(error[,fit, scen], main=ananames[fit], xlab="Error")
    abline(v=0, lwd=3, col='red')
    abline(v=mean(error[,fit, scen], na.rm=TRUE), lwd=3, col='blue', lty=2)
    if(fit==1)
      title(main=paste("Truth =", gennames[scen]), line=3, xpd=TRUE)
  }
}
