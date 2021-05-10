
# R code for the blog post at
#  https://mmeredith.net/blog/2019/MSOM_CrossVal.htm

# Load and inspect the data
load(url("http://mmeredith.net/data/fishData.RData"))
ls()
README
str(Y)
str(Psi.cov)
str(X)

modelPsi
modelP

# Running the analysis
# --------------------

# Run all 5 models without cross-validation

library(jagsUI)
wanted <- c("muPsiCoef", "sigmaPsiCoef",
            "muDetCoef", "sigmaDetCoef",
            "detCoef", "psiCoef")
modelOut <- list()
for(i in seq_along(modelPsi)) {
  inits <- function() list(Z = matrix(1, nSpecies, nSites),
       halfsigmaPsiCoef = runif(length(modelPsi[[i]]), 0, 5),
       halfsigmaDetCoef = runif(length(modelP[[i]]), 0, 5))

  jdata <- list(Y = Y, nSurveys = nSurveys,
   Psi.cov = Psi.cov[, modelPsi[[i]]], Det.cov = X[, , modelP[[i]]])

  modelOut[[i]] <- jags(jdata, inits, wanted, "MSOM_KnownSps.jags",
    n.chains=3, n.adapt=1000, n.iter= 5000, parallel=TRUE)
}
for(i in seq_along(modelOut))
  mcmcOutput::diagPlot(modelOut[[i]], param=1:18, main=names(modelPsi)[i])

# Dividing the date set into folds

nFolds <- 5
( foldSize <- nSites %/% nFolds + 1 )
tmp <- matrix(0, nFolds, foldSize)
set.seed(123) # Note that a new RNGkind is used for sample.int from R 3.6.0
for(i in 1:foldSize)
  tmp[, i] <- sample.int(nFolds)
foldID <- as.vector(tmp)[1:nSites]
table(foldID)

# Scoring the fit of the predictions to the observed data

# probs : an array of probabilities
# preds : an array of 1/0 predictions
# bool  : a matching array of observed TRUE/FALSE or 1/0 values

# output: scalar, the sum or mean score

# 0-1 loss function
score01 <- function(preds, bool)
  sum(abs(preds - bool), na.rm=TRUE)

# Area under ROC curve
auroc <- function(probs, bool) {
  ok <- is.finite(bool) & is.finite(probs)  # remove NAs
  bool <- bool[ok]
  probs <- probs[ok]
  (mean(rank(probs)[bool==1]) - (sum(bool)+1)/2) / sum(bool==0)
}

# Brier score
brier <- function(probs, bool) {
  sum(bool * (1 - probs)^2 + (1 - bool) * probs^2, na.rm=TRUE)
}

# Log(likelihood) score for Bernoulli trials
loglikBern <- function(probs, bool) {
  sum(dbinom(bool, 1, probs, log=TRUE), na.rm=TRUE)
}

# Setting up the parallel processing

library(parallel)
library(foreach)
library(doParallel)
detectCores()  # On my big desktop: 24
ncores <- 5  # 5 folds running in parallel assuming >5 cores.
  # If >= 15 cores, run 3 parallel chains per fold, otherwise do in serial.
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Running the models and estimating the scores

system.time(
# result <- foreach(fold = 1:5, .errorhandling='pass', .packages=c("jagsUI")) %dopar% {
result <- foreach(fold = 1:5, .errorhandling='pass', .packages=c("jagsUI")) %do% {

  TEST <- foldID == fold  # Just for readability
  TRAIN <- foldID != fold

  CVscores <- matrix(fold, nrow=length(modelPsi), ncol=5) # matrix of results for this fold
  rownames(CVscores) <- names(modelPsi)
  colnames(CVscores) <- c("fold", "score01", "Log_score", "AUC", "Brier")

  # Run all 5 models
  for(m in seq_along(modelPsi)) {

    inits <- function() list(Z = matrix(1, nSpecies, sum(TRAIN)),
         halfsigmaPsiCoef = runif(length(modelPsi[[m]]), 0, 5),
         halfsigmaDetCoef = runif(length(modelP[[m]]), 0, 5))

    jdata <- list(Y = Y[, TRAIN, ], nSurveys = nSurveys[TRAIN],
      Psi.cov = Psi.cov[TRAIN, modelPsi[[m]]], Det.cov = X[TRAIN, , modelP[[m]]])

    wanted <- c("psiCoef", "detCoef") # all we need for score calculation

    out <- jags(jdata, inits, wanted, "MSOM_KnownSps.jags", codaOnly=wanted,
        # n.chains=3, n.adapt=1000, n.iter= 5000, n.thin=5, parallel=TRUE)
        n.chains=3, n.adapt=100, n.iter= 500, parallel=TRUE) # only use parallel=TRUE if you have enough cores

    # Work through the MCMC chains to calculate psi, detect, and ucpd, also predictions for Z and Y.

    nIter <- length(out$sims.list$deviance)
    fPsi.cov = Psi.cov[TEST, modelPsi[[m]]] # covariate matrices for this fold and model
    fDet.cov = X[TEST, , modelP[[m]]]
    sc01 <- llk <- AUC <- br <- numeric(nIter)   # objects to hold the output
    psi <- Zpred <- matrix(NA, nSpecies, sum(TEST))
    detect <- true.detect <- Ypred <- ucpd <- array(NA, c(nSpecies, sum(TEST), maxSurveys))

    for(iter in 1:nIter) {
      psiCoef <- out$sims.list$psiCoef[iter,,]
      detCoef <- out$sims.list$detCoef[iter,,]
      # The code here mimics the JAGS code:
      for (i in 1:nSpecies) {
        # ecological model
        for (j in 1:sum(TEST)) {
          psi[i, j] <- plogis(sum(psiCoef[i, ] * fPsi.cov[j, ]))
          Zpred[i, j] <- rbinom(1, 1, psi[i, j])
          # observation model
          for (k in 1:nSurveys[TEST][j]) {
            detect[i, j, k] <- plogis(sum(detCoef[i, ] * fDet.cov[j, k, ]))
            true.detect[i, j, k] <- detect[i, j, k] * Zpred[i, j]
            Ypred[i, j, k] <- rbinom(1, 1, true.detect[i, j, k])
            ucpd[i, j, k] <- detect[i, j, k] * psi[i, j]
          }
        }
      }
      sc01[iter] <- score01(Ypred, Y[, TEST, ])
      llk[iter] <- loglikBern(ucpd, Y[, TEST, ])
      AUC[iter] <- auroc(ucpd, Y[, TEST, ])
      br[iter] <- brier(ucpd, Y[, TEST, ])
    } # end of MCMC chains

    CVscores[m, 'score01'] <- mean(sc01)
    CVscores[m, 'Log_score'] <- -2 * mean(llk)
    CVscores[m, 'AUC'] <- mean(AUC)
    CVscores[m, 'Brier'] <- mean(br)
  } # end "Run all 5 models"


  CVscores
}  )  # end of fold processing

result
tmp <- simplify2array(result)
apply(tmp, 1:2, mean) # Compare with Broms et al Table 2

