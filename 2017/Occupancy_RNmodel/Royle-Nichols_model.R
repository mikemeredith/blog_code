
# Code for the blog post at
#  https://mmeredith.net/blog/2017/Occupancy_RNmodel.htm

# Simulating data
# ---------------

# Parameters for simulated data:
nSite <- 200   # number of sites
nRep <- 10     # number of replicate surveys
lambda <- 1.5  # mean number of individuals at each site
r <- 0.1       # probability of detection of one individual

set.seed(2017)
# Biological process:
N <- rpois(nSite, lambda)  # Latent number available for detection
mean(N) # Mean N for this sample (not same as lambda)
# [1]  1.36
mean(N > 0) # True occupancy for this sample
# [1]  0.725

# Observation process:
p <- 1 - (1-r)^N
y <- rbinom(nSite, nRep, p) # number of occasions the species was detected at each site
head(y)
# [1]  3 1 2 2 2 2
mean(y > 0)  # naive estimate of occupancy
# [1] 0.555

# Maximum likelihood estimation
# -----------------------------

library(wiqid)
occSSrn0(y, n=nRep) # This needs a recent version of wiqid
# Call: occSSrn0(y = y, n = nRep)

# Real values:
             # est  lowCI  uppCI
# psiHat    0.6830 0.5755 0.7856
# lambdaHat 1.1488 0.8569 1.5401
# rHat      0.1143 0.0859 0.1506

# AIC: 1424.452

occSS0(y, n=nRep) # This also works with recent versions of wiqid
# Call: occSS0(y = y, n = nRep)

# Real values:
          # est  lowCI  uppCI
# psiHat 0.6247 0.5402 0.7023
# pHat   0.1969 0.1717 0.2247

# AIC: 1431

# The JAGS model
# --------------

cat(file="modelRN.jags", "
model {
  # Likelihood
  for (i in 1:nSite) {
    # Process model
    N[i] ~ dpois(lambda)
    z[i] <- step(N[i] - 1) # z=1 if N>0, ie. site is occupied

    # Observation model
    p[i] <- 1 - pow(1-r, N[i])
    y[i] ~ dbin(p[i], nRep)
  }
  # Priors
  lambda ~ dunif(0, 10)
  r ~ dunif(0, 1)

  # Derived quantities
  psi.sample <- mean(z[])   # Proportion of occupied sites in the sample
  psi <- 1 - exp( -lambda)  # Prob(N > 0 | lambda), probability that a random site is occupied.
}" )

# Running the model from R
# ------------------------

library(jagsUI)

# The data
str(JAGSdata <- list(y = y, nRep = nRep, nSite = nSite))
# List of 3
 # $ y    : int [1:200] 3 1 2 2 2 2 0 2 1 0 ...
 # $ n : num 10
 # $ nSite: num 200

# Starting values
Nst <- as.numeric(y > 0)  # N must be >0 if y > 0, so start at 1
inits <- function() list(N = Nst, lambda=runif(1,0,3))

# Parameters to monitor
params <- c("lambda", "r", "psi", "psi.sample", "N")

# Takes < 2 mins on my laptop:
out <- jags(JAGSdata, inits, params, "modelRN.jags", DIC=FALSE,
  n.chains = 3, n.adapt=1000, n.burnin=0, n.iter=1e5, n.thin=10, parallel=TRUE)
out
#             mean    sd  2.5%   50% 97.5% overlap0 f Rhat  n.eff
# lambda     1.183 0.184 0.887 1.161 1.602    FALSE 1    1  30000
# r          0.114 0.016 0.083 0.114 0.146    FALSE 1    1  22845
# psi        0.689 0.054 0.588 0.687 0.798    FALSE 1    1  30000
# psi.sample 0.689 0.041 0.620 0.685 0.780    FALSE 1    1  30000
# N[1]       2.096 0.922 1.000 2.000 4.000    FALSE 1    1  30000
# N[2]       1.385 0.635 1.000 1.000 3.000    FALSE 1    1  30000
# N[3]       1.705 0.808 1.000 2.000 4.000    FALSE 1    1  23728
# N[4]       1.697 0.811 1.000 2.000 4.000    FALSE 1    1  30000

Nhat <- unlist(out$mean)[-(1:4)]
plot(jitter(N), Nhat, ylim=range(0, Nhat))
abline(0, 1, col='red')
abline(v=mean(N), h=mean(Nhat), lty=3)

# Get mean(Nhat) for each value of N:
table(N)
# N
 # 0  1  2  3  4  6  7
# 55 70 40 23 10  1  1
tapply(Nhat, N, mean)
        # 0         1         2         3         4         6         7
# 0.3663258 1.0717541 1.6843473 1.9255339 2.2463930 2.1042233 4.2685067

