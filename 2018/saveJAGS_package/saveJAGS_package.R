
# Code for the blog post at
#  https://mmeredith.net/blog/2018/Intro_saveJAGS.htm

library(saveJAGS)

# The classic salamander data: the number of occasions (out of 5)
#  that salamanders were detected at each of 39 sites
sal <- rep(0:4, c(21,12,1,4,1))

# JAGS code for psi(.) p(.) occupancy model
modelText <- "
model {
  for(i in 1:nSites) {
    z[i] ~ dbern(psi)
    y[i] ~ dbin(p * z[i], n)
  }
  psi ~ dbeta(1, 1)
  p ~ dbeta(1, 1)
} "
writeLines(modelText, con = "JAGSmodel.jags")

JAGSdata <- list(y = sal, n = 5, nSites = length(sal), z=ifelse(sal > 0, 1, NA))
inits <- function(chain) list(a=NULL)
wanted <- c("p", "psi", "z")

# That's all the same; here's what's new:

dir.create("mySaves")

res1 <- saveJAGS(JAGSdata, , wanted, "JAGSmodel.jags",
        chains=3, sample2save=1000, nSaves=4, burnin=1000, thin=10,
        fileStub="mySaves/testing")

str(res1)
# List of 3
 # $ AA: chr [1:4] "mySaves/test_AA_001_210507_1432.RData" ...
 # $ AB: chr [1:4] "mySaves/test_AB_001_210507_1432.RData" ..
 # $ AC: chr [1:4] "mySaves/test_AC_001_210507_1432.RData" ...
 # - attr(*, "class")= chr [1:2] "saveJAGSfileList" "list"

summary(res1)
# File list with 3 chains, each with 4 files.
# Median run time per file less than 1 min.
# Adaptation was adequate.
# Chains already thinned by: 10
# Draws saved: 1000 per file, 4000 per chain, 12000 total.
# Total values saved: 492 thousand; expected object size: 7.87 MB
# 3 parameters with 41 nodes monitored.
#     nodes
# p       1
# psi     1
# z      39

# If the saveJAGS run was terminated before completion,
#   we can recover the file list with recoverSaves:
res2 <- recoverSaves("mySaves/testing")
all.equal(res1, res2)

# Load the results into R as an mcmc.list object
mcmc1 <- combineSaves(res2)
str(mcmc1)
object.size(mcmc1)  # only 4 MB

#... or subset by parameters, or thin:
mcmc2 <- combineSaves(res1, params=c("psi","p"), thin=4)
str(mcmc2)
object.size(mcmc2)

# Load the results into R as an mcmcOutput object
( mco1 <- mcmcOutput(res1) )

#... or subset by parameters, or thin:
( mco2 <- mcmcOutput(res1, params=c("psi","p"), thin=4) )

diagPlot(mco2)      # diagnostic plots
plot(mco2)          # plot posterior distributions
summary(mco2)       # display a summary in the Console
View(summary(mco1)  # view the summary in spreadsheet format

p0 <- mco2$p
summary(p0)

# Resuming the run
# ----------------

newRes <- resumeJAGS(fileStub="mySaves/testing", nSaves=3)
summary(newRes)
