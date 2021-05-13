
# R code for the blog post at
#  https://mmeredith.net/blog/2021/Posterior_mean_median.htm

set.seed(2021)
lpsiA <- rnorm(1e5, -4, 2)
lpsiB <- rnorm(1e5, -3, 0.7)

library(mcmcOutput)
par(mfrow=2:1)
postPlot(lpsiA, xlim=c(-10, 5), xlab="logit(psi)")
postPlot(lpsiB, xlim=c(-10, 5), xlab="logit(psi)")

# Convert to probability
psiA <- plogis(lpsiA)
psiB <- plogis(lpsiB)
postPlot(psiA, xlim=c(0, 0.4), xlab="psi")
postPlot(psiB, xlim=c(0, 0.4), xlab="psi")

# Using the median
postPlot(lpsiA, xlim=c(-10, 5), xlab="logit(psi)", center="median")
postPlot(lpsiB, xlim=c(-10, 5), xlab="logit(psi)", center="median") # not shown

postPlot(psiA, xlim=c(0, 0.4), xlab="psi", center="median")
postPlot(psiB, xlim=c(0, 0.4), xlab="psi", center="median")
