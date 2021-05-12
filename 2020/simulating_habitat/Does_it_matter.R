
# R code for the blog post at
#  https://mmeredith.net/blog/2020/simulating_habitat.htm

# Does it matter?
# ===============

# Simulating occupancy with p = 1 and a site covariate, analysis with glm.

# Load real covariates
# --------------------
library(AHMbook)
data(willowWarbler)
gdd <- willowWarbler$cells$gdd

library(unmarked)
data(Switzerland)
str(ch <- Switzerland)

area <- read.csv("Adirondack_Wetlands.csv", comment="#")$ha
larea <- log(area)

# Set up the simulation
# ---------------------
nsites <- 200
b0 <- 0
b1 <- 1

nreps <- 1e5
coefA <- array(NA, c(nreps, 6, 2))
dimnames(coefA) <- list(NULL, c("normal", "uniform", "elevation", "forest",
    "logArea", "gdd"),
    c("intrcpt", "slope"))
x <- matrix(NA, nsites, 6)
system.time(
for(i in 1:nreps) {
  x[,1] <- standardize(rnorm(nsites))
  x[,2] <- standardize(runif(nsites))
  x[,3] <- standardize(sample(ch$elevation, size=nsites))
  x[,4] <- standardize(sample(ch$forest, size=nsites))
  x[,5] <- standardize(sample(larea, size=nsites))
  x[,6] <- standardize(sample(gdd, size=nsites))
  for(j in 1:6) {
    psi <- plogis(b0 + x[,j]*b1)
    z <- rbinom(nsites, 1, psi)
    fit <- glm(z ~ x[,j], family="binomial")
    coefA[i,j,] <- coef(fit)
  }
} ) # with 100 sites, 1e4 took 85 secs, 1e5 17 mins
save(coefA, file="coefA n=200.RData")

apply(coefA, 2:3, mean)
apply(coefA, 2:3, sd)

library(beeswarm)
df <- as.data.frame(coefA[,,2])
beeswarm(df[1:300, ], method='hex', cex=0.4, col='blue', las=1,
  xlab="Covariate", ylab="Estimated slope")
abline(h=1, col='red', lwd=2)

RMSE <- sqrt(colMeans((coefA[,,2] - 1)^2))
sort(RMSE)

hist(coefA[,,2])
# Propn with Absolute Error greater than
crit <- 0.25
bad <- colMeans(abs(coefA[,,2] - 1) > crit)
sort(bad)

tst <- standardize(runif(10000))
mean(abs(range(tst)))
hist(tst)
curve(dnorm(x), -3, 3, ylim=c(0, 1.2))
rect(-1.74,0,+1.74, 1/(2*1.74), border="red")
abline(h=0)
lines(density(standardize(ch$elevation)), col='brown')
lines(density(standardize(ch$forest)), col='green')
lines(density(standardize(larea)), col='blue')

# Run larea with different sample sizes
# -------------------------------------
nsites <- 21:26*10
b0 <- 0
b1 <- 1

nreps <- 1e5
coefL <- array(NA, c(nreps, 6, 2))
dimnames(coefL) <- list(NULL, nsites, c("intrcpt", "slope"))
system.time(
for(i in 1:nreps) {
  for(j in 1:6) {
    x <- standardize(sample(larea, size=nsites[j]))
    psi <- plogis(b0 + x*b1)
    z <- rbinom(nsites[j], 1, psi)
    fit <- glm(z ~ x, family="binomial")
    coefL[i,j,] <- coef(fit)
  }
} ) # with 100 sites, 1e4 took 85 secs, 1e5 15 mins
save(coefL, file="coefL.RData")

apply(coefL, 2:3, mean)
apply(coefL, 2:3, sd)

RMSE <- sqrt(colMeans((coefL[,,2] - 1)^2))
sort(RMSE)

hist(coefL[,,2])
# Propn with Absolute Error greater than
crit <- 0.2
bad <- colMeans(abs(coefL[,,2] - 1) > crit)
sort(bad)
