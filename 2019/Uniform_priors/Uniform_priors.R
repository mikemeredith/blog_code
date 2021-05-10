
# R code for the blog post at
#  https://mmeredith.net/blog/2019/Uniform_priors.htm

# Uniform priors for logit models

library(jagsUI)
library(mcmcOutput)

set.seed(123)
x1 <- rep(1:0, c(40,60))
x2 <- rnorm(100)
y = x1   # separation
table(y)
glm(y ~ x1 + x2, family="binomial")

cat(file="logisUniform.jags", "
# uniform priors
model{
  for(i in 1:100) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1 * x1[i] # + b2 * x2[i]
  }
  b0 ~ dunif(-5, 5)
  b1 ~ dunif(-5, 5)
  b2 ~ dunif(-5, 5)
} ")

cat(file="logisNormal.jags", "
# normal priors
model{
  for(i in 1:100) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1 * x1[i] # + b2 * x2[i]
  }
  b0 ~ dnorm(0, 0.1)
  b1 ~ dnorm(0, 0.1)
  b2 ~ dnorm(0, 0.1)
} ")

bdata <- list(y = y, x1 = x1, x2 = x2)

wanted <- c("b0", "b1", "b2")

# First plot: separation with normal priors

out1 <- jags(bdata,,wanted, "logisNormal.jags", DIC=FALSE,
  n.chains=3, n.iter=5100, n.burn=100, n.thin=1)
out1
diagPlot(out1, main="Prior: dnorm(0, 0.1)")

# Second plot: separation with uniform priors

out2 <- jags(bdata,,wanted, "logisUniform.jags", DIC=FALSE,
  n.chains=3, n.iter=5100, n.burn=100, n.thin=10)
out2
diagPlot(out2, main="Prior: dunif(-5, 5)")

# Response variable is missing (or mis-coded)

bdata <- list(Y = y, x1 = x1, x2 = x2)

# Third plot: no data, normal priors

out3 <- jags(bdata,,wanted, "logisNormal.jags", DIC=FALSE,
  n.chains=3, n.iter=1100, n.burn=100)
out3
diagPlot(out3, main="Prior: dnorm(0, 0.1)")

# Fourth plot: no data, uniform priors

out4 <- jags(bdata,,wanted, "logisUniform.jags", DIC=FALSE,
  n.chains=3, n.iter=1100, n.burn=100)
out4
diagPlot(out4, main="Prior: dunif(-5, 5)")
