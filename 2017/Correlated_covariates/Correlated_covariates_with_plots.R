
# Effect of rainfall and temperature on tree size

# Code for the blog post at
#  https://mmeredith.net/blog/2017/Collinearity.htm

rain <- c(2,1,1,0,0,-1,-1,-2,-2)
temp <- c(-2, -2, -1,-1,0,0,1,1,2)
dbh <- 2.5 + 1*temp + 0.5*rain
# Add is small error term, otherwise 'lm' gives strange results
set.seed(2017)
( dbh <- dbh + rnorm(9, 0, 0.02) )

# Plot the "truth"
( gr <- cor(temp, rain) )
symbols(temp, rain, circles = dbh, inches=0.4, lwd=2, bty='l', las=1,
  xlab="Temperature", ylab="Rainfall")
points(temp, rain, pch=19)
abline(mean(rain) - gr*mean(temp), gr, col="grey", lwd=4, lty=3) # add correlation
text(-0.013, 0.7, 'Correlation = -0.93', pos=4, cex=1.2)
arrows(-0.0126, 0.6954, -0.6962, 0.267, length=0.1)
# add explanations
arrows(-2:1 +0.1, 1:-2, -1:2-0.1, 1:-2, length=0.1, col='red', lwd=3)
text(0,-0.5, "hotter = bigger", col='red', cex=1.5, pos=4)
arrows(-2:1, 1:-2+0.1, -2:1, 2:-1-0.1, length=0.1, col='blue', lwd=3)
text(-2, 1.5, "wetter = bigger", col='blue', cex=1.5, pos=4)

# Fit a linear model with both covars
coef( tr <- lm(dbh ~ temp + rain) )

# Using just one covariate
# ========================

# Ignore rainfall
symbols(temp, rep(0,9), circles = dbh, inches=0.4, lwd=2,
  ylab="", yaxt='n', bty='n', xlab="Temperature")
points(temp, rep(0, 9), pch=19)
arrows(-2:1+0.1, rep(0,4), -1:2-0.1, rep(0,4),
  length=0.1, col='red', lwd=3)
text(0,0.4, "hotter = bigger", col='red', cex=1.5)

# Ignore temperature
symbols(rep(0,9), rain, circles = dbh, inches=0.4, lwd=2, las=1,
  xlab="", ylab="Rainfall", xaxt='n', bty='n')
points(rep(0, 9), rain, pch=19)
arrows(rep(0,4), -2:1+0.1, rep(0,4), -1:2-0.1,
  length=0.1, col='blue', lwd=3)
text(0.1,0, "wetter = \nsmaller?!", col='blue', cex=1.5, pos=4)

coef(t0 <- lm(dbh ~ temp))
coef(r0 <- lm(dbh ~ rain))

# Which model should we use to  make predictions?

AIC(tr, t0, r0)

# Why is collinearity a problem?
# ==============================

cov2cor(vcov(tr))

library(car)
vif(tr)
sqrt(vif(tr))

# Does this apply to Bayesian estimation too?
# ===========================================

library(jagsUI)
# JAGS model vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
writeLines("model {
for(i in 1:length(dbh)) {
  dbh[i] ~ dnorm(mu[i], tau)
  mu[i] <- b0 + bR*rain[i] + bT*temp[i]
  }
  b0 ~ dnorm(0, 0.01)
  bR ~ dnorm(0, 0.01)
  bT ~ dnorm(0, 0.01)
  sigma ~ dunif(0, 5)
  tau <- pow(sigma, -2)
} " , con="model.jags")
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
jagsData <- list(dbh=dbh, rain=rain, temp=temp)
wanted <- c("b0", "bT", "bR", "sigma")

( out <- jags(jagsData, NULL, wanted, "model.jags", DIC=FALSE,
  n.chains=3, n.adapt=100, n.iter=5000, n.burnin=0) )
plot(out)

with(out$sims.list, plot(bT, bR))
plot(out$sims.list$bT + out$sims.list$bR, out$sims.list$bT - out$sims.list$bR)

# Prediction for rain = 1.5, temp = 0.5
mu <- with(out$sims.list, b0 + bT * 0.5 + bR * 1.5)
hist(mu)
mean(mu)
sd(mu)
quantile(mu, c(0.025, 0.5, 0.975))

# Compare CVs:
( CV_bT <- out$sd$bT / out$mean$bT )
( CV_bR <- out$sd$bR / out$mean$bR )
( CV_mu <- sd(mu) / mean(mu) )
