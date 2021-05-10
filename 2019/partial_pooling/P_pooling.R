
# R and JAGS code for the blog post at
#  https://mmeredith.net/blog/2019/partial_pooling.htm

library(jagsUI)

n <- c(60, 60, 60, 60, 7, 6, 4, 3, 0)
y <- c(42, 38, 21, 19, 3, 5, 0, 2, 0)
names(n) <- LETTERS[1:9]
cbind(n=n, y=y, R=round(y/n, 2))

sum(y)/sum(n)

# No pooling
modeltxt <- "
model{
  for(i in 1:8) {
    y[i] ~ dbin(R[i], n[i])
    R[i] ~ dbeta(1, 1)
  }
}"
writeLines(modeltxt, con="nopooling.jags")

jdata <- list(y=y, n=n)
wanted <- "R"

( outnp <- jags(jdata, NULL, wanted, "nopooling.jags", DIC=FALSE,
  n.burn=100, n.iter=1100, n.chains=3) )
wiqid::diagPlot(outnp)

# Total pooling
modeltxt <- "
model{
  for(i in 1:9) {
    y[i] ~ dbin(R, n[i])
  }
  R ~ dbeta(1, 1)
}"
writeLines(modeltxt, con="totpooling.jags")

# jdata <- list(y=y, n=n)
# wanted <- "R"

( outtp <- jags(jdata, NULL, wanted, "totpooling.jags", DIC=FALSE,
  n.burn=100, n.iter=1100, n.chains=3) )
wiqid::diagPlot(outtp)

# Partial pooling
modeltxt <- "
model{
  for(i in 1:9) {
    # observation model
    y[i] ~ dbin(R[i], n[i])
    # prior for R[i]
    logitR[i] ~ dnorm(mu, tau)
    R[i] <- ilogit(logitR[i])
  }
  # hyperpriors for mu and tau
  Rmean ~ dbeta(1, 1)
  mu <- logit(Rmean)
  sd ~ dunif(0, 10)
  tau <- sd^-2
}"
writeLines(modeltxt, con="parpooling.jags")

wanted <- c("R", "Rmean", "sd")

( outpp <- jags(jdata, NULL, wanted, "parpooling.jags", DIC=FALSE,
  n.burn=100, n.iter=10100, n.chains=3) )
wiqid::diagPlot(outpp)

res <- cbind(n=n,
             y=y,
             naive = y/n,
             np = c(outnp$mean$R, NA),
             tp = rep(outtp$mean$R, 9),
             pp = outpp$mean$R)
round(res, 2)

# Plot all that

plot(0.5,5, type='n', xlim=0:1, ylim=c(1,11), yaxt='n', xlab="RAI", ylab="")
axis(2, at=1:11, labels=c("Rmean", LETTERS[9:1], "Rpooled"), las=1)
segments(x0=outpp$q2.5$R, y0=10:2-0.2, x1=outpp$q97.5$R, y1=10:2-0.2, lwd=2)
points(outpp$mean$R, 10:2-0.2, pch=16, cex=1.5)
segments(x0=outnp$q2.5$R, y0=10:3, x1=outnp$q97.5$R, y1=10:3, col='blue', lwd=2)
points(outnp$mean$R, 10:3, col='blue', pch=17, cex=1.5)
segments(x0=outtp$q2.5$R, y0=11, x1=outtp$q97.5$R, y1=11, col='red', lwd=2)
points(outtp$mean$R, 11, col='red', pch=4, cex=1.5)
segments(x0=outpp$q2.5$Rmean, y0=1, x1=outpp$q97.5$Rmean, y1=1, lwd=2)
points(outpp$mean$Rmean, 1, pch=1, cex=1.5)
points(y/n, 10:2+0.2, pch=6, cex=1.5)
legend('topleft', c("naive","no pooling", "total pooling", "partial pooling"),
  pch=c(6, 17, 4, 16), col=c('black','blue', 'red','black'),
  lwd=c(NA, 2,2,2), bty='n')


