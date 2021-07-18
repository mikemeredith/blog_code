
# What do implied priors for data augmentation look like?
# Simulations of discrete uniform prior, Bill Link's "scale prior",
#   and informative priors.

library(jagsUI)

cat(file="DAprior.jags", "
model{
  omega ~ dbeta(DAprior[1], DAprior[2])
  for(i in 1:M){
    w[i] ~ dbern(omega)
  }
  N <- sum(w)
} ")


# The discrete uniform prior
# ==========================

jdata <- list(M=50, DAprior=c(1,1))

out1 <- jags(jdata, NULL, c("N", "omega"), "DAprior.jags", DIC=FALSE,
  n.iter=11e6, n.burn=1e6, n.chain=3, n.thin=1, parallel=TRUE)
diagPlot(out1)

Npost <- out1$sims.list$N
Ndist <- table(Npost)/length(Npost)  # takes several seconds
windows(5,5)
plot(Ndist, ylab="Prior probability", xlab="N",
    main="Discrete uniform prior")

# Bill Link's "scale" prior
# =========================

jdata2 <- list(M=50, DAprior=c(0.001,1))
out2 <- jags(jdata2, NULL, c("N", "omega"), "DAprior.jags", DIC=FALSE,
  n.iter=11e6, n.burn=1e6, n.chain=3, n.thin=1, parallel=TRUE)
diagPlot(out2)

Npost <- out2$sims.list$N
mean(Npost==0)
mean(Npost==1)
Ndist <- table(Npost)/length(Npost)
plot(Ndist, ylab="Prior probability", xlab="N",
    main="Link's scale prior")
plot(Ndist, ylab="Prior probability", xlab="N", ylim=c(0,0.001),
    main="Link's scale prior")
arrows(3,0.00101,0,0.00101, length=0.1, xpd=TRUE)
text(3, 0.00101, "Probability of N=0 is 0.995", pos=4)

# Link's prior plus 1 animal detected
# -----------------------------------
jdata2a <- list(M=50, DAprior=c(0.001,1), w=c(1, rep(NA,49)))
out2a <- jags(jdata2a, NULL, c("N", "omega"), "DAprior.jags", DIC=FALSE,
  n.iter=11e6, n.burn=1e6, n.chain=3, n.thin=1, parallel=TRUE)

hist(out2a$sims.list$omega, freq=FALSE)

Npost <- out2a$sims.list$N
mean(Npost==0)
mean(Npost==1)
Ndist <- table(Npost)/length(Npost)
plot(Ndist, ylab="Prior probability", xlab="N",
    main="Scale prior, 1 animal present")

# Link's prior plus 10 animals detected
# -------------------------------------
jdata2b <- list(M=50, DAprior=c(0.001,1), w=c(rep(1,10), rep(NA,40)))
out2b <- jags(jdata2b, NULL, c("N", "omega"), "DAprior.jags", DIC=FALSE,
  n.iter=11e6, n.burn=1e6, n.chain=3, n.thin=1, parallel=TRUE)

hist(out2b$sims.list$omega, freq=FALSE)

Npost <- out2b$sims.list$N
mean(Npost==0)
mean(Npost==1)
Ndist <- table(Npost)/length(Npost)
plot(Ndist, ylab="Prior probability", xlab="N",
    main="Scale prior, 10 animals present")


# Informative prior
# =================

library(wiqid)
jdata <- list(M=50, DAprior=c(getBeta3Par(mode=20/50, concentration=3)))
# Change the value for concentration and rerun to do the three plots

out3 <- jags(jdata, NULL, c("N", "omega"), "DAprior.jags", DIC=FALSE,
  n.iter=11e6, n.burn=1e6, n.chain=3, n.thin=1, parallel=TRUE)
Npost <- out3$sims.list$N

Ndist_3 <- table(Npost)/length(Npost)
Ndist_5 <- table(Npost)/length(Npost)
Ndist_10 <- table(Npost)/length(Npost)

windows(9,3)
par(mfrow=c(1,3))
plot(Ndist_3, ylab="Prior probability", xlab="N",
    main="Concentration = 3")
plot(Ndist_5, ylab="Prior probability", xlab="N",
    main="Concentration = 5")
plot(Ndist_10, ylab="Prior probability", xlab="N",
    main="Concentration = 10")
