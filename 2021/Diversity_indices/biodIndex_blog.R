
# Fitting a multinomial distribution to data in JAGS
#  and calculating Shannon and Gini-Simpson indices.

# Blog version

# Look at data
library(wiqid)
?KillarneyBirds
data(KillarneyBirds)
KillarneyBirds
# number of birds in each sample
colSums(KillarneyBirds)
# Oak1   Oak2   Oak3    Yew  Sitka Norway  Mixed Patchy Swampy
#  170    182    112    110     75    198     91    119    100

# number of species recorded in each block of habitat:
colSums(KillarneyBirds > 0)
# Oak1   Oak2   Oak3    Yew  Sitka Norway  Mixed Patchy Swampy
#   20     22     15     15      8     14     17     21     18

apply(KillarneyBirds, 2, biodShannon)

# Focus on Norway spruce plantation
# ---------------------------------
y <- KillarneyBirds$Norway
( y <- y[y > 0] )
# [1] 30 30  3 65 20 11  4 14  2  9  3  5  1  1

p <- y / sum(y)
# Shannon index
-sum(p*log(p))
# 2.055991
# GiniSimpson
1 - sum(p^2)
# 0.8243036

# Fitting a binomial model to the Goldcrest data
# ----------------------------------------------

# The JAGS model:
cat(file="binomial.jags", "
model{
  # Prior
  p ~ dbeta(1,1)

  # Likelihood
  y ~ dbin(p, n)
}
")

jdata <- list(y=65, n=198)

wanted = "p"

library(jagsUI)
out_bin <- jags(jdata, inits=NULL, wanted, "binomial.jags", DIC=FALSE,
    n.chains=3, n.iter=1e4)
out_bin
#   mean    sd  2.5%   50% 97.5% overlap0 f Rhat n.eff
# p 0.33 0.033 0.268 0.329 0.396    FALSE 1    1 30000
mcmcOutput::postPlot(out_bin)

# Fitting the multinomial model
# -----------------------------

# The JAGS model:
cat(file="multinomial.jags", "
model{
  # Prior
  p ~ ddirch(alpha)

  # Model fit
  y ~ dmulti(p, nsps)

  # Calculate indices
  shannon <- -sum(p * log(p))
  GS <- 1 - sum(p^2)
}
")

nsps <- length(y)
jdata <- list(y=y, nsps=nsps, alpha=rep(1, nsps))

wanted = c("shannon", "GS", "p")

out_mult <- jags(jdata, inits=NULL, wanted, "multinomial.jags", DIC=FALSE,
    n.chains=3, n.iter=1e4)
out_mult
#          mean    sd  2.5%   50% 97.5% overlap0 f Rhat n.eff
# shannon 2.106 0.061 1.984 2.107 2.221    FALSE 1    1 17850
# GS      0.834 0.015 0.801 0.835 0.860    FALSE 1    1 30000
# p[1]    0.146 0.024 0.102 0.145 0.197    FALSE 1    1 26676
# p[2]    0.146 0.024 0.102 0.145 0.197    FALSE 1    1 30000
# ...
mcmcOutput::postPlot(out_mult, c("shannon", "GS"), layout=1:2)

# Comparison of sites with violin plot
# ------------------------------------

shannon_all <- matrix(NA, 3e4, 9)
colnames(shannon_all) <- colnames(KillarneyBirds)

for(i in 1:9) {
  y <- KillarneyBirds[, i]
  y <- y[y > 0]
  nsps <- length(y)
  jdata <- list(y=y, nsps=nsps, alpha=rep(1, nsps))
  out_tmp <- jags(jdata, inits=NULL, wanted, "multinomial.jags", DIC=FALSE,
      n.chains=3, n.iter=1e4)
  shannon_all[, i] <- out_tmp$sims.list$shannon
}

library(plotrix)
violin_plot(shannon_all, col="dodgerblue", main="Shannon index",
    x_axis_labels=colnames(shannon_all))

# Check the difference between Sitka spruce and Norway spruce
# -----------------------------------------------------------
diff <- shannon_all[, "Norway"] - shannon_all[, "Sitka"]
postPlot(diff, compVal=0, xlab="Norway - Sitka", main="Shannon index")
mean(diff > 0)

# Before and after
ynew <- c(15, 15, 3, 15, 10, 6, 3, 7, 5, 2, 2)
pnew <- ynew / sum(ynew)
# Shannon index
-sum(pnew*log(pnew))  # 2.170
# Gini-Simpson index
1 - sum(pnew^2)       # 0.868

nsps <- length(ynew)
jdata <- list(y=ynew, nsps=nsps, alpha=rep(1, nsps))
out_new <- jags(jdata, inits=NULL, wanted, "multinomial.jags", DIC=FALSE,
    n.chains=3, n.iter=1e4)
out_new

par(mfrow=1:2)
diff <- out_new$sims.list$shannon - out_mult$sims.list$shannon
postPlot(diff, compVal=0, xlab="new - old", main="Shannon index")

diffGS <- out_new$sims.list$GS - out_mult$sims.list$GS
postPlot(diffGS, compVal=0, xlab="new - old", main="Gini-Simpson index")
