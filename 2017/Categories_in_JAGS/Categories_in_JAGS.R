
# Code from the blog post at
#  https://mmeredith.net/blog/2017/Categories_in_JAGS.htm

# Simulate some data
# ------------------
playerNames <- c("Andy", "Bob", "Chan", "Dora")
skill <- c(0.3, 0.6, 0.7, 0.8)
set.seed(2017)
player.char <- sample(playerNames, size=200, replace=TRUE, prob=c(1,3,2,2))
player <- factor(player.char)
summary(player)
# Andy  Bob Chan Dora 
  # 25   84   52   39 
success <- rbinom(200, 1, p=skill[as.integer(player)])
fluff <- rnorm(200)
tapply(success, player, mean) # Check: should be close to 'skill'
     # Andy       Bob      Chan      Dora 
# 0.2800000 0.5833333 0.7115385 0.8461538 

# Frequentist analysis with glm
# -----------------------------

glm(success ~ player + fluff, family='binomial')

# Call:  glm(formula = success ~ player + fluff, family = "binomial")

# Coefficients:
# (Intercept)      playerBob     playerChan     playerDora        fluff  
    # -0.9589         1.2927         1.8726         2.6850       0.1198  

# Degrees of Freedom: 199 Total (i.e. Null);  195 Residual
# Null Deviance:      263.6 
# Residual Deviance: 239.2        AIC: 249.2

# The null model - no covariates
# ------------------------------

modelText <- "
model {
  for(i in 1:length(success)) {
    success[i] ~ dbern(p)
  }

  # Priors
  p ~ dbeta(1,1)
}"
writeLines(modelText, con="modelNull.jags")

JAGSdata <- list(success=success)
wanted <- "p"

library(jagsUI)
outNull <- jags(JAGSdata, NULL, wanted, "modelNull.jags",
  n.chains=3, n.adapt=100, n.iter=2000)
outNull
# ...
            # mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# p          0.629 0.033   0.563   0.629   0.692    FALSE 1 1.000  3668
# ...
plot(outNull)

# One categorical covariate
# -------------------------

modelText <- "
model {
  for(i in 1:length(success)) {
    success[i] ~ dbern(p[player[i]])
  }

  # Priors
  for(i in 1:4) {
    p[i] ~ dbeta(1,1)
  }
}"
writeLines(modelText, con="modelPlayer.jags")

JAGSdata <- list(success=success, player=as.integer(player))
str(JAGSdata)
# List of 2
 # $ success: int [1:200] 0 1 0 0 1 1 0 1 1 1 ...
 # $ player : int [1:200] 1 3 3 2 4 4 2 3 3 2 ...
 
wanted <- "p"

outPlayer <- jags(JAGSdata, NULL, wanted, "modelPlayer.jags",
  n.chains=3, n.adapt=100, n.iter=2000)
outPlayer
# ...
            # mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# p[1]       0.299 0.087   0.144   0.294   0.487    FALSE 1 1.001  1627
# p[2]       0.582 0.053   0.479   0.582   0.684    FALSE 1 1.000  3484
# p[3]       0.705 0.061   0.578   0.708   0.817    FALSE 1 1.000  6000
# p[4]       0.830 0.058   0.706   0.835   0.927    FALSE 1 1.001  2855
# ...
# Compare with the proportion of successes for each player:
tapply(success, player, mean)
     # Andy       Bob      Chan      Dora 
# 0.2800000 0.5833333 0.7115385 0.8461538 

# More than one covariate
# -----------------------

# One intercept per category
# ''''''''''''''''''''''''''

modelText <- "
model {
  for(i in 1:length(success)) {
    logit(p[i]) <- bPlayer[player[i]] + bFluff * fluff[i]
    success[i] ~ dbern(p[i])
  }

  # Priors
  for(i in 1:4) {
    pPlayer[i] ~ dbeta(1,1)         # prior specified on probability scale...
    bPlayer[i] <- logit(pPlayer[i]) # ... then converted to logit scale.
  }
  bFluff ~ dunif(-5,5)
}"
writeLines(modelText, con="model1.jags")

JAGSdata <- list(success=success, player=as.integer(player), fluff=fluff)
str(JAGSdata)
# List of 3
 # $ success: int [1:200] 0 1 0 0 1 1 0 1 1 1 ...
 # $ player : int [1:200] 1 3 3 2 4 4 2 3 3 2 ...
 # $ fluff  : num [1:200] 0.48 0.688 -0.201 0.101 -0.577 ...
 
wanted <- c("bPlayer", "bFluff")

out1 <- jags(JAGSdata, NULL, wanted, "model1.jags",
  n.chains=3, n.adapt=100, n.iter=2000)
out1
# ...
              # mean    sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# bPlayer[1]  -1.013 0.448  -1.924  -0.998  -0.157    FALSE 0.989 1.002  5123
# bPlayer[2]   0.337 0.221  -0.095   0.338   0.764     TRUE 0.936 1.002  1195
# bPlayer[3]   0.943 0.316   0.345   0.931   1.588    FALSE 0.999 1.000  6000
# bPlayer[4]   1.835 0.466   1.011   1.811   2.827    FALSE 1.000 1.000  6000
# bFluff       0.123 0.166  -0.196   0.122   0.448     TRUE 0.771 1.000  6000
# ...

logitP1 <- out1$sims.list$bPlayer
P1 <- plogis(logitP1)
colMeans(P1)
# [1] 0.275 0.582 0.715 0.853

# One category as reference category
# ''''''''''''''''''''''''''''''''''

modelText <- "
model {
  for(i in 1:length(success)) {
    logit(p[i]) <- b0 + bPlayer[player[i]] + bFluff * fluff[i]
    success[i] ~ dbern(p[i])
  }

  # Priors
  b0 ~ dunif(-5,5)
  bPlayer[1] <- 0       # First level is reference category
  for(i in 2:4) {
    bPlayer[i] ~ dunif(-5,5) # Difference between each player and the reference player
  }
  bFluff ~ dunif(-5,5)
}"
writeLines(modelText, con="model2.jags")

wanted <- c("b0", "bPlayer", "bFluff")

out2 <- jags(JAGSdata, NULL, wanted, "model2.jags",
  n.chains=3, n.adapt=100, n.iter=2000)
out2
              # mean    sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# b0          -0.974 0.444  -1.891  -0.970  -0.117    FALSE 0.988 1.018   147
# bPlayer[1]   0.000 0.000   0.000   0.000   0.000    FALSE 1.000    NA     1
# bPlayer[2]   1.312 0.495   0.366   1.309   2.294    FALSE 0.996 1.014   172
# bPlayer[3]   1.916 0.540   0.881   1.904   3.017    FALSE 1.000 1.011   255
# bPlayer[4]   2.773 0.634   1.563   2.749   4.068    FALSE 1.000 1.008   276
# bFluff       0.124 0.165  -0.196   0.123   0.447     TRUE 0.773 1.001  2184

player_Bob <- relevel(player, ref="Bob")
summary(player_Bob)
# Bob Andy Chan Dora 
 # 84   25   52   39 
JAGSdata <- list(success=success, player=as.integer(player_Bob), fluff=fluff)

out2 <- jags(JAGSdata, NULL, wanted, "model2.jags",
  n.chains=3, n.adapt=100, n.iter=2000)
out2
            # mean    sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# b0         0.337 0.222  -0.078   0.329   0.783     TRUE 0.941 1.000  4153
# bPlay[1]   0.000 0.000   0.000   0.000   0.000    FALSE 1.000    NA     1
# bPlay[2]  -1.349 0.517  -2.431  -1.332  -0.387    FALSE 0.998 1.000  6000
# bPlay[3]   0.599 0.383  -0.151   0.593   1.364     TRUE 0.942 1.000  4634
# bPlay[4]   1.467 0.503   0.510   1.453   2.502    FALSE 0.999 1.001  2687
# bFluff     0.124 0.166  -0.199   0.123   0.450     TRUE 0.770 1.001  2864

logitP2 <- out2$sims.list$bPlayer + out2$sims.list$b0
P2 <- plogis(logitP2)
colMeans(P2)
# [1] 0.581 0.276 0.714 0.850

# Categorical coefficients sum to 0
# '''''''''''''''''''''''''''''''''

modelText <- "
model {
  for(i in 1:length(success)) {
    logit(p[i]) <- b0 + bPlayer[player[i]] + bFluff * fluff[i]
    success[i] ~ dbern(p[i])
  }

  # Priors
  b0 ~ dunif(-5,5)
  for(i in 1:4) {
    temp[i] ~ dnorm(0, 0.01)
  }
  bPlayer <- temp - mean(temp)
  bFluff ~ dunif(-5,5)
}"
writeLines(modelText, con="model3.jags")

out3 <- jags(JAGSdata, NULL, wanted, "model3.jags",
  n.chains=3, n.adapt=100, n.iter=2000)
out3
#               mean    sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# b0           0.524 0.192   0.141   0.523   0.903    FALSE 0.999 1.000  6000
# bPlayer[1]  -0.185 0.246  -0.679  -0.183   0.285     TRUE 0.771 1.000  6000
# bPlayer[2]  -1.532 0.375  -2.291  -1.522  -0.825    FALSE 1.000 1.001  2442
# bPlayer[3]   0.421 0.287  -0.141   0.417   0.996     TRUE 0.931 1.001  1380
# bPlayer[4]   1.296 0.378   0.603   1.282   2.064    FALSE 1.000 1.000  6000
# bFluff       0.126 0.168  -0.204   0.123   0.447     TRUE 0.777 1.000  6000


tst <- rowSums(out3$sims.list$bPlayer)
range(tst)
# [1] -1.065814e-14  1.065814e-14

logitP3 <- out3$sims.list$bPlayer + out3$sims.list$b0
P3 <- plogis(logitP3)
colMeans(P3)
# [1] 0.583 0.277 0.716 0.852




