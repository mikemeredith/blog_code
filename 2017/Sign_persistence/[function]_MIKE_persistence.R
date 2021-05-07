
# This is the function for sign persistence from the MIKE (Monitoring of
#   Illegal Killing of Elephants) program.

# See this blog post:
#  https://mmeredith.net/blog/2017/Sign_persistence.htm

MIKE.persistence <- function(DATA) {

  ##   Fit logistic regression model to STATE on DAYS, extract coefficients
  dung.glm <- glm(STATE ~ DAYS, data=DATA, family=binomial(link = "logit"))
  betas <- coefficients(dung.glm)
   
  ##   Calculate mean persistence time
  mean.decay <- -(1+exp(-betas[1])) * log(1+exp(betas[1])) / betas[2]

  ## Calculate the variance of the estimate
  vcovar <- vcov(dung.glm)
  var0 <- vcovar[1,1]  # variance of beta0
  var1 <- vcovar[2,2]  # variance of beta1
  covar <- vcovar[2,1] # covariance

  deriv0 <- -(1-exp(-betas[1]) * log(1+exp(betas[1])))/betas[2]
  deriv1 <- -mean.decay/betas[2]

  var.mean <- var0*deriv0^2 + 2*covar*deriv0*deriv1 + var1*deriv1^2

  ## Calculate the SE and CV and return
  se.mean <- sqrt(var.mean)
  cv.mean <- se.mean/mean.decay

  out <- c(mean.decay, se.mean, 100*cv.mean)
  names(out) <- c("Mean persistence time", "SE", "%CV")
  return(out)
}
