
# Code for the blog post at
#  https://mmeredith.net/blog/2017/UnderOverflow.htm

.Machine
?.Machine
noquote(unlist(format(.Machine)))

# Some example probabilities
# --------------------------

logit_p <- seq(-50, 50, by=10)
( p <- plogis(logit_p) )
 # [1] 1.928750e-22 4.248354e-18 9.357623e-14 2.061154e-09 4.539787e-05
 # [6] 5.000000e-01 9.999546e-01 1.000000e+00 1.000000e+00 1.000000e+00
# [11] 1.000000e+00

p == 1
# FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE TRUE TRUE

# Using log probabilities
# -----------------------

log(p)
 # [1] -5.000000e+01 -4.000000e+01 -3.000000e+01 -2.000000e+01 -1.000005e+01
 # [6] -6.931472e-01 -4.539890e-05 -2.061154e-09 -9.348078e-14 0.000000e+00
# [11] 0.000000e+00
# that doesn't work, the last 2 values are now exactly 0.

( log_p <- plogis(logit_p, log.p=TRUE) )
 # [1] -5.000000e+01 -4.000000e+01 -3.000000e+01 -2.000000e+01 -1.000005e+01
 # [6] -6.931472e-01 -4.539890e-05 -2.061154e-09 -9.357623e-14 -4.248354e-18
# [11] -1.928750e-22


# Adding probabilities
# --------------------

log_p <- -745:-760
exp(log_p) == 0
 # [1] FALSE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [12] TRUE TRUE TRUE TRUE TRUE

# A simple solution

exp(log_p - max(log_p))
 # [1] 1.000000e+00 3.678794e-01 1.353353e-01 4.978707e-02 1.831564e-02
 # [6] 6.737947e-03 2.478752e-03 9.118820e-04 3.354626e-04 1.234098e-04
# [11] 4.539993e-05 1.670170e-05 6.144212e-06 2.260329e-06 8.315287e-07
# [16] 3.059023e-07

sum(exp(log_p - max(log_p)))
# [1] 1.581977
log(sum(exp(log_p - max(log_p))))
# [1] 0.458675
log(sum(exp(log_p - max(log_p)))) + max(log_p)
# [1] -744.5413

# A better version

( p_max <- which.max(log_p) )  # find which value is the largest
# [1] 1
exp(log_p[-p_max] - max(log_p))
    # use negative index to exclude the largest (which is 1)
 # [1] 3.678794e-01 1.353353e-01 4.978707e-02 1.831564e-02 6.737947e-03
 # [6] 2.478752e-03 9.118820e-04 3.354626e-04 1.234098e-04 4.539993e-05
# [11] 1.670170e-05 6.144212e-06 2.260329e-06 8.315287e-07 3.059023e-07
sum(exp(log_p[-p_max] - max(log_p)))
# [1] 0.5819765
log1p(sum(exp(log_p[-p_max] - max(log_p))))  # log1p puts back the 1
# [1] 0.458675
log1p(sum(exp(log_p[-p_max] - max(log_p)))) + max(log_p)
# [1] -744.5413

logSumExp <- function(log_p) {
  p_max <- which.max(log_p)
  log1p(sum(exp(log_p[-p_max] - max(log_p)))) + max(log_p)
}

logSumExp(log_p)
# [1] -744.5413

# Subtracting probabilities
# -------------------------
# calculating 1 - p

logit_p <- seq(-50, 50, by=10)
( log_p <- plogis(logit_p, log.p=TRUE) )
 # [1] -5.000000e+01 -4.000000e+01 -3.000000e+01 -2.000000e+01 -1.000005e+01
 # [6] -6.931472e-01 -4.539890e-05 -2.061154e-09 -9.357623e-14 -4.248354e-18
# [11] -1.928750e-22

( log_1mp <- plogis( -logit_p, log.p=TRUE) )
 # [1] -1.928750e-22 -4.248354e-18 -9.357623e-14 -2.061154e-09 -4.539890e-05
 # [6] -6.931472e-01 -1.000005e+01 -2.000000e+01 -3.000000e+01 -4.000000e+01
# [11] -5.000000e+0

A <- log( -expm1(log_p))

B <- log1p( -exp(log_p))

 # cbind(logit_p, A, log_1mp, B)
      # logit_p             A       log_1mp             B
 # [1,]     -50  0.000000e+00 -1.928750e-22 -1.928750e-22
 # [2,]     -40  0.000000e+00 -4.248354e-18 -4.248354e-18
 # [3,]     -30 -9.359180e-14 -9.357623e-14 -9.357623e-14
 # [4,]     -20 -2.061154e-09 -2.061154e-09 -2.061154e-09
 # [5,]     -10 -4.539890e-05 -4.539890e-05 -4.539890e-05
 # [6,]       0 -6.931472e-01 -6.931472e-01 -6.931472e-01
 # [7,]      10 -1.000005e+01 -1.000005e+01 -1.000005e+01
 # [8,]      20 -2.000000e+01 -2.000000e+01 -2.000000e+01
 # [9,]      30 -3.000000e+01 -3.000000e+01 -2.999983e+01
# [10,]      40 -4.000000e+01 -4.000000e+01          -Inf
# [11,]      50 -5.000000e+01 -5.000000e+01          -Inf

log1mExp <- function(log_p)
  ifelse(log_p > -0.693, log(-expm1(log_p)), log1p(-exp(log_p)))

all.equal(log_1mp, log1mExp(log_p))
# [1] TRUE

# See more in this later post:
# https://mmeredith.net/blog/2017/UnderOverflow_2.htm
