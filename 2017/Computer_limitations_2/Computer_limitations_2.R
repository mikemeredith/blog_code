
# Code for the blog post at
#  https://mmeredith.net/blog/2017/UnderOverflow_2.htm
# which is a follow up to
#  https://mmeredith.net/blog/2017/UnderOverflow.htm


# Summing probabilities
# ---------------------

logSumExp <- function(log_p) {
  p_max <- which.max(log_p)
  log1p(sum(exp(log_p[-p_max] - max(log_p)))) + max(log_p)
}

log_p <- -745:-760
exp(log_p) == 0
 # [1] FALSE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [11] TRUE TRUE TRUE TRUE TRUE TRUE
logSumExp(log_p)
# [1] -744.5413

p0 <- rep(0, 5)
( log_p0 <- log(p0) )
# [1] -Inf -Inf -Inf -Inf -Inf
logSumExp(log_p0)
# [1] NaN

logSumExp <- function(log_p) {
  if(max(log_p) == -Inf)
    return(-Inf)
  p_max <- which.max(log_p)
  log1p(sum(exp(log_p[-p_max] - max(log_p)))) + max(log_p)
}

logSumExp(log_p0)
# [1] -Inf

# Adding two vectors
# ------------------

( p1 <- (1:9) / 10 )
# [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
( p2 <- runif(9) * (1 - p1) )
# [1] 0.13170373 0.27850206 0.57787359 0.29211162 0.08092575
# [6] 0.15465695 0.22973600 0.16997316 0.05788608
( add <- p1 + p2 )
# [1] 0.2317037 0.4785021 0.8778736 0.6921116 0.5809257
# [6] 0.7546570 0.9297360 0.9699732 0.9578861
log_p1 <- log(p1)
log_p2 <- log(p2)

logAddExp <- function(logp1, logp2) {
  bigger <- pmax(logp1, logp2)
  smaller <- pmin(logp1, logp2)
  log1p(exp(smaller - bigger)) + bigger
}

( logadd <- logAddExp(log_p1, log_p2) )
# [1] -1.46229573 -0.73709475 -0.13025268 -0.36800804 -0.54313233
# [6] -0.28149200 -0.07285460 -0.03048688 -0.04302642
all.equal(logadd, log(add))
# [1] TRUE

p1[3] <- 0
p2[3] <- 0
log_p1 <- log(p1)
log_p2 <- log(p2)
( logadd <- logAddExp(log_p1, log_p2) )
# [1] -1.46229573 -0.73709475         NaN -0.36800804 -0.54313233
# [6] -0.28149200 -0.07285460 -0.03048688 -0.04302642
all.equal(logadd, log(add))
# [1] "'is.NA' value mismatch: 0 in current 1 in target"

logAddExp <- function(logp1, logp2) {
  bigger <- pmax(logp1, logp2)
  smaller <- pmin(logp1, logp2)
  fix <- smaller > -Inf
  bigger[fix] <- log1p(exp(smaller[fix] - bigger[fix])) + bigger[fix]
  return(bigger)
}

( logadd <- logAddExp(log_p1, log_p2) )
# [1] -1.46229573 -0.73709475        -Inf -0.36800804 -0.54313233
# [6] -0.28149200 -0.07285460 -0.03048688 -0.04302642
all.equal(logadd, log(p1 + p2))
# [1] TRUE

# Matrix multiplication
# ---------------------

# Example matrices (not probabilities)
(A <- matrix(1:6, 3, 2))
     # [,1] [,2]
# [1,]    1    4
# [2,]    2    5
# [3,]    3    6

(B <- matrix(11:18, 2, 4, byrow=TRUE))
     # [,1] [,2] [,3] [,4]
# [1,]   11   12   13   14
# [2,]   15   16   17   18

A %*% B
     # [,1] [,2] [,3] [,4]
# [1,]   71   76   81   86
# [2,]   97  104  111  118
# [3,]  123  132  141  150

X <- matrix(NA, nrow(A), ncol(B))
for(i in 1:nrow(X))
  for(j in 1:ncol(X))
    X[i,j] <- sum(A[i, ] * B[, j])
all.equal(X, A %*% B)
# [1] TRUE

logMatMultExp <- function(logA, logB) {
  logX <- matrix(NA_real_, nrow(logA), ncol(logB))
  for(i in 1:nrow(logX))
    for(j in 1:ncol(logX))
      logX[i,j] <- logSumExp(logA[i, ] + logB[, j])
  return(logX)
}

( logA <- log(A) )
          # [,1]     [,2]
# [1,] 0.0000000 1.386294
# [2,] 0.6931472 1.609438
# [3,] 1.0986123 1.791759

( logB <- log(B) )
         # [,1]     [,2]     [,3]     [,4]
# [1,] 2.397895 2.484907 2.564949 2.639057
# [2,] 2.708050 2.772589 2.833213 2.890372

log(A %*% B)
         # [,1]     [,2]     [,3]     [,4]
# [1,] 4.262680 4.330733 4.394449 4.454347
# [2,] 4.574711 4.644391 4.709530 4.770685
# [3,] 4.812184 4.882802 4.948760 5.010635

all.equal(logMatMultExp(logA, logB), log(A %*% B))
# [1] TRUE

A <- matrix(sample(seq(0, 1, length=6)), 3, 2)
( logA <- log(A) )
           # [,1]       [,2]
# [1,] -1.6094379 -0.9162907
# [2,] -0.5108256 -0.2231436
# [3,]  0.0000000       -Inf

B <- matrix(sample(seq(0, 1, length=8)), 2, 4)
( logB <- log(B) )
           # [,1]       [,2]       [,3] [,4]
# [1,] -0.5596158 -0.8472979 -0.3364722    0
# [2,] -0.1541507 -1.2527630 -1.9459101 -Inf

log(A %*% B)
            # [,1]       [,2]       [,3]       [,4]
# [1,] -0.78275934 -1.6094379 -1.6094379 -1.6094379
# [2,]  0.02817088 -0.7221347 -0.6109091 -0.5108256
# [3,] -0.55961579 -0.8472979 -0.3364722  0.0000000

all.equal(logMatMultExp(logA, logB), log(A %*% B))
# [1] TRUE

A <- matrix(0, 3, 2)
( logA <- log(A) )
     # [,1] [,2]
# [1,] -Inf -Inf
# [2,] -Inf -Inf
# [3,] -Inf -Inf

A %*% B
     # [,1] [,2] [,3] [,4]
# [1,]    0    0    0    0
# [2,]    0    0    0    0
# [3,]    0    0    0    0

all.equal(logMatMultExp(logA, logB), log(A %*% B))
# [1] TRUE

# Final touches
# -------------

logMatMultExp <- function(logA, logB) {
  if(ncol(logA) != nrow(logB))
    stop("non-conformable matrices")
  logX <- matrix(NA_real_, nrow(logA), ncol(logB))
  for(i in 1:nrow(logX))
    for(j in 1:ncol(logX))
      logX[i,j] <- logSumExp(logA[i, ] + logB[, j])
  return(logX)
}

logp1 <- plogis(seq(-50, 50, length=6), log.p=TRUE)
logA <- matrix(logp1, 3, 2)
logp2 <- plogis(sample(seq(-50, 50, length=8)), log.p=TRUE)
logB <- matrix(logp2, 2, 4)
logMatMultExp(logA, logB)
              # [,1]          [,2]      [,3]          [,4]
# [1,] -4.539939e-05 -0.0008355769 -49.99926 -4.539890e-05
# [2,] -4.940512e-10 -0.0007901781 -37.14364 -9.388489e-14
# [3,] -4.939352e-10 -0.0007447453 -17.14369 -3.086479e-16
