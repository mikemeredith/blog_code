
# Code for the web page on MLE for salamanders.

y <- c(4, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
n <- 5
sum(y > 0)

PSI <- PI <- seq(0.0025, 1, 0.005)
range(PSI)
length(PSI)

LIKE <- matrix(NA, 200, 200)
for (i in 1:200) {
  for(j in 1:200) {
    LIKE[i, j] <- prod(PSI[i] * dbinom(y, n, PI[j]) + ifelse(y==0, 1 - PSI[i], 0))
  }
}
range(LIKE)

image(LIKE, xlab="occupancy, psi", ylab="detection, pi", main="Likelihood surface",
    col = terrain.colors(12), las=1, asp=1)
contour(LIKE, add=TRUE, drawlabels=FALSE)

# This is the persective plot at the top of the page
# I rotate.persp to see what rotation would be best
# library(TeachingDemos)
# rotate.persp(PSI, PI, LIKE)

persp(PSI, PI, LIKE, col='skyblue', expand=0.5,
  theta=-35, phi=25, shade=0.1,
  xlab="occupancy, psi", ylab="detection, pi",
  zlab="likelihood") #, main="Likelihood surface")

# What is the MLE of psi and pi?
( where <- arrayInd(which.max(LIKE), dim(LIKE)) )
PSI[where[1]]   # MLE of psi
# [1] 0.595
PI[where[2]]     # MLE of pi
# [1] 0.26
abline(v = PSI[where[1]], h = PI[where[2]], col='white')

# Do MLE with optim
get_negLogLike <- function(x, y, n) {
  psi <- x[1]
  pi <- x[2]
  like <- prod(psi * dbinom(y, n, pi) + ifelse(y==0, 1 - psi, 0))
  return(-log(like))
}
get_negLogLike(c(0.7, 0.3), y=y, n=5)  # check
# 49.8396

optim(c(0.5,0.5), get_negLogLike, y=y, n=5)

