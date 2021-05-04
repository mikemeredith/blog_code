
# Animation of random walk sampler with salamander data

# Code for the animated GIF on the web page at
#   https://mmeredith.net/blog/2015/RandomWalk_MCMC.htm

# The first section is the same as the script "Simple_randomWalk_sampler.R",
#   but with only 1000 iterations, and we don't remove the burn-in.
# For the animation on the web page, I used a step size which is too small,
#   you can change that if you wish.

y <- c(2, 1, 1, 4, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
n <- 5
get_like <- function(psi, pi, y, n) {
  prod(psi * dbinom(y, n, pi) + ifelse(y==0, 1 - psi, 0))
}
priorPsi <- function(psi, shape1=1, shape2=1) dbeta(psi, shape1, shape2)
priorPi <- function(pi, shape1=1, shape2=1) dbeta(psi, shape1, shape2)
# ---------------------------
nIter <- 1000               # number of iterations
stepmaxPsi <- 0.1   # maximum step size
stepmaxPi <- 0.1
psi <- pi <- numeric(nIter) # objects to hold results
psi[1] <-  pi[1] <- 0.5     # starting values
likelihood <- get_like(psi[1], pi[1], y, n)  # initial likelihood
set.seed(123) # for reproducible results
for(i in 2:nIter) {
  # 1. Update psi:
   cand <- psi[i-1] + runif(1, -stepmaxPsi, stepmaxPsi)
  if (cand < 0 | cand > 1) {
    A <- 0
  } else {
    candLlh <- get_like(cand, pi[i-1], y, n)
    jointOld <- likelihood * priorPsi(psi[i-1])
    jointCand <- candLlh * priorPsi(cand)
    A <- min(1, jointCand / jointOld)
  }
  if(A > runif(1)) {   # if accepted
    psi[i] <- cand
    likelihood <- candLlh
  } else {
    psi[i] <- psi[i-1]
  }

  # 1. Update pi:
  cand <- pi[i-1] + runif(1,  -stepmaxPi, stepmaxPi)
  if (cand < 0 | cand > 1) {
    A <- 0
  } else {
    candLlh <- get_like(psi[i], cand, y, n)
    jointOld <- likelihood * priorPi(pi[i-1])
    jointCand <- candLlh * priorPi(cand)
    A <- min(1, jointCand / jointOld)
  }
  if(A > runif(1)) {   # if accepted
    pi[i] <- cand
    likelihood <- candLlh
  } else {
    pi[i] <- pi[i-1]
  }
}

plot(psi, pi, xlim=c(0.2, 1), ylim=c(0, 0.8), type='l', col='grey',
  xlab="occupancy, psi", ylab="detection, pi")
points(psi, pi, pch=16, cex=0.1)
points(psi[1], pi[1], pch=16, col='red')  # starting point


# Do stepwise plots (first 1000 iterations)
# -----------------------------------------
# Do you want to write the plots to PNG files, or look at each one?

# EITHER run the next 2 lines AND 'dev.off()' at the end to write to files
dir.create("gfx")  # sub-directory to store individual plot files
png(file = "gfx/randomWalk%03d.png", width=400, height=400) ; png.dev <- dev.cur()

# OR run these 2 lines to display the plots one by one:
# if(interactive())
  # devAskNewPage(TRUE)

( show <- c(2, 2, 2, 2, 3, 3, 3, 4, 4,
    round(seq(2, sqrt(1000), length=90)^2), rep(1000, 3)) )
stepText <- paste0("Maximum step size for psi = ",
    stepmaxPsi, ", for pi = ", stepmaxPi)
for(i in 1:3) {
  plot(0.5, 0.5, xlim=c(0.2, 1), ylim=c(0, 0.8), type='n', asp=1, las=1,
    xlab="occupancy, psi", ylab="detection, pi",
    main="Starting values")
  points(psi[1], pi[1], pch=16, col='red')
  mtext(stepText, 3, 0.5)
}
for(i in show) {
  plot(0.5, 0.5, xlim=c(0.2, 1), ylim=c(0, 0.8), type='n', asp=1, las=1,
    xlab="occupancy, psi", ylab="detection, pi",
    main=paste("Iteration #", i))
  points(psi[1], pi[1], pch=16, col='red')
  lines(psi[1:i], pi[1:i], col='grey')
  points(psi[2:i], pi[2:i], pch=16, cex=0.5)
  mtext(stepText, 3, 0.5)
}
dev.off(png.dev) # Do this to write plots for animation

# Make the animated GIF file

files <- list.files(pattern = ".png$", recursive = TRUE)
library(gifski)
gifski(files, "Random_walk.gif", width = 400, height = 400,
  delay = 0.25)
