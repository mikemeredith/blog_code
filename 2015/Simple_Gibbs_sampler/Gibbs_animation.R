
# Animation of Gibbs sampler with salamander data

# Code for the animated GIF on the web page at
#   https://mmeredith.net/blog/2015/Gibbs_sampler.htm

# The first section is the same as the script "Simple_Gibbs_sampler.R",
#   but with only 1000 iterations, and we don't remove the burn-in.

S <- 39    # Number of sites
K <- 5     # Number of visits to each site
xObs <- 18 # Number of sites where salamanders were detected
d <- 30    # total number of detections

priPsi <- c(1, 1)
priPi <- c(1, 1)

nIter <- 1000               # number of iterations
psi <- pi <- numeric(nIter) # objects to hold results
psi[1] <-  pi[1] <- 0.5     # starting values

for(i in 2:nIter) {
  # 1. Calculate prob(occupied | not detected):
  psi0 <- (psi[i-1] * (1 - pi[i-1])^K) / (psi[i-1] * (1 - pi[i-1])^K + (1 - psi[i-1]))
  # ...and draw number of additional sites occupied
  xAdd <- rbinom(1, S - xObs, psi0)
  x <- xObs + xAdd
  # 2a. Draw new psi from beta(occupied+prior, unoccupied+prior)
  psi[i] <- rbeta(1, x + priPsi[1], S - x + priPsi[2])
  # 2b. Draw new pi from beta(detected+prior, undetected+prior).
  pi[i] <- rbeta(1, d + priPi[1], x * K - d + priPi[2])
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
png(file = "gfx/Gibbs%03d.png", width=400, height=400) ; png.dev <- dev.cur()

# OR run these 2 lines to display the plots one by one:
# if(interactive())
  # devAskNewPage(TRUE)

( show <- c(2, 2, 2, 2, 3, 3, 3, 4, 4,
    round(seq(2, sqrt(1000), length=90)^2)) )
for(i in 1:3) {
  plot(0.5, 0.5, xlim=c(0.2, 1), ylim=c(0, 0.8), type='n', asp=1, las=1,
    xlab="occupancy, psi", ylab="detection, pi",
    main="Starting values")
  points(psi[1], pi[1], pch=16, col='red')
}
for(i in show) {
  plot(0.5, 0.5, xlim=c(0.2, 1), ylim=c(0, 0.8), type='n', asp=1, las=1,
    xlab="occupancy, psi", ylab="detection, pi",
    main=paste("Iteration #", i))
  points(psi[1], pi[1], pch=16, col='red')
  lines(psi[1:i], pi[1:i], col='grey')
  points(psi[2:i], pi[2:i], pch=16, cex=0.5)
}

dev.off(png.dev) # Do this to write plots for animation

# Make the animated GIF file

files <- list.files(pattern = ".png$", recursive = TRUE)
library(gifski)
gifski(files, "Gibbs_sampler.gif", width = 400, height = 400,
  delay = 0.25)
