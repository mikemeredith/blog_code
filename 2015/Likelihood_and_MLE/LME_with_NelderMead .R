
# Comparison of MLE and MCMC mechanisms

## Data set to work on: occupancy data
## ===================================

y <- c(2, 1, 1, 4, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
n <- 5

# Function to calculate likelihood:
get_like <- function(psi, pi, y, n) {
  prod(psi * dbinom(y, n, pi) + ifelse(y==0, 1 - psi, 0))
}
# Vectorize it so that it deals with a vector of psi's and pi's:
Get_like <- Vectorize(get_like, vectorize.args=c("psi", "pi"))


## Maximum Likelihood with simplex (Nelder-Mead) algorithm
## =======================================================
target <- 0.001 # Stop when improvement in < 0.1%
maxIter <- 30   #   or after 30 iterations.
Xpsi <- Xpi <- matrix(NA, maxIter, 3)  # Triangle coordinates stored in here
move <- numeric(maxIter)              # Record the type of step
Xpsi[1, ] <- c(0.48, 0.52, 0.48)  # Arbitrary starting values
Xpi[1, ]   <- c(0.48, 0.52, 0.52)    #     close to (0.5, 0.5)
for(i in 1:(maxIter-1)) {
  # 1. Calculate likelihoods, find which is worst (lowest llh):
  lh <- Get_like(Xpsi[i, ], Xpi[i, ], y=y, n=5)
  worst <- which.min(lh)
  # 2. Reflection: Calculate centroid of remaining vertices, find
  #      reflected point and its likelihood
  cent_psi <- mean(Xpsi[i, -worst])
  cent_p   <- mean(Xpi[i, -worst])
  new_psi <- cent_psi + (cent_psi - Xpsi[i, worst])
  new_p   <- cent_p   + (cent_p   - Xpi  [i, worst])
  ref_lh <- Get_like(new_psi, new_p, y=y, n=5)
  new_move <- 1
  # 3. If reflected point is best, try Expansion
  if (ref_lh > max(lh)) {
    exp_psi <- cent_psi + 2 * (cent_psi - Xpsi[i, worst])
    exp_p   <- cent_p   + 2 * (cent_p   - Xpi  [i, worst])
    if (Get_like(exp_psi, exp_p, y=y, n=5) > ref_lh) {
      new_move <- 2
      new_psi <- exp_psi
      new_p   <- exp_p
    }
  # 4. If reflected  point is worse than other vertices, try contraction
  } else if (ref_lh < min(lh[-worst])) {
    cont_psi <- cent_psi - 0.5 * (cent_psi - Xpsi[i, worst])
    cont_p   <- cent_p   - 0.5 * (cent_p   - Xpi  [i, worst])
    if (Get_like(cont_psi, cont_p, y=y, n=5) > lh[worst]) {
      new_move <- 3
      new_psi <- cont_psi
      new_p   <- cont_p
    }
  }
  # Store results for the next iteration
  Xpsi[i+1, ] <- Xpsi[i, ]
  Xpsi[i+1, worst] <- new_psi
  Xpi  [i+1, ] <- Xpi  [i, ]
  Xpi  [i+1, worst] <- new_p
  move[i+1] <- new_move
  # 5. Check if improvement < target
  improve <- (min(Get_like(Xpsi[i+1, ], Xpi[i+1, ], y=y, n=5)) - min(lh)) / min(lh)
  if(improve < target)
    break
}
( actualIter <- i + 1)
mean(Xpsi[actualIter, ])
mean(Xpi[actualIter, ])

# Do a single plot with final result
# ----------------------------------
par(las=1)
PSI <- PI <- seq(0.0025, 1, 0.005)
LH <- outer(PSI, PI, Get_like, y=y, n=n)

contour(LH, drawlabels=FALSE, xlim=c(0.4, 1), ylim=c(0, 0.6), asp=1,
  col='grey', xlab="occupancy, psi", ylab="detection, pi")
polygon(Xpsi[1, ], Xpi[1, ], col='grey')
for(j in 2:actualIter) {
  polygon(Xpsi[j, ], Xpi[j, ], col=c('skyblue', 'yellow', 'pink')[move[j]])
}
abline(v=mean(Xpsi[actualIter, ]), h=mean(Xpi[actualIter, ]), col='grey')
legend('topright', c("Reflect", "Expand", "Contract"), pch=24, 
    pt.bg=c('skyblue', 'yellow', 'pink'), pt.cex=2, bty='n' )

# Do stepwise plots
# -----------------
# Do you want to write the plots to PNG files, or look at each one?

# Run the next 2 lines AND 'dev.off()' at the end to write to files
dir.create("gfx")  # sub-directory to store individual plot files
png(file = "gfx/NM%03d.png", width=400, height=400) ; png.dev <- dev.cur()

# Run this line to display the plots one by one:
# if(interactive())
  # devAskNewPage(TRUE)

contour(LH, drawlabels=FALSE, xlim=c(0.4, 1), ylim=c(0, 0.6),
  col='grey', main="Starting values", 
  xlab="occupancy, psi", ylab="detection, pi")
polygon(Xpsi[1, ], Xpi[1, ], col='grey')
for(i in 2:actualIter) {
  contour(LH, drawlabels=FALSE, xlim=c(0.4, 1), ylim=c(0, 0.6),
    col='grey', main=paste("Iteration", i),
    xlab="occupancy, psi", ylab="detection, pi")
  polygon(Xpsi[1, ], Xpi[1, ], col='grey')
  legend('topright', c("Reflect", "Expand", "Contract"),
    pch=24, pt.bg=c('blue', 'yellow', 'red'), pt.cex=2, bty='n')
  for(j in 2:i)
    polygon(Xpsi[j, ], Xpi[j, ], col=c('blue', 'yellow', 'red')[move[j]])
  mtext(c("Reflect", "Expand", "Contract")[move[i]], 1, 2)
}
contour(LH, drawlabels=FALSE, xlim=c(0.4, 1), ylim=c(0, 0.6),
  col='grey', main="Final result",
  xlab="occupancy, psi", ylab="detection, pi")
polygon(Xpsi[1, ], Xpi[1, ], col='grey')
legend('topright', c("Reflect", "Expand", "Contract"),
  pch=24, pt.bg=c('blue', 'yellow', 'red'), pt.cex=2, bty='n')
for(j in 2:i)
  polygon(Xpsi[j, ], Xpi[j, ], col=c('blue', 'yellow', 'red')[move[j]])
abline(v=mean(Xpsi[actualIter, ]), h=mean(Xpi[actualIter, ]), col='grey')

dev.off(png.dev) # Do this to write plots for animation

files <- list.files(pattern = ".png$", recursive = TRUE)

# For the .mp4 file:
library(av)
av_encode_video(files, "NelderMead.mp4", framerate=2)

# For the animated GIF file
library(gifski)
gifski(files, "NelderMead.gif", width = 400, height = 400,
  delay = 0.5)
