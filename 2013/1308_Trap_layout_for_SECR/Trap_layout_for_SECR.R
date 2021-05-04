
# This is the code for the plots at
# https://mmeredith.net/blog/2013/1308_Trap_layout_for_SECR.htm

# We compare a single 16 x 6 layout of traps with 3x3 clusters in a systematic
#   random layout.

# Load the packages we need (install these if necessary):
library(secr)
library(spatstat)
library(sp)
library(beeswarm)

# Load the polygon for the study block
# ====================================
# an irregular polygon with area 2000 km sq.
load("blockX.Rda")
ls()
class(win)
class(spdf)
plot(win, main="")
( area <- area.owin(win) )   # sq metres

## The tiger population
## ====================
# Small population, 5 animals in the block
n <- 5
D <- n / area * 1e4 # density per ha.
sigma <- 2500
g0 <- 0.02
territory <- 2 * sigma  # Minimum spacing of activity centres

set.seed(4)
# Generate simulated population
AClocs <- rSSI(r=territory, n=n, win=win)
popn <- data.frame(x=AClocs$x, y=AClocs$y)
class(popn) <- c("popn", "data.frame")

## Trapping design
## ===============
spacing <- 2500  # trap spacing
nOcc <- 90  # 3 months

## Plot the two trap layouts
## =========================

windows(9, 5)
par(mfrow=1:2)

# One big block layout
# --------------------
# The traps:
traps1 <- make.grid(16, 6, spacing=spacing, detector='proximity',
    originxy=c(14380, 21865))
plot(win, main="")
points(AClocs, pch=19)
plot(traps1, add=TRUE)
# Simulate captures:
simCH <- sim.capthist(traps1, 
  popn = popn,
  detectfn = 0, 
  detectpar = list(g0 = g0, sigma = sigma),
  noccasions = nOcc)
plot(simCH, add=TRUE, tracks=TRUE, rad=500)
points(AClocs, pch=19)

# Multi-cluster layout
# --------------------
# Get traps layout
cluster <- make.grid(nx=3, ny=3, spacing=spacing, detector='proximity')
pattern <- c(5, 4) # Actual number of clusters will be less than this.
for(i in 1:100)  {
  traps <- make.systematic(pattern, cluster, spdf)
  if(nrow(traps) %in% 95:97)
    break
}
plot(win, main="")
points(AClocs, pch=19)
plot(traps, add=TRUE)
# simulated capture histories for current configuration:
simCH <- sim.capthist(traps, 
      popn = popn,
      detectfn = 0, 
      detectpar = list(g0 = g0, sigma = sigma),
      noccasions = nOcc)
plot(simCH, add=TRUE, tracks=TRUE, rad=500)
points(AClocs, pch=19)


# Load the results of 300 simulations
# ===================================

Dhat <- read.csv("simulation_results.csv")
head(Dhat)  # units are animals per 100 sq km

# Look at failures
# ================
# How many failures?
round(colMeans(is.na(Dhat)) * 100, 2) # Percentage of failures

# Do some nice plots
# ==================
# beeswarm and box plots
windows(5, 5)
# Make a temporary data frame with NAs replaced by a negative number
tmp <- Dhat
tmp[is.na(Dhat)] <- -0.025
# Do the beeswarm plot
beeswarm(tmp, method='hex', xpd=TRUE, yaxs='i', las=1,
    pwcol=c('blue', 'red')[(tmp<0)+1], cex=0.7, xaxt='n',
    xlab="Trap configuration", ylab="Density estimate (animals/100 km2)",
    main="Effect of trap configuration\nTrue density = 0.25/100 sq km") 
axis(1, 1:2, c("single 6x16 block", "multiple 3x3 clusters"), tick=FALSE)
# Add the boxplot
boxplot(Dhat, add=TRUE, col="#0000ff22", yaxt='n', xaxt='n', pch=NA)
# Add horizontal lines and legend
abline(h=0)
abline(h=0.25, col='red', lwd=3, lty=3)
legend('topleft', c("Valid estimates", "Data not usable"), pch=1, cex=0.8,
  col=c('blue', 'red'), bty='n')
