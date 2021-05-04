
# Script to produce plots for the blog page
# https://mmeredith.net/blog/2013/1308_Circular_home_ranges.htm

# Do we assume that home ranges are circular when we do SECR?

library(MASS)
library(spatstat)

# Define functions used later to fit curves to histograms
# =======================================================
# half-normal curve:
ssHnorm <- function(sig, x) {
  tmp <- dnorm(x, 0, sig)
  hnorm <- tmp / sum(tmp)
  sum((hnorm - dens)^2)
}

# exponential curve:
ssExp <- function(sig, x) {
  tmp <- dexp(x, sig)
  fitted <- tmp / sum(tmp)
  sum((fitted - dens)^2)
}

# Andy Royle's hybrid function:
# NB: scale > 0 and 1 <= shape <= 2 and x >= 0
# shape = 1 is exponential
# shape = 2 is Gaussian kernel
dhybrid <- function(x, shape, scale)
  exp(-(x/scale)^shape)

ssHybrid <- function(pars, x) {
  scale <- exp(pars[1]) # must be >1
  shape <- plogis(pars[2]) + 1  # must be 1<shape<2
  tmp <- dhybrid(x, shape, scale)
  fitted <- tmp / sum(tmp)
  sum((fitted - dens)^2)
}

# Murray calls this the "w exponential":
dwex <- function(x, w, scale)
  ifelse(x < w, 1, exp(-(x-w) / scale))

ssWex <- function(pars, x) {
  scale <- exp(pars[1]) # must be >1
  w <- exp(pars[2])
  tmp <- dwex(x, w, scale)
  fitted <- tmp / sum(tmp)
  sum((fitted - dens)^2)
}

# x and y vlaues to plot circles
angle <- seq(0, 2*pi, len=300)
xx <- sin(angle)
yy <- cos(angle)

# good size of plotting window
windows(9, 5)


# plot for circular, fuzzy home range
# (symmetrical bivariate normal)
# ===================================
Sigma <- matrix(c(1,0,0,1),2,2)   # symmetric bivariate normal
set.seed(123)
simdat <- mvrnorm(2e4, c(0,0), Sigma)
simdatP <- scale(simdat, scale=FALSE) # 2e4 enough for scatter plots
simdat <- mvrnorm(1e7, c(0,0), Sigma)
simdatH <- scale(simdat, scale=FALSE) # need more for histograms

dist <- sqrt(rowSums(simdatH^2))
med <- median(dist)
cutpoints <- seq(0, 4, length=51)
binWidth <- mean(diff(cutpoints))
midpoints <- (cutpoints[-1] + cutpoints[1:50]) / 2
counts <- hist(dist, c(cutpoints, max(dist)), plot=FALSE)$counts[1:50]
binArea <- pi*(cutpoints[-1]^2 - cutpoints[1:50]^2)
dens0 <- counts / binArea
dens <- dens0 / sum(dens0)

# fit a half-normal curve to the densities:
sigFit <- optimise(ssHnorm, c(0,1000), x=midpoints)$minimum
tmp <- dnorm(midpoints, 0, sigFit)
const <- sum(tmp)
hnorm <- tmp / const
g0 <- dnorm(0, 0, sigFit) / const

# plot it:
par(mfrow=1:2)
eqscplot(simdatP, pch='.', ylim=c(-3,3), xlim=c(-3,3), xpd=TRUE, xaxt='n', yaxt='n',
  col='darkgrey', bty='n')
arrows(-1,1,0,0,,length=0.1, lwd=2)
text(-1,1, "AC", pos=3, cex=1.5)
points(xx*med, yy*med, col='blue', type='l', lwd=2) 
arrows(0, 0, med, 0, length=0.1, code=3, lwd=2, col='blue')
text(0.5, -0.2, "d", cex=1.5, col='blue')
points(0, 0, pch=19, col='red')
arrows(1,2,0,med,,length=0.1, lwd=2)
text(1,2, "trap", pos=3, cex=1.5)
points(0, med, pch=15, col='red')

plot( cutpoints[1:50], dens, type='n', ylim=c(0, g0), xaxs='i', yaxs='i',
  xaxt='n', yaxt='n', bty='l', 
  xlab="Distance between Activity Centre and trap",
  ylab="Probability of capture", )
rect(cutpoints[1:50], 0, cutpoints[-1], dens, col='grey92')
points(midpoints, hnorm, type='l', col='blue', lwd=2, xpd=TRUE)
axis(1, 0, labels = 0,)
axis(2, 0, 0, las=1, xpd=TRUE)
axis(2, g0, 'g0', las=1, xpd=TRUE, col.axis='blue', col.ticks='blue')
abline(v=med, col='blue', lwd=2, lty=3)
axis(1, med, "d", col.axis='blue', col.ticks='blue')


# plot for elongated home range
# (asymmetric bivariate normal)
# =============================
Sigma <- matrix(c(5, 0.7, 0.7, 1/5),2,2) # asymmetric bivariate normal
set.seed(123)
simdat <- mvrnorm(2e4, c(0,0), Sigma)
simdatP <- scale(simdat, scale=FALSE) # for plots
simdat <- mvrnorm(1e7, c(0,0), Sigma)
simdatH <- scale(simdat, scale=FALSE) # for histograms

dist <- sqrt(rowSums(simdatH^2))
med <- median(dist)
cutpoints <- seq(0, 8, length=51)
binWidth <- mean(diff(cutpoints))
midpoints <- (cutpoints[-1] + cutpoints[1:50]) / 2
counts <- hist(dist, c(cutpoints, max(dist)), plot=FALSE)$counts[1:50]
binArea <- pi*(cutpoints[-1]^2 - cutpoints[1:50]^2)
dens0 <- counts / binArea
dens <- dens0 / sum(dens0)

# fit a half-normal curve to the densities:
sigFit <- optimise(ssHnorm, 
c(0,1000), x=midpoints)$minimum
tmp <- dnorm(midpoints, 0, sigFit)
const <- sum(tmp)
hnorm <- tmp / const
g0hn <- dnorm(0, 0, sigFit) / const

# fit an exponential curve to the densities:
rateFit <- optimise(ssExp, c(0,1), x=midpoints)$minimum
tmp <- dexp(midpoints, rateFit)
const <- sum(tmp)
exp <- tmp / const
g0e <- dexp(0, rateFit) / const

# fit a hybrid curve:
result <- nlm(ssHybrid, c(0,0), x=midpoints)
scale <- exp(result$estimate[1])
shape <- plogis(result$estimate[2]) + 1  # must be 1<shape<2
# tmp <- exp(- midpoints^shape / scale)
tmp <- dhybrid(midpoints, shape, scale)
const <- sum(tmp)
fitted <- tmp / const
g0h <- dhybrid(0, shape, scale) / const

ymax <- max(g0hn, g0e, g0h, dens)

# plot it:
par(mfrow=1:2)
eqscplot(simdatP, pch='.', ylim=c(-2,2), xlim=c(-5,5), xpd=TRUE, xaxt='n', yaxt='n',
  col='gray', bty='n')
points(0, 0, pch=19, col='red')

points(xx*med, yy*med, col='blue', type='l', lwd=2) 
arrows(0, 0, med, 0, length=0.1, code=3, lwd=2, col='blue')
text(0.5, -0.2, "d", cex=1.5, col='blue')
points(0, 0, pch=19, col='red')
arrows(0,3,0,med,,length=0.1, lwd=2)
text(0,3, "low probability of capture", pos=3, cex=1.5)
points(0, med, pch=15, col='red')
arrows(-med,-2,-med,0,,length=0.1, lwd=2)
text(0,-2, "high probability of capture", pos=1, cex=1.5)
points(-med, 0, pch=15, col='red')

plot( cutpoints[1:50], dens, type='n', ylim=c(0, ymax), xaxs='i', yaxs='i',
  xaxt='n', yaxt='n', bty='l', 
  xlab="Distance between Activity Centre and trap", ylab="")
rect(cutpoints[1:50], 0, cutpoints[-1], dens, col='grey92')
points(c(0, midpoints), c(g0e, exp), type='l', col='red', lwd=2, xpd=TRUE)
points(c(0, midpoints), c(g0hn, hnorm), type='l', col='blue', lwd=2, xpd=TRUE)
points(c(0, midpoints), c(g0h, fitted), type='l', lty=2, lwd=2, xpd=TRUE)
axis(1, 0, labels = 0,)
axis(2, 0, 0, las=1, xpd=TRUE)
axis(2, g0hn, 'hnorm-g0', las=1, xpd=TRUE, col.axis='blue', col.ticks='blue')
axis(2, g0e, 'expo-g0', las=1, xpd=TRUE, col.axis='red', col.ticks='red')
axis(2, g0h, 'hybrid-g0', las=1, xpd=TRUE)
legend('topright', c("half-normal", "exponential", "hybrid"), col=c('blue', 'red', 'black'),
  lty=c(1,1,2), lwd=2, title="Detection function", bty='n')


# plot for irregular home range
# (mixture of previous two)
# =============================
Sigma1 <- matrix(c(1,0,0,1),2,2)   # symmetric bivariate normal
Sigma2 <- matrix(c(5, 0.7, 0.7, 1/5),2,2) # assymmetric
set.seed(123)
simdat <- rbind(mvrnorm(1e4, c(2,0), Sigma1),
                mvrnorm(1e4, c(0,0), Sigma2))
simdatP <- scale(simdat, scale=FALSE)
simdat <- rbind(mvrnorm(5e6, c(2,0), Sigma1),
                mvrnorm(5e6, c(0,0), Sigma2))
simdatH <- scale(simdat, scale=FALSE)

dist <- sqrt(rowSums(simdatH^2))
med <- median(dist)
cutpoints <- seq(0, 8, length=51)
binWidth <- mean(diff(cutpoints))
midpoints <- (cutpoints[-1] + cutpoints[1:50]) / 2
counts <- hist(dist, c(cutpoints, max(dist)), plot=FALSE)$counts[1:50]
binArea <- pi*(cutpoints[-1]^2 - cutpoints[1:50]^2)
dens0 <- counts / binArea
dens <- dens0 / sum(dens0)

# fit an exponential curve to the densities:
rateFit <- optimise(ssExp, c(0,1), x=midpoints)$minimum
tmp <- dexp(midpoints, rateFit)
const <- sum(tmp)
exp <- tmp / const
g0e <- dexp(0, rateFit) / const
ymax <- max(g0e, dens)

# fit a hybrid curve:
result <- nlm(ssHybrid, c(0,0), x=midpoints)
scale <- exp(result$estimate[1])
shape <- plogis(result$estimate[2]) + 1  # must be 1<shape<2
tmp <- dhybrid(midpoints, shape, scale)
const <- sum(tmp)
fitted <- tmp / const
g0h <- dhybrid(0, shape, scale) / const

ymax <- max(g0e, g0h, dens)

# plot it:
par(mfrow=1:2)
eqscplot(simdatP, pch='.', ylim=c(-2,2), xlim=c(-5,5), xpd=TRUE, xaxt='n', yaxt='n',
  col='grey', bty='n')
points(0, 0, pch=19, col='red')

points(xx*med, yy*med, col='blue', type='l', lwd=2) 
arrows(0, 0, med, 0, length=0.1, code=3, lwd=2, col='blue')
text(0.5, -0.2, "d", cex=1.5, col='blue')
points(0, 0, pch=19, col='red')

plot( cutpoints[1:50], dens, type='n', ylim=c(0, ymax), xaxs='i', yaxs='i',
  xaxt='n', yaxt='n', bty='l', 
  xlab="Distance between Activity Centre and trap",
  ylab="Probability of capture", )
rect(cutpoints[1:50], 0, cutpoints[-1], dens, col='grey92')
points(c(0, midpoints), c(g0e, exp), type='l', col='red', lwd=2, xpd=TRUE)
points(c(0, midpoints), c(g0h, fitted), type='l', lty=2, lwd=2, xpd=TRUE)
axis(1, 0, labels = 0,)
axis(2, 0, 0, las=1, xpd=TRUE)
axis(2, g0e, 'expo-g0', las=1, xpd=TRUE, col.axis='red', col.ticks='red')
axis(2, g0h, 'hybrid-g0', las=1, xpd=TRUE)
legend('topright', c("exponential", "hybrid"), col=c('red', 'black'),
  lty=c(1,2), lwd=2, title="Detection function", bty='n')


# plot for hard-edged irregular home range
# ========================================
HR <- data.frame(
  x = c(-6.254, -5.475, -1.860, -1.144,  0.072,  1.817,  3.063,
         3.593,  3.562,  2.596, -1.300, -5.631),
  y = c(-0.315,  0.371,  0.682,  1.181,  2.490,  2.739,  1.555,
         0.371, -0.751, -1.561, -1.437, -1.250))
win <- owin(poly=HR[nrow(HR):1, ])
set.seed(2)
simdat <- as.data.frame(runifpoint(2e4, win))
simdatP <- scale(simdat, scale=FALSE)
simdat <- as.data.frame(runifpoint(1e6, win))
simdatH <- scale(simdat, scale=FALSE)

dist <- sqrt(rowSums(simdatH^2))
med <- median(dist)
cutpoints <- seq(0, 8, length=51)
binWidth <- mean(diff(cutpoints))
midpoints <- (cutpoints[-1] + cutpoints[1:50]) / 2
counts <- hist(dist, c(cutpoints, max(dist)), plot=FALSE)$counts[1:50]
binArea <- pi*(cutpoints[-1]^2 - cutpoints[1:50]^2)
dens0 <- counts / binArea
dens <- dens0 / sum(dens0)

# fit a half-normal curve to the densities:
sigFit <- optimise(ssHnorm, c(0,1000), x=midpoints)$minimum
tmp <- dnorm(midpoints, 0, sigFit)
const <- sum(tmp)
hnorm <- tmp / const
g0hn <- dnorm(0, 0, sigFit) / const

# fit a WEX curve:
result <- nlm(ssWex, c(0,0), x=midpoints)
scale <- exp(result$estimate[1])
w <- exp(result$estimate[2])
tmp <- dwex(midpoints, w, scale)
const <- sum(tmp)
fitted <- tmp / const
g0w <- dwex(0, w, scale) / const

ymax <- max(g0hn, g0w, dens)

# plot it:
par(mfrow=1:2)
eqscplot(simdatP, pch='.', ylim=c(-2,2), xlim=c(-5,5), xpd=TRUE, xaxt='n', yaxt='n',
  col='grey', bty='n')
points(0, 0, pch=19, col='red')

points(xx*med, yy*med, col='blue', type='l', lwd=2) 
arrows(0, 0, med, 0, length=0.1, code=3, lwd=2, col='blue')
text(0.5, -0.2, "d", cex=1.5, col='blue')
points(0, 0, pch=19, col='red')

plot( cutpoints[1:50], dens, type='n', ylim=c(0, ymax), xaxs='i', yaxs='i',
  xaxt='n', yaxt='n', bty='l', 
  xlab="Distance between Activity Centre and trap",
  ylab="Probability of capture", )
rect(cutpoints[1:50], 0, cutpoints[-1], dens, col='grey92')
points(c(0, midpoints), c(g0hn, hnorm), type='l', col='blue', lwd=2, xpd=TRUE)
points(c(0, midpoints), c(g0w, fitted), type='l', lty=2, lwd=2, xpd=TRUE)
axis(1, 0, labels = 0,)
axis(2, 0, 0, las=1, xpd=TRUE)
axis(2, g0hn, 'hnorm-g0', las=1, xpd=TRUE, col.axis='blue', col.ticks='blue')
axis(2, g0w, 'wex-g0', las=1, xpd=TRUE)
legend('topright', c("half-normal", "w exponential"), col=c('blue', 'black'),
  lty=c(1,2), lwd=2, title="Detection function", bty='n')


# plot for hard-edged circular home range
# =======================================
CHR <- data.frame(x = xx * 2, y = yy * 2)
win <- owin(poly = CHR[nrow(CHR):1, ])
set.seed(123)
simdat <- as.data.frame(runifpoint(2e4, win))
simdatP <- scale(simdat, scale=FALSE)
simdat <- as.data.frame(runifpoint(1e6, win))
simdatH <- scale(simdat, scale=FALSE)

dist <- sqrt(rowSums(simdatH^2))
med <- median(dist)
cutpoints <- seq(0, 8, length=51)
binWidth <- mean(diff(cutpoints))
midpoints <- (cutpoints[-1] + cutpoints[1:50]) / 2
counts <- hist(dist, c(cutpoints, max(dist)), plot=FALSE)$counts[1:50]
binArea <- pi*(cutpoints[-1]^2 - cutpoints[1:50]^2)
dens0 <- counts / binArea
dens <- dens0 / sum(dens0)

# fit a half-normal curve to the densities:
sigFit <- optimise(ssHnorm, c(0,1000), x=midpoints)$minimum
tmp <- dnorm(midpoints, 0, sigFit)
const <- sum(tmp)
hnorm <- tmp / const
g0hn <- dnorm(0, 0, sigFit) / const

# fit a WEX curve:
result <- nlm(ssWex, c(0,0), x=midpoints)
scale <- exp(result$estimate[1])
w <- exp(result$estimate[2])
tmp <- dwex(midpoints, w, scale)
const <- sum(tmp)
fitted <- tmp / const
g0w <- dwex(0, w, scale) / const

ymax <- max(g0hn, g0w, dens)

# plot it:
par(mfrow=1:2)
eqscplot(simdatP, pch='.', ylim=c(-3,3), xlim=c(-3,3), xpd=TRUE, xaxt='n', yaxt='n',
  col='grey', bty='n')
points(0, 0, pch=19, col='red')

points(xx*med, yy*med, col='blue', type='l', lwd=2) 
arrows(0, 0, med, 0, length=0.1, code=3, lwd=2, col='blue')
text(0.5, -0.2, "d", cex=1.5, col='blue')
points(0, 0, pch=19, col='red')

plot( cutpoints[1:50], dens, type='n', ylim=c(0, ymax), xaxs='i', yaxs='i',
  xaxt='n', yaxt='n', bty='l', 
  xlab="Distance between Activity Centre and trap",
  ylab="Probability of capture", )
rect(cutpoints[1:50], 0, cutpoints[-1], dens, col='grey92')
points(c(0, midpoints), c(g0hn, hnorm), type='l', col='blue', lwd=2, xpd=TRUE)
points(c(0, midpoints), c(g0w, fitted), type='l', lty=2, lwd=2, xpd=TRUE)
axis(1, 0, labels = 0,)
axis(2, 0, 0, las=1, xpd=TRUE)
axis(2, g0hn, 'hnorm-g0', las=1, xpd=TRUE, col.axis='blue', col.ticks='blue')
axis(2, g0w, 'wex-g0', las=1, xpd=TRUE)
legend('topright', c("half-normal", "w exponential"), col=c('blue', 'black'),
  lty=c(1,2), lwd=2, title="Detection function", bty='n')
