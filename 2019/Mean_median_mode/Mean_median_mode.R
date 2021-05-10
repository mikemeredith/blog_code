
# R code for the web page at
#  https://mmeredith.net/blog/2019/Mean_median_mode.htm


# How many coins in the jar?

# Posterior probability

xx <- 500:900
dd <- pmin(xx-500, (900-xx)/3)
d <- dd / sum(dd)
plot(xx, d, type='l', lwd=2,
  xlab="Number of coins", ylab="Probability")
temp <- cumsum(d)
range(temp)
( median <- which.min(abs(temp - 0.5)) + 500 )
abline(v=600, lwd=2)
abline(v=median, lwd=2, col='blue')
( mean <- sum(xx*d) )
abline(v=mean, col='red', lwd=2)

text(562.1, 0.0045, 'mode', pos=2)
arrows(562.1, 0.0045, 600, 0.004972, length=0.1)
text(703, 0.0046, 'median', pos=4)
arrows(703, 0.0046, 656.039604, 0.0046, length=0.1)
text(703, 0.0011, 'mean', pos=4)
arrows(703, 0.0011, 666.840684, 0.0011, length=0.1)

# Which to use?


myguess <- 550:700

# 1. $100 prize for right answer, nothing if wrong
# Expected value
X1 <- d * 100
xx[which.max(X1)]

# 2. $100 prize for right answer, - $1 for each unit of error
# Expected value
X2 <- numeric(length(d))
for(i in seq_along(xx)) {
  loss <- abs(xx - xx[i]) * 1
  X2[i] <- sum(loss * d)
}
plot(xx, X2)
xx[which.min(X2)]

# 3. $100 prize for right answer, - 5c for each squared unit of error
# Expected value
X3 <- numeric(length(d))
for(i in seq_along(xx)) {
  loss <- (xx - xx[i])^2 * 0.05
  X3[i] <- sum(loss * d)
}
plot(xx, X3)
xx[which.min(X3)]

# The 3 Reward curves

error <- 0:100
plot(error, type='n', ylab="Reward ($)", xlab="Error", las=1)
segments(0, 0, 0, 100, lwd=3)
segments(0, 0, 100, 0, lwd=3)
lines(error, 100-error, lwd=3, col='blue')
lines(1:45, 100-(1:45)^2*0.05, lwd=3, col='red')

text(6.48, 27, 'Ana', pos=4)
arrows(6.48, 27, 0.57, 25.75, length=0.1)
text(62.35, 56.33, 'Beth', pos=4)
arrows(62.35, 56.33, 57.53, 43.59, length=0.1)
text(48.33, 18.11, 'Carl', pos=4)
arrows(48.33, 18.11, 41.97, 13.87, length=0.1)
