
# R code for the blog post at
#  https://mmeredith.net/blog/2021/binary_vs_count.htm


# The simulated data
# ==================

elev <- seq(20, 1980, 40)
( n <- length(elev) )
elevS <- wiqid::standardize(elev)
elevS2 <- elevS^2

beta0 <- log(2.5)
beta1 <- -1.2
beta2 <- -1

lambda <- exp(beta0 + beta1*elevS + beta2*elevS2)
psi <- 1 - exp(-lambda) # True prob of occupancy

set.seed(42)
C <- rpois(n, lambda)
B <- as.numeric(C > 0)

# Binary data
# ===========

plot(elev, B, pch=16, col='red', xlab="Elevation, m", ylab="Presence/absence")

b0 <- glm(B ~ 1, family="binomial")
b1 <- glm(B ~ elevS, family="binomial")
b2 <- glm(B ~ elevS + elevS2, family="binomial")
AIC(b0, b1, b2)
#    df      AIC
# b0  1 57.10799
# b1  2 26.11651
# b2  3 26.44560

summary(b2)

# Plot the fitted values
plot(elev, B, pch=16, col='red', xlab="Elevation (m)", ylab="Presence/absence")
points(elev, fitted(b1), cex=1.5)
points(elev, fitted(b2), cex=1.5, col='blue')
lines(elev, 1 - exp(-lambda))


# Count data
# ==========

plot(elev, C, type='h', col='red', xlab="Elevation (m)", ylab="Count")
points(elev, C, pch=16, col='red')

c0 <- glm(C ~ 1, family="poisson")
c1 <- glm(C ~ elevS, family="poisson")
c2 <- glm(C ~ elevS + elevS2, family="poisson")
AIC(c0, c1, c2)
#    df      AIC
# c0  1 214.8638
# c1  2 185.7632
# c2  3 153.9122

summary(c2)
# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.1165     0.1322   8.445  < 2e-16 ***
# elevS        -1.1250     0.2065  -5.447 5.12e-08 ***
# elevS2       -0.8953     0.1825  -4.905 9.36e-07 ***

plot(elev, C, type='h', col='red', xlab="Elevation (m)", ylab="Count")
points(elev, C, pch=16, col='red')
points(elev, fitted(c1))
points(elev, fitted(c2), col='blue')
lines(elev, lambda)
