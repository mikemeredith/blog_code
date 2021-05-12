
# R code for the blog post at
#  https://mmeredith.net/blog/2020/simulating_habitat.htm

library(fitdistrplus)
?fitdistrplus
library(wiqid)  # For AICtable

# Patch size
# ==========

# CSWA data
# ---------

library(AHMbook)
data(cswa)
str(cswa)

x <- cswa$spotMaps$parea

descdist(x)
# Out in the wilderness, high skew and kurtosis

# Fit a range of distributions:
( fitU <- fitdist(x, "unif") )
( fitN <- fitdist(x, "norm") )
( fitL <- fitdist(x, "lnorm") )
( fitG <- fitdist(x, "gamma") )
( fitE <- fitdist(x, "exp") )
( fitW <- fitdist(x, "weibull") )
# All MLE fits work
AICtable(data.frame(n=c(2,2,2,2,1,2),
    aic=c(norm=fitN$aic,
    lnorm=fitL$aic,
    uni=fitU$aic,
    gam=fitG$aic,
    exp=fitE$aic,
    wei=fitW$aic)))
# log-normal has lowest AIC

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# None can cope with the largest value

# Do the plot:
palette("R3")
denscomp(list(fitN, fitL, fitE), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="AHMbook::cswa parea",
    xlab="Patch area, ha")

# Adirondack wetlands
# -------------------

x <- read.csv("Adirondack_Wetlands.csv", comment="#")$ha
range(x)
mean(x > 200)

x <- x[x < 200]

descdist(x)
# Out in the wilderness, high skew and kurtosis

# Fit a range of distributions:
( fitU <- fitdist(x, "unif") )
( fitN <- fitdist(x, "norm") )
( fitL <- fitdist(x, "lnorm") )
( fitG <- fitdist(x, "gamma") )
( fitE <- fitdist(x, "exp") )
( fitW <- fitdist(x, "weibull") )
# All MLE fits work
AICtable(data.frame(n=c(2,2,2,2,1,2),
    aic=c(norm=fitN$aic,
    lnorm=fitL$aic,
    uni=fitU$aic,
    gam=fitG$aic,
    exp=fitE$aic,
    wei=fitW$aic)))
# log-normal has lowest AIC by far.

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# Nothing looks good, not even log-normal

# Do the plot:
palette("R3")
denscomp(list(fitN, fitL, fitE), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright", xlim=c(0, 200),
    main="Adirondack wetlands",
    xlab="Wetland area, ha")
text(150, 0.1, "0.3% greater than 200ha\nLargest 1282ha")

# Log-transformed data
# --------------------

x <- read.csv("Adirondack_Wetlands.csv", comment="#")$ha
x <- log(x)

descdist(x)
# Close to normal, somewhat higher kurtosis

# Fit a range of distributions:
( fitU <- fitdist(x, "unif") )
( fitN <- fitdist(x, "norm") )
( fitL <- fitdist(x, "lnorm") )    # failed
( fitG <- fitdist(x, "gamma") )    # failed
( fitE <- fitdist(x, "exp") )      # failed
( fitW <- fitdist(x, "weibull") )  # failed

# try MME fits, but no SEs or AIC
# ( fitU <- fitdist(x, "unif", method="mme") )
# ( fitN <- fitdist(x, "norm", method="mme") )
( fitL <- fitdist(x, "lnorm", method="mme") )    # failed
( fitG <- fitdist(x, "gamma", method="mme") )    # failed
( fitE <- fitdist(x, "exp", method="mme") )      # failed
( fitW <- fitdist(x, "weibull", method="mme") )  # failed

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
# Normal looks good, except for tails in the Q-Q plot

# Do the plot:
palette("R3")
denscomp(list(fitN), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="Adirondack wetlands",
    xlab="Log(wetland area)")
