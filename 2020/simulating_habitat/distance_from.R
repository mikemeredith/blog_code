
# R code for the blog post at
#  https://mmeredith.net/blog/2020/simulating_habitat.htm

library(fitdistrplus)
?fitdistrplus
library(wiqid)  # For AICtable

# Distance from ...
# =================

# Distance from road
# ------------------

source("distance_covar.R")

x <- x / 1000  # convert to km

descdist(x)
# Near normal, but high skew

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
# Weibull has lowest AIC, gamma second

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# Weibull not too bad

# Do the plot:
palette("R3")
denscomp(list(fitN, fitW), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="Distance from road",
    xlab="Distance from road, km")

# Log-transformed data
# --------------------

x <- log(x)

descdist(x)
# Close to normal, but higher skew and kurtosis

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
    demp=TRUE, xlegend = "topleft",
    main="Distance from road",
    xlab="Log(distance from road)")

# square-root-transformed data
# ----------------------------

source("distance_covar.R")
x <- x / 1000  # convert to km
x <- sqrt(x)

descdist(x)
# Between uniform and normal, no skew

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
# Weibull has lowest AIC

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# Weibull not too bad

# Do the plot:
palette("R3")
denscomp(list(fitN, fitW), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="Distance from road",
    xlab="Sqrt(distance from road)")

# Yellowstone roads
# =================

library(raster)
rddist <- raster("Yellowstone_roadDist.tif")
x <- values(rddist)
x <- x[!is.na(x)]   # pixels outside the park are NA

descdist(x)
# Far from anywhere

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
# Weibull has lowest AIC

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# Weibull Q-Q plot not good

# Do the plot:
palette("R3")
denscomp(list(fitN, fitW), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="Yellowstone NP roads",
    xlab="Distance from road, km")

# square-root-transformed data
# ----------------------------

# x <- values(rddist) / 1000 # km
# x <- x[!is.na(x)]   # pixels outside the park are NA
x <- sqrt(x)

descdist(x)
# Close to normal

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
# Weibull has lowest AIC

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# Weibull not too bad

# Do the plot:
palette("R3")
denscomp(list(fitN, fitW), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="Yellowstone NP roads",
    xlab="Sqrt(distance from road)")


# Distance from boundary
# ======================

library(makeJAGSmask) # version 0.1.1.9005 or later
data(simSCR)
# help(simSCR)

x <- values(simSCR$patchR) / 1000 # km
x <- x[!is.na(x)]

descdist(x)
# Close to uniform

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
# Weibull has lowest AIC

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# None very good

# Do the plot:
palette("R3")
denscomp(list(fitN, fitW), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="Mask boundary",
    xlab="Distance from boundary")

# Square-root transform
# ---------------------
x <- sqrt(x)

descdist(x)
# Very close to uniform

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
# Normal has lowest AIC, Weibull second

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# None very good

# Do the plot:
palette("R3")
denscomp(list(fitN, fitW), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="Mask boundary",
    xlab="sqrt(distance from boundary)")
