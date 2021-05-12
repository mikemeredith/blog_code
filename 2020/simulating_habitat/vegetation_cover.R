
# R code for the blog post at
#  https://mmeredith.net/blog/2020/simulating_habitat.htm

library(fitdistrplus)
?fitdistrplus
library(wiqid)  # For AICtable

# Vegetation cover
# ================

# MHB data
# --------

library(AHMbook)
data(MHB2014)
str(MHB2014)

x <- MHB2014$sites$forest / 100
hist(x, main="AHMbook::MHB2014 forest", xlab="Proportion of forest cover")
# Failed to fit most distributions to these data

# Check for proportion of zeros
mean(x == 0) # 12%

# Remove zeros and try fitting things again
x <- x[x > 0]

descdist(x)
# Close to uniform in the Cullen-Frey plot

# Fit a range of distributions:
( fitU <- fitdist(x, "unif") )
( fitN <- fitdist(x, "norm") )
( fitL <- fitdist(x, "lnorm") )
( fitG <- fitdist(x, "gamma") )
( fitE <- fitdist(x, "exp") )
( fitW <- fitdist(x, "weibull") )
( fitB <- fitdist(x, "beta") )
# All MLE fits work
AICtable(data.frame(n=c(2,2,2,2,2,1,2),
    aic=c(norm=fitN$aic,
    beta=fitB$aic,
    lnorm=fitL$aic,
    uni=fitU$aic,
    gam=fitG$aic,
    exp=fitE$aic,
    wei=fitW$aic)))
# beta has lowest AIC, uniform second.

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
plot(fitB) ; title(main="Beta", outer=TRUE, line=-1)
# Beta best in GOF plots, Weibull second, better than uniform

# Do the plot:
palette("R3")
denscomp(list(fitN, fitB, fitW), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="AHMbook::MHB2014 forest less zeros",
    xlab="Proportion of forest cover")

# Switzerland data
# ----------------

library(unmarked)
data(Switzerland)
str(Switzerland)

x <- Switzerland$forest / 100
hist(x, main="unmarked::Switzerland forest", xlab="Proportion of forest cover")
# Failed to fit most distributions to these data

# Check for proportion of zeros
mean(x == 0) # 26%

# Remove zeros and try fitting things again
x <- x[x > 0]

descdist(x)
# Close to uniform but higher skew


# Fit a range of distributions:
( fitU <- fitdist(x, "unif") )
( fitN <- fitdist(x, "norm") )
( fitL <- fitdist(x, "lnorm") )
( fitG <- fitdist(x, "gamma") )
( fitE <- fitdist(x, "exp") )
( fitW <- fitdist(x, "weibull") )
( fitB <- fitdist(x, "beta") )  # failed
# All MLE fits work except Beta
AICtable(data.frame(n=c(2,2,2,2,1,2),
    aic=c(norm=fitN$aic,
    lnorm=fitL$aic,
    uni=fitU$aic,
    gam=fitG$aic,
    exp=fitE$aic,
    wei=fitW$aic)))
# Weibull has highest AIC by far.

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# Nothing looks good, not even Weibull

# Do the plot:
palette("R3")
denscomp(list(fitN, fitW), fitlwd=2, datacol="gray90", demp=TRUE, xlegend = "topright",
    main="unmarked::Switzerland forest less zeros",
    xlab="Proportion of forest cover")

# wagtail data
# ------------

library(AHMbook)
data(wagtail)
str(wagtail)

x <- wagtail$grass / 100
hist(x, main="AHMbook::wagtail grass", xlab="Proportion of permanent grass")

mean(x == 0)  # 63%
x <- x[x > 0]
mean(x == 1)  # 11%
