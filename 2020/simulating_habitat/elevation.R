
# R code for the blog post at
#  https://mmeredith.net/blog/2020/simulating_habitat.htm

library(fitdistrplus)
?fitdistrplus

# Elevation
# =========

# MHB data
# --------

library(AHMbook)
data(MHB2014)
str(MHB2014)

x <- MHB2014$sites$ele / 1000

descdist(x)
# Close to uniform in the Cullen-Frey plot

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
# Gamma has highest AIC.

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# Weibull best in GOF plots, but none good

# Do the plot:
palette("R3")
denscomp(list(fitN, fitG), fitlwd=2, datacol="gray90", demp=TRUE, xlegend = "topright",
    main="AHMbook::MHB2014 elevation",
    xlab="Elevation, m/1000")

# Switzerland data
# ----------------

library(unmarked)
data(Switzerland)
str(Switzerland)

x <- Switzerland$elevation / 1000

descdist(x)
# Close to uniform but higher skew

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
# log-normal has highest AIC by far.

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
denscomp(list(fitN, fitL), fitlwd=2, datacol="gray90", demp=TRUE, xlegend = "topright",
    main="unmarked::Switzerland elevation",
    xlab="Elevation, m/1000")

# Hubbard Brook data
# ------------------

library(AHMbook)
data(HubbardBrook)
str(HubbardBrook)

x <- HubbardBrook$sitecov$Elev / 1000

descdist(x)
# Between normal and uniform

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
# Weibull has highest AIC.

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# Weibull looks good

# Do the plot:
palette("R3")
denscomp(list(fitN, fitW), fitlwd=2, datacol="gray90", demp=TRUE, xlegend = "topright",
    main="AHMbook::HubbardBrook elevation",
    xlab="Elevation, m/1000")

