
# R code for the blog post at
#  https://mmeredith.net/blog/2020/simulating_habitat.htm

library(fitdistrplus)
?fitdistrplus
library(wiqid)  # For AICtable

# Other covariates
# ================

# Housing density
# ---------------

library(AHMbook)
data(MesoCarnivores)
str(MesoCarnivores)


x <- MesoCarnivores$sitecovs$HDens_5km

descdist(x)
# Off the planet

# Fit a range of distributions:
( fitU <- fitdist(x, "unif") )
( fitN <- fitdist(x, "norm") )
( fitL <- fitdist(x, "lnorm") )  # failed
( fitG <- fitdist(x, "gamma") )  # failed
( fitE <- fitdist(x, "exp") )
( fitW <- fitdist(x, "weibull") )# failed
# All MLE fits work
AICtable(data.frame(n=c(2,2,1),
    aic=c(norm=fitN$aic,
    uni=fitU$aic,
    exp=fitE$aic)))
# Exponential has lowest AIC
# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
# All disasterous

# Do the plot:
palette("R3")
denscomp(list(fitN, fitE), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="AHMbook::MesoCarnivores HDens",
    xlab="Housing density")

# Log-transformed data
# --------------------

x <- log(x + 1)

descdist(x)
# Still out in the wilderness

# Fit a range of distributions:
( fitU <- fitdist(x, "unif") )
( fitN <- fitdist(x, "norm") )
( fitL <- fitdist(x, "lnorm") )    # failed
( fitG <- fitdist(x, "gamma") )    # failed
( fitE <- fitdist(x, "exp") )
( fitW <- fitdist(x, "weibull") )  # failed
AICtable(data.frame(n=c(2,2,1),
    aic=c(norm=fitN$aic,
    uni=fitU$aic,
    exp=fitE$aic)))
# Exponential has lowest AIC

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitE) ; title(main="Exponential", outer=TRUE, line=-1)
# Exponential not too bad

# Do the plot:
palette("R3")
denscomp(list(fitN, fitE), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="AHMbook::MesoCarnivores HDens",
    xlab="log(Housing density)")

# Hubbard Brook slope
# ===================

data(HubbardBrook)
str(HubbardBrook)

x <- HubbardBrook$sitecov$Slope

descdist(x)
# Close to normal but skewed

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
# Weibull looks good

# Do the plot:
palette("R3")
denscomp(list(fitN, fitW), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="AHMbook::HubbardBrook slope",
    xlab="Slope, degrees")

# Hubbard Brook aspect
# ====================

data(HubbardBrook)
str(HubbardBrook)

x <- HubbardBrook$sitecov$Aspect

# A circular distribution
# Do a density plot
xRad <- x / 360 * 2*pi
overlap::densityPlot(xRad, xscale=360, extend=NULL, lwd=2,
    main="AHMbook::HubbardBrook aspect", xlab="Aspect, degrees")
# Strong modes around north and south


# Willow warbler GDD
# ==================

data(willowWarbler )
str(willowWarbler)

x <- willowWarbler$cells$gdd # mean growing degree days

# Plot gdd data
library(raster)
mapPalette <- colorRampPalette(c("gray", "yellow", "orange", "red"))
r1 <- with(willowWarbler$cells, rasterFromXYZ(data.frame(x = lon, y = lat, z = gdd)))
par(mar = c(1,1,3,0.5))
plot(r1, col = mapPalette(100), axes = FALSE, box = FALSE,
    main ="Map of GDD covariate")
box()

descdist(x)
# Close to uniform and normal, skewed

# Fit a range of distributions:
( fitU <- fitdist(x, "unif") )
( fitN <- fitdist(x, "norm") )
( fitL <- fitdist(x, "lnorm") )
( fitG <- fitdist(x, "gamma") )
( fitE <- fitdist(x, "exp") )     # failed
( fitW <- fitdist(x, "weibull") )
# All MLE fits work except exponential
AICtable(data.frame(n=c(2,2,2,2,2),
    aic=c(norm=fitN$aic,
    lnorm=fitL$aic,
    uni=fitU$aic,
    gam=fitG$aic,
    wei=fitW$aic)))
# Weibull has lowest AIC

# GOF plots:
plot(fitU) ; title(main="Uniform", outer=TRUE, line=-1)
plot(fitN) ; title(main="Normal", outer=TRUE, line=-1)
plot(fitL) ; title(main="Log-Normal", outer=TRUE, line=-1)
plot(fitG) ; title(main="Gamma", outer=TRUE, line=-1)
plot(fitW) ; title(main="Weibull", outer=TRUE, line=-1)
# None very good

# Do the plot:
palette("R3")
denscomp(list(fitN, fitW), fitlwd=2, datacol="gray90",
    demp=TRUE, xlegend = "topright",
    main="AHMbook::willowWarbler gdd",
    xlab="Mean growing degree days")

