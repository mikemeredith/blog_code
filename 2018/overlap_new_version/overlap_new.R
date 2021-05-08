
# Code for the web page at
#  https://mmeredith.net/blog/2018/overlap_sunTime.htm

library(overlap)
data(simCalls)
str(simCalls)
# 'data.frame':   100 obs. of  2 variables:
 # $ time : num  2.38 4.07 1.17 2.26 1.86 ...
 # $ dates: chr  "2017-01-05" "2017-01-09" "2017-01-13" "2017-01-14" ...

Dates <- as.POSIXct(simCalls$dates, tz="GMT")

coords <- matrix(c(-3, 56), nrow=1)
Coords <- sp::SpatialPoints(coords,
   proj4string=sp::CRS("+proj=longlat +datum=WGS84"))

st <- sunTime(simCalls$time, Dates, Coords)

par(mfrow=2:1)
densityPlot(st, col='red', lwd=2, xaxt='n', main="Sun time")
axis(1, at=c(0, 6, 12, 18, 24),
    labels=c("midnight", "sunrise", "noon", "sunset", "midnight"))
densityPlot(simCalls$time, lwd=2, main="Clock time")
par(mfrow=c(1,1))
