
# R code for the blog post at
#  https://mmeredith.net/blog/2019/plotting_rasters.htm

library(raster)

# Get nice examples to play with:
library(AHMbook)
data(BerneseOberland)
str(BerneseOberland)
rbo <- rasterFromXYZ(BerneseOberland[1:3])
plot(rbo)

library(unmarked)
data(Switzerland)
str(Switzerland)
rche <- rasterFromXYZ(Switzerland[c(1,2,3)]) # elevation
rchf <- rasterFromXYZ(Switzerland[c(1,2,4)]) # forest

# Remove axes
par(mfrow=1:2)
plot(rchf, col=rampGreens(), axes = FALSE, main="Forest cover")
plot(rche, axes = FALSE, main="Elevation")

# Reduce the margins (also manually resize the plot window)
par(mfrow=1:2, mar=c(0,0,3,1)+0.1)
plot(rchf, col=rampGreens(), axes = FALSE, main="Forest cover")
plot(rche, axes = FALSE, main="Elevation")

# Increase legend.mar, remove boxes
plot(rchf, col=rampGreens(), axes = FALSE, main="Forest cover", legend.mar=8.1, box=FALSE)
plot(rche, axes = FALSE, main="Elevation", legend.mar=8.1, box=FALSE)

# Horizontal legend bar
# ---------------------
par(mar=c(0,0,0,0)+0.1)
plot(rchf, col=rampGreens(256), axes = FALSE, horizontal = TRUE, box=FALSE)

# Zooming in
# ----------

par(mfrow=1:2, mar=c(0,0,3,0)+0.1)
plot(rchf, axes = FALSE, main="With ext", ext=extent(rbo), legend=FALSE)
plot(rchf, axes = FALSE, main="With xlim, ylim", colNA = 'black',
  xlim=c(910442,960442), ylim=c(53776,103776), legend=FALSE)

# Interpolation
# -------------

par(mfrow=1:2, mar=c(0,0,3,0)+0.1)
plot(rbo, axes = FALSE, main="Original", legend=FALSE)
plot(rbo, axes = FALSE, main="Interpolated", interpolate=TRUE, legend=FALSE)


# Plots with the same scale and a single legend
# ---------------------------------------------

par(mfrow=1:2, mar=c(0,0,3,1)+0.1)
plot(rche, axes = FALSE, main="Switzerland")
plot(rbo, axes = FALSE, main="Bernese Oberland")

# Same scale for both, going up to 4000+. Check the range for each
rche ; rbo
zlim <- range(values(rche), na.rm=TRUE)
plot(rche, axes = FALSE, main="Switzerland", zlim=zlim)
plot(rbo, axes = FALSE, main="Bernese Oberland", zlim=zlim)

# Just one legend
plot(rche, axes = FALSE, main="Switzerland", zlim=zlim, legend=FALSE, box=FALSE)
plot(rbo, axes = FALSE, main="Bernese Oberland", zlim=zlim, legend.mar=8.1, box=FALSE)

# Multiple plots, same legend
# ---------------------------
zlim <- range(values(rbo))
par(mfrow=1:2, mar=c(0,0,3,0.5)+0.1, oma=c(0,0,0,0))
plot(rbo, axes = FALSE, main="Before", zlim=zlim, legend=FALSE)
plot(rbo, axes = FALSE, main="After", zlim=zlim, legend.mar=8)

plot(rbo, axes = FALSE, main="Before", zlim=zlim, legend=FALSE)
plot(rbo, axes = FALSE, main="After", zlim=zlim, legend=FALSE)
plot(rbo, axes = FALSE, main="After", zlim=zlim, legend.only=TRUE)
# plot(rbo, axes = FALSE, main="After", zlim=zlim, legend.only=TRUE, legend.mar=8)

# Use outer margin, 'oma', to make space for the legend
par(mfrow=1:2, mar=c(0,0,3,0.5)+0.1, oma=c(0.5,0.5,0,2))
plot(rbo, axes = FALSE, main="Before", zlim=zlim, legend=FALSE)
plot(rbo, axes = FALSE, main="After", zlim=zlim, legend=FALSE)
plot(rbo, zlim=zlim, legend.only=TRUE)
