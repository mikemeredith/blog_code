# Code from
# https://mmeredith.net/blog/2012/1212_GIS_layer_for_Distance_from.htm
# downloaded 2021-05-02

library(raster)
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal)

# Produce a layer with distance-from-nearest-road
# ===============================================

# Read in the roads shape file and plot it.
# .v or .r on the end of the name indicate vector or raster.

roads.v <- readOGR("toy_example/roads.shp")
plot(roads.v)

LU.r <- raster("toy_example/LandUse.asc")
LU.r
# class      : RasterLayer 
# dimensions : 330, 520, 171600  (nrow, ncol, ncell)
# resolution : 1000, 1000  (x, y)
# extent     : -230000, 290000, -2e+05, 130000  (xmin, xmax, ymin, ymax)
# crs        : NA 
# source     : C:/.../toy_example/LandUse.asc 
# names      : LandUse 
# values     : 1, 5  (min, max)



# Note that 'crs' is NA
crs(LU.r) <- CRS("+init=epsg:3395")

proj4string(LU.r)
# "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
proj4string(roads.v)  # Not the same
# "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs +type=crs"

roads.v <- spTransform(roads.v, crs(LU.r))
plot(LU.r)
lines(roads.v)

# Now we use rasterize to create a roads raster matching the land use raster.
roads.r <- rasterize(roads.v, LU.r, field=1)
summary(roads.r)          # pixels crossed by a road have "1"
plot(roads.r, add=TRUE, legend=FALSE)

# Create the distance-to-road raster
system.time(
roaddist.r <- distance(roads.r) )
plot(roaddist.r)
lines(roads.v)

writeRaster(roaddist.r, "DistFromRoad.tiff", "GTiff")


# Get distance-from-nearest-road for each of our cameras
# ======================================================

# Read in the cameras shape file and check the CRS
cams.v <- readOGR("toy_example/cameras.shp")
proj4string(cams.v)  # Not same
cams.v <- spTransform(cams.v, crs(LU.r))

# Add cameras to the plot:
points(cams.v, pch=16)

# Now extract the distance-from-road information from roaddist.r for the camera locations:
distFromRoad <- extract(roaddist.r, cams.v)
cams.v$roadDist <- distFromRoad
head(cams.v)
#   id    brand height  roadDist
# 1  1 cantrack     10  27658.63
# 2  2 cantrack     12  79000.00
# 3  3 cantrack      9 118270.88
# 4  4    cutty      7  51971.15
# 5  5    cutty     14  80056.23
# 6  6    cutty     13  71000.00

# Getting distance-from-water
# ===========================

water.r <- LU.r
water.r[LU.r != 5] <- NA
summary(water.r)
plot(water.r)

streams.v <- readOGR("toy_example/streams.shp")
lines(streams.v) # not ok!

proj4string(streams.v) 
streams.v <- spTransform(streams.v, crs(LU.r))
lines(streams.v)

# Create a raster from the stream layer just as we did for the roads earlier, and plot it to check:
stream.r <- rasterize(streams.v, LU.r, field=1)
plot(water.r, col='blue')
plot(stream.r, add=TRUE, col='blue')

water.r[!is.na(stream.r)] <- 1
summary(water.r)
#         LandUse
# Min.          1
# 1st Qu.       5
# Median        5
# 3rd Qu.       5
# Max.          5
# NA's     166736
plot(water.r, col='blue')

# The rest of the procedure is the same as for the distance-from-road information.

# For the distance-from-water layer, calculate distance from each pixel to the nearest non-NA pixel in water.r:
waterdist.r <- distance(water.r)
plot(waterdist.r)
plot(water.r, add=TRUE, col='blue', legend=FALSE)
writeRaster(waterdist.r, "DistFromWater.tiff", "GTiff")

# To get the distance-from-water for each camera and add it to the cameras attribute table
distFromWater <- extract(waterdist.r, cams.v)
cams.v$waterDist <- distFromWater
head(cams.v)
#   id    brand height  roadDist waterDist
# 1  1 cantrack     10  27658.63  52611.79
# 2  2 cantrack     12  79000.00  13416.41
# 3  3 cantrack      9 118270.88  11045.36
# 4  4    cutty      7  51971.15  40311.29
# 5  5    cutty     14  80056.23  31384.71
# 6  6    cutty     13  71000.00  35902.65

writeOGR(cams.v,dsn="toy_example/cameras.gpkg", layer="cameras", driver="GPKG")
