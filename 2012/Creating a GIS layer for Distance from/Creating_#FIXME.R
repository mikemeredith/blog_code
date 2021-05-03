# Code from
# https://mmeredith.net/blog/2012/1212_GIS_layer_for_Distance_from.htm
# downloaded 2021-05-02

library(raster)
library(rgdal)

# Produce a layer with distance-from-nearest-road
# ===============================================

# Read in the roads shape file and plot it.
# .v or .r on the end of the name indicate vector or raster.

roads.v <- readOGR("toy_example/roads.shp")
plot(roads.v)

LU.r <- raster("toy_example/LandUse.asc")
LU.r
# Note that 'crs' is NA
crs(LU.r) <- CRS("+init=epsg:3395")
plot(LU.r)

proj4string(LU.r)
proj4string(roads.v)  # Not the same

roads.v <- spTransform(roads.v, crs(LU.r))

lines(roads.v)

r <- LU.r  # this will be the template
r[] <- NA  # assigns all values as NA
summary(r) # shows you what you have: all NA's

# Now we use rasterize to create a roads raster matching the template.
roads.r <- rasterize(roads.v, r, field=1)
summary(roads.r)          # pixels crossed by a road have "1"
plot(roads.r, add=TRUE)

# Create the distance-to-road raster
roaddist.r <- distance(roads.r)
class(roaddist.r)
# Check:
plot(roaddist.r)
lines(roads.v)

writeRaster(roaddist.r, "DistFromRoad.tiff", "GTiff")


# Get distance-from-nearest-road for each of our cameras
# ======================================================

# Read in the cameras shape file and plot cameras and roads together to check:
cams.v <- readOGR("toy_example/cameras.shp")
proj4string(cams.v)  # Not same
cams.v <- spTransform(cams.v, crs(LU.r))

plot(roads.v)
points(cams.v)

# Now extract the distance-from-road information from roaddist.r for the camera locations:
distFromRoad <- extract(roaddist.r, cams.v)
distFromRoad

# You can add this information to the attribute table for the cameras and save it.
# We'll use the modern GeoPackage format.
cams2.v <- cams.v
cams2.v$roadDist <- distFromRoad
head(cams2.v)
writeOGR(cams2.v, dsn="toy_example/cameras.gpkg", layer="cameras2", driver="GPKG")

# Getting distance-from-water
# ===========================

water.r <- LU.r
water.r[LU.r != 5] <- NA
summary(water.r)
plot(water.r)

streams.v <- readOGR("toy_example/streams.shp")  ## FIXME
lines(streams.v) # not ok!

proj4string(streams.v) 
streams2.v <- spTransform(streams.v, crs(water.r))
lines(streams2.v)

# Create a raster from the stream layer just as we did for the roads earlier, and plot it to check:
stream2.r <- rasterize(streams2.v, r, field=1)
summary(stream2.r)
plot(water.r)
plot(stream2.r, add=TRUE)

water2.r <- water.r
water2.r[!is.na(stream2.r)] <- 1
summary(water2.r)
plot(water2.r, col=rainbow(5))

# The rest of the procedure is the same as for the distance-from-road information.

# For the distance-from-water layer, calculate distance from each pixel to the nearest non-NA pixel in water2.r:
waterdist.r <- distance(water2.r)
plot(waterdist.r)
lines(streams2.v, col='blue')
plot(water.r, add=TRUE, col='blue', legend=FALSE)

writeRaster(waterdist.r, "DistFromWater.tiff", "GTiff")

# To get the distance-from-water for each camera and add it to the cameras attribute table
distFromWater <- extract(waterdist.r, cams.v)
distFromWater

cams3.v <- cams2.v
cams2.v$waterDist <- distFromWater
head(cams3.v)

writeSpatialShape(cams3.v, "cameras3")

Load cameras3.shp into QGIS and check the attribute table.

These general methods could be extended to other distance-from-... variables, such as distance from dwellings (shown on a vector layer of points) or distance from forest edge (polygons), either outwards from the forest or inwards into the forest.

