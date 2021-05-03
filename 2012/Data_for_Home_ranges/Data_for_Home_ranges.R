# Code from
# https://mmeredith.net/blog/1212_Data_for_home_range_analysis_in_R.htm
# downloaded 2021-05-02


# 1. Data formats needed by adehabitatHR
# --------------------------------------

library(adehabitatHR)
data(puechabonsp)
str(puechabonsp, max=2)
# List of 2
#  $ map   :Formal class 'SpatialPixelsDataFrame' [package "sp"] with 7 slots
#  $ relocs:Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots

# 2. Choosing your Coordinate Reference System (CRS) and resolution (no code)

# 3. Importing and formatting your tracking data
# ----------------------------------------------

track <- read.csv("data/tracking_WGS84.csv", stringsAsFactors = TRUE)
summary(track)
head(track)
#   Name Age Sex   Date        X        Y
# 1 Chou   3   2 920729 3.571935 43.72938
# 2 Chou   3   2 920802 3.573262 43.72154
# 3 Chou   3   2 920803 3.573899 43.72293
# 4 Chou   3   2 920804 3.573661 43.72281
# 5 Chou   3   2 920805 3.575234 43.72204
# 6 Chou   3   2 920809 3.574405 43.72215

library(raster)
library(rgdal)
# library(maptools)

coordinates(track) <- c("X", "Y")
class(track)
# [1] "SpatialPointsDataFrame"
plot(track, col=track$Name)

# track has no CRS information:
proj4string(track)
# NA

proj4string(track) <- CRS("+init=epsg:4326")

track <- spTransform(track, CRS("+init=epsg:27573"))
summary(track)
# Object of class SpatialPointsDataFrame
# Coordinates:
# min max
# X 698626 701410
# Y 3157848 3161678
# Is projected: TRUE
# proj4string :
# [+init=epsg:27573 +proj=lcc +lat_1=44.10000000000001  ## FIXME
# [+proj=lcc +lat_1=44.10000000000001
# +lat_0=44.10000000000001 +lon_0=0 +k_0=0.999877499 +x_0=600000
# +y_0=3200000 +a=6378249.2 +b=6356515 +towgs84=-168,-60,320,0,0,0,0
# +pm=paris +units=m +no_defs]
# Number of points: 119
# Data attributes:
# ...

# 4. Importing raster layers

herb <- raster("data/herbaceous.asc")
plot(herb)
herb
# class      : RasterLayer 
# dimensions : 84, 79, 6636  (nrow, ncol, ncell)
# resolution : 100, 100  (x, y)
# extent     : 697850, 705750, 3156750, 3165150  (xmin, xmax, ymin, ymax)
# crs        : NA 
# source     : C:/.../data/herbaceous.asc 
# names      : herbaceous 
# values     : -2147483648, 2147483647  (min, max)


proj4string(herb) <- CRS("+init=epsg:27573")

# 5. Importing shape files and converting to rasters
# --------------------------------------------------

legal.v <- readOGR("data/Legal_status.shp", stringsAsFactors=TRUE)
plot(legal.v)
summary(legal.v)
# Object of class SpatialPolygonsDataFrame
# Coordinates:
#         min       max
# x  395030.7  403156.1
# y 5420247.3 5427317.2
# Is projected: TRUE 
# proj4string :
# [+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1
# +units=m +nadgrids=@null +wktext +no_defs +type=crs]
# Data attributes:
 # Type 
 # A:1  
 # B:2  
 # C:2  
 # X:1  

legal.v <- spTransform(legal.v, crs(herb))
plot(herb)
plot(legal.v, add=TRUE) # looks good.

head(legal.v)
 # Type
# 0 A
# 1 B
# 2 C
# 3 X
# 4 B
# 5 C

legal.v$typeCode <- as.numeric(legal.v$Type)
head(legal.v)
#   Type typeCode
# 0    A        1
# 1    B        2
# 2    C        3
# 3    X        4
# 4    B        2
# 5    C        3

legal.r <- rasterize(legal.v, herb, field="typeCode")
summary(legal.r)
#         layer
# Min.        1
# 1st Qu.     1
# Median      2
# 3rd Qu.     2
# Max.        4
# NA's     4982

legal.r[is.na(legal.r[])] <- 0
table(legal.r[])
#    0     1     2     3     4
# 4982   526   853   274     1

plot(legal.r)

# 6. Downloading data from the internet, and extracting the data you want.
# ------------------------------------------------------------------------

?getData

elev.big <- getData('SRTM', lon=4, lat=44)
elev.mid <- crop(elev.big, extent(3.5, 3.7, 43.65, 43.8))

# This loads elev.mid from the data folder:
elev.mid <- raster("data/elevation_mid")

elev.mid
# class      : RasterLayer 
# dimensions : 180, 240, 43200  (nrow, ncol, ncell)
# resolution : 0.0008333333, 0.0008333333  (x, y)
# extent     : 3.5, 3.7, 43.65, 43.8  (xmin, xmax, ymin, ymax)
# crs        : +proj=longlat +datum=WGS84 +no_defs 
# source     : data/elevation_mid.grd 
# names      : srtm_37_04 
# values     : 32, 814  (min, max)


plot(elev.mid)

herb.padded <- extend(herb, 1)
herb
# dimensions : 84, 79, 6636  (nrow, ncol, ncell)
herb.padded
# dimensions : 86, 81, 6966  (nrow, ncol, ncell)

elev <- projectRaster(elev.mid, herb.padded)
plot(elev)

# 7. Calculating derived layers: getting slope and aspect from an elevation layer.
# --------------------------------------------------------------------------------

?terrain
slope <- terrain(elev, opt='slope', unit='degrees')
summary(slope)
#                slope
# Min.      0.08353635
# 1st Qu.   3.44000787
# Median    7.60628675
# 3rd Qu.  14.94807435
# Max.     40.11498634
# NA's    330.00000000
plot(slope)

aspect.deg <- terrain(elev, opt='aspect', unit='degrees')
summary(aspect.deg)
#               aspect
# Min.      0.09894173
# 1st Qu. 124.12124560
# Median  178.17135391
# 3rd Qu. 255.01922499
# Max.    359.70694748
# NA's    330.00000000

?reclassify

# To use this, we need the rcl matrix with columns "from", "to", and "becomes":
rcl <- matrix(c(
    0,  45, 1,
   45, 135, 2,
  135, 225, 3,
  225, 315, 4,
  315, 360, 1), ncol=3, byrow=TRUE)
rcl
#      [,1] [,2] [,3]
# [1,]    0   45    1
# [2,]   45  135    2
# [3,]  135  225    3
# [4,]  225  315    4
# [5,]  315  360    1

aspect <- reclassify(aspect.deg, rcl=rcl)
summary(aspect)
#         aspect
# Min.         1
# 1st Qu.      2
# Median       3
# 3rd Qu.      4
# Max.         4
# NA's       330
plot(aspect)

# 8. Putting the habitat rasters together and formatting for adehabitatHR
# -----------------------------------------------------------------------
# we put all our raster layers together in a rasterStack:
habitat.rs <- stack(elev, slope, aspect, herb.padded, legal.r)
# Error in compareRaster(x) : Different extents

# Oops! I forgot that legal2.r matches herb, not herb.padded; we'll pad it the same way:
legal.padded <- extend(legal.r, 1)
habitat.rs <- stack(elev, slope, aspect, herb.padded, legal.padded)
plot(habitat.rs)

names(habitat.rs) <- c("elevation", "slope", "aspect", "herbaceous", "legalStatus")
plot(habitat.rs)

habitat.rs <- mask(habitat.rs, herb.padded)
plot(habitat.rs)

habitat <- as(habitat.rs, "SpatialPixelsDataFrame")
mimage(habitat)

bbox(habitat)
      # min     max
# x  697850  705750
# y 3156750 3165150
bbox(herb)
       # min     max
# s1  697850  705750
# s2 3156750 3165150

# Finally...
# ----------

proj4string(track) == proj4string(habitat)
# [1] TRUE

my.homerange.data <- list(map = habitat, relocs = track)

# 9. Try it out!
# --------------

cp <- mcp(my.homerange.data$relocs[,1], percent=95)

cprast <- hr.rast(cp, my.homerange.data$map)

# image(my.homerange.data$map[,5])
names(my.homerange.data$map)  # Check the layer names!
image(my.homerange.data$map[, "herbaceous"])
image(cprast[,2], add=TRUE, col="black", useRasterImage=FALSE)
