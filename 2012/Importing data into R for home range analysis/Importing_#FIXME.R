# Code from
# https://mmeredith.net/blog/1212_Data_for_home_range_analysis_in_R.htm
# downloaded 2021-05-02


# 1. Data formats needed by adehabitatHR
# --------------------------------------

library(adehabitatHR)
data(puechabonsp)
str(puechabonsp, max=2)

# 2. Choosing your Coordinate Reference System (CRS) and resolution (no code)

# 3. Importing and formatting your tracking data

# track <- read.csv("data/tracking_WGS84.csv")
track <- read.csv("data/tracking_WGS84.csv", stringsAsFactors = TRUE)  ## FIXME
summary(track)
head(track)

library(raster)
library(rgdal)
library(maptools)

coordinates(track) <- c("X", "Y")
class(track)
plot(track)

plot(track, col=track$Name)

# track has no CRS information:
proj4string(track)

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

proj4string(herb) <- CRS("+init=epsg:27573")

# 5. Importing shape files and converting to rasters

# legal.v <- readShapeSpatial("Legal_status.shp")  ## FIXME
legal.v <- readOGR("data/Legal_status.shp", stringsAsFactors=TRUE)
plot(legal.v)
summary(legal.v)

# proj4string(legal.v) <- CRS("+init=epsg:3857") # FIXME no longer needed

Now we transform this layer to our chosen CRS and check that it matches the herb raster layer:
legal2.v <- spTransform(legal.v, CRS("+init=epsg:27573"))
plot(herb)
plot(legal2.v, add=TRUE) # looks good.

head(legal2.v)
 # Type
# 0 A
# 1 B
# 2 C
# 3 X
# 4 B
# 5 C

legal2.v$typeCode <- as.numeric(legal2.v$Type)
head(legal2.v)
  Type typeCode
0    A        1
1    B        2
2    C        3
3    X        4
4    B        2
5    C        3

legal2.r <- rasterize(legal2.v, herb, field="typeCode")
summary(legal2.r)
# Min. 1.000
# 1st Qu. 1.000
# Median 2.000
# Mean 1.849
# 3rd Qu. 2.000
# Max. 4.000
# NA's 4982.000

legal2.r[is.na(legal2.r[])] <- 0
summary(legal2.r)
table(legal2.r[])
   # 0     1     2     3     4
# 4982   526   853   274     1

# 6. Downloading data from the internet, and extracting the data you want.

?getData

elev.big <- getData('SRTM', lon=4, lat=44)
elev.mid <- crop(elev.big, extent(3.5, 3.7, 43.65, 43.8))
writeRaster(elev.mid, "elevation_mid")

# This loads elev.mid from the data folder:
elev.mid <- raster("data/elevation_mid")

elev.mid
# class : RasterLayer
# dimensions : 180, 240, 43200 (nrow, ncol, ncell)
# resolution : 0.0008333333, 0.0008333333 (x, y)
# extent : 3.5, 3.7, 43.65, 43.8 (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0
# data source : C:/Documents and Settings/Mike/My Documents/aadocs/Blog_bits/Data for adehabitatHR/data/elevation_mid.grd
# names : srtm_37_04
# values : 32, 814 (min, max)

plot(elev.mid)

herb.padded <- extend(herb, 1)
herb
herb.padded

elev <- projectRaster(elev.mid, herb.padded)
plot(elev)
elev

# 7. Calculating derived layers: getting slope and aspect from an elevation layer.
# --------------------------------------------------------------------------------

?terrain
slope <- terrain(elev, opt='slope', unit='degrees')
summary(slope)
plot(slope)

aspect.deg <- terrain(elev, opt='aspect', unit='degrees')
summary(aspect.deg)

?reclassify

# To use this, we need the rcl matrix with columns "from", "to", and "becomes":
rcl <- matrix(c(
  0,  45, 1,
 45, 135, 2,
135, 225, 3,
225, 315, 4,
315, 360, 1), ncol=3, byrow=TRUE)
rcl
aspect <- reclassify(aspect.deg, rcl=rcl)
summary(aspect)
plot(aspect)

# 8. Putting the habitat rasters together and formatting for adehabitatHR
# -----------------------------------------------------------------------
# we put all our raster layers together in a rasterStack:
habitat.rs <- stack(elev, slope, aspect, herb.padded, legal2.r)
# Error in compareRaster(x) : Different extent

# Oops! I forgot that legal2.r matches herb, not herb.padded; we'll pad it the same way:
legal.padded <- extend(legal2.r, 1)
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

### FIXME need to rewrite!!!

cp <- mcp(my.homerange.data$relocs[,1], percent=95)

cprast <- hr.rast(cp, my.homerange.data$map)

# image(my.homerange.data$map[,5])
names(my.homerange.data$map)  # Check the layer names!
image(my.homerange.data$map[, "herbaceous"])
image(cprast[,2], add=TRUE, col="black", useRasterImage=FALSE)


