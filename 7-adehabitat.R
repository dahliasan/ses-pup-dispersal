## Try adehabitat

rm(list= ls())

# Load packages -----------------------------------------------------------
library(tidyverse)
library(adehabitatHR)
source('convert2polarsf.R')


# Set common graphical elements ------------------------------------------
crs <- "+proj=stere +lon_0=170 +lat_0=-90 +units=m +datum=WGS84"
proj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# world map
world_sf <- st_as_sf(rworldmap::getMap(resolution = "li")) %>% st_transform(., crs = crs)

# mq island
mq <- data.frame(lon = 158.95, lat = -54.5)
mq_sf <- convert2polarsf(mq)


# Load data ---------------------------------------------------------------

load('./Output/foragingTracks_filtered.RData')
load("./Output/foragingTracks_with_survival.RData")
# raw <- read_csv('./Data/mq-ellie-weaners-argos.csv')

# seal <- loc %>% filter(id == "mq4-Alice-00") %>% 
#   mutate(x = lon, y = lat) %>% 
#   dplyr::select(id, date, x, y)

# KDE ---------------------------------------------------------------------
## adehabitat takes sp objects (SpatialPoints - one animal or SpatialPointsDataFrame - many animals)

## Set projection
loc_sf <- loc_s %>% 
  st_as_sf(coords = 7:8, crs = proj) %>%
  st_transform(crs = crs)
loc_sf$lon = st_coordinates(loc_sf)[, 1]
loc_sf$lat = st_coordinates(loc_sf)[, 2]
st_geometry(loc_sf) <- NULL

## Convert to a SpatialPointsDataFrame
d1 <- loc_sf %>% filter(trip == 1, SUS == FALSE, n_seen > 0) %>% dplyr::select(id, date, lon, lat)
coordinates(d1) = c("lon", "lat") # specify column names
class(d1)
head(d1)

## kde
kud <- kernelUD(d1[,1], h = "href", extent = 2)
image(kud)

## Estimate home range
hr <- getverticeshr(kud)
plot(hr, col=1:length(unique(d1$id)))

hr <- getverticeshr(kud, percent = 50)
plot(hr, col=1:length(unique(d1$id)))

elev <- puechabonsp$map
image(elev, 1)
plot(ver, add=TRUE, col=rainbow(4))
legend(699000, 3165000, legend = names(ud), fill = rainbow(4))

image(kud[[1]])
xyz <- as.image.SpatialGridDataFrame(kud[[1]])
contour(xyz, add=TRUE)

vud <- getvolumeUD(kud)
image(vud[[1]])
xyzv <- as.image.SpatialGridDataFrame(vud[[1]])
contour(xyzv, add=TRUE)

## Example plot of one animal
fud <- vud[[1]]
hr50 <- as.data.frame(fud)[,1]
hr50 <- as.numeric(hr50 <= 50)
hr50 <- as.data.frame(hr50)
coordinates(hr50) <- coordinates(vud[[1]])
gridded(hr50) <- TRUE
image(hr50)
map(add = T, fill=T, col="grey")

## Determining home range UD cell use for all animals (50% or 95%)

y <- estUDm2spixdf(vud)
fud_all <- as.data.frame(y@data)

for (i in 1:ncol(y@data)) {
  fudall[i] <- as.data.frame(y@data[i])
  fudall[i] <- as.numeric(y@data[i] <= 95)
}

coordinates(fudall) <- coordinates(y@coords)
gridded(fudall) <- TRUE
image(fudall[1])
