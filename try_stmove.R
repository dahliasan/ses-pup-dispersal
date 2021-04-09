## Trying out stmove r package

rm(list = ls())

# Load packages -----------------------------------------------------------

library(tidyverse)
# install.packages("remotes")
# remotes::install_github("dpseidel/stmove")
library(stmove)
library(ctmm)
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


# Pre-processing ----------------------------------------------------------
## convert to correct projection
# loc_sf <- convert2polarsf(loc)

loc_sf <- loc %>%
  st_as_sf(coords = 7:8, crs = proj) %>%
  st_transform(crs = crs)


loc_sf %>% 
  ggplot() + 
  geom_sf()

## stmove requires locations to be at regular intervals, in this case, output from foiegras already made it so at 4h intervals
## check for NAs in location
table(is.na(loc$lat))

## check # suspicious locations
table(loc$SUS)

loc_final <- loc_sf

loc_final$x = st_coordinates(loc_sf)[, 1]
loc_final$y = st_coordinates(loc_sf)[, 2]
loc_final <- loc_final %>% 
  dplyr::select(x, y, date, id)
st_geometry(loc_final) <- NULL

# Build report ------------------------------------------------------------

seals <- loc_final$id %>% unique

build_report(loc_final %>% filter(id == "mq4-Alice-00"),
             path = '~/OneDrive - University of Tasmania/Elephant-Weaners',
             proj4 = crs,
             stats = c("rolling", "diurnal", "lunar"),
             construct = "akde")


