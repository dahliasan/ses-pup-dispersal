#############################################
# Extract chla (copernicus) data to tracks
# Do this after extracting enviro data from raadtools 
# Output: add chl column to original raadtools dataset
#############################################
rm(list=ls(all=TRUE))

# Install packages
# install.packages(c("raster", "ncdf4", "ggplot2", "lubridate", "sf", "leaflet",
#                    "rasterVis", "rnaturalearth", "rnaturalearthdata"))

# Load package
library(ncdf4)
library(raster)
library(sf)
library(leaflet)
library(lubridate)
library(rasterVis)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)


#############################################
# Inspect NetCDF files
#############################################

# Set the path for the NetCDF file
ncfiles <- dir('./Data/chla-copernicus', full.names = T, pattern = '.nc')


#############################################
# Import NetCDF as Raster
#############################################

# import multi-band NetCDF file
chl <- stack(brick(ncfiles[1]),
             brick(ncfiles[2]),
             brick(ncfiles[3]))

# Load location dataset
load("./Output/raadtools_extract_sim_tracks_12h.rdata")

# Extract chla

ddates <- unique(d$day)
rdates <- as.POSIXct(names(chl), format = "X%Y.%m.%d.%H.%M.%S", tz = "UTC") %>% as.Date %>% sort()
r <- chl[[which(rdates %in% ddates)]]
rdates <- as.POSIXct(names(r), format = "X%Y.%m.%d.%H.%M.%S", tz = "UTC") %>% as.Date %>% sort()

d1 <- NULL
if(all.equal(ddates, rdates)) 
  for(i in 1:length(ddates)){
    df1 <- d[d$day==ddates[i],]
    r1 <- r[[i]]
    # extract chl
    df1$chl <- raster::extract(r1, df1[, c("lon", "lat")])
    # extract chlgrad
    df1$chlgrad <- raster::extract(terrain(r1, opt = 'slope'), df1[, c("lon", "lat")])
    # combine
    d1 <- rbind(d1, df1)
    print(paste(i,'of',length(ddates)))
  }

d <- d1
d <- d %>% arrange(id, sim, date)
# save(d, file = "./Output/raadtools_extract_sim_tracks_12h_chla.rdata")
beepr::beep(4)
