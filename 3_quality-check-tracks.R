# Quality check weaner tracks

# Objectives:
# 1)  Plot each seal individually to check the tracks are sensible – i.e. that
# the seal actually left the island, or doesn’t have weird gaps.
# 2) save each of these maps for future reference

## load libraries
library(tidyverse)
library(lubridate)
library(foieGras)
library(sf)
library(raster)
library(viridis)

## load ssm filtered tracks
load("./Output/weaner_ssm_4h.Rdata")

## extract the predicted values - ie estimated location every day - in lat/lon format
ssm <- grab(fit_all, "predicted", as_sf=FALSE)
dmp1 <-  ssm

## plot the ssm data
coordinates(dmp1) <- c("lon", "lat")
projection(dmp1) <-  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

ssmp <- dmp1 %>%
  st_as_sf(dmp1)

## project data to Antarctic Stereographic
crs <- "+proj=stere +lon_0=170 +lat_0=-90 +units=km +datum=WGS84"
ssmp_sf <- ssmp %>% st_transform(crs = crs)


## get high-res land & apply stereo projection
world_sf <- st_as_sf(rworldmap::getMap(resolution = "high")) %>%
  st_transform(., crs = crs)

bb <- extent(ssmp_sf)
ssm_plot <- ggplot() +
  geom_sf(data = ssmp_sf, size=.1) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  xlim(bb[1]-1000, bb[2]+500) + ylim(bb[3]-1000, bb[4]+1000)

ssm_plot

## plot each seal individually
by_id <- ssmp_sf %>% split(ssmp$id)

ssm_plots <- purrr::map(by_id, function(df){
  
  # create equidistant sequence of dates to use as labels
  lab_dates <- pretty(df$date)
  
  plot <- ggplot() +
    geom_sf(data = df, size=.1, aes(colour = date)) +
    geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
    # xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) +
    xlim(bb[1]-1000, bb[2]+500) + ylim(bb[3]-1000, bb[4]+1000) +
    labs(title = df$id[1]) + 
    scale_colour_viridis(breaks = as.numeric(lab_dates), 
                         labels = lab_dates %>% as.character())
  
  plot
  
})
ssm_plots[[3]]

ids <- names(ssm_plots)

for(i in ids){
  png(paste('./Output/ssm plots/', i, '.png', sep = ''), units = 'in', width = 4, height = 3, res = 240)
  print(ssm_plots[i])
  dev.off()
  
}

table(ssmp_sf$id)