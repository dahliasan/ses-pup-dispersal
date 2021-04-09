## Try ctmm

rm(list= ls())

# Load packages -----------------------------------------------------------

library(tidyverse)
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
raw <- read_csv('./Data/mq-ellie-weaners-argos.csv')

# Pre-processing ----------------------------------------------------------

seal <- loc %>% filter(id == "mq4-Alice-00") %>% 
  ungroup() %>% 
  mutate(timestamp = as.character(date), 
         individual.local.identifier = id, 
         location.long = lon, 
         location.lat = lat) %>% 
  select(individual.local.identifier, timestamp, location.long, location.lat)

## Convert to telemetry object
seal_tel <- as.telemetry(seal, timeformat = '%Y-%m-%d %H:%M:%S', projection = crs)


## Calculate variogram
vg <- variogram(seal_tel)

## Plot up to 50% of the maximum lag in the data
plot(vg)

## Zoom in on the shortest lags
plot(vg, fraction=0.005)

## Use the sliders provided by variogram.fit to specify starting values.
## The default choices are usually acceptable.
variogram.fit(vg)

## Automatically fit the range-resident models via maximum likelihood
## using the initial parameter values obtained from variogram.fit()
fitted.mods <- ctmm.select(seal_tel, CTMM=GUESS, verbose=TRUE)
summary(fitted.mods)

## Examine model results
ouf <- fitted.mods$`OUF anisotropic`
plot(vg, CTMM = ouf, col.CTMM="#1b9e77")
plot(vg, CTMM = ouf, col.CTMM="#1b9e77", fraction=0.005)

summary(ouf)

# Autocorrelated KDE estimate
akde <- akde(seal_tel, CTMM = ouf)

summary(akde)
plot(seal_tel, UD = akde)
