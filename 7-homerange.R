# try amt 

rm(list= ls())

# Load packages -----------------------------------------------------------
library(tidyverse)
# devtools::install_github("jmsigner/amt")
library(amt)
source('convert2polarsf.R')


# Set common graphical elements ------------------------------------------
crs <- "+proj=stere +lon_0=170 +lat_0=-90 +units=km +datum=WGS84"
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


# Start -------------------------------------------------------------------

# Make track object
trk <- loc_s %>% 
  filter(trip == 1, SUS == FALSE) %>% # first trip only
  ungroup() %>% 
  make_track(lon, lat, date, id = id, crs = sp::CRS(proj), n_seen = n_seen)

# Transform projection
trk_tr <- transform_coords(trk, sp::CRS(crs))

# Make template raster for KDE
trast <- make_trast(trk_tr, res = 5)

# KDE
# hr <- hr_kde(trk_tr, trast = trast, levels = c(0.5, 0.9))
# plot(hr)


## KDE by individual
dat <- trk_tr %>% 
  nest(data = -id) %>% 
  mutate(kde = map(data, hr_kde, trast = trast, levels = c(0.5, 0.95)))

dat <- dat %>% 
  mutate(n = map_int(data, nrow))

# Convert hr to sf
hr_sf <- dat %>% hr_to_sf(kde, id, n)

# Plot
ggplot() + 
  geom_sf(data = hr_sf %>% filter(level == 0.95)) + 
  geom_sf(data = hr_sf %>% filter(level == 0.5), fill = 'red') + 
  facet_wrap(~id)



## KDE by survivor
dat <- trk_tr %>% 
  mutate(survive = n_seen > 1) %>% 
  nest(data = -survive) %>% 
  filter(!is.na(survive)) %>% 
  mutate(kde = map(data, hr_kde, trast = trast, levels = c(0.5, 0.95))) %>% 
  mutate(n = map_int(data, nrow))


# Convert hr to sf
hr_sf <- dat %>% hr_to_sf(kde, survive, n)
loc_sf <- convert2polarsf(loc_s)
bb <- extent(loc_sf) 

# Plot
png('UD50_survivors.png', width=25, height=20, units='cm',
    res=500)
ggplot() + 
  # geom_sf(data = hr_sf %>% filter(level == 0.95)) + 
  geom_sf(data = hr_sf %>% filter(level == 0.5), aes(fill = survive), alpha = 0.5) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  geom_sf(data = mq_sf, size = 0.5, colour = 'red') + 
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) + 
  labs(title = 'UD 50%')
dev.off()

ggplot() + 
  geom_sf(data = hr_sf %>% filter(level == 0.95), aes(fill = survive), alpha = 0.5) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  geom_sf(data = mq_sf, size = 0.5, colour = 'red') + 
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4])
