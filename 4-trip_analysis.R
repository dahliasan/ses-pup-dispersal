# descriptive stats – duration of the transmissions, total distance travelled,
# mean distance travelled, lat an lon of most distant point. Perhaps divide
# these into the complete tracks (those that return to MI) and incomplete. 
# ** trim any on-land locations before the seals head off to sea, otherwise
# these will bias the trip durations.


library(tidyverse)
library(sf)
library(SGAT)
library(raster)
library(lubridate)
library(foieGras)
library(viridis)
library(tmap)
library(ggpubr)
library(plotly)

## load data
load('./Output/weaner_ssm_4h.Rdata') #ssm locations
raw <- grab(fit_all, 'data', as_sf = F) #raw data
load('./Data/mq-ellie-weaners-haulout.RData') #haulout data
source('convert2polarsf.R')

haul <- haul %>% 
  mutate(startdate = S.DATE %>% mdy_hms(),
         enddate = E.DATE %>% mdy_hms())

## get names of seals with complete and incomplete tracks
id_keep <- dir('./Output/ssm plots/keep') %>% str_replace('.png', '')

ssm <- ssm %>% 
  mutate(type = 'never left', 
         type = replace(type, id %in% id_keep, 'keep'))

## get only seals that left MI
ssm1 <- ssm %>% filter(!type == 'never left')


# Set common graphical elements ------------------------------------------
crs <- "+proj=stere +lon_0=170 +lat_0=-90 +units=km +datum=WGS84"
proj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# world map
world_sf <- st_as_sf(rworldmap::getMap(resolution = "li")) %>% st_transform(., crs = crs)

# Explore starting locations ----------------------------------------------

dmp1 <- ssm1 %>% 
  group_by(id) %>% 
  summarise(lon = first(lon), lat = first(lat))

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

bb <- extent(ssmp_sf) # this bb will be used as the standard for individual plots later
ssm_plot <- ggplot() +
  geom_sf(data = ssmp_sf, colour = 'orange') +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4])

ssm_plot




# Step 1: trim locations when seal on land ----------------------------------------
## create sf object of macquarie island only
mqbb <- extent(c(158, 159, -55, -54)) %>% st_bbox()
mq_sf <- st_as_sf(rworldmap::getMap(resolution = "high")) %>% 
  st_crop(mqbb) %>% 
  dplyr::select(geometry)

# plot(mq_sf)

## create sf object for seal tracks
dmp1 <- ssm1 %>% dplyr::select(id, date, lon, lat,type)
coordinates(dmp1) <- c("lon", "lat")
projection(dmp1) <-  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
ssm_sf <- dmp1 %>% st_as_sf(dmp1)

## group into foraging trips
# find when locations intersect with macquarie island
ssm_sf$land <- st_intersects(ssm_sf, mq_sf, sparse = FALSE) %>% as.vector()
ssm_sf <- ssm_sf %>% 
  group_by(id) %>% 
  mutate(trip = cumsum(land)) 

## categorise complete and incomplete trips
a <- ssm_sf
st_geometry(a) <- NULL
b <- a %>% 
  group_by(id, trip) %>% 
  summarise(date = max(date) + 60*60*4) %>% # add 4 hours (sample) to each observation
  ungroup()

# get only locations when seal in on land
b <- left_join(b, ssm_sf %>% ungroup %>% dplyr::select(id, date, land)) 
b <- b %>% filter(land == TRUE) %>% dplyr::select(-geometry, -date) %>% rename(complete = land)

# get only locations when seal is at sea
tm <- left_join(ssm_sf, ssm1 %>% dplyr::select(id, date, lon, lat)) %>% 
  filter(land == FALSE) 

tm <- left_join(tm, b)
tm <- tm %>% mutate(complete = ifelse(is.na(complete), FALSE, TRUE), 
                    type = ifelse(complete == TRUE, 'complete', 'incomplete'))

plot(tm)

# Step 2: trip metrics ------------------------------------------------------------
# duration of the transmissions, total distance travelled,
# mean distance travelled, lat an lon of most distant point.
mqcol <- cbind(158.95, -54.5)
tm1 <- tm %>%
  group_by(id, trip) %>% 
  mutate(dur = difftime(max(date), min(date), units = 'days'),
         dist2next = gcDist(cbind(lon, lat), cbind(lead(lon), lead(lat))),
         dist2col = gcDist(mqcol, cbind(lon, lat)),
         lon.max = lon[dist2col == max(dist2col)],
         lat.max = lat[dist2col == max(dist2col)])


## get trip summary
trips_all <- tm1 %>% 
  group_by(id, trip) %>% 
  summarise(type = first(type),
            startdate = min(date), 
            enddate = max(date),
            dur = mean(dur), 
            totaldist = sum(dist2next, na.rm = TRUE), 
            maxdist2col = max(dist2col), 
            lon.max = mean(lon.max), 
            lat.max = mean(lat.max))

trips1 <- trips_all %>% 
  filter(dur > 5) %>% 
  mutate(dur = dur %>% as.numeric)
  
## summary metrics
summ1 <- trips1 %>% 
  group_by(id) %>% 
  mutate(n= n()) %>% 
  group_by(type) %>% 
  summarise(nseals = unique(id) %>% length(), 
            ntrips = n(),
            meanntrips = mean(n), 
            meandist = mean(totaldist), 
            maxdist = max(totaldist),
            meandur = mean(dur), 
            maxdur = max(dur)) 
  

## summary: transmission duration
summ2 <- ssm1 %>% 
  group_by(id) %>% 
  summarise(startdate = min(date), enddate = max(date),
            dur = difftime(enddate, startdate, units = 'days') %>% as.numeric())  


## save summaries
# save(trips1, file = './Output/trips.Rdata')
# write_csv(trips1 %>% mutate_if(is.numeric, ~signif(.,3)), file = './Output/trip_analysis.csv')
# write_csv(summ1 %>% mutate_if(is.numeric, ~signif(.,3)), file = './Output/incomplete_vs_complete_trips_summary.csv')
# write_csv(summ2 %>% mutate_if(is.numeric, ~signif(.,3)), file = './Output/seal_transmission_duration.csv')

## Example plot of 1 seal
seal1 <- trips1 %>% filter(id == 'mq3-22488-99')
seal1 <- seal1 %>% st_transform(crs = crs)
ggplot() +
  geom_sf(data = seal1, aes(colour = trip))


summ1 %>% 
  st_transform(crs = crs) %>% 
  ggplot() + 
  geom_sf(aes(colour = type))
  

# Update: 3 March 2021 (Wed) ----------------------------------------------
# Identify individual trips -----------------------------------------------
## Identify on land locations
# Create sf object of MI only
mqbb <- extent(c(158, 159, -55, -54)) %>% st_bbox()
mq_sf <- st_as_sf(rworldmap::getMap(resolution = "high")) %>% 
  st_crop(mqbb) %>% 
  dplyr::select(geometry)

# plot(mq_sf)

# Create sf object for seal tracks
dmp1 <- ssm1 %>% dplyr::select(id, date, lon, lat,type)
coordinates(dmp1) <- c("lon", "lat")
projection(dmp1) <-  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
ssm_sf <- dmp1 %>% st_as_sf(dmp1)

# Find seal locations that intersect with MI
ssm_sf$land <- st_intersects(ssm_sf, mq_sf, sparse = FALSE) %>% as.vector()
ssm_sf <- ssm_sf %>% 
  group_by(id) %>% 
  mutate(trip = cumsum(land)) 

trip <- ssm_sf 
st_geometry(trip) <- NULL
trip <- trip %>% left_join(ssm1 %>% dplyr::select(id, date, lon, lat))

# / Distance from colony over time ------------------------------------------
mqcol <- cbind(158.95, -54.5) #MI location

## Calculate distance to colony 
trip <- trip %>% mutate(dist2col = gcDist(mqcol, cbind(lon, lat)))
raw1 <- raw %>% filter(keep == TRUE & id %in% trip$id ) %>% mutate(dist2col = gcDist(mqcol, cbind(lon, lat)))

## Calculate days since departure
trip <- trip %>% 
  group_by(id) %>% 
  mutate(departDays = difftime(date, min(date), units = 'days') %>% as.numeric())

raw1 <- raw1 %>% 
  group_by(id) %>% 
  mutate(departDays = difftime(date, min(date), units = 'days') %>% as.numeric())

## Plot distance to colony over time
# png('dist2col_all.png', width=20, height=16, units='cm',
#     res=500)
# trip %>% 
#   ggplot(aes(x = departDays, y = dist2col)) + 
#   geom_point(size = 0.1, colour = 'black') + 
#   geom_point(data = raw1, size = 0.1, colour = 'red') + 
#   facet_wrap(~id) + 
#   # scale_x_datetime(breaks = scales::breaks_pretty(3)) + 
#   labs(title = 'Dist to Macquarie Is. Black = ssm, red = raw', y = 'Dist2col (km)', x = 'Days after departure') +
#   theme(text = element_text(size = 7),
#         axis.text.x = element_text( size = 4)) 
# dev.off()

## Plot raw vs ssm tracks
# trip_sf <- convert2polarsf(trip)
# raw1_sf <- raw1 %>% dplyr::select(id, date, lc, lon, lat) %>% convert2polarsf()
# bb <- extent(dat1_sf)

# png('tracks_ssm_vs_raw_all.png', width=30, height=20, units='cm',
#     res=500)
# ggplot() +
#   geom_sf(data = dat1_sf, size = 0.05, colour = 'black') +
#   geom_sf(data = raw1_sf, size = 0.05, colour = 'red') +
#   geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
#   xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) + 
#   facet_wrap(~id) + 
#   theme(text = element_text(size = 8)) +
#   labs(title = 'black = ssm locations, red = raw locations')
# dev.off()

# / Explore trip duration threshold ---------------------------------------
# reset trip id to start from 1 = 1st trip, 2 = 2nd trip etc
tmp <- trip %>% 
  group_by(id, trip) %>% 
  summarise(n = n()) %>% 
  filter(n > 1) 

tmp %>% ggplot(aes(n)) + geom_histogram(binwidth = 6) 


## Check histogram of # trip locations to determine threshold
# png('unfiltered_trip_dur_histogram_all.png', width=20, height=16, units='cm',
#     res=500)
# tmp %>% ggplot(aes(n)) +
#   geom_histogram(binwidth = 10, fill = 'red') +
#   facet_wrap(~id) + 
#   labs(title = '# of 4-h locations per trip (filtering out short trips) - binwidth = 10 locations') + 
#   theme_pubr(border = TRUE) + 
#   theme(text = element_text(size = 7))
# 
# dev.off()


# / Filter real foraging  trips ------------------------------------------
## Criteria for real foraging trip = > 4 days 4d*24h/4h

trip1 <- trip %>% 
  group_by(id, trip) %>% 
  summarise(maxdist = max(dist2col), n = n()) %>% 
  filter(n > 1) #remove land locations

trip1 %>%
  ggplot(aes(maxdist)) + 
  geom_histogram(binwidth = 20) #40-60km threshold

trip1 %>%
  ggplot(aes(n)) + 
  geom_histogram(binwidth = 6) #6 locations = 24h = 1 day, 3 day threshold

trip1 <- trip1 %>% 
  filter(maxdist > 60, n > 24) %>% 
  mutate(trip1 = 1) %>% 
  group_by(id) %>% 
  mutate(trip1 = cumsum(trip1))

# tmp <- tmp %>% 
#   filter(n > 34) %>% #set trip duration threshold
#   mutate(trip1 = 1) %>% 
#   group_by(id) %>% 
#   mutate(trip1 = cumsum(trip1))

trip2 <- left_join(ssm_sf, trip1) %>% 
  mutate(trip1 = factor(trip1))

## Plot
trip2_sf <- trip2 %>% st_transform(crs = crs)
bb <- extent(trip2_sf)

trip2_sf <- trip2_sf %>% mutate(dataset = str_split(id, '-')[[1]][1])
# png('trips_all_threshold=60km,4d.png', width=30, height=20, units='cm',
#     res=500)
# ggplot() +
#   geom_sf(data = trip2_sf, aes(colour = trip1), size = 0.0001) +
#   geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
#   xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) + 
#   facet_wrap(~id) + 
#   scale_color_viridis_d() + 
#   theme(text = element_text(size = 8))
# dev.off()

 ## clean dataset for saving
loc <- trip2_sf %>% 
  mutate(trip = ifelse(is.na(trip1), 0, trip1)) %>% 
  dplyr::select(-type, -trip1, -land)

st_geometry(loc) <- NULL

loc <- left_join(loc, ssm1 %>% dplyr::select(id, date, lon, lat))

# save(loc, file = './Output/foraging_tracks_filtered.RData')


# Check seals with long gaps ----------------------------------------------

load('./Output/tracks_filtered.RData')

seals_sus <- c('mq3-26627-99', 'mq3-2846-99', 'mq3-28496-99', 'mq4-FirstOne-00')
loc_sus <- loc %>% 
  filter(id %in% seals_sus)

raw_sus <- raw1 %>% filter(id %in% seals_sus)

# check interval between each location
loc_sus %>% group_by(id) %>% mutate(interval = date - lag(date)) %>% pull(interval) %>% unique
raw_sus <- raw_sus %>% group_by(id) %>% mutate(interval = date - lag(date), interval = interval %>% as.numeric())


p1 <- raw_sus %>% 
  ggplot(aes(x = date, y = 1)) + 
  geom_point(size = 0.5) +  
  facet_wrap(~id)

ggplotly(p1) 

loc_sus_sf <- loc_sus %>% convert2polarsf()
raw_sus_sf <- raw_sus %>% dplyr::select(id, date, lc, lon, lat) %>% convert2polarsf()
  
bb <- extent(raw_sus_sf)

p1 <- ggplot() +
  geom_sf(data = loc_sus_sf, size = 0.05, colour = 'black') +
  geom_sf(data = raw_sus_sf, size = 0.05, colour = 'red', aes(text = date)) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) +
  facet_wrap(~id)

ggplotly(p1)


# Spurious locations
# mq3-28496-99    > 2000-04-19 17:29:10
# mq4-FirstOne-00    > 2001-01-30 16:26:14
# mq3-2846-99    1999-12-26 20:20:30 – 2000-03-01 10:49:51

loc$SUS <- FALSE
loc$SUS[loc$id == 'mq3-28496-99' & loc$date > ymd_hms('2000-04-19 17:29:10')] <- TRUE
loc$SUS[loc$id == 'mq4-FirstOne-00' & loc$date > ymd_hms('2001-01-30 16:26:14')] <- TRUE
loc$SUS[loc$id == 'mq3-2846-99' & 
          loc$date > ymd_hms('1999-12-26 20:20:30') & 
          loc$date < ymd_hms('2000-03-01 10:49:51')] <- TRUE

# save(loc, file = './Output/tracks_filtered.RData')

          # Haulout -----------------------------------------------------------------
# haul1 <- haul %>% 
#   mutate(haulout.num = HAULOUT.NUMBER) %>% 
#   dplyr::select(ref, haulout.num, lat, lon, startdate, enddate) %>% 
#   filter(ref %in% id_keep) %>% 
#   mutate(hauldur = difftime(enddate, startdate, units = 'days') %>% as.numeric())
#   
# haul1 %>% 
#   ggplot(aes(x = hauldur)) +
#   geom_histogram(binwidth = 1/24)
# 
# ## Plot haulouts
# haul1 %>% 
#   ggplot(aes(x = lon, y = lat)) + 
#   geom_point(aes(colour = hauldur)) + 
#   facet_wrap(~ref, scales = 'free')
# 
# 
# dmp1 <- haul1 %>% filter(!is.na(lon))
# 
# coordinates(dmp1) <- c("lon", "lat")
# projection(dmp1) <-  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# ssmp <- dmp1 %>%
#   st_as_sf(dmp1)
# 
# ## project data to Antarctic Stereographic
# crs <- "+proj=stere +lon_0=170 +lat_0=-90 +units=km +datum=WGS84"
# ssmp_sf <- ssmp %>% st_transform(crs = crs)
# 
# ## get high-res land & apply stereo projection
# world_sf <- st_as_sf(rworldmap::getMap(resolution = "high")) %>%
#   st_transform(., crs = crs)
# 
# bb <- extent(ssmp_sf)
# ssm_plot <- ggplot() +
#   geom_sf(data = ssmp_sf, aes(colour = hauldur)) +
#   geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
#   xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) +
#   scale_color_viridis()
#   
# ssm_plot


