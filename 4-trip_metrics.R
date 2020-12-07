# descriptive stats – duration of the transmissions, total distance travelled,
# mean distance travelled, lat an lon of most distant point. Perhaps divide
# these into the complete tracks (those that return to MI) and incomplete. 
# ** trim any on-land locations before the seals head off to sea, otherwise
# these’ll bias the trip durations.


library(tidyverse)
library(sf)
library(SGAT)

## load data
load('./Output/weaner_ssm_4h.Rdata')

## get names of seals with complete and incomplete tracks
id_keep <- dir('./Output/ssm plots/keep') %>% str_replace('.png', '')

ssm <- ssm %>% 
  mutate(type = 'never left', 
         type = replace(type, id %in% id_keep, 'keep'))

## filter only complete and incomplete tracks
ssm1 <- ssm %>% filter(!type == 'never left')

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

## create sf object for seal tracks
dmp1 <- ssm1 %>% dplyr::select(id, date, lon, lat,type)
coordinates(dmp1) <- c("lon", "lat")
projection(dmp1) <-  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
ssm_sf <- dmp1 %>% st_as_sf(dmp1)

## group into foraging trips
ssm_sf$land <- st_intersects(ssm_sf, mq_sf, sparse = FALSE) %>% as.vector()
ssm_sf <- ssm_sf %>% 
  group_by(id) %>% 
  mutate(trip = cumsum(land)) 

## categorise complete and incomplete trips
a <- ssm_sf
st_geometry(a) <- NULL
b <- a %>% 
  group_by(id, trip) %>% 
  summarise(date = max(date) + 60*60*4) %>% 
  ungroup()

b <- left_join(b, ssm_sf %>% ungroup %>% dplyr::select(id, date, land)) 
b <- b %>% filter(land == TRUE) %>% dplyr::select(-geometry, -date) %>% rename(complete = land)

tm <- left_join(ssm_sf, ssm1 %>% dplyr::select(id, date, lon, lat)) %>% 
  filter(land == FALSE)

tm <- left_join(tm, b)
tm <- tm %>% mutate(complete = ifelse(is.na(complete), FALSE, TRUE), 
                    type = ifelse(complete == TRUE, 'complete', 'incomplete'))

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
save(trips1, file = './Output/trips.Rdata')
write_csv(trips1 %>% mutate_if(is.numeric, ~signif(.,3)), file = './Output/trip_analysis.csv')
write_csv(summ1 %>% mutate_if(is.numeric, ~signif(.,3)), file = './Output/incomplete_vs_complete_trips_summary.csv')
write_csv(summ2 %>% mutate_if(is.numeric, ~signif(.,3)), file = './Output/seal_transmission_duration.csv')

 ## Example plot of 1 seal
seal1 <- trips1 %>% filter(id == 'mq3-22488-99')
seal1 <- seal1 %>% st_transform(crs = crs)
ggplot() +
  geom_sf(data = seal1, aes(colour = trip))


summ1 %>% 
  st_transform(crs = crs) %>% 
  ggplot() + 
  geom_sf(aes(colour = type))
  