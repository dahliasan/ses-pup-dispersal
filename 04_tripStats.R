# descriptive stats – duration of the transmissions, total distance travelled,
# mean distance travelled, lat an lon of most distant point. Perhaps divide
# these into the complete tracks (those that return to MI) and incomplete. 
# ** trim any on-land locations before the seals head off to sea, otherwise
# these will bias the trip durations.
rm(list = ls())
options(tibble.width = Inf)


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
library(rworldmap)
sf::sf_use_s2(FALSE)
load("baseInfo.Rdata")


# Load data ---------------------------------------------------------------
# load('./Output/weaner_ssm_4h.Rdata') # ssm locations (loads ssm + fit_all)
# load('./Output/weaner_ssm_24h.Rdata') # ssm locations (loads ssm + fit_all)
# raw <- grab(fit_all, 'data', as_sf = F) # raw data
fitssm <- readRDS("./Output/foiegras_fitssm_vmax4_12h.rds")
ssm <- grab(fitssm, 'predicted', as_sf = F)
raw <- grab(fitssm, 'data', as_sf = F) # raw data
load('./Data/mq-ellie-weaners-haulout.RData') # haulout data (loads haul)
load("./Data/mq-ellie-weaners-dive.RData") # dive data
source('convert2polarsf.R')


ssm1 <- ssm

# Visualise haulout data --------------------------------------------------
## visualise haulout data
haul <- haul %>% 
  filter(ref %in% ssm1$id) %>% 
  mutate(startdate = S.DATE %>% mdy_hms(),
         enddate = E.DATE %>% mdy_hms(), 
         cohort = substr(ref, 1,3),
         ref = as.character(ref))

haul %>% 
  filter(cohort == 'mq1') %>% 
  ggplot(aes(y = ref)) + 
  geom_segment(aes(x = startdate, xend = enddate, yend = ref), 
               size = 3, 
               colour = 'red',
               linetype = 2) +
  facet_wrap(~cohort, scales = 'free') +
  scale_x_datetime(date_breaks = '1 month', date_labels = '%m/%y')
  

# Set common graphical elements ------------------------------------------
crs <- "+proj=stere +lon_0=170 +lat_0=-90 +units=km +datum=WGS84"
proj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# world map
world_sf <- st_as_sf(rworldmap::getMap(resolution = "li")) %>% 
  st_transform(., crs = crs)

# # Explore starting locations ----------------------------------------------
# 
# dmp1 <- ssm1 %>% 
#   group_by(id) %>% 
#   summarise(lon = first(lon), lat = first(lat))
# 
# ## plot the ssm data
# coordinates(dmp1) <- c("lon", "lat")
# projection(dmp1) <-  proj
# ssmp <- dmp1 %>% st_as_sf(dmp1)
# 
# ## project data to Antarctic Stereographic
# ssmp_sf <- ssmp %>% st_transform(crs = crs)
# 
# bb <- extent(ssmp_sf)
# ssm_plot <- ggplot() +
#   geom_sf(data = ssmp_sf, colour = 'red') +
#   geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
#   xlim(bb[1], bb[2]) + ylim(bb[3], bb[4])
# 
# ssm_plot # 1 outlier location


# Identify trips ----------------------------------------------------------

# + Find land locations ---------------------------------------------------
## Defined as within x km to island, cross-checked with dive and haulout data

BUFFER_DIST <- 50 # buffer in km. Chose 50km same as mcconell et al 2002

# Create sf object of MI only
mqbb <- extent(c(158, 159, -55, -54)) %>% st_bbox()
mq_sf <- st_as_sf(rworldmap::getMap(resolution = "high")) %>% 
  st_crop(mqbb) %>% 
  dplyr::select(geometry) %>% 
  st_transform(crs = crs)
  
# Create buffer around MI sf
mq_sfbuffer <- st_buffer(mq_sf, BUFFER_DIST) 

plot(mq_sfbuffer)
plot(mq_sf, add = TRUE)

# Create sf object for seal tracks
dmp1 <- ssm1 %>% dplyr::select(id, date, lon, lat)
coordinates(dmp1) <- c("lon", "lat")
projection(dmp1) <-  proj
ssm_sf <- dmp1 %>% st_as_sf(dmp1) %>% st_transform(crs=crs)

# Find seal locations that intersect with MI
ssm_sf$land <- st_intersects(ssm_sf, mq_sfbuffer, sparse = FALSE) %>% as.vector()

# Identify trips
ssm_sf <- ssm_sf %>% 
  group_by(id) %>% 
  mutate(trip = cumsum(land))

plot(mq_sfbuffer)
plot(mq_sf, add = TRUE)
plot(ssm_sf %>% filter(land == TRUE), add = TRUE)

trip <- ssm_sf 
# st_geometry(trip) <- NULL
trip <- trip %>% left_join(ssm1 %>% dplyr::select(id, date, lon, lat))

## Calculate days since tag deployment (to standardise plotting scale later)
trip <- trip %>% 
  group_by(id) %>% 
  mutate(daysFromDeployment = difftime(date, min(date), units = 'days') %>% as.numeric())

# + Match haulouts to locations -------------------------------------------
i <- survival::neardate(haul$ref, trip$id, as.Date(haul$startdate), as.Date(trip$date))
i <- ifelse(abs(haul$startdate - trip$date[i]) > 86400, NA, i)
trip$haulout <- NA
trip$haulout[i] <- TRUE

ggplot() +
  geom_sf(data = trip, size = 0.3) +
  geom_sf(data = trip %>% filter(haulout == TRUE), colour = 'red', size = 0.3) + 
  facet_wrap(~id)


# + Calculate dist from colony --------------------------------------------
mqcol <- cbind(158.95, -54.5) #MI location
trip <- trip %>% mutate(dist2col = gcDist(mqcol, cbind(lon, lat)))

# Visualise dist2col over time + colour haulouts
# trip %>% 
#   ggplot(aes(x = daysFromDeployment, y = dist2col)) + 
#   geom_path() +
#   geom_point(data = trip %>% filter(haulout == TRUE), colour = 'red') + 
#   facet_wrap(~id)


# + Filter real foraging trips (NEW) --------------------------------------------
trip1 <- trip
trip1 <- trip1 %>% 
  group_by(id, trip) %>% 
  mutate(tripdur = difftime(last(date), first(date), units = 'days'))


# Min trip duration threshold (in days)
TRIPDUR_THRESHOLD <- 1

# Visualise dist2col over time + colour haulouts
trip1 %>% 
  filter(tripdur > TRIPDUR_THRESHOLD) %>% 
  ggplot(aes(x = daysFromDeployment, y = dist2col)) + 
  geom_path(aes(colour = factor(trip))) +
  geom_point(data = trip %>% filter(haulout == TRUE), colour = 'red') +
  facet_wrap(~id)


# Create new trip numbers based on actual foraging trips + 
# identify completed trips
tmp <- trip1 
st_geometry(tmp) <- NULL
tmp <- tmp %>% 
  group_by(id, trip) %>% 
  summarise(trip2 = 1, tripdur = first(tripdur)) %>% 
  group_by(id) %>% 
  mutate(is_trip_complete = ifelse(is.na(trip - lead(trip)), FALSE, TRUE)) %>% 
  filter(tripdur > TRIPDUR_THRESHOLD) %>%
  mutate(trip2 = cumsum(trip2)) %>% 
  dplyr::select(-tripdur)

trip1 <- trip1 %>% left_join(tmp)

bb <- extent(trip1)
ggplot() +
  geom_sf(data = trip1, aes(colour = factor(trip2)), size = 0.1) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) +
  facet_wrap(~id) +
  theme(text = element_text(size = 8))


beepr::beep(4)
# Min trip dur threshold set to > 1 day, 
# However for mq4-Billie-00's first trip is ignored since it's complete and it
# is clear it wasn't the real first long trip. 
trip1$trip2[trip1$id == 'mq4-Billie-00'] <- trip1$trip2[trip1$id == 'mq4-Billie-00'] - 1
trip1$trip2[which(trip1$trip2 == 0)] <- NA

# png('allTrips.png', width=20, height=16, units='cm', res=500)
ggplot() +
  geom_sf(data = trip1, aes(colour = factor(trip2)), size = 0.01) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) +
  facet_wrap(~id) +
  theme(text = element_text(size = 8))
# dev.off()


## Visual check if trip completed info is correct
# Check first trips vs if trip completed vs haulout data
# png('firstTrips_ifcomplete.png', width=24, height=16, units='cm', res=500)
trip1 %>% 
  filter(trip2 == 1) %>% 
  ggplot() +
  geom_sf(size = 0.05, aes(colour = is_trip_complete)) +
  geom_sf(data = trip1 %>% filter(trip2 == 1, haulout == TRUE), colour = 'blue', size = 0.05) + 
  labs(title = "blue = haul out locations") +
  facet_wrap(~id) +
  theme(text = element_text(size = 10))  
# dev.off()

# png('dist2col_allTrips_ifcomplete.png', width=24, height=16, units='cm', res=500)
trip1 %>% 
  filter(!is.na(is_trip_complete)) %>% 
  ggplot(aes(x = daysFromDeployment, y = dist2col)) + 
  geom_path(aes(colour = is_trip_complete)) +
  geom_point(data = trip %>% filter(haulout == TRUE), colour = 'blue', size = 0.05) + 
  facet_wrap(~id)+
  labs(title = "blue = haul out locations") +
  theme(text = element_text(size = 10))  
# dev.off()


## clean dataset for saving
out <- trip1 %>%
  dplyr::select(-trip) %>% 
  rename(trip = trip2)

saveRDS(out, file = './Output/tracks_processed_12h.rds')
# save(loc, file = './Output/tracks_processed.RData')


# Plot raw vs ssm tracks ------------------------------------------------
raw1 <- raw %>% 
  filter(keep == TRUE) %>% 
  mutate(dist2col = gcDist(mq, cbind(lon, lat)))

raw1 <- raw1 %>% 
  group_by(id) %>% 
  mutate(daysFromDeployment = difftime(date, min(date), units = 'days') %>% as.numeric())

## Plot raw vs ssm tracks
raw1_sf <- raw1 %>% dplyr::select(id, date, lc, lon, lat) %>% convert2polarsf()
bb <- extent(trip)

png('tracks_ssm_vs_raw_all.png', width=30, height=20, units='cm', res=500)
ggplot() +
  geom_sf(data = trip, size = 0.05, colour = 'black') +
  geom_sf(data = raw1_sf, size = 0.05, colour = 'red') +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) +
  facet_wrap(~id) +
  theme(text = element_text(size = 8)) +
  labs(title = 'black = ssm locations, red = raw locations')
dev.off()

# Check suspicious locations ----------------------------------------------
# loc <- readRDS("./Output/tracks_processed_12h.rds")

seals_sus <- c('mq3-22484-99',
               'mq3-2846-99', 
               'mq3-28496-99',
               'mq3-28494-99')

loc_sus <- loc %>% 
  filter(id %in% seals_sus)

raw <- grab(fitssm, 'data', as_sf = F) # raw data
raw_sus <- raw %>% filter(id %in% seals_sus)

# check interval between each location
loc_sus %>% 
  group_by(id) %>% 
  mutate(interval = date - lag(date)) %>% 
  pull(interval) %>% #in hours
  unique

raw_sus <- raw_sus %>% 
  group_by(id) %>% 
  mutate(interval = date - lag(date), 
         interval = interval %>% as.numeric())


p1 <- raw_sus %>% 
  ggplot(aes(x = date, y = 1)) + 
  geom_point(size = 0.5) +  
  facet_wrap(~id)

ggplotly(p1) 

loc_sus_sf <- loc_sus
raw_sus_sf <- raw_sus %>% dplyr::select(id, date, lc, lon, lat) %>% convert2polarsf()

bb <- extent(loc_sus_sf)

p1 <- ggplot() +
  geom_sf(data = loc_sus_sf, size = 0.1, colour = 'black') +
  geom_sf(data = raw_sus_sf, size = 0.1, colour = 'red', aes(text = date)) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) +
  facet_wrap(~id)

ggplotly(p1)


# Spurious locations
# mq1-5811-95 1995-12-25 04:43:08 - 1996-01-05 14:17:45
# mq3-28496-99    2000-04-19 17:29:10 - 2000-11-05 03:24:04
# mq3-2846-99    1999-12-26 20:20:30 – 2000-03-01 10:49:51

loc$SUS <- FALSE

# loc$SUS[loc$id == 'mq1-5811-95' & 
#           loc$date > ymd_hms('1995-12-25 04:43:08') &
#           loc$date < ymd_hms('1996-01-05 14:17:45')] <- TRUE

loc$SUS[loc$id == 'mq3-28494-99' & 
          loc$date > ymd_hms('2000-01-18 16:01:53')] <- TRUE
loc$SUS[loc$id == 'mq3-28496-99' & 
          loc$date > ymd_hms('2000-04-19 17:29:10')] <- TRUE
loc$SUS[loc$id == 'mq3-2846-99' & 
          loc$date > ymd_hms('1999-12-26 20:20:30')] <- TRUE

# visual check
bb <- extent(loc)
ggplot() +
  geom_sf(data = loc, size = 0.05, aes(colour = SUS)) +
  geom_sf(data = raw1_sf, size = 0.05, colour = 'red') +
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) +
  scale_colour_manual(values = c('black', 'blue')) + 
  facet_wrap(~id) +
  theme(text = element_text(size = 8))
beepr::beep(4)
# saveRDS(loc, file = './Output/tracks_processed_12h.rds')
# 





