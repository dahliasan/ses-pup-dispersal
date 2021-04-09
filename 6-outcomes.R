## Comparing first trips to outcome metrics (i.e. body condition, survival)

rm(list = ls())
library(tidyverse)
library(sf)
library(SGAT)
library(raster)
library(lubridate)
library(foieGras)
library(viridis)
library(tmap)
library(ggpubr)
library(see)



# Load Datasets -----------------------------------------------------------
load('./Output/allFit_prefilter=0.3.RData')
load('./Output/foragingtracks_filtered.RData')
source('convert2polarsf.R')


# Set common graphical elements ------------------------------------------
crs <- "+proj=stere +lon_0=170 +lat_0=-90 +units=km +datum=WGS84"
proj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# world map
world_sf <- st_as_sf(rworldmap::getMap(resolution = "li")) %>% st_transform(., crs = crs)

# mq island
mq <- data.frame(lon = 158.95, lat = -54.5)
mq_sf <- convert2polarsf(mq)

# Drift Dives ---------------------------------------------------------------
## Extract drift dives from augFit
drift <- purrr::map(allFit, function(x) {
  xx <- x$augFit
  if(xx %>% length > 1) xx$data
  
  })

drift <- drift %>% 
  reduce(bind_rows) %>%
  dplyr::select(ref, DE.DATE, SURF.DUR, DIVE.DUR, MAX.DEP, drate) %>% 
  rename(date = DE.DATE, surfdur = SURF.DUR, divedur = DIVE.DUR, maxdep = MAX.DEP) %>% 
  as_tibble()

save(drift, file = './Output/driftRate_augFit_prefilter=0.3.RData')

## Get daily average drift rates
drift1 <- drift %>%
  mutate(day = as.Date(date), ref = ref %>% as.character()) %>%
  group_by(ref, day) %>%
  summarise(drate.avg = median(drate), drate.sd = sd(drate), drate.max = max(drate)) %>%
  rename(id = ref)

## Get daily average drift rates
drift1 <- drift %>%
  mutate(day = as.Date(date), ref = ref %>% as.character()) %>%
  group_by(ref, day) %>%
  summarise(drate.avg = median(drate), drate.sd = sd(drate), drate.max = max(drate)) %>%
  rename(id = ref)

# ## Split drift1 dates into 4 hour intervals to match locations
# drift1 <- drift %>% 
#   mutate(date1 = cut(date, '4 hours') %>% ymd_hms(), ref = ref %>% as.character()) %>% 
#   rename(id = ref)

## Combine location and drift rates datasets
d <- left_join(loc %>% mutate(day = date %>% as.Date), drift1) 
d <- d %>% filter(id %in% drift1$id)

## Calculate daily drift rate change
tmp <- d %>% 
  dplyr::select(id, day, trip, drate.avg) %>% 
  unique() %>% 
  group_by(id, trip) %>% 
  mutate(dRateAvgLog = log(drate.avg + abs(min(d$drate.avg, na.rm = T)) + 0.001),
         deltaDRate = drate.avg - lag(drate.avg), deltaDRatePerc = deltaDRate/lag(drate.avg))

d <- left_join(d, tmp)

## Plot tracks with drift rates
d_sf <- convert2polarsf(d)
bb <- extent(d_sf) 

png('tracks_augFitpf0.3_drates.png', width=25, height=20, units='cm',
    res=500)
ggplot() +
  geom_sf(data = d_sf %>% filter(trip ==1), size = 0.1, colour = 'lightgrey') +
  geom_sf(data = d_sf %>% filter(trip ==1, !is.na(drate.max)), aes(colour = drate.max), size = 0.2) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  geom_sf(data = mq_sf, size = 0.5, colour = 'red') + 
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) + 
  facet_wrap(~id) +
  scale_colour_viridis(name = 'daily drift rate', na.value = 'lightgrey') +
  # scale_colour_distiller(palette = 'RdYlBu') + 
  # theme(text = element_text(size = 8)) +
  theme_pubr(border = TRUE, base_size = 8)
dev.off()

## Daily change in drift rate
png('tracks_augFitpf0.3_deltadrates.png', width=25, height=20, units='cm',
    res=500)
ggplot() +
  geom_sf(data = d_sf %>% filter(trip == 1), size = 0.1, colour = 'darkgrey', alpha = 0.1) +
  geom_sf(data = d_sf %>% filter(trip == 1, !is.na(deltaDRate)), aes(colour = deltaDRate), size = 0.05) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  geom_sf(data = mq_sf, size = 0.5, colour = 'red') + 
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) + 
  facet_wrap(~id) +
  scale_colour_viridis(name = 'daily change drift rate', na.value = 'lightgrey') +
  # scale_colour_gradient2(high = '#EB5300', low = '#01F4EB') +
  # theme(text = element_text(size = 8)) +
  theme_pubr(border = TRUE, base_size = 8) + 
  theme_abyss()
dev.off()

# Survival Data -----------------------------------------------------------

s <- read_csv('./Data/CMC/TAGS/6mth_Macca 180810.csv', col_types = cols('WeaningLength' = 'd', 
                                                                        'WeaningGirth' = 'd'))

s <- s %>% select(-Index)
colnames(s) <- colnames(s) %>% tolower()
s <- s %>%
  mutate(birthdate = birthdate %>% dmy(), weandate = weandate %>% dmy()) %>% 
  rename(seal_ID = sealid)

s1 <- s %>% gather("year_season", "seen", 2:25)

## Create separate year and season columns
s1$season <- 'summer'
s1$season[str_detect(s1$year_season, 'w')] <- 'winter'
s1$year <- s1$year_season %>% str_split(pattern = '[:alpha:]') %>% map(~.x[1]) %>% unlist()

## Summarise survival data...
# some metrics to consider: What is the year/season seal was last seen?
s2 <- s1 %>% 
  group_by(seal_ID, birthdate, birthyear, birthlocation, weandate, sex, blackmass, weanmass) %>% 
  summarise(n_seen = sum(seen, na.rm = T))

s2 %>% 
  ggplot(aes(x = n_seen)) + 
  geom_histogram(binwidth = 1) # majority only seen once (i.e. when they were tagged)


## Match seal_ID to seal track df
# Load seal IDs matching
load('./Output/seal_ID.RData')
loc_s <- left_join(loc, sealID_all %>% dplyr::select(ref, SEAL_ID) %>% rename(id = ref, seal_ID = SEAL_ID))
loc_s <- loc_s %>% left_join(s2) # ***WHY SOME SEALS NO RESIGHT DATA AT ALL?

## Plot tracks colour coded to n_seen
loc_s_sf <- convert2polarsf(loc_s)
bb <- extent(loc_s_sf) 

png('tracks_seen_again_1map.png', width=25, height=20, units='cm',
res=500)
ggplot() +
  geom_sf(data = loc_s_sf %>% filter(trip ==1, !is.na(n_seen), SUS == FALSE), aes(colour = n_seen > 1), size = 0.1) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  geom_sf(data = mq_sf, size = 0.5, colour = 'red') + 
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) + 
  labs(title = 'Survival: if seal was seen again after the first time')
  # facet_wrap(~id) +
  # facet_wrap(~n_seen > 1) +
  # scale_colour_viridis() +
  # scale_colour_distiller(palette = 'RdYlBu') + 
  theme(text = element_text(size = 8)) 
dev.off()

save(loc_s, file = "./Output/foragingTracks_with_survival.RData")

# Combine DD and Survival Datasets ----------------------------------------

dim(d)
dim(loc_s)

loc1 <- left_join(loc_s, d)

save(loc1, file = './Output/foragingtracks_drift_survival_combined.RData')



