## Matching and tabulating survival data
## Output: cleaned survival data

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
library(see)
library(orsifronts)
library(pals)
library(lme4)
library(nlme)


# Load Datasets -----------------------------------------------------------
load('./Output/seal_ID.RData') # all the different types of seal IDs
load('./Output/allFit_prefilter=0.3.RData')
# load('./Output/tracks_processed.RData')
loc <- readRDS("./Output/tracks_processed_12h.rds")
load('baseInfo.RData')
source('convert2polarsf.R')

# Create master reference object for study seals' IDs ---------------------
st_geometry(loc) <- NULL
sseals <- loc
sseals <- left_join(sseals %>% dplyr::select(id), 
                    sealID_all %>% dplyr::select(ref, SEAL_ID), 
                    by = c('id' = 'ref')) %>% 
  rename(seal_id = SEAL_ID) %>% 
  filter(!duplicated(id)) %>%
  filter(id %in% id_keep)

loc <- loc %>% filter(id %in% sseals$id) %>% left_join(sseals)

# 1) Clive's 6 mth Survival Data -----------------------------------------------------------

# Clive's 6mth survival matrix
ssix <- read_csv('./Data/CMC/TAGS/6mth_Macca 180810.csv', 
              col_types = cols('WeaningLength' = 'd', 'WeaningGirth' = 'd'))

## Clean data - transpose data
ssix <- ssix %>% dplyr::select(-Index)
colnames(ssix) <- colnames(ssix) %>% tolower()
ssix <- ssix %>%
  mutate(birthdate = birthdate %>% dmy(), weandate = weandate %>% dmy()) %>% 
  rename(seal_id = sealid)

ssix <- ssix %>% gather("year_season", "is_seen", 2:25)

## Create separate year and season columns
ssix$season <- 'summer'
ssix$season[str_detect(ssix$year_season, 'w')] <- 'winter'
ssix$year <- ssix$year_season %>% str_split(pattern = '[:alpha:]') %>% purrr::map(~.x[1]) %>% unlist()

## Summarise survival data
# How many times was each seal seen at 6 month intervals?
ssix1 <- ssix %>%
  group_by(seal_id, birthdate, birthyear, birthlocation, weandate, sex, blackmass, weanmass, weaninglength, weaninggirth, lactationperiod) %>%
  filter(seal_id %in% sseals$seal_id) %>% 
  summarise(n_6mth_seen = sum(is_seen, na.rm = T),
            trip1_survive = ifelse(n_6mth_seen > 1, T, F),
            year1_survive = ifelse(n_6mth_seen > 2, T, F),
            last_seen_year = max(year[which(is_seen == 1)] %>% as.numeric(), na.rm = T))

# ssix1 %>%
#   ggplot(aes(x = n_6mth_seen)) +
#   geom_histogram(binwidth = 1) # majority only seen once (i.e. when they were tagged)

# ## Match seal_ID to seal track df
# loc_s <- loc %>% left_join(ssix1) # *** 3 seals with no 6mth survival data (26605, 28775, 28783)
# # table(loc_s$seal_id, loc_s$n_6mth_seen > 0)
# 
# # IDs of seals with no 6mth survival data
# loc_s %>% filter(last_seen_year %>% is.na()) %>% pull(id) %>% unique()


# 1) All sightings -----------------------------------------------

# Load all sightings data
sightings <- read_csv('./Data/CMC/TAGS/TAGS_ELLIE_SIGHTINGS.csv') %>%
  mutate(DATE_OBSERVED = mdy_hms(DATE_OBSERVED) %>% date)
colnames(sightings) <- colnames(sightings) %>% tolower()
sightings <- sightings %>% 
  left_join(sseals) %>% 
  filter(id %in% sseals$id) %>% 
  dplyr::select(id, seal_id, date_observed, weight, standard_length, girth)  %>% 
  mutate(date_observed = mdy_hms(date_observed) %>% date) 

# Truncate all sightings to 6 month periods using all sightings dataset
# winter = march - aug
sightings_6mon <- sightings %>% 
  mutate(year = year(date_observed), 
         month = month(date_observed), 
         season = ifelse(month > 2 & month < 9, 'winter', 'summer')) %>% 
  group_by(id, seal_id, year, season) %>% 
  summarise()

# Summarise survival history
survival <- sightings_6mon %>% 
  group_by(id, seal_id) %>% 
  summarise(n_6m_seen = n(), 
            seen_6m = ifelse(n_6m_seen > 1, T, F),
            seen_1y = ifelse(n_6m_seen > 2, T, F),
            last_seen_year = max(year))

# Add other life history metadata from Clive's 6mth survival dataset
survival <- left_join(survival, 
                      ssix1 %>% dplyr::select(seal_id, blackmass, weanmass, birthlocation, birthdate, weandate, lactationperiod, weaninglength, weaninggirth, birthyear, sex))

# Combine study seal all sightings to their survival data
sightings1 <- sightings %>% 
  left_join(survival) %>% 
  filter(!duplicated(.)) %>% 
  mutate(birthyear = min(year(date_observed)), last_seen_age = last_seen_year - birthyear)

sightings1_summary <- sightings1 %>% 
  group_by(id) %>% 
  summarise(weanmass = first(weanmass),  last_seen_age = first(last_seen_age))

sightings1_summary %>% 
  ggplot(aes(x = weanmass, y = last_seen_age)) + 
  geom_point() + 
  geom_smooth()
  

# Model: age ~ weanmass ---------------------------------------------------

m1 <- glm(last_seen_age ~ weanmass, family = poisson(), data = sightings1_summary)
summary(m1)
performance::model_performance(m1)
plot(m1)
effects::allEffects(m1) %>% plot


hist(sightings1_summary$last_seen_age, breaks = seq(0, max(sightings1_summary$last_seen_age), 1))


save(sightings1, survival, file = "./Output/survival.RData")
save.image("11_getSurvival_WORKSPACE.RData")


# + plot all sightings timeline -------------------------------------------
# Plot timeline of sightings
winter_dates <- tibble(startdate = seq(ymd(19950301), ymd(20090301), by = '1 year'),
                       enddate = seq(ymd(19950831), ymd(20090831), by = '1 year'))
winter_dates$group <- seq_len(nrow(winter_dates))

# png('all_sightings_timeline.png', width=30, height=20, units='cm', res=500)

ggplot() +
  geom_rect(data = winter_dates,
            aes(xmin = startdate,
                xmax = enddate,
                group = group,
                ymin = min(sightings1$id),
                ymax = max(sightings1$id)),
            fill = 'lightblue', alpha = 0.4) +
  geom_point(data = sightings1, aes(x = date_observed, y = id, colour = seen_6m)) +
  scale_x_date(date_breaks = '12 months', date_labels = '%y') +
  labs(title = 'All sightings (1995-2010). Colour = If seal was seen again after 6 months? Blue shade = winter months')+
  theme_linedraw()

# dev.off()




# Plot all measured weights
# (p1 <- sightings %>% 
#   filter(!is.na(weight)) %>% 
#   # mutate(date = ifelse(date_observed %>% month > 8, 
#   #                      format(date_observed, '0000-%m-%d'), 
#   #                      format(date_observed, '0001-%m-%d')),
#     #        date = date %>% ymd) %>% 
#     ggplot(aes(x = date_observed, y = weight)) + 
#     geom_point() +
#     scale_x_date(date_labels = '%m/%y') +
#     facet_wrap(~paste(id, seal_id, sep = ' | '), scales = 'free_x') +
#     theme_linedraw(base_size = 7) + 
#     theme(panel.spacing.x = unit(1, "lines")))
# 
# png('all_weighings.png', width=30, height=20, units='cm', res=500)
# p1
# dev.off()


# Reference to clive's definition of weaning mass
# s %>% dplyr::select(seal_ID, blackmass, weanmass, lactationperiod, weaninglength, weaninggirth) %>% filter(seal_ID %in% loc_s$seal_ID) %>% View()


# 3) All summary (deployment / trip / survival) --------------------------------------------

# Load raw argos trips to get full deployment info
fit_all <- readRDS("./Output/foiegras_fitssm_vmax4_12h.rds") #ssm locations
loc_raw <- grab(fit_all, 'data', as_sf = F) #raw data

# Get deployment start and end dates
deploy <- loc_raw %>% 
  filter(id %in% loc$id) %>% 
  group_by(id) %>% 
  summarise(deploystart = min(date), deployend = max(date))

# Generate deployment summary
ds <- loc %>% 
  filter(trip == 1) %>% 
  group_by(id, seal_id, trip) %>% 
  summarise(trip1_start = min(date), 
            trip1_end = max(date), 
            trip1_complete = first(is_trip_complete),
            trip1_dur = difftime(trip1_end, trip1_start, units = 'days'))

ds$trip1_complete[is.na(ds$trip1_complete)] <- FALSE

# Add deployment start and end dates
ds <- left_join(ds, deploy)

# Add life history and sightings info
ds <- left_join(ds, survival)

# Add if seal seen alive after 6 months (survived 1st trip)
ds <- ds %>% mutate(trip1_alive = seen_6m)

# Which seal wasn't seen after 6 months but had complete trip1 track?
ds$trip1_alive[which(ds$seen_6m == FALSE & ds$trip1_complete == TRUE)] <- TRUE

# Calculate time spent on land before first foraging trip
ds <- ds %>%
  mutate(birth2trip1_dur = difftime(trip1_start, birthdate, units = 'days') %>% round(), 
         wean2trip1_dur = difftime(trip1_start, weandate, units = 'days') %>% round())

# Join deployment summary to seal locations as a column list
ds <- left_join(ds, loc %>% 
  group_by(id, seal_id) %>% 
  nest() %>% 
  rename(track = data))

# Join deployment summary to drift info as a column list.
load('./Output/drift_rate.RData')

ds <- ds %>% 
  left_join(drift %>% 
              group_by(id) %>% 
              summarise(ddrate_max = max(drate_delta, na.rm = T), 
                        ddrate_min = min(drate_delta, na.rm = T),
                        ddrate_sd = sd(drate_delta, na.rm = T), 
                        ddrate_mean = mean(drate_delta, na.rm = T))) %>% 
  left_join(drift %>% 
              group_by(id) %>% 
              nest() %>% 
              rename(drift = data))

# Add post trip 1 weights (closest observation to trip end date)
ds <- ds %>% left_join(
  left_join(sightings, ds %>% 
              dplyr::select(id, trip1_end) %>%
              mutate(trip1_enddate = trip1_end %>% date)) %>% 
    group_by(id) %>% 
    filter(date_observed > trip1_enddate, !is.na(weight)) %>% 
    mutate(date_diff = difftime(date_observed, trip1_enddate)) %>% 
    filter(date_diff == min(date_diff)) %>% 
    rename(postTrip1_date_observed = date_observed,
           postTrip1_weight = weight,
           postTrip1_slength = standard_length,
           postTrip1_girth = girth,
           postTrip1_dateDiff = date_diff) %>% 
    dplyr::select(-trip1_enddate)
)

write.csv(ds %>% dplyr::select(-drift, -track), file = './Output/overall_summary_table(1).csv', na = '-', row.names = FALSE)

# save.image("11_survival_WORKSPACE.RData")
# 


# Summary Tables ----------------------------------------------------------

# seen again at least after 1 year 
table(ds$seen_1y) 
# seen again at least after 6 months
table(ds$seen_6m) 
# had complete foraging tracks
table(ds$trip1_complete) 
# either completed trip 1 or was seen again after 6m
table(ds$trip1_alive) 

