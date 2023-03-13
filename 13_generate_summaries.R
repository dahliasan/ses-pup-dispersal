# Generate data summaries

rm(list = ls())
options(tibble.width = Inf)
library(tidyverse)


# Load data ---------------------------------------------------------------

locw <- readRDS("./Output/tracks_processed_12h.rds")
a <- readRDS("./Output/all_data_combined.rds")
load('baseInfo.rdata')
source('convert2polarsf.R')
source('functions.R')

a1 <- a %>% 
  filter(sim == 0, trip == 1, land == FALSE, !is.na(bearing.pt), !is.na(weanmass))

a2 <- a1 %>% 
  ungroup() %>% 
  select(1:12, 33:48, 51:56 )

indiv_summ <- a2 %>% 
  group_by(id) %>% 
  summarise(tripdur = first(tripdur), 
            dist2col_max = max(dist2col),
            survived_trip1 = ifelse(first(is_trip_complete) == T | first(seen_6m) == T, T, F),
            survived_year1 = ifelse(first(seen_1y) == TRUE, T, F),
            last_seen_year = first(last_seen_year),
            birthlocation = first(birthlocation),
            birthdate = first(birthdate),
            weandate = first(weandate),
            blackmass = first(blackmass),
            weanmass = first(weanmass),
            birthyear = first(birthyear),
            bearing = first(bearing.w),
            bearing_particleTrace = first(bearing.pt),
            compass_direction = first(compass_zone)
            )

group_summ <- indiv_summ %>% 
  ungroup() %>% 
  summarise(n = n(),
            tripdur_mean = mean(tripdur),
            tripdur_se = sd(tripdur)/sqrt(n), 
            tripdur_min = min(tripdur), 
            tripdur_max = max(tripdur), 
            dist2col_mean = mean(dist2col_max),
            dist2col_se = sd(dist2col_max)/sqrt(n),
            dist2col_min = min(dist2col_max),
            dist2col_max = max(dist2col_max),
            n_survived_trip1 = sum(survived_trip1),
            perc_survived_trip1 = n_survived_trip1/n,
            n_survived_year1 = sum(survived_year1),
            perc_survived_year1 = n_survived_year1/n,
            birthdate_mean = yday(birthdate) %>% mean %>% as_date(origin = "2000-01-01"),
            birthdate_se = sd(yday(birthdate))/sqrt(n),
            weandate_mean = yday(weandate) %>% mean %>% as_date(origin = "2000-01-01"),
            weandate_se = sd(yday(weandate))/sqrt(n),
            )


  

write_csv(indiv_summ, "./output/individual_summary.csv")
write_csv(group_summ, "./output/group_summary.csv")


# Adult female summaries --------------------------------------------------
REF_DATE = 5 
load("./Data/macca_winter_locs.Rdata")  # adult females

# plot deployment timeline of each female
g <- loc %>% 
  tibble() %>% 
  group_by(seal) %>% 
  summarise(startdate = min(gmt), 
            enddate = max(gmt), season = first(season))

g$enddate %>% month() %>% hist()
g$startdate %>% month() %>% hist()

loc %>% 
  ggplot(aes(x = gmt, y = seal)) +
  geom_point(aes(color = season)) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%b-%y")

locf <- loc %>% as_tibble() %>% 
  arrange(seal, gmt) # time range "2000-02-03 00:00:00 GMT" "2010-10-24 19:31:03 GMT"

locf_sf <- convert2polarsf(locf) 
locf_sf <- locf_sf %>% left_join(locf) %>% 
  rename(id = seal)

locf_sf %>% 
  ggplot() + 
  geom_sf()

# create days since start of trip, and dist2col variables
locf_sf <- locf_sf %>% 
  # filter(season == 'pm1') %>% 
  group_by(id) %>% 
  mutate(daysSinceStart = difftime(gmt, min(gmt), units = 'day') %>% round(),
         dist2col = SGAT::gcDist(mq, cbind(lon, lat)))

locf_sf %>% 
  group_by(id) %>% 
  summarise(dist2col = first(dist2col)) %>% 
  filter(dist2col > 200) %>% pull(id) ->  bad_females



f <- locf_sf %>% 
  filter(daysSinceStart <= REF_DATE) %>%
  filter(!id %in% bad_females) %>% 
  group_by(id) %>% 
  summarise(startdate = min(gmt),
            enddate = max(gmt), 
            bearing = bearing(mq, c(last(lon), last(lat))),
            type = last(type), 
            dist2col = last(dist2col))
