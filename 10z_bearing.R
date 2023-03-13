
rm(list = ls())
options(tibble.width = Inf)
# devtools::install_github("tidyverse/reprex")
library(tidyverse)
library(sf)
library(lubridate)
library(geosphere)
library(circular)
library(ggpubr)
library(broom)
library(rstatix)
library(effectsize)
library(reprex)

# Load data ---------------------------------------------------------------

locw <- readRDS("./Output/tracks_processed_12h.rds")
masterData <- readRDS("./Output/all_data_combined.rds")
load('baseInfo.rdata')
source('convert2polarsf.R')
source('functions.R')

load("10z_bearing_WORKSPACE.Rdata")


# Identify and remove seals with no weanmass data -------------------------

masterData %>% 
  group_by(id) %>% 
  summarise(weanmass = mean(weanmass, na.rm = T),
            weanmass_first = first(weanmass)) %>% 
  filter(is.na(weanmass)) %>% pull(id) -> sealsWithNoWeanmass


# Set up table outputs ----------------------------------------------------
departureDates_results <- tibble()

# 1) Weaners ----------------------------------------------------------
# the location animal is at on this date to calculate bearing from the colony
REF_DATE <- 5

# Get only first trips
locw1 <- locw %>% 
  group_by(id) %>% 
  filter(trip == 1, SUS == FALSE) %>% 
  filter(!id %in% sealsWithNoWeanmass) %>% 
  # create variable days since start of first trip 
  mutate(daysSinceStart = difftime(date, min(date), units = 'day') %>% round())

# # visualise tracks
# locw1 %>% 
#   ggplot() +
#   geom_sf(size = 0.2, aes(colour = daysSinceStart <= REF_DATE)) 
# 
# # Get bearing for the early part (6 days) of the first trip
loc_early <- locw1 %>%
  filter(daysSinceStart == REF_DATE) %>%
  group_by(id) %>%
  summarise(startdate = min(date),
            enddate = max(date),
            bearing = bearing(mq, c(last(lon), last(lat))),
            type = 'argos',
            dist2col = last(dist2col))

# loc_early$enddate %>% month %>% unique
# loc_early$enddate %>% year %>% unique

# Get mean departure date
loc_early %>% 
  group_by(id) %>% 
  summarise(startdate = min(startdate)) %>% 
  pull(startdate) -> jdates

# Add to output table
departureDates_results <- saveDepartureDateResult(departureDates_results, jdates, "weaners")


# 2) Adult females ---------------------------------------------------
# 10.1111/gcb.13776 Full winter foraging trips were obtained from a total of 67
# seals in six years: 2000 (n = 7), 2001 (n = 8), 2002 (n = 10), 2004 (n = 15),
# 2005 (n = 12) and 2010 (n = 15). The mean winter (post-moult) trip duration
# was 211.0 +- 51.5 days, with a median departure date of 4th of February and a
# median return date of 23rd of September, arriving for the subsequent
# breeding season. The seals had a mean departure mass in February of 370 +- 66
# kg, and a mean arrival mass of 514 +- 43 kg after their time at sea.
load("./Data/macca_winter_locs.Rdata") # adult females
locf <- loc %>% as_tibble() %>% 
  arrange(seal, gmt) # time range "2000-02-03 00:00:00 GMT" "2010-10-24 19:31:03 GMT"

locf_sf <- convert2polarsf(locf) 
locf_sf <- locf_sf %>% left_join(locf) %>% 
  rename(id = seal)

# inspect and visualise data
# relationship between month and season? 
# locf %>% 
#   ggplot(aes(x = month, y = season)) + 
#   geom_point()

# create days since start of trip, and dist2col variables
locf_sf <- locf_sf %>% 
  # filter(season == 'pm1') %>% 
  group_by(id) %>% 
  mutate(daysSinceStart = difftime(gmt, min(gmt), units = 'day') %>% round(),
         dist2col = SGAT::gcDist(mq, cbind(lon, lat)))



# some seal tracks don't seem to start from the colony... remove those
locf_sf %>%
  ggplot() + 
  geom_sf(aes(colour = daysSinceStart <= REF_DATE), size = 0.01)

# # which seals first dist2col > 200km from colony?
# locf_sf %>% 
#   group_by(id) %>% 
#   summarise(dist2col = first(dist2col)) %>% 
#   ggplot(aes(x = dist2col)) + 
#   geom_histogram() 

locf_sf %>% 
  group_by(id) %>% 
  summarise(dist2col = first(dist2col)) %>% 
  filter(dist2col > 200) %>% pull(id) ->  bad_females

# Visualise again without bad females
# locf_sf %>%
#   filter(!id %in% bad_females) %>%
#   ggplot() + 
#   geom_sf(aes(colour = daysSinceStart <= REF_DATE), size = 0.01)


locf_sf <- locf_sf %>% filter(!id %in% bad_females)

# Create new dataframe for just early locations and calculate bearing
locf_early <-locf_sf %>% 
  mutate(date = as_date(gmt)) %>% 
  group_by(id, date) %>%
  arrange(id, gmt) %>% 
  filter(gmt == last(gmt)) %>% 
  filter(daysSinceStart == REF_DATE) %>% 
  group_by(id) %>% 
  summarise(startdate = min(gmt),
            enddate = max(gmt), 
            bearing = bearing(mq, c(last(lon), last(lat))),
            type = last(type), 
            dist2col = last(dist2col))

# save mean departure date
locf_early %>% 
  group_by(id) %>% 
  summarise(startdate = min(startdate)) %>% 
  pull(startdate) -> jdates

departureDates_results <- saveDepartureDateResult(departureDates_results, jdates, "adult females")

# write.csv(departureDates_results, "./output/dispersal bearing/median departure dates.csv", row.names = FALSE)

# 3) Particle trace ----------------------------------------------------------
pt <- readRDS("./Output/currently_particleTrace.rds")

pt <- pt %>% group_by(id) %>% 
  filter(!is.na(x)) %>% 
  mutate(daysSinceDeparture = difftime(date, min(date), units = 'days'))

# plot trace
pt %>% 
  ggplot() + 
  geom_point(size = 0.1, aes(lon, lat, col = daysSinceDeparture %>% as.numeric, fill = NULL, group = group)) + 
  annotate(geom = 'point', x =158.95, y = -54.5, shape = 17, col = 'red', size = 5) +
  labs(caption = 'Particle trace of 44 seals from MI', col = 'Days since departure') +
  viridis::scale_color_viridis() + 
  lims(x = c(155, 180))  +
  theme_bw()

# convert trace to sf
pt_sf <- convert2polarsf(pt) %>% 
  left_join(pt) %>% 
  mutate(dist2col = SGAT::gcDist(mq, cbind(lon, lat)))

pt_sf %>% 
  group_by(group) %>% 
  summarise(min(daysSinceDeparture))


# ~ calc bearing ------------------------------------

## METHOD 1: determine where to calculate bearing from for particle trace based on
## corresponding seal dist2col at REF_DATE
pt_sf <- pt_sf %>% 
  filter(id %in% loc_early$id) %>% 
  left_join(as.data.frame(loc_early) %>% 
              dplyr::select(id, dist2col) %>% 
              rename(dist2col.w = dist2col))

pt_sf %>%
  # filter(dist2col <= dist2col.w) %>%
  ggplot() +
  geom_sf(size = 0.1, aes(color = dist2col <= dist2col.w))

pt_early <- pt_sf %>% 
  group_by(id) %>% 
  mutate(dist2col_diff = abs(dist2col.w - dist2col), 
         min_dist2col_diff = min(dist2col_diff)) %>% 
  filter(dist2col_diff == min_dist2col_diff) %>% 
  mutate(bearing = bearing(mq, c(lon, lat)),
         type = NA) %>% 
  dplyr::select(id, bearing, type, dist2col, dist2col.w)

# OR METHOD 2: determine where to calculate bearing from for particle trace based on the mean
# distance of all weaners from colony at REF_DATE into their trip
# pt_sf %>%
#   filter(dist2col <= mean(loc_early$dist2col)) %>%
#   ggplot() +
#   geom_sf(size = 0.1)
# 
# pt_early <- pt_sf %>%
#   group_by(id) %>%
#   filter(dist2col <= mean(loc_early$dist2col)) %>%
#   summarise(bearing = bearing(mq, c(last(lon), last(lat))),
#             type = NA,
#             dist2col = last(dist2col))


## OR METHOD 3: determine reference point based solely on days since start of track
# pt_sf %>%
#   filter(daysSinceDeparture <= REF_DATE) %>%
#   ggplot() +
#   geom_sf(size = 0.1)
# 
# pt_early <- pt_sf %>%
#   group_by(id) %>%
#   filter(daysSinceDeparture == REF_DATE) %>%
#   summarise(bearing = bearing(mq, c(last(lon), last(lat))),
#             type = NA,
#             dist2col = last(dist2col))


# ~ plot weaner tracks / particle trace -----------------------------------
bb <- extent(locw1)
orsi_sf$nudge_x <- c(-1000, -4800, -4500, -4000)
orsi_sf$nudge_y <- c(0, 5300, 4800, 3000)

# all in one map
png('./output/weaner_particleTrace_adultFemales_1map.png', width=20, height=13, units='cm', res=500)
bb <- extent(locf_sf)
ggplot() + 
  geom_sf(data = locf_sf, size = .1, col = 'grey') +
  geom_sf(data = locw1 %>% filter(id %in% loc_early$id), size = 0.1, col = 'firebrick1') + 
  geom_sf(data = pt_sf, col = 'dodgerblue', size = 0.1, alpha = .3) + 
  labs(caption = 'blue = particle trace, red = weaners, grey = adult females', x = 'lon', y = 'lat') +
  geom_sf(data = orsi_sf, linetype = 'dashed', colour = 'grey10') + 
  geom_sf_text(data = orsi_sf, aes(label = front), nudge_x = orsi_sf$nudge_x, nudge_y = orsi_sf$nudge_y) + 
  geom_sf(data = mq_sf, shape = 17, fill = 'black' ) +
  xlim(bb[1], bb[2]) + ylim(bb[3]-200, bb[4]+500) +
  theme_bw() 
dev.off()


# facet by individual
png('./output/weaner_particleTrace_by_seal.png', width=20, height=18, units='cm', res=500)
ggplot() + 
  geom_sf(data = locw1 %>% filter(id %in% loc_early$id), size = 0.1, col = 'firebrick1') + 
  geom_sf(data = pt_sf, col = 'dodgerblue', size = 0.1, alpha = .3) + 
  labs(caption = 'blue = particle trace, red = weaner', x = 'lon', y = 'lat') +
  # add reference points for bearing calculation
  geom_sf(data = loc_early, shape = 10) + 
  geom_sf(data = pt_early, shape = 10) +
  theme_bw() + 
  facet_wrap(~id)
dev.off()

# Adult females only with reference locations
png('./output/adultFemale_tracks_by_seal.png', width=20, height=25, units='cm', res=500)
ggplot() + 
  geom_sf(data = locf_sf, size = 0.1, alpha = .3) + 
  labs(x = 'lon', y = 'lat') +
  # add reference points for bearing calculation
  geom_sf(data = locf_early, shape = 10, size = .8, col = 'red') +
  theme_bw() + 
  facet_wrap(~id)
dev.off()





# 5) Circular stats ----------------------------------------------------------
circ.f <- circular(locf_early$bearing%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
plot(circ.f, stack = TRUE, shrink = 1, main = 'females')
arrows.circular(mean(circ.f))

circ.w <- circular(loc_early$bearing%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
plot(circ.w, stack = TRUE, shrink = 1.6, main = 'weaners')
arrows.circular(mean(circ.w))

circ.p <- circular(pt_early$bearing%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
plot(circ.p, stack = TRUE, shrink = 1.6, main = 'particle trace')
arrows.circular(mean(circ.p))


# ~~ !generate summaries! --------------------------------------------------

prepIndividualDf <- function(x, group_name, bearing, y) {
  newDf <- tibble(id = x$id, 
                  group = group_name, 
                  type = x$type, 
                  bearing = bearing, 
                  refLocation_dist2col = x$dist2col)
  
  if(exists("y")) {
    y1 <- y %>% 
      group_by(id) %>% 
      summarise(startdate = min(date) %>% as_date,
                enddate = max(date) %>% as_date, 
                trip_dur = difftime(enddate, startdate))
  }
  
  out <- left_join(newDf, y1) %>% select(-geometry)
  return(out)
}

indiv_summ <- prepIndividualDf(loc_early, "weaner", circ.w, locw) %>% 
  bind_rows(prepIndividualDf(locf_early, "adult female", circ.f, locf_sf %>% rename(date = gmt))) %>% 
  bind_rows(prepIndividualDf(pt_early, "particle trace", circ.p, pt_sf))

indiv_summ <- indiv_summ %>% mutate(compass_zone = sapply(bearing, whichZone))


group_summ <- indiv_summ %>%
  group_by(group) %>% 
  summarise(n = n(), 
            bearing_mean = mean(bearing),
            bearing_se = sd(bearing)/sqrt(n),
            bearing_min = min(bearing), 
            bearing_max = max(bearing),
            refLocation_dist2col_mean = mean(refLocation_dist2col),
            refLocation_dist2col_se = sd(refLocation_dist2col)/sqrt(n),
            refLocation_dist2col_min = min(refLocation_dist2col),
            refLocation_dist2col_max = max(refLocation_dist2col),
            tripDur_mean = mean(trip_dur),
            tripDur_se = sd(trip_dur)/sqrt(n),
            tripDur_min = min(trip_dur),
            tripDur_max = max(trip_dur)
            )

write_csv(group_summ, "./output/dispersal bearing/group_summary.csv")
write_csv(indiv_summ, "./output/dispersal bearing/allBearings.csv")


# 6) Plot circular histogram -------------------------------------------------
# adult female
p1 <- ggplot() +
  geom_histogram(data = data.frame(circ.f), aes(x = circ.f), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.f), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.f), y = 20, label = mean(circ.f) %>% round(1)) +
  labs(subtitle = paste('adult females (n=', length(circ.f), ')', sep = ''), y = "") + 
  theme_bw() 


# weaner
p2 <- ggplot() +
  geom_histogram(data = data.frame(circ.w), aes(x = circ.w), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.w), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.w), y = 20, label = mean(circ.w) %>% round(1)) + 
  labs(subtitle = paste('weaners (n=', length(circ.w), ')',sep = ''), 
       y = "count")  + 
  theme_bw()

# particle trace
p3 <- ggplot() +
  geom_histogram(data = data.frame(circ.p), aes(x = circ.p), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.p), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.p), y = 20, label = mean(circ.p) %>% round(1)) + 
  labs(subtitle = paste('particle trace (n=', length(circ.p), ')',sep = ''),
       y = "") + 
  theme_bw()

# ~ save plots ------------------------------------------------------------
png('dispersal_direction.png', width=20, height=10, units='cm', res=500)
ggarrange(p2, p3, p1, ncol = 3, labels = c('a', 'b', 'c'))
dev.off()

# png('dispersal_direction_weaner_juvenile_adult.png', width=50, height=12, units='cm', res=500)
# ggarrange(p2, p4, p5, p6, p7, p1, ncol = 6)
# dev.off()


# ~ tests -------------------------------------------------------------------

### test if data follows a von Mises distribution (p < 0.05)
watson.test(circ.f)
watson.test(circ.w)


### compare means of females and weaners
# https://www.biorxiv.org/content/10.1101/2021.03.25.436932v1.full
mean(circ.f) %>% as.numeric(); mean(circ.w) %>% as.numeric()
watson.two.test(circ.f, circ.w) # s. p < 0.001
watson.two.test(circ.p, circ.w) # n.s. p > 0.10
watson.two.test(circ.p, circ.f) # s. p < 0.001

### Test for direction uniformity 
# Hermans-Rasson test for uniformity (alternative to Rayleigh test). Better for multimodel distributions
# https://bmcecol.biomedcentral.com/articles/10.1186/s12898-019-0246-8
# p < 0.05 = directions are not uniform ie animals have a preferred direction.
library(CircMLE)
HR_test(circ.f) 
HR_test(circ.w)

# Our data is unimodal - so rayleigh test is fine. Either way both tests p < 0.05
rayleigh.test(circ.f)
rayleigh.test(circ.w)


# # ~ manova approach (14 Jun 2022) -----------------------------------------
# # Try this approach for comparing bearings of multiple groups (> 2)
# # Code from: https://www.nature.com/articles/s41598-021-99299-5
# g1 <- c("w", "f")
# g2 <- list(circ.w, circ.f)
# 
# all_circ_df <- purrr::map2(g1, g2, function(x, y) {
#   return(bind_cols(bearing = y, group = x))
# }) %>% reduce(bind_rows)
# # all_circ_df <- all_circ_df %>% filter( group %in% c("w", "f"))
# 
# # tutorial: https://www.reneshbedre.com/blog/manova.html
# dep_vars <- cbind(cos(all_circ_df$bearing), sin(all_circ_df$bearing))
# age_groups <- all_circ_df$group
# 
# man_fit <- manova(dep_vars ~ age_groups)
# tidy(man_fit)
# # write.csv(tidy(man_fit), file = "./output/dispersal bearing/manova_fit.csv", row.names = FALSE)
# 
# # effect size
# effectsize::eta_squared(man_fit)
# 
# # check anova homogeneity assumption
# all_circ_df2 <- tibble(
#   cos_bearing = cos(all_circ_df$bearing) %>% as.numeric(), 
#   sin_bearing = sin(all_circ_df$bearing) %>% as.numeric(), 
#   age = age_groups)
# 
# grouped_data <- all_circ_df2 %>% 
#   gather(key = "variable", value = "value", cos_bearing, sin_bearing) %>%
#   group_by(variable) 
# 
# grouped_data %>% levene_test(value ~ age) # homogeneity assumption met for both response variables
# 
# # one-way anova test
# grouped_data %>% anova_test(value ~ age) # only cos(bearing) is different between groups ie EAST <> WEST direction (x-component)
# 
# # post-hoc tests
# tukey_hsd(all_circ_df2, cos_bearing ~ age)


# Create output dataframe -------------------------------------------------

out <- tibble(loc_early) %>% 
  dplyr::select(id, bearing) %>% 
  rename(bearing.w = bearing) %>% 
  left_join(tibble(pt_early) %>% dplyr::select(id, bearing)) %>% 
  rename(bearing.pt = bearing) %>% 
  group_by(id) %>% 
  mutate(bearing_diff = angle_diff(bearing.w, bearing.pt), 
         ew_zone = eastOrWest(bearing.w), 
         compass_zone = whichZone(bearing.w))

# saveRDS(out, file = './output/dispersal_bearing___REF_DATE=6d.rds')

# Survive vs dead weaners -------------------------------------------------

# ~ bearing ---------------------------------------------------------------

d <- readRDS('./Output/all_data_combined.rds')
ds <- d %>% 
  group_by(id) %>% 
  filter(trip == 1, sim == 0) %>% 
  summarise(birthyear = first(birthyear),
            surviveTrip1 = ifelse(first(is_trip_complete) == T | first(seen_6m) == T, T, F),
            surviveYear1 = ifelse(first(seen_1y) == TRUE, T, F),
            weanmass = first(weanmass))

b <- out %>% 
  ungroup() %>% 
  left_join(ds)

# Add size category variable
b$sizeCat <- 'avg'
b$sizeCat[b$weanmass > 135] <- 'heavy'
b$sizeCat[b$weanmass < 96] <- 'light'

table(b$surviveTrip1)
table(b$surviveYear1)


# ~~ plot 1st trip survival -----------------------------------------------


# set which response variable to use
b <- b %>% mutate(survive = surviveTrip1)

circ.live <- circular(b %>% filter(survive == TRUE) %>% pull(bearing.w)%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
circ.dead <- circular(b %>% filter(survive == FALSE) %>% pull(bearing.w)%%360, units = 'degrees', zero = pi/2, rotation = 'clock')

p1 <- ggplot() +
  geom_histogram(data = data.frame(circ.live), aes(x = circ.live), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  # add mean direction annotation
  geom_vline(xintercept = mean(circ.live), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.live), y = 20, label = mean(circ.live) %>% round(1)) +
  # labs
  labs(subtitle = paste0('1st trip survived, n = ', circ.live %>% length), y = "count") + 
  theme_bw() 

p2 <- ggplot() +
  geom_histogram(data = data.frame(circ.dead), aes(x = circ.dead), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  # add mean direction annotation
  geom_vline(xintercept = mean(circ.dead), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.dead), y = 4, label = mean(circ.dead) %>% round(1)) +
  # labs
  labs(subtitle = paste0('1st trip died, n = ', circ.dead %>% length), y = "") + 
  theme_bw() 

# png('./output/dispersal_direction_by_trip1survival.png', width=20, height=10, units='cm', res=500)
ggarrange(p1, p2, ncol = 2, labels = c('a', 'b'), align = 'hv')
# dev.off()



# ~~ plot 1st year survival -----------------------------------------------

b <- b %>% mutate(survive = surviveYear1)

circ.live <- circular(b %>% filter(survive == TRUE) %>% pull(bearing.w)%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
circ.dead <- circular(b %>% filter(survive == FALSE) %>% pull(bearing.w)%%360, units = 'degrees', zero = pi/2, rotation = 'clock')

p3 <- ggplot() +
  geom_histogram(data = data.frame(circ.live), aes(x = circ.live), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  # add mean direction annotation
  geom_vline(xintercept = mean(circ.live), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.live), y = 20, label = mean(circ.live) %>% round(1)) +
  # labs
  labs(subtitle = paste0('1st year survived, n = ', circ.live %>% length), y = "count") + 
  theme_bw() 

p4 <- ggplot() +
  geom_histogram(data = data.frame(circ.dead), aes(x = circ.dead), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  # add mean direction annotation
  geom_vline(xintercept = mean(circ.dead), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.dead), y = 4, label = mean(circ.dead) %>% round(1)) +
  # labs
  labs(subtitle = paste0('1st year died, n = ', circ.dead %>% length), y = "") + 
  theme_bw() 

# png('./output/dispersal_direction_by_year1survival.png', width=20, height=10, units='cm', res=500)
ggarrange(p3, p4, ncol = 2, labels = c('a', 'b'), align = 'hv')
# dev.off()



# ~~ combined survival plot -----------------------------------------------

# png('./output/dispersal_direction_by_year1 & trip1 survival.png', width=20, height=20, units='cm', res=500)
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = c('a', 'b', 'c', 'd'), align = 'hv')
# dev.off()


# 
# # ~ tests -----------------------------------------------------------------
# 
# ## Test if directing means of survive vs dead are different
# watson.two.test(circ.live, circ.dead)
# 
# ## Mean test of particle trace vs survive and dead seals
# watson.two.test(circ.p, circ.dead)
# watson.two.test(circ.p, circ.live)
# 
# ## test against adult females
# watson.two.test(circ.f, circ.dead)
# watson.two.test(circ.f, circ.live)
# 
# 
# ## Test each group for uniformity of directions
# CircMLE::HR_test(circ.live) # p = 0.0001
# CircMLE::HR_test(circ.dead) # p = 0.0015
# 
# 
# ## do weaners that go with the flow but DON'T survive lighter in weanmass?
# b_flow <- b %>% 
#   filter(compass_zone == 'E-SE')
# 
# table(b_flow$survive) 
# t.test(weanmass ~ survive, data = b_flow) # p>0.05 not significantly different
# 
# ## survival rate by weight group
# table(b$sizeCat, b$survive) 
# 6/10 #avg
# 15/20 #light
# 11/(7+11) #heavy
# 
# 
# ## closer look at flow seals
# d_flow <- d %>% 
#   filter(id %in% b_flow$id, sim == 0) %>% 
#   filter(trip == 1, land == FALSE, SUS == FALSE) %>% 
#   left_join(b_flow %>% dplyr::select(id, survive))
# 
# d_flow <- convert2polarsf(d_flow)
# 
# bb <- extent(d_flow)
# survive.labs <- c('died n=6', "survived n=23")
# names(survive.labs) <- c(FALSE, TRUE)
# 
# d_flow %>% ggplot() +
#   geom_sf(size = 0.5, aes(col = g)) + 
#   # orsi fronts
#   geom_sf(data = orsi_sf, linetype = 'dashed', colour = 'grey10') + 
#   geom_sf_text(data = orsi_sf, aes(label = front), nudge_x = orsi_sf$nudge_x, nudge_y = orsi_sf$nudge_y) +
#   xlim(bb[1], bb[2]) + ylim(bb[3]-200, bb[4]+500) +
#   facet_wrap(~survive, labeller = labeller(survive = survive.labs)) + 
#   viridis::scale_color_viridis() + 
#   facet_wrap(~id) +
#   theme_bw()
# 
# d_flow %>% ggplot() +
#   geom_sf(size = 0.5, aes(col = g < 0.5)) + 
#   # orsi fronts
#   geom_sf(data = orsi_sf, linetype = 'dashed', colour = 'grey10') + 
#   geom_sf_text(data = orsi_sf, aes(label = front), nudge_x = orsi_sf$nudge_x, nudge_y = orsi_sf$nudge_y) +
#   xlim(bb[1], bb[2]) + ylim(bb[3]-200, bb[4]+500) +
#   facet_wrap(survive~id) +
#   theme_bw()
# 
# t.test(g ~ survive, data = d_flow)
# 
# # ~ timing of departure effect ---------------------------------------------
# b <- b %>% left_join(loc_early %>% dplyr::select(id, startdate, enddate))
# b$geometry <- NULL
# 
# ## create dummy date with same year for plotting
# bb <- b %>% 
#   mutate(startdate2 = startdate %>% format('2000-%m-%d') %>% as.Date())
# bb$startdate2[bb$id == 'mq3-17217-99'] <- bb$startdate[bb$id == 'mq3-17217-99'] %>% format('2001-%m-%d') %>% as.Date()
# bb$startdate2[bb$id == 'mq4-Flora-00'] <- bb$startdate[bb$id == 'mq4-Flora-00'] %>% format('2001-%m-%d') %>% as.Date()
# bb$startdate2[bb$id == 'mq4-Ella-00'] <- bb$startdate[bb$id == 'mq4-Ella-00'] %>% format('2001-%m-%d') %>% as.Date()
# 
# 
# bb %>% 
#   ggplot(aes(x = startdate2)) + 
#   geom_histogram(aes(fill = survive), binwidth = 1) + 
#   facet_wrap(~birthyear) +
#   labs(x = 'departure date') +
#   theme_bw()
# 
# 
# bb %>% 
#   ggplot(aes(x = startdate2)) + 
#   geom_histogram(aes(fill = sizeCat), binwidth = 1) + 
#   facet_wrap(~birthyear) +
#   labs(x = 'departure date') +
#   theme_bw()
# 
# # Weanmass vs dispersal bearing -------------------------------------------
# 
# 
# circ.avgSize <- circular(b %>% filter(sizeCat == 'avg') %>% pull(bearing.w)%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
# circ.light <- circular(b %>% filter(sizeCat == 'light') %>% pull(bearing.w)%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
# circ.heavy <- circular(b %>% filter(sizeCat == 'heavy') %>% pull(bearing.w)%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
# 
# 
# p1 <- b %>% filter(sizeCat == 'heavy')  %>% 
#   ggplot() +
#   geom_histogram(aes(x = bearing.w%%360, fill = survive), 
#                  breaks = seq(0, 360, 45), 
#                  colour = "black") + 
#   coord_polar() +
#   scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
#   # add mean direction annotation
#   geom_vline(xintercept = mean(circ.heavy), color = "black", linetype = 2, size = 1) +
#   annotate("label", x = mean(circ.heavy), y = 6, label = mean(circ.heavy) %>% round(1)) +
#   # labs
#   labs(title = 'weaners heavy (>135kg)', y = "number of individuals") + 
#   theme_bw()+
#   theme(legend.position = 'none')
# 
# p2 <- b %>% filter(sizeCat == 'avg')  %>% 
#   ggplot() +
#   geom_histogram(aes(x = bearing.w%%360, fill = survive), 
#                  breaks = seq(0, 360, 45), 
#                  colour = "black") + 
#   coord_polar() +
#   scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
#   # add mean direction annotation
#   geom_vline(xintercept = mean(circ.avgSize), color = "black", linetype = 2, size = 1) +
#   annotate("label", x = mean(circ.avgSize), y = 4, label = mean(circ.avgSize) %>% round(1)) +
#   # labs
#   labs(title = 'weaners average', y = "") + 
#   theme_bw() +
#   theme(legend.position = 'bottom')
# 
# 
# 
# p3 <- b %>% filter(sizeCat == 'light')  %>% 
#   ggplot() +
#   geom_histogram(aes(x = bearing.w%%360, fill = survive), 
#                  breaks = seq(0, 360, 45), 
#                  colour = "black") + 
#   coord_polar() +
#   scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
#   # add mean direction annotation
#   geom_vline(xintercept = mean(circ.light), color = "black", linetype = 2, size = 1) +
#   annotate("label", x = mean(circ.light), y = 4, label = mean(circ.light) %>% round(1)) +
#   # labs
#   labs(title = 'weaners light (<96kg)', y = "") + 
#   theme_bw() +
#   theme(legend.position = 'none')
# 
# 
# # png('dispersal_direction_by_survival.png', width=20, height=10, units='cm', res=500)
# ggarrange(p1, p2,p3,  ncol = 3, labels = c('a', 'b', 'c'), align = 'hv')
# # dev.off()

# Save Image --------------------------------------------------------------
save.image(file = "10z_bearing_WORKSPACE.Rdata")

