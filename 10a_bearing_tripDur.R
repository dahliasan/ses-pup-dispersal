### Testing different trip duration lengths to calculate bearing


rm(list = ls())
options(tibble.width = Inf)
library(tidyverse)
library(sf)
library(lubridate)
library(geosphere)
library(circular)
library(ggpubr)

# Load data ---------------------------------------------------------------

locw <- readRDS("./Output/tracks_processed_12h.rds")
load('baseInfo.rdata')
source('convert2polarsf.R')

lastDayThreshold <- 60

# 1) Analyse weaners ----------------------------------------------------------

# Get only first trips
locw1 <- locw %>% 
  group_by(id) %>% 
  filter(trip == 1, SUS == FALSE) %>% 
  # create variable days since start of first trip 
  mutate(daysSinceStart = difftime(date, min(date), units = 'day') %>% round())

# # visualise tracks
# locw1 %>% 
#   ggplot() +
#   geom_sf(size = 0.2, aes(colour = daysSinceStart <= lastDayThreshold)) + 
#   facet_wrap(~id) 

# Get bearing for the early part  of the first trip
loc_early <- locw1 %>% 
  filter(daysSinceStart <= lastDayThreshold) %>% 
  group_by(id) %>% 
  summarise(bearing = bearing(mq, c(last(lon), last(lat))),
            type = 'argos', 
            dist2col = last(dist2col))


# 2) Adult female analysis ---------------------------------------------------
# 10.1111/gcb.13776 Full winter foraging trips were obtained from a total of 67
# seals in six years: 2000 (n = 7), 2001 (n = 8), 2002 (n = 10), 2004 (n = 15),
# 2005 (n = 12) and 2010 (n = 15). The mean winter (post-moult) trip duration
# was 211.0 +- 51.5 days, with a median departure date of 4th of February and a
# median return date of 23rd of September, arriving for the subsequent
# breeding season. The seals had a mean departure mass in February of 370 +- 66
# kg, and a mean arrival mass of 514 +- 43 kg after their time at sea.
load("./Data/macca_winter_locs.Rdata") # adult females
locf <- loc %>% as_tibble() %>% 
  arrange(seal, gmt)

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
  filter(season == 'pm1') %>% 
  group_by(id) %>% 
  mutate(daysSinceStart = difftime(gmt, min(gmt), units = 'day') %>% round(),
         dist2col = SGAT::gcDist(mq, cbind(lon, lat)))

# # some seal tracks don't seem to start from the colony... remove those
# locf_sf %>%
#   ggplot() + 
#   geom_sf(aes(colour = daysSinceStart <= lastDayThreshold), size = 0.01)

# which seals first dist2col > 200km from colony?
# locf_sf %>%
#   group_by(id) %>%
#   summarise(dist2col = first(dist2col)) %>%
#   ggplot(aes(x = dist2col)) +
#   geom_histogram()

locf_sf %>% 
  group_by(id) %>% 
  summarise(dist2col = first(dist2col)) %>% 
  filter(dist2col > 200) %>% pull(id) ->  bad_females

locf_sf <- locf_sf %>% filter(!id %in% bad_females)

# Create new dataframe for just early locations
locf_early <- locf_sf %>% 
  filter(daysSinceStart <= lastDayThreshold)

# Get bearing for the early part of the first trip
locf_early <- locf_early %>% 
  group_by(id) %>% 
  summarise(bearing = bearing(mq, c(last(lon), last(lat))),
            type = last(type), 
            dist2col = last(dist2col))

## plot compass
# ggplot(data = locf_early, aes(bearing, dist2col)) +
#   geom_segment(aes(xend = bearing, yend = 0, colour = type)) +
#   scale_x_continuous(limits = c(-180, 180),
#                      breaks =  seq(-180, 180, 90)) +
#   coord_polar(start = pi) +
#   theme_bw()


# 3) Particle trace ----------------------------------------------------------
pt <- readRDS("./Output/currently_particleTrace.rds")

pt <- pt %>% group_by(id) %>% 
  filter(!is.na(x)) %>% 
  mutate(daysSinceDeparture = difftime(date, min(date), units = 'days'))

# # plot trace
# pt %>% 
#   ggplot() + 
#   geom_point(size = 0.1, aes(lon, lat, col = daysSinceDeparture %>% as.numeric, fill = NULL, group = group)) + 
#   annotate(geom = 'point', x =158.95, y = -54.5, shape = 17, col = 'red', size = 5) +
#   labs(caption = 'Particle trace of 44 seals from MI', col = 'Days since departure') + 
#   viridis::scale_color_viridis() + 
#   lims(x = c(155, 180))  +
#   theme_bw()

# convert trace to sf
pt_sf <- convert2polarsf(pt) %>% 
  left_join(pt) %>% 
  mutate(dist2col = SGAT::gcDist(mq, cbind(lon, lat)))

pt_sf %>% 
  group_by(group) %>% 
  summarise(min(daysSinceDeparture))


# ~ plot weaner tracks / particle trace -----------------------------------
bb <- extent(locw1)
orsi_sf$nudge_x <- c(-1000, -4800, -4500, -4000)
orsi_sf$nudge_y <- c(0, 5300, 4800, 3000)
# png('weaner_particleTrace_track.png', width=20, height=15, units='cm', res=500)
# ggplot() + 
#   geom_sf(data = pt_sf, col = 'blue', size = 0.1) + 
#   geom_sf(data = locw1, size = 0.1, col = 'red') + 
#   labs(caption = 'blue = particle trace, red = weaner. complete 1st trip', x = 'long', y = 'lat') +
#   geom_sf(data = orsi_sf, linetype = 'dashed', colour = 'grey10') + 
#   geom_sf_text(data = orsi_sf, aes(label = front), nudge_x = orsi_sf$nudge_x, nudge_y = orsi_sf$nudge_y) + 
#   xlim(bb[1], bb[2]) + ylim(bb[3]-200, bb[4]+500) +
#   theme_bw() #+ facet_wrap(~id)
# dev.off()


# ~ calc bearing ------------------------------------

## Determine where to calculate bearing from for particle trace based on
## corresponding seal dist2col lastDayThreshold
pt_sf <- pt_sf %>% left_join(as.data.frame(loc_early) %>% 
                               dplyr::select(id, dist2col) %>% 
                               rename(dist2col.w = dist2col))

# pt_sf %>% 
#   filter(dist2col <= dist2col.w) %>% 
#   ggplot() + 
#   geom_sf(size = 0.1)

pt_early <- pt_sf %>% 
  group_by(id) %>% 
  filter(dist2col <= dist2col.w) %>% 
  summarise(bearing = bearing(mq, c(last(lon), last(lat))),
            type = NA, 
            dist2col = last(dist2col))

# Circular stats ----------------------------------------------------------
circ.f <- circular(locf_early$bearing%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
# plot(circ.f, stack = TRUE, shrink = 1, main = 'females')
# arrows.circular(mean(circ.f))
# mean(circ.f)

circ.w <- circular(loc_early$bearing%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
# plot(circ.w, stack = TRUE, shrink = 1.6, main = 'weaners')
# arrows.circular(mean(circ.w))

circ.p <- circular(pt_early$bearing%%360, units = 'degrees', zero = pi/2, rotation = 'clock')
# plot(circ.p, stack = TRUE, shrink = 1.6, main = 'particle trace')
# arrows.circular(mean(circ.p))


# ~ plot circular histogram -------------------------------------------------
p1 <- ggplot() +
  geom_histogram(data = data.frame(circ.f), aes(x = circ.f), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.f), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.f), y = 20, label = mean(circ.f) %>% round(1)) +
  labs(title = 'adult females', y = "number of individuals") + 
  theme_bw() 



p2 <- ggplot() +
  geom_histogram(data = data.frame(circ.w), aes(x = circ.w), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "red") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.w), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.w), y = 20, label = mean(circ.w) %>% round(1)) + 
  labs(title = 'weaners', y = "") + 
  theme_bw()


p3 <- ggplot() +
  geom_histogram(data = data.frame(circ.p), aes(x = circ.p), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "red") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.p), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.p), y = 20, label = mean(circ.p) %>% round(1)) + 
  labs(title = 'particle trace', y = "", 
       caption = paste('bearing calculated based on\nlast location of weaner track day', lastDayThreshold, sep = ' ')) + 
  theme_bw()

# png('dispersal_direction.png', width=20, height=10, units='cm', res=500)
ggarrange(p1, p2, p3, ncol = 3, labels = c('a', 'b', 'c'), align = 'hv')
# dev.off()




# # ~ tests -------------------------------------------------------------------
# 
# ### test if data follows a von Mises distribution (p < 0.05)
# watson.test(circ.f)
# watson.test(circ.w)
# 
# 
# ### compare means of females and weaners
# # https://www.biorxiv.org/content/10.1101/2021.03.25.436932v1.full
# mean(circ.f); mean(circ.w)
# # watson.williams.test(list(females = circ.f, weaners = circ.w))
# # watson.wheeler.test(list(females = circ.f, weaners = circ.w))
watson.two.test(circ.f, circ.w)
watson.two.test(circ.p, circ.w)
# watson.wheeler.test(list(circ.w, circ.p))
# 
# 
# 
# ### Test for direction uniformity 
# # Hermans-Rasson test for uniformity (alternative to Rayleigh test). Better for multimodel distributions
# # https://bmcecol.biomedcentral.com/articles/10.1186/s12898-019-0246-8
# # p < 0.05 = directions are not uniform ie animals have a preferred direction.
# library(CircMLE)
# HR_test(circ.f) 
# HR_test(circ.w)
# 
# # Our data is unimodal - so rayleigh test is fine. Either way both tests p < 0.05
# rayleigh.test(circ.f)
# rayleigh.test(circ.w)
# 
# # save.image(file = "09_bearing_workspace_18Dec21.Rdata")
