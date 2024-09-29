

#   -----------------------------------------------------------------------
#   Dispersal bearing analysis for juvenile seals
#   -----------------------------------------------------------------------

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


# 4) Juveniles ------------------------------------------------------------
# load juvenile data
juv <- readRDS('./output/macca_juvenile_locs.rds')
juv_sf <-  convert2polarsf(juv) %>% 
  left_join(juv)
juv_sf <- juv_sf %>% arrange(id, date)

# create days since start of trip, and dist2col variables
juv_sf <- juv_sf %>% 
  group_by(id) %>% 
  mutate(daysSinceStart = difftime(date, min(date), units = 'day') %>% round(),
         dist2col = SGAT::gcDist(mq, cbind(lon, lat)))

# plot tracks
juv_sf %>%
  ggplot() + 
  geom_sf(aes(colour = daysSinceStart <= REF_DATE), size = 0.01) + 
  geom_sf(data = mq_sf, col = 'red')

# closer look at starting locations
# juv_sf %>% 
#   filter(daysSinceStart <= 4) %>% 
#   ggplot() + 
#   geom_sf(aes(col = daysSinceStart <= 0))


# which seals first dist2col is very far from colony?
juv_sf %>% 
  group_by(id) %>% 
  summarise(dist2col = first(dist2col)) %>% 
  ggplot(aes(x = dist2col)) + 
  geom_histogram() 

# Create new dataframe for just early locations
juv_early <- juv_sf %>% 
  filter(daysSinceStart <= REF_DATE)

# Get bearing for the early part (4 days) of the first trip
juv_early <- juv_early %>% 
  group_by(id) %>% 
  summarise(startdate = min(date), 
            enddate = max(date), 
            bearing = bearing(mq, c(last(lon), last(lat))),
            type = last(type), 
            dist2col = last(dist2col),
            year_born = first(year_born),
            age = first(age))

juv_early$enddate %>% month %>% unique
juv_early$enddate %>% year %>% unique

## summary of methods
juv_early %>% 
  mutate(year = year(startdate)) %>% 
  group_by(age) %>% 
  summarise(n = length(unique(id)))

## mean departure date
juv_early %>% 
  group_by(id) %>% 
  summarise(startdate = min(startdate)) %>% 
  pull(startdate) -> jdates

purrr::map(jdates, function(x) {
  ifelse(month(x) != 12, format(x, '2000-%m-%d'), format(x, '1999-%m-%d'))
}) %>% unlist() %>% as.Date() %>% median()

purrr::map(jdates, function(x) {
  ifelse(month(x) != 12, format(x, '2000-%m-%d'), format(x, '1999-%m-%d'))
}) %>% unlist() %>% as.Date() %>% sd()



# circular stats ----------------------------------------------------------

circ.j1 <- circular(juv_early[juv_early$age == 1,]$bearing%%360, units = 'degrees', zero = pi/2, rotation = 'clock')

circ.j2 <- circular(juv_early[juv_early$age == 2,]$bearing%%360, units = 'degrees', zero = pi/2, rotation = 'clock')

circ.j3 <- circular(juv_early[juv_early$age == 3,]$bearing%%360, units = 'degrees', zero = pi/2, rotation = 'clock')

circ.j4 <- circular(juv_early[juv_early$age == 4,]$bearing%%360, units = 'degrees', zero = pi/2, rotation = 'clock')

# # plot with circular method
# plot(circ.f, stack = TRUE,
#      main = 'Dispersal from Macquarie Island',
#      sub = 'black = adult females, red = weaners')
# arrows.circular(mean(circ.f))
# points(circ.w, stack = TRUE,col = 'red')
# arrows.circular(mean(circ.w), col = 'red')



# Plots -------------------------------------------------------------------

# juvenile age 1
p4 <- ggplot() +
  geom_histogram(data = data.frame(circ.j1), aes(x = circ.j1), 
                 breaks = seq(0, 360, 45), 
                 fill = 'grey',
                 colour = "black") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.j1), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.j1), y = 1, label = mean(circ.j1) %>% round(1)) +
  labs(title = paste('juvenile age 1 n=', length(circ.j1), sep = ''), y = "") + 
  theme_bw()

# juvenile age 2
p5 <- ggplot() +
  geom_histogram(data = data.frame(circ.j2), aes(x = circ.j2), 
                 breaks = seq(0, 360, 45), 
                 fill = 'grey',
                 colour = "black") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.j2)%%360, color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.j2)%%360, y = 3, 
           label = mean(circ.j2)%%360 %>% round(1)) +
  labs(title = paste('juvenile age 2 n=', length(circ.j2), sep = ''), y = "")+ 
  theme_bw()

# juvenile age 3
p6 <- ggplot() +
  geom_histogram(data = data.frame(circ.j3), aes(x = circ.j3), 
                 breaks = seq(0, 360, 45), 
                 fill = 'grey',
                 colour = "black") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.j3)%%360, color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.j3)%%360, y = 3, 
           label = mean(circ.j3)%%360 %>% round(1)) +
  labs(title = paste('juvenile age 3 n=', length(circ.j3), sep = ''), y = "") + 
  theme_bw()

# juvenile age 4
p7 <- ggplot() +
  geom_histogram(data = data.frame(circ.j4), aes(x = circ.j4), 
                 breaks = seq(0, 360, 45), 
                 fill = 'grey',
                 colour = "black") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.j4)%%360, color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.j4)%%360, y = 3, 
           label = mean(circ.j4)%%360 %>% round(1)) +
  labs(title = paste('juvenile age 4 n=', length(circ.j4), sep = ''), y = "") +
  theme_bw()




