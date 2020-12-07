library(tidyverse)
library(lubridate)
library(foieGras)
library(sf)
library(sp)
library(raster)


# Step 1: Prepare data ------------------------------------------------------------
w1 <- read_csv("./Data/mq-ellie-weaners-argos.csv")
w1 <- w1 %>% 
  mutate(date = d.date %>% mdy_hms())

w2 <-  w1 %>%
  dplyr::select(ref, date, lq, lon, lat) %>%
  dplyr::rename(id=ref, lc=lq)

## convert to foisGras format
w2 <- w2 %>% mutate(lc = replace(lc, lc == -9, "Z"), lc = replace(lc, lc == -2, "B"), lc = replace(lc, lc == -1, "A"))

## only keep seals with 10 days of data
dur <- w2 %>%
  group_by(id) %>%
  summarise(first = min(date), last = max(date)) %>%
  mutate(dur = difftime(last, first)) %>%
  filter(dur > 240) %>%
  pull(id)

##ensure the data are sorted by id and date and remove duplicates and locations too far north
d1 <- w2 %>%
  filter(id %in% dur) %>%
  arrange(id, date) %>%
  group_by(id) %>%
  distinct(date, .keep_all = TRUE) %>%
  filter(lat < -20) %>%
  dplyr::select(1:5)


# Plot raw locations ------------------------------------------------------
############################################################
## initial plots of tracks
crs <- "+proj=stere +lon_0=170 +lat_0=-90 +units=km +datum=WGS84"
coordinates(d1) <- c("lon", "lat")
projection(d1) <-  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

## Set up as a "simple feature" for plotting in ggplot
d <- d1 %>%
  st_as_sf(d1)

## project data to Antarctic Stereographic
d_sf <- d %>% st_transform(crs = crs)

## get high-res land & apply stereo projection
world_sf <- st_as_sf(rworldmap::getMap(resolution = "high")) %>%
  st_transform(., crs = crs)

bb <- extent(d_sf)
raw_plot <- ggplot() +
  geom_sf(data = d_sf, size=.1) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  # xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) + # show extent of all raw locations
  xlim(bb[1]+4500, bb[2]-5000) + ylim(bb[3]+4000, bb[4]-500) # crop out outlier locations

raw_plot


# Step 2. Fit the SSM -----------------------------------------------------
## uses a speed max of 7.0 m/s (approx 25 km/h)
d3 <- as_tibble(d1) 

table(d3$id)
range(d3$date)

##time step of 24 = 1 per day
fit_all <- foieGras::fit_ssm(d3,
                             vmax = 7,
                             map = list(psi = factor(NA)),
                             model = "crw",
                             optim = c("nlminb"),
                             time.step = 4)

## extract the predicted values - ie estimated location every day - in lat/lon format
ssm <- grab(fit_all, "predicted", as_sf=FALSE)
dmp1 <-  ssm

## plot the ssm data
coordinates(dmp1) <- c("lon", "lat")
projection(dmp1) <-  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

ssmp <- dmp1 %>%
  st_as_sf(dmp1)

## project data to Antarctic Stereographic
ssmp_sf <- ssmp %>% st_transform(crs = crs)

ssm_plot <- ggplot() +
  geom_sf(data = ssmp_sf, size=.1) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
  xlim(bb[1]+3500, bb[2]-5000) + ylim(bb[3]+4000, bb[4]-500)
ssm_plot

## save the output
save(fit_all, ssm, file = "./Output/weaner_ssm_4h.Rdata")
