
rm(list = ls())
library(tidyverse)
library(lubridate)
library(foieGras)
library(sf)
library(sp)
library(raster)
library(viridis)


# Step 1: Prepare data ------------------------------------------------------------
load('baseInfo.RData')

w1 <- read_csv("./Data/mq-ellie-weaners-argos.csv")
w1 <- w1 %>% 
  mutate(date = d.date %>% mdy_hms())

w1 <- w1 %>% filter(ref %in% id_keep)

w2 <-  w1 %>%
  dplyr::select(ref, date, lq, lon, lat) %>%
  dplyr::rename(id=ref, lc=lq)

## convert to foisGras format
w2 <- w2 %>% mutate(lc = replace(lc, lc == -9, "Z"), lc = replace(lc, lc == -2, "B"), lc = replace(lc, lc == -1, "A"))

## only keep seals with at least 10 days of data
dur <- w2 %>%
  group_by(id) %>%
  summarise(first = min(date), last = max(date)) %>%
  mutate(dur = difftime(last, first,units = 'days')) %>%
  filter(dur >= 10) %>% 
  pull(id)

##ensure the data are sorted by id and date and remove duplicates and locations too far north
d1 <- w2 %>%
  # keep only seals with >= 10 days deployment
  filter(id %in% dur) %>%
  arrange(id, date) %>%
  group_by(id) %>%
  # remove near duplicates 
  mutate(timediff = difftime(lead(date), date, units = 'min')) %>%
  filter(timediff > 2) %>% 
  # remove locations too far north
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
d3 <- as_tibble(d1) 

table(d3$id)
range(d3$date)

##time step of 24 = 1 per day
fit_all <- foieGras::fit_ssm(d3,
                             vmax = 4, 
                             map = list(psi = factor(NA)),
                             model = "crw",
                             time.step = 12)

beepr::beep(4)

plot(fit_all, what = "fitted", type = 2)

## save the output
# saveRDS(fit_all, file = './Output/foiegras_fitssm_vmax4_12h.rds')

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



# Step 3. Infer behavioural index -----------------------------------------

# Encountered issues fitting individuals all at once as a jmpm. 
# Solution is to fit them individually as mpm. See https://github.com/ianjonsen/foieGras/issues/29

fit_all <- readRDS("./Output/foiegras_fitssm_vmax4_12h.rds")

fmp_all <- purrr::map(split(fit_all, 1:nrow(fit_all)), function(x) {
  fit_mpm(x, what = "predicted", model = "mpm", control = mpm_control(verbose = 0))
}) %>% 
  reduce(bind_rows)

beepr::beep(4)

# which failed?
fmp_all %>% filter(converged == FALSE) 


# ~ Compare fmp fitted vs predicted -----------------------------------------
# Fit fmp to fitted locations

fmp_fitted <- purrr::map(split(fit_all, 1:nrow(fit_all)), function(x) {
  fit_mpm(x, what = "fitted", model = "mpm", control = mpm_control(verbose = 0))
})

beepr::beep(4)

fmp_fitted <- fmp_fitted %>% reduce(bind_rows)

beepr::beep(4)
# plot(fmp_fitted, rev = TRUE)

f <- grab(fmp_fitted, "fitted", as_sf=FALSE) %>% 
  mutate(type = 'fitted')

tmp <- grab(fmp_all, "fitted", as_sf=FALSE) %>% 
  mutate(type = 'vmax4')
f <- f %>% bind_rows(tmp)

# tmp <- grab(fmp_all, "fitted", as_sf=FALSE) %>% 
#   mutate(type = 'vmax10')
# f <- f %>% bind_rows(tmp)
# 
# tmp <- grab(fmp_all, "fitted", as_sf=FALSE) %>% 
#   mutate(type = 'vmax3.5')
# f <- f %>% bind_rows(tmp)

f %>% 
  ggplot(aes(x = date, y = g)) + 
  geom_line(aes(colour = type)) + 
  facet_wrap(~id, scale = 'free') # vmax = 4 is best so far. 


# save(fmp_all, file = './output/foiegras_fitmpm_vmax4_12h.Rdata')


# Plot fmp values ---------------------------------------------------------

f <- grab(fmp_all, "fitted", as_sf=FALSE) %>% 
  group_by(id) %>% 
  mutate(daysSinceDeparture = difftime(date, min(date), units = 'days'))
f %>% 
  ggplot(aes(x = daysSinceDeparture, y = g)) + 
  geom_line(aes(colour = g), size = 1) + 
  facet_wrap(~id) + 
  scale_colour_viridis(option = 'D') +
  # scale_colour_see_c() +
  theme_bw()


