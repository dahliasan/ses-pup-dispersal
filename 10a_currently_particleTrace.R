# particle dispersion
rm(list=ls(all=TRUE))
# devtools::install_github("AustralianAntarcticDivision/raadtools")
# devtools::install_github("mdsumner/currently")


library(raadtools)  ## assuming version >= 0.6.0.9005
library(currently) 
library(tidyverse)
library(sf)

d <- readRDS('./Output/tracks_processed_12h.rds')
d <- d %>% filter(sim == 0) # only first trips
mq <- cbind(158.95, -54.5)

d1 <- d %>% group_by(id) %>% 
  summarise(startdate = first(date), startlon = first(lon), startlat = first(lat), tripdur = first(tripdur))

out <- NULL
for(i in 1:nrow(d1)){
  print(paste(i, "out of", nrow(d1), sep = ' '))
  pts <- cbind(d1$startlon[i], d1$startlat[i])
  sdate <- d1$startdate[i] %>% as.Date()
  duration <- d1$tripdur[i]
  id <- d1$id[i]
  pt <- particle_trace(pts, time_step = 12 * 3600, start_date = sdate, end_date = sdate + duration)
  pt$group <- i
  pt$id <- id
  out <- out %>% bind_rows(pt)
}


# saveRDS(out, file = 'currently_particleTrace.rds')
# save.image('currently_workspace.Rdata')

out <- out %>% 
  group_by(id) %>% 
  mutate(daysSinceDeparture = difftime(date, min(date), units = 'days'))

#ddt <- as.data.frame(setNames(topo, "z"), xy = TRUE)


ggplot() + #geom_raster(data = ddt, aes(x, y, fill = z))
  geom_path(data = out, aes(lon, lat, col = daysSinceDeparture %>% as.numeric, fill = NULL, group = group)) + 
  annotate(geom = 'point', x =158.95, y = -54.5, shape = 17, col = 'red', size = 5) +
  labs(caption = 'Particle trace of 44 seals from MI', col = 'Days since departure') + 
  viridis::scale_color_viridis() + 
  lims(x = c(155, 180))  +
  theme_bw()



# Calculate particle trace for adult females ------------------------------
load("./Data/macca_winter_locs.Rdata") # adult females
loc <- as_tibble(loc)
loc <- loc %>% rename(id = seal, date = gmt)
loc <- loc %>% mutate(dist2col = SGAT::gcDist(mq, cbind(lon, lat)))

# Check if adult female track is one continuous track (foraging trip)
loc %>% 
  ggplot(aes(x = date, y = dist2col)) + 
  geom_point() +
  facet_wrap(~id, scales = 'free')

d1 <- loc %>% group_by(id) %>% 
  summarise(startdate = first(date), startlon = first(lon), startlat = first(lat))

out <- NULL
for(i in 1:nrow(d1)){
  print(paste(i, "out of", nrow(d1), sep = ' '))
  pts <- cbind(d1$startlon[i], d1$startlat[i])
  sdate <- d1$startdate[i] %>% as.Date()
  duration <- 30 # 30 days
  id <- d1$id[i]
  pt <- particle_trace(pts, time_step = 12 * 3600, start_date = sdate, end_date = sdate + duration)
  pt$group <- i
  pt$id <- id
  out <- out %>% bind_rows(pt)
}

