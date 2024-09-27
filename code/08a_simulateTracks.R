rm(list=ls(all=TRUE))

############################# 
# STEP 1: track simulations #
#############################

library(remotes)
# install_github("AustralianAntarcticDivision/availability")
library(availability)
library(tidyverse)
library(crawl)
library(geosphere)
library(sf)

# Using the vector-AR method (see [1]):

# estimated animal locations at regular time intervals, 
# while accounting for errors around Argos locations by 
# fitting a continuous-time correlated random walk model 
# (Johnson et al. 2008) to each track using the ‘crawl’ 
# 1.5 package (Johnson 2015). We used Jonsen's bsam approach.

# set the time interval as the whole number nearest the median 
# time interval in the unprocessed tracking data

# load('./Output/tracks_processed.RData')
d1 <- readRDS("./Output/tracks_processed_12h.rds")

## Get only first trips and remove suspicious locations
d1 <- d1 %>% 
  filter(trip == 1, SUS == FALSE) 

# Complete as well as incomplete tracks were retained
st_geometry(d1) <- NULL

# fitting a first-order vector autoregressive model (~15 mins)
d <- NULL
dp <- NULL

id <- unique(d1$id)
for(i in 1:length(id)) {
  
  wd <- "./output/pseudo tracks/track simulation/"
  png(paste(wd,id[i],".png",sep=''), res=300,height=14,width=23,units="cm")
  
  df1 <- d1[d1$id==id[i],]
  with(df1, plot(lon, lat, pch=19, cex=0.6, type='b', col='royalblue'))
  
  # To characterize the environment potentially available to individuals, 
  # and thus allowing a case-control design for habitat preference modelling 
  # (Aarts et al. 2008), we simulated random or pseudo-tracks.  For each 
  # observed track we simulated 20 pseudo-tracks by fitting a first-order 
  # vector autoregressive model characterized by the step lengths and turning
  # characteristics of the real track, as detailed in Raymond et al. (2015). 
  # Pseudo-locations falling on land were rejected and re-sampled.
  for(j in 1:20) {
    realtrack <- as.matrix(df1[,c('lon','lat')]) ## 2-column matrix of longitude and latitude
    arf <- surrogateARModel(realtrack) ## fit AR model to track
    st <- surrogateAR(arf, realtrack, point.check = gshhsMask()) ## simulate new track
    time <- df1$date
    df2 <- data.frame(id=id[i], date=time, lon=st$xs[,1], lat=st$xs[,2], sim=j)
    dp <- rbind(dp, df2) # too slow
  }
  
  library(maps)
  with(dp[dp$id==id[i],], points(lon, lat, pch=19, cex=0.6, type='p', col='lightblue'))
  with(df1, points(lon, lat, pch=19, cex=0.6, type='b', col='royalblue'))
  map("world", add=T, fill=TRUE, col="lightgrey", lwd=1.5)	
  print(paste(i,'of',length(id)))
  dev.off() # scan()
}

# bind original and simulated tracks 
d <- d1 %>% mutate(sim = 0)
d <- bind_rows(d,dp) %>% arrange(id, sim)

d %>% 
  group_by(id) %>% 
  summarise(n(), min(sim), max(sim))

# Inspect real and sim speeds
time_step <- 4
step_len <- function(z) distVincentyEllipsoid(z[-nrow(z), 1:2], z[-1, 1:2]) / 1e3
step_len_wrap <- function(xx, bysim = F) {
  if(bysim == FALSE) {
    x_split <- split(xx, xx$id)
  } else {
      x_split <- split(xx, paste(xx$id, xx$sim, sep = '_'))
    }
  purrr::map(x_split, function(x) {
    yx <- as.matrix(x[,c('lon','lat')])
    step_len(yx) / time_step
  }) %>%
    do.call("c", .)
}

beepr::beep(4)

## For multiple seals
yx <- step_len_wrap(d1)
temp <- data.frame(speed = yx, track = "Observed")
yx <- step_len_wrap(dp, bysim = T)
temp <- rbind(temp, data.frame(speed = yx, track = "AR"))
ggplot(data = temp, aes(x = speed, fill = track)) + 
  geom_histogram(binwidth = 0.5) +
  facet_wrap(~track, ncol = 1, scales = 'free') + 
  scale_x_continuous(breaks = seq(from = 0, to = 60, by = 5)) +
  xlab("Speed (km/h)")

## For single seals 
# yx <- as.matrix(d1[,c('lon','lat', 'id')])
# temp <- data.frame(speed = step_len(yx) / time_step, track = "Observed")
# yx <- as.matrix(dp[,c('lon','lat')])
# temp <- rbind(temp, data.frame(speed = step_len(yx) / time_step, track = "AR"))
# ggplot(data = temp, aes(x = speed, fill = track)) + 
#   geom_histogram(binwidth = 0.5) +
#   facet_wrap(~track, ncol = 1, scales = 'free') + 
#   scale_x_continuous(breaks = seq(from = 0, to = 60, by = 5)) +
#   xlab("Speed (km/h)")

temp %>% 
  group_by(track) %>% 
  summarise(min(speed), max(speed))


# save output
# save(d, dp, file='./output/pseudo tracks/simulated_tracks_12h.RData')


source('convert2polarsf.R')
crs <- "+proj=stere +lon_0=170 +lat_0=-90 +units=km +datum=WGS84"
proj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# visually inspect tracks and sims --------------------------------------
# world map
world_sf <- st_as_sf(rworldmap::getMap(resolution = "li")) %>% 
  st_transform(., crs = crs)
dd <- d
dd <- convert2polarsf(dd)
bb <- extent(dd)
png('./output/pseudo tracks/real_pseudo_tracks.png', width=30, height=20, units='cm', res=500)
ggplot() +
  geom_sf(data = dd, aes(colour = sim < 1), size = 0.01) +
  geom_sf(data = world_sf, fill = grey(0.4), colour = NA) + 
  facet_wrap(~id) + 
  xlim(bb[1], bb[2]) + ylim(bb[3], bb[4])
dev.off()
beepr::beep(4)


# test if any points on land
mask  <- gshhsMask() ## initialize land mask function
x <- split(d, d$id)[[2]]
purrr::map2(x$lon, x$lat, ~mask(0,c(.x, .y)))
