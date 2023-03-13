library(raster)
library(raadtools)
library(lubridate)
library(rasterVis)
library(rts)
library(RNetCDF)
library(rgeos)
library(ggplot2)
library(maptools)
library(zoo)

# library(help = raadtools)


load('simulated_tracks_12h.RData')

# Extraction code for 17 environmental covariates 
# - mainly remotely-sensed by satellites - 
# using ‘raadtools’ (Sumner 2016) and ‘raster’ (Hijmans 2016) 
coords <- d[ ,c("lon", "lat")]
lims <- c(min(d$lon)-2, max(d$lon)+2, min(d$lat)-2, max(d$lat)+2) # extract environmental raster for tracks (obs and sims)

# first and last day of month to cover entire track set
startmonth <- as.Date(as.yearmon(min(d$date)))
endmonth <- as.Date(as.yearmon(max(d$date)))
endmonth <- endmonth %m+% months(1) #add one month with lubridate to cover period
dates_month <- seq(startmonth, endmonth, by = "1 month")
dates_day <- seq(startmonth, endmonth, by = "1 day")


# now extract variables
# --------------------------------------------------------------

#### 1. etopo1 = 1 arcminute (0.016 degrees)
d$topo <- extract(readtopo("etopo1"), d[, c("lon", "lat")]) 


#### 2. sst
# https://rdrr.io/github/AustralianAntarcticDivision/raadtools/man/readsst.html
d$sst <- extract(readsst, d[, c("lon", "lat", "date")])


#### 3. ssh - delayed-time, multi altimeter satellite, gridded ssh computed with respect to a 20-year mean (default=daily)
d$ssha <- extract(readssh, d[, c("lon", "lat", "date")], ssha=TRUE) 


#### 4. curru- horizontal geostrophic velocity - delayed-time, multi-altimeter satellite, gridded vertical geostrophic velocity (default=daily)
#### 5. currv - vertical geostrophic velocity - delayed-time, multi-altimeter satellite, gridded vertical geostrophic velocity (default=daily)
#### 6. eddy kenetic energy
d$curru <- extract(readcurr, d[, c("lon", "lat", "date")], uonly=TRUE) 
d$currv <- extract(readcurr, d[, c("lon", "lat", "date")], vonly=TRUE)
d$eke <- 0.5*((d$currv)^2 + (d$curru)^2)


#### 7. mixed layer depth - only monthly climatology available (there's only 12 layers for that)
# d$mld <- extract(readmld, d[, c("lon", "lat", "date")])


#### 8. windu - horizontal wind (default=6h)
#### 9. windv - vertical wind (default=6h)
d$windu <- extract(readwind, d[, c("lon", "lat", "date")], uonly=TRUE) 
d$windv <- extract(readwind, d[, c("lon", "lat", "date")], vonly=TRUE)


#### 10. TRI - terrain ruggedness index DOI 10.1080/01490410701295962
bath <- readbathy(topo="etopo1", xylim=lims)
tri <- terrain(bath, opt="TRI")
d$tri <- extract(x=tri, y=coords)
rm(tri)


#### 11. SLOPE
slope <- terrain(bath, opt='slope', unit='degrees') #bath read in above (for TRI)
d$slope <- extract(x=slope, y=coords)
rm(slope, bath)


#### 12. SSTgrad
readsstgrad <- function(date, returnfiles = FALSE, ...) {
  xx <- readsst(date, returnfiles = returnfiles, ...)
  if (returnfiles) return(xx)
  terrain(xx, opt = "slope")
}
d$SSTgrad <- extract(readsstgrad, d[, c("lon", "lat", "date")])
plot(d$SSTgrad)


#### 13. SSHgrad
readsshgrad <- function(date, returnfiles = FALSE, ...) {
  xx <- readssh(date, returnfiles = returnfiles, ...)
  if (returnfiles) return(xx)
  terrain(xx, opt = "slope", unit = "degrees")
}
d$SSHgrad <- extract(readsshgrad, d[, c("lon", "lat", "date")])
plot(d$SSHgrad)

save.image("raadtools_extract_enviro_WORKSPACE.RData")


#### 14. ICE
# https://rdrr.io/github/AustralianAntarcticDivision/raadtools/man/readice.html
#d$ice <- extract(readice_monthly, d[, c("lon", "lat", "date")])
d$ice <- extract(readice, d[, c("lon", "lat", "date")])
d$ice[is.na(d$ice)] <- 0
plot(d$ice)


#### 15. DISTice distance to ice
#d$DISTice <- extract(distance_to_ice_edge, d[, c("lon", "lat", "date")], threshold=15, hemisphere = "south")
#plot(d$DISTice)
## d1 is the data frame of track data
library(sospatial)
library(dplyr)
library(sf)
## expand this summary data set into the longitude/latitude of daily
## sea ice edge
ice <- tibble(lon = rep_len(seq(-180, 179), length.out = length(packed_lats)), 
              lat = packed_lats/10,
              day = as.Date(rep(ice_dates, each = 360)))
range(ice$day)
#[1] "1978-10-26 10:00:00 AEST" "2017-11-09 11:00:00 AEDT"

## split our track data set into whole days
d1 <- d
d1$day <- as.Date(round(d1$date, "day"))
list_d1 <- split(d1, d1$day)

for (iday in seq_along(list_d1)) {
  day1 <- list_d1[[iday]]
  ## this day's ice summary
  ice1 <- ice %>% dplyr::filter(day == day1$day[1])
  ## workaround for missing lat
  ice1$lat[1] <- mean(ice1$lat[c(2, nrow(ice1))])
  ice_lats <- st_as_sf(ice1, coords = c("lon", "lat"), crs = 4326)
  trk_lats <- st_as_sf(day1, coords = c("lon", "lat"), crs = 4326)
  list_d1[[iday]]$dist_to_ice_m <- c(apply(st_distance(trk_lats, ice_lats), 1, min))
}
d2 <- bind_rows(list_d1)
d <- d2
rm(d2, d1)
save.image("raadtools_extract_enviro_WORKSPACE.RData")
#library(ggplot2)
#ggplot(d1, aes(lon, lat, colour = dist_to_ice_m)) + geom_point()


#### 16. CHL - talk to Mike - read monthly mean by first giving every date in the month (doesn't work like other funcions)
# raadtools works in this form but will be deprecated...
# Also note, using seawifs not ideal??
# d$weeks <- as.character(cut(d$date, "1 week"))
# dates <- sort(as.POSIXct(as.character(unique(d$weeks)), tz="GMT"))
# chl <- readchla(dates[c(200)], product='MODISA', algorithm='johnson')
# d$chl <- extract(x=chl, y=d[, c("lon", "lat")])
# 
# d1 <- NULL
# date_list <- substring(as.character(dates),1,10)
# for(i in 142:length(dates)) { # data available from Sep 2002 onwards
#   df1 <- d[d$weeks==date_list[i],]
#   chl <- readchla(dates[i], product='MODISA', algorithm='johnson')
#   df1$chl <- extract(x=chl, y=df1[, c("lon", "lat")])
#   d1 <- rbind(d1, df1)
#   print(paste(i,'of',length(dates)))
# }
# plot(d$chl)
# ###


save.image("raadtools_extract_enviro_WORKSPACE.RData")

#### 17. distance to colony
#derivaadcproducts() # list all derived data products available through readderivaadc
#d1$DISTcolony <- extract(readderivaadc, products='distance_colony', d[, c("lon", "lat", "date")])
####


save(d, file='raadtools_extract_sim_tracks.RData')