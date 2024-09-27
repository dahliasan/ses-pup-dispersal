## download and process data from copernicus
library(lubridate)
library(dplyr)
library(raster)
library(stringr)
rm(list = ls())


# Download data -----------------------------------------------------------

pw = "'soAFO9UK1HllL22%fa%j'"
user = "'dfoo1'"
dates = c(seq(ymd_hms("1999-11-30 12:00:00"), ymd_hms("2002-01-06 12:00:00"), '60 days'), ymd_hms("2002-01-06 12:00:00")) %>% 
  as.character() %>% 
  paste("'", ., "'", sep = '')



dir = "'/Users/dahliafoo/Downloads'"

## Longitude 95-180 (east)

for(i in 4:length(dates)-1) {
  filename = paste("'seapodym_1999-2002", '(', i, ").nc'", sep = '')
  command = paste("/Users/dahliafoo/opt/anaconda3/bin/python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id GLOBAL_MULTIYEAR_BGC_001_033-TDS --product-id cmems_mod_glo_bgc_my_0.083deg-lmtl_PT1D-i --longitude-min 95 --longitude-max 180 --latitude-min -79 --latitude-max -30 --date-min", dates[i], "--date-max", dates[i+1],  "--variable mnkc_hmlmeso --variable mnkc_lmeso --variable mnkc_mlmeso --variable mnkc_mumeso --variable mnkc_umeso --out-dir", dir,  "--out-name", filename, "--user", user, "--pwd", pw)
  
  system(command)
  
  print(paste(i, 'of', length(dates), sep = ' '))
  
}

## Longitude -180 to -120 (west)

for(i in 7:(length(dates)-1)) {
  filename = paste("'seapodym_1999-2002_2", '(', i, ").nc'", sep = '')
  command = paste("/Users/dahliafoo/opt/anaconda3/bin/python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id GLOBAL_MULTIYEAR_BGC_001_033-TDS --product-id cmems_mod_glo_bgc_my_0.083deg-lmtl_PT1D-i --longitude-min -180 --longitude-max -120 --latitude-min -79 --latitude-max -30 --date-min", dates[i], "--date-max", dates[i+1],  "--variable mnkc_hmlmeso --variable mnkc_lmeso --variable mnkc_mlmeso --variable mnkc_mumeso --variable mnkc_umeso --out-dir", dir,  "--out-name", filename, "--user", user, "--pwd", pw)
  
  system(command)
  
  print(paste(i, 'of', length(dates), sep = ' '))
  
}


#  Process data -----------------------------------------------------------

## Set the path for the NetCDF file
ncfiles <- dir('./Data/seapodym-copernicus', full.names = T, pattern = '.nc')
ncfiles1 <- ncfiles[1:13]
ncfiles2 <- ncfiles[14:length(ncfiles)]

q <- list()
varnames <- c('mnkc_umeso', 'mnkc_mumeso', 'mnkc_hmlmeso', 'mnkc_mlmeso', 'mnkc_lmeso')

## merge all nc files

for(j in 2:length(varnames)){
  # stack one section of nc files first (13 files)
  varname <- varnames[j]
  s <- list()
  for(i in 1:length(ncfiles1)){
    s[[i]] <- brick(ncfiles[i], varname = varname)
  }
  
  s <- stack(s)
  rdates <- as.POSIXct(names(s), format = "X%Y.%m.%d.%H.%M.%S", tz = "UTC") %>% as.Date()
  dup <- which(rdates %>% duplicated)
  s <- s[[-dup]]
  rdates <- as.POSIXct(names(s), format = "X%Y.%m.%d.%H.%M.%S", tz = "UTC") %>% as.Date()
  s <- s[[order(rdates)]]
  
  out <- list(s)
  
  # stack other section of ncfiles
  s <- list()
  for(i in 1:length(ncfiles2)){
    s[[i]] <- brick(ncfiles2[i], varname = varname)
  }
  
  s <- stack(s)
  rdates <- as.POSIXct(names(s), format = "X%Y.%m.%d.%H.%M.%S", tz = "UTC") %>% as.Date()
  dup <- which(rdates %>% duplicated)
  s <- s[[-dup]]
  rdates <- as.POSIXct(names(s), format = "X%Y.%m.%d.%H.%M.%S", tz = "UTC") %>% as.Date()
  s <- s[[order(rdates)]]
  
  out[[2]] <- s
  
  ## merge rasterbricks
  s <- merge(out[[1]], out[[2]])
  rdates <- as.POSIXct(names(out[[1]]), format = "X%Y.%m.%d.%H.%M.%S", tz = "UTC") %>% as.Date()
  names(s) <- rdates
  
  ## save to final output list
  q[[varname]] <- s
  
  print(paste(j, 'of', length(varnames)))
  
}

saveRDS(q, file = './output/seapodym_list.rds')











# Load location dataset
load("./Output/raadtools_extract_sim_tracks_12h.rdata")

# Extract chla

ddates <- unique(d$day)
rdates <- as.POSIXct(names(chl), format = "X%Y.%m.%d.%H.%M.%S", tz = "UTC") %>% as.Date %>% sort()
r <- chl[[which(rdates %in% ddates)]]
rdates <- as.POSIXct(names(r), format = "X%Y.%m.%d.%H.%M.%S", tz = "UTC") %>% as.Date %>% sort()

d1 <- NULL
if(all.equal(ddates, rdates)) 
  for(i in 1:length(ddates)){
    df1 <- d[d$day==ddates[i],]
    r1 <- r[[i]]
    # extract chl
    df1$chl <- raster::extract(r1, df1[, c("lon", "lat")])
    # extract chlgrad
    df1$chlgrad <- raster::extract(terrain(r1, opt = 'slope'), df1[, c("lon", "lat")])
    # combine
    d1 <- rbind(d1, df1)
    print(paste(i,'of',length(ddates)))
  }

d <- d1
d <- d %>% arrange(id, sim, date)
# save(d, file = "./Output/raadtools_extract_sim_tracks_12h_chla.rdata")




