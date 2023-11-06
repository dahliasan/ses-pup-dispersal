# Extract mean environmental conditions for prediction

##### Needs to be run on a VM with raadtools data ########

library(raster)
library(raadtools)
library(lubridate)
library(rasterVis)
library(RNetCDF)
library(rgeos)
library(ggplot2)
library(maptools)
library(zoo)

# Define dates to use

# First and last day of month to cover entire track set
startmonth <- as.Date(as.yearmon("1995-11-29"))
endmonth <- as.Date(as.yearmon("2001-10-03"))
endmonth <- endmonth %m+% months(1) # Add one month with lubridate to cover period
dates_month <- seq(startmonth, endmonth, by = "1 month")
dates_day <- seq(startmonth, endmonth, by = "1 day")

# Select only months there seals have gone the furthest.
smr <- c(12, 1, 2, 3, 4, 5, 6)
summer_dates_month <- dates_month[month(dates_month) %in% smr]
summer_dates_day <- dates_day[month(dates_day) %in% smr]

# Define spatial extent of the study area
# Generate a template raster
slim <- c(95, 180, -79, -32) # Use pre-defined limits
slim2 <- c(-180, -120, -79, -32) # Use pre-defined limits

# save.image('envdata_workspace.RData')
wgs <- CRS("+proj=longlat +ellps=WGS84")
temp <- raster(ext = extent(slim), res = 0.1, crs = wgs)

# ---------------------------------------------------
# Get the variables

#### 1. etopo1 = 1 arcminute (0.016 degrees)
topo <- readtopo("etopo1", xylim = slim)
topo <- resample(topo, temp, method = "bilinear")

#### 2. sst
sst <- readsst(summer_dates_month, xylim=extent(temp), lon180 = TRUE, time.resolution = "monthly")
sst <- mean(sst, na.rm = T)
sst <- resample(sst, temp, method = "bilinear")

#### 3. ssha
ssha <- readssh(summer_dates_day, xylim=extent(temp), lon180 = TRUE, ssha = TRUE)
ssha <- mean(ssha, na.rm = T)
ssha <- resample(ssha, temp, method = "bilinear")

#### 4. curru - horizontal geostrophic velocity - delayed-time, multi-altimeter satellite, gridded vertical geostrophic velocity (default=daily)
#### 5. currv - vertical geostrophic velocity - delayed-time, multi-altimeter satellite, gridded vertical geostrophic velocity (default=daily)
#### 6. eke - eddy kenetic energy
curru <- readcurr(summer_dates_day, xylim=extent(temp), lon180 = TRUE, uonly = TRUE)
currv <- readcurr(summer_dates_day, xylim=extent(temp), lon180 = TRUE, vonly = TRUE)

curru <- mean(curru, na.rm = T)
curru <- resample(curru, temp, method = "bilinear")

currv <- mean(currv, na.rm = T)
currv <- resample(currv, temp, method = "bilinear")

eke <- 0.5*((currv)^2 + (curru)^2)

#### 7. mixed layer depth - only monthly climatology available (there's only 12 layers for that)
mld <- readmld(summer_dates_day, xylim=extent(temp), lon180 = TRUE)
mld <- mean(mld, na.rm = T)
mld <- resample(mld, temp, method = "bilinear")

#### 8. windu - horizontal wind (default=6h)
#### 9. windv - vertical wind (default=6h)
windu <- readwind(summer_dates_day, xylim=extent(temp), lon180 = TRUE, uonly = TRUE)
windv <- readwind(summer_dates_day, xylim=extent(temp), lon180 = TRUE, vonly = TRUE)

windu <- mean(windu, na.rm = T)
windu <- resample(windu, temp, method = "bilinear")

windv <- mean(windv, na.rm = T)
windv <- resample(windv, temp, method = "bilinear")

#### 10. tri - terrain ruggedness index DOI 10.1080/01490410701295962
bath <- readbathy(topo="etopo1", xylim=slim)
tri <- terrain(bath, opt="TRI")
tri <- resample(tri, temp, method = "bilinear")

### 11. slope
slope <- terrain(bath, opt='slope', unit='degrees') #bath read in above (for tri)
slope <- crop(slope, temp)
slope <- resample(slope, temp, method = "bilinear")

#### 12. SSTgrad
readsstgrad <- function(date, returnfiles = FALSE, ...) {
  xx <- readsst(date, returnfiles = returnfiles, ...)
  if (returnfiles) return(xx)
  terrain(xx, opt = "slope")
}

SSTgrad <- readsstgrad(summer_dates_month, xylim=extent(temp), lon180 = TRUE, time.resolution = "monthly")
SSTgrad <- mean(SSTgrad, na.rm = T)
SSTgrad <- resample(SSTgrad, temp, method = "bilinear")

#### 13. SSHgrad
readsshgrad <- function(date, returnfiles = FALSE, ...) {
  xx <- readssh(date, returnfiles = returnfiles, ...)
  if (returnfiles) return(xx)
  terrain(xx, opt = "slope", unit = "degrees")
}

SSHgrad <- readsshgrad(summer_dates_day, xylim=extent(temp), lon180 = TRUE)
SSHgrad <- mean(SSHgrad, na.rm = T)
SSHgrad <- resample(SSHgrad, temp, method = "bilinear")

### 14. ice
ice <- readice(summer_dates_day, xylim=extent(temp), lon180 = TRUE, icea = TRUE)
ice <- mean(ice, na.rm = T)
ice <- resample(ice, temp, method = "bilinear")

### 15. chl
chl <- readchla(summer_dates_month, product='MODISA', algorithm='johnson', xylim=extent(temp))
chl <- mean(chl, na.rm = T)
chl <- resample(chl, temp, method = "bilinear")

## Combine into a single stack

# If some rasters were set to NULL
vars <- c("topo",
          "sst",
          "ssha",
          "curru",
          "currv",
          "eke",
          "mld",
          "windu",
          "windv",
          "tri",
          "slope",
          "SSTgrad",
          "SSHgrad",
          "ice",
          "chl")

envdata <- stack()

for (i in 1:length(vars)) {
  print(vars[i])
  if(exists(vars[i])) {
    envdata <- stack(envdata, get0(vars[i]))
  }
}

names(envdata) <- vars

## Mask based on land
# msk <- bath
# rm(bath)
# msk[msk > 0] <- NA
# envdata <- mask(x = envdata, mask = msk)

plot(envdata, col=viridis(125))

save(envdata, file="./output/envdata.Rdata")
#writeRaster(envdata, "envdata.grd", format = "raster")
