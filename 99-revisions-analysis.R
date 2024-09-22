Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/terra/proj") # fix Cannot find proj.db warning
# Read in nc file ---------------------------------------------------------
library(ncdf4)
# library(raster)
library(sf)
library(tidyverse)
library(terra)

# Get file
dir(full.names = T, pattern = "cmems")
file <- dir(full.names = T, pattern = "cmems")[3]
file
# Create a raster brick from the nc file with uo and vo variables
vo <- rast(file, subds = "vo")
uo <- rast(file, subds = "uo")

# Resample to 0.25deg or load existing files
template <- rast(file, subds = "vo")
res(template) <- 0.25

vo1_filename <- "vo_186m_0.25deg.tif"
uo1_filename <- "uo_186m_0.25deg.tif"

if (file.exists(vo1_filename)) {
  vo1 <- rast(vo1_filename)
} else {
  vo1 <- resample(vo, template, filename = vo1_filename)
}

if (file.exists(uo1_filename)) {
  uo1 <- rast(uo1_filename)
} else {
  uo1 <- resample(uo, template, filename = uo1_filename, overwrite = TRUE)
}

source("particle_trace_custom.R")

d <- readRDS("./Output/tracks_processed_12h.rds")
st_geometry(d) <- NULL
d <- d %>% dplyr::filter(trip == 1, SUS == FALSE)
mq <- cbind(158.95, -54.5)

d1 <- d %>%
  group_by(id) %>%
  summarise(startdate = first(date), startlon = first(lon), startlat = first(lat), tripdur = first(tripdur))
d1
out <- NULL

# i <- 1
for (i in 2:2) {
  print(paste(i, "out of", nrow(d1), sep = " "))
  pts <- cbind(d1$startlon[i], d1$startlat[i])
  sdate <- d1$startdate[i] %>% as.Date()
  duration <- d1$tripdur[i]
  id <- d1$id[i]

  # custom particle trace function to accept my own curr_u and curr_v as input
  pt <- particle_trace_custom(
    xy = pts,
    time_step = 24 * 3600,
    start_date = sdate,
    end_date = sdate + duration,
    curr_u = uo1,
    curr_v = vo1,
  )
  pt$group <- i
  pt$id <- id
  pt$date_processed <- Sys.Date()
  out <- out %>% bind_rows(pt)
}


# save to rds
# saveRDS(out, "output/currently_particleTrace_0.49m_0.25deg.rds")

# load rds
old <- readRDS("output/currently_particleTrace.rds") # this was my original anlaysis results using raadtools in currently.
old <- old %>% dplyr::filter(id %in% out$id)

out <- out %>%
  group_by(id) %>%
  mutate(daysSinceDeparture = difftime(date, min(date), units = "days")) %>%
  rename(lon = x, lat = y)

p <- ggplot() +
  geom_path(
    data = out, aes(lon, lat, fill = NULL, group = group),
    col = "black"
  ) +
  geom_path(
    data = old, aes(lon, lat, fill = NULL, group = group),
    col = "blue"
  ) +
  annotate(geom = "point", x = 158.95, y = -54.5, shape = 17, col = "red", size = 5) +
  viridis::scale_color_viridis() +
  # lims(x = c(155, 180))  +
  theme_bw()

p
