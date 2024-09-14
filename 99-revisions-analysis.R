# Read in nc file ---------------------------------------------------------
library(ncdf4)
library(raster)
library(currently)

# Get file
dir(full.names = T, pattern = "cmems")
file <- dir(full.names = T, pattern = "cmems")[1]

# Load the nc file
nc_data <- nc_open(file)

# Get var names
names(nc_data$var)

# Create a raster brick from the nc file with uo and vo variables
vo <- brick(file, varname = "vo")
uo <- brick(file, varname = "uo")

# plot(uo[[1]])
# plot(vo[[1]])
#
# r <- stack(uo[[1]], vo[[1]])
# plot(r)
#
# extent(uo)

Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/terra/proj") # fix Cannot find proj.db warning
source("particle_trace_subsurface.R")

library(currently)
library(tidyverse)
library(sf)

d <- readRDS("./Output/tracks_processed_12h.rds")
st_geometry(d) <- NULL
d <- d %>% filter(trip == 1, SUS == FALSE)
mq <- cbind(158.95, -54.5)


d1 <- d %>%
  group_by(id) %>%
  summarise(startdate = first(date), startlon = first(lon), startlat = first(lat), tripdur = first(tripdur))

d1

out <- NULL
for (i in 1:nrow(d1)) {
  print(paste(i, "out of", nrow(d1), sep = " "))
  pts <- cbind(d1$startlon[i], d1$startlat[i])
  sdate <- d1$startdate[i] %>% as.Date()
  duration <- d1$tripdur[i]
  id <- d1$id[i]
  pt <- particle_trace_custom(
    xy = pts,
    uo = uo,
    vo = vo,
    time_step = 24 * 3600,
    start_date = sdate,
    end_date = sdate + duration,
  )
  pt$group <- i
  pt$id <- id
  pt$date_processed <- Sys.Date()
  out <- out %>% bind_rows(pt)
}

# saveRDS(out, file = 'currently_particleTrace_200m.rds')
# save.image('currently_workspace.Rdata')

# load rds
old <- readRDS("output/currently_particleTrace.rds")

old <- old %>% filter(id == out$id[1])

out <- out %>%
  group_by(id) %>%
  mutate(daysSinceDeparture = difftime(date, min(date), units = "days"))

ggplot() + # geom_raster(data = ddt, aes(x, y, fill = z))
  # geom_path(data = out, aes(lon, lat, col = daysSinceDeparture %>% as.numeric, fill = NULL, group = group)) +
  geom_path(
    data = out, aes(lon, lat, fill = NULL, group = group),
    col = "black"
  ) +
  geom_path(
    data = old, aes(lon, lat, fill = NULL, group = group),
    col = "blue"
  ) +
  # annotate(geom = 'point', x =158.95, y = -54.5, shape = 17, col = 'red', size = 5) +
  labs(caption = "Particle trace of 44 seals from MI", col = "Days since departure") +
  viridis::scale_color_viridis() +
  # lims(x = c(155, 180))  +
  theme_bw()
