Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/terra/proj") # fix Cannot find proj.db warning
# devtools::install_github("mdsumner/currently")
source("99-revisions-particle_trace_custom.R")
# Read in nc file ---------------------------------------------------------
library(ncdf4)
# library(raster)
library(sf)
library(tidyverse)
library(terra)
library(conflicted)

conflicts_prefer(terra::extract, terra::rotate, dplyr::select, dplyr::filter, purrr::map)


# Get file
# dir(full.names = T, pattern = "cmems_mod_glo")
files <- dir(full.names = T, pattern = "cmems")[c(3)]
# "./cmems_mod_glo_phy_my_0.083deg_P1D-m_uo-vo_110.00E-234.00E_76.92S-29.67S_0.49m_1995-11-15-2001-10-01.nc", "./cmems_mod_glo_phy_my_0.083deg_P1D-m_uo-vo_110.00E-234.00E_76.92S-29.67S_186.13m_1995-11-15-2001-10-01.nc"

# extract the depth string

for (file in files) {
  print(file)

  # get current depth
  depth <- str_extract(file, "\\d+\\.\\d+m")

  # # Load the nc file
  # nc_data <- nc_open(file)

  # # Get var names
  # names(nc_data$var)

  vo <- rast(file, subds = "vo")
  uo <- rast(file, subds = "uo")


  # Function to check and rotate files
  check_and_rotate <- function(var, depth) {
    file_name <- paste0(var, "180_", depth, ".tif")
    if (!file.exists(file_name)) {
      rotated <- terra::rotate(get(var), filename = file_name)
    } else {
      rotated <- rast(file_name)
    }
    return(rotated)
  }

  check_and_resample <- function(var, depth) {
    file_name <- paste0(var, "180_", depth, "_25.tif")
    if (!file.exists(file_name)) {
      resampled <- resample(get(var), template, threads = TRUE, filename = file_name)
    } else {
      resampled <- rast(file_name)
    }
    return(resampled)
  }

  # Rotate or load the rotated files
  vo180 <- check_and_rotate("vo", depth)
  uo180 <- check_and_rotate("uo", depth)

  # Reproject to 25 deg resolution
  # create template spatRaster
  template <- rast(vo180)
  res(template) <- 0.25
  vo25 <- check_and_resample("vo", depth)
  uo25 <- check_and_resample("uo", depth)



  d <- readRDS("Output/tracks_processed_12h.rds")
  st_geometry(d) <- NULL
  d <- d %>% dplyr::filter(trip == 1, SUS == FALSE)
  mq <- cbind(158.95, -54.5)

  d1 <- d %>%
    group_by(id) %>%
    summarise(startdate = first(date), startlon = first(lon), startlat = first(lat), tripdur = first(tripdur))

  out <- NULL

  for (i in 2:2) {
    message(paste(i, "out of", nrow(d1), sep = " "))
    pts <- cbind(d1$startlon[i], d1$startlat[i])
    sdate <- d1$startdate[i] %>% as.Date()
    duration <- d1$tripdur[i]
    id <- d1$id[i]

    # custom particle trace function to accept my own curr_u and curr_v as input
    pt <- tryCatch(
      {
        particle_trace_ll_custom(
          xy = pts,
          time_step = 24 * 3600,
          start_date = sdate,
          end_date = sdate + duration,
          curr_u = uo25,
          curr_v = vo25,
        )
      },
      error = function(e) {
        message(paste("Error in iteration", i, ":", e$message))
        return(NULL)
      }
    )

    if (!is.null(pt)) {
      pt$group <- i
      pt$id <- id
      pt$date_processed <- Sys.Date()
      out <- out %>% bind_rows(pt)
    }
  }

  # save rds
  # saveRDS(out, paste0("output/particle-trace-rk", depth, ".rds"))
  out <- readRDS(paste0("output/particle-trace-", depth, ".rds"))

  # load rds
  old <- readRDS("output/currently_particleTrace.rds") # this was my original anlaysis results using raadtools in currently.
  old <- old %>% dplyr::filter(id %in% out$id)

  length(out$id %>% unique())

  out <- out %>%
    group_by(id) %>%
    mutate(daysSinceDeparture = difftime(date, min(date), units = "days"))

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
    lims(x = c(155, 180)) +
    theme_bw()

  p

  # save plot
  ggsave(paste0("output/particle-trace-", depth, ".png"))
}
