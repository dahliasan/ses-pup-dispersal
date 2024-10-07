# Revisions for reviewers comments

rm(list = ls())
library(tidyverse)
library(sf)
# devtools::install_github("ianjonsen/ssmTMB")
# devtools::install_github("SCAR/RAATD/R/duckConfit")
library(duckConfit)


# source("code/functions/functions.R")

# Read animal tracks
locw <- readRDS("./Output/tracks_processed_12h.rds") %>%
  filter(trip == 1) %>%
  mutate(
    group = str_sub(id, 1, 3),
    lon = duckConfit::unwrapLon(lon),
  ) %>%
  st_drop_geometry()

datasets_info <- locw %>%
  group_by(group) %>%
  summarise(
    start_date = min(date) %>% as.Date(),
    end_date = max(date) %>% as.Date(),
    min_lon = min(lon),
    max_lon = max(lon),
    min_lat = min(lat),
    max_lat = max(lat)
  )

# Download subsurface data from copernicus --------------------------------

library(reticulate)


# Create virtual env if it doesn't exist
if (!virtualenv_exists("CopernicusMarine")) {
  virtualenv_create(envname = "CopernicusMarine")
  virtualenv_install("CopernicusMarine", packages = c("copernicusmarine"))
}


# Use try-catch to handle potential errors when activating the virtual environment
tryCatch(
  {
    use_virtualenv("CopernicusMarine", required = TRUE)
  },
  error = function(e) {
    stop("Error activating virtual environment: ", e$message)
  }
)

download_dataset <- function(group_data) {
  tryCatch(
    {
      cmt <- import("copernicusmarine")

      # Set login credentials (uncomment and replace with actual credentials)
      # cmt$login("your_username", "your_password")
      # cmt$login("dfoo1", "soAFO9UK1HllL22%fa%j")

      cmt$subset(
        dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
        variables = list("uo", "vo"),
        minimum_latitude = group_data$min_lat - 1,
        maximum_latitude = group_data$max_lat + 1,
        start_datetime = format(group_data$start_date, "%Y-%m-%dT%H:%M:%S"),
        end_datetime = format(group_data$end_date, "%Y-%m-%dT%H:%M:%S"),
        minimum_depth = 200,
        maximum_depth = 200,
        output_filename = paste0("data/copernicus_", group_data$group, "_", group_data$start_date, "_", group_data$end_date, "_200m.nc")
      )
    },
    error = function(e) {
      warning("Error processing group ", group_data$group, ": ", e$message)
    }
  )
}

datasets_info[2:nrow(datasets_info), ] %>%
  group_split(group) %>%
  purrr::walk(download_dataset)


# import copernicusmarine
#
# copernicusmarine.subset(
#   dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
#   variables=["uo", "vo"],
#   minimum_longitude=110,
#   maximum_longitude=234,
#   minimum_latitude=-76.943,
#   maximum_latitude=-29.64961,
#   start_datetime="1995-11-15T00:00:00",
#   end_datetime="2001-10-01T00:00:00",
#   minimum_depth=200,
#   maximum_depth=200,
# )
