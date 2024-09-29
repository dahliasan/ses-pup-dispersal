library(tidyverse)
library(lubridate)
library(sf)

preprocess_seal_particle_data <- function(seal_data_path = "./Output/tracks_processed_12h.rds", particle_data_path = "./Output/currently_particleTrace.rds") {
    
    # Load data
    cat("Loading seal data\n")
    locw <- readRDS(seal_data_path)

    cat("Loading particle data\n")
    locpt <- readRDS(particle_data_path)

    cat("Preprocessing particle data\n")
    locpt <- locpt %>% dplyr::mutate(date = with_tz(date, tzone = "UTC"))

    cat("Preprocessing seal data\n")
    # Add this line to convert the spatial object to a regular dataframe
    locw <- locw %>% st_drop_geometry()

    # Data preprocessing
    locpt <- locpt %>% dplyr::mutate(date = with_tz(date, tzone = "UTC"))

    # Identify and remove seals with no weanmass data
    all_data_weaners <- readRDS("./Output/all_data_combined.rds")
    seals_with_no_weanmass <- all_data_weaners %>%
        dplyr::group_by(id) %>%
        dplyr::summarise(weanmass = mean(weanmass, na.rm = TRUE)) %>%
        filter(is.na(weanmass)) %>%
        dplyr::pull(id)

    # Prepare locw data
    locw <- locw %>%
        filter(!id %in% seals_with_no_weanmass) %>%
        filter(trip == 1, SUS == FALSE) # only keep first trip

    # Resample data function
    resample_data <- function(df, interval = "1 day") {
        if ("sf" %in% class(df)) {
            cat("Converting spatial object to dataframe\n")
            df <- df %>% st_drop_geometry()
        }
        df %>%
            dplyr::mutate(date = floor_date(date, unit = interval)) %>%
            dplyr::group_by(id, date) %>%
            dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
            dplyr::ungroup()
    }

    # Resample both datasets
    locpt_resampled <- resample_data(locpt, interval = "1 day")
    cat("Resampling seal data\n")
    locw_resampled <- locw %>%
        dplyr::select(id, date, lat, lon) %>%
        resample_data(interval = "1 day")

    # Match date range function
    match_date_range <- function(df1, df2) {
        df1 %>%
            inner_join(df2 %>% select(id, date), by = c("id", "date"))
    }

    locpt_matched <- match_date_range(locpt_resampled, locw_resampled)
    locw_matched <- match_date_range(locw_resampled, locpt_resampled)

    # Calculate days since start
    locpt_matched <- locpt_matched %>%
        dplyr::group_by(id) %>%
        dplyr::mutate(days_since_start = difftime(date, min(date), units = "day"))
    locw_matched <- locw_matched %>%
        dplyr::group_by(id) %>%
        dplyr::mutate(days_since_start = difftime(date, min(date), units = "day"))

    # Ensure both dataframes have the same number of rows per id
    stopifnot(all(table(locpt_matched$id) == table(locw_matched$id)))

    # Sort both dataframes by id and date to ensure alignment
    locpt_matched <- locpt_matched %>% dplyr::arrange(id, date)
    locw_matched <- locw_matched %>% dplyr::arrange(id, date)

    return(list(locpt_matched = locpt_matched, locw_matched = locw_matched))
}
