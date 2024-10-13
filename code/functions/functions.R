require(tidyverse)
require(lubridate)
require(conflicted)

conflicts_prefer(dplyr::summarise, dplyr::filter, dplyr::lag, purrr::map, dplyr::select, .quiet = TRUE)

source("code/functions/bearing_utils.R")
source("code/functions/convert2polarsf.R")
source("code/functions/preprocess_seal_particle_data.R")
source("code/functions/print_and_save_results.R")
source("code/functions/unwrap_lon.R")
load("baseInfo.Rdata")


# Custom Functions ---------------------------------------------------------
calcMeanDepartureDate <- function(dates, type = "median") {
  dates_formatted <- purrr::map(dates, function(x) {
    ifelse(!month(x) %in% c(11, 12), format(x, "2000-%m-%d"), format(x, "1999-%m-%d"))
  }) %>%
    unlist() %>%
    as.Date()

  if (type == "median") {
    return(dates_formatted %>% median())
  }
  if (type == "mean") {
    return(dates_formatted %>% mean())
  }
  if (type == "sd") {
    return(dates_formatted %>% sd())
  }
}

saveDepartureDateResult <- function(results_table, dates, age_group) {
  return(
    bind_rows(
      results_table,
      tibble(
        age_group = age_group,
        n = dates %>% length(),
        median = calcMeanDepartureDate(dates = dates),
        mean = calcMeanDepartureDate(dates, "mean"),
        sd = calcMeanDepartureDate(dates, "sd"),
        median_yday = calcMeanDepartureDate(dates) %>% yday()
      )
    )
  )
}

# survival
# create function to determine survival
get_survival_data <- function() {
  require(tidyverse)
  all_data_weaners <- readRDS("./Output/all_data_combined.rds")
  all_data_weaners %>%
    filter(sim == 0, trip == 1) %>%
    group_by(id) %>%
    summarise(
      seen_6m = if (any(seen_6m == TRUE)) TRUE else FALSE,
      seen_1y = if (any(seen_1y == TRUE)) TRUE else FALSE,
      is_trip_complete = if (all(is_trip_complete == TRUE)) TRUE else FALSE,
      weanmass = first(weanmass),
      birthyear = first(birthyear),
      survive_trip_1 = ifelse(is_trip_complete == TRUE | seen_6m == TRUE, TRUE, FALSE),
      survive_year_1 = ifelse(seen_1y == TRUE, TRUE, FALSE),
    ) %>%
    mutate(
      size = case_when(
        weanmass > 135 ~ "heavy",
        weanmass < 96 ~ "light",
        TRUE ~ "avg"
      ),
      weanmass = ifelse(is.na(weanmass), NA, weanmass)
    )
}

load_all_data_weaners <- function() {
  all_data_weaners <- readRDS("./Output/all_data_combined.rds")
  return(all_data_weaners)
}
