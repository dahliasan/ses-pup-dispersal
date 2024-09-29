require(lubridate)
require(dplyr)
require(purrr)


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


# dispersal bearing -------------------------------------------------------

angle_diff <- function(theta1, theta2) {
  theta <- abs(theta1 - theta2) %% 360
  return(ifelse(theta > 180, 360 - theta, theta))
}

eastOrWest <- function(theta) {
  theta <- theta %% 360
  return(ifelse(theta >= 0 & theta <= 180, "east", "west"))
}

whichZone <- function(theta) {
  theta <- theta %% 360
  if (theta >= 0 & theta < 45) {
    return("N-NE")
  } else if (theta >= 45 & theta < 90) {
    return("NE-E")
  } else if (theta >= 90 & theta < 135) {
    return("E-SE")
  } else if (theta >= 135 & theta < 180) {
    return("SE-S")
  } else if (theta >= 180 & theta < 225) {
    return("S-SW")
  } else if (theta >= 225 & theta < 270) {
    return("SW-W")
  } else if (theta >= 270 & theta < 315) {
    return("W-NW")
  } else {
    return("NW-N")
  }
}


# survival
# create function to determine survival
get_survival <- function() {
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
      )
    )
}
