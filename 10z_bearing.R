# Load necessary libraries
library(tidyverse)
library(sf)
library(lubridate)
library(geosphere)
library(circular)
library(ggpubr)
library(broom)
library(rstatix)
library(effectsize)
library(reprex)
library(conflicted)

conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate, dplyr::summarise, dplyr::group_by, dplyr::rename, dplyr::left_join, stats::sd)

# Load data
locw <- readRDS("./Output/tracks_processed_12h.rds")
masterData <- readRDS("./Output/all_data_combined.rds")
load("baseInfo.rdata")
source("convert2polarsf.R")
source("functions.R")
load("./Data/macca_winter_locs.Rdata") # adult females
locf <- loc %>%
  as_tibble() %>%
  arrange(seal, gmt)
pt <- readRDS("./Output/currently_particleTrace.rds")

# Identify and remove seals with no weanmass data
sealsWithNoWeanmass <- masterData %>%
  group_by(id) %>%
  summarise(weanmass = mean(weanmass, na.rm = TRUE)) %>%
  filter(is.na(weanmass)) %>%
  pull(id)

# Set up table outputs
departureDates_results <- tibble()

# Function to process weaners
process_weaners <- function(locw, sealsWithNoWeanmass, REF_DATE) {
  locw1 <- locw %>%
    group_by(id) %>%
    filter(trip == 1, SUS == FALSE) %>%
    filter(!id %in% sealsWithNoWeanmass) %>%
    mutate(daysSinceStart = difftime(date, min(date), units = "day") %>% round())

  loc_early <- locw1 %>%
    filter(daysSinceStart == REF_DATE) %>%
    group_by(id) %>%
    summarise(
      startdate = min(date),
      enddate = max(date),
      bearing = bearing(mq, c(last(lon), last(lat))),
      type = "argos",
      dist2col = last(dist2col)
    )

  jdates <- loc_early %>%
    group_by(id) %>%
    summarise(startdate = min(startdate)) %>%
    pull(startdate)

  departureDates_results <- saveDepartureDateResult(departureDates_results, jdates, "weaners")
  return(list(loc_early = loc_early, departureDates_results = departureDates_results))
}

# Function to process adult females
process_adult_females <- function(locf, REF_DATE) {
  locf_sf <- convert2polarsf(locf) %>%
    left_join(locf) %>%
    rename(id = seal) %>%
    group_by(id) %>%
    mutate(
      daysSinceStart = difftime(gmt, min(gmt), units = "day") %>% round(),
      dist2col = SGAT::gcDist(mq, cbind(lon, lat))
    )

  bad_females <- locf_sf %>%
    group_by(id) %>%
    summarise(dist2col = first(dist2col)) %>%
    filter(dist2col > 200) %>%
    pull(id)

  locf_sf <- locf_sf %>% filter(!id %in% bad_females)

  locf_early <- locf_sf %>%
    mutate(date = as_date(gmt)) %>%
    group_by(id, date) %>%
    arrange(id, gmt) %>%
    filter(gmt == last(gmt)) %>%
    filter(daysSinceStart == REF_DATE) %>%
    group_by(id) %>%
    summarise(
      startdate = min(gmt),
      enddate = max(gmt),
      bearing = bearing(mq, c(last(lon), last(lat))),
      type = last(type),
      dist2col = last(dist2col)
    )

  jdates <- locf_early %>%
    group_by(id) %>%
    summarise(startdate = min(startdate)) %>%
    pull(startdate)

  departureDates_results <- saveDepartureDateResult(departureDates_results, jdates, "adult females")
  return(list(locf_early = locf_early, departureDates_results = departureDates_results))
}

# Function to process particle trace
process_particle_trace <- function(pt, loc_early, REF_DATE) {
  pt <- pt %>%
    group_by(id) %>%
    filter(!is.na(x)) %>%
    mutate(daysSinceDeparture = difftime(date, min(date), units = "days"))

  pt_sf <- convert2polarsf(pt) %>%
    left_join(pt) %>%
    mutate(dist2col = SGAT::gcDist(mq, cbind(lon, lat)))

  pt_sf <- pt_sf %>%
    filter(id %in% loc_early$id) %>%
    left_join(as.data.frame(loc_early) %>%
      dplyr::select(id, dist2col) %>%
      rename(dist2col.w = dist2col))

  pt_early <- pt_sf %>%
    group_by(id) %>%
    mutate(
      dist2col_diff = abs(dist2col.w - dist2col),
      min_dist2col_diff = min(dist2col_diff)
    ) %>%
    filter(dist2col_diff == min_dist2col_diff) %>%
    mutate(
      bearing = bearing(mq, c(lon, lat)),
      type = NA
    ) %>%
    dplyr::select(id, bearing, type, dist2col, dist2col.w)

  return(pt_early)
}

# Function to plot circular histograms
plot_circular_histogram <- function(circ_data, title, subtitle) {
  ggplot() +
    geom_histogram(
      data = data.frame(circ_data), aes(x = circ_data),
      breaks = seq(0, 360, 45),
      colour = "black",
      fill = "grey"
    ) +
    coord_polar() +
    scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
    geom_vline(xintercept = mean(circ_data), color = "black", linetype = 2, size = 1) +
    annotate("label", x = mean(circ_data), y = 20, label = mean(circ_data) %>% round(1)) +
    labs(subtitle = subtitle, y = "") +
    theme_bw()
}

# Main script
REF_DATE <- 5
weaner_results <- process_weaners(locw, sealsWithNoWeanmass, REF_DATE)
loc_early <- weaner_results$loc_early
departureDates_results <- weaner_results$departureDates_results

adult_females_results <- process_adult_females(locf, REF_DATE)
locf_early <- adult_females_results$locf_early
departureDates_results <- adult_females_results$departureDates_results

pt_early <- process_particle_trace(pt, loc_early, REF_DATE)

# Plot circular histograms
circ_f <- circular(locf_early$bearing %% 360, units = "degrees", zero = pi / 2, rotation = "clock")
circ_w <- circular(loc_early$bearing %% 360, units = "degrees", zero = pi / 2, rotation = "clock")
circ_p <- circular(pt_early$bearing %% 360, units = "degrees", zero = pi / 2, rotation = "clock")

p1 <- plot_circular_histogram(circ_f, "females", paste("adult females (n=", length(circ_f), ")", sep = ""))
p2 <- plot_circular_histogram(circ_w, "weaners", paste("weaners (n=", length(circ_w), ")", sep = ""))
p3 <- plot_circular_histogram(circ_p, "particle trace", paste("particle trace (n=", length(circ_p), ")", sep = ""))

# Save plots
png("/output/dispersal_direction_v2.png", width = 20, height = 15, units = "cm", res = 500)
ggarrange(p2, p1, ncol = 2, labels = c("a", "b"))
dev.off()

# Save workspace
save.image(file = "10z_bearing_WORKSPACE.Rdata")
