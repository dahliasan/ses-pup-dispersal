{
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
  library(CircMLE)

  conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate, dplyr::summarise, dplyr::group_by, dplyr::rename, dplyr::left_join, stats::sd, .quiet = TRUE)

  # Set output path
  output_path <- paste0("./output/dispersal_bearing/", Sys.Date(), "/")
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  output_path

  # Load data
  locw <- readRDS("./Output/tracks_processed_12h.rds")
  masterData <- readRDS("./Output/all_data_combined.rds")
  load("baseInfo.rdata")
  source("code/convert2polarsf.R")
  source("code/functions.R")
  load("./Data/macca_winter_locs.Rdata") # adult females
  source("code/10z_generate_summaries_functions.R")
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

  # Function to process weaners
  process_weaners <- function(locw, sealsWithNoWeanmass, days_after_departure) {
    locw1 <- locw %>%
      group_by(id) %>%
      filter(trip == 1, SUS == FALSE) %>%
      filter(!id %in% sealsWithNoWeanmass) %>%
      mutate(daysSinceStart = difftime(date, min(date), units = "day") %>% round())

    loc_early <- locw1 %>%
      filter(daysSinceStart == days_after_departure) %>%
      group_by(id) %>%
      summarise(
        startdate = min(date),
        enddate = max(date),
        bearing = bearing(mq, c(last(lon), last(lat))),
        type = "argos",
        dist2col = last(dist2col)
      )

    return(list(
      loc_early = loc_early
    ))
  }

  # Function to process adult females
  process_adult_females <- function(locf, days_after_departure) {
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
      filter(daysSinceStart == days_after_departure) %>%
      group_by(id) %>%
      summarise(
        startdate = min(gmt),
        enddate = max(gmt),
        bearing = bearing(mq, c(last(lon), last(lat))),
        type = last(type),
        dist2col = last(dist2col)
      )

    return(list(
      locf_early = locf_early,
      locf_sf = locf_sf
    ))
  }

  # Function to process particle trace
  process_particle_trace <- function(pt, loc_early, days_after_departure) {
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

    return(
      list(
        pt_early = pt_early,
        pt_sf = pt_sf
      )
    )
  }

  # Function to plot circular histograms
  plot_circular_histogram <- function(circ_data, title, subtitle, filename) {
    p <- ggplot() +
      geom_histogram(
        data = data.frame(circ_data), aes(x = circ_data),
        breaks = seq(0, 360, 45),
        colour = "black",
        fill = "grey"
      ) +
      coord_polar() +
      scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
      geom_vline(xintercept = mean(circ_data), color = "black", linetype = 2, linewidth = 1) +
      annotate("label", x = mean(circ_data), y = 20, label = mean(circ_data) %>% round(1)) +
      labs(subtitle = subtitle, y = "") +
      theme_bw()

    ggsave(filename, plot = p, width = 20, height = 15, units = "cm", dpi = 500)

    return(p)
  }

  # Function to perform circular statistics
  perform_circular_stats <- function(loc_early, locf_early, pt_early) {
    circ_f <- circular(locf_early$bearing %% 360, units = "degrees", zero = pi / 2, rotation = "clock")
    circ_w <- circular(loc_early$bearing %% 360, units = "degrees", zero = pi / 2, rotation = "clock")
    circ_p <- circular(pt_early$bearing %% 360, units = "degrees", zero = pi / 2, rotation = "clock")

    # Plot circular histograms
    plot_circular_histogram(circ_f, "females", paste("adult females (n=", length(circ_f), ")", sep = ""), paste0(output_path, "females_histogram.png"))
    plot_circular_histogram(circ_w, "weaners", paste("weaners (n=", length(circ_w), ")", sep = ""), paste0(output_path, "weaners_histogram.png"))
    plot_circular_histogram(circ_p, "particle trace", paste("particle trace (n=", length(circ_p), ")", sep = ""), paste0(output_path, "particle_trace_histogram.png"))

    # perform statistical tests
    tests <- list(
      watson_f = watson.test(circ_f),
      watson_w = watson.test(circ_w),
      watson_f_w = watson.two.test(circ_f, circ_w),
      watson_p_w = watson.two.test(circ_p, circ_w),
      watson_p_f = watson.two.test(circ_p, circ_f),
      HR_f = HR_test(circ_f),
      HR_w = HR_test(circ_w),
      rayleigh_f = rayleigh.test(circ_f),
      rayleigh_w = rayleigh.test(circ_w),
      circ_f = circ_f,
      circ_w = circ_w,
      circ_p = circ_p
    )

    return(tests)
  }

  # Function to generate plots and perform tests based on survival
  plot_survival_analysis <- function(indiv_summ, masterData, loc_early, locf_early, pt_early) {
    # Purpose: This function analyzes and visualizes the relationship between dispersal direction and survival rates for seal pups

    # Step 1: Merge survival data with individual summaries
    d <- readRDS("./Output/all_data_combined.rds")
    ds <- d %>%
      group_by(id) %>%
      filter(trip == 1, sim == 0) %>%
      summarise(
        birthyear = first(birthyear),
        surviveTrip1 = ifelse(first(is_trip_complete) == TRUE | first(seen_6m) == TRUE, TRUE, FALSE),
        surviveYear1 = ifelse(first(seen_1y) == TRUE, TRUE, FALSE),
        weanmass = first(weanmass)
      )

    b <- indiv_summ %>%
      ungroup() %>%
      left_join(ds)

    # Step 2: Categorize seals by weight
    b <- b %>% mutate(
      sizeCat = case_when(
        weanmass > 135 ~ "heavy",
        weanmass < 96 ~ "light",
        TRUE ~ "avg"
      )
    )

    # Step 3: Analyze survival for the first trip
    b <- b %>% mutate(survive = surviveTrip1)

    # Create circular data for seals that survived and didn't survive the first trip
    circ_live <- circular(b %>% filter(survive == TRUE) %>% pull(bearing) %% 360, units = "degrees", zero = pi / 2, rotation = "clock")
    circ_dead <- circular(b %>% filter(survive == FALSE) %>% pull(bearing) %% 360, units = "degrees", zero = pi / 2, rotation = "clock")

    # Plot histograms for 1st trip survival
    p1 <- plot_circular_histogram(circ_live, "Survived Trip 1", paste0("1st trip survived, n = ", length(circ_live)), paste0(output_path, "survived_trip1_histogram.png"))
    p2 <- plot_circular_histogram(circ_dead, "Died Trip 1", paste0("1st trip died, n = ", length(circ_dead)), paste0(output_path, "died_trip1_histogram.png"))

    # Step 4: Analyze survival for the first year
    b <- b %>% mutate(survive = surviveYear1)

    # Create circular data for seals that survived and didn't survive the first year
    circ_live_y1 <- circular(b %>% filter(survive == TRUE) %>% pull(bearing) %% 360, units = "degrees", zero = pi / 2, rotation = "clock")
    circ_dead_y1 <- circular(b %>% filter(survive == FALSE) %>% pull(bearing) %% 360, units = "degrees", zero = pi / 2, rotation = "clock")

    # Plot histograms for 1st year survival
    p3 <- plot_circular_histogram(circ_live_y1, "Survived Year 1", paste0("1st year survived, n = ", length(circ_live_y1)), paste0(output_path, "survived_year1_histogram.png"))
    p4 <- plot_circular_histogram(circ_dead_y1, "Died Year 1", paste0("1st year died, n = ", length(circ_dead_y1)), paste0(output_path, "died_year1_histogram.png"))

    # Step 5: Create and save a combined plot of all survival analyses
    combined_plot <- ggarrange(
      p1,
      p2,
      p3,
      p4,
      ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), align = "hv"
    )
    ggsave(paste0(output_path, "dispersal_direction_by_survival.png"), plot = combined_plot, width = 20, height = 20, units = "cm", dpi = 500)

    # Interpretation: These circular histograms show the distribution of dispersal directions for different groups.
    # Compare the patterns between survivors and non-survivors to see if there are any preferred directions associated with survival.
    # Look for peaks in the histograms that might indicate favorable directions for survival.

    # Step 6: Perform statistical tests to compare dispersal directions
    # Return the results of the statistical tests
    list(
      watson_live_dead_trip1 = watson.two.test(circ_live, circ_dead), watson_p_live_dead_trip1 = watson.two.test(circ_p, circ_dead),
      watson_p_female = watson.two.test(circ_p, circ_f),
      HR_live = HR_test(circ_live),
      HR_dead = HR_test(circ_dead)
    )
  }

  # New function to handle departure dates
  process_departure_dates <- function(locw, locf, sealsWithNoWeanmass, days_after_departure) {
    departureDates_results <- tibble()

    # Process weaners
    locw1 <- locw %>%
      group_by(id) %>%
      filter(trip == 1, SUS == FALSE) %>%
      filter(!id %in% sealsWithNoWeanmass) %>%
      mutate(daysSinceStart = difftime(date, min(date), units = "day") %>% round())

    weaner_dates <- locw1 %>%
      filter(daysSinceStart == days_after_departure) %>%
      group_by(id) %>%
      summarise(startdate = min(date)) %>%
      pull(startdate)

    departureDates_results <- saveDepartureDateResult(departureDates_results, weaner_dates, "weaners")

    # Process adult females
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

    female_dates <- locf_sf %>%
      mutate(date = as_date(gmt)) %>%
      group_by(id, date) %>%
      arrange(id, gmt) %>%
      filter(gmt == last(gmt)) %>%
      filter(daysSinceStart == days_after_departure) %>%
      group_by(id) %>%
      summarise(startdate = min(gmt)) %>%
      pull(startdate)

    departureDates_results <- saveDepartureDateResult(departureDates_results, female_dates, "adult females")

    return(departureDates_results)
  }

  run_analysis <- function(days_after_departure) {
    # Process departure dates
    departureDates_results <- process_departure_dates(locw, locf, sealsWithNoWeanmass, days_after_departure)

    # Process Weaners
    weaner_results <- process_weaners(locw, sealsWithNoWeanmass, days_after_departure)
    loc_early <- weaner_results$loc_early

    # Process Adult Females
    adult_females_results <- process_adult_females(locf, days_after_departure)
    locf_early <- adult_females_results$locf_early

    # Process Particle Trace
    particle_trace_results <- process_particle_trace(pt, loc_early, days_after_departure)
    pt_early <- particle_trace_results$pt_early

    # Perform Circular Statistics
    circular_stats <- perform_circular_stats(loc_early, locf_early, pt_early)

    # Generate Summaries
    circ.w <- circular_stats$circ_w
    circ.f <- circular_stats$circ_f
    circ.p <- circular_stats$circ_p
    locf_sf <- adult_females_results$locf_sf
    pt_sf <- particle_trace_results$pt_sf

    summary_results <- generateSummaries(loc_early, locf_early, pt_early, circ.w, circ.f, circ.p, locw, locf_sf, pt_sf, output_path, days_after_departure)

    indiv_summ <- summary_results$individual_summary

    # Plot Survival Analysis and Perform Tests
    survival_tests <- plot_survival_analysis(
      indiv_summ,
      masterData,
      loc_early,
      locf_early,
      pt_early
    )

    # return results
    return(list(
      days_after_departure = days_after_departure,
      circular_stats = circular_stats,
      adult_females_results = adult_females_results,
      particle_trace_results = particle_trace_results,
      weaner_results = weaner_results,
      departureDates_results = departureDates_results
    ))
  }
}

# Test different reference dates and save results
days <- 5
days_after_departure <- 5
results <- map(days, run_analysis)

results[[1]]$circular_stats$watson_p_w

results

# Save the results to a CSV file
write.csv(results, "circular_stats_results.csv", row.names = FALSE)

# Save Workspace
save.image(file = "10z_bearing_WORKSPACE.Rdata")
