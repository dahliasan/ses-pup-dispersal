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

  conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate, dplyr::summarise, dplyr::group_by, dplyr::rename, dplyr::left_join, stats::sd)

  # Set output path
  output_path <- paste0("./output/dispersal_bearing/", Sys.Date(), "/")
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  output_path

  # Load data
  locw <- readRDS("./Output/tracks_processed_12h.rds")
  masterData <- readRDS("./Output/all_data_combined.rds")
  load("baseInfo.rdata")
  source("convert2polarsf.R")
  source("functions.R")
  load("./Data/macca_winter_locs.Rdata") # adult females
  load("10z_generate_summaries_functions.R")
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
  process_weaners <- function(locw, sealsWithNoWeanmass, bearing_ref_day) {
    locw1 <- locw %>%
      group_by(id) %>%
      filter(trip == 1, SUS == FALSE) %>%
      filter(!id %in% sealsWithNoWeanmass) %>%
      mutate(daysSinceStart = difftime(date, min(date), units = "day") %>% round())

    loc_early <- locw1 %>%
      filter(daysSinceStart == bearing_ref_day) %>%
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
    return(list(
      loc_early = loc_early,
      departureDates_results = departureDates_results
    ))
  }

  # Function to process adult females
  process_adult_females <- function(locf, bearing_ref_day) {
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
      filter(daysSinceStart == bearing_ref_day) %>%
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

    return(list(
      locf_early = locf_early,
      departureDates_results = departureDates_results,
      locf_sf = locf_sf
    ))
  }

  # Function to process particle trace
  process_particle_trace <- function(pt, loc_early, bearing_ref_day) {
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
      geom_vline(xintercept = mean(circ_data), color = "black", linetype = 2, size = 1) +
      annotate("label", x = mean(circ_data), y = 20, label = mean(circ_data) %>% round(1)) +
      labs(subtitle = subtitle, y = "") +
      theme_bw()

    ggsave(filename, plot = p, width = 20, height = 15, units = "cm", dpi = 500)
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
    plot_circular_histogram(circ_live, "Survived Trip 1", paste0("1st trip survived, n = ", length(circ_live)), paste0(output_path, "survived_trip1_histogram.png"))
    plot_circular_histogram(circ_dead, "Died Trip 1", paste0("1st trip died, n = ", length(circ_dead)), paste0(output_path, "died_trip1_histogram.png"))

    # Step 4: Analyze survival for the first year
    b <- b %>% mutate(survive = surviveYear1)

    # Create circular data for seals that survived and didn't survive the first year
    circ_live_y1 <- circular(b %>% filter(survive == TRUE) %>% pull(bearing) %% 360, units = "degrees", zero = pi / 2, rotation = "clock")
    circ_dead_y1 <- circular(b %>% filter(survive == FALSE) %>% pull(bearing) %% 360, units = "degrees", zero = pi / 2, rotation = "clock")

    # Plot histograms for 1st year survival
    plot_circular_histogram(circ_live_y1, "Survived Year 1", paste0("1st year survived, n = ", length(circ_live_y1)), paste0(output_path, "survived_year1_histogram.png"))
    plot_circular_histogram(circ_dead_y1, "Died Year 1", paste0("1st year died, n = ", length(circ_dead_y1)), paste0(output_path, "died_year1_histogram.png"))

    # Step 5: Create and save a combined plot of all survival analyses
    combined_plot <- ggarrange(
      plot_circular_histogram(circ_live, "Survived Trip 1", paste0("1st trip survived, n = ", length(circ_live))),
      plot_circular_histogram(circ_dead, "Died Trip 1", paste0("1st trip died, n = ", length(circ_dead))),
      plot_circular_histogram(circ_live_y1, "Survived Year 1", paste0("1st year survived, n = ", length(circ_live_y1))),
      plot_circular_histogram(circ_dead_y1, "Died Year 1", paste0("1st year died, n = ", length(circ_dead_y1))),
      ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), align = "hv"
    )
    ggsave(paste0(output_path, "dispersal_direction_by_survival.png"), plot = combined_plot, width = 20, height = 20, units = "cm", dpi = 500)

    # Interpretation: These circular histograms show the distribution of dispersal directions for different groups.
    # Compare the patterns between survivors and non-survivors to see if there are any preferred directions associated with survival.
    # Look for peaks in the histograms that might indicate favorable directions for survival.

    # Step 6: Perform statistical tests to compare dispersal directions
    tests_survival <- list(
      watson_live_dead_trip1 = watson.two.test(circ_live, circ_dead),
      watson_p_live_dead_trip1 = watson.two.test(circ_p, circ_dead),
      watson_p_female = watson.two.test(circ_p, circ_f)
    )

    HR_tests_survival <- list(
      HR_live = HR_test(circ_live),
      HR_dead = HR_test(circ_dead)
    )

    # Interpretation of statistical tests:
    # 1. Watson's two-sample test (watson_live_dead_trip1):
    #    - Null hypothesis: The two samples (survivors and non-survivors) have the same distribution.
    #    - If p-value < 0.05, reject the null hypothesis, suggesting a significant difference in dispersal directions between survivors and non-survivors.

    # 2. Watson's two-sample test (watson_p_live_dead_trip1):
    #    - Compares the particle trace directions with the directions of seals that didn't survive.
    #    - A significant result (p < 0.05) suggests that the modeled particle directions differ from the actual directions of seals that didn't survive.

    # 3. Watson's two-sample test (watson_p_female):
    #    - Compares the particle trace directions with adult female directions.
    #    - A significant result indicates a difference between modeled dispersal and actual adult female movement patterns.

    # 4. HR test (Hermans-Rasson test):
    #    - Tests for circular uniformity (i.e., no preferred direction) in each group.
    #    - Null hypothesis: The distribution is uniform (no preferred direction).
    #    - If p-value < 0.05, reject the null hypothesis, suggesting a preferred direction in that group.

    # Overall interpretation:
    # - Compare the test results between survivors and non-survivors to see if there are significant differences in dispersal patterns.
    # - Look for any preferred directions in either group that might be associated with survival.
    # - Consider how the modeled particle directions compare to actual seal movements, both for pups and adult females.

    # Return the results of the statistical tests
    return(list(tests_survival = tests_survival, HR_tests_survival = HR_tests_survival))
  }


  run_analysis <- function(bearing_ref_day) {
    # Process Weaners
    weaner_results <- process_weaners(locw, sealsWithNoWeanmass, bearing_ref_day)
    loc_early <- weaner_results$loc_early
    departureDates_results <- weaner_results$departureDates_results

    # Process Adult Females
    adult_females_results <- process_adult_females(locf, bearing_ref_day)
    locf_early <- adult_females_results$locf_early
    departureDates_results <- adult_females_results$departureDates_results

    # Process Particle Trace
    particle_trace_results <- process_particle_trace(pt, loc_early, bearing_ref_day)
    pt_early <- particle_trace_results$pt_early

    # Perform Circular Statistics
    circular_stats <- perform_circular_stats(loc_early, locf_early, pt_early)

    # Generate Summaries
    circ.w <- circular_stats$circ_w
    circ.f <- circular_stats$circ_f
    circ.p <- circular_stats$circ_p
    locf_sf <- adult_females_results$locf_sf
    pt_sf <- particle_trace_results$pt_sf

    generateSummaries(loc_early, locf_early, pt_early, circ.w, circ.f, circ.p, locw, locf_sf, pt_sf, output_path, bearing_ref_day)

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
      bearing_ref_day = bearing_ref_day,
      circular_stats = circular_stats,
      adult_females_results = adult_females_results,
      particle_trace_results = particle_trace_results,
      weaner_results = weaner_results
    ))
  }
}
# Test different reference dates and save results
days <- 5
results <- map(days, run_analysis)

results[[1]]$circular_stats$watson_p_w

results

# Save the results to a CSV file
write.csv(results, "circular_stats_results.csv", row.names = FALSE)

# Save Workspace
save.image(file = "10z_bearing_WORKSPACE.Rdata")
