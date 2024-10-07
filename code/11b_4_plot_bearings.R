library(tidyverse)
library(circular)
library(sf)
library(move)
library(conflicted)
library(patchwork)
# install.packages("move")

conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate, dplyr::group_by)

seal_data_path <- "./Output/tracks_processed_12h.rds"
female_data_path <- "data/adult_female_locs.rds"
particle_data_path <- "./Output/particle-trace-186.13m.rds"
surface_particle_data_path <- "./Output/currently_particleTrace.rds"

source("code/functions/functions.R")

colony <- cbind(158.95, -54.5)

calculate_distance_from_start <- function(data, colony = cbind(158.95, -54.5)) {
    data %>%
        group_by(id) %>%
        mutate(
            distance_from_start = geosphere::distGeo(
                cbind(colony[1], colony[2]),
                cbind(lon, lat)
            )
        )
}

identify_outbound_trip <- function(data) {
    data %>%
        calculate_distance_from_start() %>%
        group_by(id) %>%
        arrange(id, date) %>%
        mutate(
            max_distance = max(distance_from_start),
            is_outbound = distance_from_start <= max_distance & row_number() <= which.max(distance_from_start)
        ) %>%
        ungroup()
}

normalize_longitude <- function(lon) {
    ((lon + 180) %% 360) - 180
}

interpolate_track <- function(data) {
    data <- data %>% as.data.frame()

    # Print the number of rows before processing
    cat("Total rows before processing:", nrow(data), "\n")

    processed_data <- data %>%
        group_split(id) %>%
        map_dfr(~ {
            id <- .$id %>% unique()
            cat("Processing ID:", id, "\n")
            cat("Rows for this ID before processing:", nrow(.), "\n")

            # Create move object
            move_obj <- move(
                x = .$lon,
                y = .$lat,
                time = .$date,
                proj = CRS("+proj=longlat +datum=WGS84"),
                data = .,
            )

            # Only interpolate if there are gaps larger than 1 day
            time_diffs <- difftime(.$date, lag(.$date), units = "days") %>% as.numeric()
            print(time_diffs %>% unique())
            if (all(time_diffs %>% na.omit() != 1)) {
                cat("Time diffs not 1 day. Interpolating...\n")
                interpolated <- interpolateTime(move_obj, time = as.difftime(24, units = "hours"), spaceMethod = "greatcircle")
            } else {
                cat("Time diffs are 1 day. Skipping interpolation.\n")
                interpolated <- move_obj
            }

            result <- interpolated %>%
                as_tibble() %>%
                mutate(id = id, date = as.Date(date)) %>%
                dplyr::select(-lon, -lat) %>%
                rename(lon = coords.x1, lat = coords.x2) %>%
                dplyr::select(id, date, lon, lat, everything())

            cat("Rows for this ID after processing:", nrow(result), "\n\n")
            result
        })

    # Print the number of rows after processing
    cat("Total rows after processing:", nrow(processed_data), "\n")

    return(processed_data)
}

preprocess_track <- function(data) {
    data %>%
        mutate(lon = normalize_longitude(lon)) %>%
        # drop_na(lon, lat) %>%
        arrange(id, date) %>%
        group_by(id) %>%
        interpolate_track() %>%
        calculate_bearings() %>%
        mutate(bearing = make_circular(bearing)) %>%
        identify_outbound_trip()
}

# Calculate days since start for each id
calculate_days_since_start <- function(data) {
    data %>%
        group_by(id) %>%
        arrange(id, date) %>%
        mutate(days_since_start = as.numeric(difftime(date, first(date), units = "days"))) %>%
        ungroup()
}

# Read and preprocess data
locw <- readRDS(seal_data_path) %>%
    filter(trip == 1, SUS == FALSE) %>%
    dplyr::select(-daysFromDeployment, -SUS, -trip, -dist2col, -haulout) %>%
    preprocess_track() %>%
    convert2polarsf(remove_coords = FALSE) %>%
    filter(is_outbound)

locp <- readRDS(particle_data_path) %>%
    preprocess_track() %>%
    convert2polarsf(remove_coords = FALSE)

locps <- readRDS(surface_particle_data_path) %>%
    preprocess_track() %>%
    convert2polarsf(remove_coords = FALSE)

locf <- readRDS(female_data_path) %>%
    preprocess_track() %>%
    convert2polarsf(remove_coords = FALSE) %>%
    filter(is_outbound)

## Plot all tracks
ggplot() +
    geom_sf(data = locf, aes(color = "Female"), alpha = 0.5, size = .2) +
    geom_sf(data = locw, aes(color = "Weaner"), alpha = 0.5, size = .2) +
    geom_sf(data = locps, aes(color = "Particle 0m"), alpha = 0.5, size = .2) +
    # geom_sf(data = locp, aes(color = "Particle 200m")) +
    theme_bw()

## Create circular plot of all bearings for each data set
plot_bearing_histogram <- function(data, title = "") {
    # Calculate mean bearing
    # mean_bearing <- mean.circular(data$bearing, na.rm = TRUE)
    mean_bearing <- data$bearing %>% median.circular(na.rm = TRUE)

    data %>%
        ggplot(aes(x = bearing)) +
        geom_histogram(binwidth = 10, boundary = 0) +
        coord_polar(start = 0, direction = 1) +
        scale_x_continuous(
            breaks = seq(0, 360, by = 45),
            limits = c(0, 360),
            labels = c("N", "", "E", "", "S", "", "W", "", "N")
        ) +
        theme_bw() +
        labs(x = NULL, y = "Count", title = title) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        ) +
        # Add arrow for mean bearing
        geom_segment(aes(x = mean_bearing, y = 0, xend = mean_bearing, yend = Inf),
            arrow = arrow(length = unit(0.5, "cm")), color = "red", size = 1
        ) +
        # Add text label for mean bearing
        annotate("text",
            x = mean_bearing, y = Inf, label = sprintf("%.1f°", as.numeric(mean_bearing)),
            color = "red", vjust = 1.1, hjust = -0.1
        )
}

# All sequential outbound bearings
p_list <- list(locw, locps, locf)
p_names <- c("Weaner", "Particle 0m", "Female")

map2(p_list, p_names, plot_bearing_histogram) %>%
    wrap_plots(ncol = 2)

# Mean outbound bearings by id
process_by_id <- function(data) {
    data %>%
        group_by(id) %>%
        summarize(mean_bearing = mean.circular(bearing, na.rm = TRUE), sd_bearing = sd.circular(bearing, na.rm = TRUE)) %>%
        rename(bearing = mean_bearing)
}

locw_by_id <- locw %>% process_by_id()
locps_by_id <- locps %>% process_by_id()
locf_by_id <- locf %>% process_by_id()

p_list <- list(locw_by_id, locps_by_id, locf_by_id)
p_names <- c("Weaner", "Particle 0m", "Female")

map2(p_list, p_names, plot_bearing_histogram) %>%
    wrap_plots(ncol = 2)


## Plot just particle traces
ggplot() +
    geom_sf(data = locps, aes(color = "Particle 0m"), alpha = 0.5) +
    geom_sf(data = locp, aes(color = "Particle 200m"), alpha = 0.5) +
    theme_bw()

## Function to compare bearings between two groups
compare_bearings <- function(data1, data2, group1_name, group2_name, bearing_varname) {
    # Load required library
    library(circular)

    data1 <- data1 %>% rename(bearing = !!bearing_varname)
    data2 <- data2 %>% rename(bearing = !!bearing_varname)

    # Calculate mean bearings
    mean_bearing_1 <- mean.circular(data1$bearing, na.rm = TRUE)
    mean_bearing_2 <- mean.circular(data2$bearing, na.rm = TRUE)

    # Perform Watson-Williams test
    watson_test <- watson.two.test(data1$bearing, data2$bearing)

    # Print results
    cat("Mean bearing for", group1_name, ":", as.numeric(mean_bearing_1), "degrees sd ", sd.circular(data1$bearing, na.rm = TRUE), "\n")
    cat("Mean bearing for", group2_name, ":", as.numeric(mean_bearing_2), "degrees sd ", sd.circular(data2$bearing, na.rm = TRUE), "\n")
    cat("\nWatson-Williams test results:\n")
    print(watson_test)

    # Visualize the difference
    ggplot() +
        geom_histogram(data = data1, aes(x = bearing, fill = group1_name), alpha = 0.5, binwidth = 10, boundary = 0) +
        geom_histogram(data = data2, aes(x = bearing, fill = group2_name), alpha = 0.5, binwidth = 10, boundary = 0) +
        coord_polar(start = 0, direction = 1) +
        scale_x_continuous(
            breaks = seq(0, 360, by = 45),
            limits = c(0, 360),
            labels = c("N", "", "E", "", "S", "", "W", "", "N")
        ) +
        scale_fill_manual(values = c(setNames(c("blue", "red"), c(group1_name, group2_name)))) +
        theme_bw() +
        labs(x = NULL, y = "Count", fill = "Bearings", title = paste("Comparison of", group1_name, "and", group2_name, "Bearings")) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        ) +
        geom_segment(aes(x = mean_bearing_1, y = 0, xend = mean_bearing_1, yend = Inf),
            arrow = arrow(length = unit(0.5, "cm")), color = "blue", size = 1
        ) +
        geom_segment(aes(x = mean_bearing_2, y = 0, xend = mean_bearing_2, yend = Inf),
            arrow = arrow(length = unit(0.5, "cm")), color = "red", size = 1
        ) +
        annotate("text",
            x = mean_bearing_1, y = Inf,
            label = sprintf("%.1f°", as.numeric(mean_bearing_1)),
            color = "blue", vjust = 1.1, hjust = -0.1
        ) +
        annotate("text",
            x = mean_bearing_2, y = Inf,
            label = sprintf("%.1f°", as.numeric(mean_bearing_2)),
            color = "red", vjust = 1.1, hjust = -0.1
        )
}

## Overall group comparisons
compare_bearings(locw, locf, "Weaner", "Female", "bearing")
compare_bearings(locw, locps, "Weaner", "Particle 0m", "bearing")

## Compare summarised bearings
compare_bearings(locw_by_id, locf_by_id, "Weaner", "Female", "bearing")
compare_bearings(locw_by_id, locps_by_id, "Weaner", "Particle 0m", "bearing")

## Individual comparisons

compare_individual_bearings <- function(group1, group2, group1_name, group2_name, bearing_varname = "bearing") {
    # testing purposes
    id <- locw_filtered$id %>%
        unique() %>%
        sample(1)
    group1 <- locw_filtered %>% filter(id == id)
    group2 <- locp_filtered %>% filter(id == id)
    group1_name <- "Weaner"
    group2_name <- "Particle 200m"
    bearing_varname <- "bearing"

    group1 <- group1 %>% rename(bearing = !!bearing_varname)
    group2 <- group2 %>% rename(bearing = !!bearing_varname)

    # Calculate mean bearings
    mean_bearing_1 <- mean.circular(group1$bearing, na.rm = TRUE)
    mean_bearing_2 <- mean.circular(group2$bearing, na.rm = TRUE)

    # Perform Watson's two-sample test
    results <- watson.two.test(group1$bearing, group2$bearing)

    cat("Watson-Williams test results for", group1_name, "and", group2_name, ":\n")
    print(results)

    # plot
    p <- ggplot() +
        geom_histogram(data = group1, aes(x = bearing, fill = group1_name), alpha = 0.5, binwidth = 10, boundary = 0) +
        geom_histogram(data = group2, aes(x = bearing, fill = group2_name), alpha = 0.5, binwidth = 10, boundary = 0) +
        coord_polar(start = 0, direction = 1) +
        scale_x_continuous(
            breaks = seq(0, 360, by = 45),
            limits = c(0, 360),
            labels = c("N", "", "E", "", "S", "", "W", "", "N")
        ) +
        theme_bw() +
        theme(
            axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()
        ) +
        geom_segment(aes(x = mean_bearing_1, y = 0, xend = mean_bearing_1, yend = Inf),
            arrow = arrow(length = unit(0.5, "cm")), color = "blue", size = 1
        ) +
        annotate("text",
            x = mean_bearing_1, y = Inf,
            label = sprintf("%.1f°", as.numeric(mean_bearing_1)),
            color = "blue", vjust = 1.1, hjust = -0.1
        ) +
        annotate("text",
            x = mean_bearing_2, y = Inf,
            label = sprintf("%.1f°", as.numeric(mean_bearing_2)),
            color = "red", vjust = 1.1, hjust = -0.1
        )

    # Combine results into a tibble
    tibble(
        id = group1$id %>% unique(),
        mean_1 = group1$bearing %>% mean.circular(na.rm = TRUE),
        sd_1 = group1$bearing %>% sd.circular(na.rm = TRUE),
        mean_2 = group2$bearing %>% mean.circular(na.rm = TRUE),
        sd_2 = group2$bearing %>% sd.circular(na.rm = TRUE),
        watson_test = list(results),
        plot = list(p)
    )
}

w_ps <- locw_filtered %>%
    group_split(id) %>%
    map(~ compare_individual_bearings(
        .,
        locp_filtered,
        "Weaner",
        "Particle 200m",
        "bearing"
    )) %>%
    bind_rows()


p <- w_ps$plot %>% wrap_plots() + plot_layout(ncol = 4, guides = "collect") & theme(legend.position = "bottom", aspect.ratio = 1)

p
ggsave(p, filename = "./Output/figures/bearings_individual_comparison.png", width = 10, height = 10)

w_ps$plot[1:nrow(w_ps)] %>%
    map(~.x) %>%
    wrap_plots(ncol = 8, guides = "collect") & theme(legend.position = "bottom")

# 3. Calculate mean angular difference
calculate_circular_cummean <- function(bearings) {
    n <- length(bearings)

    sapply(1:n, function(i) {
        bearings_subset <- bearings[1:i]
        mean_bearing <- mean.circular(bearings_subset)
    })
}

w_pt <- locw %>%
    left_join(locps %>% st_drop_geometry() %>% select(id, date, bearing), by = c("id", "date")) %>%
    group_by(id) %>%
    mutate(
        cummean_bearing.x = calculate_circular_cummean(bearing.x),
        cummean_bearing.y = calculate_circular_cummean(bearing.y),
        days_since_start = as.numeric(date - first(date)),
        angle_diff = angle_diff(cummean_bearing.x, cummean_bearing.y),
    )

w_pt %>%
    ggplot() +
    geom_histogram(aes(x = angle_diff), fill = "blue", alpha = 0.5) +
    theme_bw()


w_pt %>%
    ggplot() +
    geom_point(aes(x = cummean_bearing.x, y = days_since_start, color = "Weaner"), size = 0.1, alpha = 0.25) +
    geom_point(aes(x = cummean_bearing.y, y = days_since_start, color = "Particle 0m"), size = 0.1, alpha = 0.25) +
    facet_wrap(~id, scales = "free") +
    coord_polar() +
    scale_x_continuous(limits = c(0, 360)) +
    scale_color_manual(values = c("Weaner" = "red", "Particle 0m" = "#6ebeff")) +
    theme_bw()

w_pt %>%
    ggplot() +
    geom_point(aes(x = days_since_start, y = angle_diff, color = "Particle 0m"), size = 0.2) +
    facet_wrap(~id) +
    theme_bw()


w_pt %>%
    group_by(id) %>%
    summarize(
        min_angle_diff = min(angle_diff, na.rm = TRUE),
        day_at_min_angle_diff = days_since_start[which.min(angle_diff)],
    ) -> w_pt_summary

w_pt_summary %>%
    ggplot() +
    geom_histogram(aes(x = day_at_min_angle_diff), fill = "blue", alpha = 0.5, binwidth = 1) +
    theme_bw()

w_pt_summary %>%
    st_drop_geometry() %>%
    ungroup() %>%
    summarise(mean_day = mean(day_at_min_angle_diff, na.rm = TRUE)) %>%
    pull(mean_day) -> early_period_day

## Filter datasets to the median days since start

locw_ep <- locw %>%
    calculate_days_since_start() %>%
    filter(days_since_start <= early_period_day) %>%
    process_by_id()
locps_ep <- locps %>%
    calculate_days_since_start() %>%
    filter(days_since_start <= early_period_day) %>%
    process_by_id()
locf_ep <- locf %>%
    calculate_days_since_start() %>%
    filter(days_since_start <= early_period_day) %>%
    process_by_id()

## Compare bearings in the early period
compare_bearings(locw_ep, locps_ep, "Weaner", "Particle 0m", "bearing")
compare_bearings(locw_ep, locf_ep, "Weaner", "Female", "bearing") # Filtering by early period shows worst results than just using full trip data that's available for each indvidual.


### Categorise seals as following current if they are within 45 deg of the mean trace direction
locw_by_id <- locw_by_id %>%
    left_join(locps_by_id %>% st_drop_geometry() %>% select(id, bearing), by = "id", suffix = c("_seal", "_trace")) %>%
    mutate(
        angle_diff = angle_diff(bearing_seal, bearing_trace),
        is_following = angle_diff < 45
    )

by_following <- locw_by_id %>%
    group_by(is_following) %>%
    summarise(mean_bearing = mean.circular(bearing_seal, na.rm = TRUE), sd_bearing = sd.circular(bearing_seal, na.rm = TRUE), n = n())
