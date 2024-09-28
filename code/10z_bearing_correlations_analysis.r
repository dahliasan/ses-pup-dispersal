{ # Load required libraries
    library(tidyverse)
    library(circular)
    library(geosphere)
    library(lubridate)
    library(conflicted)
    conflicts_prefer(dplyr::filter, dplyr::lag, purrr::map, .quiet = TRUE)

    calculate_bearings <- function(data) {
        data %>%
            mutate(bearing = bearing(cbind(lon, lat), cbind(lead(lon), lead(lat))))
    }

    # Main analysis function
    run_analysis <- function(seal_data, particle_data) {
        # calculate bearing between each location for each id for seal_data and particle_data
        seal_bearings <- calculate_bearings(seal_data)
        particle_bearings <- calculate_bearings(particle_data)

        # Plot bearings over time
        ggplot() +
            geom_point(data = seal_bearings, aes(x = days_since_start, y = bearing), color = "red", size = 0.5) +
            geom_point(data = particle_bearings, aes(x = days_since_start, y = bearing), color = "blue", size = 0.5) +
            theme_minimal() +
            facet_wrap(~id, scales = "free") +
            theme(legend.position = "none") +
            labs(
                title = "Bearings over Time for Seals and Particles",
                x = "Date",
                y = "Bearing (degrees)"
            )

        # plot lat and lon, with bearing vector as an arrow for each point to check if calculation is correct

        ggplot() +
            geom_point(data = seal_bearings, aes(x = lon, y = lat), size = 0.5) +
            geom_segment(
                data = seal_bearings,
                aes(
                    x = lon, y = lat,
                    xend = lead(lon), yend = lead(lat)
                ),
                color = "red", arrow = arrow(length = unit(0.2, "cm"))
            ) +
            theme_minimal() +
            facet_wrap(~id, scales = "free") +
            theme(legend.position = "none") +
            labs(
                title = "Lat and Lon for Seals with Path Arrows",
                x = "Longitude",
                y = "Latitude"
            )


        # Calculate cumulative circular correlation between seal and particle bearings
        cumulative_correlations <- seal_bearings %>%
            left_join(particle_bearings, by = c("id", "days_since_start"), suffix = c("_seal", "_particle")) %>%
            group_by(id) %>%
            arrange(id, days_since_start) %>%
            mutate(
                bearings_used_for_correlation_seal = map(seq_along(bearing_seal), function(i) {
                    paste(bearing_seal[1:i] %>% round(2), collapse = ", ")
                }) %>% unlist(),
                bearings_used_for_correlation_particle = map(seq_along(bearing_particle), function(i) {
                    paste(bearing_particle[1:i] %>% round(2), collapse = ", ")
                }) %>% unlist(),
                cumulative_correlation = map(seq_along(bearing_seal), function(i) {
                    circular::cor.circular(
                        circular::circular(bearing_seal[1:i], units = "degrees", type = "directions", template = "geographics", modulo = "2pi", rotation = "clock"),
                        circular::circular(bearing_particle[1:i], units = "degrees", type = "directions", template = "geographics", modulo = "2pi", rotation = "clock"),
                        test = TRUE
                    )
                })
            ) %>%
            rowwise() %>%
            mutate(
                cor = cumulative_correlation$cor,
                p_value = cumulative_correlation$p.value,
                statistic = cumulative_correlation$statistic
            ) %>%
            select(-cumulative_correlation) %>%
            ungroup()

        cumulative_correlations %>% View()

        # Plot cumulative correlations
        ggplot(cumulative_correlations, aes(x = days_since_start, y = cumulative_correlation)) +
            geom_line() +
            theme_minimal() +
            lims(y = c(-1, 1)) +
            facet_wrap(~id, scales = "free") +
            labs(
                title = "Cumulative Correlation of Bearings between Seals and Particles",
                x = "Days Since Start",
                y = "Cumulative Correlation"
            )

        # Plot histogram
        ggplot(cumulative_correlations, aes(x = cor)) +
            geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
            theme_minimal() +
            labs(
                title = "Distribution of Correlations",
                x = "Correlation", y = "Count"
            )

        # Find critical period
        critical_period <- cumulative_correlations %>%
            # cannot be too early
            filter(days_since_start > 1) %>%
            group_by(id) %>%
            summarize(
                max_correlation = max(cor, na.rm = TRUE),
                days_to_max = days_since_start[which.max(cor)] %>% as.numeric()
            )

        psych::describe(critical_period)
        Hmisc::describe(critical_period)

        # Plot distribution of maximum correlations
        ggplot(critical_period, aes(x = max_correlation)) +
            geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
            theme_minimal() +
            labs(
                title = "Distribution of Maximum Correlations",
                x = "Maximum Correlation", y = "Count"
            )

        # Test if the mean maximum correlation is significantly different from 0
        t_test_result <- t.test(critical_period$max_correlation)
        print(t_test_result)

        #  relationship between maximum correlation and time to max correlation:
        ggplot(critical_period, aes(x = days_to_max, y = max_correlation)) +
            geom_point() +
            geom_smooth(method = "lm") +
            theme_minimal() +
            labs(
                title = "Maximum Correlation vs. Time to Max Correlation",
                x = "Days to Max Correlation", y = "Maximum Correlation"
            )

        correlation_test <- cor.test(critical_period$days_to_max, critical_period$max_correlation)
        print(correlation_test) # negative correlation


        library(mgcv)

        glimpse(cumulative_correlations)

        mod_dat <- cumulative_correlations %>%
            mutate(
                days_since_start = as.numeric(days_since_start),
                id = as.factor(id)
            ) %>%
            filter(!is.na(cor))

        gamm_model <- gamm(cor ~ s(days_since_start, k = 10) + s(id, bs = "re"),
            data = mod_dat,
            method = "REML"
        )
        plot(gamm_model$gam)
        summary(gamm_model$gam)
    }
}

# ------------------------------------------------------------------------------
# Load data
locw <- readRDS("./Output/tracks_processed_12h.rds")
all_data_weaners <- readRDS("./Output/all_data_combined.rds")
load("baseInfo.rdata")
source("code/convert2polarsf.R")
source("code/functions.R")
load("./Data/macca_winter_locs.Rdata") # adult females
locf <- loc %>%
    as_tibble() %>%
    arrange(seal, gmt)
rm(loc)
locpt <- readRDS("./Output/currently_particleTrace.rds") %>%
    # set tz to utc
    mutate(date = with_tz(date, tzone = "UTC"))
# Identify and remove seals with no weanmass data
seals_with_no_weanmass <- all_data_weaners %>%
    group_by(id) %>%
    summarise(weanmass = mean(weanmass, na.rm = TRUE)) %>%
    filter(is.na(weanmass)) %>%
    pull(id)

# Prepare locw data
locw <- locw %>%
    filter(!id %in% seals_with_no_weanmass) %>%
    filter(trip == 1, SUS == FALSE) # only keep first trip

# resample locpt and locw to have same time intervals to make them compare
# we want same number of rows

resample_data <- function(df, interval = "1 day") {
    df %>%
        mutate(date = floor_date(date, unit = interval)) %>%
        group_by(id, date) %>%
        summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
        ungroup()
}

# Resample both datasets
locpt_resampled <- resample_data(locpt, interval = "1 day")
locw_resampled <- locw %>%
    select(id, date, lat, lon) %>%
    # remove geometry
    st_drop_geometry() %>%
    resample_data(interval = "1 day")

# Match date range per id to make sure df have same number of rows
match_date_range <- function(df1, df2) {
    df1 %>%
        inner_join(df2 %>% select(id, date), by = c("id", "date"))
}

locpt_matched <- match_date_range(locpt_resampled, locw_resampled)
locw_matched <- match_date_range(locw_resampled, locpt_resampled)

# calculate days since start
locpt_matched <- locpt_matched %>%
    group_by(id) %>%
    mutate(days_since_start = difftime(date, min(date), units = "day"))
locw_matched <- locw_matched %>%
    group_by(id) %>%
    mutate(days_since_start = difftime(date, min(date), units = "day"))

# plot to check
ggplot() +
    geom_point(data = locpt_matched, aes(x = days_since_start, y = lat), color = "blue", size = 0.5) +
    geom_point(data = locw_matched, aes(x = days_since_start, y = lat), color = "red", size = 0.5) +
    theme_minimal() +
    facet_wrap(~id, scales = "free") +
    # remove legend
    theme(legend.position = "none") +
    labs(
        title = "Latitude over Time for Particles and Seals",
        x = "Days since start",
        y = "Latitude",
        color = "ID"
    )

# Ensure both dataframes have the same number of rows per id
stopifnot(all(table(locpt_matched$id) == table(locw_matched$id)))

# Sort both dataframes by id and date to ensure alignment
locpt_matched <- locpt_matched %>% arrange(id, date)
locw_matched <- locw_matched %>% arrange(id, date)

# Now locpt_resampled and locw_resampled should have comparable time intervals

results <- run_analysis(seal_data = locw_matched, particle_data = locpt_matched)
