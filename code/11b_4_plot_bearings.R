options(width = 10000)

library(tidyverse)
library(circular)
library(sf)
library(move)
library(conflicted)
library(patchwork)


# Define output folder
output_folder <- "output/11b_4_plot_bearings"

# Create the output folder if it doesn't exist
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Start redirecting console output to a file


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

sink(file.path(output_folder, "console_output.txt"))

cat("Number of unique ids in locf:", locf$id %>% unique() %>% length(), "\n")
cat("Number of unique ids in locw:", locw$id %>% unique() %>% length(), "\n")
cat("Number of unique ids in locps:", locps$id %>% unique() %>% length(), "\n")


## Plot all tracks
all_tracks_plot <- ggplot() +
    geom_sf(data = locf, aes(color = "Female"), alpha = 0.5, size = .2) +
    geom_sf(data = locw, aes(color = "Weaner"), alpha = 0.5, size = .2) +
    geom_sf(data = locps, aes(color = "Particle 0m"), alpha = 0.5, size = .2) +
    theme_bw()

# ggsave(file.path(output_folder, "all_tracks_plot.png"), all_tracks_plot, width = 10, height = 8)

## Create circular plot of all bearings for each data set
plot_bearing_histogram <- function(data, title = NULL) {
    # Calculate mean bearing
    # mean_bearing <- mean.circular(data$bearing, na.rm = TRUE)
    mean_bearing <- data$bearing %>% median.circular(na.rm = TRUE)

    p <- data %>%
        ggplot(aes(x = bearing)) +
        geom_histogram(binwidth = 10, boundary = 0) +
        coord_polar(start = 0, direction = 1) +
        scale_x_continuous(
            breaks = seq(0, 360, by = 45),
            limits = c(0, 360),
            labels = c("N", "", "E", "", "S", "", "W", "", "N")
        ) +
        theme_bw() +
        labs(x = NULL, y = "Count") +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        ) +
        # Add arrow for mean bearing
        annotate("segment",
            x = mean_bearing, y = 0, xend = mean_bearing, yend = Inf,
            arrow = arrow(length = unit(0.5, "cm")), color = "red", linewidth = 1
        ) +
        # Add text label for mean bearing
        annotate("text",
            x = mean_bearing, y = Inf, label = sprintf("%.1f°", as.numeric(mean_bearing)),
            color = "red", vjust = 1.1, hjust = -0.1
        )

    if (!is.null(title)) {
        p <- p + labs(title = title)
    } else {
        p <- p + theme(plot.title = element_blank())
    }

    return(p)
}

# All sequential outbound bearings
p_list <- list(locw, locps, locf)
p_names <- c("Weaner", "Particle 0m", "Female")

# Map of all bearings for each data set
all_bearings_plot <- map2(p_list, p_names, plot_bearing_histogram) %>%
    wrap_plots(ncol = 2)

ggsave(file.path(output_folder, "all_bearings_plot.png"), all_bearings_plot, width = 12, height = 10)

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
p_names <- c("Weaner", "Particle", "Female")

# Map of mean bearings by id
mean_bearings_plot <- map2(p_list, p_names, plot_bearing_histogram) %>%
    wrap_plots(ncol = 2, guides = "collect")

mean_bearings_plot

ggsave(file.path(output_folder, "mean_bearings_plot.png"), mean_bearings_plot, width = 12, height = 10)


## Plot just particle traces
particle_traces_plot <- ggplot() +
    geom_sf(data = locps, aes(color = "Particle 0m"), alpha = 0.5) +
    geom_sf(data = locp, aes(color = "Particle 200m"), alpha = 0.5) +
    theme_bw()

ggsave(file.path(output_folder, "particle_traces_plot.png"), particle_traces_plot, width = 10, height = 8)

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

    cb_palette <- c("#D55E00", "#0072B2")

    ggplot() +
        geom_histogram(data = data1, aes(x = bearing, fill = group1_name), alpha = 1, binwidth = 10, boundary = 0) +
        geom_histogram(data = data2, aes(x = bearing, fill = group2_name), alpha = 0.6, binwidth = 10, boundary = 0) +
        coord_polar(start = 0, direction = 1) +
        scale_x_continuous(
            breaks = seq(0, 360, by = 45),
            limits = c(0, 360),
            labels = c("N", "", "E", "", "S", "", "W", "", "N")
        ) +
        scale_fill_manual(values = setNames(cb_palette, c(group1_name, group2_name))) +
        theme_bw() +
        labs(x = NULL, y = "Count", fill = "Bearings") +
        theme(
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank(),
            legend.position = "bottom"
        ) +
        annotate("segment",
            x = mean_bearing_1, y = 0, xend = mean_bearing_1, yend = Inf,
            arrow = arrow(length = unit(0.5, "cm")), color = cb_palette[1], linewidth = 1
        ) +
        annotate("segment",
            x = mean_bearing_2, y = 0, xend = mean_bearing_2, yend = Inf,
            arrow = arrow(length = unit(0.5, "cm")), color = cb_palette[2], linewidth = 1
        ) +
        annotate("text",
            x = mean_bearing_1, y = Inf,
            label = sprintf("%.1f°", as.numeric(mean_bearing_1)),
            color = cb_palette[1], vjust = 1.1, hjust = -0.1
        ) +
        annotate("text",
            x = mean_bearing_2, y = Inf,
            label = sprintf("%.1f°", as.numeric(mean_bearing_2)),
            color = cb_palette[2], vjust = 1.1, hjust = -0.1
        )
}

## Overall group comparisons
# weaner_female_comparison <- compare_bearings(locw, locf, "Weaner", "Female", "bearing")
# ggsave(file.path(output_folder, "weaner_female_comparison.png"), weaner_female_comparison, width = 10, height = 8)

# weaner_particle_comparison <- compare_bearings(locw, locps, "Weaner", "Particle 0m", "bearing")
# ggsave(file.path(output_folder, "weaner_particle_comparison.png"), weaner_particle_comparison, width = 10, height = 8)

## Compare summarised bearings
summarised_weaner_female_comparison <- compare_bearings(locw_by_id, locf_by_id, "Pups", "Adult Females", "bearing")

summarised_weaner_female_comparison
ggsave(file.path(output_folder, "summarised_weaner_female_comparison.png"), summarised_weaner_female_comparison, width = 10, height = 8)

summarised_weaner_particle_comparison <- compare_bearings(locw_by_id, locps_by_id, "Pups", "Particle", "bearing")

summarised_weaner_particle_comparison
ggsave(file.path(output_folder, "summarised_weaner_particle_comparison.png"), summarised_weaner_particle_comparison, width = 10, height = 8)

summarised_weaner_particle_comparison + summarised_weaner_female_comparison + plot_layout(ncol = 2) + plot_annotation(tag_levels = "a")

ggsave(file.path(output_folder, "combined_mean_bearings_comparison.png"), width = 10, height = 6)


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

### Categorise seals as following current if they are within 45 deg of the mean trace direction
locw_by_id <- locw_by_id %>%
    left_join(locps_by_id %>% st_drop_geometry() %>% select(id, bearing, sd_bearing), by = "id", suffix = c("_seal", "_trace")) %>%
    mutate(
        angle_diff = angle_diff(bearing_seal, bearing_trace),
        is_following = angle_diff < 45
    )

cat("\n\n--- Summary of locw_by_id ---\n")
locw_by_id %>%
    st_drop_geometry() %>%
    mutate(across(where(is.numeric), ~ round(., 1))) %>%
    as.data.frame() %>%
    print()

locw_by_id %>%
    st_drop_geometry() %>%
    psych::describe() %>%
    print()

by_following <- locw_by_id %>%
    group_by(is_following) %>%
    summarise(mean_bearing = mean.circular(bearing_seal, na.rm = TRUE), sd_bearing = sd.circular(bearing_seal, na.rm = TRUE), n = n())


cat("\n\n--- Summary of by_following ---\n")
by_following %>%
    st_drop_geometry() %>%
    mutate(across(where(is.numeric), ~ round(., 1))) %>%
    as.data.frame() %>%
    print()

w_pt <- w_pt %>%
    left_join(locw_by_id %>% st_drop_geometry() %>% select(id, is_following), by = "id")


w_pt %>%
    ggplot(aes(x = days_since_start, y = angle_diff, color = is_following)) +
    geom_line() +
    facet_wrap(~id, scales = "free_x") +
    theme_bw() +
    labs(x = "Days since start", y = "Angle difference (degrees)")

ggsave(paste0(output_folder, "/angle_diff_line_plot.png"), width = 10, height = 8, dpi = 300)


w_pt %>%
    ggplot(aes(x = days_since_start, y = id)) +
    geom_tile(aes(fill = angle_diff)) +
    scale_fill_viridis_c(option = "B") +
    theme_bw() +
    facet_wrap(~is_following, scales = "free") +
    labs(x = "Days since start", y = "Seal ID", fill = "Angle difference")


ggsave(paste0(output_folder, "/angle_diff_tile_plot.png"), width = 10, height = 8, dpi = 300)

library(patchwork)

plot_cummean_bearings <- function(data) {
    data %>%
        ggplot() +
        geom_point(aes(x = cummean_bearing.x, y = days_since_start, color = "Pups"), size = 0.1, alpha = 0.25) +
        geom_point(aes(x = cummean_bearing.y, y = days_since_start, color = "Particle"), size = 0.1, alpha = 0.25) +
        facet_wrap(~id, scales = "free", ncol = 4) +
        coord_polar() +
        scale_x_continuous(limits = c(0, 360)) +
        scale_color_manual(values = c("Pups" = "#D55E00", "Particle" = "#0072B2")) +
        theme_bw() +
        theme(legend.position = "bottom", strip.text = element_text(size = 8)) +
        labs(x = NULL, y = "Days since start", color = "Group")
}

p1 <- w_pt %>%
    filter(is_following) %>%
    plot_cummean_bearings()

w_pt %>%
    filter(is_following) %>%
    distinct(id) %>%
    nrow() -> n1


p2 <- w_pt %>%
    filter(!is_following) %>%
    plot_cummean_bearings()

w_pt %>%
    filter(!is_following) %>%
    distinct(id) %>%
    nrow() -> n2

p1 + p2 + plot_layout(ncol = 2, guides = "collect", widths = c(n1, n2)) + plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom")

ggsave(paste0(output_folder, "/cummean_bearings_plot.png"), width = 16, height = 10, dpi = 300)

w_pt %>%
    ggplot() +
    geom_histogram(aes(x = angle_diff, fill = is_following), alpha = 0.5, binwidth = 10, position = "identity") +
    theme_bw()

ggsave(paste0(output_folder, "/angle_diff_histogram_plot.png"), width = 10, height = 8, dpi = 300)


{
    # Plot sliding mean bearings
    #     library(slider)

    #     calculate_sliding_mean <- function(bearings, before = 2, after = 2) {
    #         slide(bearings, ~ mean.circular(., na.rm = TRUE), .before = before, .after = after, .complete = TRUE) %>%
    #             # convert NULL to NA
    #             map(~ if (is.null(.)) NA else .) %>%
    #             unlist()
    #     }

    #     w_pt$sliding_bearing.x <- calculate_sliding_mean(w_pt$bearing.x, 5, 5)
    #     w_pt$sliding_bearing.y <- calculate_sliding_mean(w_pt$bearing.y, 5, 5)
    #     w_pt$angle_diff.sliding <- angle_diff(w_pt$sliding_bearing.x, w_pt$sliding_bearing.y)

    #     w_pt %>%
    #         ggplot() +
    #         geom_path(aes(x = sliding_bearing.x, y = days_since_start, color = "Weaner"), size = .5, alpha = 0.5) +
    #         geom_path(aes(x = sliding_bearing.y, y = days_since_start, color = "Particle 0m"), size = .5, alpha = 0.5) +
    #         facet_wrap(~id, scales = "free") +
    #         coord_polar() +
    #         scale_x_continuous(limits = c(0, 360)) +
    #         scale_color_manual(values = c("Weaner" = "red", "Particle 0m" = "#6ebeff")) +
    #         theme_bw()

    #     w_pt %>%
    #         ggplot(aes(x = days_since_start, y = angle_diff.sliding)) +
    #         geom_line() +
    #         facet_wrap(~id, scales = "free_x") +
    #         theme_bw() +
    #         labs(x = "Days since start", y = "Angle difference (degrees)")
}


## Survival model analysis of following vs not following
library(tidyverse)
library(lubridate)
library(sf)
library(MuMIn)
library(broom)
library(DHARMa)
all_data_weaners <- load_all_data_weaners() # includes environment data

predictor_vars <- c("sst", "ssha", "eke", "tri", "slope", "SSTgrad", "ice", "dist_to_ice_m", "chl", "chlgrad")

model_data <- all_data_weaners %>%
    filter(id %in% w_pt$id, trip == 1, SUS == FALSE) %>%
    group_by(id) %>%
    summarise(across(any_of(predictor_vars), ~ mean(., na.rm = TRUE))) %>%
    left_join(get_survival_data(), by = "id") %>%
    left_join(
        w_pt %>%
            select(id, is_following) %>%
            st_drop_geometry() %>%
            group_by(id) %>%
            summarise(is_following = last(is_following)),
        by = "id"
    )

cat("\n\nmodel data summary:\n")
print(model_data %>% summary())

print(model_data %>% Hmisc::describe())
print(model_data %>% psych::describe())

# Create contingency tables for is_following vs first trip survival and first year survival

# First trip survival contingency table
trip_contingency <- model_data %>%
    group_by(is_following, survive_trip_1) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = survive_trip_1, values_from = count, values_fill = 0) %>%
    mutate(
        total = `FALSE` + `TRUE`,
        prop_not_survived = `FALSE` / total,
        prop_survived = `TRUE` / total
    )

colnames(trip_contingency)[2:3] <- c("died first trip", "survived first trip")

cat("\n\n--- Contingency table: Is following vs First trip survival ---\n")
print(trip_contingency)

# First year survival contingency table
year_contingency <- model_data %>%
    group_by(is_following, survive_year_1) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = survive_year_1, values_from = count, values_fill = 0) %>%
    mutate(
        total = `FALSE` + `TRUE`,
        prop_not_survived = `FALSE` / total,
        prop_survived = `TRUE` / total
    )

colnames(year_contingency)[2:3] <- c("died first year", "survived first year")

cat("\n\n--- Contingency table: Is following vs First year survival ---\n")
print(year_contingency)

# Calculate odds ratios for survival based on following behavior

# Function to calculate odds ratio and confidence interval
calculate_odds_ratio <- function(contingency_table) {
    odds_ratio <- (contingency_table[2, 2] * contingency_table[1, 1]) /
        (contingency_table[1, 2] * contingency_table[2, 1])

    log_odds_ratio <- log(odds_ratio)
    se_log_odds_ratio <- sqrt(sum(1 / contingency_table))
    ci_lower <- exp(log_odds_ratio - 1.96 * se_log_odds_ratio)
    ci_upper <- exp(log_odds_ratio + 1.96 * se_log_odds_ratio)

    return(list(odds_ratio = odds_ratio, ci_lower = ci_lower, ci_upper = ci_upper))
}

# For first trip survival
trip_cont_matrix <- as.matrix(trip_contingency[, c("died first trip", "survived first trip")])
trip_odds <- calculate_odds_ratio(trip_cont_matrix)

cat("\nOdds of surviving first trip if following vs not following:\n")
cat(sprintf(
    "Odds ratio: %.2f (95%% CI: %.2f - %.2f)\n",
    trip_odds$odds_ratio, trip_odds$ci_lower, trip_odds$ci_upper
))
cat(sprintf(
    "Seals following the current have %.1f times higher odds of surviving their first trip.\n\n",
    trip_odds$odds_ratio
))

# For first year survival
year_cont_matrix <- as.matrix(year_contingency[, c("died first year", "survived first year")])
year_odds <- calculate_odds_ratio(year_cont_matrix)

cat("Odds of surviving first year if following vs not following:\n")
cat(sprintf(
    "Odds ratio: %.2f (95%% CI: %.2f - %.2f)\n",
    year_odds$odds_ratio, year_odds$ci_lower, year_odds$ci_upper
))
cat(sprintf(
    "Seals following the current have %.1f times higher odds of surviving their first year.\n",
    year_odds$odds_ratio
))



# source("code/11b_3.1_get_high_cor_vars.R")
# highly correlated variables: "tri" "dist_to_ice_m" "chl"

model_data <- model_data %>% drop_na(weanmass)

model_trip <- glm(survive_trip_1 ~ is_following * weanmass + birthyear + sst + ssha + eke + slope + SSTgrad + ice + chlgrad,
    data = model_data, family = binomial(link = "logit"), na.action = "na.fail"
)

cat("\n\n--- Summary of model_trip (global model) ---\n")
print(summary(model_trip))

# Add DHARMa model checking for trip survival model
cat("\n\n--- DHARMa model checking for trip survival model ---\n")
dharma_trip <- simulateResiduals(fittedModel = model_trip, n = 1000)
plot(dharma_trip)
testDispersion(dharma_trip)
testZeroInflation(dharma_trip)
testOutliers(dharma_trip)

dredge_trip <- dredge(model_trip)

best_models <- get.models(dredge_trip, subset = delta <= 2)

# Create a function to extract important metrics
extract_model_info <- function(model, model_name) {
    glance_info <- glance(model)
    tibble(
        model_name = model_name,
        formula = paste(format(formula(model)), collapse = " "),
        AICc = AICc(model), # Use AICc instead of AIC
        BIC = glance_info$BIC,
        deviance = glance_info$deviance,
        df.residual = glance_info$df.residual,
        null.deviance = glance_info$null.deviance,
        df.null = glance_info$df.null
    )
}

# Apply the function to each model
model_summary <- map2_dfr(best_models, names(best_models), extract_model_info) %>%
    mutate(delta_AICc = AICc - min(AICc)) %>% # Use AICc instead of AIC
    arrange(delta_AICc)

cat("\n\n--- Model selection summary for trip survival ---\n")
print(model_summary %>% as.data.frame())

avg_model <- model.avg(best_models)
cat("\n\n--- Model averaging summary for trip survival ---\n")
print(summary(avg_model))

importance <- sw(dredge_trip)
cat("\n\n--- Variable importance for trip survival ---\n")
print(importance)

## Year survival
model_year <- glm(survive_year_1 ~ is_following * weanmass + birthyear + sst + ssha + eke + slope + SSTgrad + ice + chlgrad,
    data = model_data, family = binomial(link = "logit"), na.action = "na.fail"
)

cat("\n\n--- Summary of model_year (global model) ---\n")
print(summary(model_year))

# Add DHARMa model checking for year survival model
cat("\n\n--- DHARMa model checking for year survival model ---\n")
dharma_year <- simulateResiduals(fittedModel = model_year, n = 1000)
plot(dharma_year)
testDispersion(dharma_year)
testZeroInflation(dharma_year)
testOutliers(dharma_year)

dredge_year <- dredge(model_year)

best_models_year <- get.models(dredge_year, subset = delta <= 2)

model_summary_year <- map2_dfr(best_models_year, names(best_models_year), extract_model_info) %>%
    mutate(delta_AICc = AICc - min(AICc)) %>%
    arrange(delta_AICc)

cat("\n\n--- Model selection summary for year survival ---\n")
print(model_summary_year %>% as.data.frame())

avg_model_year <- model.avg(best_models_year)
cat("\n\n--- Model averaging summary for year survival ---\n")
print(summary(avg_model_year))

importance_year <- sw(dredge_year)
cat("\n\n--- Variable importance for year survival ---\n")
print(importance_year)

# Generate general summary table


cat("\n\n--- Summary table for all weaners ---\n")

all_data_weaners %>%
    filter(id %in% w_pt$id, trip == 1, SUS == FALSE) %>%
    left_join(get_survival_data()) %>%
    group_by(id, tripdur, last_seen_year, birthdate, blackmass, weanmass, survive_trip_1, survive_year_1) %>%
    summarise(
        max_distance = max(dist2col, na.rm = TRUE) %>% round(0)
    ) %>%
    left_join(
        locw_by_id %>%
            st_drop_geometry(),
        by = "id"
    ) %>%
    ungroup() %>%
    mutate(across(all_of(c("survive_trip_1", "survive_year_1", "is_following")), ~ ifelse(. == TRUE, "yes", "no"))) %>%
    mutate(across(where(is.numeric), ~ round(., 1))) %>%
    mutate(across(any_of(c("tripdur", "max_distance")), ~ round(., 0))) %>%
    mutate(
        bearing_trace = paste0(bearing_trace, " ± ", sd_bearing_trace),
        bearing_seal = paste0(bearing_seal, " ± ", sd_bearing_seal)
    ) %>%
    select(-c(sd_bearing_trace, sd_bearing_seal)) %>%
    rename(
        "Trip duration (days)" = tripdur,
        "Last seen (year)" = last_seen_year,
        "Weaning mass (kg)" = weanmass,
        "Birthdate" = birthdate,
        "Black mass (kg)" = blackmass,
        "Survived trip 1" = survive_trip_1,
        "Survived year 1" = survive_year_1,
        "Max distance to colony (km)" = max_distance,
        "ID" = id,
        "Seal bearing (°)" = bearing_seal,
        "PT bearing (°)" = bearing_trace,
        "Followed PT" = is_following,
        "Δ bearing (°)" = angle_diff
    ) %>%
    as.data.frame() %>%
    print() %>%
    write_csv(paste0(output_folder, "/summary_table_weaners.csv"))



cat("\n\n--- Analysis complete. Output saved. ---\n")

sink()

# Create a new sink for model outputs in Markdown format
sink(file.path(output_folder, "model_outputs.md"))

cat("# Model Analysis Results\n\n")

# First Trip Survival Model
cat("## 1. First Trip Survival Model\n\n")

cat("### 1.1 Global Model Summary\n")
cat("This summary shows the coefficients, standard errors, z-values, and p-values for all variables in the full model.\n\n")
cat("```\n")
print(summary(model_trip))
cat("```\n\n")

cat("### 1.2 DHARMa Model Diagnostics\n")
cat("These diagnostics check the model's assumptions and fit. The plots show residual patterns, while the tests check for specific issues.\n\n")
cat("```\n")
cat("#### Dispersion test (checks if the model is over- or under-dispersed)\n")
print(testDispersion(dharma_trip))
cat("\n#### Zero-inflation test (checks if there are more zeros than expected)\n")
print(testZeroInflation(dharma_trip))
cat("\n#### Outliers test (identifies potential outliers)\n")
print(testOutliers(dharma_trip))
cat("```\n\n")

cat("### 1.3 Model Selection Summary\n")
cat("This table shows the top models based on AICc. Lower AICc values indicate better model fit, considering both goodness of fit and model complexity.\n\n")
cat("```\n")
print(model_summary %>%
    mutate(across(where(is.numeric), ~ signif(., 3))) %>%
    as.data.frame())
cat("```\n\n")

cat("### 1.4 Model Averaging Summary\n")
cat("This summary shows the averaged coefficients across the top models, accounting for model uncertainty.\n\n")
cat("```\n")
print(summary(avg_model))
cat("```\n\n")

cat("### 1.5 Variable Importance\n")
cat("These values indicate the relative importance of each variable across all considered models. Higher values suggest greater importance.\n\n")
cat("```\n")
print(importance)
cat("```\n\n")

# First Year Survival Model
cat("## 2. First Year Survival Model\n\n")

cat("### 2.1 Global Model Summary\n")
cat("This summary shows the coefficients, standard errors, z-values, and p-values for all variables in the full model for first-year survival.\n\n")
cat("```\n")
print(summary(model_year))
cat("```\n\n")

cat("### 2.2 DHARMa Model Diagnostics\n")
cat("These diagnostics check the model's assumptions and fit for the first-year survival model.\n\n")
cat("```\n")
cat("#### Dispersion test\n")
print(testDispersion(dharma_year))
cat("\n#### Zero-inflation test\n")
print(testZeroInflation(dharma_year))
cat("\n#### Outliers test\n")
print(testOutliers(dharma_year))
cat("```\n\n")

cat("### 2.3 Model Selection Summary\n")
cat("This table shows the top models for first-year survival based on AICc.\n\n")
cat("```\n")
print(model_summary_year %>%
    mutate(across(where(is.numeric), ~ signif(., 3))) %>%
    as.data.frame())
cat("```\n\n")

cat("### 2.4 Model Averaging Summary\n")
cat("This summary shows the averaged coefficients across the top models for first-year survival.\n\n")
cat("```\n")
print(summary(avg_model_year))
cat("```\n\n")

cat("### 2.5 Variable Importance\n")
cat("These values indicate the relative importance of each variable across all considered models for first-year survival.\n\n")
cat("```\n")
print(importance_year)
cat("```\n\n")

# Close the model outputs sink
sink()

cat("\n\nModel analysis results have been saved to 'model_outputs.md'\n")

cat("Analysis complete. Output saved to:", output_folder, "\n")
