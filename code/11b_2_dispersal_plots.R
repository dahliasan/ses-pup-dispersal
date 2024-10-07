{
    library(tidyverse)
    library(sf)
    library(patchwork)
    library(rnaturalearth)
    library(rnaturalearthdata)

    source("code/functions/functions.R")

    # read in the data all_results
    file <- list.files(output_path, pattern = ".rds", full.names = TRUE)
    dispersal_results <- readRDS(file)

    seal_following_particle <- dispersal_results$analysis_results$seal_following_particle %>%
        select(-contains("watson"))

    seal_dispersal_data <- seal_following_particle %>%
        select(data) %>%
        map_dfr(~.x)

    weaner_first_trip_locs <- seal_dispersal_data %>%
        ungroup() %>%
        left_join(seal_following_particle %>% select(-data), by = "id")

    all_data_weaners <- load_all_data_weaners()

    weaner_locs <- all_data_weaners %>%
        ungroup() %>%
        filter(id %in% seal_following_particle$id, trip == 1, SUS == FALSE) %>%
        select(id, lon, lat) %>%
        left_join(seal_following_particle, by = "id") %>%
        select(-data, -contains("seen"), -is_trip_complete)

    names(weaner_locs)

    weaner_locs_sf <- convert2polarsf(weaner_locs, remove_coords = FALSE)

    # Function to calculate endpoint in stereographic projection
    calculate_endpoint_stere <- function(start_x, start_y, bearing_deg, distance_km) {
        bearing_rad <- bearing_deg * pi / 180
        end_x <- start_x + distance_km * sin(bearing_rad)
        end_y <- start_y + distance_km * cos(bearing_rad)
        return(c(end_x, end_y))
    }

    # Calculate endpoints for each seal
    endpoints <- weaner_locs_sf %>%
        group_by(id, survive_trip_1) %>%
        summarise(
            start_x = first(st_coordinates(.)[, 1]),
            start_y = first(st_coordinates(.)[, 2]),
            seal_mean_bearing = first(seal_mean_bearing),
            end_coords = list(calculate_endpoint_stere(
                first(st_coordinates(.)[, 1]),
                first(st_coordinates(.)[, 2]),
                first(seal_mean_bearing),
                1000
            ))
        ) %>%
        unnest_wider(end_coords, names_sep = "_") %>%
        mutate(
            survive_trip_1 = if_else(survive_trip_1 == TRUE, "survived", "died"),
        ) %>%
        left_join(seal_following_particle %>% select(id, is_following), by = "id")

    particle_endpoints <- weaner_locs_sf %>%
        group_by(id) %>%
        summarise(
            start_x = first(st_coordinates(.)[, 1]),
            start_y = first(st_coordinates(.)[, 2]),
            seal_mean_bearing = first(overall_particle_mean),
            end_coords = list(calculate_endpoint_stere(
                start_x, start_y, first(overall_particle_mean), 1000
            ))
        ) %>%
        unnest_wider(end_coords, names_sep = "_")

    # map
    world <- ne_countries(scale = "medium", returnclass = "sf")

    bbox <- st_bbox(weaner_locs_sf)

    # arramge data by survival first
    weaner_locs_sf <- weaner_locs_sf %>%
        arrange(survive_trip_1, id)

    caption <- paste("Mean particle bearing", seal_following_particle$overall_particle_mean[1] %>% round())

    # seals following
    seals_following <- seal_following_particle %>%
        filter(is_following) %>%
        pull(id)

    seals_not_following <- seal_following_particle %>%
        filter(!is_following) %>%
        pull(id)

    create_dispersal_plot <- function(seal_data, filter_following, title, ncol = 4) {
        ggplot() +
            geom_sf(data = world, fill = "gray80") +
            geom_sf(data = seal_data, size = 0.1, color = "gray", alpha = 0.5) +
            geom_segment(
                data = endpoints %>% filter(id %in% filter_following),
                aes(x = start_x, y = start_y, xend = end_coords_1, yend = end_coords_2, color = survive_trip_1),
                arrow = arrow(length = unit(2, "mm")),
            ) +
            geom_segment(
                data = particle_endpoints %>% filter(id %in% filter_following),
                aes(x = start_x, y = start_y, xend = end_coords_1, yend = end_coords_2),
                arrow = arrow(length = unit(2, "mm")),
                color = "black"
            ) +
            theme_bw() +
            facet_wrap(~id, ncol = ncol) +
            coord_sf(crs = st_crs(seal_data)) +
            labs(x = "Longitude", y = "Latitude", title = title)
    }

    # Create the plot
    p1 <- weaner_locs_sf %>%
        filter(id %in% seals_following) %>%
        create_dispersal_plot(seals_following, "Following", ncol = 5)

    p2 <- weaner_locs_sf %>%
        filter(id %in% seals_not_following) %>%
        create_dispersal_plot(seals_not_following, "Not following", ncol = 3)

    p <- p1 + p2 +
        labs(caption = caption) +
        plot_layout(ncol = 2, guides = "collect", widths = c(1.5, 1)) &
        theme(
            legend.position = "bottom",
            # make facet text small
            strip.text = element_text(size = 10),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),
            # remove grid lines
            panel.grid = element_blank()
        ) &
        lims(x = bbox[c("xmin", "xmax")] + c(-0.1, 0.1), y = bbox[c("ymin", "ymax")] + c(-0.1, 0.1))

    ggsave(file.path(output_path, "dispersal_plot.png"), plot = p, width = 15, height = 10, dpi = 300, bg = "white")
}

{
    ## Plot particle trace
    particle_trace <- readRDS(particle_data_path) %>%
        na.omit()

    particle_trace_sf <- convert2polarsf(particle_trace, remove_coords = FALSE)

    # Group the data by particle ID and convert to lines
    particle_trace_lines <- particle_trace_sf %>%
        group_by(id) %>%
        summarise(do_union = FALSE) %>%
        st_cast("LINESTRING")

    bbox <- st_bbox(particle_trace_sf)

    p <- ggplot() +
        geom_sf(data = particle_trace_lines, size = 0.5, alpha = 0.5, color = "blue") +
        geom_sf(data = world, fill = "gray80") +
        theme_bw() +
        coord_sf(crs = st_crs(particle_trace_sf)) +
        theme(
            panel.grid = element_blank()
        ) +
        lims(x = bbox[c("xmin", "xmax")] + c(-0.1, 0.1), y = bbox[c("ymin", "ymax")] + c(-0.1, 0.1)) +
        labs(x = "Longitude", y = "Latitude", title = "Particle trace")

    ggsave(file.path(output_path, "particle_trace.png"), plot = p, width = 10, height = 10, dpi = 300, bg = "white")
}
