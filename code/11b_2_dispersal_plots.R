{
    library(tidyverse)
    library(sf)
    library(patchwork)
    library(rnaturalearth)
    library(rnaturalearthdata)

    source("code/functions/functions.R")

    # read in the data all_results
    file <- list.files("output/dispersal_analysis_2", pattern = ".rds", recursive = TRUE, full.names = TRUE) %>% last()
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
        )

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

    # Create the plot
    p <- ggplot() +
        geom_sf(data = world, fill = "gray80") +
        geom_sf(data = weaner_locs_sf, size = 0.1, color = "gray", alpha = 0.5) +
        geom_segment(
            data = endpoints,
            aes(x = start_x, y = start_y, xend = end_coords_1, yend = end_coords_2, color = survive_trip_1),
            arrow = arrow(length = unit(2, "mm")),
        ) +
        geom_segment(
            data = particle_endpoints,
            aes(x = start_x, y = start_y, xend = end_coords_1, yend = end_coords_2),
            arrow = arrow(length = unit(2, "mm")),
            color = "black"
        ) +
        theme_bw() +
        facet_wrap(~id) +
        coord_sf(crs = st_crs(weaner_locs_sf)) +
        theme(
            # make facet text small
            strip.text = element_text(size = 8),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 8),
            # remove grid lines
            panel.grid = element_blank()
        ) +
        lims(x = bbox[c("xmin", "xmax")] + c(-0.5, 0.5), y = bbox[c("ymin", "ymax")] + c(-0.5, 0.5)) +
        labs(x = "Longitude", y = "Latitude")

    ggsave(file.path(output_path, "dispersal_plot.png"), plot = p, width = 12, height = 10, dpi = 300, bg = "white")
}
