Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/terra/proj") # fix Cannot find proj.db warning
# devtools::install_github("mdsumner/currently")
source("code/functions/particle_trace_ll_custom.R")
source("code/functions/functions.R")
# Read in nc file ---------------------------------------------------------
library(ncdf4)
# library(raster)
library(sf)
library(tidyverse)
library(terra)
library(conflicted)
library(parallel)
library(patchwork)

conflicts_prefer(terra::extract, terra::rotate, dplyr::select, dplyr::filter, purrr::map)

get_depths <- function(nc_file) {
    nc_open(file)$dim$depth$vals %>% round()
}

# Get .nc files
files <- dir("data", full.names = T, pattern = ".nc")

# Get weaner locs
d <- readRDS("output/tracks_processed_12h.rds") %>%
    st_drop_geometry() %>%
    filter(trip == 1, SUS == FALSE) %>%
    mutate(group = str_sub(id, 1, 3))

d_list <- d %>%
    group_by(group) %>%
    group_split()

# mq <- cbind(158.95, -54.5)

run_particle_trace <- function(d, uo, vo) {
    d1 <- d %>%
        group_by(id) %>%
        summarise(startdate = first(date), startlon = first(lon), startlat = first(lat), tripdur = first(tripdur))

    traces <- mclapply(
        split(d1, d1$id),
        function(x) {
            particle_trace_ll_custom(
                xy = cbind(x$startlon, x$startlat),
                time_step = 24 * 3600,
                start_date = x$startdate,
                end_date = x$startdate + x$tripdur,
                curr_u = uo,
                curr_v = vo
            )
        },
        mc.cores = detectCores() - 1
    )

    return(traces)
}

process_nc_file <- function(files) {
    uo <- rast(file, subds = "uo")
    vo <- rast(file, subds = "vo")
    return(list(uo = uo, vo = vo))
}

get_cohort_nc_file <- function(cohort, files) {
    files[str_detect(files, cohort)]
}

# Run particle trace for each cohort
out <- d_list %>%
    map(~ {
        cat(paste("Processing cohort:", .x$group[1], "\n"))
        cohort <- .x$group[1]
        file <- get_cohort_nc_file(cohort, files)
        depths <- get_depths(file)
        cat(paste("depths:", depths, "\n"))
        uv <- process_nc_file(file)
        run_particle_trace(.x, uv$uo, uv$vo)
    })

out <- out %>% flatten()

for (i in 1:length(out)) {
    id <- names(out)[i]
    out[[i]] <- out[[i]] %>%
        mutate(id = id) %>%
        select(id, everything())
}

saveRDS(out, "output/weaner_particle_traces_at_200m.rds")


########################################################
out <- readRDS("output/weaner_particle_traces_at_200m.rds")
out <- out %>% bind_rows()
dim(out)
# Check proportions of NA locs (lon and lat)

summarise_missing_locs <- function(data) {
    data %>%
        group_by(id) %>%
        summarise(
            n = n(),
            n_na_loc = sum(is.na(lon) | is.na(lat)),
            prop_na_loc = n_na_loc / n
        )
}


out %>%
    summarise_missing_locs() %>%
    View()

## Compare to previous particle traces
out1 <- readRDS("output/particle-trace-186.13m.rds")
out0 <- readRDS("output/currently_particleTrace.rds")

out1 %>%
    summarise_missing_locs() %>%
    pull(prop_na_loc) %>%
    round(2) %>%
    hist(breaks = seq(0, 1, 0.05))

p1 <- out0 %>%
    summarise_missing_locs() %>%
    ggplot(aes(x = prop_na_loc)) +
    geom_histogram(binwidth = 0.05) +
    labs(title = "Particle trace at 0m") +
    theme_bw()

p2 <- out %>%
    summarise_missing_locs() %>%
    ggplot(aes(x = prop_na_loc)) +
    geom_histogram(binwidth = 0.05) +
    labs(title = "Particle trace at 186m") +
    theme_bw()

p1 + p2 + plot_layout(ncol = 2) + plot_annotation(title = "Proportion of missing locations in particle trace")

ggsave("output/figures/particle_trace_missing_locs.png", width = 8, height = 4)

## Plot traces
out_sf <- out %>%
    drop_na(lon, lat) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    convert2polarsf()

out1_sf <- out1 %>%
    drop_na(lon, lat) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    convert2polarsf()

d_sf <- d %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    convert2polarsf()

ggplot() +
    geom_sf(data = d_sf, size = 0.1, alpha = 0.2, aes(color = "seal")) +
    geom_sf(data = out_sf, size = 0.1, aes(color = "particle")) +
    facet_wrap(~id) +
    labs(title = "Particle trace at 186m") +
    theme_bw()

ggsave("output/figures/particle_trace_at_186m.png", width = 12, height = 8)

out0_sf <- out0 %>%
    drop_na(lon, lat) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    convert2polarsf()

ggplot() +
    geom_sf(data = d_sf, size = 0.1, alpha = 0.2, aes(color = "seal")) +
    geom_sf(data = out0_sf, size = 0.1, aes(color = "particle")) +
    facet_wrap(~id) +
    labs(title = "Particle trace at 0m") +
    theme_bw()

ggsave("output/figures/particle_trace_at_0m.png", width = 12, height = 8)
