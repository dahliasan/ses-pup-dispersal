---
title: Dispersal Analysis
layout: default
---

---
title: "Southern Elephant Seal Pup Dispersal Analysis"
author: "Dahlia Foo"
date: today
format:
    html:
        embed-resources: true
        code-fold: true
        code-line-numbers: true
---




``` r
library(tidyverse)
library(circular)
```

```
## 
## Attaching package: 'circular'
```

```
## The following objects are masked from 'package:stats':
## 
##     sd, var
```

``` r
library(sf)
library(move)
```

```
## Warning: package 'move' was built under R version 4.4.1
```

```
## Loading required package: geosphere
```

```
## Warning: package 'geosphere' was built under R version 4.4.1
```

```
## Loading required package: sp
```

```
## Loading required package: raster
```

```
## 
## Attaching package: 'raster'
```

```
## The following object is masked from 'package:dplyr':
## 
##     select
```

```
## A successor to `move` has been developed and is available on cran (`move2`, https://bartk.gitlab.io/move2/). This brings speed improvements and is based on `sf`. Feedback (including missing functionality) on this new package is welcome, for new projects it might be worth considering starting directly with `move2`.
```

``` r
library(conflicted)
library(patchwork)
library(skimr)
library(gtsummary)
```

```
## Warning: package 'gtsummary' was built under R version 4.4.1
```

``` r
library(gt)
```

```
## Warning: package 'gt' was built under R version 4.4.1
```

``` r
library(gtExtras)
library(summarytools)
library(GGally)
```

```
## Registered S3 method overwritten by 'GGally':
##   method from   
##   +.gg   ggplot2
```

``` r
library(broom)
library(DHARMa)
```

```
## This is DHARMa 0.4.6. For overview type '?DHARMa'. For recent changes, type news(package = 'DHARMa')
```

``` r
library(MuMIn)
library(here)

conflicts_prefer(dplyr::filter, dplyr::select, dplyr::mutate, dplyr::group_by)
```

```
## [conflicted] Will prefer dplyr::filter over any other package.
```

```
## [conflicted] Will prefer dplyr::select over any other package.
## [conflicted] Will prefer dplyr::mutate over any other package.
## [conflicted] Will prefer dplyr::group_by over any other package.
```

``` r
theme_set(theme_bw())

# Define output folder
output_folder <- here("output", "dispersal_analysis_3")

# Create the output folder if it doesn't exist
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
```

## Data Loading and Preprocessing


``` r
seal_data_path <- here("output", "tracks_processed_12h.rds")
female_data_path <- here("data", "adult_female_locs.rds")
particle_data_path <- here("output", "particle-trace-186.13m.rds")
surface_particle_data_path <- here("output", "currently_particleTrace.rds")

source(here("code", "functions", "functions.R"))

colony <- cbind(158.95, -54.5)

# Include your function definitions

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
            cat("Processing ", id, "...\n")
            initial_rows <- nrow(.)
            do_interpolate <- FALSE

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
            if (all(time_diffs %>% na.omit() != 1)) {
                cat("time interval not 1 day. interpolating...\n")
                do_interpolate <- TRUE
                interpolated <- interpolateTime(move_obj, time = as.difftime(24, units = "hours"), spaceMethod = "greatcircle")
            } else {
                interpolated <- move_obj
            }

            result <- interpolated %>%
                as_tibble() %>%
                mutate(id = id, date = as.Date(date)) %>%
                dplyr::select(-lon, -lat) %>%
                rename(lon = coords.x1, lat = coords.x2) %>%
                dplyr::select(id, date, lon, lat, everything())

            if (do_interpolate) {
                cat("interpolated from", initial_rows, "to", nrow(result), "rows\n")
            }
            return(result)
        })

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
```


``` r
locw <- readRDS(seal_data_path) %>%
    filter(trip == 1, SUS == FALSE) %>%
    dplyr::select(-daysFromDeployment, -SUS, -trip, -dist2col, -haulout) %>%
    preprocess_track() %>%
    convert2polarsf(remove_coords = FALSE) %>%
    filter(is_outbound)
```

```
## Total rows before processing: 10690 
## Processing  mq1-17215-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 23 to 12 rows
## Processing  mq1-17217-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 242 to 121 rows
## Processing  mq1-17219-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 65 to 33 rows
## Processing  mq1-20916-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 21 to 11 rows
## Processing  mq1-20918-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 270 to 135 rows
## Processing  mq1-22483-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 242 to 121 rows
## Processing  mq1-22484-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 45 to 23 rows
## Processing  mq1-22486-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 231 to 116 rows
## Processing  mq1-22490-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 276 to 138 rows
## Processing  mq1-22499-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 351 to 176 rows
## Processing  mq1-22500-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 291 to 146 rows
## Processing  mq1-22501-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 74 to 37 rows
## Processing  mq1-26625-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 217 to 109 rows
## Processing  mq1-26627-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 128 to 64 rows
## Processing  mq1-26628-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 237 to 119 rows
## Processing  mq1-26629-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 21 to 11 rows
## Processing  mq1-26633-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 20 to 10 rows
## Processing  mq1-26635-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 288 to 144 rows
## Processing  mq1-2849-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 259 to 130 rows
## Processing  mq1-5811-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 70 to 35 rows
## Processing  mq1-5814-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 294 to 147 rows
## Processing  mq2-20916-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 22 to 11 rows
## Processing  mq2-20917-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 10 to 5 rows
## Processing  mq2-22500-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 13 to 7 rows
## Processing  mq2-26623-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 349 to 175 rows
## Processing  mq2-28482-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 91 to 46 rows
## Processing  mq3-17217-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 353 to 177 rows
## Processing  mq3-20918-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 300 to 150 rows
## Processing  mq3-22484-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 312 to 156 rows
## Processing  mq3-22488-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 483 to 242 rows
## Processing  mq3-22498-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 93 to 47 rows
## Processing  mq3-26627-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 355 to 178 rows
## Processing  mq3-26629-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 219 to 110 rows
## Processing  mq3-2849-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 366 to 183 rows
## Processing  mq3-28494-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 54 to 27 rows
## Processing  mq3-28496-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 235 to 118 rows
## Processing  mq3-28497-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 272 to 136 rows
## Processing  mq3-28500-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 19 to 10 rows
## Processing  mq3-28504-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 350 to 175 rows
## Processing  mq3-5812-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 335 to 168 rows
## Processing  mq4-20918-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 293 to 147 rows
## Processing  mq4-Alice-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 410 to 205 rows
## Processing  mq4-Billie-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 401 to 201 rows
## Processing  mq4-Cleo-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 462 to 231 rows
## Processing  mq4-Doris-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 555 to 278 rows
## Processing  mq4-Ella-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 186 to 93 rows
## Processing  mq4-FirstOne-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 135 to 68 rows
## Processing  mq4-Flora-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 352 to 176 rows
```

``` r
locf <- readRDS(female_data_path) %>%
    preprocess_track() %>%
    convert2polarsf(remove_coords = FALSE) %>%
    filter(is_outbound)
```

```
## Total rows before processing: 101905 
## Processing  1 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 926 to 232 rows
## Processing  2 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 934 to 234 rows
## Processing  3 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 926 to 232 rows
## Processing  4 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 988 to 247 rows
## Processing  7 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 914 to 229 rows
## Processing  10 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 950 to 238 rows
## Processing  11 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 1044 to 261 rows
## Processing  13 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 871 to 218 rows
## Processing  14 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 909 to 228 rows
## Processing  17 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 617 to 155 rows
## Processing  18 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 948 to 237 rows
## Processing  22 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 955 to 239 rows
## Processing  25 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 912 to 228 rows
## Processing  26 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 754 to 189 rows
## Processing  27 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 935 to 234 rows
## Processing  31 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 936 to 234 rows
## Processing  34 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 932 to 233 rows
## Processing  35 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 916 to 229 rows
## Processing  36 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 939 to 235 rows
## Processing  39 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 931 to 233 rows
## Processing  42 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 929 to 233 rows
## Processing  43 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 939 to 235 rows
## Processing  45 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 991 to 248 rows
## Processing  46 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 991 to 248 rows
## Processing  47 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 991 to 248 rows
## Processing  48 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 916 to 229 rows
## Processing  49 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 981 to 246 rows
## Processing  50 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 936 to 234 rows
## Processing  53 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 482 to 121 rows
## Processing  55 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 953 to 239 rows
## Processing  57 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 930 to 233 rows
## Processing  58 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 944 to 236 rows
## Processing  68 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 964 to 241 rows
## Processing  70 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 1846 to 154 rows
## Processing  71 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 1875 to 157 rows
## Processing  72 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 1754 to 147 rows
## Processing  73 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2161 to 181 rows
## Processing  74 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 986 to 83 rows
## Processing  75 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 842 to 71 rows
## Processing  76 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2629 to 220 rows
## Processing  77 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2763 to 231 rows
## Processing  78 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 308 to 26 rows
## Processing  79 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2534 to 212 rows
## Processing  80 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2615 to 218 rows
## Processing  81 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2824 to 236 rows
## Processing  82 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2625 to 219 rows
## Processing  83 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 650 to 55 rows
## Processing  85 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2229 to 186 rows
## Processing  86 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 913 to 229 rows
## Processing  87 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 966 to 242 rows
## Processing  89 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 939 to 235 rows
## Processing  90 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 937 to 235 rows
## Processing  94 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2795 to 233 rows
## Processing  95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2833 to 237 rows
## Processing  96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2746 to 229 rows
## Processing  97 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2203 to 184 rows
## Processing  98 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 3190 to 266 rows
## Processing  99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2877 to 240 rows
## Processing  100 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2091 to 175 rows
## Processing  101 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2765 to 231 rows
## Processing  102 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2763 to 231 rows
## Processing  103 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 839 to 70 rows
## Processing  104 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2746 to 229 rows
## Processing  105 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2806 to 234 rows
## Processing  106 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2772 to 232 rows
## Processing  107 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2907 to 243 rows
## Processing  108 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 2992 to 250 rows
```

``` r
locp <- readRDS(particle_data_path) %>%
    preprocess_track() %>%
    convert2polarsf(remove_coords = FALSE)
```

```
## Total rows before processing: 5603 
## Processing  mq1-17215-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-17217-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-17219-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-20916-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-20918-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-22483-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-22484-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-22486-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-22490-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-22499-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-22500-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-22501-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-26625-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-26627-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-26628-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-26629-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-26633-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-26635-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-2849-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-5811-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq1-5814-95 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq2-20916-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq2-20917-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq2-22500-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq2-26623-96 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq2-28482-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-17217-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-20918-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-22484-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-22488-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-22498-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-26627-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-26629-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-2849-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-28494-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-28496-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-28497-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-28500-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-28504-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq3-5812-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq4-20918-00 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq4-Alice-00 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq4-Billie-00 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq4-Cleo-00 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq4-Doris-00 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq4-Ella-00 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq4-FirstOne-00 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## Processing  mq4-Flora-00 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

``` r
locps <- readRDS(surface_particle_data_path) %>%
    preprocess_track() %>%
    convert2polarsf(remove_coords = FALSE)
```

```
## Total rows before processing: 11164 
## Processing  mq1-17215-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 23 to 12 rows
## Processing  mq1-17217-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 242 to 121 rows
## Processing  mq1-17219-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 65 to 33 rows
## Processing  mq1-20916-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 21 to 11 rows
## Processing  mq1-20918-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 270 to 135 rows
## Processing  mq1-22483-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 242 to 121 rows
## Processing  mq1-22484-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 45 to 23 rows
## Processing  mq1-22486-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 231 to 116 rows
## Processing  mq1-22490-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 276 to 138 rows
## Processing  mq1-22499-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 351 to 176 rows
## Processing  mq1-22500-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 291 to 146 rows
## Processing  mq1-22501-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 74 to 37 rows
## Processing  mq1-26625-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 217 to 109 rows
## Processing  mq1-26627-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 128 to 64 rows
## Processing  mq1-26628-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 237 to 119 rows
## Processing  mq1-26629-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 21 to 11 rows
## Processing  mq1-26633-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 20 to 10 rows
## Processing  mq1-26635-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 288 to 144 rows
## Processing  mq1-2849-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 259 to 130 rows
## Processing  mq1-5811-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 70 to 35 rows
## Processing  mq1-5814-95 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 294 to 147 rows
## Processing  mq2-20916-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 22 to 11 rows
## Processing  mq2-20917-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 10 to 5 rows
## Processing  mq2-22500-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 13 to 7 rows
## Processing  mq2-26623-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 349 to 175 rows
## Processing  mq2-28482-96 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 91 to 46 rows
## Processing  mq3-17217-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 353 to 177 rows
## Processing  mq3-20918-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 300 to 150 rows
## Processing  mq3-22484-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 312 to 156 rows
## Processing  mq3-22488-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 483 to 242 rows
## Processing  mq3-22498-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 93 to 47 rows
## Processing  mq3-26627-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 355 to 178 rows
## Processing  mq3-26629-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 219 to 110 rows
## Processing  mq3-2849-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 366 to 183 rows
## Processing  mq3-28494-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 82 to 41 rows
## Processing  mq3-28496-99 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 681 to 234 rows
## Processing  mq3-28497-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 272 to 136 rows
## Processing  mq3-28500-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 19 to 10 rows
## Processing  mq3-28504-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 350 to 175 rows
## Processing  mq3-5812-99 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 335 to 168 rows
## Processing  mq4-20918-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 293 to 147 rows
## Processing  mq4-Alice-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 410 to 205 rows
## Processing  mq4-Billie-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 401 to 201 rows
## Processing  mq4-Cleo-00 ...
```

```
## Warning in .local(x, y, time, data, proj, ...): There were NA locations detected and omitted. Currently they are not stored in unusedrecords
```

```
## Warning: Unknown or uninitialised column: `citation`.
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 462 to 213 rows
## Processing  mq4-Doris-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 555 to 278 rows
## Processing  mq4-Ella-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 186 to 93 rows
## Processing  mq4-FirstOne-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 135 to 68 rows
## Processing  mq4-Flora-00 ...
```

```
## Warning: Unknown or uninitialised column: `citation`.
## Setting row names on a tibble is deprecated.
```

```
## time interval not 1 day. interpolating...
## interpolated from 352 to 176 rows
```

## Data Summary


```
## Number of unique ids in locf: 67
```

```
## Number of unique ids in locw: 48
```

```
## Number of unique ids in locps: 48
```

```
## Number of unique ids in locp: 48
```

## Plotting All Tracks


``` r
# Option 1: Colorblind-friendly palette with muted tones
palette1 <- c(
    "Female" = "#999999", # Muted grey
    "Weaner" = "#E69F00", # Orange
    "Particle 0m" = "#56B4E9", # Light blue
    "Particle 186m" = "#0072B2" # Dark blue
)

# Option 2: Palette with richer colors
palette2 <- c(
    "Female" = "#7F7F7F", # Warm grey
    "Weaner" = "#D55E00", # Rich orange-red
    "Particle 0m" = "#69B3E7", # Sky blue
    "Particle 186m" = "#1B4B8C" # Deep navy
)

# Option 3: Palette with more contrast
palette3 <- c(
    "Female" = "#8C8C8C", # Medium grey
    "Weaner" = "#CC3311", # Bright red
    "Particle 0m" = "#33B2FF", # Bright blue
    "Particle 186m" = "#004C99" # Dark blue
)


# Get the bounding box of all tracks combined
bb <- sf::st_bbox(
    sf::st_union(
        sf::st_union(locf),
        sf::st_union(locw),
        sf::st_union(locps),
        sf::st_union(locp)
    )
)

# Add nudge values for the labels
orsi_sf <- orsi_sf %>%
    mutate(
        nudge_x = c(-1000, -4800, -4500, -4000),
        nudge_y = c(0, 5300, 4800, 3000),
        front = toupper(front)
    )

# Convert Macquarie Island coordinates to sf object
mq_sf <- st_as_sf(data.frame(lon = 158.95, lat = -54.5),
    coords = c("lon", "lat"),
    crs = 4326
) # assuming WGS84


all_tracks_plot <- ggplot() +
    geom_sf(data = orsi_sf, linetype = "dashed", colour = "grey30") +
    geom_sf_text(
        data = orsi_sf,
        aes(label = front),
        nudge_x = orsi_sf$nudge_x,
        nudge_y = orsi_sf$nudge_y
    ) +
    geom_sf(data = locf, aes(color = "Female"), alpha = 0.4, size = .5) +
    geom_sf(data = locw, aes(color = "Weaner"), alpha = 0.4, size = .5) +
    geom_sf(data = locps, aes(color = "Particle 0m"), alpha = 0.4, size = .5) +
    geom_sf(data = locp, aes(color = "Particle 186m"), alpha = 0.4, size = .5) +
    geom_sf(data = mq_sf, shape = 17, size = 3, color = "black") +
    theme_bw() +
    # Set the x and y limits using the bounding box
    coord_sf(
        xlim = c(bb["xmin"], bb["xmax"]),
        ylim = c(bb["ymin"], bb["ymax"] + 1000)
    ) +
    scale_color_manual(values = palette3) +
    labs(color = "Group", x = "Longitude", y = "Latitude")

all_tracks_plot
```

![All tracks plot with ORSI fronts](figure/plot-all-tracks-1.png)

``` r
ggsave(file.path(output_folder, "all_tracks_plot.png"), all_tracks_plot, width = 8, height = 4)
```


## Bearing Analysis


``` r
## Create circular plot of all bearings for each data set
plot_bearing_histogram <- function(data, title = NULL) {
    # Calculate mean bearing
    mean_bearing <- mean.circular(data$bearing, na.rm = TRUE)
    # mean_bearing <- data$bearing %>% median.circular(na.rm = TRUE)

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
p_list <- list(locw, locps, locf, locp)
p_names <- c("Weaner", "Particle 0m", "Female", "Particle 186m")

# Map of all bearings for each data set
all_bearings_plot <- map2(p_list, p_names, plot_bearing_histogram) %>%
    wrap_plots(ncol = 2)

all_bearings_plot
```

![All sequential outbound bearings](figure/bearing-analysis-1.png)

``` r
ggsave(file.path(output_folder, "all_bearings_plot.png"), all_bearings_plot, width = 12, height = 10)
```

## Create mean bearings by id


``` r
# Calculate mean bearings by id
process_by_id <- function(data) {
    data %>%
        group_by(id) %>%
        summarize(mean_bearing = mean.circular(bearing, na.rm = TRUE), sd_bearing = sd.circular(bearing, na.rm = TRUE)) %>%
        rename(bearing = mean_bearing)
}


locw_by_id <- locw %>% process_by_id()
locps_by_id <- locps %>% process_by_id()
locf_by_id <- locf %>% process_by_id()
locp_by_id <- locp %>% process_by_id()

p_list <- list(locw_by_id, locps_by_id, locf_by_id, locp_by_id)
p_names <- c("Weaner", "Particle 0m", "Female", "Particle 186m")

# Map of mean bearings by id
mean_bearings_plot <- map2(p_list, p_names, plot_bearing_histogram) %>%
    wrap_plots(ncol = 2, guides = "collect")

mean_bearings_plot
```

![Mean bearings by id](figure/mean-bearings-by-id-1.png)

``` r
ggsave(file.path(output_folder, "mean_bearings_plot.png"), mean_bearings_plot, width = 12, height = 10)
```

## Comparing Bearings


``` r
compare_bearings <- function(data1, data2, group1_name, group2_name, bearing_varname) {
    data1 <- data1 %>% rename(bearing = !!bearing_varname)
    data2 <- data2 %>% rename(bearing = !!bearing_varname)

    # Calculate mean bearings
    mean_bearing_1 <- mean.circular(data1$bearing, na.rm = TRUE)
    mean_bearing_2 <- mean.circular(data2$bearing, na.rm = TRUE)

    # Perform Watson-Williams test
    watson_test <- watson.two.test(data1$bearing, data2$bearing)

    # Create tidy table for results
    results_table <- tibble(
        Group = c(group1_name, group2_name),
        Mean_Bearing = c(as.numeric(mean_bearing_1), as.numeric(mean_bearing_2)),
        SD_Bearing = c(sd.circular(data1$bearing, na.rm = TRUE), sd.circular(data2$bearing, na.rm = TRUE))
    )

    # Print results
    print(results_table)
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
            color = cb_palette[1], vjust = 1.1, hjust = -0.1, size = 3
        ) +
        annotate("text",
            x = mean_bearing_2, y = Inf,
            label = sprintf("%.1f°", as.numeric(mean_bearing_2)),
            color = cb_palette[2], vjust = 1.1, hjust = -0.1, size = 3
        )
}

## Compare individual mean bearings
p1 <- compare_bearings(locw_by_id, locf_by_id, "Pups", "Adult Females", "bearing")
```

```
## # A tibble: 2 × 3
##   Group         Mean_Bearing SD_Bearing
##   <chr>                <dbl>      <dbl>
## 1 Pups                  120.      0.969
## 2 Adult Females         138.      1.23 
## 
##       Watson's Two-Sample Test of Homogeneity 
## 
## Test Statistic: 0.221 
## 0.01 < P-value < 0.05 
## 
```

``` r
p1
```

![plot of chunk compare-bearings](figure/compare-bearings-1.png)

``` r
p2 <- compare_bearings(locw_by_id, locps_by_id, "Pups", "Particle 0m", "bearing")
```

```
## # A tibble: 2 × 3
##   Group       Mean_Bearing SD_Bearing
##   <chr>              <dbl>      <dbl>
## 1 Pups                120.      0.969
## 2 Particle 0m         115.      0.629
## 
##       Watson's Two-Sample Test of Homogeneity 
## 
## Test Statistic: 0.1768 
## 0.05 < P-value < 0.10 
## 
```

``` r
p2
```

![plot of chunk compare-bearings](figure/compare-bearings-2.png)

``` r
p3 <- compare_bearings(locw_by_id, locp_by_id, "Pups", "Particle 186m", "bearing")
```

```
## # A tibble: 2 × 3
##   Group         Mean_Bearing SD_Bearing
##   <chr>                <dbl>      <dbl>
## 1 Pups                 120.       0.969
## 2 Particle 186m         92.8      0.471
## 
##       Watson's Two-Sample Test of Homogeneity 
## 
## Test Statistic: 0.4125 
## P-value < 0.001 
## 
```

``` r
p3
```

![plot of chunk compare-bearings](figure/compare-bearings-3.png)

![plot of chunk save-combined-mean-bearings-comparison](figure/save-combined-mean-bearings-comparison-1.png)![plot of chunk save-combined-mean-bearings-comparison](figure/save-combined-mean-bearings-comparison-2.png)

## Categorise seals as following current if they are within 45 deg of the mean trace direction


``` r
### Categorise seals as following current if they are within 45 deg of the mean trace direction
locw_by_id <- locw_by_id %>%
    left_join(locps_by_id %>% st_drop_geometry() %>% select(id, bearing, sd_bearing), by = "id", suffix = c("_seal", "_pt0m")) %>%
    mutate(
        angle_diff_seal_pt0m = angle_diff(bearing_seal, bearing_pt0m),
        is_following_0m = angle_diff_seal_pt0m < 45
    )

# Add 186m particle data
locw_by_id <- locw_by_id %>%
    left_join(locp_by_id %>% st_drop_geometry() %>% select(id, bearing, sd_bearing) %>% rename(bearing_pt186m = bearing, sd_bearing_pt186m = sd_bearing), by = "id") %>%
    mutate(
        angle_diff_seal_pt186m = angle_diff(bearing_seal, bearing_pt186m),
        is_following_186m = angle_diff_seal_pt186m < 45
    )

locw_by_id %>%
    st_drop_geometry() %>%
    left_join(get_survival_data(), by = "id") %>%
    select(id, is_following_0m, is_following_186m, survive_trip_1, survive_year_1) %>%
    pivot_longer(
        cols = c(is_following_0m, is_following_186m),
        names_to = "particle_depth",
        values_to = "is_following"
    ) %>%
    mutate(particle_depth = if_else(particle_depth == "is_following_0m", "0m", "186m")) %>%
    group_by(particle_depth, is_following) %>%
    summarise(
        n = n(),
        prop = (n / nrow(locw_by_id)) %>% round(2),
        n_survived_trip_1 = sum(survive_trip_1, na.rm = TRUE),
        prop_survived_trip_1 = (n_survived_trip_1 / n()) %>% round(2),
        n_survived_year_1 = sum(survive_year_1, na.rm = TRUE),
        prop_survived_year_1 = (n_survived_year_1 / n()) %>% round(2),
        .groups = "drop"
    ) %>%
    mutate(is_following = if_else(is_following, "Following", "Not following")) %>%
    gt(groupname_col = "particle_depth") %>%
    fmt_percent(columns = c(prop, prop_survived_trip_1, prop_survived_year_1), decimals = 0) %>%
    cols_label(
        is_following = "Following Particle",
        n = "Count",
        prop = "Proportion",
        n_survived_trip_1 = "Survived Trip 1",
        prop_survived_trip_1 = "% Survived Trip 1",
        n_survived_year_1 = "Survived Year 1",
        prop_survived_year_1 = "% Survived Year 1"
    ) %>%
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_row_groups()
    ) %>%
    cols_merge(
        columns = c(n, prop),
        pattern = "{1} ({2})"
    ) %>%
    cols_merge(
        columns = c(n_survived_trip_1, prop_survived_trip_1),
        pattern = "{1} ({2})"
    ) %>%
    cols_merge(
        columns = c(n_survived_year_1, prop_survived_year_1),
        pattern = "{1} ({2})"
    ) %>%
    cols_label(
        n = "Count (%)",
        n_survived_trip_1 = "Survived Trip 1 (%)",
        n_survived_year_1 = "Survived Year 1 (%)"
    ) %>%
    gtsave(file.path(output_folder, "following_particle_survival.docx"))
```

## Calculate cumulative mean bearings


``` r
# Calculate cumulative mean bearings
calculate_circular_cummean <- function(bearings) {
    n <- length(bearings)

    sapply(1:n, function(i) {
        bearings_subset <- bearings[1:i]
        mean_bearing <- mean.circular(bearings_subset)
    })
}

# Combine weaner and particle data
w_pt <- locw %>%
    left_join(locps %>% st_drop_geometry() %>% select(id, date, bearing), by = c("id", "date"), suffix = c("_seal", "_pt")) %>%
    left_join(locp %>% st_drop_geometry() %>% select(id, date, bearing) %>% rename(bearing_pt186m = bearing), by = c("id", "date")) %>%
    group_by(id) %>%
    mutate(
        cummean_bearing_seal = calculate_circular_cummean(bearing_seal),
        cummean_bearing_pt = calculate_circular_cummean(bearing_pt),
        cummean_bearing_pt186m = calculate_circular_cummean(bearing_pt186m),
        days_since_start = as.numeric(date - first(date)),
        angle_diff_seal_pt = angle_diff(cummean_bearing_seal, cummean_bearing_pt),
        angle_diff_seal_pt186m = angle_diff(cummean_bearing_seal, cummean_bearing_pt186m),
    ) %>%
    left_join(locw_by_id %>% st_drop_geometry() %>% select(id, is_following_0m, is_following_186m), by = "id")
```

## Explore difference in cumulative bearings between seal and particle over time


``` r
# Plot angle difference over time
p1 <- w_pt %>%
    ggplot(aes(x = days_since_start, y = id)) +
    geom_tile(aes(fill = angle_diff_seal_pt)) +
    scale_fill_viridis_c(option = "B") +
    theme_bw() +
    facet_wrap(~is_following_0m, scales = "free", labeller = as_labeller(c(`TRUE` = "following", `FALSE` = "not following"))) +
    labs(x = "Days since start", y = "Seal ID", fill = "Δ Bearing - 0 m (°)") +
    theme(legend.position = "bottom", text = element_text(size = 8))

p1
```

![plot of chunk plot-angle-difference-over-time](figure/plot-angle-difference-over-time-1.png)

``` r
p2 <- w_pt %>%
    ggplot(aes(x = days_since_start, y = id)) +
    geom_tile(aes(fill = angle_diff_seal_pt186m)) +
    facet_wrap(~is_following_186m, scales = "free", labeller = as_labeller(c(`TRUE` = "following", `FALSE` = "not following"))) +
    scale_fill_viridis_c(option = "B") +
    labs(x = "Days since start", y = "Seal ID", fill = "Δ Bearing - 200 m (°)") +
    theme(legend.position = "bottom", text = element_text(size = 8))

p2
```

![plot of chunk plot-angle-difference-over-time](figure/plot-angle-difference-over-time-2.png)

![plot of chunk save-angle-difference-tile-plot](figure/save-angle-difference-tile-plot-1.png)

## Plot angle difference tracks

``` r
## Plot angle difference tracks

plots1 <- w_pt %>%
    rename(x1 = cummean_bearing_seal, x2 = cummean_bearing_pt, y = days_since_start) %>%
    select(id, is_following_0m, x1, x2, y) %>%
    group_split(id) %>%
    map(~ ggplot(data = .x) +
        geom_point(aes(x = x1, y = y, color = "Pups"), size = 0.2, alpha = .5) +
        geom_point(aes(x = x2, y = y, color = "Particle"), size = 0.2, alpha = .5) +
        coord_polar() +
        scale_x_continuous(limits = c(0, 360)) +
        scale_color_manual(values = c("Pups" = "#D55E00", "Particle" = "gray20")) +
        theme_bw() +
        theme(
            legend.position = "bottom", text = element_text(size = 8),
            panel.border = element_rect(color = ifelse(.x$is_following_0m[1], "#0072B2", "grey70"), linewidth = 1)
        ) +
        labs(x = NULL, y = "Days since start", color = "Group", subtitle = .x$id[1]))

o <- order(w_pt %>% group_by(id, is_following_0m) %>% nest() %>% pull(is_following_0m))

plots1 <- plots1[o]

p1 <- wrap_plots(plots1, guides = "collect") +
    plot_annotation(title = "0 m") &
    theme(legend.position = "bottom", axis.title = element_blank())

p1
```

![Angle difference tracks](figure/plot-angle-difference-tracks-1.png)

``` r
ggsave(file.path(output_folder, "angle_diff_tracks_0m.png"), p1, width = 10, height = 10, dpi = 300, bg = "white")


plots2 <- w_pt %>%
    rename(x1 = cummean_bearing_seal, x2 = cummean_bearing_pt186m, y = days_since_start) %>%
    select(id, is_following_186m, x1, x2, y) %>%
    group_split(id) %>%
    map(~ ggplot(data = .x) +
        geom_point(aes(x = x1, y = y, color = "Pups"), size = 0.2, alpha = .5) +
        geom_point(aes(x = x2, y = y, color = "Particle"), size = 0.2, alpha = .5) +
        coord_polar() +
        scale_x_continuous(limits = c(0, 360)) +
        scale_color_manual(values = c("Pups" = "#D55E00", "Particle" = "gray20")) +
        theme_bw() +
        theme(
            legend.position = "bottom", text = element_text(size = 8),
            panel.border = element_rect(color = ifelse(.x$is_following_186m[1], "#0072B2", "grey70"), linewidth = 1)
        ) +
        labs(x = NULL, y = "Days since start", color = "Group", subtitle = .x$id[1]))

o <- order(w_pt %>% group_by(id, is_following_186m) %>% nest() %>% pull(is_following_186m))

plots2 <- plots2[o]

p2 <- wrap_plots(plots2, guides = "collect") +
    plot_annotation(title = "200 m") &
    theme(legend.position = "bottom", axis.title = element_blank())

p2
```

![Angle difference tracks](figure/plot-angle-difference-tracks-2.png)

``` r
ggsave(file.path(output_folder, "angle_diff_tracks_186m.png"), p2, width = 10, height = 10, dpi = 300, bg = "white")
```

# Survival Model Analysis


``` r
# Prepare data for survival model
all_data_weaners <- load_all_data_weaners() # includes environment data

predictor_vars <- c("sst", "ssha", "eke", "tri", "slope", "SSTgrad", "ice", "dist_to_ice_m", "chl", "chlgrad")

model_data_all <- all_data_weaners %>%
    filter(id %in% w_pt$id, trip == 1, SUS == FALSE) %>%
    group_by(id) %>%
    summarise(across(any_of(predictor_vars), ~ mean(., na.rm = TRUE))) %>%
    left_join(get_survival_data(), by = "id") %>%
    left_join(
        w_pt %>%
            select(id, is_following_0m, is_following_186m) %>%
            st_drop_geometry() %>%
            group_by(id) %>%
            summarise(is_following_0m = last(is_following_0m), is_following_186m = last(is_following_186m)),
        by = "id"
    ) %>%
    select(-contains("seen"))

## Summarise model data
print(model_data_all %>% summary())
```

```
##       id                 sst              ssha                eke               tri            slope          SSTgrad               ice           dist_to_ice_m          chl            chlgrad          is_trip_complete    weanmass       birthyear    survive_trip_1  survive_year_1      size           is_following_0m is_following_186m
##  Length:48          Min.   :0.3145   Min.   :-0.131524   Min.   :0.01063   Min.   :20.05   Min.   :1.093   Min.   :6.895e-06   Min.   :0.000000   Min.   : 404874   Min.   :0.2639   Min.   :9.745e-08   Mode :logical    Min.   : 56.0   Min.   :1995   Mode :logical   Mode :logical   Length:48          Mode :logical   Mode :logical    
##  Class :character   1st Qu.:1.9070   1st Qu.:-0.021707   1st Qu.:0.02138   1st Qu.:37.58   1st Qu.:1.882   1st Qu.:9.499e-06   1st Qu.:0.000000   1st Qu.: 632734   1st Qu.:0.2839   1st Qu.:2.294e-07   FALSE:29         1st Qu.: 89.0   1st Qu.:1995   FALSE:16        FALSE:29        Class :character   FALSE:17        FALSE:16         
##  Mode  :character   Median :2.6798   Median : 0.006362   Median :0.02953   Median :41.99   Median :2.132   Median :1.046e-05   Median :0.000000   Median : 725079   Median :0.2951   Median :3.387e-07   TRUE :19         Median :102.0   Median :1996   TRUE :32        TRUE :19        Mode  :character   TRUE :31        TRUE :32         
##                     Mean   :3.1342   Mean   :-0.001090   Mean   :0.03943   Mean   :40.88   Mean   :2.099   Mean   :1.112e-05   Mean   :0.004963   Mean   : 738672   Mean   :0.3005   Mean   :4.657e-07                    Mean   :118.5   Mean   :1997                                                                                       
##                     3rd Qu.:4.2578   3rd Qu.: 0.016178   3rd Qu.:0.04157   3rd Qu.:46.71   3rd Qu.:2.388   3rd Qu.:1.199e-05   3rd Qu.:0.000000   3rd Qu.: 845446   3rd Qu.:0.3153   3rd Qu.:6.335e-07                    3rd Qu.:145.0   3rd Qu.:1999                                                                                       
##                     Max.   :7.2130   Max.   : 0.167460   Max.   :0.15998   Max.   :57.08   Max.   :2.960   Max.   :2.003e-05   Max.   :0.114613   Max.   :1153773   Max.   :0.3694   Max.   :2.597e-06                    Max.   :197.0   Max.   :2000                                                                                       
##                                                                                                                                                                                                                           NA's   :3
```

``` r
print(model_data_all %>% Hmisc::describe())
```

```
## . 
## 
##  19  Variables      48  Observations
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## id 
##        n  missing distinct 
##       48        0       48 
## 
## lowest : mq1-17215-95    mq1-17217-95    mq1-17219-95    mq1-20916-95    mq1-20918-95   , highest: mq4-Cleo-00     mq4-Doris-00    mq4-Ella-00     mq4-FirstOne-00 mq4-Flora-00   
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## sst 
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
##       48        0       48        1    3.134    1.843   0.9798   1.3815   1.9070   2.6798   4.2578   5.2382   6.0392 
## 
## lowest : 0.314476 0.692268 0.824426 1.26849  1.34853 , highest: 5.42762  5.50923  6.32455  6.88435  7.213   
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## ssha 
##         n   missing  distinct      Info      Mean       Gmd       .05       .10       .25       .50       .75       .90       .95 
##        48         0        48         1  -0.00109    0.0453 -0.082855 -0.058873 -0.021707  0.006362  0.016178  0.030506  0.045553 
## 
## lowest : -0.131524  -0.0839622 -0.083435  -0.0817769 -0.0760478, highest: 0.0317246  0.0410802  0.0479608  0.0707955  0.16746   
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## eke 
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
##       48        0       48        1  0.03943  0.02928  0.01326  0.01519  0.02138  0.02953  0.04157  0.07682  0.11338 
## 
## lowest : 0.0106258 0.0119708 0.0131632 0.0134405 0.0142879, highest: 0.0972405 0.110141  0.11513   0.123047  0.159983 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## tri 
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
##       48        0       48        1    40.88    9.463    24.35    27.91    37.58    41.99    46.71    49.10    53.72 
## 
## lowest : 20.049  22.6509 23.9863 25.0326 27.5868, highest: 49.3846 53.3999 53.8981 55.7982 57.084 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## slope 
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
##       48        0       48        1    2.099    0.487    1.239    1.512    1.882    2.132    2.388    2.506    2.896 
## 
## lowest : 1.09263 1.19855 1.21547 1.28326 1.38115, highest: 2.52443 2.85913 2.91525 2.92215 2.95994
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## SSTgrad 
##         n   missing  distinct      Info      Mean       Gmd       .05       .10       .25       .50       .75       .90       .95 
##        48         0        48         1 1.112e-05 2.894e-06 7.992e-06 8.521e-06 9.499e-06 1.046e-05 1.199e-05 1.519e-05 1.612e-05 
## 
## lowest : 6.89491e-06 6.93019e-06 7.94706e-06 8.0751e-06  8.34814e-06, highest: 1.5231e-05  1.5533e-05  1.6439e-05  1.68904e-05 2.00338e-05
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## ice 
##        n  missing distinct     Info     Mean      Gmd 
##       48        0        5     0.23 0.004963 0.009545 
##                                                                  
## Value      0.00000000 0.02597403 0.03966942 0.05797101 0.11461318
## Frequency          44          1          1          1          1
## Proportion      0.917      0.021      0.021      0.021      0.021
## 
## For the frequency table, variable is rounded to the nearest 0
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## dist_to_ice_m 
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
##       48        0       48        1   738672   182002   469810   557889   632734   725079   845446   935091   979479 
## 
## lowest : 404874  409554  450919  504892  549607 , highest: 955431  965330  987097  1074360 1153770
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## chl 
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
##       48        0       48        1   0.3005  0.02782   0.2701   0.2714   0.2839   0.2951   0.3153   0.3381   0.3459 
## 
## lowest : 0.263859 0.268478 0.269682 0.270794 0.271071, highest: 0.340619 0.343853 0.346996 0.363428 0.369352
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## chlgrad 
##         n   missing  distinct      Info      Mean       Gmd       .05       .10       .25       .50       .75       .90       .95 
##        48         0        48         1 4.657e-07 3.568e-07 1.306e-07 1.589e-07 2.294e-07 3.387e-07 6.335e-07 7.982e-07 8.768e-07 
## 
## lowest : 9.74526e-08 1.21839e-07 1.27012e-07 1.37174e-07 1.51993e-07, highest: 8.28255e-07 8.3738e-07  8.98033e-07 1.13201e-06 2.5967e-06 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## is_trip_complete 
##        n  missing distinct 
##       48        0        2 
##                       
## Value      FALSE  TRUE
## Frequency     29    19
## Proportion 0.604 0.396
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## weanmass 
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
##       45        3       37    0.999    118.5    41.58     78.4     82.2     89.0    102.0    145.0    171.4    179.0 
## 
## lowest :  56  78  80  81  84, highest: 173 175 180 195 197
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## birthyear 
##        n  missing distinct     Info     Mean      Gmd 
##       48        0        4    0.886     1997    2.308 
##                                   
## Value       1995  1996  1999  2000
## Frequency     21     5    14     8
## Proportion 0.438 0.104 0.292 0.167
## 
## For the frequency table, variable is rounded to the nearest 0
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## survive_trip_1 
##        n  missing distinct 
##       48        0        2 
##                       
## Value      FALSE  TRUE
## Frequency     16    32
## Proportion 0.333 0.667
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## survive_year_1 
##        n  missing distinct 
##       48        0        2 
##                       
## Value      FALSE  TRUE
## Frequency     29    19
## Proportion 0.604 0.396
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## size 
##        n  missing distinct 
##       48        0        3 
##                             
## Value        avg heavy light
## Frequency     10    18    20
## Proportion 0.208 0.375 0.417
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## is_following_0m 
##        n  missing distinct 
##       48        0        2 
##                       
## Value      FALSE  TRUE
## Frequency     17    31
## Proportion 0.354 0.646
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## is_following_186m 
##        n  missing distinct 
##       48        0        2 
##                       
## Value      FALSE  TRUE
## Frequency     16    32
## Proportion 0.333 0.667
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```

``` r
print(model_data_all %>% psych::describe())
```

```
##                   vars  n      mean        sd    median   trimmed       mad       min        max     range  skew kurtosis       se
## id*                  1 48     24.50     14.00     24.50     24.50     17.79      1.00      48.00     47.00  0.00    -1.28     2.02
## sst                  2 48      3.13      1.63      2.68      3.04      1.51      0.31       7.21      6.90  0.59    -0.37     0.24
## ssha                 3 48      0.00      0.05      0.01      0.00      0.02     -0.13       0.17      0.30  0.19     3.45     0.01
## eke                  4 48      0.04      0.03      0.03      0.03      0.01      0.01       0.16      0.15  2.09     4.02     0.00
## tri                  5 48     40.88      8.49     41.99     41.26      6.90     20.05      57.08     37.04 -0.53    -0.08     1.23
## slope                6 48      2.10      0.44      2.13      2.11      0.37      1.09       2.96      1.87 -0.29    -0.04     0.06
## SSTgrad              7 48      0.00      0.00      0.00      0.00      0.00      0.00       0.00      0.00  1.11     1.16     0.00
## ice                  8 48      0.00      0.02      0.00      0.00      0.00      0.00       0.11      0.11  4.36    20.00     0.00
## dist_to_ice_m        9 48 738671.61 160821.71 725078.76 737636.08 147581.23 404873.77 1153772.78 748899.01  0.18    -0.07 23212.61
## chl                 10 48      0.30      0.03      0.30      0.30      0.02      0.26       0.37      0.11  0.89     0.17     0.00
## chlgrad             11 48      0.00      0.00      0.00      0.00      0.00      0.00       0.00      0.00  3.23    14.33     0.00
## is_trip_complete    12 48       NaN        NA        NA       NaN        NA       Inf       -Inf      -Inf    NA       NA       NA
## weanmass            13 45    118.51     36.96    102.00    116.05     31.13     56.00     197.00    141.00  0.50    -1.05     5.51
## birthyear           14 48   1997.10      2.15   1996.00   1997.03      1.48   1995.00    2000.00      5.00  0.18    -1.86     0.31
## survive_trip_1      15 48       NaN        NA        NA       NaN        NA       Inf       -Inf      -Inf    NA       NA       NA
## survive_year_1      16 48       NaN        NA        NA       NaN        NA       Inf       -Inf      -Inf    NA       NA       NA
## size*               17 48      2.21      0.77      2.00      2.25      1.48      1.00       3.00      2.00 -0.36    -1.28     0.11
## is_following_0m     18 48       NaN        NA        NA       NaN        NA       Inf       -Inf      -Inf    NA       NA       NA
## is_following_186m   19 48       NaN        NA        NA       NaN        NA       Inf       -Inf      -Inf    NA       NA       NA
```



``` r
## Explore explanatory variables and their relationships with survival
pp <- ggpairs(all_data_weaners %>% select(all_of(c(predictor_vars, "survive_trip_1", "is_following_0m", "is_following_186m"))), aes(color = factor(survive_trip_1)))

ggsave(file.path(output_folder, "model_data_pairs_plot.png"), pp, width = 12, height = 10, scale = 1.5)
```



``` r
evaluate_model <- function(model) {
    cat("\n\n--- Summary of model_trip (global model) ---\n\n")
    print(model$call)
    tidy(model) %>% print()
    glance(model) %>% print()

    # Add DHARMa model checking for trip survival model
    cat("\n\n--- DHARMa model checking for trip survival model ---\n\n")
    dharma <- simulateResiduals(fittedModel = model, n = 1000)
    plot(dharma)
    testDispersion(dharma) %>% print()
    testZeroInflation(dharma) %>% print()
    testOutliers(dharma) %>% print()

    # model selection
    dredge <- dredge(model)

    best_models <- get.models(dredge, subset = delta <= 2)

    # summarise model selection
    .extract_model_info <- function(model, model_name) {
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

    model_summary <- map2_dfr(best_models, names(best_models), .extract_model_info) %>%
        mutate(delta_AICc = AICc - min(AICc)) %>% # Use AICc instead of AIC
        arrange(delta_AICc)

    cat("\n\n--- Dredge (model selection) summary ---\n\n")
    print(model_summary %>% as.data.frame())

    # model averaging
    avg_model <- model.avg(best_models)
    cat("\n\n--- Model averaging summary ---\n\n")
    print(summary(avg_model))

    # variable importance
    importance <- sw(dredge)
    cat("\n\n--- Variable importance ---\n\n")
    print(importance)
}
```

## Fit first trip survival model

``` r
model_data <- model_data_all %>%
    drop_na(weanmass) %>%
    mutate(across(where(is.numeric), ~ scale(.)))

model_trip <- glm(survive_trip_1 ~ is_following_0m + is_following_186m + weanmass + birthyear + sst + ssha + eke + slope + SSTgrad + ice + chlgrad,
    data = model_data, family = binomial(link = "logit"), na.action = "na.fail"
)

model <- model_trip

evaluate_model(model_trip)
```

```
## 
## 
## --- Summary of model_trip (global model) ---
## 
## glm(formula = survive_trip_1 ~ is_following_0m + is_following_186m + 
##     weanmass + birthyear + sst + ssha + eke + slope + SSTgrad + 
##     ice + chlgrad, family = binomial(link = "logit"), data = model_data, 
##     na.action = "na.fail")
## # A tibble: 12 × 5
##    term                  estimate std.error statistic p.value
##    <chr>                    <dbl>     <dbl>     <dbl>   <dbl>
##  1 (Intercept)             0.313      0.934    0.335   0.737 
##  2 is_following_0mTRUE     0.704      1.09     0.647   0.518 
##  3 is_following_186mTRUE   0.168      1.04     0.161   0.872 
##  4 weanmass                0.0228     0.476    0.0479  0.962 
##  5 birthyear               0.467      0.577    0.810   0.418 
##  6 sst                     0.931      0.832    1.12    0.263 
##  7 ssha                    1.32       0.849    1.56    0.120 
##  8 eke                     0.693      0.992    0.699   0.485 
##  9 slope                   0.588      0.528    1.11    0.265 
## 10 SSTgrad                -1.12       0.660   -1.70    0.0884
## 11 ice                    -0.110      0.433   -0.253   0.800 
## 12 chlgrad                 0.182      0.393    0.463   0.644 
## # A tibble: 1 × 8
##   null.deviance df.null logLik   AIC   BIC deviance df.residual  nobs
##           <dbl>   <int>  <dbl> <dbl> <dbl>    <dbl>       <int> <int>
## 1          58.6      44  -23.1  70.1  91.8     46.1          33    45
## 
## 
## --- DHARMa model checking for trip survival model ---
```

![plot of chunk fit-trip-survival-model](figure/fit-trip-survival-model-1.png)![plot of chunk fit-trip-survival-model](figure/fit-trip-survival-model-2.png)

```
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
## 
## data:  simulationOutput
## dispersion = 1.0257, p-value = 0.892
## alternative hypothesis: two.sided
```

![plot of chunk fit-trip-survival-model](figure/fit-trip-survival-model-3.png)

```
## 
## 	DHARMa zero-inflation test via comparison to expected zeros with simulation under H0 = fitted model
## 
## data:  simulationOutput
## ratioObsSim = 1.0057, p-value = 1
## alternative hypothesis: two.sided
```

![plot of chunk fit-trip-survival-model](figure/fit-trip-survival-model-4.png)

```
## 
## 	DHARMa bootstrapped outlier test
## 
## data:  dharma
## outliers at both margin(s) = 0, observations = 45, p-value = 1
## alternative hypothesis: two.sided
##  percent confidence interval:
##  0 0
## sample estimates:
## outlier frequency (expected: 0 ) 
##                                0
```

```
## Fixed term is "(Intercept)"
```

```
## 
## 
## --- Dredge (model selection) summary ---
## 
##    model_name                                                     formula     AICc      BIC deviance df.residual null.deviance df.null delta_AICc
## 1          17                        survive_trip_1 ~ is_following_0m + 1 58.27032 61.59793 53.98460          43      58.57363      44  0.0000000
## 2          65                                  survive_trip_1 ~ slope + 1 59.20669 62.53430 54.92098          43      58.57363      44  0.9363759
## 3          81                survive_trip_1 ~ is_following_0m + slope + 1 59.66646 64.50109 53.08110          42      58.57363      44  1.3961479
## 4         145                 survive_trip_1 ~ is_following_0m + ssha + 1 59.73652 64.57114 53.15116          42      58.57363      44  1.4662045
## 5          33                      survive_trip_1 ~ is_following_186m + 1 59.80832 63.13593 55.52261          43      58.57363      44  1.5380068
## 6         533        survive_trip_1 ~ eke + is_following_0m + SSTgrad + 1 59.88602 66.11267 50.88602          41      58.57363      44  1.6157009
## 7         529              survive_trip_1 ~ is_following_0m + SSTgrad + 1 59.90034 64.73497 53.31498          42      58.57363      44  1.6300270
## 8          21                  survive_trip_1 ~ eke + is_following_0m + 1 59.95555 64.79018 53.37019          42      58.57363      44  1.6852376
## 9         577                        survive_trip_1 ~ slope + SSTgrad + 1 60.02630 64.86092 53.44093          42      58.57363      44  1.7559792
## 10        193                           survive_trip_1 ~ slope + ssha + 1 60.08203 64.91665 53.49666          42      58.57363      44  1.8117116
## 11        273                  survive_trip_1 ~ is_following_0m + sst + 1 60.09842 64.93304 53.51305          42      58.57363      44  1.8280988
## 12         49    survive_trip_1 ~ is_following_0m + is_following_186m + 1 60.09844 64.93306 53.51307          42      58.57363      44  1.8281220
## 13        661 survive_trip_1 ~ eke + is_following_0m + ssha + SSTgrad + 1 60.10013 67.59498 48.56167          40      58.57363      44  1.8298135
## 
## 
## --- Model averaging summary ---
## 
## 
## Call:
## model.avg(object = best_models)
## 
## Component model call: 
## glm(formula = survive_trip_1 ~ <13 unique rhs>, family = binomial(link = "logit"), data = model_data, na.action = na.fail)
## 
## Component models: 
##      df logLik  AICc delta weight
## 2     2 -26.99 58.27  0.00   0.16
## 4     2 -27.46 59.21  0.94   0.10
## 24    3 -26.54 59.67  1.40   0.08
## 25    3 -26.58 59.74  1.47   0.07
## 3     2 -27.76 59.81  1.54   0.07
## 127   4 -25.44 59.89  1.62   0.07
## 27    3 -26.66 59.90  1.63   0.07
## 12    3 -26.69 59.96  1.69   0.07
## 47    3 -26.72 60.03  1.76   0.06
## 45    3 -26.75 60.08  1.81   0.06
## 26    3 -26.76 60.10  1.83   0.06
## 23    3 -26.76 60.10  1.83   0.06
## 1257  5 -24.28 60.10  1.83   0.06
## 
## Term codes: 
##               eke   is_following_0m is_following_186m             slope              ssha               sst           SSTgrad 
##                 1                 2                 3                 4                 5                 6                 7 
## 
## Model-averaged coefficients:  
## (full average) 
##                       Estimate Std. Error Adjusted SE z value Pr(>|z|)
## (Intercept)           -0.02991    0.62232     0.63392   0.047    0.962
## is_following_0mTRUE    0.95279    0.86611     0.87816   1.085    0.278
## slope                  0.17749    0.34363     0.34714   0.511    0.609
## ssha                   0.09413    0.27919     0.28320   0.332    0.740
## is_following_186mTRUE  0.11723    0.41371     0.41879   0.280    0.780
## eke                    0.14139    0.41796     0.42316   0.334    0.738
## SSTgrad               -0.15009    0.35643     0.36046   0.416    0.677
## sst                    0.01465    0.10313     0.10527   0.139    0.889
##  
## (conditional average) 
##                       Estimate Std. Error Adjusted SE z value Pr(>|z|)  
## (Intercept)           -0.02991    0.62232     0.63392   0.047   0.9624  
## is_following_0mTRUE    1.35713    0.72093     0.74140   1.830   0.0672 .
## slope                  0.58527    0.38821     0.39839   1.469   0.1418  
## ssha                   0.46941    0.46102     0.47306   0.992   0.3211  
## is_following_186mTRUE  0.86920    0.78447     0.80423   1.081   0.2798  
## eke                    0.70958    0.68816     0.70391   1.008   0.3134  
## SSTgrad               -0.56427    0.49387     0.50476   1.118   0.2636  
## sst                    0.23415    0.34442     0.35464   0.660   0.5091  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## 
## 
## --- Variable importance ---
## 
##                      is_following_0m SSTgrad slope ssha is_following_186m eke  sst  ice  chlgrad birthyear weanmass
## Sum of weights:      0.49            0.44    0.42  0.40 0.33              0.31 0.28 0.24 0.23    0.23      0.22    
## N containing models: 1024            1024    1024  1024 1024              1024 1024 1024 1024    1024      1024
```

## Fit first year survival model

``` r
## Year survival
model_year <- glm(survive_year_1 ~ is_following_0m + is_following_186m + weanmass + birthyear + sst + ssha + eke + slope + SSTgrad + ice + chlgrad,
    data = model_data, family = binomial(link = "logit"), na.action = "na.fail"
)

evaluate_model(model_year)
```

```
## 
## 
## --- Summary of model_trip (global model) ---
## 
## glm(formula = survive_year_1 ~ is_following_0m + is_following_186m + 
##     weanmass + birthyear + sst + ssha + eke + slope + SSTgrad + 
##     ice + chlgrad, family = binomial(link = "logit"), data = model_data, 
##     na.action = "na.fail")
## # A tibble: 12 × 5
##    term                  estimate std.error statistic p.value
##    <chr>                    <dbl>     <dbl>     <dbl>   <dbl>
##  1 (Intercept)            -1.81       1.15    -1.58    0.114 
##  2 is_following_0mTRUE     1.62       1.28     1.27    0.205 
##  3 is_following_186mTRUE   0.0464     1.28     0.0364  0.971 
##  4 weanmass                1.14       0.499    2.29    0.0219
##  5 birthyear              -0.374      0.552   -0.677   0.498 
##  6 sst                     0.279      0.756    0.369   0.712 
##  7 ssha                    0.717      0.521    1.38    0.169 
##  8 eke                     0.364      0.704    0.517   0.605 
##  9 slope                   0.325      0.489    0.665   0.506 
## 10 SSTgrad                -0.0907     0.536   -0.169   0.866 
## 11 ice                    -0.114      0.557   -0.204   0.838 
## 12 chlgrad                -0.727      0.789   -0.921   0.357 
## # A tibble: 1 × 8
##   null.deviance df.null logLik   AIC   BIC deviance df.residual  nobs
##           <dbl>   <int>  <dbl> <dbl> <dbl>    <dbl>       <int> <int>
## 1          59.7      44  -22.5  68.9  90.6     44.9          33    45
## 
## 
## --- DHARMa model checking for trip survival model ---
```

![plot of chunk fit-year-survival-model](figure/fit-year-survival-model-1.png)![plot of chunk fit-year-survival-model](figure/fit-year-survival-model-2.png)

```
## 
## 	DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
## 
## data:  simulationOutput
## dispersion = 1.031, p-value = 0.85
## alternative hypothesis: two.sided
```

![plot of chunk fit-year-survival-model](figure/fit-year-survival-model-3.png)

```
## 
## 	DHARMa zero-inflation test via comparison to expected zeros with simulation under H0 = fitted model
## 
## data:  simulationOutput
## ratioObsSim = 0.99312, p-value = 1
## alternative hypothesis: two.sided
```

![plot of chunk fit-year-survival-model](figure/fit-year-survival-model-4.png)

```
## 
## 	DHARMa bootstrapped outlier test
## 
## data:  dharma
## outliers at both margin(s) = 0, observations = 45, p-value = 1
## alternative hypothesis: two.sided
##  percent confidence interval:
##  0 0
## sample estimates:
## outlier frequency (expected: 0 ) 
##                                0
```

```
## Fixed term is "(Intercept)"
```

```
## 
## 
## --- Dredge (model selection) summary ---
## 
##    model_name                                                                       formula     AICc      BIC deviance df.residual null.deviance df.null delta_AICc
## 1        1217                                  survive_year_1 ~ slope + ssha + weanmass + 1 59.30001 65.52666 50.30001          41      59.66692      44 0.00000000
## 2        1425                  survive_year_1 ~ is_following_0m + ssha + sst + weanmass + 1 59.32833 66.82318 47.78987          40      59.66692      44 0.02831627
## 3        1173                  survive_year_1 ~ eke + is_following_0m + ssha + weanmass + 1 59.47861 66.97346 47.94014          40      59.66692      44 0.17859472
## 4        1169                        survive_year_1 ~ is_following_0m + ssha + weanmass + 1 59.70655 65.93320 50.70655          41      59.66692      44 0.40654175
## 5        1047          survive_year_1 ~ chlgrad + eke + is_following_0m + weanmass +      1 59.83252 67.32737 48.29406          40      59.66692      44 0.53251014
## 6        1473                            survive_year_1 ~ slope + ssha + sst + weanmass + 1 59.99098 67.48584 48.45252          40      59.66692      44 0.69097360
## 7        1045                         survive_year_1 ~ eke + is_following_0m + weanmass + 1 60.24620 66.47285 51.24620          41      59.66692      44 0.94618574
## 8        1043                     survive_year_1 ~ chlgrad + is_following_0m + weanmass + 1 60.29998 66.52663 51.29998          41      59.66692      44 0.99997170
## 9        1153                                          survive_year_1 ~ ssha + weanmass + 1 60.31085 65.14548 53.72549          42      59.66692      44 1.01084363
## 10       1041                               survive_year_1 ~ is_following_0m + weanmass + 1 60.59700 65.43162 54.01164          42      59.66692      44 1.29699211
## 11       1299          survive_year_1 ~ chlgrad + is_following_0m + sst + weanmass +      1 60.66698 68.16183 49.12852          40      59.66692      44 1.36696714
## 12       1427   survive_year_1 ~ chlgrad + is_following_0m + ssha + sst + weanmass +      1 60.68489 69.31434 46.47437          39      59.66692      44 1.38488277
## 13       1218                      survive_year_1 ~ birthyear + slope + ssha + weanmass + 1 60.69941 68.19426 49.16095          40      59.66692      44 1.39939574
## 14       1175   survive_year_1 ~ chlgrad + eke + is_following_0m + ssha + weanmass +      1 60.73187 69.36132 46.52134          39      59.66692      44 1.43186001
## 15       1089                                         survive_year_1 ~ slope + weanmass + 1 60.76038 65.59500 54.17501          42      59.66692      44 1.46036529
## 16       1221                            survive_year_1 ~ eke + slope + ssha + weanmass + 1 60.76448 68.25934 49.22602          40      59.66692      44 1.46447312
## 17       1219                        survive_year_1 ~ chlgrad + slope + ssha + weanmass + 1 60.86178 68.35663 49.32332          40      59.66692      44 1.56176762
## 18       1171         survive_year_1 ~ chlgrad + is_following_0m + ssha + weanmass +      1 60.87568 68.37053 49.33721          40      59.66692      44 1.57566526
## 19        193                                             survive_year_1 ~ slope + ssha + 1 60.89231 65.72694 54.30695          42      59.66692      44 1.59230372
## 20       1170       survive_year_1 ~ birthyear + is_following_0m + ssha + weanmass +      1 60.91938 68.41423 49.38092          40      59.66692      44 1.61936829
## 21       1489     survive_year_1 ~ is_following_0m + slope + ssha + sst + weanmass +      1 60.95446 69.58391 46.74394          39      59.66692      44 1.65445241
## 22       1091                               survive_year_1 ~ chlgrad + slope + weanmass + 1 60.96738 67.19403 51.96738          41      59.66692      44 1.66736874
## 23       1157                                    survive_year_1 ~ eke + ssha + weanmass + 1 60.97137 67.19802 51.97137          41      59.66692      44 1.67136193
## 24       1233           survive_year_1 ~ is_following_0m + slope + ssha + weanmass +      1 60.98930 68.48415 49.45084          40      59.66692      44 1.68929191
## 25       1025                                                 survive_year_1 ~ weanmass + 1 61.22367 64.55128 56.93796          43      59.66692      44 1.92365836
## 26       1174 survive_year_1 ~ birthyear + eke + is_following_0m + ssha + weanmass +      1 61.28456 69.91401 47.07403          39      59.66692      44 1.98454661
## 
## 
## --- Model averaging summary ---
## 
## 
## Call:
## model.avg(object = best_models)
## 
## Component model call: 
## glm(formula = survive_year_1 ~ <26 unique rhs>, family = binomial(link = "logit"), data = model_data, na.action = na.fail)
## 
## Component models: 
##       df logLik  AICc delta weight
## 568    4 -25.15 59.30  0.00   0.07
## 4678   5 -23.89 59.33  0.03   0.07
## 3468   5 -23.97 59.48  0.18   0.06
## 468    4 -25.35 59.71  0.41   0.06
## 2348   5 -24.15 59.83  0.53   0.05
## 5678   5 -24.23 59.99  0.69   0.05
## 348    4 -25.62 60.25  0.95   0.04
## 248    4 -25.65 60.30  1.00   0.04
## 68     3 -26.86 60.31  1.01   0.04
## 48     3 -27.01 60.60  1.30   0.04
## 2478   5 -24.56 60.67  1.37   0.03
## 24678  6 -23.24 60.68  1.38   0.03
## 1568   5 -24.58 60.70  1.40   0.03
## 23468  6 -23.26 60.73  1.43   0.03
## 58     3 -27.09 60.76  1.46   0.03
## 3568   5 -24.61 60.76  1.46   0.03
## 2568   5 -24.66 60.86  1.56   0.03
## 2468   5 -24.67 60.88  1.58   0.03
## 56     3 -27.15 60.89  1.59   0.03
## 1468   5 -24.69 60.92  1.62   0.03
## 45678  6 -23.37 60.95  1.65   0.03
## 258    4 -25.98 60.97  1.67   0.03
## 368    4 -25.99 60.97  1.67   0.03
## 4568   5 -24.73 60.99  1.69   0.03
## 8      2 -28.47 61.22  1.92   0.03
## 13468  6 -23.54 61.28  1.98   0.03
## 
## Term codes: 
##       birthyear         chlgrad             eke is_following_0m           slope            ssha             sst        weanmass 
##               1               2               3               4               5               6               7               8 
## 
## Model-averaged coefficients:  
## (full average) 
##                     Estimate Std. Error Adjusted SE z value Pr(>|z|)  
## (Intercept)         -1.22976    0.82217     0.83601   1.471   0.1413  
## slope                0.22988    0.39109     0.39567   0.581   0.5612  
## ssha                 0.50887    0.50065     0.50901   1.000   0.3174  
## weanmass             0.78189    0.42738     0.43754   1.787   0.0739 .
## is_following_0mTRUE  0.95217    1.05545     1.06887   0.891   0.3730  
## sst                  0.14218    0.34562     0.34963   0.407   0.6843  
## eke                  0.15434    0.32251     0.32637   0.473   0.6363  
## chlgrad             -0.21778    0.49360     0.50131   0.434   0.6640  
## birthyear           -0.04073    0.18777     0.19083   0.213   0.8310  
##  
## (conditional average) 
##                     Estimate Std. Error Adjusted SE z value Pr(>|z|)  
## (Intercept)          -1.2298     0.8222      0.8360   1.471   0.1413  
## slope                 0.6325     0.4076      0.4196   1.507   0.1317  
## ssha                  0.7188     0.4507      0.4638   1.550   0.1211  
## weanmass              0.8065     0.4106      0.4215   1.913   0.0557 .
## is_following_0mTRUE   1.5890     0.9204      0.9459   1.680   0.0930 .
## sst                   0.6708     0.4572      0.4713   1.423   0.1547  
## eke                   0.5600     0.3876      0.3991   1.403   0.1606  
## chlgrad              -0.7647     0.6613      0.6813   1.122   0.2617  
## birthyear            -0.4592     0.4531      0.4673   0.983   0.3258  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## 
## 
## --- Variable importance ---
## 
##                      weanmass ssha is_following_0m slope chlgrad eke  sst  birthyear is_following_186m ice  SSTgrad
## Sum of weights:      0.73     0.55 0.49            0.39  0.38    0.33 0.30 0.25      0.25              0.24 0.23   
## N containing models: 1024     1024 1024            1024  1024    1024 1024 1024      1024              1024 1024
```


## Create summary table for all pups

```
## Joining with `by = join_by(id, is_trip_complete, seen_6m, seen_1y, weanmass, birthyear)`
## `summarise()` has grouped output by 'id', 'tripdur', 'birthdate', 'weanmass', 'survive_trip_1'. You can override using the `.groups` argument.
```
