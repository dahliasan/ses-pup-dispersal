---
title: Data Export
layout: default
---

---
title: "Export data for publication"
format: html
---


``` r
library(tidyverse)
library(here)
library(sf)

# Load data
files <- list(
    seal_data_path = here("output", "tracks_processed_12h.rds"),
    female_data_path = here("data", "adult_female_locs.rds"),
    particle_data_path = here("output", "particle-trace-186.13m.rds"),
    surface_particle_data_path = here("output", "currently_particleTrace.rds")
)
```



``` r
load(here("output", "seal_id.rdata"))
ls()
```

```
##  [1] "add_frontmatter"               "all_bearings_plot"             "all_data_weaners"              "all_tracks_plot"               "angle_diff"                    "bb"                            "calcMeanDepartureDate"         "calculate_bearings"            "calculate_circular_cummean"    "calculate_days_since_start"    "calculate_distance_from_start" "colony"                        "compare_bearings"              "convert2polarsf"               "crs"                           "eastOrWest"                    "evaluate_model"                "female_data_path"              "files"                         "get_survival_data"             "id_keep"                       "identify_outbound_trip"        "interpolate_track"             "load_all_data_weaners"         "locf"                          "locf_by_id"                    "locp"                          "locp_by_id"                    "locps"                         "locps_by_id"                   "locw"                          "locw_by_id"                    "make_circular"                 "mean_bearings_plot"            "model"                         "model_data"                    "model_data_all"                "model_trip"                    "model_year"                    "mq"                            "mq_sf"                         "normalize_longitude"           "o"                             "orsi_sf"                       "output_folder"                 "p_list"                        "p_names"                       "p1"                            "p2"                            "p3"                            "palette1"                      "palette2"                      "palette3"                      "particle_data_path"            "plot_bearing_histogram"        "plots1"                        "plots2"                        "predictor_vars"                "preprocess_seal_particle_data" "preprocess_track"              "print_and_save_results"        "process_by_id"                 "proj"                          "saveDepartureDateResult"       "seal_data_path"                "sealID_all"                    "surface_particle_data_path"    "unwrapLon"                     "w_pt"                          "whichZone"                     "world_sf"                      "wrapLon"
```

``` r
sealID_all %>% print(n = 100)
```

```
## # A tibble: 98 × 6
##    ref               ptt brand tag1  tag2  SEAL_ID
##    <chr>           <dbl> <chr> <chr> <chr>   <dbl>
##  1 mq1-2849-95      2849 J373  9788  <NA>     8108
##  2 mq1-5811-95      5811 J032  8048  8049     7254
##  3 mq1-5814-95      5814 J239  8392  8393     7375
##  4 mq1-17215-95    17215 J622  9500  9501     7932
##  5 mq1-20916-95    20916 J723  8728  8729     7582
##  6 mq1-20918-95    20918 J503  8203  <NA>     7308
##  7 mq1-22483-95    22483 J806  9538  9539     7947
##  8 mq1-22484-95    22484 J557  8874  8875     7627
##  9 mq1-22486-95    22486 J781  8870  8871     7625
## 10 mq1-22490-95    22490 J429  8927  8928     7664
## 11 mq1-22500-95    22500 J803  8340  8341     7392
## 12 mq1-22501-95    22501 J012  8258  8259     7357
## 13 mq1-26625-95    26625 J766  8114  8115     7279
## 14 mq1-26628-95    26628 J483  9230  9231     7855
## 15 mq1-26629-95    26629 J708  8598  8599     7497
## 16 mq1-26633-95    26633 J591  9017  9018     7710
## 17 mq1-26635-95    26635 J375  8742  8743     7542
## 18 mq1-17217-95    17217 J226  9346  9347     7915
## 19 mq1-17219-95    17219 J615  8040  8041     7244
## 20 mq1-22499-95    22499 J415  9711  9712     8055
## 21 mq1-26627-95    26627 J091  8607  <NA>     7501
## 22 mq1-2846-95      2846 J100  9027  9028     7698
## 23 mq1-5810-95      5810 J219  8310  8311     7353
## 24 mq1-5815-95      5815 J332  9069  9070     7761
## 25 mq1-5816-95      5816 J604  8957  8958     7726
## 26 mq1-17213-95    17213 J627  9173  9174     7780
## 27 mq1-17216-95    17216 J336  9240  9241     7823
## 28 mq1-17221-95    17221 J387  8722  8723     7579
## 29 mq1-20915-95    20915 J026  8282  8283     7369
## 30 mq1-20917-95    20917 J239  8392  8393     7375
## 31 mq1-22482-95    22482 J622  9500  9501     7932
## 32 mq1-22485-95    22485 J575  8696  8697     7570
## 33 mq1-22487-95    22487 J243  8292  8293     7326
## 34 mq1-22488-95    22488 J693  8896  8897     7671
## 35 mq1-22497-95    22497 J726  8618  8619     7522
## 36 mq1-26624-95    26624 J362  8935  8936     7668
## 37 mq1-26626-95    26626 J356  8453  8454     7445
## 38 mq1-26631-95    26631 J708  8598  8599     7497
## 39 mq1-26632-95    26632 J803  8340  8341     7392
## 40 mq2-20916-96    20916 K887  <NA>  <NA>     9774
## 41 mq2-20917-96    20917 K354  966   967      9748
## 42 mq2-22500-96    22500 K339  669   670      9574
## 43 mq2-26623-96    26623 K397  795   796      9637
## 44 mq2-28482-96    28482 K223  725   726      9602
## 45 mq2-28479-96    28479 K875  2359  2360    10450
## 46 mq3-2846-99      2846 T232  <NA>  <NA>    25783
## 47 mq3-2849-99      2849 T825  <NA>  <NA>    26468
## 48 mq3-5812-99      5812 T887  <NA>  <NA>    26605
## 49 mq3-17217-99    17217 T088  <NA>  <NA>    26415
## 50 mq3-20918-99    20918 T875  <NA>  <NA>    26290
## 51 mq3-22484-99    22484 T806  <NA>  <NA>    25910
## 52 mq3-22488-99    22488 T899  <NA>  <NA>    25966
## 53 mq3-22498-99    22498 T857  <NA>  <NA>    25859
## 54 mq3-26627-99    26627 T867  <NA>  <NA>    26244
## 55 mq3-26629-99    26629 T773  <NA>  <NA>    26535
## 56 mq3-28494-99    28494 T807  <NA>  <NA>    26342
## 57 mq3-28496-99    28496 T823  <NA>  <NA>    26139
## 58 mq3-28497-99    28497 T839  <NA>  <NA>    25949
## 59 mq3-28500-99    28500 T772  <NA>  <NA>    26560
## 60 mq3-28504-99    28504 T719  <NA>  <NA>    25673
## 61 mq3-2849a-99     2849 T825  <NA>  <NA>    26468
## 62 mq3-5812a-99     5812 T887  <NA>  <NA>    26605
## 63 mq3-17217a-99   17217 T088  <NA>  <NA>    26415
## 64 mq3-28496a-99   28496 T823  <NA>  <NA>    26139
## 65 mq3-28504a-99   28504 T719  <NA>  <NA>    25673
## 66 mq2-22483-96    22483 <NA>  <NA>  <NA>       NA
## 67 mq2-22490-96    22490 <NA>  <NA>  <NA>       NA
## 68 mq2-26624-96    26624 <NA>  <NA>  <NA>       NA
## 69 mq2-22486-96    22486 <NA>  <NA>  <NA>       NA
## 70 mq2-22488-96    22488 <NA>  <NA>  <NA>       NA
## 71 mq2-22497-96    22497 <NA>  <NA>  <NA>       NA
## 72 mq2-22499-96    22499 <NA>  <NA>  <NA>       NA
## 73 mq2-28483-96    28483 <NA>  <NA>  <NA>       NA
## 74 mq3-2841-99      2841 <NA>  <NA>  <NA>       NA
## 75 mq3-2845-99      2845 <NA>  <NA>  <NA>       NA
## 76 mq3-26635-99    26635 <NA>  <NA>  <NA>       NA
## 77 mq3-2442-99      2442 <NA>  <NA>  <NA>       NA
## 78 mq3-2841a-99     2841 <NA>  <NA>  <NA>       NA
## 79 mq3-2841b-99     2841 <NA>  <NA>  <NA>       NA
## 80 mq3-2841c-99     2841 <NA>  <NA>  <NA>       NA
## 81 mq3-2845a-99     2845 <NA>  <NA>  <NA>       NA
## 82 mq3-2845b-99     2845 <NA>  <NA>  <NA>       NA
## 83 mq3-2845c-99     2845 <NA>  <NA>  <NA>       NA
## 84 mq3-22500-99    22500 <NA>  <NA>  <NA>       NA
## 85 mq3-26635a-99   26635 <NA>  <NA>  <NA>       NA
## 86 mq3-26635b-99   26635 <NA>  <NA>  <NA>       NA
## 87 mq3-26635c-99   26635 <NA>  <NA>  <NA>       NA
## 88 mq3-28502-99    28502 <NA>  <NA>  <NA>       NA
## 89 mq3-28505a-99   28505 <NA>  <NA>  <NA>       NA
## 90 mq3-28505-99    28505 <NA>  <NA>  <NA>       NA
## 91 mq4-Alice-00     1547 <NA>  W3762 W3763   28562
## 92 mq4-Billie-00   17214 <NA>  W3998 W3999   28680
## 93 mq4-Cleo-00     17216 <NA>  W2706 W2707   28029
## 94 mq4-FirstOne-00 17217 <NA>  W3000 W3001   28181
## 95 mq4-Doris-00    20916 <NA>  W4197 W4198   28775
## 96 mq4-20918-00    20918 <NA>  W4191 W4231   28783
## 97 mq4-Ella-00     22489 <NA>  W3325 W3326   28345
## 98 mq4-Flora-00    26624 <NA>  W3843 W3844   28603
```



``` r
get_dataset_info <- function(file_path, dataset_name) {
    data <- read_rds(file_path) %>%
        ungroup() %>%
        st_drop_geometry()

    # Time coverage using reframe instead of summarise
    time_range <- data %>%
        reframe(
            start_date = min(date, na.rm = TRUE),
            end_date = max(date, na.rm = TRUE),
            n_records = n()
        )

    # Spatial coverage using reframe
    spatial_range <- data %>%
        reframe(
            lon_min = min(lon, na.rm = TRUE),
            lon_max = max(lon, na.rm = TRUE),
            lat_min = min(lat, na.rm = TRUE),
            lat_max = max(lat, na.rm = TRUE)
        )

    # Variable names and classes (removing duplicates)
    var_info <- data %>%
        summarise(across(everything(), ~ class(.x)[1])) %>% # Take only first class
        pivot_longer(everything(),
            names_to = "variable",
            values_to = "class"
        ) %>%
        distinct() # Remove any remaining duplicates

    # Number of unique IDs
    n_ids <- data %>%
        pull(id) %>%
        n_distinct()

    # Print summary
    cat("\n=== Dataset:", dataset_name, "===\n")
    cat("\nTime Coverage:\n")
    print(time_range)

    cat("\nSpatial Coverage:\n")
    print(spatial_range)

    cat("\nNumber of unique IDs:", n_ids, "\n")

    cat("\nVariables and their classes:\n")
    print(var_info)

    cat("\nSample size:", nrow(data), "records\n")

    # Return invisibly for potential further use
    invisible(list(
        time_range = time_range,
        spatial_range = spatial_range,
        var_info = var_info,
        n_ids = n_ids,
        data = data # Include the data for CSV export
    ))
}
```


``` r
# Get metadata for each dataset
seal_info <- get_dataset_info(files$seal_data_path, "seal_data")
```

```
## 
## === Dataset: seal_data ===
## 
## Time Coverage:
## # A tibble: 1 × 3
##   start_date          end_date            n_records
##   <dttm>              <dttm>                  <int>
## 1 1995-11-27 02:00:00 2002-01-05 18:00:00     14162
## 
## Spatial Coverage:
## # A tibble: 1 × 4
##   lon_min lon_max lat_min lat_max
##     <dbl>   <dbl>   <dbl>   <dbl>
## 1   -180.    180.   -67.0   -43.8
## 
## Number of unique IDs: 49 
## 
## Variables and their classes:
## # A tibble: 12 × 2
##    variable           class    
##    <chr>              <chr>    
##  1 id                 character
##  2 date               POSIXct  
##  3 land               logical  
##  4 lon                numeric  
##  5 lat                numeric  
##  6 daysFromDeployment numeric  
##  7 haulout            logical  
##  8 dist2col           numeric  
##  9 tripdur            difftime 
## 10 trip               numeric  
## 11 is_trip_complete   logical  
## 12 SUS                logical  
## 
## Sample size: 14162 records
```

``` r
female_info <- get_dataset_info(files$female_data_path, "female_data")
```

```
## 
## === Dataset: female_data ===
## 
## Time Coverage:
## # A tibble: 1 × 3
##   start_date          end_date            n_records
##   <dttm>              <dttm>                  <int>
## 1 2000-02-03 00:00:00 2010-10-24 19:31:03    101905
## 
## Spatial Coverage:
## # A tibble: 1 × 4
##   lon_min lon_max lat_min lat_max
##     <dbl>   <dbl>   <dbl>   <dbl>
## 1    110.    241.   -76.5   -43.1
## 
## Number of unique IDs: 67 
## 
## Variables and their classes:
## # A tibble: 7 × 2
##   variable class    
##   <chr>    <chr>    
## 1 id       factor   
## 2 date     POSIXct  
## 3 lon      numeric  
## 4 lat      numeric  
## 5 type     character
## 6 month    numeric  
## 7 season   character
## 
## Sample size: 101905 records
```

``` r
particle_info <- get_dataset_info(files$particle_data_path, "particle_data")
```

```
## 
## === Dataset: particle_data ===
## 
## Time Coverage:
## # A tibble: 1 × 3
##   start_date          end_date            n_records
##   <dttm>              <dttm>                  <int>
## 1 1995-11-29 00:00:00 2001-10-03 00:00:00      5603
## 
## Spatial Coverage:
## # A tibble: 1 × 4
##   lon_min lon_max lat_min lat_max
##     <dbl>   <dbl>   <dbl>   <dbl>
## 1    158.    173.   -61.4   -52.3
## 
## Number of unique IDs: 48 
## 
## Variables and their classes:
## # A tibble: 6 × 2
##   variable       class    
##   <chr>          <chr>    
## 1 lon            numeric  
## 2 lat            numeric  
## 3 group          integer  
## 4 date           POSIXct  
## 5 id             character
## 6 date_processed Date     
## 
## Sample size: 5603 records
```

``` r
surface_particle_info <- get_dataset_info(files$surface_particle_data_path, "surface_particle_data")
```

```
## 
## === Dataset: surface_particle_data ===
## 
## Time Coverage:
## # A tibble: 1 × 3
##   start_date          end_date            n_records
##   <dttm>              <dttm>                  <int>
## 1 1995-11-29 11:00:00 2001-10-03 10:00:00     11164
## 
## Spatial Coverage:
## # A tibble: 1 × 4
##   lon_min lon_max lat_min lat_max
##     <dbl>   <dbl>   <dbl>   <dbl>
## 1   -180.    180.   -65.3   -47.3
## 
## Number of unique IDs: 48 
## 
## Variables and their classes:
## # A tibble: 7 × 2
##   variable class    
##   <chr>    <chr>    
## 1 x        numeric  
## 2 y        numeric  
## 3 group    integer  
## 4 lon      numeric  
## 5 lat      numeric  
## 6 date     POSIXct  
## 7 id       character
## 
## Sample size: 11164 records
```


``` r
# Export CSVs
seal_info$data %>%
    ungroup() %>%
    st_drop_geometry() %>%
    rename(flag = SUS) %>%
    select(id, date, lon, lat, flag) %>%
    write_csv(paste0(tools::file_path_sans_ext(files$seal_data_path), ".csv"))

female_info$data %>%
    select(id, date, lon, lat, type, month, season) %>%
    write_csv(paste0(tools::file_path_sans_ext(files$female_data_path), ".csv"))

particle_info$data %>%
    select(id, date, lon, lat) %>%
    write_csv(paste0(tools::file_path_sans_ext(files$particle_data_path), ".csv"))

surface_particle_info$data %>%
    select(id, date, lon, lat) %>%
    write_csv(paste0(tools::file_path_sans_ext(files$surface_particle_data_path), ".csv"))
```


``` r
# Combine all metadata into a list
all_metadata <- list(
    seal_data = seal_info,
    female_data = female_info,
    particle_data = particle_info,
    surface_particle_data = surface_particle_info
)

# Save metadata as RDS for potential future use
write_rds(all_metadata, here("output", "dataset_metadata.rds"))

# Create a simplified metadata summary as CSV
metadata_summary <- bind_rows(
    seal_info$time_range %>% mutate(dataset = "seal_data"),
    female_info$time_range %>% mutate(dataset = "female_data"),
    particle_info$time_range %>% mutate(dataset = "particle_data"),
    surface_particle_info$time_range %>% mutate(dataset = "surface_particle_data")
) %>%
    select(dataset, everything())

write_csv(metadata_summary, here("output", "dataset_metadata_summary.csv"))
#
```


