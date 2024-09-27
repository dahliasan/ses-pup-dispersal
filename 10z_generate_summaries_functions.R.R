require(dplyr)
require(lubridate)
require(readr)
require(sf)

### Generate Summaries -----------------------------------------------------

prepIndividualDf <- function(x, group_name, bearing, y) {
    newDf <- tibble(
        id = x$id,
        group = group_name,
        type = x$type,
        bearing = bearing,
        refLocation_dist2col = x$dist2col
    )

    if (exists("y")) {
        y1 <- y %>%
            group_by(id) %>%
            summarise(
                startdate = min(date) %>% as_date(),
                enddate = max(date) %>% as_date(),
                trip_dur = difftime(enddate, startdate)
            )
    }

    out <- left_join(newDf, y1) %>%
        select(-geometry) %>%
        # round to 2 decimal places
        mutate(across(where(is.numeric), round, 2))

    return(out)
}

generateSummaries <- function(loc_early, locf_early, pt_early, circ.w, circ.f, circ.p, locw, locf_sf, pt_sf, output_path, bearing_ref_day) {
    indiv_summ <- prepIndividualDf(loc_early, "weaner", circ.w, locw) %>%
        bind_rows(prepIndividualDf(locf_early, "adult female", circ.f, locf_sf %>% rename(date = gmt))) %>%
        bind_rows(prepIndividualDf(pt_early, "particle trace", circ.p, pt_sf))

    indiv_summ <- indiv_summ %>% mutate(compass_zone = sapply(bearing, whichZone))

    group_summ <- indiv_summ %>%
        group_by(group) %>%
        summarise(
            n = n(),
            bearing_mean = mean(bearing),
            bearing_se = sd(bearing) / sqrt(n),
            bearing_min = min(bearing),
            bearing_max = max(bearing),
            refLocation_dist2col_mean = mean(refLocation_dist2col),
            refLocation_dist2col_se = sd(refLocation_dist2col) / sqrt(n),
            refLocation_dist2col_min = min(refLocation_dist2col),
            refLocation_dist2col_max = max(refLocation_dist2col),
            tripDur_mean = mean(trip_dur),
            tripDur_se = sd(trip_dur) / sqrt(n),
            tripDur_min = min(trip_dur),
            tripDur_max = max(trip_dur)
        ) %>%
        # round to 2 decimal places
        mutate(across(where(is.numeric), round, 2))

    write_csv(group_summ, paste0(output_path, "group_summary-bearing_ref_day=", bearing_ref_day, ".csv"))
    write_csv(indiv_summ, paste0(output_path, "allBearings-bearing_ref_day=", bearing_ref_day, ".csv"))

    return(list(individual_summary = indiv_summ, group_summary = group_summ))
}
