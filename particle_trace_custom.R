# xy = in longlat +proj=longlat +datum=WGS84, resolution in degrees, lon in -180/180
# curr_u and curr_v are in longlat +proj=longlat +datum=WGS84, resolution in degrees, lon in 0-360
# output should be in longlat +proj=longlat +datum=WGS84, resolution in degrees, lon in -180/180

# # Helper function to convert velocity from m/s to degrees/second
# velocity_to_degree <- function(uv, lat, earth_radius) {
#     # Convert u (longitude) velocity
#     uv[, 1] <- uv[, 1] / (earth_radius * cos(lat * pi / 180)) * (180 / pi)

#     # Convert v (latitude) velocity
#     uv[, 2] <- uv[, 2] / earth_radius * (180 / pi)

#     return(uv)
# }

xyuv2l <- function(xy, x) {
    xscale <- 1 / cos(xy[, 2] * pi / 180)
    cbind(x[, 1] * xscale, x[, 2]) / (1852 * 60)
}

particle_trace_custom <- function(
    xy, time_step = 24 * 3600, start_date = NULL, end_date = NULL,
    method = "bilinear", plot = FALSE, silent = FALSE, rk = FALSE, curr_u, curr_v) {
    if (time_step < 0) {
        stop("'time_step' must be positive (use 'start_date' and 'end_date' to specify direction)")
    }
    if (is.null(start_date)) {
        stop("must specify start_date e.g. as.Date(\"2021-02-15\")")
    }
    start_date <- as.POSIXct(start_date, tz = "UTC")
    if (is.null(end_date)) {
        stop("must specify end_date e.g. as.Date(\"2021-02-15\")")
    }
    end_date <- as.POSIXct(end_date, tz = "UTC")
    if (end_date == start_date) {
        stop("'end_date' must not be equal to 'start_date'")
    }
    if (end_date < start_date) {
        time_step <- -time_step
    }

    dates <- seq(start_date, end_date, by = time_step)

    N <- length(dates)

    if (!silent) {
        message(sprintf(
            "proceeding with %s time step from %s to %s with %i increments",
            c("positive", "negative")[(sign(time_step) < 0) +
                1], as.character(start_date), as.character(end_date),
            N
        ))
    }

    l <- vector("list", N)

    if (!silent) {
        pb <- progress::progress_bar$new(total = N)
    }

    # Get the right date from uo and vo
    curr_dates <- time(curr_u)

    cat("before xy", xy, "\n")
    # Convert input longitudes from -180/180 to 0-360 range
    xy[, 1] <- (xy[, 1] + 360) %% 360
    cat("after xy", xy, "\n")

    # Earth's radius in meters
    earth_radius <- 6371000

    for (jj in seq_along(dates)) {
        model_time <- dates[jj] # eg. "1995-12-15 UTC"
        # model_time %>% class() # POSIXct

        # Which date is closest to the model_time
        curr_date <- curr_dates[which.min(abs(curr_dates - model_time))]

        # Get the index of the date
        curr_idx <- which(curr_dates == curr_date)

        # Get the current
        uo_curr <- curr_u[[curr_idx]]
        vo_curr <- curr_v[[curr_idx]]

        # Combine the current
        curr <- c(uo_curr, vo_curr) %>% raster::brick()

        uv0 <- raster::extract(curr, xy, method = method)

        # Handle NAs
        valid_points <- complete.cases(uv0)
        if (sum(!valid_points) > 0) {
            warning(sprintf("%d points had NA velocities and were removed", sum(!valid_points)))
            xy <- xy[valid_points, ]
            uv0 <- uv0[valid_points, ]
        }

        if (rk) {
            # 4th-order Runge-Kutta method
            k1 <- xyuv2l(xy, uv0)
            k2 <- xyuv2l(xy, raster::extract(curr, xy + (k1 * time_step) / 2, method = method))
            k3 <- xyuv2l(xy, raster::extract(curr, xy + (k2 * time_step) / 2, method = method))
            k4 <- xyuv2l(xy, raster::extract(curr, xy + (k3 * time_step), method = method))

            uv <- (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)
        } else {
            # Simple Euler method
            uv <- xyuv2l(xy, uv0)
        }

        pti <- xy + uv * time_step

        cat("pti", pti, "\n")

        # Ensure longitudes stay in 0-360 range
        pti[, 1] <- pti[, 1] %% 360
        cat("after pti", pti, "\n")

        xy <- pti
        l[[jj]] <- xy

        if (!silent) {
            pb$tick()
        }
        if (plot) {
            points(xy, pch = ".")
        }
    }
    pts <- setNames(tibble::as_tibble(do.call(rbind, l)), c(
        "x",
        "y"
    ))
    pts$group <- rep(1:nrow(l[[1]]), length.out = nrow(pts))
    # Remove or adjust the final reprojection
    # pts[c("lon", "lat")] <- reproj::reproj_xy(cbind(pts$x, pts$y),
    #     "+proj=longlat +datum=WGS84",
    #     source = projection(curr)
    # )
    pts$date <- rep(dates, lengths(l) / ncol(l[[1]]))
    # Convert final longitudes back to -180/180 range
    pts$x <- ((pts$x + 180) %% 360) - 180
    pts
}
