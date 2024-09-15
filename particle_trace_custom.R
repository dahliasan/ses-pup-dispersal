# library(conflicted)
# conflict_prefer("extract", "raster")

# xy <- pts
# time_step <- 24 * 3600
# method <- "bilinear"
# plot <- FALSE
# silent <- FALSE
# rk <- FALSE

# curr_u <- uo1
# curr_v <- vo1
library(logger)
# install.packages("logger")

particle_trace_custom <- function(
    xy, time_step = 24 * 3600, start_date = NULL, end_date = NULL,
    method = "bilinear", plot = FALSE, silent = FALSE, rk = FALSE, curr_u, curr_v) {
    # readcurr_POLAR <- memoise::memoize(raadtools:::readcurr_polar)


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

    # files <- raadfiles::altimetry_currents_polar_files()

    jj <- 1

    for (jj in seq_along(dates)) {
        model_time <- dates[jj] # eg. "1995-12-15 UTC"
        # model_time %>% class() # POSIXct

        # Get the right date from uo and vo
        curr_dates <- raster::getZ(curr_u)
        curr_dates <- lubridate::ymd(curr_dates)
        # convert to POSIXct
        curr_dates <- as.POSIXct(curr_dates, tz = "UTC")

        # Which date is closest to the model_time
        curr_date <- curr_dates[which.min(abs(curr_dates - model_time))]

        # Get the index of the date
        curr_idx <- which(curr_dates == curr_date)

        # Get the current
        uo_curr <- curr_u[[curr_idx]]
        vo_curr <- curr_v[[curr_idx]]

        # Combine the current
        curr <- raster::brick(uo_curr, vo_curr)
        # curr <- readcurr_POLAR(model_time, inputfiles = files)

        if (jj == 1L) {
            xy <- reproj::reproj_xy(xy, raster::projection(curr),
                source = "+proj=longlat +datum=WGS84"
            )
            if (plot) {
                plot(xy, asp = 1)
            }
        }
        uv0 <- raster::extract(curr, xy, method = method)
        if (rk) {
            uv1 <- raster::extract(curr, xy + (uv0 * time_step) / 2, method = method)
            uv2 <- raster::extract(curr, xy + (uv1 * time_step) / 2, method = method)
            uv3 <- raster::extract(curr, xy + (uv2 * time_step), method = method)
            pti <- xy + (uv0 / 6 + uv1 / 3 + uv2 / 3 + uv3 / 6) * time_step
        } else {
            pti <- xy + uv0 * time_step
        }
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
    pts[c("lon", "lat")] <- reproj::reproj_xy(cbind(pts$x, pts$y),
        "+proj=longlat +datum=WGS84",
        source = projection(curr)
    )
    pts$date <- rep(dates, lengths(l) / ncol(l[[1]]))
    pts
}
