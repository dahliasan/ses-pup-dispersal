particle_trace_custom <- function(xy,
                                  uo, vo,
                                  time_step = 24 * 3600,
                                  start_date = NULL,
                                  end_date = NULL,
                                  method = "bilinear",
                                  plot = FALSE, 
                                  silent = FALSE, 
                                  rk = FALSE,
                                  lon180 = TRUE) {
  
  source("mds_rotate.R")
  
  # Validate the input dates
  if (time_step < 0) stop("'time_step' must be positive (use 'start_date' and 'end_date' to specify direction)")
  if (is.null(start_date)) stop("must specify start_date e.g. as.Date('2021-02-15')")
  start_date <- as.POSIXct(start_date, tz = "UTC")
  if (is.null(end_date)) stop("must specify end_date e.g. as.Date('2021-02-15')")
  end_date <- as.POSIXct(end_date, tz = "UTC")
  if (end_date == start_date) stop("'end_date' must not be equal to 'start_date'")
  if (end_date < start_date) time_step <- -time_step
  
  # Generate sequence of dates
  dates <- seq(start_date, end_date, by = time_step) 
  N <- length(dates)
  
  # Progress message
  if (!silent) {
    message(sprintf("Proceeding with %s time step from %s to %s with %i increments",
                    ifelse(sign(time_step) < 0, "negative", "positive"),
                    as.character(start_date), as.character(end_date), N))
  }
  
  # Initialize list to store particle positions
  l <- vector("list", N)
  
  # Initialize progress bar if not silent
  if (!silent) pb <- progress::progress_bar$new(total = N)
  
  # Loop through the time steps
  for (jj in seq_along(dates)) {
    model_time <- dates[jj]
    message(model_time)
    
    # Get the right date layer (z) from the uo and vo datasets.
    z <- which.min(abs(as.numeric(difftime(model_time, getZ(uo), units = "days"))))
    
    if (lon180) {
      uo1 <- .rotate(uo[[z]])
      vo1 <- .rotate(vo[[z]])
    } else {
      uo1 <- uo[[z]]
      vo1 <- vo[[z]]
    }
    
    target_crs <- "+proj=laea +lat_0=-90 +lon_0=0 +datum=WGS84 +units=m +no_defs"
    
    uo1 <- projectRaster(uo1, crs = target_crs, res = 25000)
    vo1 <- projectRaster(vo1, crs = target_crs, res = 25000)

    uo1
    vo1
    
    # Ensure the order is u-component first (uo), v-component second (vo)
    curr <- brick(uo1, vo1)
    
    # Convert the starting points to the CRS of the current dataset
    if (jj == 1L) {
      xy <- reproj::reproj_xy(xy, raster::projection(curr),
                              source = "+proj=longlat +datum=WGS84")
      message(xy)
      if (plot) graphics::plot(xy, asp = 1)
    }

    # Extract current velocity vectors (with method specified)
    uv0 <- raster::extract(curr, xy, method = method)
    message(uv0)
    
    
    
    # Apply Runge-Kutta method if required
    if (rk) {
      uv1 <- raster::extract(curr, xy + (uv0 * time_step)/2)
      uv2 <- raster::extract(curr, xy + (uv1 * time_step)/2)
      uv3 <- raster::extract(curr, xy + (uv2 * time_step))
      
      pti <- xy + (uv0/6 + uv1/3 + uv2/3 + uv3/6) * time_step
    } else {
      pti <- xy + uv0 * time_step
    }
  
    # Update particle positions
    xy <- pti
    l[[jj]] <- xy
    
    if (!silent) pb$tick()
    if (plot) graphics::points(xy, pch = ".")
  }
  
  # Prepare output as tibble with projected lon/lat and timestamps
  pts <- stats::setNames(tibble::as_tibble(do.call(rbind, l)), c("x", "y"))
  pts$group <- rep(1:nrow(l[[1]]), length.out = nrow(pts))
  
  message(raster::projection(curr))
  pts[c("lon", "lat")] <- reproj::reproj_xy(cbind(pts$x, pts$y), "+proj=longlat +datum=WGS84", source = raster::projection(curr))
  pts$date <- rep(dates, lengths(l)/ncol(l[[1]]))
  
  return(pts)
}
