library(raster)
library(dplyr)
library(nabor)
library(tabularaster)
library(reproj)

# Define the custom target grid with Lambert Azimuthal Equal Area projection
default_grid <- function() {
    prjj <- "+proj=laea +lat_0=-90 +datum=WGS84"
    raster(
        spex::buffer_extent(projectExtent(
            raster(extent(-180, 180, -90, -30),
                crs = "+init=epsg:4326"
            ),
            prjj
        ), 25000),
        res = 25000, crs = prjj
    )
}

# Assuming you've already loaded the Copernicus data as `uo` and `vo`
# Example:
# vo <- brick(file, varname = "vgos")  # Meridional (north-south) component
# uo <- brick(file, varname = "ugos")  # Zonal (east-west) component

# Load the target grid for reprojected output
target <- default_grid()

# Extract the coordinates of the target grid
xy <- coordinates(target)

# Reproject the coordinates from target grid (LAEA) to original (lon-lat) using reproj
xyi <- reproj::reproj(xy, target = "OGC:CRS84", source = projection(target))

# Initialize the index for the target grid
index <- tibble(cellindex = seq_len(ncell(uo)))

# Project the source (Copernicus) grid coordinates to the target CRS and store for nearest-neighbor
source_crs <- projection(uo) # Assuming the source data is in Equidistant Cylindrical projection
index[c("qx", "qy")] <- reproj::reproj(xyFromCell(uo, index$cellindex), source = source_crs, target = "OGC:CRS84")
index[c("sx", "sy")] <- reproj::reproj(xyFromCell(uo, index$cellindex), source = source_crs, target = projection(target))

# KNN for nearest-neighbor search
knn <- WKNNF(as.matrix(index[c("sx", "sy")]))

# Initialize empty lists to store processed rasters
uo_processed <- list()
vo_processed <- list()

# Loop through each layer of the `uo` and `vo` rasters (representing different time steps)
for (i in seq_len(nlayers(uo))) {
    # Extract the u (east-west) and v (north-south) components for the current time step
    U <- raster::subset(uo, i)
    V <- raster::subset(vo, i)

    # Compute the new x and y coordinates after adding the u and v components
    index$x1 <- index$qx + values(U)
    index$y1 <- index$qy + values(V)

    # Reproject the calculated coordinates back to the target CRS
    index[c("ex", "ey")] <- reproj::reproj(as.matrix(index[c("x1", "y1")]), target = projection(target), source = source_crs)

    # Calculate the shifts in the x and y directions (for proper alignment of currents)
    index$pu <- index$ex - index$sx
    index$pv <- index$ey - index$sy

    # Extract valid values
    ee <- raster::extract(U, xyi)

    # Handle non-NA values for re-projection
    xyq <- xy[!is.na(c(ee)), ]
    idx <- knn$query(xyq, k = 1, eps = 0, radius = 0)

    # Create rasters for the reprojected u and v components
    uu <- vv <- target
    uu[!is.na(ee)] <- index$pu[idx$nn.idx]
    vv[!is.na(ee)] <- index$pv[idx$nn.idx]

    # Set time information in the output rasters
    uu <- setZ(uu, names(uo)[i])
    vv <- setZ(vv, names(vo)[i])

    # Store the processed layers in memory
    uo_processed[[i]] <- uu
    vo_processed[[i]] <- vv

    # Print progress
    print(paste("Processed layer", i, "out of", nlayers(uo)))
}

# Combine the processed layers into RasterBricks
uo_final <- brick(uo_processed)
vo_final <- brick(vo_processed)

# Return the final processed bricks
uo_final
vo_final
plot(uo_final[[1]])
