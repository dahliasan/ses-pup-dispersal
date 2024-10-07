## function for converting dataset for plotting with polar projection
# 3 requires lon and lat columns in that col name.
## returns an sf object that is transformed with polar projection

convert2polarsf <- function(data,
                            crs = "+proj=laea +lat_0=-90 +lon_0=170 +datum=WGS84 +units=km +ellps=WGS84 +no_defs",
                            proj = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                            remove_coords = TRUE) {
  require(sf)

  data <- data %>%
    st_as_sf(remove = remove_coords, coords = c("lon", "lat")) %>%
    st_set_crs(proj) %>%
    st_transform(crs = crs)

  return(data)
}
