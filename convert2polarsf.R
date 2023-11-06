## function for converting dataset for plotting with polar projection
#3 requires lon and lat columns in that col name.
## returns an sf object that is transformed with polar projection

convert2polarsf <- function(data, 
         crs = "+proj=stere +lon_0=170 +lat_0=-90 +units=km +datum=WGS84", 
         proj = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"){
  
  require(raster)
  require(sf)

  coordinates(data) <- c("lon", "lat")
  projection(data) <-  proj
  data <- data %>% st_as_sf(data) %>% st_transform(crs = crs)
  
  return(data)
  
}
