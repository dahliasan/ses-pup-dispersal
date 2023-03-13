## ARS plot from foiegras::mpm result


rm(list = ls())
library(tidyverse)
library(raster)
library(sf)
library(viridis)
library(orsifronts)
load("baseInfo.Rdata")
source('convert2polarsf.R')



# Load data ---------------------------------------------------------------
d <- readRDS('./Output/all_data_combined.rds')

## Get only real tracks (remove pseudo locations)
d <- d %>% filter(sim == 0, land == FALSE, trip == 1)



# Convert dataframe to raster ---------------------------------------------
# https://irapoenya.wordpress.com/2018/10/11/r-studio-create-raster-by-rasterize-function/

df <- d %>% ungroup() %>% dplyr::select(lon, lat, g) %>% 
  rename(x = lon, y = lat) %>% 
  filter(g < 0.5)


#create the extent and the raster
wgs <- CRS("+proj=longlat +ellps=WGS84")
e <- extent(df[,(1:2)])
r <- raster(e, resolution = 0.5, crs= wgs)
r_new <- rasterize(df[,1:2], r, df[,3], fun=mean)
r_polar <- projectRaster(r_new, crs = crs)
r_polar <- trim(r_polar)

r <- r_polar

# Plot raster
# set up orsifronts
ofp <- spTransform(orsifronts, CRS(projection(crs)))
ofp <- crop(ofp, extent(r))

plot(r, col = viridis(125, direction = -1))
plot(ofp, add=T, col="black")



df_new <- as.data.frame(r_polar, xy=TRUE) %>% 
  tibble() %>% 
  filter(!is.na(layer)) %>% 
  rename(g = layer)


# Plot individual ARS -----------------------------------------------------

df <- d %>% ungroup() %>% dplyr::select(lon, lat, g, id) %>% 
  filter(g < 0.5)


df %>% 
  convert2polarsf() %>% 
  ggplot() +
  geom_sf(data = convert2polarsf(d), color = 'black', alpha = 0.1, size = 0.1) + 
  geom_sf(aes(colour = g), size = 0.2) +
  facet_wrap(~id) + 
  theme_bw() + 
  scale_colour_viridis(direction = -1)

# Smoothr -----------------------------------------------------------------

library(smoothr)
library(units)
r <- r_polar < 0.25
plot(rasterToPolygons(r), col = NA, border = NA) # set up plot extent
plot(r, col = c("white", "#4DAF4A"), legend = FALSE, add = TRUE, box = FALSE)

# polygonize
r_poly <- rasterToPolygons(r, function(x){x == 1}, dissolve = TRUE) %>% 
  st_as_sf()

# drop cells
r_poly_dropped <- drop_crumbs(r_poly, set_units(2086.38, km^2))
# plot
plot(rasterToPolygons(r), col = NA, border = NA) # set up plot extent
plot(r_poly_dropped, col = "#4DAF4A", border = "grey20", lwd = 1.5, add = TRUE)

# fill holes
r_poly_filled <- fill_holes(r_poly_dropped, set_units(2086.38 * 10, km^2))
# plot
plot(rasterToPolygons(r), col = NA, border = NA) # set up plot extent
plot(r_poly_filled, col = "#4DAF4A", border = "grey20", lwd = 1.5, add = TRUE)

# smooth
r_poly_smooth <- smooth(r_poly_filled, method = "ksmooth", smoothness = 2)
# plot
plot(rasterToPolygons(r), col = NA, border = NA) # set up plot extent
plot(r_poly_smooth, col = "#4DAF4A", border = "grey20", lwd = 1.5, add = TRUE)

# Plot
bb <- extent(r_polar)

convert2polarsf(d) %>% 
  ggplot() +
  geom_sf(size = 0.1, alpha = 0.2, color = 'grey') + 
  geom_sf(data = orsi_sf, linetype = 'dashed', colour = 'grey10') +
  geom_sf(data = r_poly_smooth, fill= 'red') + 
  xlim(bb[1], bb[2]) + ylim(bb[3]-200, bb[4]+500) + 
  theme_bw()
  

df_new %>% 
  st_as_sf(coords = c('x', 'y'), crs = crs) %>% 
  ggplot() + 
  geom_sf(aes(color = g)) + 
  viridis::scale_color_viridis() +
  geom_sf(data = orsi_sf, linetype = 'dashed', colour = 'grey10') +
  xlim(bb[1], bb[2]) + ylim(bb[3]-200, bb[4]+500) +
  theme_bw()


tm_shape(r_polar) +
  tm_raster(palette = "viridis")
