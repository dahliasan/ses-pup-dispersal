# Create and export workspace (.RData) with all essential data

rm(list = ls())
options(tibble.width = Inf)
library(tidyverse)
library(sf)
library(raster)
library(viridis)
library(tmap)
library(ggpubr)
library(see)
library(orsifronts)
library(pals)
source('convert2polarsf.R')
load('./output/seal_id.rdata')


# Seal IDs to keep --------------------------------------------------------

id_keep <- dir('./Output/ssm plots/keep') %>% str_replace('.png', '')

id_keep <- id_keep[!id_keep %in% c("mq3-26635-99", "mq3-2841-99", "mq3-2845-99")] # Non-weaner seals

# missing_ids <- c('mq2-22483-96',
#                  'mq2-22490-96',
#                  'mq2-26624-96')

seals_with_ids <- sealID_all %>% filter(!is.na(SEAL_ID)) %>% pull(ref)
id_keep <- id_keep[id_keep %in% seals_with_ids]

# Set common graphical elements ------------------------------------------
crs <- "+proj=stere +lon_0=170 +lat_0=-90 +units=km +datum=WGS84"
proj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# world map
world_sf <- st_as_sf(rworldmap::getMap(resolution = "li")) %>% st_transform(., crs = crs)

# mq island
mq <- data.frame(lon = 158.95, lat = -54.5)
mq_sf <- convert2polarsf(mq)

# orsi fronts
orsi_sf <- st_as_sf(orsifronts) %>% 
  filter(front %in% c('stf', 'saf', 'pf', 'saccf')) %>%  
  st_transform(., crs = crs) %>% 
  mutate()

save(id_keep, crs, proj, world_sf, mq, mq_sf, orsi_sf, file = 'baseInfo.RData')
