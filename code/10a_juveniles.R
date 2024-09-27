# Prepare and take a look at juvenile SES tracks sent by Mark
rm(list = ls())

library(tidyverse)
library(foieGras)
source('convert2polarsf.R')

# Load data ---------------------------------------------------------------
d <- read_csv('./Data/MI Juvenile Seals data.csv') %>% 
  rename(lon = Long, lat = Lat, id = SealID)
names(d) <- tolower(names(d))
d <- d %>% mutate(gmt = lubridate::dmy_hms(paste(date, time)))
d$lat <- d$lat * -1




# Fit ssm -----------------------------------------------------------------
d1 <- d %>% 
  mutate(lc = 'GL', lonerr = 1, laterr = 1) %>% 
  dplyr::select(id, gmt, lc, lon, lat, lonerr, laterr) %>% 
  rename(date = gmt)

fit_all <- foieGras::fit_ssm(d1,
                             vmax = 4, 
                             map = list(psi = factor(NA)),
                             model = "crw",
                             time.step = 12)
ssm <- grab(fit_all, "predicted", as_sf=FALSE)
saveRDS(ssm, file = './output/foiegras_fitssm_juveniles.rds')

# Plots -------------------------------------------------------------------

# Basic summary of data
d$id %>% unique

# Visualise tracks
d_sf <- convert2polarsf(d)
d_sf %>% ggplot() +
  geom_sf(size = 0.1, aes(col = id)) + 
  theme(legend.position = "none")

ssm_sf <- convert2polarsf(ssm)
ssm_sf %>% ggplot() +
  geom_sf(size = 0.1, aes(col = id)) + 
  theme(legend.position = "none")


  ssm %>% 
  ggplot(aes(x = date, y = id)) + 
  geom_point(size = .5)


# Assign year born --------------------------------------------------------
# SES Tag identification year born
# F & H - 1994
# J & L - 1995
# K - 1996
# D - 1997
# P -  1998
# T - 1999


year_born <- function(x) {
  letter <- substr(x, 1, 1) %>% tolower()
  if(letter == 'f' | letter == 'h') return(1994) 
  else if(letter == 'j' | letter == 'l') return(1995)
  else if(letter == 'k') return(1996)
  else if(letter == 'd') return(1997)
  else if(letter == 'p') return(1998)
  else if(letter == 't') return(1999)
  else return(NA)
}

ssm$year_born <- lapply(ssm$id, year_born) %>% unlist()
ssm <- ssm %>% mutate(age = lubridate::year(date) - year_born)
ssm$type <- 'gls'
out <- ssm %>% dplyr::select(id, date, lon, lat, year_born, age, type)

saveRDS(out, file = './output/macca_juvenile_locs.rds')
