## Combine all variables / data together to create master dataset for modelling etc.

rm(list = ls())
options(tibble.width = Inf)
library(tidyverse)
library(lubridate)
library(foieGras)

# Load Datasets -----------------------------------------------------------
load('./Output/allFit_prefilter=0.3.RData')
load("./Output/raadtools_extract_sim_tracks_12h_chla.rdata")
load("./Output/survival.RData")
source('convert2polarsf.R')


# 1) Calculate delta drift rate ----------------------------------------------
## Extract drift dives from augFit
drift <- purrr::map(allFit, function(x) {
  xx <- x$augFit
  if(xx %>% length > 1) {
    data <- xx$data # original timestamps
    pred <- xx$pred # daily predicted timestamps
    pred$id <- data$ref[1]
    pred %>% 
      mutate(date = DE.DATE %>% as.Date) %>% 
      dplyr::select(id, date, fit, sd.fit, lwr, upr, segment) %>% 
      as_tibble()
  }
}) %>% reduce(bind_rows)

## Calculate daily drift rate change
drift <- drift %>% 
  filter(!duplicated(drift)) %>% # remove duplicate rows
  group_by(id) %>% 
  mutate(drate_chg = lead(fit) - fit) %>% 
  rename(drate = fit) %>% 
  ungroup()

## Check duplication - all ok
# drift <- drift %>% 
#   group_by(id) %>% 
#   mutate(dup = duplicated(date))
# 
# dup <- drift %>% filter(dup == TRUE)

## Combine location and drift rate data
d1 <- d %>% ungroup() %>% 
  mutate(day = as.Date(date))
drift <- drift %>% rename(day = date)
d1 <- left_join(d1, drift %>% dplyr::select(id, day, drate, drate_chg)) 
d1 <- d1 %>% filter(!duplicated(d1))
dim(d1)
dim(d)

# Checking things...
# drift %>% filter(id == 'mq2-20916-96') %>% ggplot(aes(date, drate_delta)) + geom_point() +
# d %>% filter(id == 'mq2-20916-96') %>% ggplot(aes(day, drate_delta)) + geom_point()
# d %>% filter(drate_delta < 0.08) %>% ggplot(aes(dist2col, drate_delta)) + geom_point() +facet_wrap(~id)

# save(d, file = './Output/driftRate_locations_1d.RData')

# 2) Combine all datasets -------------------------------------------------
d1 <- d1 %>% left_join(survival)

d1 <- d1 %>% 
  group_by(id, trip, sim) %>% 
  arrange(id, trip, date, sim) %>% 
  mutate(daysFromDepart =  difftime(date, min(date), units = 'days') %>% as.numeric())

# saveRDS(d1, file = './Output/all_data_combined.rds')


# 3) add/combine mpm model + bearing info ----------------------------------
# last edit: 20 Dec 2021
d1 <- readRDS('./Output/all_data_combined.rds')
load('./output/foiegras_fitmpm_vmax4_12h.Rdata')
bear <- readRDS('./output/dispersal_bearing___REF_DATE=6d.rds')

fmp <- grab(fmp_all, "fitted")
d11 <- d1 %>% filter(sim == 0)
d11 <- d11 %>% left_join(fmp) %>% left_join(bear)

out <- d11 %>% bind_rows(d1 %>% filter(sim != 0))

# saveRDS(out, file = './Output/all_data_combined.rds')



# 4) Fix birth year column ------------------------------------------------

d1 <- readRDS('./Output/all_data_combined.rds')


d1$fieldseason <- substr(d1$id, 1, 3)
d1$birthyear[d1$fieldseason == 'mq1'] <- 1995
d1$birthyear[d1$fieldseason == 'mq2'] <- 1996
d1$birthyear[d1$fieldseason == 'mq3'] <- 1999
d1$birthyear[d1$fieldseason == 'mq4'] <- 2000

# saveRDS(d1, file = './Output/all_data_combined.rds')
