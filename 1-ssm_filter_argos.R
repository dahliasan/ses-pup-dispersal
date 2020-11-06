library(Hmisc)
library(tidyverse)
library(lubridate)
library(trip)


## Install foieGras
# remotes::install_github("ianjonsen/foieGras")
library(foieGras)



# Read MDB files (Microsoft access, ew.) ----------------------------------

mdb <-  mdb.get('./Data/AccessDB/mq1.mdb')
mdb2 <-  mdb.get('./Data/AccessDB/mq2.mdb')

df <- mdb$diag %>% as_tibble() %>% select(-PTT) %>% mutate(ACTUAL.PTT = ACTUAL.PTT %>% as.character)
df2 <- mdb2$diag %>% as_tibble() %>% select(-PTT) %>% mutate(ACTUAL.PTT = ACTUAL.PTT %>% as.character)
df <- bind_rows(df, df2)

# write_csv(df, './Data/mq-ellie-weaners-argos.csv')

# Tidy up data ------------------------------------------------------------

## Parse Date Time 
df <- df %>% 
  mutate(D.DATE = D.DATE %>% as.character() %>% mdy_hms())

## Get rid of the "labelled" class attribute for numeric columns
df <- df %>% mutate_if(is.numeric, as.numeric)


## Select only relevant columns and rename them
weaners <- df %>% 
  select(ACTUAL.PTT, D.DATE, LQ, LON, LAT)

colnames(weaners) <- c('id', 'date', 'lc', 'lon', 'lat')

weaners$id <- weaners$id %>% as.character()
weaners$lc <- weaners$lc %>% factor()


# Fit SSM -----------------------------------------------------------------
library(argosfilter)
w <- weaners %>% 
  group_by(id) %>% 
  nest()

w$data <- w %>% 
  .$data %>% 
  map(function(x) {
    x$filter <- with(x, vmask(lat, lon, date, 4))
    x %>% filter(filter == 'not')
  })

w <- w %>% 
  unnest(cols = c(data)) %>% 
  select(-filter) %>% 
  ungroup()

w <- w %>% filter(!id == '17213')

fit <- fit_ssm(weaners %>% filter(id == '2849'), model= 'crw',  vmax = 5, time.step = 48, verbose = 0, 
               map = list(rho_o = factor(NA)))


# Plot raw locations ------------------------------------------------------


w %>% 
  ggplot(aes(lon, lat, colour = id)) + 
  geom_point() + 
  geom_path()
