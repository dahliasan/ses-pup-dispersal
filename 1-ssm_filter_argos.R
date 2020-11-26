library(Hmisc)
library(tidyverse)
library(lubridate)


## Install foieGras
# remotes::install_github("ianjonsen/foieGras")
library(foieGras)



# Read MDB files (Microsoft access, ew.) ----------------------------------

mdb <-  mdb.get('./Data/AccessDB/mq1.mdb')
mdb2 <-  mdb.get('./Data/AccessDB/mq2.mdb')

df <- mdb$diag %>% as_tibble() %>% select(-PTT, -ACTUAL.PTT)
df2 <- mdb2$diag %>% as_tibble() %>% select(-PTT, -ACTUAL.PTT)
df <- bind_rows(df, df2)

write.csv(df, './Data/mq_ellie_weaners_argos.csv')


# Tidy up data ------------------------------------------------------------

## Parse Date Time 
df <- df %>% 
  mutate(D.DATE = D.DATE %>% as.character() %>% mdy_hms())

## Get rid of the "labelled" class attribute for numeric columns
df <- df %>% mutate_if(is.numeric, as.numeric)


## Select only relevant columns and rename them
weaners <- df %>% 
  select(ref, D.DATE, LQ, LON, LAT)

colnames(weaners) <- c('id', 'date', 'lc', 'lon', 'lat')

weaners$id <- weaners$id %>% as.character()


# Fit SSM -----------------------------------------------------------------
weaners <- weaners %>% filter(!id == 'mq1-22488-95') # this seal had < 3 locations

fit <- fit_ssm(weaners, model= 'crw', vmax = 2, time.step = NA, verbose = 0, map = list(psi = factor(NA)))