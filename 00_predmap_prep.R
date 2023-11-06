## Determining which months to use as input data for prediction map

library(tidyverse)
library(lubridate)
d1 <- readRDS("./output/all_data_combined.rds")
names(d1)
tmp <- d1 %>% 
  filter(trip == 1) %>% 
  group_by(id) %>% 
  summarise(startdate = min(date), 
            enddate = max(date), 
            tripdur = mean(tripdur),
            furthest_loc_date = date[ which.max(dist2col)],
            furthest_month = month(furthest_loc_date))
tmp$furthest_month %>% unique
