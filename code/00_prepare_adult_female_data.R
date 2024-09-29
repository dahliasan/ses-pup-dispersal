library(tidyverse)
load("data/macca_winter_locs.Rdata") # adult females
loc <- loc %>%
    as_tibble() %>%
    rename(id = seal, date = gmt)
saveRDS(loc, "data/adult_female_locs.rds")
