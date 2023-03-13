### Match seal IDs ###
rm(list = ls())

library(tidyverse)
library(lubridate)


w1 <- read_csv("./Data/mq-ellie-weaners-argos.csv")
w2 <- w1 %>% select(ref, ptt) %>% distinct()

## identify seals that left colony
# keep = if seal left colony or not (determined from visual checking)
id_keep <- dir('./Output/ssm plots/keep') %>% str_replace('.png', '')
w2$keep <- 0
w2$keep[w2$ref %in% id_keep] <- 1
w2 <- w2 %>% arrange(desc(keep))

## add year column
# mq1 = 1995, mq2 = 1996, mq3 = 1999, mq4 = 2000
w2$cohort <- substr(w2$ref, 1, 3)
w2$year <- 1995
w2$year[w2$cohort == 'mq2'] <- 1996
w2$year[w2$cohort == 'mq3'] <- 1999
w2$year[w2$cohort == 'mq4'] <- 2000


# IDs ---------------------------------------------------------------------

# + Match to brand + flipper tag id ---------------------------------------
# NOTE: some seals have same PTT number but are different individuals! Maybe PTT reused?
tag1 <- read_csv('./Data/TagID/mqPTT_95_96.csv') %>% 
  mutate(Tag1 = Tag1 %>% as.character(), Tag2 = Tag2 %>% as.character())
tag2 <- read_csv('./Data/TagID/PTT_seals_1999_2000.csv') %>% 
  mutate(PTT = PTT %>% as.character())


## Convert tag1 (1995-96) seal code to ptt number
tag1 <- tag1 %>% 
  mutate(PTTx = str_split(SRDLbody, '-'))

tag1$PTT <- purrr::map(tag1$PTTx, function(x){
  if(!is.na(x)) x[2] else NA
}) %>% unlist()

tag1$Year <- purrr::map(tag1$PTTx, function(x){
  if(!is.na(x)) x[3] else NA
}) %>% unlist()

tag1$Year <- strptime(tag1$Year, '%y') %>% year()

tag <- bind_rows(tag1 %>% select(-PTTx), tag2)
tag <- tag %>% 
  select(Year, PTT, Brand, Tag1, Tag2)
colnames(tag) <- tolower(colnames(tag))
tag_final <- tag %>% filter(!is.na(ptt))

## check for duplications
# tag_final <- tag_final %>% mutate(dups_ptt = duplicated(ptt))
# tag_final %>% filter(dups_ptt == TRUE)
tag_final %>% select(year, ptt) %>% duplicated() %>% sum() # ptts may have be reused but on different seals

## combine ptt id to brand + flipper tag id
w2 <- w2 %>% mutate(year_ptt = paste(year, ptt, sep = '_'))
tag_final <- tag_final %>% mutate(year_ptt = paste(year, ptt, sep = '_')) %>% 
  mutate(ptt = ptt %>% as.double())
w2 <- left_join(w2, tag_final)


# + Add in missing but now discovered seal ids ----------------------------
# Feb 2022: these were missing, but then discovered from Martin Biuw and George.

missing <- w2 %>% filter(year == '1995', is.na(brand), is.na(tag1), is.na(tag2)) %>% 
  select(-brand, -tag1, -tag2)

wean95 <- read_csv('./Data/Satellite tagged weaners 1995 full data.csv') %>% 
  mutate(tag1 = tag1 %>% as.character, tag2 = tag2 %>% as.character()) %>% 
  mutate(year_ptt = paste(1995, ptt, sep = '_'),
                            brand = paste('J', brand, sep ='')) %>% 
  select(ptt, brand, tag1, tag2, year_ptt)

missing <- left_join(missing, wean95)

w2 <- w2 %>% filter(!year_ptt %in% missing$year_ptt)
w2 <- bind_rows(w2, missing) %>% arrange(cohort)


## Fix 1995 J91 and J26 weaners brand ids 
w2$brand[w2$brand == 'J91'] <- 'J091'
w2$brand[w2$brand == 'J26'] <- 'J026'

# + Match to database seal ids --------------------------------------------
## Read in database seal ids
sealID <- read_csv('./Data/CMC/TAGS/TAGS_ELLIE_TAGS.csv', col_types = 'dddcccc')
sealID <- sealID %>% dplyr::select(SEAL_ID, SEAL_TAG, SEAL_TAG_TYPE) %>% unique()

tmp <- w2

## Join by brand ID first 
tmp <- left_join(tmp, sealID %>% filter(SEAL_TAG_TYPE == 'B'), by = c('brand' = 'SEAL_TAG'))

# get only remaining with no SEAL_ID yet then join by tag
tmp1 <- tmp %>% filter(is.na(SEAL_ID)) %>% select(-SEAL_ID, -SEAL_TAG_TYPE) %>%  
  left_join(sealID %>% filter(SEAL_TAG_TYPE == 'T'), by = c('tag1' = 'SEAL_TAG'))

w3 <- bind_rows(tmp %>% filter(!is.na(SEAL_ID)), tmp1)


## Save as RDATA
sealID_all <- w3 %>% dplyr::select(ref, ptt, brand, tag1, tag2, SEAL_ID)
save(sealID_all, file = './Output/seal_ID.RData')
