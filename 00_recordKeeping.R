### Tracking Data Output ###

library(tidyverse)
library(lubridate)


w1 <- read_csv("./Data/mq-ellie-weaners-argos.csv")
w2 <- w1 %>% 
  group_by(ref, ptt) %>% 
  mutate(d.date = d.date %>% mdy_hms()) %>% 
  summarise(deployment.start = min(d.date), deployment.end = max(d.date))

# Calculate deployment duration
w2 <- w2 %>% mutate(deployment.dur = difftime(deployment.end, deployment.start, units = 'days') %>% round(1))


id_keep <- dir('./Output/ssm plots/keep') %>% str_replace('.png', '')

w2$left.colony <- 0
w2$left.colony[w2$ref %in% id_keep] <- 1


## Load drift dive data
load('./Output/allFit.RData')

# combFit
j <- purrr::map(allFit, function(x) {
  xx <- x$combFit
  class(xx)
})

tmp <- data.frame(ref = names(allFit), combFit = unlist(j))

# augFit
j <- purrr::map(allFit, function(x) {
  xx <- x$augFit
  class(xx)
  
})
tmp$augFit = unlist(j)

# augFit - prefilter = 0.3
load('./Output/allFit_prefilter=0.3.RData')

j <- purrr::map(allFit, function(x) {
  xx <- x$augFit
  class(xx)
  
})
tmp$augFit_pf0.3 = unlist(j)

tmp <- tmp %>% 
  mutate(combFit = ifelse(combFit == 'list', 1, 0),
         augFit = ifelse(augFit == 'list', 1, 0),
         augFit_pf0.3 = ifelse(augFit_pf0.3 == 'list', 1, 0))


w2 <- left_join(w2, tmp)

## Arrange by left colony or not
w2 <- w2 %>% arrange(left.colony)


# IDs ---------------------------------------------------------------------
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


# check for duplications
# tag_final <- tag_final %>% mutate(dups_ptt = duplicated(ptt))
# tag_final %>% filter(dups_ptt == TRUE)
tag_final %>% select(year, ptt) %>% duplicated() %>% sum() # ptts may have be reused but on different seals

## Combine with full record df
w3 <- w2 %>% mutate(year_ptt = paste(year(deployment.start), ptt, sep = '_'))
tag_final <- tag_final %>% mutate(year_ptt = paste(year, ptt, sep = '_'))
w3 <- left_join(w3 %>% mutate(ptt = ptt %>% as.character()), tag_final)

sealID <- read_csv('./Data/CMC/TAGS/TAGS_ELLIE_TAGS.csv', col_types = 'dddcccc')
sealID <- sealID %>% dplyr::select(SEAL_ID, SEAL_TAG, SEAL_TAG_TYPE) %>% unique()

wtmp <- w3 %>% dplyr::select(ref, ptt, year_ptt, brand, tag1, tag2)

# join by brand ID first 
wtmp <- left_join(wtmp, sealID %>% filter(SEAL_TAG_TYPE == 'B'), by = c('brand' = 'SEAL_TAG'))

# get only remaining with no SEAL_ID yet then join by tag
wtmp1 <- wtmp %>% filter(is.na(SEAL_ID)) %>% select(-SEAL_ID, -SEAL_TAG_TYPE) %>%  
  left_join(sealID %>% filter(SEAL_TAG_TYPE == 'T'), by = c('tag1' = 'SEAL_TAG'))

w4 <- bind_rows(wtmp %>% filter(!is.na(SEAL_ID)), wtmp1)

w5 <- left_join(w3,w4)

# write.csv(w5, './Output/ellies_analysisRecords_db.csv',  na = "-")

## Save as RDATA
sealID_all <- w5 %>% dplyr::select(ref, ptt, brand, tag1, tag2, SEAL_ID)
save(sealID_all, file = './Output/seal_ID.RData')
# 
# 
# # Try something else for ID -----------------------------------------------
# wtmp <- w3 %>% dplyr::select(ref, ptt, year_ptt, brand, tag1, tag2)
# sealID <- read_csv('./Data/CMC/TAGS/TAGS_ELLIE_TAGS.csv', col_types = 'dddcccc')
# 
# sealID <- sealID %>% filter(SEAL_TAG %in% na.omit(c(wtmp$brand, wtmp$tag1, wtmp$tag2)))
# sealID <- sealID %>% arrange(sealID)
# 
# d <- data.frame(table(sealID$SEAL_ID, sealID$SEAL_TAG))
# colnames(d) <- c("id", "tag", "freq")
# d <- d %>% filter(freq != 0) %>% 
#   arrange(id)
# 
# sealID %>% 
#   ggplot(aes(x = as.factor(SEAL_ID), y = SEAL_TAG)) + 
#   geom_point()
# 
# 
# 
# # join by brand ID first 
# sealID1 <- sealID %>% dplyr::select(SEAL_ID, SEAL_TAG, SEAL_TAG_TYPE) %>% unique()
# wtmp <- left_join(wtmp, sealID1 %>% filter(SEAL_TAG_TYPE == 'B'), by = c('brand' = 'SEAL_TAG'))
# 
# # Next join by flipper tag
# wtmp1 <- wtmp %>% filter(is.na(SEAL_ID)) %>% select(-SEAL_ID, -SEAL_TAG_TYPE) %>%  
#   left_join(sealID1 %>% filter(SEAL_TAG_TYPE == 'T'), by = c('tag1' = 'SEAL_TAG'))
# 
# w4 <- bind_rows(wtmp %>% filter(!is.na(SEAL_ID)), wtmp1)
# 
# w5 <- left_join(w3,w4)
# 
# sid <-  w5 %>% dplyr::select(ref, ptt, brand, tag1, tag2, SEAL_ID)
# 
# dups <- sid %>% 
#   filter(!is.na(SEAL_ID)) %>% 
#   mutate(dup = duplicated(SEAL_ID)) %>% 
#   filter(dup == TRUE)
# 
# dup <- na.omit(sid$SEAL_ID)[sid$SEAL_ID %>% na.omit %>% duplicated]
