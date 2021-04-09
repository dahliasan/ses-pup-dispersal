library(Hmisc)
library(tidyverse)
library(labelled)


# Read MDB files (Microsoft access, ew.) ----------------------------------

# Extract ARGOS locations -------------------------------------------------

mdb_flist <- dir('./Data/AccessDB', full.names = T)

mdb_list <- purrr::map(mdb_flist, function(x){
  mdb <-  mdb.get(x)
  df <- mdb$diag %>% as_tibble() %>% select(1:6)
  colnames(df) <- tolower(colnames(df)) # change all colnames to lowercase
  df <- remove_labels(df)
  df
})

mdb <- mdb_list %>% bind_rows()

write_csv(mdb, './Data/mq-ellie-weaners-argos.csv')


# Extract dive data -------------------------------------------------------

mdb_flist <- dir('./Data/AccessDB', full.names = T)

mdb_list <- purrr::map(mdb_flist, function(x){
  mdb <-  mdb.get(x)
  df <- mdb$dive %>% as_tibble() 
  # # remove columns with all NA
  # cols_na <- purrr::map(df, function(x){
  #   is.na(x) %>% all()
  # }) %>% unlist()
  # df <- df[,!cols_na]
  df <- remove_labels(df)
  df
})

## Fix different column names of mq3 and mq4 data before combining all
names(mdb_list[[3]])[1] <- 'ref' 
names(mdb_list[[4]])[1] <- 'ref'

# Combine all mq data
# dive <- mdb_list %>% bind_rows()
dive <- mdb_list %>% plyr::rbind.fill()
dive <- dive %>% as_tibble()

save(dive, file = './Data/mq-ellie-weaners-dive.RData')

# write_csv(dive, './Data/mq-ellie-weaners-dive.csv')


# Extract haulout data -------------------------------------------------------

mdb_flist <- dir('./Data/AccessDB', full.names = T)

mdb_list <- purrr::map(mdb_flist, function(x){
  mdb <-  mdb.get(x)
  df <- mdb$haulout_orig %>% as_tibble() 
  df <- remove_labels(df)
  df
})

## Fix different column names of mq3 and mq4 data before combining all
names(mdb_list[[3]])[1] <- 'ref' 
names(mdb_list[[4]])[1] <- 'ref'

# Combine all mq data
# dive <- mdb_list %>% bind_rows()
haul <- mdb_list %>% plyr::rbind.fill()
haul <- haul %>% as_tibble()

save(haul, file = './Data/mq-ellie-weaners-haulout.RData')



# Capture Mark Recapture --------------------------------------------------

cmc_file <- dir('./Data/CMC/TAGS', pattern = '.mdb', full.name = T)
mdb <-  mdb.get(cmc_file)
