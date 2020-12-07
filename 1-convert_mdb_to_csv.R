library(Hmisc)
library(tidyverse)
library(labelled)


# Read MDB files (Microsoft access, ew.) ----------------------------------

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
