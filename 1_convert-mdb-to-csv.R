library(Hmisc)
library(tidyverse)


# Read MDB files (Microsoft access, ew.) ----------------------------------

mdb <-  mdb.get('./Data/AccessDB/mq1.mdb')
mdb2 <-  mdb.get('./Data/AccessDB/mq2.mdb')

df <- mdb$diag %>% as_tibble() %>% select(-PTT) %>% mutate(ACTUAL.PTT = ACTUAL.PTT %>% as.character)
df2 <- mdb2$diag %>% as_tibble() %>% select(-PTT) %>% mutate(ACTUAL.PTT = ACTUAL.PTT %>% as.character)
df <- bind_rows(df, df2)

# write_csv(df, './Data/mq-ellie-weaners-argos.csv')