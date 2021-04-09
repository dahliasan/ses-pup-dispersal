TAGS_sightings <- function(seal_tag='P001',
                           db='D:/emb7/My Documents/Weaner project/TAGS/TAGS/Elies_db.accdb') {
  require(RODBC)
  con <- odbcDriverConnect(paste("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=",db))
  seal_id <- sqlQuery(con,
                      paste("SELECT DISTINCT SEAL_ID FROM TAGS_ELLIE_TAGS WHERE SEAL_TAG=",
                            paste("'", seal_tag, "'", sep=''), sep=''))
  seal_sex <- as.character(sqlQuery(con,
                                    paste("SELECT DISTINCT SEX FROM TAGS_ELLIE_SEALS WHERE SEAL_ID=",
                                          seal_id, sep=''), as.is=T))
  sightings <- sqlQuery(con,
                        paste("SELECT * FROM TAGS_ELLIE_SIGHTINGS WHERE SEAL_ID=",
                              seal_id, sep=''))
  odbcClose(con)
  sightings$COHORT <- format(min(sightings$DATE_OBSERVED), '%Y')
  sightings$SEX <- rep(seal_sex, nrow(sightings))
  sightings
}
