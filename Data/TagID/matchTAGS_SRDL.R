mqPTT <-tabulizer::extract_tables('C:/Users/a5406/Documents/Articles/McConnell et al. JAE2002.pdf', 
                                  5, encoding='utf-8') 

nam <- trimws(apply(mqPTT[[1]][c(1:3),], 2, function(x) paste(x, collapse=' ')))

mqPTT <- as.data.frame(mqPTT[[1]][-c(1:3),])

names(mqPTT) <- nam

mqPTT$`Seal code` <- gsub('â\200“', '-', mqPTT$`Seal code`)
mqPTT$`Seal code` <- paste('mq', mqPTT$`Seal code`, sep='-')
mqPTT <- mqPTT[-match(c('HF', 'HM', 'H', 'LF'), mqPTT$Seal),]

mqPTT$Cohort <- paste('19', substr(mqPTT$`Seal code`, nchar(mqPTT$`Seal code`)-1, nchar(mqPTT$`Seal code`)), sep='')

mqPTT$`Weaning mass (kg)` <- as.numeric(mqPTT$`Weaning mass (kg)`)


mqTAGS <- read.csv("D:/emb7/My Documents/Weaner project/Brand+tags 95-96 PTT's.csv", header=F)
names(mqTAGS) <- c('Brand', 'Tag1', 'Tag2')

mqTAGS$Brand[which(nchar(mqTAGS$Brand)==2)] <- paste(substr(mqTAGS$Brand[which(nchar(mqTAGS$Brand)==2)], 1, 1),
                                                     '00', substr(mqTAGS$Brand[which(nchar(mqTAGS$Brand)==2)], 2, 2), sep='')
mqTAGS$Brand[which(nchar(mqTAGS$Brand)==3)] <- paste(substr(mqTAGS$Brand[which(nchar(mqTAGS$Brand)==3)], 1, 1),
                                                     '0', substr(mqTAGS$Brand[which(nchar(mqTAGS$Brand)==3)], 2, 3), sep='')

mqTAGS$PTT <- mqTAGS$SMRUmass <- mqTAGS$TAGSmass <- rep(NA, nrow(mqTAGS))


for(i in 1:nrow(mqTAGS)) {
  tmp <- TAGS_sightings(mqTAGS$Brand[i])
  mqTAGS$TAGSmass[i] <- tmp$WEIGHT[which(tmp$AGE_CLASS_ID==2)]
  candidates <- which(mqPTT$Cohort==tmp$COHORT & substr(mqPTT$Seal, 2, 2)==tmp$SEX[which(tmp$AGE_CLASS_ID==2)])
  candidate.mass <- mqPTT$`Weaning mass (kg)`[candidates]
  closest.mass <- which.min(abs(candidate.mass-tmp$WEIGHT[which(tmp$AGE_CLASS_ID==2)]))
  mqTAGS$PTT[i] <- mqPTT$`Seal code`[candidates[closest.mass]] 
  mqTAGS$SMRUmass[i] <- mqPTT$`Weaning mass (kg)`[candidates[closest.mass]]
}

duplicates <- which(table(mqTAGS$PTT)>1)

for(i in 1:length(duplicates)) {
  which.name <- which(mqTAGS$PTT==names(duplicates)[i])
  closest.mass <- which.min(abs(mqTAGS$TAGSmass-mqTAGS$SMRUmass))
  which.not <- which.name[-closest.mass]
  mqTAGS$SMRUmass[which.not] <- NA
  mqTAGS$PTT[which.not] <- NA
}

write.csv(mqTAGS, 'mqPTT_95_96.csv', row.names=F)
