## Code from Martin Biuw to analyse drift dives

rm(list = ls())
# install.packages("remotes")
# remotes::install_github("farcego/slimmingDive")
# devtools::install_github("embiuw/drifteR")  
# devtools::install_github("embiuw/rSRDL")  

library(rSRDL)
library(drifteR)
library(tidyverse) 
library(lubridate)

# mq1 <- get.all.SRDLdb('mq1')
# mq2 <- get.all.SRDLdb('mq2')
# mq3 <- get.all.SRDLdb('mq3')
# mq4 <- get.all.SRDLdb('mq4')

# dive <- plyr::rbind.fill(mq1$dive, mq2$dive, mq3$dive, mq4$dive)
# names(dive)[1] <- 'ref'

load("./Data/mq-ellie-weaners-dive.RData")
dive <-  dive %>% mutate(DE.DATE = DE.DATE %>% mdy_hms())

id_keep <- dir('./Output/ssm plots/keep') %>% str_replace('.png', '')

dive$ref <- dive$ref %>% factor()
dive <- dive %>% filter(ref %in% id_keep)

# Prepare data ------------------------------------------------------------
## Calculate local date
# dive$DE.DATE <- as.POSIXct(strptime(dive$DE.DATE, '%d/%m/%y %H:%M:%S'), tz='UTC')
dive$L.DATE <- dive$DE.DATE+((dive$lon/360)*3600)
dive$L.TIME <- as.POSIXct(strptime(paste('1970-01-01',
                                         format(dive$L.DATE, '%H:%M:%S')),
                                   '%Y-%m-%d %H:%M:%S'), tz='GMT')
dive$L.DAY <- as.POSIXct(strptime(format(dive$L.DATE, '%Y-%m-%d'),
                                  '%Y-%m-%d'), tz='GMT')

## Filter only dives with inflection data
detailed <- unlist(lapply(dive$DEPTH.STR, function(x) length(unlist(strsplit(x, ',')))))
sdive <- dive %>% filter(detailed==max(detailed, na.rm=T))

## Remove tibble class so script works
sdive <- data.frame(sdive)

## Quality check dives
sdive <- qcDive(sdive)
## sdive <- sdive %>% filter(D.MASK>0)
sdive$ref <- factor(sdive$ref)



# Get Drift Dives -----------------------------------------------------------
## Get drift rates and probability weights of each dive being a drift dive
sdive <- rbs(sdive, num=NA)
sdive <- pDrift(sdive, depth.wpar = 50)
sdive$sqDrate <- sdive$drate
sdive$drate <- sign(sdive$drate) * sdive$drate^2


# + Test fit ----------------------------------------------------------------
tstFit <- fitDrate(sdive %>% filter(ref=='mq4-Doris-00'),
##                   weight=c('rRes', 'dWt', 'drWt', 'skew', 'span'), 
                   prefilter=0.5)
fitPlot(tstFit, yLims=c(-0.1, 0.1))

# tstFit.nn <- fitDrate(sdive %>% filter(ref=='mq4-Doris-00'),
#                    ##                   weight=c('rRes', 'dWt', 'drWt', 'skew', 'span'), 
#                    mixture=2, prefilter=0.4)
# 
# fitPlot(tstFit.nn, yLims=c(-0.1, 0.1))

folder <- './Output/drift dive/processDive/temp/'
for(i in 1:nlevels(sdive$ref)) {
  print(i)
  tmp <- try(augmentWt(sdive %>% filter(ref==levels(sdive$ref)[i]), yLims=c(-0.1, 0.1)), silent=T)
  if(class(tmp)!='try-error') {
    tmpFit <- try(fitDrate(tmp, weight='combWeight', prefilter=0.2), silent=T)
    if(class(tmpFit)!='try-error') {
      png(paste(folder, levels(sdive$ref)[i], '.comb.png', sep=''), width=20, height=15, units='cm', 
          res=500, pointsize=10)
      fitPlot(tmpFit, yLims=c(-0.1, 0.1))
      mtext(side=3, levels(sdive$ref)[i])
      dev.off()
    }
  }
}


# + Real fit --------------------------------------------------------------
allFit <- list()
folder <- './Output/drift dive/processDive/allFit/'
for(i in 1:length(unique(sdive$ref))) {
  print(i)
  tmp <- try(augmentWt(sdive %>% filter(ref==levels(sdive$ref)[i]), yLims=c(-0.1, 0.1)), silent=T)
  if(class(tmp)[1]!='try-error') {
    combFit <- try(fitDrate(tmp, weight='combWeight', prefilter=0.3), silent=T)
    augFit <- try(fitDrate(tmp, weight='augWeight', prefilter=0.3), silent=T)
    if(class(combFit)!='try-error' | class(augFit)!='try-error') {
      png(paste(folder, levels(sdive$ref)[i], '.png', sep=''), width=20, height=15, units='cm',
          res=500, pointsize=10)
      par(mfrow=c(2,1), mar=c(2,4,2,1))
      if(class(combFit)!='try-error') {
        fitPlot(combFit, yLims=c(-0.1, 0.1))
      } else {
        plot.new()
      }
      mtext(side=3, paste(levels(sdive$ref)[i], ' - combFit', sep = ''))

      if(class(augFit)!='try-error') {
        fitPlot(augFit, yLims=c(-0.1, 0.1))
      } else {
        plot.new()
      }
      mtext(side=3, paste(levels(sdive$ref)[i], ' - augFit', sep = ''))
      dev.off()
    }
  }
  allFit[[i]] <- list(combFit=combFit, augFit=augFit)
}

names(allFit) <- levels(sdive$ref)
save(allFit, file = './Output/allFit_prefilter=0.3.RData')

# Plot all fits -----------------------------------------------------------
# 
# png('combFit_all.png', width=20, height=16, units='cm',
#     res=500, pointsize=10)
# par(mfrow = c(7, 5), mar = c(2, 4, 2, 1))
# for(i in 1:length(names(allFit))){
#  
#   if(class(allFit[[i]]$combFit)!='try-error') {
#     fitPlot(allFit[[i]]$combFit, yLims=c(-0.1, 0.1))
#     mtext(side=3, names(allFit)[i], cex = .7)
#   } else {
#     plot.new()
#     mtext(side=3, names(allFit)[i], cex = .7)
#   }
# 
# }
# dev.off()
#    

png('augFit_prefilter=0.3_all.png', width=20, height=16, units='cm',
    res=500, pointsize=10)
par(mfrow = c(7, 5), mar = c(2, 4, 2, 1))
for(i in 1:length(names(allFit))){
  
  if(class(allFit[[i]]$augFit)!='try-error') {
    fitPlot(allFit[[i]]$augFit, yLims=c(-0.1, 0.1))
    mtext(side=3, names(allFit)[i], cex = .7)
  } else {
    plot.new()
    mtext(side=3, names(allFit)[i], cex = .7)
  }
  
}
dev.off()
