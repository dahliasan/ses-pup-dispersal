## Kernel Density

library(tidyverse)
library(adehabitatHR)
# library(maps)
# library(mapdata)
library(sp)
library(sf)

# Load Dataset ------------------------------------------------------------

load("./Output/foragingtracks_drift_survival_combined.RData")


# Not my code -------------------------------------------------------------

## Convert to a SpatialPointsDataFrame, necessary for home range analysis
d1 <- loc1
class(d1)
coordinates(d1) = c("lon", "lat") # specify column names
class(d1)
head(d1)


# Step 1: Calculating the Utilisation Distribution (UD) -------------------
# Add small amount of noise using jitter to get LSCV to converge
dd <- jitter(d1@coords, amount=1.2)
d1@coords <- dd

##### Remove any ids that don't converge
#d1 <- d1[d1@data$id_trips20!="4794_bi09_t1_EM",]

############################################ Calculating UD passing a raster map (envt variables)
load("\\\\psf\\Home\\Documents\\PhD\\Data\\Multi-Year Animals\\site fidelity\\envt data\\HomeRange_grid.Rdata") # raster grid of area from sst data (0.25 degree)
str(sst.g)

gridx <- sst.g@grid
gridx <- SpatialPoints(gridx)# Convert first to class SpatialPoints
gridx <- SpatialPixels(gridx)# Then convert to class SpatialPixels
str(gridx)

kudm <- kernelUD(d1[,1], grid=gridx, h="LSCV")
image(kudm)

image(kudm[[1]])  ##e.g. of one animal
map(add = T, fill=T, col="grey")
x <- as.image.SpatialGridDataFrame(kudm[[1]])
contour(x, add=TRUE)

###########################################
## Compare the estimation with ad hoc and LSCV method
#(cuicui1 <- kernel.area(kud)) ## ad hoc
#plot(cuicui1)
#(cuicui2 <- kernel.area(kudl)) ## LSCV
#plot(cuicui2)

## Diagnostic of the cross-validation
plotLSCV(kudl)

####                          STEP 2                         ####
##############################################################################################################################################
########################       Estimating home range from the UD          ###################################################################
##############################################################################################################################################

############################################ Raster Mode                       
vud <- getvolumeUD(kudm)
vud  ## IDENTIFY ANY THAT DON"T CONVERGE UNDER LSCV, then remove from d1 and make kudl again
image(vud)

image(vud[[1]]) #example animal 
xy <- as.image.SpatialGridDataFrame(vud[[1]])
contour(xy, add=TRUE)
map(add = T, fill=T, col="grey")

par(mfrow=c(2,1))
########### 95% Home Range animal 1
fud <- vud[[1]]
hr95 <- as.data.frame(fud)[,1]
hr95 <- as.numeric(hr95 <= 95)
hr95 <- as.data.frame(hr95)
coordinates(hr95) <- coordinates(vud[[1]])
gridded(hr95) <- TRUE
image(hr95)
map(add = T, fill=T, col="grey")

########### 50% Home Range animal 1
fud <- vud[[1]]
hr50 <- as.data.frame(fud)[,1]
hr50 <- as.numeric(hr50 <= 50)
hr50 <- as.data.frame(hr50)
coordinates(hr50) <- coordinates(vud[[1]])
gridded(hr50) <- TRUE
image(hr50)
map(add = T, fill=T, col="grey")

################################ Vector Mode(results in a SpatialPolygonsDataFrame)
#homerange <- getverticeshr(kud)
#class(homerange)
#as.data.frame(homerange)
#plot(homerange, col=1:13)
#map(add = T, fill=T, col="grey")

#image(kud[[1]])
#x <- as.image.SpatialGridDataFrame(kud[[1]])
#contour(x, add=TRUE)

################################################################################
## Determining home range UD cell use for all animals (50% or 95%)

y <- estUDm2spixdf(vud)
fudall <- as.data.frame(y@data)

for (i in 1:ncol(y@data)) {
  fudall[i] <- as.data.frame(y@data[i])
  fudall[i] <- as.numeric(y@data[i] <= 95)
}

coordinates(fudall) <- coordinates(y@coords)
gridded(fudall) <- TRUE
image(fudall[1])

#-------------------------------------------------------------------------------
# create a data frame of cell use for all animals with lat/lons
celluse <- cbind(as.data.frame(fudall@coords), as.data.frame(fudall@data))
colnames(celluse)[1] <- "lon"
colnames(celluse)[2] <- "lat"
head(celluse)

# Estimate overlapped cells and non-overlapped cells for each individual
mc1 <- as.data.frame(apply(celluse[,3:4],1,sum))
mc2 <- as.data.frame(apply(celluse[,5:7],1,sum))
mc3 <- as.data.frame(apply(celluse[,8:9],1,sum))
mc4 <- as.data.frame(apply(celluse[,10:12],1,sum))
mc5 <- as.data.frame(apply(celluse[,13:14],1,sum))
mc6 <- as.data.frame(apply(celluse[,15:16],1,sum))
mc7 <- as.data.frame(apply(celluse[,17:18],1,sum))
mc8 <- as.data.frame(apply(celluse[,19:21],1,sum))
multicell <- cbind(mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8)

a <- colnames(celluse[3:length(celluse)])
nom <- vector()

for (k in 1:length(a)) {
  nom[k] <- substring(a, 1,5)[k]
}

nom <- unique(nom)

colnames(multicell) <- nom
x <- cbind(celluse[,1:2],multicell) 
multicell <- x
head(multicell)

### Save overlap
save(fudall, multicell, file="overlap_95UD_ARS_between_yrs.Rdata")

## Plot overlap
multilap <- multicell[,3:ncol(multicell)]
coordinates(multilap) <- multicell[,1:2]
gridded(multilap) <- TRUE

image(multilap[3], xlim=c(0,80), ylim=c(-70,-30))
map(add = T, fill=T, col="grey")
title(main = list("PP620 o'lap 95% Hrange"))

####                          STEP 3                         ####
##############################################################################################################################################
########################       Calculating home range size                ###################################################################
##############################################################################################################################################
#From raster estimate
hrsize <- kernel.area(kudm, percent=95, unin="km", unout="km2")
hrsize

hrsize <- t(hrsize)
mean(hrsize)
sd(hrsize)
se.hrsize <- sqrt(var(hrsize)/length(hrsize))
se.hrsize
range(hrsize)

####                          STEP 4                         ####
##############################################################################################################################################
########################       Estimating home range overlap              ####################################################################
##############################################################################################################################################

## overlap estimated from the UDs
udoverlap <- kerneloverlaphr(kudm, method="BA", lev=95, conditional=FALSE)  # BA is Bhattacharrya's affinity
udoverlap

save.image("HrangeOverlap_ARS_betweenyrs.Rdata")

#### Pull out overlap values for animals
#ov <- diag(udoverlap)
ov <- as.data.frame(udoverlap)
ov$id <- row.names(ov)    
x <- row.names(ov)
x <- substring(x,1,5)
y <- as.factor(x)

ids <- data.frame()
for (i in 1:nrow(ov)) {
  ov.dat <- ov[ov$id==y[i],]
  ids[i,1] <- y[i] 
}

head(ids)
ids <- as.data.frame(ids[,1])
colnames(ids) <- c("seal")
ov <- cbind(ov,ids)

############################# tabulate overlap for each gls_site level
setwd("\\\\psf/Home/Documents/PhD/Data/Multi-Year Animals/site fidelity/home range/between year/ARS") #change working directory as each trip is going to save out a .csv in the following loop
individ <- data.frame()
z <- levels(ov$seal)

for (j in 1:nlevels(ov$seal))  {
  ind.dat <- ov[ov$seal==z[j],]
  sealtr <- z[j] 
  l <- rownames(ind.dat)  
  individ <- ind.dat[l]  
  meanolp <- matrix(sapply(individ[,1:ncol(individ)], sum)-1) #   ncol(individ)-1) #-1/ncol(individ)-1
  meanolps <- meanolp[1:length(meanolp)]/(nrow(meanolp)-1)          
  #meanolp <- matrix(sapply(individ[,1:ncol(individ)], mean))  
  individ <- cbind(individ, meanolps)
  write.csv(individ, paste(sealtr, "olap_betweenyr_ARS.csv", sep="_"))    ## Change
}

#merge all together
fs3 <- list.files(pattern = "csv")
fs3     

eee <- array()
fff <- array()
for (ki in 1:length(fs3)) {
  fff <- read.csv(fs3[ki])
  j <- nrow(fff)+3
  fff[[j]] <- 0
  colnames(fff)[j] <- (substring(fs3[ki], 1,nchar(fs3[ki])-19))  
  fff[[j]] <- fff[1:j-1]
  fff <- as.list(fff)
  fff <- fff[j]
  eee <- rbind(eee, fff)
}

# after loop
dim(eee)   # number of rows
eee <- eee[2:length(eee)]
eee


# Save out using sink()
sink("Betweenyr_olap_ARS.txt")
eee
sink()


### Some basic stats ###
# Mean overlap 
molp <- c(eee[[1]]$meanolp, eee[[2]]$meanolp, eee[[3]]$meanolp, eee[[4]]$meanolp, 
          eee[[5]]$meanolp, eee[[6]]$meanolp, eee[[7]]$meanolp, eee[[8]]$meanolp)

mean(molp) #MEAN
sd(molp) #SD
se.molp <- sqrt(var(molp)/length(molp))
se.molp #SE
range(molp) #RANGE

save(molp, file="mean_overlap_betweenARS.Rdata")















