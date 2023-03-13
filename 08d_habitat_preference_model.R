## Modelling

rm(list = ls())
options(tibble.width = Inf)
library(tidyverse)
library(sf)
library(lubridate)
library(foieGras)
library(viridis)
library(tmap)
library(ggpubr)
library(see)
library(orsifronts)
library(pals)
library(lme4)
library(DHARMa)
library(MuMIn)
library(mgcv)


# Load Datasets -----------------------------------------------------------

source('convert2polarsf.R')
load('baseInfo.Rdata')
d1 <- readRDS("./output/all_data_combined.rds")

# + pseudo v real ----------------------------------------------------------

dat <- d1 %>% 
  ungroup() %>% 
  mutate(across(topo:SSHgrad, ~scale(.)[,1]),
         chl = log(chl) %>% scale(.),
         chlgrad = log(chlgrad) %>% scale(.),
         dist_to_ice_m = scale(dist_to_ice_m)[,1], 
         wind = sqrt(windu^2 + windv^2) %>% scale) %>% 
  mutate(real = ifelse(sim == 0, 1, 0)) %>% 
  dplyr::select(id, date, lon, lat, topo:chlgrad, wind, real, -ice)

dat <- na.omit(dat)

## look at correlation between variables 
r <- dat %>% dplyr::select_if(is.numeric) %>% 
  dplyr::select(-lon, -lat, -curru, -currv, -windv, -windu)

rr <- cor(r, use="complete.obs")
round(rr,2) #tri-slope,  sshgrad-curr, disttoice-sst, 

# r <- dat %>% dplyr::select_if(is.numeric) %>% 
#   dplyr::select(-lon, -lat, -curru, -currv, -windv, -windu) %>% 
#   filter(real == 1)
# 
# GGally::ggpairs(r)


# ~ GLMM ------------------------------------------------------------------

options(na.action = "na.omit")
m.full <- glmer(real ~ topo + ssha + eke + slope + SSTgrad + dist_to_ice_m + chl +  (1|id), 
                data = dat, 
                family = binomial)
summary(m.full)

options(na.action = "na.fail")
dd <- dredge(m.full)
dd

# save(dd, file = './Output/dredge1.Rdata')

#models with delta.aicc < 2
summary(model.avg(dd, subset = delta < 2))

#or as a 95% confidence set:
avgmod.95p <- model.avg(dd, cumsum(weight) <= .95)
summary(avgmod.95p)
confint(avgmod.95p)

# save.image('09_habitat_preference_model_workspace.Rdata')


# ~ Generation prediction maps --------------------------------------------
# To generate prediction maps, we calculated the mean of each environmental covariate for the study period (summer and winter separately) based on input data at the same spatio- temporal resolution as that used to model habitat pref- erence. We used ordinary Kriging to interpolate any missing values. Values of these variables were then sampled on a new 0.1° × 0.1° grid which was used for prediction. By predicting to the mean environmen- tal conditions over the study period while matching the locations used in the models to the spatio- temporally nearest environmental covari- ates, the interannual variation typical of this dynamic marine system should be accounted for to some extent. (Reisinger et al 2018)


## Quick and dirty way (not using raster?)
preddat <- dat %>% 
  group_by(lon, lat) %>% 
  summarise(across(c(topo, ssha, eke, slope, SSTgrad, dist_to_ice_m, chl), mean, na.rm = T))

dat %>% 
  ggplot(aes(x = lon)) + 
  geom_histogram(binwidth = .1)




# ~ GAMM --------------------------------------------------------------------

gamm.global <- gamm(real ~
                      # s(topo) +
                      # s(slope) +
                      # s(SSTgrad) +
                      # s(eke) +
                      # s(ssha) +
                      # s(chl) +
                      s(dist_to_ice_m),
                    random = list(id = ~1),
                    family = binomial,
                    # method = "REML",
                    data = dat,
                    # correlation = corAR1(),
                    niterPQL = 1000)


summary(gamm.global$lme)
plot(gamm.global$gam)
