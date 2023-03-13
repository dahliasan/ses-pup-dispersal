## Drift dive analysis

# install.packages("remotes")
# remotes::install_github("farcego/slimmingDive")
library(tidyverse)
library(lubridate)
library(slimmingDive)
library(ggpubr)


# Step 1: Prep data ----------------------------------------------------------
# w1 <- read_csv("./Data/mq-ellie-weaners-dive.csv")
load("./Data/mq-ellie-weaners-dive.RData")
w1 <- dive
w1 <- w1 %>% 
  mutate(DE.DATE = DE.DATE %>% mdy_hms())

## Filter only rows with inflection points data
w2 <- w1 %>% filter(!is.na(D1))

## Format dataset for slimmingDive
w2 <-  w2 %>%
  dplyr::select(ref, DE.DATE, SURF.DUR, DIVE.DUR, MAX.DEP, D1, D2, D3, D4, T1, T2, T3, T4, lat, lon)
colnames(w2) <- str_replace(colnames(w2), '\\.', '_')
Data <- formatDives(w2) %>% data.frame()


# Step 2: Drift dive analysis --------------------------------------------------
Data[c('order','minresid')] <- t(apply(Data,1,RBSM,retrieve='both'))
Data$minresid <- as.numeric(Data$minresid)
D1 <- newVarsVect(Data)
D1[c('NDE','ds')] <- t(apply(D1,1,NDE, extract='both'))

## Remove non-drift dives
plotDrift(D1, xlab = 'foraging days', ylab = 'Drift rate',
          cex.lab = 1.5, cex = 2, main = '1. Raw drift rates (mq2 only)')
D2 <- driftFilter(D1)
plotDrift(D2, xlab = 'foraging days', ylab = 'Drift rate', cex.lab = 1.5, cex = 2, main = '2. rbsm filter (mq2 only)')

## Kalman Filter (optional? - further filtering of non-drift dives)
D3 <- kalman(D2,400000, 10000, n.adap = 10000)
class(D3) <- append(class(D3), 'kalman')
D4 <- postKalman(D3)
D4 <- D4[D4$zetas > .5, ]
plotDrift(D4, xlab = 'foraging days',cex.lab = 1.5, cex = 2,
          ylab = 'Drift rate', main = '3. kalman filter (mq2 only)') # Drift rate is column NDE


## Make the gam
D5 <- makeTheGam(D4) 
plotDrift(D4, xlab = 'foraging days',cex.lab = 1.5, cex = 2,
          ylab = 'Drift rate', ) 

DD <- D4 %>% mutate(season = purrr::map(ref, function(x) str_split(x, '-')[[1]][1]) %>% unlist()) %>% 
  select(ref, season, Date, NDE) %>% 
  as_tibble()

fo <- purrr::map(D5, ~as_tibble(.))
fo[[1]]$season <- 'mq3'
fo[[2]]$season <- 'mq3'
fo[[3]]$season <- 'mq4'
fo[[4]]$season <- 'mq4'

DD %>% 
  ggplot(aes(x = Date, y = NDE)) + 
  geom_point(size = 2) + 
  ylim(c(-0.45, 0.3)) + 
  scale_x_datetime(date_breaks = '2 months', date_labels = '%b') +
  theme_pubr(border = T) + 
  labs(x = 'Foraging Days', y = 'Drift Rate', title = '4. Final with fitted GAMs') +
  facet_wrap(~season, scales = 'free') +
  geom_line(data = fo[[1]], aes(y = pred), color = 'red', size = 1.5) + 
  geom_line(data = fo[[2]], aes(y = pred), color = 'red', size = 1.5) + 
  geom_line(data = fo[[3]], aes(y = pred), color = 'red', size = 1.5) +
  geom_line(data = fo[[4]], aes(y = pred), color = 'red', size = 1.5) 




# Individual drift rate time series ---------------------------------------
d2 <- D4 %>% as_tibble() %>% 
  select(ref, Date, NDE)

d2 %>% 
  ggplot(aes(x = Date, y = NDE)) + 
  geom_point() + 
  geom_smooth() + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  scale_x_datetime(date_breaks = '2 months', date_labels = '%b') +
  facet_wrap(~ref, scales = 'free') + 
  theme_pubr(border = T) +
  ylim(c(min(d2$NDE), max(d2$NDE)))

  