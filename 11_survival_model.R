# survival model

rm(list = ls())
options(tibble.width = Inf)
library(tidyverse)
library(sf)
library(lubridate)
# library(foieGras)
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
library(effects)
library(ggeffects)
library(rstatix)
# library(sjPlot)

load('11_survival_model_workspace.RData')


# Load Datasets -----------------------------------------------------------

source('convert2polarsf.R')
load('baseInfo.Rdata')
source("functions.R")
d1 <- readRDS("./Output/all_data_combined.rds")
# d1 <- readRDS("./Output/all_data_combined__REF_DATE=6d.rds")
dr <- d1 %>% filter(sim == 0, trip == 1, land == FALSE, !is.na(bearing.pt)) # only real locations

# Determine if seal swimming with the current -----------------------------
dr$compass_zone.pt <- sapply(dr$bearing.pt, whichZone)

# Some exploratory plots --------------------------------------------------
# look at tracks of seals that wasn't seen again after 6m but track was completed.
tmp <- dr %>% 
  filter(seen_6m == FALSE, is_trip_complete == TRUE) 

tmp <- convert2polarsf(tmp)
tmp %>% ggplot() + 
  geom_sf(aes(col = daysFromDeployment), size = .1) + 
  geom_sf(data = convert2polarsf(mq), col = 'red', size = 0.1) + 
  facet_wrap(~id) +
  scale_color_viridis() +
  labs(caption = "seals that were not seen at the 6 month resight, but their tracks were complete.")



# Model -------------------------------------------------------------------
# + 1st trip survival -----------------------------------------------------
## prepare model data
dat <-dr %>% 
  filter(trip == 1, land == FALSE) %>% 
  group_by(id) %>% 
  arrange(id, trip) %>% 
  summarise(
    year = first(birthyear),
    surviveTrip1 = ifelse(first(seen_6m) == TRUE | first(is_trip_complete) == TRUE, 1, 0),
    surviveYear1 = ifelse(first(seen_1y) == TRUE, 1, 0),
    # across(c(topo:ice, chl, chlgrad, dist_to_ice_m), mean, na.rm = T),
    across(c(blackmass, weanmass, bearing_diff, ew_zone, compass_zone, compass_zone.pt), first),
    # g_sd = sd(g),
    # g_mean = mean(g),
    # g_propARS = sum(g<=0.5)/n(),
    is_ESE = ifelse(compass_zone == "E-SE", TRUE, FALSE))

dat <- dat %>% filter(!is.na(weanmass), !is.na(is_ESE))
# dat <- dat %>% mutate(across(topo:bearing_diff, ~scale(.)[,1]))

# seals that survived 1st trip
table(dat$surviveTrip1 == TRUE)
table(dat$year)

options(na.action = 'na.omit')

m.global <- glm(surviveTrip1 ~ weanmass*is_ESE + year, data = dat, family = binomial)
options(na.action = 'na.fail')

# m.global <- glmer(surviveTrip1 ~ weanmass*is_ESE + (1|year), data = dat, family = binomial) #isSingular warning therefore did not include as a random term. Furthermore does it make biological sense?

# # does weaning mass vary between year?
# rstatix::anova_test(dat, weanmass ~ factor(year))
# rstatix::tukey_hsd(dat, weanmass ~ factor(year))
# dat %>% 
#   ggplot(aes(x = factor(year), y = weanmass, group = year)) + 
#   geom_boxplot()

summary(m.global)
dd <- dredge(m.global)
dd

# write.csv(dd %>% mutate(across(where(is.numeric), ~signif(., 3))),
          # "./output/surivalTrip1__dredge__model_selection.csv", na = "")



avgmod.2delta <- model.avg(dd, subset = delta < 2, fit = TRUE)
summary(avgmod.2delta)
importance(avgmod.2delta) 

# save model summary
# write.csv(summary(avgmod.2delta)$coefmat.subset %>% data.frame(), "./output/1st trip survival model table.csv")

m.final <- glm(surviveTrip1 ~  is_ESE, data = dat, family = binomial)
summary(m.final)
effects::allEffects(m.final) %>% plot()

simulationOutput <- simulateResiduals(fittedModel = m.global, n = 1000)
plot(simulationOutput)

# ~ plot avg model --------------------------------------------------------
a <- tibble(is_ESE = TRUE, weanmass = dat$weanmass)
b <- tibble(is_ESE = FALSE, weanmass = dat$weanmass)

pred <- predict(avgmod.2delta, newdata = a, se.fit = TRUE) %>% as_tibble() %>% 
  bind_cols(tibble(is_ESE = TRUE, weanmass = dat$weanmass))

pred <- pred %>% bind_rows(
  predict(avgmod.2delta, newdata = b, se.fit = TRUE) %>% as_tibble() %>% bind_cols(tibble(is_ESE = FALSE, weanmass = dat$weanmass))
)

# png('./output/survival model/1st trip survival -- avg model effects.png', width=15, height=15, units='cm', res=500)
  ggplot(data = pred, aes(x = is_ESE, group = is_ESE)) + 
    geom_point(aes(y = fit)) + 
    geom_errorbar(aes(ymin = fit-se.fit, ymax = fit+se.fit), width = .1) +
    labs(y = "1st trip survival", x = "travelled with the flow") + 
    theme_pubr(border = T)
# dev.off()




# + 1st year survival -----------------------------------------------------
# seals that survived 1st year
table(dat$surviveYear1 == TRUE)

options(na.action = 'na.omit')

m.global <- glm(surviveYear1 ~ weanmass*is_ESE + year, data = dat, family = binomial)
options(na.action = 'na.fail')
dd <- dredge(m.global)
dd

# write.csv(dd %>% mutate(across(where(is.numeric), ~signif(., 3))),
          # "./output/surivalYear1__dredge__model_selection.csv", na = "")

avgmod.2delta <- model.avg(dd, subset = delta < 2, fit = TRUE)
summary(avgmod.2delta)
importance(avgmod.2delta)

# write.csv(summary(avgmod.2delta)$coefmat.subset %>% data.frame(), "./output/1st year survival model table.csv")

# Test assumptions
m.final <- glm(surviveYear1 ~  weanmass*is_ESE, data = dat, family = binomial)
summary(m.final)
effects::allEffects(m.final) %>% plot()

simulationOutput <- simulateResiduals(fittedModel = m.global, n = 1000)
plot(simulationOutput)

# ~ plot avg model --------------------------------------------------------
a <- tibble(is_ESE = TRUE, weanmass = dat$weanmass)
b <- tibble(is_ESE = FALSE, weanmass = dat$weanmass)

pred <- predict(avgmod.2delta, newdata = a, se.fit = TRUE) %>% as_tibble() %>% 
  bind_cols(tibble(is_ESE = TRUE, weanmass = dat$weanmass))

pred <- pred %>% bind_rows(
  predict(avgmod.2delta, newdata = b, se.fit = TRUE) %>% as_tibble() %>% bind_cols(tibble(is_ESE = FALSE, weanmass = dat$weanmass))
)

# config facet labels
withFlow_labs <- c("is_ESE = FALSE", "is_ESE = TRUE")
names(withFlow_labs) <- c("FALSE", "TRUE")

# png('./output/survival model/1st year survival -- avg model effects.png', width=20, height=13, units='cm', res=500)
  ggplot(data = pred, aes(x = weanmass)) + 
    geom_line(aes(y = fit)) + 
    # geom_point(aes(y = fit)) + 
    geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = 0.1) + 
    facet_wrap(~is_ESE, labeller = labeller(is_ESE = withFlow_labs)) + 
    labs(y = "1st year survival", x = "weanmass (kg)") + 
    theme_pubr(border = T) +
    theme(panel.spacing = unit(2, "lines"))
# dev.off()

ggplot(data = pred, aes(x = weanmass, group = is_ESE, colour = is_ESE)) + 
  geom_line(aes(y = fit)) + 
  # geom_line(aes(y = fit - se.fit), linetype = "dashed") + 
  # geom_line(aes(y = fit + se.fit), linetype = "dashed") + 
  labs(y = "1st year survival", x = "weanmass (kg)") + 
  theme_pubr(border = T)


# Weanmass vs going with flow ---------------------------------------------

anova_test(dat, weanmass ~ is_ESE) # n.s.

dat_filtered <- dat %>% filter(surviveTrip1 == 1)


png('./output/weanmass_vs_withFlow_allSurvivalOutcomes.png', width=15, height=15, units='cm', res=500)
dat %>% 
  ggplot(aes(x = is_ESE, y = weanmass, group = is_ESE)) + 
  geom_point(aes(color = ifelse(surviveTrip1 == 1, T, F))) + 
  geom_point(data = dat_filtered, 
             aes(color = ifelse(surviveYear1 == 1, T, F)), 
             position = position_nudge(x = .1), 
             shape = 17) + 
  geom_segment(data = dat_filtered, aes(x = is_ESE + 1.02, xend = is_ESE + 1.08, yend = weanmass),
               arrow = arrow(length = unit(0.01, "npc")), color = 'grey') +
  labs(color = "survived", y = "weaning mass (kg)", x = "went with flow",
       caption = "circle = 1st trip survival, triangle = 1st year survival") + 
  scale_color_manual(values = c("firebrick1", "dodgerblue")) +
  theme_bw()
dev.off()



# # 1st trip to 1st year survival -------------------------------------------
# # create post trip survival dataframe
# sur <- dr %>% 
#   group_by(id) %>% 
#   summarise(trip_survive = ifelse(first(seen_6m) == TRUE | first(is_trip_complete == TRUE), 1, 0), 
#             seen_1y = first(seen_1y), 
#             weanmass = first(weanmass), 
#             ew_zone = first(ew_zone),
#             birthlocation = first(birthlocation),
#             fieldseason = first(fieldseason)) %>% 
#   mutate(post_trip_died = trip_survive - seen_1y,
#          post_trip_died = as.logical(post_trip_died),)
# 
# sur_trip1 <- sur %>% 
#   filter(trip_survive == 1)
# 
# t.test(sur_trip1$weanmass[sur_trip1$post_trip_died == TRUE], sur_trip1$weanmass[sur_trip1$post_trip_died == FALSE])
# 
# 
# sur_trip1 %>% 
#   ggplot(aes(y = weanmass, x = post_trip_died)) + 
#   geom_boxplot()
# 
# ## Lighter seals died after their first trip.
# 
# 

# # Weanmass vs birth location ----------------------------------------------
# sur %>% 
#   ggplot(aes(y = weanmass, x = birthlocation)) + 
#   geom_boxplot(aes(fill = seen_1y), position=position_dodge(.9))
# 
# sur %>% 
#   ggplot(aes(x = birthlocation)) + 
#   geom_bar(aes(fill = seen_1y)) + 
#   facet_wrap(~ew_zone)
# 
# sur %>% 
#   ggplot(aes(x = fieldseason)) + 
#   geom_bar(aes(fill = seen_1y)) + 
#   facet_wrap(~ew_zone)
#   
# sur %>% 
#   ggplot(aes(y = weanmass, x = fieldseason)) + 
#   geom_boxplot(aes(fill = seen_1y), position=position_dodge(.9))
# 
# sur %>% 
#   ggplot(aes(x = birthlocation)) + 
#   geom_bar(position = 'stack', aes(fill = fieldseason)) + 
#   scale_fill_bluebrown()


# Save Workspace ----------------------------------------------------------


save.image('11_survival_model_workspace.RData')
