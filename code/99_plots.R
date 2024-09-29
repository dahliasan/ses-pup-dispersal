## Make plots



d1 <- readRDS("./Output/all_data_combined.rds")
dr <- d1 %>% filter(sim == 0) # only real locations

dr %>% 
  ggplot(aes(x = drate_chg)) + 
  geom_histogram() + 
  facet_wrap(~id, scale = 'free')

# Plots -------------------------------------------------------------------

d1_sf <- convert2polarsf(d1)
dr_sf <- convert2polarsf(dr)
bb <- extent(dr_sf) 


# + timeline of trips -----------------------------------------------------
# png('full_deployment_timelines.png', width=30, height=20, units='cm', res=500)
dr %>%  
  ggplot(aes(x = date, y = id)) + 
  geom_point(aes(colour = trip %>% factor), size = .5) +
  geom_point(data = dr %>% filter(land == TRUE), size = .5, colour = 'red') +
  geom_point(data = dr, aes(x = birthdate ), colour = 'black', size = 1) +
  geom_point(data = dr, aes(x = weandate), colour = 'black', size = 1) +
  scale_x_date(date_breaks = '3 months', date_labels = '%m/%y', date_minor_breaks = '1 month') +
  scale_color_viridis_d(name = 'trip') +
  labs(title = 'Trip Timeline (black dots = birth and wean date; red dot = at colony)') + 
  theme_linedraw()
# dev.off()


# + map of tracks by 1st trip survival ---------------------------------------
orsi_sf$nudge_x <- c(-1000, -4800, -4500, -4000)
orsi_sf$nudge_y <- c(0, 5300, 4800, 3000)
orsi_sf <- orsi_sf %>% 
  mutate(front = front %>% toupper())

orsi_sf %>% 
  ggplot() +
  geom_sf() + 
  geom_sf_text(aes(label = front), 
               nudge_x = orsi_sf$nudge_x,
               nudge_y = orsi_sf$nudge_y)

(p1 <- ggplot() +
    geom_sf(data = dr_sf %>% filter(trip == 1, !is.na(seen_6m), SUS == FALSE), 
            aes(colour = seen_6m), size = 0.1) +
    # map components
    geom_sf(data = orsi_sf, linetype = 'dashed', colour = 'grey10') + 
    geom_sf_text(data = orsi_sf, aes(label = front), nudge_x = orsi_sf$nudge_x, nudge_y = orsi_sf$nudge_y) + 
    geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
    geom_sf(data = mq_sf, size = 1, colour = 'yellow', shape = 17) + 
    xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) + 
    labs(title = 'If seal seen after first trip') +
    scale_color_see() + 
    theme_linedraw(base_size = 8))

png('tracks_survival_1map.png', width=25, height=20, units='cm', res=500)
p1
dev.off()


# + Plot drift rates -----------------------------------------------------------

dr %>% 
  filter(!is.na(drate), trip == 1) %>% 
  mutate(drate = drate*100) %>% 
  ggplot(aes(x = daysFromDeparture, y = drate)) +
  geom_point() +
  geom_smooth() + 
  facet_wrap(~id)


# + Tracks vs environmental variables (TODO) -------------------------------------

# + Assess tag loss -------------------------------------------------------
# Plot track x first trip seen x first trip complete
(p2 <- ggplot() +
   geom_sf(data = dr_sf %>% filter(trip == 1), 
           aes(colour = is_complete), size = 0.2) +
   geom_sf(data = dr_sf %>% filter(trip == 1, SUS == TRUE), 
           size = 0.2, 
           color = 'grey40') + 
   geom_sf(data = world_sf, fill = grey(0.4), colour = NA) +
   geom_sf(data = mq_sf, size = 1, colour = 'black') + 
   xlim(bb[1], bb[2]) + ylim(bb[3], bb[4]) + 
   labs(title = 'Seal seen after trip 1 + complete trip 1 track. Tag loss = seen again + incomplete track') + 
   facet_wrap(seen_6m~id) +
   theme_linedraw(base_size = 8))

png('tracks_tagloss.png', width=25, height=20, units='cm', res=500)
p2
dev.off()

tagloss <- dr %>% filter(trip == 1) %>% 
  group_by(id) %>% 
  summarise(is_complete = first(is_complete), 
            seen_6m = first(seen_6m))

table(tagloss$is_complete, tagloss$seen_6m) # tag loss n = 13

