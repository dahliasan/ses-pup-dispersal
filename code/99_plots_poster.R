# For poster
# need to load 10z_bearing.R first.

library(showtext)
font_add_google("Inter")
showtext_auto()
# Set the plot dimensions in inches
width <- 9 # 3 inches x 3
height <- 6 # 2 inches x 3

# Set the resolution in dots per inch
res <- 300

poster_theme <-   theme(text = element_text(family = "Inter", size = 18))

p1 <- ggplot() + 
  geom_sf(data = locf_sf, size = .1, col = 'grey') +
  geom_sf(data = locw1 %>% filter(id %in% loc_early$id), size = 0.1, col = 'firebrick1') + 
  geom_sf(data = pt_sf, col = 'dodgerblue', size = 0.1, alpha = .3) + 
  # labs(caption = 'blue = particle trace, red = weaners, grey = adult females', x = 'lon', y = 'lat') +
  geom_sf(data = orsi_sf, linetype = 'dashed', colour = 'grey10') + 
  # geom_sf_text(data = orsi_sf, aes(label = front), nudge_x = orsi_sf$nudge_x, nudge_y = orsi_sf$nudge_y) + 
  geom_sf(data = mq_sf, shape = 17, fill = 'black', size = 2 ) +
  xlim(bb[1], bb[2]) + ylim(bb[3]-200, bb[4]+500) +
  theme_bw() +
  poster_theme


png("./output/poster/weaner_particleTrace_adultFemales_1map.png", width = width, height = height, units = "in", res = res)
print(p1)
dev.off()


# Weaner and Female Direction ---------------------------------------------

p1 <- ggplot() +
  geom_histogram(data = data.frame(circ.f), aes(x = circ.f), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  geom_vline(xintercept = mean(circ.f), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.f), y = 20, label = mean(circ.f) %>% round(1)) +
  labs(subtitle = paste('adult females (n=', length(circ.f), ')', sep = ''), y = "") + 
  theme_bw() + 
  poster_theme



# weaner
p2 <- ggplot() +
  # add weaner 
  geom_histogram(data = data.frame(circ.w), aes(x = circ.w), 
                 breaks = seq(0, 360, 45), 
                 colour = NA, 
                 fill = "firebrick3", 
                 alpha = 1) + 
  # add particle trace 
  geom_histogram(data = data.frame(circ.p), aes(x = circ.p), 
                 breaks = seq(0, 360, 45), 
                 fill = NA,
                 color = "black",
                 alpha = 1,
                 size = 1) + 
  
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  # weaner mean direction
  geom_vline(xintercept = mean(circ.w), color = "black", linetype = 1, size = 1) +
  annotate("label", x = mean(circ.w), y = 20, label = mean(circ.w) %>% round(1)) + 
  # particle mean direction
  geom_vline(xintercept = mean(circ.p), color = "pink", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.p), y = 20, label = mean(circ.p) %>% round(1)) +
  labs(subtitle = paste('weaners (n=', length(circ.w), ')',sep = ''), 
       y = "count")  + 
  theme_bw() +
  poster_theme

png('./output/poster/dispersal_direction_v2.png', width=12, height=6, units='in', res=300)
ggarrange(p2, p1, ncol = 2, labels = c('a', 'b'))
dev.off()



# ~~ plot 1st year survival -----------------------------------------------


p3 <- ggplot() +
  geom_histogram(data = data.frame(circ.live), aes(x = circ.live), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  # add mean direction annotation
  geom_vline(xintercept = mean(circ.live), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.live), y = 20, label = mean(circ.live) %>% round(1)) +
  # labs
  labs(subtitle = paste0('1st year survived, n = ', circ.live %>% length), y = "count") + 
  theme_bw() + poster_theme

p4 <- ggplot() +
  geom_histogram(data = data.frame(circ.dead), aes(x = circ.dead), 
                 breaks = seq(0, 360, 45), 
                 colour = "black", 
                 fill = "grey") + 
  coord_polar() +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 45)) +
  # add mean direction annotation
  geom_vline(xintercept = mean(circ.dead), color = "black", linetype = 2, size = 1) +
  annotate("label", x = mean(circ.dead), y = 4, label = mean(circ.dead) %>% round(1)) +
  # labs
  labs(subtitle = paste0('1st year died, n = ', circ.dead %>% length), y = "") + 
  theme_bw() + poster_theme

png('./output/poster/dispersal_direction_by_year1survival.png', width=12, height=6, units='in', res=res)
ggarrange(p3, p4, ncol = 2, labels = c('a', 'b'), align = 'hv')
dev.off()



# Survival models - need to load 11_survival_model.R ----------------------

png('./output/poster/1st trip survival -- avg model effects.png', width=6, height=6, units='in', res=res)
ggplot(data = pred, aes(x = is_ESE, group = is_ESE)) + 
  geom_point(aes(y = fit)) + 
  geom_errorbar(aes(ymin = fit-se.fit, ymax = fit+se.fit), width = .1) +
  labs(y = "1st trip survival", x = "Travelled with the flow") + 
  theme_pubr(border = T) + poster_theme
dev.off()

withFlow_labs <- c("swam against current", "swam with current")
names(withFlow_labs) <- c("FALSE", "TRUE")
png('./output/poster/1st year survival -- avg model effects.png', width=12, height=6, units='in', res=res)
ggplot(data = pred, aes(x = weanmass)) + 
  geom_line(aes(y = fit)) + 
  # geom_point(aes(y = fit)) + 
  geom_ribbon(aes(ymin = fit-se.fit, ymax = fit+se.fit), alpha = 0.1) + 
  facet_wrap(~is_ESE, labeller = labeller(is_ESE = withFlow_labs)) + 
  labs(y = "1st year survival", x = "Weanmass (kg)") + 
  theme_pubr(border = T) +
  theme(panel.spacing = unit(2, "lines")) + poster_theme
dev.off()

