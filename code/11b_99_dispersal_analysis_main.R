# Run all scripts for dispersal analysis

current_date <- format(Sys.Date(), "%Y-%m-%d")
output_path <- file.path("./Output/dispersal_analysis_2/200m", current_date)
seal_data_path <- "./Output/tracks_processed_12h.rds"
# particle_data_path <- "./Output/currently_particleTrace.rds"
# pt1 <- readRDS(particle_data_path)
particle_data_path <- "output/particle-trace-186.13m.rds"
# pt2 <- readRDS(particle_data_path)

source("code/11b_1_dispersal_analysis.R")
source("code/11b_2_dispersal_plots.R")
source("code/11b_3_survival_models.R")


# # plot pt
# source("code/functions/functions.R")
# pt2_sf <- pt2 %>%
#     na.omit() %>%
#     convert2polarsf()

# pt1_sf <- pt1 %>%
#     na.omit() %>%
#     convert2polarsf()

# ggplot() +
#     geom_sf(data = pt1_sf, col = "blue") +
#     geom_sf(data = pt2_sf, col = "red") +
#     theme_bw()



# ggplot(pt1, aes(x, y)) +
#     geom_point() +
#     geom_point(data = pt2, col = "red") +
#     theme_bw()
