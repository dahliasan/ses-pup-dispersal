# Run all scripts for dispersal analysis

current_date <- format(Sys.Date(), "%Y-%m-%d")
output_path <- file.path("./Output/dispersal_analysis_2/200m", current_date)
seal_data_path <- "./Output/tracks_processed_12h.rds"
# particle_data_path <- "./Output/currently_particleTrace.rds"
particle_data_path <- "output/particle-trace-186.13m.rds"

source("code/11b_1_dispersal_analysis.R")
source("code/11b_2_dispersal_plots.R")
source("code/11b_3_survival_models.R")
