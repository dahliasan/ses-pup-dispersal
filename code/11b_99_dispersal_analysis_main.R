# Run all scripts for dispersal analysis

current_date <- format(Sys.Date(), "%Y-%m-%d")
output_path <- file.path("./Output/dispersal_analysis_2", current_date)
seal_data_path <- "./Output/tracks_processed_12h.rds"
particle_data_path <- "./Output/currently_particleTrace.rds"

source("code/11b_1_dispersal_analysis.R")
source("code/11b_2_dispersal_plots.R")
source("code/11b_3_survival_models.R")
