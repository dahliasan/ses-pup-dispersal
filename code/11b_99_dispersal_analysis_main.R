# Run all scripts for dispersal analysis

use_critical_period <- FALSE # set to FALSE to use all data for the seal
filter_overall_particle_mean <- FALSE # set to FALSE to use all data for the particle

# current_date <- format(Sys.Date(), "%Y-%m-%d")

subfolder <- sprintf("use_critical_period_%s_filter_overall_particle_mean_%s", use_critical_period, filter_overall_particle_mean)
output_path <- file.path("./Output/dispersal_analysis_2/0m", subfolder)
seal_data_path <- "./Output/tracks_processed_12h.rds"
particle_data_path <- "./Output/currently_particleTrace.rds"
# particle_data_path <- "output/particle-trace-186.13m.rds"

source("code/11b_1_dispersal_analysis.R")
source("code/11b_2_dispersal_plots.R")
source("code/11b_3_survival_models.R")


list.files("Output/dispersal_analysis_2/", pattern = ".rds", recursive = TRUE, full.names = TRUE) %>%
    file.info() %>%
    arrange(desc(mtime))

results_file <- "Output/dispersal_analysis_2//200m/2024-09-30/all_results.rds"

all_results <- readRDS(results_file)

print_and_save_results(all_results, output_path)
