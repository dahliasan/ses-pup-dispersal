# Check if output_path exists
if (!exists("output_path")) {
    cat("Output directory does not exist.\n")
    user_input <- readline(prompt = "Enter a valid output path: ")
    output_path <- paste0(user_input, Sys.Date())
}

# Create the directory if it doesn't exist
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

sink(paste0(output_path, "/model_output.txt"))

print("Starting survival model analysis...")

library(tidyverse)
library(lubridate)
library(sf)
library(MuMIn)
library(broom)

# Do binomial GLM for survival
modal_data <- weaner_locs_sf %>%
    select(id, seal_mean_bearing, is_following, weanmass, birthyear, size, survive_trip_1, survive_year_1) %>%
    st_drop_geometry() %>%
    distinct()

# names(all_data_weaners)

predictor_vars <- c("sst", "ssha", "eke", "tri", "slope", "SSTgrad", "ice", "dist_to_ice_m", "chl", "chlgrad")

predictor_data <- all_data_weaners %>%
    group_by(id) %>%
    summarise(across(all_of(predictor_vars), ~ mean(., na.rm = TRUE)))

modal_data <- modal_data %>% left_join(predictor_data, by = "id")

cat("\n\nModal data summary:\n")
print(modal_data %>% summary())

source("code/11b_3.1_get_high_cor_vars.R")
# highly correlated variables: "tri" "dist_to_ice_m" "chl"

model_trip <- glm(survive_trip_1 ~ is_following * weanmass + birthyear + sst + ssha + eke + slope + SSTgrad + ice + chlgrad,
    data = modal_data, family = binomial(link = "logit"), na.action = "na.fail"
)

cat("\n\n--- Summary of model_trip (global model) ---\n")
print(summary(model_trip))

dredge_trip <- dredge(model_trip)

best_models <- get.models(dredge_trip, subset = delta <= 2)

# Create a function to extract important metrics
extract_model_info <- function(model, model_name) {
    glance_info <- glance(model)
    tibble(
        model_name = model_name,
        formula = paste(format(formula(model)), collapse = " "),
        AICc = AICc(model), # Use AICc instead of AIC
        BIC = glance_info$BIC,
        deviance = glance_info$deviance,
        df.residual = glance_info$df.residual,
        null.deviance = glance_info$null.deviance,
        df.null = glance_info$df.null
    )
}

# Apply the function to each model
model_summary <- map2_dfr(best_models, names(best_models), extract_model_info) %>%
    mutate(delta_AICc = AICc - min(AICc)) %>% # Use AICc instead of AIC
    arrange(delta_AICc)

cat("\n\n--- Model selection summary for trip survival ---\n")
print(model_summary)

avg_model <- model.avg(best_models)
cat("\n\n--- Model averaging summary for trip survival ---\n")
print(summary(avg_model))

importance <- sw(dredge_trip)
cat("\n\n--- Variable importance for trip survival ---\n")
print(importance)

## Year survival
model_year <- glm(survive_year_1 ~ is_following * weanmass + birthyear + sst + ssha + eke + slope + SSTgrad + ice + chlgrad,
    data = modal_data, family = binomial(link = "logit"), na.action = "na.fail"
)

cat("\n\n--- Summary of model_year (global model) ---\n")
print(summary(model_year))

dredge_year <- dredge(model_year)

best_models_year <- get.models(dredge_year, subset = delta <= 2)

model_summary_year <- map2_dfr(best_models_year, names(best_models_year), extract_model_info) %>%
    mutate(delta_AICc = AICc - min(AICc)) %>%
    arrange(delta_AICc)

cat("\n\n--- Model selection summary for year survival ---\n")
print(model_summary_year)

avg_model_year <- model.avg(best_models_year)
cat("\n\n--- Model averaging summary for year survival ---\n")
print(summary(avg_model_year))

importance_year <- sw(dredge_year)
cat("\n\n--- Variable importance for year survival ---\n")
print(importance_year)

cat("\n\n--- Analysis complete. Output saved. ---\n")

sink()
