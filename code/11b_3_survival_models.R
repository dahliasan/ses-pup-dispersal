output_path <- paste0("output/dispersal_analysis_2/", Sys.Date())
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

sink(paste0(output_path, "/model_output.txt"))

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

names(all_data_weaners)

predictor_vars <- c("sst", "ssha", "eke", "tri", "slope", "SSTgrad", "ice", "dist_to_ice_m", "chl", "chlgrad")

predictor_data <- all_data_weaners %>%
    group_by(id) %>%
    summarise(across(all_of(predictor_vars), ~ mean(., na.rm = TRUE)))

modal_data <- modal_data %>% left_join(predictor_data, by = "id")

modal_data %>% summary()

model_trip <- glm(survive_trip_1 ~ is_following * weanmass + birthyear + sst + ssha + eke + tri + slope + SSTgrad + ice + dist_to_ice_m + chl + chlgrad,
    data = modal_data, family = binomial(link = "logit"), na.action = "na.fail"
)

dredge_trip <- dredge(model_trip)

best_models <- get.models(dredge_trip, subset = delta <= 2)

best_models %>% summary()

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
model_summary <- map2_dfr(best_models, names(best_models), extract_model_info)

# Add delta AICc
model_summary <- model_summary %>%
    mutate(delta_AICc = AICc - min(AICc)) %>% # Use AICc instead of AIC
    arrange(delta_AICc)

# Print the summary
print(model_summary)

avg_model <- model.avg(best_models)
summary(avg_model)

importance <- sw(dredge_trip)
print(importance)

# Interpretation of variable importance:
# - sst (0.71) and is_following (0.58) show the strongest evidence of importance,
#   appearing in many of the best models.
# - dist_to_ice_m (0.46) shows moderate importance.
# - eke, ssha, ice, SSTgrad, slope, tri, weanmass, birthyear, chlgrad, and chl
#   (all between 0.24-0.35) show weaker evidence of importance.
# - The interaction is_following:weanmass (0.03) shows very little evidence of importance.
#
# This suggests that sea surface temperature and whether the seal is following
# currents are the most influential factors in the model, while the interaction
# between following currents and wean mass is least important.

model_year <- glm(survive_year_1 ~ is_following * weanmass + birthyear + sst + ssha + eke + tri + slope + SSTgrad + ice + dist_to_ice_m + chl + chlgrad,
    data = modal_data, family = binomial(link = "logit"), na.action = "na.fail"
)

dredge_year <- dredge(model_year)

best_models_year <- get.models(dredge_year, subset = delta <= 2)

best_models_year %>% summary()

model_summary_year <- map2_dfr(best_models_year, names(best_models_year), extract_model_info)

model_summary_year <- model_summary_year %>%
    mutate(delta_AICc = AICc - min(AICc)) %>%
    arrange(delta_AICc)

print(model_summary_year)

avg_model_year <- model.avg(best_models_year)
summary(avg_model_year)

importance_year <- sw(dredge_year)
print(importance_year)

sink()
