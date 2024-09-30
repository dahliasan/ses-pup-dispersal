Methods
Data Collection and Preprocessing
We analyzed the dispersal patterns and survival rates of weaned elephant seals using tracking data collected from [insert data collection method and timeframe]. The seal tracking data was processed to provide locations at 12-hour intervals. Concurrently, we utilized particle tracking data to represent ocean currents in the study area. Both datasets were preprocessed and matched temporally and spatially using custom R functions.

Dispersal Analysis

We conducted a comprehensive dispersal analysis to investigate the relationship between seal movements and ocean currents. For each seal and corresponding particle track, we calculated bearings between each sequential location using the `geosphere` package in R. This provided us with a series of bearings representing the direction of movement for both seals and particles over time. We then computed cumulative circular correlations between these seal and particle bearings using the `circular` package to assess the similarity in movement patterns.

To identify critical periods in the seals' dispersal, we determined the maximum correlation and the time to reach this maximum for each individual. We used the median of these times to maximum correlation across all seals as the critical period. For each seal, we calculated the mean bearing during this critical period. We tested whether the mean maximum correlation differed significantly from zero using a one-sample t-test. Additionally, we examined the relationship between the time to maximum correlation and the maximum correlation value using a correlation test.

To account for individual variability and temporal patterns, we fitted a Generalized Additive Mixed Model (GAMM) using the `mgcv` package. The model included a smooth term for days since the start of the journey and a random effect for individual seals.

We compared each seal's mean bearing during the critical period to the overall mean particle bearing to determine if seals were following the prevailing currents. The overall mean particle bearing was calculated using all particle traces in the study area, providing a representation of the general current direction across the entire region. This approach was chosen to capture the broader oceanic flow patterns that might influence seal dispersal, rather than focusing on localized current variations. A seal was considered to be "following" if its mean bearing was within 45 degrees of this overall particle mean bearing. This classification method allowed us to assess whether seals were generally aligning their movement with the dominant current patterns in the region, and was later used as a factor in our survival analysis.

Survival Analysis

Survival data for the seals was obtained from [insert data source]. We conducted survival analyses for two periods: the first trip and the first year post-weaning. We used binomial Generalized Linear Models (GLMs) to investigate the relationship between survival and various predictors, including whether the seal was following currents, weaning mass, birth year, and environmental variables such as sea surface temperature, sea surface height anomaly, eddy kinetic energy, terrain ruggedness index, slope, sea surface temperature gradient, ice coverage, distance to ice edge, chlorophyll concentration, and chlorophyll gradient.

Model selection was performed using the dredge function from the MuMIn package, which fits all possible combinations of predictor variables. We considered models with a delta AICc (Akaike Information Criterion corrected for small sample sizes) of 2 or less to be the best-fitting models. We then used model averaging to account for model selection uncertainty and to obtain robust parameter estimates and relative variable importance.

Statistical Analysis and Visualization

All statistical analyses were conducted in R (version [insert version]). We used the `tidyverse` suite of packages for data manipulation and visualization, `sf` for spatial data handling, `patchwork` for combining multiple plots, and `rnaturalearth` for mapping. We generated various plots to visualize the dispersal patterns, including bearings over time, spatial trajectories, and cumulative correlations.

For the survival analysis, we used the `broom` package to extract and summarize model results. We calculated variable importance using the `sw` function from the MuMIn package.

All results, including summary statistics, test results, model outputs, and plots, were saved to a date-stamped output directory for reproducibility and further analysis.
